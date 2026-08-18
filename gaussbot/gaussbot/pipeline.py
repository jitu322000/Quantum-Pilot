"""
pipeline.py

Chains the low-level pieces (geometry, input_builder, local_runner,
parser) into the resilient optimize-and-verify loop described by the
project's actual working intuition: run an Opt Freq, and if it comes
back as anything other than a clean minimum, try to repair it rather
than just giving up on the first bad log --

  - didn't converge to a stationary point -> retry with more Opt
    cycles, continuing from wherever the optimizer actually got to
    (plus a small random jitter) rather than restarting from the
    exact same starting guess every time -- just adding cycles on an
    unchanged geometry can walk right back into the same wall.
  - converged, but has an imaginary frequency -> nudge the geometry
    along that mode and reoptimize from there.

and, one level up, if PM6 can't be coaxed into a minimum at all from
the geometry it was given -> refetch a geometry from PubChem and
retry PM6 -> if that still doesn't work -> escalate to HF/STO-3G.
"""

from __future__ import annotations

import dataclasses
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from .executor import run_com
from .geometry import Structure, displace_structure, jitter_structure, from_pubchem
from .input_builder import build_input, write_com, BasisGroups
from .local_runner import check_normal_termination, GaussianRunError
from .parser import parse_log, imaginary_mode_displacement, last_geometry, last_convergence_status, LogResult, ConvergenceStatus

# Small nudge applied whenever a cycle-bump retry continues from the
# last geometry an unconverged optimization actually reached, rather
# than restarting from the exact same guess -- enough to break an
# oscillation/stall, small enough not to disturb genuine progress. Used
# when _repair_strategy() says "continue" (most convergence criteria
# already met); "escalate" scales this up instead, see below.
CYCLE_BUMP_JITTER = 0.02

# Cap on how far the escalating jitter is allowed to grow -- enough to
# genuinely leave a stuck region without being violent enough to blow
# up a reasonable starting structure.
MAX_ESCALATE_JITTER = 0.5


def _repair_strategy(status: Optional[ConvergenceStatus]) -> str:
    """
    Decide how to react to a non-converged optimization based on
    Gaussian's own convergence table at the point it gave up, instead
    of always reacting the same way regardless of how close it
    actually got:

      - "continue": at least 3 of the 4 criteria (Maximum/RMS Force,
        Maximum/RMS Displacement) are already satisfied -- probably
        just needs more cycles from where it got to, not a different
        starting point.
      - "escalate": fewer than that are met -- the small nudge this
        loop already tried clearly isn't enough on its own, so use a
        bigger one instead of repeating the same thing and expecting a
        different result.

    Falls back to "continue" (today's old behavior) if there's no
    convergence data to go on at all -- e.g. a log that crashed before
    any optimization step printed one.
    """
    if status is None or status.n_total == 0:
        return "continue"
    return "continue" if status.n_converged >= status.n_total - 1 else "escalate"


def _convergence_summary(status: Optional[ConvergenceStatus]) -> str:
    """Short human-readable summary for the trail log, e.g.
    '3/4 criteria met (only RMS Displacement short: 0.0013 vs 0.0012 threshold)'
    or '0/4 criteria met, nothing trending toward convergence'."""
    if status is None or status.n_total == 0:
        return "no convergence data available"
    unmet = [(name, item) for name, item in status.items if not item.converged]
    if not unmet:
        return f"{status.n_converged}/{status.n_total} criteria met"
    if len(unmet) == 1:
        name, item = unmet[0]
        return f"{status.n_converged}/{status.n_total} criteria met (only {name} short: {item.value:g} vs {item.threshold:g} threshold)"
    return f"{status.n_converged}/{status.n_total} criteria met, nothing trending toward convergence"


@dataclass
class OptOutcome:
    """What repair_and_optimize()/preopt_with_escalation() came back with."""

    success: bool
    structure: Optional[Structure]
    result: Optional[LogResult]
    method: str
    basis: str
    log: List[str] = field(default_factory=list)  # human-readable trail of what was tried
    log_path: Optional[str] = None  # the actual .log that produced `structure`/`result` (may be a _tryN attempt)


def repair_and_optimize(
    structure: Structure,
    com_path: str,
    method: str = "PM6",
    basis: str = "",
    basis_groups: Optional[BasisGroups] = None,
    ecp_groups: Optional[BasisGroups] = None,
    nprocs: int = 2,
    mem_gb: int = 2,
    initial_maxcycles: Optional[int] = None,
    max_repair_attempts: int = 4,
    displacement_step: float = 0.3,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> OptOutcome:
    """
    Run a {method}/{basis} Opt Freq on `structure`, repairing in place
    when the result isn't a clean minimum instead of just failing:

      - non-convergence -> looks at Gaussian's own convergence table
        (parser.last_convergence_status) at the point it gave up, not
        just the bare fact that it didn't converge. If most of the 4
        criteria (Force/Displacement, max/RMS) are already met, that's
        a "close, just needs more cycles" situation -- continue from
        the last geometry reached with MaxCycles doubled and a small
        jitter, same as before. If few or none are met, the small
        nudge clearly isn't enough on its own -- use a bigger one
        instead of repeating something that already didn't work (see
        _repair_strategy/_convergence_summary). Up to
        `max_repair_attempts` times either way.
      - an imaginary frequency -> nudge the geometry along that mode
        (imaginary_mode_displacement + displace_structure) and
        reoptimize, up to `max_repair_attempts` times.

    Each retry is written to its own .com/.log (com_path with a
    "_tryN" suffix) so every attempt stays on disk for inspection.
    Gives up once both repair budgets are spent, or on a g09/error
    termination there's no way to repair from.

    `executor` ("local" or "pbs", default "local") picks how each
    attempt actually runs -- see executor.run_com().
    """
    trail: List[str] = []
    current = structure
    maxcycles = initial_maxcycles

    attempt = 0
    cycle_bumps_used = 0
    escalates_used = 0
    displacements_used = 0

    def _bump(log_path_for_status: str) -> tuple:
        """Decide strategy from the log's convergence table and return
        (jitter_scale, strategy_note) -- also bumps escalates_used."""
        nonlocal escalates_used
        status = last_convergence_status(log_path_for_status)
        strategy = _repair_strategy(status)
        summary = _convergence_summary(status)
        if strategy == "continue":
            return CYCLE_BUMP_JITTER, f"{summary} -- continuing with more cycles"
        jitter = min(CYCLE_BUMP_JITTER * (2 ** (escalates_used + 1)), MAX_ESCALATE_JITTER)
        escalates_used += 1
        return jitter, f"{summary} -- using a bigger kick ({jitter:.2f} A) instead of repeating the same nudge"

    while True:
        attempt += 1
        path = com_path if attempt == 1 else com_path.replace(".com", f"_try{attempt}.com")
        write_com(
            build_input(
                "pm6_opt_freq",
                [current],
                method=method,
                basis=basis,
                basis_groups=basis_groups,
                ecp_groups=ecp_groups,
                nprocs=nprocs,
                mem_gb=mem_gb,
                opt_maxcycles=maxcycles,
                chk=path[:-4],  # strip ".com" -- keeps the checkpoint unique per attempt and next to its .com/.log, not a shared generic name in cwd
            ),
            path,
        )

        try:
            log_path = run_com(path, executor=executor, cancel_event=cancel_event)
        except GaussianRunError as e:
            # Hitting MaxCycles is itself an *error* termination in
            # Gaussian ("Optimization stopped -- Number of steps
            # exceeded"), not a normal one -- that's a recoverable
            # "needs more cycles" signal, not a real failure.
            log_path = str(Path(path).with_suffix(".log"))
            log_text = ""
            try:
                with open(log_path) as f:
                    log_text = f.read()
            except OSError:
                pass

            if "Number of steps exceeded" in log_text:
                if cycle_bumps_used >= max_repair_attempts:
                    trail.append(
                        f"attempt {attempt}: still hasn't converged after "
                        f"{cycle_bumps_used} repair attempts -- giving up"
                    )
                    return OptOutcome(False, current, None, method, basis, trail, log_path)
                cycle_bumps_used += 1
                maxcycles = (maxcycles or 100) * 2
                jitter_scale, note = _bump(log_path)
                partial = last_geometry(log_path)
                if partial is not None:
                    current = jitter_structure(dataclasses.replace(current, atoms=partial), scale=jitter_scale)
                    trail.append(
                        f"attempt {attempt}: hit the step limit without converging -- {note}, "
                        f"retrying with MaxCycles={maxcycles}"
                    )
                else:
                    trail.append(
                        f"attempt {attempt}: hit the step limit without converging -- couldn't "
                        f"read the partial geometry, retrying the same guess with MaxCycles={maxcycles}"
                    )
                continue

            trail.append(f"attempt {attempt} ({method}/{basis or 'no basis'}): g09 error -- {e}")
            return OptOutcome(False, None, None, method, basis, trail, log_path)

        if not check_normal_termination(log_path):
            trail.append(f"attempt {attempt}: {log_path} did not terminate normally -- giving up")
            return OptOutcome(False, current, None, method, basis, trail, log_path)

        result = parse_log(log_path)

        if not result.stationary_point_found:
            if cycle_bumps_used >= max_repair_attempts:
                trail.append(
                    f"attempt {attempt}: still hasn't converged after "
                    f"{cycle_bumps_used} repair attempts -- giving up"
                )
                return OptOutcome(False, current, result, method, basis, trail, log_path)
            cycle_bumps_used += 1
            maxcycles = (maxcycles or 100) * 2
            jitter_scale, note = _bump(log_path)
            if result.final_geometry is not None:
                current = jitter_structure(dataclasses.replace(current, atoms=result.final_geometry), scale=jitter_scale)
                trail.append(
                    f"attempt {attempt}: didn't converge -- {note}, retrying with MaxCycles={maxcycles}"
                )
            else:
                trail.append(f"attempt {attempt}: didn't converge -- retrying with MaxCycles={maxcycles}")
            continue

        if not result.frequencies:
            trail.append(f"attempt {attempt}: no frequency data in {log_path} -- giving up")
            return OptOutcome(False, current, result, method, basis, trail, log_path)

        if result.imaginary_freqs:
            if displacements_used >= max_repair_attempts:
                trail.append(
                    f"attempt {attempt}: still imaginary after {displacements_used} "
                    f"displacement attempts (freqs={result.imaginary_freqs}) -- giving up"
                )
                return OptOutcome(False, current, result, method, basis, trail, log_path)
            vectors = imaginary_mode_displacement(log_path)
            if vectors is None:
                trail.append(
                    f"attempt {attempt}: imaginary frequency but couldn't read its "
                    "displacement vector -- giving up"
                )
                return OptOutcome(False, current, result, method, basis, trail, log_path)
            displacements_used += 1
            current = displace_structure(current, vectors, scale=displacement_step)
            trail.append(
                f"attempt {attempt}: imaginary frequency {result.imaginary_freqs} -- "
                "nudging geometry along that mode and reoptimizing"
            )
            continue

        current = dataclasses.replace(current, atoms=result.final_geometry)
        trail.append(f"attempt {attempt}: converged, clean minimum, SCF energy = {result.scf_energy}")
        return OptOutcome(True, current, result, method, basis, trail, log_path)


def preopt_with_escalation(
    structure: Structure,
    out_dir: str,
    tag: str,
    pubchem_query: Optional[str] = None,
    nprocs: int = 2,
    mem_gb: int = 2,
    max_repair_attempts: int = 4,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> OptOutcome:
    """
    Get *some* trustworthy stationary-point geometry to hand off to the
    higher-level stage, trying progressively more drastic fallbacks:

      1. PM6 on the geometry as given.
      2. If that fails and `pubchem_query` was given, refetch the
         geometry from PubChem (same charge/multiplicity) and retry
         PM6 from there.
      3. If PM6 still can't do it, escalate to HF/STO-3G on the best
         geometry found so far.
    """
    trail: List[str] = []

    outcome = repair_and_optimize(
        structure, f"{out_dir}/{tag}_pm6opt.com", method="PM6", nprocs=nprocs, mem_gb=mem_gb,
        max_repair_attempts=max_repair_attempts, cancel_event=cancel_event, executor=executor,
    )
    trail += [f"[PM6, original geometry] {line}" for line in outcome.log]
    if outcome.success:
        outcome.log = trail
        return outcome

    best_geometry = outcome.structure or structure

    if pubchem_query:
        try:
            pubchem_structure = from_pubchem(pubchem_query, label=structure.label)
            pubchem_structure.charge = structure.charge
            pubchem_structure.multiplicity = structure.multiplicity
        except Exception as e:
            trail.append(f"[PubChem lookup for {pubchem_query!r}] failed: {e}")
            pubchem_structure = None

        if pubchem_structure is not None:
            outcome = repair_and_optimize(
                pubchem_structure, f"{out_dir}/{tag}_pm6opt_pubchem.com", method="PM6",
                nprocs=nprocs, mem_gb=mem_gb, max_repair_attempts=max_repair_attempts,
                cancel_event=cancel_event, executor=executor,
            )
            trail += [f"[PM6, PubChem geometry] {line}" for line in outcome.log]
            if outcome.success:
                outcome.log = trail
                return outcome
            best_geometry = outcome.structure or best_geometry

    outcome = repair_and_optimize(
        best_geometry, f"{out_dir}/{tag}_hf_sto3g_opt.com", method="HF", basis="STO-3G",
        nprocs=nprocs, mem_gb=mem_gb, max_repair_attempts=max_repair_attempts,
        cancel_event=cancel_event, executor=executor,
    )
    trail += [f"[HF/STO-3G] {line}" for line in outcome.log]
    outcome.log = trail
    return outcome
