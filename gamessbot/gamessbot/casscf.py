"""
casscf.py

Runs the CASSCF stage on top of a converged RHF and CIS: the CIS
excited states drive active_space.suggest_active_space() (NMCC/NDOC/
NVAL), which gamess_input.build_casscf_input() turns into a full
FORS-CASSCF input (state-averaged over `nstate` states) run through the
same local_runner.run_gamess() everything else uses.

Convergence isn't guaranteed on the first try (FULLNR can struggle,
especially for small/near-degenerate active spaces) -- run_casscf_staged()
automates the recovery ladder you described: retry with more iterations,
and if that still doesn't converge, fall back to a smaller/easier active
space first and use its converged orbitals as a much better starting
guess to reach the originally requested (larger) active space.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .active_space import ActiveSpaceSuggestion, regrow_active_space, shrink_active_space, suggest_active_space
from .cis import CISOutcome
from .executor import run_gms
from .gamess_input import build_casscf_input, extract_optimized_mcscf_vec_annotation, extract_optimized_mcscf_vec_block
from .local_runner import check_normal_termination
from .parser import CASSCFResult, parse_casscf_log, parse_rhf_log
from .rhf import RHFOutcome


@dataclass
class CASSCFOutcome:
    """What run_casscf() came back with. `success` means the MCSCF
    orbital optimization actually converged (result.converged), not
    just that GAMESS exited without crashing -- a run can terminate
    normally after exhausting MAXIT without converging, printing a
    placeholder "FINAL MCSCF ENERGY IS 0.0000000000" that isn't real;
    see parser.parse_casscf_log().

    data_block/charge/mult/gbasis_line/norb are carried over from the
    RHF outcome and optimized_vec_block (the punch's "OPTIMIZED MCSCF
    MO-S" block, only set on success) is added here, so xmcqdpt.run_xmcqdpt()
    (or a later CASSCF attempt restarting from this one) can build on
    top of this outcome alone, without needing the original rhf_outcome
    passed again."""

    success: bool
    log_path: str
    result: Optional[CASSCFResult]
    active_space: ActiveSpaceSuggestion
    data_block: str
    charge: int
    mult: int
    gbasis_line: str
    norb: int
    optimized_vec_block: Optional[str] = None
    optimized_vec_annotation: str = ""  # the "--- OPTIMIZED MCSCF MO-S --- ..." heading, if any


@dataclass
class CASSCFStagedOutcome:
    """What run_casscf_staged() came back with."""

    outcome: CASSCFOutcome
    trail: List[str] = field(default_factory=list)
    # True if every automatic attempt (MAXIT=120, then MAXIT=200) failed to converge --
    # the caller may offer the smaller-active-space recovery path (run_casscf_with_smaller_active_space_recovery())
    exhausted: bool = False


def mo_source_from_casscf_outcome(previous: CASSCFOutcome) -> RHFOutcome:
    """Adapts a converged CASSCFOutcome into the RHFOutcome-shaped
    object run_casscf()/run_casscf_staged() expect as a starting-
    orbitals source -- its optimized_vec_block/optimized_vec_annotation
    stand in for vec_block/vec_annotation, so a later combination can
    start its own orbital optimization from these converged (possibly
    MO-reordered) orbitals instead of the original closed-shell ones,
    per your request to make that a choice rather than always
    restarting from the closed-shell orbitals."""
    return RHFOutcome(
        success=True, log_path=previous.log_path, energy=None, norb=previous.norb,
        vec_block=previous.optimized_vec_block, data_block=previous.data_block,
        charge=previous.charge, mult=previous.mult, gbasis_line=previous.gbasis_line,
        vec_annotation=previous.optimized_vec_annotation,
    )


def suggest_active_space_from_cis(
    cis_outcome: CISOutcome,
    rhf_outcome: RHFOutcome,
    threshold: float = 0.20,
    max_electrons: int = 16,
    max_orbitals: int = 16,
) -> ActiveSpaceSuggestion:
    """
    Convenience wrapper: reads n_occ (the number of doubly-occupied
    orbitals -- the NMCC/NDOC boundary) straight off the RHF log, so
    callers don't have to parse it themselves before calling
    active_space.suggest_active_space().
    """
    rhf_result = parse_rhf_log(rhf_outcome.log_path)
    return suggest_active_space(
        cis_outcome.states, n_occ=rhf_result.n_occ, norb=rhf_outcome.norb,
        threshold=threshold, max_electrons=max_electrons, max_orbitals=max_orbitals,
    )


def _run_casscf_job(
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    vec_block: str,
    norb: int,
    active_space: ActiveSpaceSuggestion,
    out_dir: str,
    name: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int,
    mem_mwords: int,
    wstate: Optional[List[float]],
    maxit: int,
    cancel_event: Optional[threading.Event],
    executor: str = "local",
    vec_annotation: str = "",
) -> CASSCFOutcome:
    """The actual mechanics of one CASSCF attempt -- builds+runs+parses
    a single .inp named `<name>.inp` (so multiple attempts in the same
    out_dir, e.g. from run_casscf_staged(), don't overwrite each
    other's files)."""
    inp_text = build_casscf_input(
        data_block, charge=charge, mult=mult, gbasis_line=gbasis_line, vec_block=vec_block, norb=norb,
        active_space=active_space, nstate=nstate, wstate=wstate, mem_mwords=mem_mwords, maxit=maxit,
        vec_annotation=vec_annotation,
    )
    inp_path = os.path.join(out_dir, f"{name}.inp")
    with open(inp_path, "w") as f:
        f.write(inp_text)

    log_path = run_gms(inp_path, rungms_path, scratch_dir, ncpus=ncpus, executor=executor, cancel_event=cancel_event)

    ran_cleanly = check_normal_termination(log_path)
    result = parse_casscf_log(log_path, nstate=nstate) if ran_cleanly else None
    success = ran_cleanly and result.converged

    optimized_vec_block = None
    optimized_vec_annotation = ""
    if success:
        dat_path = os.path.join(out_dir, f"{name}.dat")
        with open(dat_path) as f:
            dat_text = f.read()
        optimized_vec_block = extract_optimized_mcscf_vec_block(dat_text)
        optimized_vec_annotation = extract_optimized_mcscf_vec_annotation(dat_text)

    return CASSCFOutcome(
        success=success, log_path=log_path, result=result, active_space=active_space,
        data_block=data_block, charge=charge, mult=mult, gbasis_line=gbasis_line, norb=norb,
        optimized_vec_block=optimized_vec_block, optimized_vec_annotation=optimized_vec_annotation,
    )


def run_casscf(
    rhf_outcome: RHFOutcome,
    active_space: ActiveSpaceSuggestion,
    out_dir: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    mem_mwords: int = 1,
    wstate: Optional[List[float]] = None,
    maxit: int = 120,
    name: str = "casscf",
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> CASSCFOutcome:
    """
    A single CASSCF attempt on top of `rhf_outcome` -- its data_block/
    charge/mult/gbasis_line/norb/vec_block feed the CASSCF input's
    $DATA/$BASIS/$GUESS directly, same as cis.run_cis(). `active_space`
    is typically suggest_active_space_from_cis()'s result, possibly
    adjusted by the user first (see cli.py/webapp.py's confirm-or-edit
    step). For the full retry/recovery ladder, see run_casscf_staged().
    """
    if not rhf_outcome.success:
        raise ValueError("run_casscf() needs a successful RHF outcome to read orbitals from")
    return _run_casscf_job(
        rhf_outcome.data_block, rhf_outcome.charge, rhf_outcome.mult, rhf_outcome.gbasis_line,
        rhf_outcome.vec_block, rhf_outcome.norb, active_space, out_dir, name, nstate,
        rungms_path, scratch_dir, ncpus, mem_mwords, wstate, maxit, cancel_event, executor,
        vec_annotation=rhf_outcome.vec_annotation,
    )


def run_casscf_staged(
    rhf_outcome: RHFOutcome,
    active_space: ActiveSpaceSuggestion,
    out_dir: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    mem_mwords: int = 1,
    wstate: Optional[List[float]] = None,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
    name_prefix: str = "casscf",
) -> CASSCFStagedOutcome:
    """
    First attempt at MAXIT=120 (the default); if that doesn't converge,
    a second attempt at MAXIT=200, logging both in the trail either way.
    If BOTH fail, `exhausted` is set so the caller can offer the
    smaller-active-space recovery path (run_casscf_with_smaller_active_space_recovery())
    instead of just giving up.

    `name_prefix` (default "casscf", giving "casscf_try1.inp"/
    "casscf_try2.inp") lets multiple active-space/state combinations
    coexist in the same out_dir without overwriting each other's files
    -- e.g. "cas-44-sa3" for a CAS(4,4), 3-state combo.
    """
    trail: List[str] = []

    outcome = run_casscf(
        rhf_outcome, active_space, out_dir, nstate, rungms_path, scratch_dir,
        ncpus=ncpus, mem_mwords=mem_mwords, wstate=wstate, maxit=120, name=f"{name_prefix}_try1",
        cancel_event=cancel_event, executor=executor,
    )
    if outcome.success:
        trail.append("CASSCF (MAXIT=120): converged")
        return CASSCFStagedOutcome(outcome=outcome, trail=trail, exhausted=False)
    trail.append("CASSCF (MAXIT=120): did not converge -- retrying with MAXIT=200")

    outcome2 = run_casscf(
        rhf_outcome, active_space, out_dir, nstate, rungms_path, scratch_dir,
        ncpus=ncpus, mem_mwords=mem_mwords, wstate=wstate, maxit=200, name=f"{name_prefix}_try2",
        cancel_event=cancel_event, executor=executor,
    )
    if outcome2.success:
        trail.append("CASSCF (MAXIT=200): converged")
        return CASSCFStagedOutcome(outcome=outcome2, trail=trail, exhausted=False)
    trail.append("CASSCF (MAXIT=200): still did not converge")
    return CASSCFStagedOutcome(outcome=outcome2, trail=trail, exhausted=True)


def run_casscf_with_smaller_active_space_recovery(
    rhf_outcome: RHFOutcome,
    original_active_space: ActiveSpaceSuggestion,
    out_dir: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    smaller_max_electrons: int,
    smaller_max_orbitals: int,
    ncpus: int = 1,
    mem_mwords: int = 1,
    wstate: Optional[List[float]] = None,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
    name_prefix: str = "casscf",
) -> Tuple[CASSCFOutcome, ActiveSpaceSuggestion, List[str]]:
    """
    The recovery path once run_casscf_staged() is exhausted: shrinks
    `original_active_space` down to (smaller_max_electrons,
    smaller_max_orbitals) using the same CIS scores it was originally
    ranked by, converges CASSCF there first (usually much easier), then
    restarts the ORIGINAL (larger) active space using those converged
    orbitals as a starting guess (active_space.regrow_active_space()) --
    per your own description of this recovery strategy.

    Returns (final_outcome, active_space_used_for_final_outcome, trail).
    If the smaller pass itself doesn't converge either, `final_outcome`
    is that failed smaller-pass outcome (nothing further to try) and
    the returned active space is the smaller one, not the original.

    `name_prefix` -- see run_casscf_staged().
    """
    trail: List[str] = []

    smaller_active_space = shrink_active_space(original_active_space, smaller_max_electrons, smaller_max_orbitals)
    trail.append(
        f"Trying a smaller active space first: NMCC={smaller_active_space.nmcc} "
        f"NDOC={smaller_active_space.ndoc} NVAL={smaller_active_space.nval} "
        f"(CAS({2 * smaller_active_space.ndoc},{smaller_active_space.ndoc + smaller_active_space.nval}))"
    )
    smaller_outcome = run_casscf(
        rhf_outcome, smaller_active_space, out_dir, nstate, rungms_path, scratch_dir,
        ncpus=ncpus, mem_mwords=mem_mwords, wstate=wstate, maxit=120, name=f"{name_prefix}_smaller",
        cancel_event=cancel_event, executor=executor,
    )
    if not smaller_outcome.success:
        trail.append("Smaller active space did not converge either -- no further automatic recovery to try")
        return smaller_outcome, smaller_active_space, trail
    trail.append("Smaller active space converged -- using its orbitals to restart the original active space")

    regrown_active_space = regrow_active_space(original_active_space, smaller_active_space)
    regrow_outcome = _run_casscf_job(
        smaller_outcome.data_block, smaller_outcome.charge, smaller_outcome.mult, smaller_outcome.gbasis_line,
        smaller_outcome.optimized_vec_block, smaller_outcome.norb, regrown_active_space, out_dir,
        f"{name_prefix}_regrown", nstate, rungms_path, scratch_dir, ncpus, mem_mwords, wstate, 120, cancel_event, executor,
        vec_annotation=smaller_outcome.optimized_vec_annotation,
    )
    if regrow_outcome.success:
        trail.append("Regrown to the original active space: converged")
    else:
        trail.append("Regrown to the original active space: still did not converge")
    return regrow_outcome, regrown_active_space, trail
