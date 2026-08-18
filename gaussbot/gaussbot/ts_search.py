"""
ts_search.py

Guesses a transition-state geometry -- interpolating between the
optimized reactant and product (geometry.interpolate_structures) if
the user doesn't supply one -- runs the TS optimization
(Opt=(TS,CalcFC,NoEigenTest) Freq, the "ts_highlevel" job type), and
checks that the resulting single imaginary frequency actually
corresponds to the reactant<->product interconversion rather than
some unrelated motion (e.g. a methyl rotation) before trusting it.

That check compares the imaginary mode's per-atom displacement vector
against the reactant->product atomic displacement vector, both
expressed in the *TS job's own frame* (align_to the reactant and
product onto the TS geometry itself -- not onto each other -- since
the mode vectors Gaussian prints are in whatever frame that
particular job used).
"""

from __future__ import annotations

import dataclasses
import math
import threading
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .executor import run_com
from .geometry import Structure, align_to, interpolate_structures
from .input_builder import build_input, write_com, BasisGroups
from .local_runner import check_normal_termination, GaussianRunError
from .parser import parse_log, imaginary_mode_displacement, LogResult


@dataclass
class TSOutcome:
    """What run_ts_search() came back with."""

    success: bool
    structure: Optional[Structure]
    result: Optional[LogResult]
    match_overlap: Optional[float]  # cosine similarity, mode vs. reaction vector
    log: List[str] = field(default_factory=list)
    imaginary_mode: Optional[List[Tuple[float, float, float]]] = None  # the TS's own imaginary-mode displacement vector
    log_path: Optional[str] = None  # the actual .log the accepted TS came from


def _reaction_vector(reactant: Structure, product: Structure, ts_frame: Structure):
    """Per-atom (dx, dy, dz) from reactant to product, both aligned onto
    `ts_frame` -- the frame the TS job's imaginary-mode vectors are in."""
    r = align_to(reactant, reference=ts_frame)
    p = align_to(product, reference=ts_frame)
    return [
        (px - rx, py - ry, pz - rz)
        for (_, rx, ry, rz), (_, px, py, pz) in zip(r.atoms, p.atoms)
    ]


def _min_pairwise_distance(structure: Structure) -> float:
    coords = [(x, y, z) for _, x, y, z in structure.atoms]
    if len(coords) < 2:
        return math.inf
    return min(
        math.dist(coords[i], coords[j]) for i in range(len(coords)) for j in range(i + 1, len(coords))
    )


def _cosine_similarity(a: List[tuple], b: List[tuple]) -> float:
    dot = sum(ax * bx + ay * by + az * bz for (ax, ay, az), (bx, by, bz) in zip(a, b))
    norm_a = math.sqrt(sum(ax**2 + ay**2 + az**2 for ax, ay, az in a))
    norm_b = math.sqrt(sum(bx**2 + by**2 + bz**2 for bx, by, bz in b))
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return dot / (norm_a * norm_b)


def run_ts_search(
    reactant: Structure,
    product: Structure,
    com_path: str,
    method: str,
    basis: str,
    basis_groups: Optional[BasisGroups] = None,
    ecp_groups: Optional[BasisGroups] = None,
    nprocs: int = 2,
    mem_gb: int = 2,
    ts_guess: Optional[Structure] = None,
    match_threshold: float = 0.4,
    max_guess_attempts: int = 3,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> TSOutcome:
    """
    Try to find the TS connecting reactant and product.

    If `ts_guess` is given, only that one guess is tried (the user
    presumably has a specific structure in mind). Otherwise, tries a
    few points along the naive linear-interpolation path between
    reactant and product (0.5, then 0.35, then 0.65) since a guess
    that lands on the wrong saddle point is a real possibility and
    the interpolation fraction is the one easy knob to turn.

    Each attempt: run the TS optimization, then require (a) it
    converged, (b) it has *exactly one* imaginary frequency (a true
    first-order saddle -- zero means it relaxed to a minimum, more
    than one means a higher-order saddle), and (c) that mode's
    displacement vector actually overlaps with the reactant->product
    motion above `match_threshold` -- otherwise it's a real TS, just
    for the wrong process (e.g. a methyl rotation), and the next
    guess is tried.
    """
    trail: List[str] = []

    if ts_guess is not None:
        guesses = [ts_guess]
    else:
        fractions = [0.5, 0.35, 0.65][:max_guess_attempts]
        guesses = [interpolate_structures(reactant, product, fraction=f) for f in fractions]

    for attempt, guess in enumerate(guesses, start=1):
        min_dist = _min_pairwise_distance(guess)
        if min_dist < 0.4:
            trail.append(
                f"attempt {attempt}: interpolated guess has two atoms only {min_dist:.2f} A "
                "apart -- a degenerate alignment, which shows up for near-linear "
                "reactant/product (the rotation about the molecular axis is undetermined, "
                "so Kabsch can pick a bad one). Skipping this guess; a hand-supplied "
                "ts_guess will do better for a reaction like this."
            )
            continue

        path = com_path if attempt == 1 else com_path.replace(".com", f"_try{attempt}.com")
        write_com(
            build_input(
                "ts_highlevel", [guess], method=method, basis=basis,
                basis_groups=basis_groups, ecp_groups=ecp_groups,
                nprocs=nprocs, mem_gb=mem_gb, chk=path[:-4],
            ),
            path,
        )

        try:
            log_path = run_com(path, executor=executor, cancel_event=cancel_event)
        except GaussianRunError as e:
            trail.append(f"attempt {attempt}: g09 error -- {e}")
            continue

        if not check_normal_termination(log_path):
            trail.append(f"attempt {attempt}: {log_path} did not terminate normally")
            continue

        result = parse_log(log_path)

        if not result.stationary_point_found:
            trail.append(f"attempt {attempt}: TS optimization didn't converge")
            continue
        if len(result.imaginary_freqs) == 0:
            trail.append(f"attempt {attempt}: relaxed to a minimum, not a TS (no imaginary frequency)")
            continue
        if len(result.imaginary_freqs) > 1:
            trail.append(
                f"attempt {attempt}: higher-order saddle point, not a true TS "
                f"({len(result.imaginary_freqs)} imaginary frequencies: {result.imaginary_freqs})"
            )
            continue

        ts_frame = dataclasses.replace(guess, atoms=result.final_geometry)
        mode = imaginary_mode_displacement(log_path)
        overlap = _cosine_similarity(mode, _reaction_vector(reactant, product, ts_frame)) if mode else 0.0

        if abs(overlap) < match_threshold:
            trail.append(
                f"attempt {attempt}: found a real TS (imaginary freq {result.imaginary_freqs[0]:.1f} cm^-1) "
                f"but its motion doesn't match the reactant->product conversion (overlap={overlap:.2f}) "
                "-- likely a spurious mode (e.g. a rotation), not the reaction coordinate. Trying a different guess."
            )
            continue

        trail.append(
            f"attempt {attempt}: TS found, imaginary freq {result.imaginary_freqs[0]:.1f} cm^-1, "
            f"overlap with reactant->product motion = {overlap:.2f} -- matches."
        )
        return TSOutcome(True, ts_frame, result, overlap, trail, imaginary_mode=mode, log_path=log_path)

    trail.append("Ran out of guesses without finding a TS whose imaginary mode matches the reaction.")
    return TSOutcome(False, None, None, None, trail)


def _attempt_one_ts_job(
    guess: Structure,
    com_path: str,
    job_type: str,
    method: str,
    basis: str,
    basis_groups: Optional[BasisGroups],
    ecp_groups: Optional[BasisGroups],
    nprocs: int,
    mem_gb: int,
    reactant: Optional[Structure],
    product: Optional[Structure],
    match_threshold: float,
    cancel_event: Optional[threading.Event],
    trail: List[str],
    attempt_label: str,
    executor: str,
) -> Optional[TSOutcome]:
    """
    One TS optimization attempt (job_type is "ts_highlevel" for CalcFC
    or "ts_calcall" for the brute-force fallback). Returns a validated
    TSOutcome, or None (with an explanatory trail line already
    appended) if this attempt didn't produce one. The reactant->product
    mode-overlap check only runs when both `reactant` and `product` are
    given -- with a TS-guess-only search there's nothing to compare
    the imaginary mode against, so this only requires a genuine
    first-order saddle (exactly one imaginary frequency).
    """
    min_dist = _min_pairwise_distance(guess)
    if min_dist < 0.4:
        trail.append(
            f"{attempt_label}: guess has two atoms only {min_dist:.2f} A apart -- "
            "a degenerate geometry, skipping."
        )
        return None

    write_com(
        build_input(
            job_type, [guess], method=method, basis=basis,
            basis_groups=basis_groups, ecp_groups=ecp_groups,
            nprocs=nprocs, mem_gb=mem_gb, chk=com_path[:-4],
        ),
        com_path,
    )

    try:
        log_path = run_com(com_path, executor=executor, cancel_event=cancel_event)
    except GaussianRunError as e:
        trail.append(f"{attempt_label}: g09 error -- {e}")
        return None

    if not check_normal_termination(log_path):
        trail.append(f"{attempt_label}: {log_path} did not terminate normally")
        return None

    result = parse_log(log_path)

    if not result.stationary_point_found:
        trail.append(f"{attempt_label}: TS optimization didn't converge")
        return None
    if len(result.imaginary_freqs) == 0:
        trail.append(f"{attempt_label}: relaxed to a minimum, not a TS (no imaginary frequency)")
        return None
    if len(result.imaginary_freqs) > 1:
        trail.append(
            f"{attempt_label}: higher-order saddle point, not a true TS "
            f"({len(result.imaginary_freqs)} imaginary frequencies: {result.imaginary_freqs})"
        )
        return None

    ts_frame = dataclasses.replace(guess, atoms=result.final_geometry)
    mode = imaginary_mode_displacement(log_path)

    overlap: Optional[float] = None
    if reactant is not None and product is not None:
        overlap = _cosine_similarity(mode, _reaction_vector(reactant, product, ts_frame)) if mode else 0.0
        if abs(overlap) < match_threshold:
            trail.append(
                f"{attempt_label}: found a real TS (imaginary freq {result.imaginary_freqs[0]:.1f} cm^-1) "
                f"but its motion doesn't match the reactant->product conversion (overlap={overlap:.2f}) "
                "-- likely a spurious mode (e.g. a rotation), not the reaction coordinate."
            )
            return None

    overlap_msg = (
        f", overlap with reactant->product motion = {overlap:.2f}" if overlap is not None
        else " (no reactant/product given -- motion not checked against a reaction vector)"
    )
    trail.append(
        f"{attempt_label}: TS found, imaginary freq {result.imaginary_freqs[0]:.1f} cm^-1{overlap_msg}"
    )
    return TSOutcome(True, ts_frame, result, overlap, list(trail), imaginary_mode=mode, log_path=log_path)


def _search_ts_at_level(
    guesses: List[Structure],
    com_path_base: str,
    method: str,
    basis: str,
    basis_groups: Optional[BasisGroups],
    ecp_groups: Optional[BasisGroups],
    nprocs: int,
    mem_gb: int,
    reactant: Optional[Structure],
    product: Optional[Structure],
    match_threshold: float,
    cancel_event: Optional[threading.Event],
    trail: List[str],
    use_calcall_fallback: bool,
    executor: str,
) -> Optional[TSOutcome]:
    """
    Try each guess in turn with the CalcFC route (ts_highlevel). If
    none of them lands a validated TS and `use_calcall_fallback` is
    set, retry once more on the primary guess (the first one -- either
    the fixed ts_guess, or the 0.5 interpolation fraction) with the
    CalcAll route (ts_calcall) -- much more likely to actually locate
    the saddle point when CalcFC can't, but also much more expensive
    (full Hessian every step instead of once), so it's opt-in rather
    than automatic.
    """
    for i, guess in enumerate(guesses, start=1):
        com_path = f"{com_path_base}.com" if i == 1 else f"{com_path_base}_try{i}.com"
        outcome = _attempt_one_ts_job(
            guess, com_path, "ts_highlevel", method, basis, basis_groups, ecp_groups,
            nprocs, mem_gb, reactant, product, match_threshold, cancel_event, trail,
            f"CalcFC attempt {i}", executor,
        )
        if outcome is not None:
            return outcome

    if not use_calcall_fallback:
        trail.append("Every CalcFC attempt failed -- CalcAll fallback is off, not trying it.")
        return None

    calcall_path = f"{com_path_base}_calcall.com"
    return _attempt_one_ts_job(
        guesses[0], calcall_path, "ts_calcall", method, basis, basis_groups, ecp_groups,
        nprocs, mem_gb, reactant, product, match_threshold, cancel_event, trail,
        "CalcAll fallback", executor,
    )


def search_ts_staged(
    out_dir: str,
    method: str,
    basis: str,
    basis_groups: Optional[BasisGroups] = None,
    ecp_groups: Optional[BasisGroups] = None,
    nprocs: int = 2,
    mem_gb: int = 2,
    reactant: Optional[Structure] = None,
    product: Optional[Structure] = None,
    ts_guess: Optional[Structure] = None,
    skip_pm6: bool = False,
    use_calcall_fallback: bool = False,
    calcfc_attempts: int = 4,
    match_threshold: float = 0.4,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> TSOutcome:
    """
    A more resilient TS search than run_ts_search(): tries a cheap PM6
    stage first (unless `skip_pm6`) before ever touching the expensive
    final method/basis, and -- only if `use_calcall_fallback` is set,
    off by default since it's a much more expensive full-Hessian-every-
    step calculation -- falls back to a CalcAll brute-force attempt at
    each stage if the usual CalcFC route can't locate a saddle point.
    Kept entirely separate from run_ts_search() (used by the existing
    reaction-mechanism study) rather than refactoring it -- this is
    additive, for the dedicated TS Search section only.

    At least one of `ts_guess`, or both `reactant` and `product`,
    must be given. With `ts_guess`: that's the one fixed starting
    geometry tried at each stage. Without it: the same 0.5/0.35/0.65
    interpolation-fraction ladder run_ts_search() uses, generated fresh
    from `reactant`/`product` (capped to `calcfc_attempts` fractions).
    The reactant->product mode-overlap validation only runs when both
    `reactant` and `product` are given -- with a guess-only search,
    "exactly one imaginary frequency" (a genuine first-order saddle) is
    the only check available.

    If the PM6 stage finds a validated TS, its converged structure
    becomes the single fixed starting guess for the final level
    (one less thing for the expensive stage to get wrong); otherwise
    the final level searches from the original guess(es) itself.
    """
    if ts_guess is None and (reactant is None or product is None):
        raise ValueError("search_ts_staged needs either ts_guess, or both reactant and product")
    calcfc_attempts = max(1, calcfc_attempts)

    trail: List[str] = []

    def _build_guesses() -> List[Structure]:
        if ts_guess is not None:
            return [ts_guess]
        fractions = [0.5, 0.35, 0.65][:calcfc_attempts]
        return [interpolate_structures(reactant, product, fraction=f) for f in fractions]

    fixed_guess: Optional[Structure] = None
    if not skip_pm6:
        trail.append("--- PM6 stage ---")
        pm6_outcome = _search_ts_at_level(
            _build_guesses(), f"{out_dir}/ts_pm6", "PM6", "", None, None,
            nprocs, mem_gb, reactant, product, match_threshold, cancel_event, trail,
            use_calcall_fallback, executor,
        )
        if pm6_outcome is not None:
            trail.append(
                "PM6 stage found a validated TS -- using its converged geometry as the "
                "starting guess for the final level."
            )
            fixed_guess = pm6_outcome.structure
        else:
            trail.append(
                "PM6 stage couldn't find a validated TS -- the final level will search "
                "from the original guess(es) instead."
            )
    else:
        trail.append("Skipping the PM6 stage -- searching directly at the final level.")

    trail.append("--- Final level ---")
    final_guesses = [fixed_guess] if fixed_guess is not None else _build_guesses()
    final_outcome = _search_ts_at_level(
        final_guesses, f"{out_dir}/ts_final", method, basis, basis_groups, ecp_groups,
        nprocs, mem_gb, reactant, product, match_threshold, cancel_event, trail,
        use_calcall_fallback, executor,
    )
    if final_outcome is not None:
        return final_outcome

    fallback_desc = "CalcFC/CalcAll" if use_calcall_fallback else "CalcFC (CalcAll fallback was off)"
    trail.append(
        f"Ran out of guesses and fallbacks (PM6 {fallback_desc}, final-level {fallback_desc}) "
        "without finding a validated TS."
    )
    return TSOutcome(False, None, None, None, trail)
