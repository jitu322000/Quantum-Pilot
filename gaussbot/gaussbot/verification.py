"""
verification.py

Two independent, opt-in cross-checks that regenerate the reactant and
product straight from the verified TS, instead of trusting the
original guess-based reactant/product optimization alone:

  - TS-mode distortion: displace the TS's converged geometry along its
    own imaginary vibrational mode, in the + and - directions, and
    reoptimize each to a nearby ground-state minimum. One side should
    relax toward the reactant, the other toward the product -- the
    same technique as the attached dis.f90 (Gaussian normal-mode
    distortion for GAMESS inputs), applied here to reoptimize in
    Gaussian directly rather than hand it off to another package. The
    underlying displacement math is the same thing pipeline.py's
    repair loop already uses (parser.imaginary_mode_displacement() +
    geometry.displace_structure()) -- Gaussian's own per-mode
    "Atom AN X Y Z" table is already a ready-to-add Cartesian
    displacement, not mass-weighted, so no manual unit conversion is
    needed the way dis.f90 does for GAMESS's raw normal-coordinate
    output.
  - IRC-endpoint reoptimization: an IRC path's forward/reverse
    endpoints are just the last point of a fixed-step path, not
    necessarily a fully converged minimum -- reoptimizing each gives a
    third candidate per side.

Either way, the point is the same: collect every candidate geometry
that exists for the reactant (and, separately, the product), and let
select_best_candidate() pick whichever is actually lowest in energy
for the final barrier/reaction-energy report -- a TS is only as
trustworthy as the reactant/product it's shown to actually connect.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

from .geometry import Structure, displace_structure, rmsd
from .irc import IRCOutcome
from .level_select import LevelChoice
from .pipeline import repair_and_optimize, OptOutcome
from .ts_search import TSOutcome

# Kept close together (< this) counts as "clearly this side" when
# classifying a reoptimized candidate as reactant- or product-derived.
# Unlike irc.py's RMSD_MATCH_THRESHOLD (which gates pass/fail), this
# only picks a *label* -- a candidate that's ambiguous either way is
# still reported, just flagged as such rather than silently guessed.
CLASSIFY_AMBIGUOUS_MARGIN = 0.3  # Angstrom


@dataclass
class SideCandidate:
    """One reoptimized candidate structure, tagged with how confidently
    it was classified reactant-side vs. product-side."""

    outcome: OptOutcome
    rmsd_to_reactant: float
    rmsd_to_product: float


@dataclass
class DistortionOutcome:
    success: bool
    reactant_candidate: Optional[SideCandidate]
    product_candidate: Optional[SideCandidate]
    log: List[str] = field(default_factory=list)


def distort_along_imaginary_mode(
    ts: Structure, mode: List[Tuple[float, float, float]], factor: float = 1.0
) -> Tuple[Structure, Structure]:
    """
    (plus, minus) structures displaced off `ts` along `mode` by
    `factor` -- x = x0 +/- factor*mode, same convention as dis.f90.
    `mode` is already a ready-to-add Cartesian displacement vector
    (see parser.imaginary_mode_displacement), so factor=1.0 is a full,
    undamped step -- deliberately, since this is displacing a
    converged TS to find the minima on either side of it, not nudging
    a stalled optimization the way pipeline.py's repair loop does.
    """
    return (
        displace_structure(ts, mode, scale=factor),
        displace_structure(ts, mode, scale=-factor),
    )


def _classify(
    outcome: OptOutcome, reactant_ref: Structure, product_ref: Structure
) -> Optional[Tuple[float, float]]:
    """RMSD of a converged candidate to each reference, or None if the
    candidate didn't converge (nothing to classify)."""
    if not outcome.success or outcome.structure is None:
        return None
    return rmsd(outcome.structure, reactant_ref), rmsd(outcome.structure, product_ref)


def _assign_sides(
    a: OptOutcome, b: OptOutcome, reactant_ref: Structure, product_ref: Structure, log: List[str], tag_a: str, tag_b: str,
    reactant_trusted: bool = True, product_trusted: bool = True,
) -> Tuple[Optional[SideCandidate], Optional[SideCandidate]]:
    """
    Given two reoptimized candidates (from the +/- distortion, or the
    forward/reverse IRC endpoints), figure out which one is
    reactant-side and which is product-side by RMSD -- same
    either-pairing comparison irc.py uses to handle the fact that
    which physical direction is "+"/"forward" is arbitrary. Returns
    (reactant_candidate, product_candidate); either can be None if
    that candidate didn't converge, and a log line is added either way
    (including when the assignment is ambiguous).

    `reactant_trusted`/`product_trusted`: set False when that side's
    reference is only a stand-in (its own optimization failed -- see
    webapp.py/cli.py's endpoint-recovery path), not a real converged
    structure to compare against. With one side untrusted, the
    RMSD-to-that-side numbers are informative only -- classification
    goes by elimination against the side that *is* trusted instead of
    the usual paired-cost comparison.
    """
    a_rmsds = _classify(a, reactant_ref, product_ref)
    b_rmsds = _classify(b, reactant_ref, product_ref)

    if a_rmsds is None:
        log.append(f"{tag_a}: didn't reach a clean minimum -- not usable as a candidate")
    if b_rmsds is None:
        log.append(f"{tag_b}: didn't reach a clean minimum -- not usable as a candidate")
    if a_rmsds is None and b_rmsds is None:
        return None, None

    if a_rmsds is None or b_rmsds is None:
        outcome, (r, p), tag = (b, b_rmsds, tag_b) if a_rmsds is None else (a, a_rmsds, tag_a)
        if product_trusted and not reactant_trusted:
            side = "product" if p < r else "reactant"
        else:
            side = "reactant" if r < p else "product"
        log.append(f"{tag}: RMSD to reactant={r:.3f} A, to product={p:.3f} A -- classified as {side}")
        cand = SideCandidate(outcome, r, p)
        return (cand, None) if side == "reactant" else (None, cand)

    a_r, a_p = a_rmsds
    b_r, b_p = b_rmsds
    log.append(f"{tag_a}: RMSD to reactant={a_r:.3f} A, to product={a_p:.3f} A")
    log.append(f"{tag_b}: RMSD to reactant={b_r:.3f} A, to product={b_p:.3f} A")

    if not reactant_trusted or not product_trusted:
        # No real reference on one side -- classify by elimination
        # against the side that *is* trusted (default to "reactant" if
        # somehow neither is, so this always resolves to something).
        trusted_is_reactant = reactant_trusted or not product_trusted
        a_trusted_dist, b_trusted_dist = (a_r, b_r) if trusted_is_reactant else (a_p, b_p)
        close, far = (a, b) if a_trusted_dist < b_trusted_dist else (b, a)
        close_r, close_p = (a_r, a_p) if close is a else (b_r, b_p)
        far_r, far_p = (a_r, a_p) if far is a else (b_r, b_p)

        recovered_label = "product" if trusted_is_reactant else "reactant"
        log.append(
            f"Classifying by elimination against the trusted side (RMSDs to the untrusted "
            f"side above are informative only) -- the closer candidate is the trusted side, "
            f"the other is the recovered {recovered_label}, not independently confirmed."
        )
        close_cand = SideCandidate(close, close_r, close_p)
        far_cand = SideCandidate(far, far_r, far_p)
        return (close_cand, far_cand) if trusted_is_reactant else (far_cand, close_cand)

    # Two pairings: (a=reactant, b=product) or (a=product, b=reactant).
    # Pick whichever minimizes total mismatch -- same logic as irc.py's
    # clean_a/clean_b, just picking the better one instead of gating
    # pass/fail on a fixed threshold.
    cost_a_reactant = a_r + b_p
    cost_a_product = a_p + b_r

    if abs(cost_a_reactant - cost_a_product) < CLASSIFY_AMBIGUOUS_MARGIN:
        log.append(
            "Couldn't confidently tell which side is which (both pairings look "
            "similarly close) -- neither candidate offered for this run. Check "
            "the .log files by hand if you want to use one anyway."
        )
        return None, None

    if cost_a_reactant < cost_a_product:
        return SideCandidate(a, a_r, a_p), SideCandidate(b, b_r, b_p)
    return SideCandidate(b, b_r, b_p), SideCandidate(a, a_r, a_p)


def run_ts_distortion_check(
    ts_outcome: TSOutcome,
    reactant_ref: Structure,
    product_ref: Structure,
    out_dir: str,
    level: LevelChoice,
    nprocs: int = 2,
    mem_gb: int = 2,
    factor: float = 1.0,
    cancel_event: Optional[threading.Event] = None,
    reactant_trusted: bool = True,
    product_trusted: bool = True,
    executor: str = "local",
) -> DistortionOutcome:
    """
    Displace the converged TS along its own imaginary mode in both
    directions, reoptimize each (ground-state Opt Freq, through the
    same repair_and_optimize() the final-level reactant/product
    reopt already uses -- non-convergence and any new imaginary
    frequency get repaired the same way), and classify which side is
    reactant-derived and which is product-derived by RMSD.
    """
    log: List[str] = []

    if ts_outcome.structure is None or ts_outcome.imaginary_mode is None:
        log.append("No TS structure/imaginary mode available -- can't run the distortion check.")
        return DistortionOutcome(False, None, None, log)

    plus, minus = distort_along_imaginary_mode(ts_outcome.structure, ts_outcome.imaginary_mode, factor=factor)

    log.append(f"Distorting the TS along its imaginary mode by factor={factor:+.2f} in both directions:")
    plus_outcome = repair_and_optimize(
        plus, f"{out_dir}/ts_distort_plus.com",
        method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    log += [f"[distort +{factor:.2f}] {line}" for line in plus_outcome.log]

    minus_outcome = repair_and_optimize(
        minus, f"{out_dir}/ts_distort_minus.com",
        method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    log += [f"[distort -{factor:.2f}] {line}" for line in minus_outcome.log]

    reactant_candidate, product_candidate = _assign_sides(
        plus_outcome, minus_outcome, reactant_ref, product_ref, log,
        tag_a=f"distort+{factor:.2f}", tag_b=f"distort-{factor:.2f}",
        reactant_trusted=reactant_trusted, product_trusted=product_trusted,
    )
    success = reactant_candidate is not None or product_candidate is not None
    return DistortionOutcome(success, reactant_candidate, product_candidate, log)


def run_irc_endpoint_reopt(
    irc_outcome: IRCOutcome,
    reactant_ref: Structure,
    product_ref: Structure,
    out_dir: str,
    level: LevelChoice,
    nprocs: int = 2,
    mem_gb: int = 2,
    cancel_event: Optional[threading.Event] = None,
    reactant_trusted: bool = True,
    product_trusted: bool = True,
    executor: str = "local",
) -> DistortionOutcome:
    """
    Reoptimize the IRC's forward/reverse endpoints to proper minima --
    the raw endpoint is just the last point of a fixed-step path, not
    necessarily fully converged -- and classify each by RMSD, same as
    run_ts_distortion_check(). Any RMSDs irc_outcome already carries
    are from *before* this reoptimization, so fresh ones are computed
    here rather than reused.
    """
    log: List[str] = []

    if irc_outcome.forward_structure is None or irc_outcome.reverse_structure is None:
        log.append("IRC has no forward/reverse endpoint to reoptimize.")
        return DistortionOutcome(False, None, None, log)

    log.append("Reoptimizing the IRC forward/reverse endpoints to proper minima:")
    forward_outcome = repair_and_optimize(
        irc_outcome.forward_structure, f"{out_dir}/irc_forward_reopt.com",
        method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    log += [f"[IRC forward reopt] {line}" for line in forward_outcome.log]

    reverse_outcome = repair_and_optimize(
        irc_outcome.reverse_structure, f"{out_dir}/irc_reverse_reopt.com",
        method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    log += [f"[IRC reverse reopt] {line}" for line in reverse_outcome.log]

    reactant_candidate, product_candidate = _assign_sides(
        forward_outcome, reverse_outcome, reactant_ref, product_ref, log,
        tag_a="IRC forward reopt", tag_b="IRC reverse reopt",
        reactant_trusted=reactant_trusted, product_trusted=product_trusted,
    )
    success = reactant_candidate is not None or product_candidate is not None
    return DistortionOutcome(success, reactant_candidate, product_candidate, log)


def select_best_candidate(candidates: Dict[str, Optional[OptOutcome]]) -> Tuple[str, OptOutcome]:
    """
    Given e.g. {"guess": reactant_outcome, "ts_distortion": ...,
    "irc_reopt": ...} (any value may be None -- skipped), pick
    whichever converged candidate has the lowest (most stable) energy
    -- ZPE-corrected via the LogResult.energy property, the same
    number energetics.py already builds its report from. Returns
    (winning_label, winning_outcome). At least one candidate must be
    usable (non-None, converged, with an energy) or this raises.
    """
    usable = {
        label: outcome
        for label, outcome in candidates.items()
        if outcome is not None and outcome.success and outcome.result is not None and outcome.result.energy is not None
    }
    if not usable:
        raise ValueError(f"No usable candidate among {list(candidates)} -- every one is missing/didn't converge")

    winner_label = min(usable, key=lambda label: usable[label].result.energy)
    return winner_label, usable[winner_label]


@dataclass
class UnclassifiedRecovery:
    """Like DistortionOutcome, but for when there's no reactant/product
    reference to classify the two distorted-and-reoptimized sides
    against -- so they come back unclassified rather than guessed."""

    success: bool
    side_a: Optional[OptOutcome]
    side_b: Optional[OptOutcome]
    log: List[str] = field(default_factory=list)


def recover_endpoints_from_ts(
    ts_outcome: TSOutcome,
    out_dir: str,
    level: LevelChoice,
    nprocs: int = 2,
    mem_gb: int = 2,
    factor: float = 1.0,
    reactant_ref: Optional[Structure] = None,
    product_ref: Optional[Structure] = None,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> Union[DistortionOutcome, UnclassifiedRecovery]:
    """
    Regenerate the reactant/product straight from a found TS -- for the
    dedicated TS Search section, where (unlike a full mechanism study)
    there may be no independently-optimized reactant/product at all.

    If both `reactant_ref`/`product_ref` are given, this is exactly
    run_ts_distortion_check() (their RMSD to each distorted-and-
    reoptimized side is what classifies which is which). If neither is
    given -- a TS-guess-only search has no reference to compare
    against -- the two sides are returned unclassified: there's no way
    to tell which one is "reactant" vs. "product" from the TS alone,
    so this doesn't guess; it's on the user to tell them apart
    chemically.
    """
    if reactant_ref is not None and product_ref is not None:
        return run_ts_distortion_check(
            ts_outcome, reactant_ref, product_ref, out_dir, level,
            nprocs=nprocs, mem_gb=mem_gb, factor=factor, cancel_event=cancel_event, executor=executor,
        )

    log: List[str] = []
    if ts_outcome.structure is None or ts_outcome.imaginary_mode is None:
        log.append("No TS structure/imaginary mode available -- can't recover endpoints.")
        return UnclassifiedRecovery(False, None, None, log)

    plus, minus = distort_along_imaginary_mode(ts_outcome.structure, ts_outcome.imaginary_mode, factor=factor)
    log.append(
        f"No reactant/product reference given -- distorting the TS along its imaginary mode "
        f"by factor={factor:+.2f} in both directions and reoptimizing each side. Which side is "
        "actually the reactant vs. the product can't be determined from the TS alone -- that's "
        "on you to tell apart chemically."
    )

    side_a = repair_and_optimize(
        plus, f"{out_dir}/ts_recover_a.com",
        method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    log += [f"[side A] {line}" for line in side_a.log]

    side_b = repair_and_optimize(
        minus, f"{out_dir}/ts_recover_b.com",
        method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    log += [f"[side B] {line}" for line in side_b.log]

    success = side_a.success or side_b.success
    return UnclassifiedRecovery(
        success, side_a if side_a.success else None, side_b if side_b.success else None, log
    )


def format_candidate_comparison(side: str, candidates: Dict[str, Optional[OptOutcome]], winner_label: str) -> str:
    """Human-readable one-line trail entry, e.g.:
    'reactant candidates: guess=-152.401123, ts_distortion=-152.410877 (winner)'"""
    parts = []
    for label, outcome in candidates.items():
        if outcome is None or not outcome.success or outcome.result is None or outcome.result.energy is None:
            parts.append(f"{label}=n/a")
            continue
        marker = " (winner)" if label == winner_label else ""
        parts.append(f"{label}={outcome.result.energy:.6f}{marker}")
    return f"{side} candidates: " + ", ".join(parts)
