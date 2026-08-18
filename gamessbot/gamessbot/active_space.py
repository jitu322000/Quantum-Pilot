"""
active_space.py

Suggests a CASSCF active space (NMCC/NDOC/NVAL) from a CIS calculation's
excited-state transitions -- the workflow you described: scan CIS states
for SAP coefficients with |c| > 0.2, the occupied MOs those transitions
originate from become NDOC (active doubly-occupied), the virtual MOs
they land on become NVAL (active virtual), and everything else occupied
stays frozen as NMCC. Confirmed against theochem.cc's own CASSCF
tutorial (same |c| > 0.20 threshold) and your real iron-complex CASSCF
input ($DRT NMCC=67 NDOC=7 NVAL=2), including its use of $GUESS
NORDER/IORDER to bring the selected orbitals -- which aren't always
contiguous with HOMO/LUMO -- into the single contiguous window $DRT's
NMCC/NDOC/NVAL scheme requires.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

from .parser import ExcitedState

DEFAULT_THRESHOLD = 0.20
DEFAULT_MAX_ELECTRONS = 16
DEFAULT_MAX_ORBITALS = 16


@dataclass
class ActiveSpaceSuggestion:
    """What suggest_active_space() came back with -- ready to hand to
    gamess_input.build_casscf_input() as-is, or to adjust first (see
    cli.py's/webapp.py's confirm-or-edit step, per your request that the
    system propose an active space and ask before running with it)."""

    nmcc: int
    ndoc: int
    nval: int
    n_occ: int
    norb: int
    occ_selected: List[int]  # active docc MOs, ascending
    virt_selected: List[int]  # active virtual MOs, ascending
    occ_dropped: List[int] = field(default_factory=list)  # candidates cut by the electron/orbital cap
    virt_dropped: List[int] = field(default_factory=list)
    scores: Dict[int, float] = field(default_factory=dict)  # MO -> best |SAP coefficient| that nominated it
    iorder: List[int] = field(default_factory=list)  # full NORB-length permutation for $GUESS NORDER/IORDER

    @property
    def capped(self) -> bool:
        return bool(self.occ_dropped or self.virt_dropped)


def _build_iorder(occ_selected: List[int], virt_selected: List[int], n_occ: int, norb: int) -> List[int]:
    """
    Builds the NORB-length permutation for $GUESS NORDER=1/IORDER using
    only symmetric pairwise swaps (IORDER(a)=b together with
    IORDER(b)=a) -- simpler and more conservative than a longer
    reordering cycle, per your preference. IORDER(i) = the original MO
    that becomes new orbital i; positions left at IORDER(i)=i (identity)
    are omitted entirely from the actual input (that's GAMESS's own
    default), so this only ever touches the MOs that actually move.

    For each active MO not already sitting in its target window (NDOC
    slots right after NMCC, then NVAL slots right after that), it's
    swapped directly with whichever MO currently occupies that window
    slot. This is always exactly resolvable with disjoint pairwise
    swaps here (never a longer cycle) because occ_selected is always
    <= n_occ and virt_selected is always > n_occ (the CIS occupied/
    virtual split), which matches the docc/val sub-windows' own
    boundary exactly -- so the two swap passes below can never collide
    or displace each other's placements.
    """
    iorder = list(range(1, norb + 1))
    nmcc = n_occ - len(occ_selected)

    def _swap_into_window(selected: List[int], window_start: int) -> None:
        window = set(range(window_start, window_start + len(selected)))
        already_placed = set(selected) & window
        outside = sorted(set(selected) - already_placed)
        vacated = sorted(window - already_placed)
        for mo, position in zip(outside, vacated):
            iorder[position - 1] = mo
            iorder[mo - 1] = position

    _swap_into_window(occ_selected, nmcc + 1)
    _swap_into_window(virt_selected, nmcc + len(occ_selected) + 1)
    return iorder


def suggest_active_space(
    states: List[ExcitedState],
    n_occ: int,
    norb: int,
    threshold: float = DEFAULT_THRESHOLD,
    max_electrons: int = DEFAULT_MAX_ELECTRONS,
    max_orbitals: int = DEFAULT_MAX_ORBITALS,
) -> ActiveSpaceSuggestion:
    """
    Scans every transition in `states` (typically cis.CISOutcome.states)
    for |coefficient| > threshold. The occupied MOs those transitions
    originate from become active-docc candidates; the virtual MOs they
    land on become active-virtual candidates. Ranked by each MO's best
    (largest-magnitude) contributing coefficient -- not by distance from
    HOMO/LUMO, since the important transitions aren't always near the
    frontier orbitals (per your own iron-complex experience, where this
    happens and $GUESS's NORDER/IORDER is what makes it workable).

    Capped to `max_electrons` (NDOC <= max_electrons//2, since each
    active-docc orbital contributes 2 electrons in the reference
    configuration) and `max_orbitals` (NDOC+NVAL <= max_orbitals) total,
    keeping the highest-scoring candidates when a cap trims the list --
    default 16/16, already a large CASSCF.

    Raises ValueError if no transition anywhere exceeds `threshold` on
    both the occupied and virtual side -- nothing to build an active
    space from; lower the threshold, look at more CIS states, or pick
    the active space by hand.
    """
    occ_scores: Dict[int, float] = {}
    virt_scores: Dict[int, float] = {}
    for state in states:
        for t in state.transitions:
            c = abs(t.coefficient)
            if c > threshold:
                occ_scores[t.from_mo] = max(occ_scores.get(t.from_mo, 0.0), c)
                virt_scores[t.to_mo] = max(virt_scores.get(t.to_mo, 0.0), c)

    if not occ_scores or not virt_scores:
        raise ValueError(
            f"No CIS transition has |coefficient| > {threshold} on both the occupied and "
            "virtual side -- nothing to build an active space from. Try a lower threshold, "
            "more CIS states, or pick the active space by hand."
        )

    occ_ranked = sorted(occ_scores, key=lambda mo: -occ_scores[mo])
    virt_ranked = sorted(virt_scores, key=lambda mo: -virt_scores[mo])

    max_ndoc = max_electrons // 2
    occ_kept, occ_dropped = occ_ranked[:max_ndoc], occ_ranked[max_ndoc:]

    remaining_orbitals = max(0, max_orbitals - len(occ_kept))
    virt_kept, virt_dropped = virt_ranked[:remaining_orbitals], virt_ranked[remaining_orbitals:]

    occ_selected = sorted(occ_kept)
    virt_selected = sorted(virt_kept)

    return ActiveSpaceSuggestion(
        nmcc=n_occ - len(occ_selected), ndoc=len(occ_selected), nval=len(virt_selected),
        n_occ=n_occ, norb=norb, occ_selected=occ_selected, virt_selected=virt_selected,
        occ_dropped=sorted(occ_dropped), virt_dropped=sorted(virt_dropped),
        scores={**occ_scores, **virt_scores},
        iorder=_build_iorder(occ_selected, virt_selected, n_occ, norb),
    )


def build_active_space(n_occ: int, norb: int, occ_selected: List[int], virt_selected: List[int]) -> ActiveSpaceSuggestion:
    """
    Builds an ActiveSpaceSuggestion directly from an explicit MO
    selection -- no CIS-driven scoring -- for when the user overrides
    the automatic suggestion with their own choice of active orbitals
    (per your own description: sometimes the important orbitals need to
    be picked by hand, by looking at the orbitals' eigenvectors, not
    just the CIS coefficients).

    occ_selected must all be occupied MOs (<= n_occ) and virt_selected
    must all be virtual MOs (> n_occ) -- this is what keeps
    _build_iorder()'s pairwise swaps collision-free; a MO on the wrong
    side (e.g. a virtual index typed into the occupied field) is a
    real, easy-to-make input mistake, so this is checked rather than
    silently building a wrong active space.
    """
    occ_selected = sorted(occ_selected)
    virt_selected = sorted(virt_selected)
    if any(mo > n_occ for mo in occ_selected):
        raise ValueError(f"Active occupied MOs must all be <= {n_occ} (n_occ): {occ_selected}")
    if any(mo <= n_occ for mo in virt_selected):
        raise ValueError(f"Active virtual MOs must all be > {n_occ} (n_occ): {virt_selected}")
    return ActiveSpaceSuggestion(
        nmcc=n_occ - len(occ_selected), ndoc=len(occ_selected), nval=len(virt_selected),
        n_occ=n_occ, norb=norb, occ_selected=occ_selected, virt_selected=virt_selected,
        iorder=_build_iorder(occ_selected, virt_selected, n_occ, norb),
    )


def shrink_active_space(
    suggestion: ActiveSpaceSuggestion, max_electrons: int, max_orbitals: int
) -> ActiveSpaceSuggestion:
    """
    Re-caps an already-built active space down to a smaller one, using
    the same scores it was originally ranked by -- for the CASSCF
    convergence-recovery path (casscf.run_casscf_staged()): when even a
    higher-MAXIT retry doesn't converge, try a smaller/easier active
    space first, then use ITS converged orbitals as a much better
    starting guess to reach the original, larger active space
    (regrow_active_space()).
    """
    max_ndoc = max_electrons // 2
    occ_ranked = sorted(suggestion.occ_selected, key=lambda mo: -suggestion.scores.get(mo, 0.0))
    occ_kept = sorted(occ_ranked[:max_ndoc])

    remaining_orbitals = max(0, max_orbitals - len(occ_kept))
    virt_ranked = sorted(suggestion.virt_selected, key=lambda mo: -suggestion.scores.get(mo, 0.0))
    virt_kept = sorted(virt_ranked[:remaining_orbitals])

    return ActiveSpaceSuggestion(
        nmcc=suggestion.n_occ - len(occ_kept), ndoc=len(occ_kept), nval=len(virt_kept),
        n_occ=suggestion.n_occ, norb=suggestion.norb, occ_selected=occ_kept, virt_selected=virt_kept,
        scores=suggestion.scores,
        iorder=_build_iorder(occ_kept, virt_kept, suggestion.n_occ, suggestion.norb),
    )


def regrow_active_space(
    original: ActiveSpaceSuggestion, converged_smaller: ActiveSpaceSuggestion
) -> ActiveSpaceSuggestion:
    """
    Builds the IORDER needed to restart `original`'s (larger) active
    space using orbitals already optimized -- and thus already
    reordered once -- by a smaller, easier CASSCF pass run first as a
    better starting guess (per your convergence-recovery request:
    smaller active space -> converge -> grow back to the originally
    requested one using those orbitals).

    The smaller pass's punched $VEC is what this run's own $GUESS reads,
    so this new IORDER is expressed relative to THAT file's own column
    order (starting from identity again, not from converged_smaller's
    permutation) -- it composes converged_smaller's reordering with
    original's implicitly, by first finding where converged_smaller's
    own IORDER actually put each of original's target MOs, then doing
    the same symmetric-pairwise-swap _build_iorder() does, using those
    current column positions as the pool instead of the raw original
    RHF MO indices.
    """
    # inverse of converged_smaller.iorder: original RHF MO index -> its current column in that punch file
    current_position = {mo: pos + 1 for pos, mo in enumerate(converged_smaller.iorder)}

    iorder = list(range(1, original.norb + 1))
    nmcc = original.n_occ - len(original.occ_selected)

    def _swap_into_window(selected: List[int], window_start: int) -> None:
        window = set(range(window_start, window_start + len(selected)))
        current_cols = [current_position[mo] for mo in selected]
        already_placed = set(current_cols) & window
        outside = sorted(set(current_cols) - already_placed)
        vacated = sorted(window - already_placed)
        for col, position in zip(outside, vacated):
            iorder[position - 1] = col
            iorder[col - 1] = position

    _swap_into_window(original.occ_selected, nmcc + 1)
    _swap_into_window(original.virt_selected, nmcc + len(original.occ_selected) + 1)

    return ActiveSpaceSuggestion(
        nmcc=original.nmcc, ndoc=original.ndoc, nval=original.nval,
        n_occ=original.n_occ, norb=original.norb,
        occ_selected=original.occ_selected, virt_selected=original.virt_selected,
        scores=original.scores, iorder=iorder,
    )


def format_active_space_summary(suggestion: ActiveSpaceSuggestion) -> str:
    """Human-readable explanation of what was picked and why -- meant to
    be shown to the user before running CASSCF, so they can confirm or
    adjust NMCC/NDOC/NVAL by hand."""
    lines = [
        f"Suggested active space: NMCC={suggestion.nmcc}  NDOC={suggestion.ndoc}  NVAL={suggestion.nval}"
        f"  (CAS({2 * suggestion.ndoc} electrons, {suggestion.ndoc + suggestion.nval} orbitals))",
        f"  Active occupied MOs: {suggestion.occ_selected}",
        f"  Active virtual MOs:  {suggestion.virt_selected}",
    ]
    for mo in suggestion.occ_selected + suggestion.virt_selected:
        lines.append(f"    MO {mo}: best |SAP coefficient| = {suggestion.scores.get(mo, 0):.4f}")
    if suggestion.capped:
        lines.append(
            "  Active space was capped at the default max (16 electrons, 16 orbitals) -- "
            "dropped lower-scoring candidates:"
        )
        if suggestion.occ_dropped:
            lines.append(f"    occupied MOs not included: {suggestion.occ_dropped}")
        if suggestion.virt_dropped:
            lines.append(f"    virtual MOs not included:  {suggestion.virt_dropped}")
    return "\n".join(lines)
