import os

import pytest

from gamessbot.active_space import (
    ActiveSpaceSuggestion,
    build_active_space,
    regrow_active_space,
    shrink_active_space,
    suggest_active_space,
)
from gamessbot.parser import parse_cis_log

IRON_CIS_LOG = "/home/srr/jitendra/iron/neutral/iron_complex_cis.log"

pytestmark = pytest.mark.skipif(
    not os.path.exists(IRON_CIS_LOG), reason="real reference log not present on this machine"
)


def test_suggest_active_space_real_iron_complex():
    states = parse_cis_log(IRON_CIS_LOG)
    # n_occ=74 is this real (ECP) log's own working occupied-orbital count
    # (the ECP-adjusted "KEPT" value) -- confirmed against the log itself.
    suggestion = suggest_active_space(states, n_occ=74, norb=363)

    assert suggestion.capped  # more than 16/16 worth of candidates exist above threshold
    assert suggestion.ndoc <= 8  # capped at max_electrons(16)//2
    assert suggestion.ndoc + suggestion.nval <= 16  # capped at max_orbitals
    assert suggestion.nmcc == 74 - suggestion.ndoc

    # the permutation is a true bijection over 1..norb
    assert len(suggestion.iorder) == 363
    assert sorted(suggestion.iorder) == list(range(1, 364))

    # the active window is contiguous and starts right after NMCC
    active = suggestion.iorder[suggestion.nmcc:suggestion.nmcc + suggestion.ndoc + suggestion.nval]
    assert set(active) == set(suggestion.occ_selected) | set(suggestion.virt_selected)

    # every reordering is a symmetric pairwise swap (an involution): if
    # position a maps to value b, position b must map back to a -- per
    # your preference for simple swaps over a longer reordering cycle
    for position, value in enumerate(suggestion.iorder, start=1):
        assert suggestion.iorder[value - 1] == position


def test_suggest_active_space_no_transitions_above_threshold():
    with pytest.raises(ValueError):
        suggest_active_space([], n_occ=5, norb=7)


def test_suggest_active_space_respects_smaller_cap():
    states = parse_cis_log(IRON_CIS_LOG)
    suggestion = suggest_active_space(states, n_occ=74, norb=363, max_electrons=4, max_orbitals=4)
    assert suggestion.ndoc == 2  # 4 electrons // 2
    assert suggestion.ndoc + suggestion.nval == 4
    assert suggestion.capped


def test_build_active_space_manual_selection():
    active_space = build_active_space(n_occ=5, norb=7, occ_selected=[4, 5], virt_selected=[6, 7])
    assert active_space.nmcc == 3
    assert active_space.ndoc == 2
    assert active_space.nval == 2
    assert active_space.iorder == [1, 2, 3, 4, 5, 6, 7]  # already contiguous, no reordering needed


def test_build_active_space_noncontiguous_reordering():
    # active MOs scattered, not adjacent to each other or to the frozen/virtual boundary
    active_space = build_active_space(n_occ=10, norb=20, occ_selected=[2, 7], virt_selected=[15, 18])
    assert active_space.nmcc == 8
    assert active_space.ndoc == 2
    assert active_space.nval == 2
    # frozen (8) + active occ (2) + active virt (2) + remaining virtuals (8) = 20
    assert len(active_space.iorder) == 20
    assert sorted(active_space.iorder) == list(range(1, 21))
    # active window starts right after NMCC=8, contiguous
    assert active_space.iorder[8:12] == [2, 7, 15, 18]

    # a symmetric pairwise swap: MO 2 <-> position 9, MO 7 <-> position 10, etc.
    assert active_space.iorder[1] == 9  # position 2 now holds what was at 9
    assert active_space.iorder[8] == 2  # and position 9 holds MO 2, symmetrically


def test_build_active_space_rejects_wrong_side_mo():
    with pytest.raises(ValueError):
        build_active_space(n_occ=5, norb=7, occ_selected=[4, 6], virt_selected=[7])  # 6 is virtual, not occupied
    with pytest.raises(ValueError):
        build_active_space(n_occ=5, norb=7, occ_selected=[4], virt_selected=[3, 7])  # 3 is occupied, not virtual


def _make_suggestion(n_occ, norb, occ_selected, virt_selected, scores):
    from gamessbot.active_space import _build_iorder  # noqa: PLC0415 -- test-only, internal helper reuse

    return ActiveSpaceSuggestion(
        nmcc=n_occ - len(occ_selected), ndoc=len(occ_selected), nval=len(virt_selected),
        n_occ=n_occ, norb=norb, occ_selected=sorted(occ_selected), virt_selected=sorted(virt_selected),
        scores=scores, iorder=_build_iorder(sorted(occ_selected), sorted(virt_selected), n_occ, norb),
    )


def test_shrink_active_space_keeps_highest_scoring():
    original = _make_suggestion(
        n_occ=10, norb=20, occ_selected=[3, 4, 7, 8], virt_selected=[12, 13, 15, 16],
        scores={3: 0.9, 4: 0.5, 7: 0.3, 8: 0.2, 12: 0.8, 13: 0.6, 15: 0.4, 16: 0.25},
    )
    smaller = shrink_active_space(original, max_electrons=4, max_orbitals=4)
    assert smaller.ndoc == 2
    assert smaller.nval == 2
    assert smaller.occ_selected == [3, 4]  # highest-scored occ MOs kept
    assert smaller.virt_selected == [12, 13]  # highest-scored virt MOs kept
    assert sorted(smaller.iorder) == list(range(1, 21))


def test_regrow_active_space_restores_original_selection():
    original = _make_suggestion(
        n_occ=10, norb=20, occ_selected=[3, 4, 7, 8], virt_selected=[12, 13, 15, 16],
        scores={3: 0.9, 4: 0.5, 7: 0.3, 8: 0.2, 12: 0.8, 13: 0.6, 15: 0.4, 16: 0.25},
    )
    smaller = shrink_active_space(original, max_electrons=4, max_orbitals=4)
    regrown = regrow_active_space(original, smaller)

    # regrown keeps original's own NMCC/NDOC/NVAL and MO selection
    assert (regrown.nmcc, regrown.ndoc, regrown.nval) == (original.nmcc, original.ndoc, original.nval)
    assert regrown.occ_selected == original.occ_selected
    assert regrown.virt_selected == original.virt_selected

    # it's a true bijection over 1..norb
    assert sorted(regrown.iorder) == list(range(1, 21))

    # tracing through BOTH reorderings (smaller's own, then regrown's,
    # applied to the smaller pass's punched $VEC), the physical MOs
    # that land in original's active window are exactly original's
    # own full (larger) selection -- i.e. this correctly "grows" the
    # active space back up using the smaller pass's optimized orbitals
    window = range(regrown.nmcc + 1, regrown.nmcc + regrown.ndoc + regrown.nval + 1)
    physical_mos_in_window = {smaller.iorder[regrown.iorder[p - 1] - 1] for p in window}
    assert physical_mos_in_window == set(original.occ_selected) | set(original.virt_selected)
