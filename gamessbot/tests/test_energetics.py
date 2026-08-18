import os

import pytest

from gamessbot.active_space import ActiveSpaceSuggestion
from gamessbot.energetics import HARTREE_TO_EV, build_energy_table, format_latex_table
from gamessbot.parser import parse_casscf_log, parse_transitn_log, parse_xmcqdpt_log

IRON_MCQDPT_LOG = "/home/srr/jitendra/iron/neutral/uhf/iron_complex_mcqdpt2.log"

BUT_TEST_5_DIR = "/home/srr/BOT/gamessbot/jobs/but-test-5"
BUT_TEST_5_CASSCF_LOG = f"{BUT_TEST_5_DIR}/casscf_try1.log"
BUT_TEST_5_XMCQDPT_LOG = f"{BUT_TEST_5_DIR}/xmcqdpt.log"
BUT_TEST_5_TRANSITN_LOG = f"{BUT_TEST_5_DIR}/transitn.log"

pytestmark = pytest.mark.skipif(
    not os.path.exists(IRON_MCQDPT_LOG), reason="real reference log not present on this machine"
)


class _FakeCASSCFResult:
    def __init__(self, state_energies):
        self.state_energies = state_energies


def test_build_energy_table_real_iron_complex():
    xmcqdpt_result = parse_xmcqdpt_log(IRON_MCQDPT_LOG, nstate=6)
    casscf_result = _FakeCASSCFResult(xmcqdpt_result.casscf_state_energies)

    rows = build_energy_table(casscf_result, xmcqdpt_result)
    assert len(rows) == 6
    assert rows[0].label == "S0"
    assert rows[0].casscf_ev == 0.0
    assert rows[0].xmcqdpt_ev == 0.0

    expected_s1_casscf_ev = (
        xmcqdpt_result.casscf_state_energies[1][1] - xmcqdpt_result.casscf_state_energies[0][1]
    ) * HARTREE_TO_EV
    assert abs(rows[1].casscf_ev - expected_s1_casscf_ev) < 1e-9

    expected_s1_xmcqdpt_ev = (
        xmcqdpt_result.mcqdpt_state_energies[1][1] - xmcqdpt_result.mcqdpt_state_energies[0][1]
    ) * HARTREE_TO_EV
    assert abs(rows[1].xmcqdpt_ev - expected_s1_xmcqdpt_ev) < 1e-9

    # excitation energies should be non-negative for higher states (S0 is the minimum)
    assert all(row.casscf_ev >= 0 for row in rows)


def test_build_energy_table_casscf_only():
    casscf_result = _FakeCASSCFResult([(1, -76.0), (2, -75.8)])
    rows = build_energy_table(casscf_result, xmcqdpt_result=None)
    assert len(rows) == 2
    assert rows[0].xmcqdpt_ev is None
    assert rows[1].xmcqdpt_ev is None
    assert abs(rows[1].casscf_ev - (0.2 * HARTREE_TO_EV)) < 1e-9


def test_build_energy_table_empty():
    assert build_energy_table(_FakeCASSCFResult([])) == []


def test_format_latex_table_caption_and_structure():
    casscf_result = _FakeCASSCFResult([(1, -76.0), (2, -75.8), (3, -75.7)])
    rows = build_energy_table(casscf_result)
    active_space = ActiveSpaceSuggestion(
        nmcc=2, ndoc=4, nval=4, n_occ=6, norb=20,
        occ_selected=[3, 4, 5, 6], virt_selected=[7, 8, 9, 10], iorder=list(range(1, 21)),
    )
    tex = format_latex_table(rows, nstate=3, active_space=active_space)

    assert r"\begin{table}" in tex
    assert r"\end{table}" in tex
    assert "SA3 (8,8)" in tex  # 2*ndoc electrons, ndoc+nval orbitals
    assert "State & CASSCF (eV) & XMCQDPT (eV)" in tex
    assert r"S$_{0}$ & 0.00 & -" in tex  # no XMCQDPT result given -> single "-"
    assert tex.count(r"\\") == 4  # header row + 3 state rows


# ------------------------------------------------------- oscillator strength column

@pytest.mark.skipif(
    not (os.path.exists(BUT_TEST_5_CASSCF_LOG) and os.path.exists(BUT_TEST_5_XMCQDPT_LOG)
         and os.path.exists(BUT_TEST_5_TRANSITN_LOG)),
    reason="real reference logs not present on this machine",
)
def test_build_energy_table_and_latex_real_butadiene_with_oscillator_strengths():
    casscf_result = parse_casscf_log(BUT_TEST_5_CASSCF_LOG, nstate=3)
    xmcqdpt_result = parse_xmcqdpt_log(BUT_TEST_5_XMCQDPT_LOG, nstate=3)
    oscillator_strengths = parse_transitn_log(BUT_TEST_5_TRANSITN_LOG, nstate=3)
    assert oscillator_strengths == {2: 1.31516, 3: 0.0}

    rows = build_energy_table(casscf_result, xmcqdpt_result, oscillator_strengths)
    assert rows[0].oscillator_strength is None  # S0
    assert rows[1].oscillator_strength == 1.31516  # S1
    assert rows[2].oscillator_strength == 0.0  # S2

    active_space = ActiveSpaceSuggestion(
        nmcc=13, ndoc=2, nval=2, n_occ=15, norb=86,
        occ_selected=[14, 15], virt_selected=[16, 17], iorder=list(range(1, 87)),
    )
    tex = format_latex_table(rows, nstate=3, active_space=active_space)
    assert "State & CASSCF (eV) & XMCQDPT (eV) & $f$" in tex
    assert r"S$_{0}$ & 0.00 & 0.00 & -" in tex
    assert r"S$_{1}$" in tex and "1.3152" in tex
    assert r"S$_{2}$" in tex and "0.0000" in tex


def test_build_energy_table_omits_oscillator_strength_when_not_given():
    casscf_result = _FakeCASSCFResult([(1, -76.0), (2, -75.8)])
    rows = build_energy_table(casscf_result)
    assert all(row.oscillator_strength is None for row in rows)


def test_format_latex_table_omits_f_column_when_no_oscillator_strengths():
    casscf_result = _FakeCASSCFResult([(1, -76.0), (2, -75.8)])
    rows = build_energy_table(casscf_result)
    active_space = ActiveSpaceSuggestion(
        nmcc=2, ndoc=4, nval=4, n_occ=6, norb=20,
        occ_selected=[3, 4, 5, 6], virt_selected=[7, 8, 9, 10], iorder=list(range(1, 21)),
    )
    tex = format_latex_table(rows, nstate=2, active_space=active_space)
    assert "$f$" not in tex
    assert r"\begin{tabular}{lcc}" in tex


def test_format_latex_table_includes_f_column_when_oscillator_strengths_given():
    casscf_result = _FakeCASSCFResult([(1, -76.0), (2, -75.8), (3, -75.7)])
    rows = build_energy_table(casscf_result, oscillator_strengths={2: 0.5})
    active_space = ActiveSpaceSuggestion(
        nmcc=2, ndoc=4, nval=4, n_occ=6, norb=20,
        occ_selected=[3, 4, 5, 6], virt_selected=[7, 8, 9, 10], iorder=list(range(1, 21)),
    )
    tex = format_latex_table(rows, nstate=3, active_space=active_space)
    assert r"\begin{tabular}{lccc}" in tex
    assert r"S$_{0}$ & 0.00 & - & -" in tex
    assert "0.5000" in tex
    assert r"S$_{2}$" in tex and "- \\\\" in tex  # S2 has no oscillator strength given
