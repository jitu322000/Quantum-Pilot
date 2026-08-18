import os

import pytest

from gamessbot.parser import parse_casscf_log, parse_cis_log, parse_rhf_log, parse_transitn_log

# Real GAMESS logs already on this machine -- ground truth for the
# exact fixed-format lines these regexes match.
EXAM01_LOG = "/home/srr/jitendra/gamess/exam01.log"
IRON_CIS_LOG = "/home/srr/jitendra/iron/neutral/iron_complex_cis.log"
IRON_CAS_LOG = "/home/srr/jitendra/iron/neutral/iron_complex_cas.log"

pytestmark = pytest.mark.skipif(
    not (os.path.exists(EXAM01_LOG) and os.path.exists(IRON_CIS_LOG) and os.path.exists(IRON_CAS_LOG)),
    reason="real reference logs not present on this machine",
)


def test_parse_rhf_log_real_example():
    result = parse_rhf_log(EXAM01_LOG)
    assert result.normal_termination is True
    assert result.norb == 7
    assert result.n_occ == 4
    assert abs(result.energy - (-37.2322678015)) < 1e-6


ETHENE_CCD_RHF_LOG = "/home/srr/BOT/gamessbot/jobs/ethene-test-1/rhf_soscf.log"


@pytest.mark.skipif(not os.path.exists(ETHENE_CCD_RHF_LOG), reason="real reference log not present on this machine")
def test_parse_rhf_log_norb_uses_variation_space_not_cartesian_count():
    # Real ethene/cc-pVDZ (GBASIS=CCD) run with ISPHER=1: 40 cartesian
    # basis functions, but 2 spherical contaminants dropped, leaving 38
    # MOs in the variation space -- norb must be 38 (what the punched
    # $VEC block and any later $GUESS ... NORB= actually need), not the
    # cartesian 40, which caused a real rungms crash ("PREMATURE END OF
    # ORBITAL INPUT ... LOOKING FOR ORBITAL 39") before this was fixed.
    result = parse_rhf_log(ETHENE_CCD_RHF_LOG)
    assert result.normal_termination is True
    assert result.norb == 38


def test_parse_rhf_log_missing_termination(tmp_path):
    log = tmp_path / "job.log"
    log.write_text("some output\nstill running...\n")
    result = parse_rhf_log(str(log))
    assert result.normal_termination is False


def test_parse_cis_log_real_iron_complex():
    states = parse_cis_log(IRON_CIS_LOG)
    assert len(states) == 6

    s1 = states[0]
    assert s1.index == 1
    assert abs(s1.energy - (-985.5877297824)) < 1e-6
    assert s1.spin == 0.0
    assert s1.space_sym == "A"
    assert len(s1.transitions) == 5
    assert s1.transitions[0].from_mo == 72
    assert s1.transitions[0].to_mo == 78
    assert abs(s1.transitions[0].coefficient - (-0.07321866)) < 1e-8
    assert abs(s1.transitions[-1].coefficient - 0.95105506) < 1e-8

    # states are in ascending order, matching the log's own ordering
    assert [s.index for s in states] == [1, 2, 3, 4, 5, 6]


def test_parse_casscf_log_real_converged_run():
    result = parse_casscf_log(IRON_CAS_LOG)
    assert result.normal_termination is True
    assert result.converged is True
    assert abs(result.final_energy - (-985.6235289481)) < 1e-6
    assert result.iterations == 28
    assert len(result.state_energies) == 6
    assert result.state_energies[0] == (1, -985.711059110)
    assert result.state_energies[-1] == (6, -985.585559377)


def test_parse_casscf_log_not_converged(tmp_path):
    log = tmp_path / "job.log"
    log.write_text(
        "some MCSCF output\n"
        " EXCESSIVE NUMBER OF ITERATIONS...\n"
        " MCSCF IS NOT CONVERGED!\n"
        " FINAL MCSCF ENERGY IS        0.0000000000 AFTER 120 ITERATIONS\n"
        " STATE #    1  ENERGY =     -96.655474557\n"
        " EXECUTION OF GAMESS TERMINATED NORMALLY\n"
    )
    result = parse_casscf_log(str(log))
    assert result.normal_termination is True
    assert result.converged is False
    assert result.final_energy is None  # the 0.0 placeholder is discarded when not converged
    assert result.iterations == 120
    assert result.state_energies == [(1, -96.655474557)]


# ------------------------------------------------------------- parse_transitn_log

# Excerpted directly from a real production RUNTYP=TRANSITN (OPERAT=DM)
# log this user ran -- no real .log file exists on this machine for it
# (this is a brand-new feature this session), so the excerpt itself is
# persisted as the reference fixture rather than a hand-fabricated one.
TRANSITN_FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "transitn_oscillator_strengths.txt")


def test_parse_transitn_log_real_example():
    result = parse_transitn_log(TRANSITN_FIXTURE, nstate=3)
    # S0 (state 1) has no oscillator strength of its own -- only the
    # "1 vs itself" expectation-value block, not a transition moment.
    assert 1 not in result
    assert result[2] == 0.369622  # S1, relative to S0
    assert result[3] == 0.034274  # S2, relative to S0


def test_parse_transitn_log_stops_after_nstate_minus_one_entries(tmp_path):
    # nstate=2 should only look for the "1 2" pair (S1), not read ahead
    # into "1 3" (S2) even though it's present right after in the file.
    result = parse_transitn_log(TRANSITN_FIXTURE, nstate=2)
    assert result == {2: 0.369622}


def test_parse_transitn_log_no_oscillator_strengths_present(tmp_path):
    log = tmp_path / "job.log"
    log.write_text("EXECUTION OF GAMESS TERMINATED NORMALLY\n")
    assert parse_transitn_log(str(log), nstate=3) == {}
