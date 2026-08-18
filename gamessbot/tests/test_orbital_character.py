import os

import pytest

from gamessbot.orbital_character import parse_orbital_character

# Real water STO-3G RHF log from an actual end-to-end run -- ground truth
# for the EIGENVECTORS table format, same file used to ad-hoc verify
# parse_orbital_character() during development.
WATER_RHF_LOG = "/home/srr/BOT/gamessbot/jobs/water_casscf_gui/rhf_soscf.log"

pytestmark = pytest.mark.skipif(
    not os.path.exists(WATER_RHF_LOG), reason="real reference log not present on this machine"
)


def test_parse_orbital_character_real_water_log():
    character = parse_orbital_character(WATER_RHF_LOG, [1, 2, 3, 4, 5, 6, 7])

    assert set(character.keys()) == {1, 2, 3, 4, 5, 6, 7}
    for mo, c in character.items():
        assert c.atom_symbol == "O"
        assert c.atom_number == 1
        assert 0.0 < abs(c.coefficient) <= 1.0001

    # MO 5 is the classic O 2pz lone-pair HOMO for water/STO-3G.
    assert character[5].shell == "Z"
    assert abs(character[5].coefficient - 1.0) < 1e-3
    assert character[5].label == "O1 Z"


def test_parse_orbital_character_missing_mo_returns_partial():
    character = parse_orbital_character(WATER_RHF_LOG, [5, 999])
    assert 5 in character
    assert 999 not in character


def test_parse_orbital_character_no_eigenvectors_section(tmp_path):
    log = tmp_path / "job.log"
    log.write_text("some output\nno eigenvectors here\n")
    assert parse_orbital_character(str(log), [1, 2, 3]) == {}
