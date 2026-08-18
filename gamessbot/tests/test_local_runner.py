import os
import stat
from pathlib import Path

import pytest

from gamessbot.local_runner import GamessRunError, check_normal_termination, run_gamess

# A tiny stand-in for rungms: writes a fake .log to stdout (redirected by
# run_gamess to <stem>.log, matching the real script's own `> $input.log`)
# and drops a .dat file in the scratch dir, exercising the subprocess/
# move-the-.dat plumbing without a real GAMESS install.
FAKE_RUNGMS_OK = """#!/bin/sh
# $1=stem $2=verno $3=ncpus $4=ncpus
echo " EXECUTION OF GAMESS TERMINATED NORMALLY Mon Aug 10 00:00:00 2026"
echo "fake punch data" > "$SCRATCH_DIR/$1.dat"
"""

FAKE_RUNGMS_FAIL = """#!/bin/sh
echo " EXECUTION OF GAMESS TERMINATED -ABNORMALLY- Mon Aug 10 00:00:00 2026"
exit 1
"""


def _install_fake_rungms(tmp_path, script_text):
    script = tmp_path / "rungms"
    script.write_text(script_text)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def test_run_gamess_success_moves_dat(tmp_path, monkeypatch):
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    monkeypatch.setenv("SCRATCH_DIR", str(scratch_dir))
    rungms = _install_fake_rungms(tmp_path, FAKE_RUNGMS_OK)

    job_dir = tmp_path / "job"
    job_dir.mkdir()
    inp = job_dir / "water.inp"
    inp.write_text(" $CONTRL SCFTYP=RHF $END\n")

    log_path = run_gamess(str(inp), str(rungms), str(scratch_dir), ncpus=1)

    assert Path(log_path) == job_dir / "water.log"
    assert check_normal_termination(log_path) is True
    assert (job_dir / "water.dat").exists()  # moved out of scratch
    assert not (scratch_dir / "water.dat").exists()


def test_run_gamess_clears_stale_scratch_files(tmp_path, monkeypatch):
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    monkeypatch.setenv("SCRATCH_DIR", str(scratch_dir))
    (scratch_dir / "water.F10").write_text("stale")
    rungms = _install_fake_rungms(tmp_path, FAKE_RUNGMS_OK)

    job_dir = tmp_path / "job"
    job_dir.mkdir()
    inp = job_dir / "water.inp"
    inp.write_text(" $CONTRL SCFTYP=RHF $END\n")

    run_gamess(str(inp), str(rungms), str(scratch_dir), ncpus=1)

    assert not (scratch_dir / "water.F10").exists()


def test_run_gamess_missing_rungms_raises(tmp_path):
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    inp = job_dir / "water.inp"
    inp.write_text(" $CONTRL SCFTYP=RHF $END\n")

    with pytest.raises(GamessRunError):
        run_gamess(str(inp), str(tmp_path / "no-such-rungms"), str(tmp_path))


def test_run_gamess_nonzero_exit_raises(tmp_path, monkeypatch):
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    rungms = _install_fake_rungms(tmp_path, FAKE_RUNGMS_FAIL)

    job_dir = tmp_path / "job"
    job_dir.mkdir()
    inp = job_dir / "water.inp"
    inp.write_text(" $CONTRL SCFTYP=RHF $END\n")

    with pytest.raises(GamessRunError):
        run_gamess(str(inp), str(rungms), str(scratch_dir), ncpus=1)


def test_check_normal_termination_false_on_incomplete_log(tmp_path):
    log = tmp_path / "job.log"
    log.write_text("some output\nstill running...\n")
    assert check_normal_termination(str(log)) is False
