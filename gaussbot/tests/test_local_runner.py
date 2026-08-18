import os
import stat
from pathlib import Path

import pytest

from gaussbot.local_runner import (
    GaussianRunError,
    check_normal_termination,
    find_gaussian_executable,
    run_local,
)

# A tiny stand-in for g09: reads the .com from stdin (so run_local's
# plumbing is exercised for real) and writes a fake log tail. This
# lets us test the subprocess wiring without a real Gaussian install.
FAKE_G09_OK = """#!/bin/sh
cat
echo ""
echo " Job cpu time: 0 days 0 hours 0 minutes 1.0 seconds."
echo " Normal termination of Gaussian 09 at Mon Aug 10 00:00:00 2026."
"""

FAKE_G09_FAIL = """#!/bin/sh
cat
echo "Error termination via Lnk1e"
exit 1
"""

SAMPLE_COM = "%chk=job.chk\n#p PM6 Opt\n\ntitle\n\n0 1\nH 0.0 0.0 0.0\n\n"


def _install_fake_g09(tmp_path, monkeypatch, script_text):
    script = tmp_path / "g09"
    script.write_text(script_text)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    return script


def test_find_gaussian_executable_missing_raises(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))  # empty dir, nothing on PATH
    with pytest.raises(GaussianRunError):
        find_gaussian_executable("g09")


def test_find_gaussian_executable_found(tmp_path, monkeypatch):
    script = _install_fake_g09(tmp_path, monkeypatch, FAKE_G09_OK)
    assert find_gaussian_executable("g09") == str(script)


def test_run_local_success_and_normal_termination(tmp_path, monkeypatch):
    _install_fake_g09(tmp_path, monkeypatch, FAKE_G09_OK)
    com = tmp_path / "job.com"
    com.write_text(SAMPLE_COM)

    log_path = run_local(str(com))

    assert Path(log_path).exists()
    assert Path(log_path).read_text().startswith("%chk=job.chk")  # stdin was piped through
    assert check_normal_termination(log_path) is True


def test_run_local_nonzero_exit_raises(tmp_path, monkeypatch):
    _install_fake_g09(tmp_path, monkeypatch, FAKE_G09_FAIL)
    com = tmp_path / "job.com"
    com.write_text(SAMPLE_COM)

    with pytest.raises(GaussianRunError):
        run_local(str(com))


def test_check_normal_termination_false_on_incomplete_log(tmp_path):
    log = tmp_path / "job.log"
    log.write_text("some output\nstill running...\n")
    assert check_normal_termination(str(log)) is False


def test_run_local_default_log_path_next_to_com(tmp_path, monkeypatch):
    _install_fake_g09(tmp_path, monkeypatch, FAKE_G09_OK)
    com = tmp_path / "myjob.com"
    com.write_text(SAMPLE_COM)

    log_path = run_local(str(com))

    assert log_path == str(tmp_path / "myjob.log")
