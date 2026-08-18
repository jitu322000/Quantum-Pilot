"""
local_runner.py

Runs a Gaussian .com file directly on this machine with g09 (or g16)
-- no PBS involved. This is the "just run it here" path, for testing
the pipeline on a login/workstation node before job_runner.py wires
up queue submission.

    g09 < file.com > file.log
"""

from __future__ import annotations

import shutil
import subprocess
import threading
from pathlib import Path
from typing import Optional


class GaussianRunError(RuntimeError):
    """Raised when g09/g16 can't be found or exits with an error."""


class GaussianCancelledError(RuntimeError):
    """Raised when a run_local() call is cancelled via `cancel_event`
    before Gaussian finished on its own -- the log has whatever partial
    output Gaussian had written up to the point it was terminated."""


def find_gaussian_executable(preferred: str = "g09") -> str:
    """
    Locate the Gaussian executable on PATH. Raises with a clear
    message if it's missing -- most commonly because the site's
    Gaussian environment script hasn't been sourced in this shell.
    """
    exe = shutil.which(preferred)
    if exe is None:
        raise GaussianRunError(
            f"Could not find '{preferred}' on PATH. On most clusters you "
            "need to source Gaussian's setup script first, e.g.:\n"
            "  source /path/to/g09/bsd/g09.profile\n"
            "(the same line your PBS scripts already use)."
        )
    return exe


def run_local(
    com_path: str,
    log_path: Optional[str] = None,
    executable: str = "g09",
    timeout: Optional[int] = None,
    cancel_event: Optional[threading.Event] = None,
) -> str:
    """
    Run a .com file: `<executable> < com_path > log_path`. Blocks
    until Gaussian exits (or `timeout` seconds elapse).

    If `cancel_event` is given and gets set while Gaussian is running,
    the process is terminated (SIGTERM, then SIGKILL after a short
    grace period if it hasn't exited) and GaussianCancelledError is
    raised instead of returning -- `log_path` will have whatever
    partial output Gaussian had written up to that point, same as any
    other in-progress log.

    Returns the log path. Raises GaussianRunError if the executable
    is missing or exits with a non-zero return code -- a zero exit
    code does NOT by itself mean the calculation converged, see
    check_normal_termination().
    """
    com_path = str(com_path)
    log_path = log_path or str(Path(com_path).with_suffix(".log"))

    exe = find_gaussian_executable(executable)

    with open(com_path) as infile, open(log_path, "w") as outfile:
        proc = subprocess.Popen([exe], stdin=infile, stdout=outfile, stderr=subprocess.STDOUT)
        elapsed = 0.0
        while True:
            try:
                returncode = proc.wait(timeout=0.5)
                break
            except subprocess.TimeoutExpired:
                if cancel_event is not None and cancel_event.is_set():
                    proc.terminate()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    raise GaussianCancelledError(
                        f"Cancelled -- {exe} was terminated (see {log_path} for partial output)"
                    )
                elapsed += 0.5
                if timeout is not None and elapsed > timeout:
                    proc.kill()
                    proc.wait()
                    raise subprocess.TimeoutExpired(exe, timeout)

    if returncode != 0:
        raise GaussianRunError(f"{exe} exited with code {returncode} -- see {log_path}")

    return log_path


def check_normal_termination(log_path: str) -> bool:
    """
    The real ground-truth check: does the log end with Gaussian's own
    success marker? A clean process exit code isn't sufficient --
    Gaussian can exit 0 on some failure modes too.
    """
    with open(log_path) as f:
        lines = f.readlines()
    tail = lines[-5:]
    return any("Normal termination of Gaussian" in line for line in tail)
