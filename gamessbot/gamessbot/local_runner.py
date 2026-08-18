"""
local_runner.py

Runs a GAMESS .inp file directly on this machine via `rungms` -- no PBS
involved (a PBS equivalent, following the same pattern gaussbot uses for
Gaussian, is future work once this round is validated).

    rungms <stem> 00 <ncpus> <ncpus> > <stem>.log

reproducing the exact shape of your own script: clear any stale scratch
files for this job name first, run, then move the resulting .dat (PUNCH
file, containing the converged orbitals) back next to the .inp/.log.
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import threading
from pathlib import Path
from typing import Optional


class GamessRunError(RuntimeError):
    """Raised when rungms can't be found or exits with a non-zero code.
    A zero exit code does NOT by itself mean the calculation converged --
    see check_normal_termination()."""


class GamessCancelledError(RuntimeError):
    """Raised when a run_gamess() call is cancelled via `cancel_event`
    before GAMESS finished on its own -- the log has whatever partial
    output GAMESS had written up to the point it was terminated."""


def run_gamess(
    inp_path: str,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    log_path: Optional[str] = None,
    timeout: Optional[int] = None,
    cancel_event: Optional[threading.Event] = None,
) -> str:
    """
    Run a .inp file: clears any stale `<scratch_dir>/<stem>.*` files for
    this job name, runs `<rungms_path> <stem> 00 <ncpus> <ncpus>` with
    stdout redirected to the log, then moves the resulting
    `<scratch_dir>/<stem>.dat` (the punched orbitals) back next to the
    .inp/.log. Blocks until GAMESS exits (or `timeout` seconds elapse).

    Returns the log path.
    """
    inp_path = str(inp_path)
    inp_dir = os.path.dirname(os.path.abspath(inp_path)) or "."
    stem = Path(inp_path).stem
    log_path = log_path or str(Path(inp_path).with_suffix(".log"))

    if not os.path.exists(rungms_path):
        raise GamessRunError(
            f"Could not find rungms at {rungms_path!r}. Check the GAMESS settings "
            "(rungms path) -- this needs to point at your GAMESS installation's rungms script."
        )

    for stale in glob.glob(os.path.join(scratch_dir, f"{stem}.*")):
        os.remove(stale)

    with open(log_path, "w") as outfile:
        proc = subprocess.Popen(
            [rungms_path, stem, "00", str(ncpus), str(ncpus)],
            cwd=inp_dir, stdout=outfile, stderr=subprocess.STDOUT,
        )
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
                    raise GamessCancelledError(
                        f"Cancelled -- rungms was terminated (see {log_path} for partial output)"
                    )
                elapsed += 0.5
                if timeout is not None and elapsed > timeout:
                    proc.kill()
                    proc.wait()
                    raise subprocess.TimeoutExpired(rungms_path, timeout)

    if returncode != 0:
        raise GamessRunError(f"rungms exited with code {returncode} -- see {log_path}")

    dat_src = os.path.join(scratch_dir, f"{stem}.dat")
    if os.path.exists(dat_src):
        shutil.move(dat_src, os.path.join(inp_dir, f"{stem}.dat"))

    return log_path


def check_normal_termination(log_path: str) -> bool:
    """
    Confirmed exact string against a real log: "EXECUTION OF GAMESS
    TERMINATED NORMALLY <timestamp>". Searched over the whole file rather
    than just the tail, since rungms wrappers can append extra
    housekeeping output after GAMESS's own last printed line.
    """
    with open(log_path) as f:
        text = f.read()
    return "TERMINATED NORMALLY" in text
