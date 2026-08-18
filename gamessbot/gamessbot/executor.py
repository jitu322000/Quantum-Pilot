"""
executor.py

One tiny dispatch point: run_gms() picks between running a .inp file
directly on this machine (local_runner.run_gamess -- the default) or
submitting it to a PBS queue (job_runner.run_pbs) via the same
"local"/"pbs" string every rhf.py/cis.py/casscf.py/xmcqdpt.py function
that runs a GAMESS job now takes as an `executor` parameter. Mirrors
gaussbot/executor.py's design. Swapping one for the other is a single
flag threaded down from webapp.py's/cli.py's session-wide PBS toggle,
not a code change anywhere a GAMESS job actually runs.
"""

from __future__ import annotations

import threading
from typing import Optional

from .local_runner import run_gamess

VALID_EXECUTORS = ("local", "pbs")


def run_gms(
    inp_path: str,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    executor: str = "local",
    cancel_event: Optional[threading.Event] = None,
) -> str:
    """Run `inp_path` with whichever executor is selected, both sharing
    the same contract: blocks until done, returns the log path, raises
    GamessRunError/GamessCancelledError on failure/cancellation.

    In "pbs" mode, `rungms_path`/`scratch_dir` (this machine's local
    rungms settings) are not used -- the PBS script has its own
    rungms-dev path (job_runner.DEFAULT_RUNGMS_DEV_PATH), edited via
    the session's PBS script template instead."""
    if executor == "local":
        return run_gamess(inp_path, rungms_path, scratch_dir, ncpus=ncpus, cancel_event=cancel_event)
    if executor == "pbs":
        from .job_runner import run_pbs

        return run_pbs(inp_path, ncpus=ncpus, cancel_event=cancel_event)
    raise ValueError(f"Unknown executor {executor!r} -- expected one of {VALID_EXECUTORS}")
