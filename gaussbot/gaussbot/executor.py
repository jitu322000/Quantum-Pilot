"""
executor.py

One tiny dispatch point: run_com() picks between running a .com file
directly on this machine (local_runner.run_local -- the default) or
submitting it to a PBS queue (job_runner.run_pbs) via the same
"local"/"pbs" string every pipeline.py/ts_search.py/irc.py/
verification.py function that runs a Gaussian job now takes as an
`executor` parameter. Swapping one for the other is a single flag
threaded down from webapp.py's/cli.py's session-wide PBS toggle, not a
code change anywhere a Gaussian job actually runs.
"""

from __future__ import annotations

import threading
from typing import Optional

from .local_runner import run_local

VALID_EXECUTORS = ("local", "pbs")


def run_com(com_path: str, executor: str = "local", cancel_event: Optional[threading.Event] = None) -> str:
    """Run `com_path` with whichever executor is selected, both sharing
    the same contract: blocks until done, returns the log path, raises
    GaussianRunError/GaussianCancelledError on failure/cancellation."""
    if executor == "local":
        return run_local(com_path, cancel_event=cancel_event)
    if executor == "pbs":
        from .job_runner import run_pbs

        return run_pbs(com_path, cancel_event=cancel_event)
    raise ValueError(f"Unknown executor {executor!r} -- expected one of {VALID_EXECUTORS}")
