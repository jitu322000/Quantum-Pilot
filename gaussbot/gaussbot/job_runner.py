"""
job_runner.py

Generates the PBS submission script for a .com file -- matching the
site's existing template exactly (Intel setvars, g09.profile, a
per-job GAUSS_SCRDIR under /scratch/$PBS_JOBID, cleaned up after) --
and submits it with qsub, as the queue-based alternative to running
g09 directly on this machine (local_runner.py).

No real PBS system is reachable from here to test submission/polling
against, so this is validated with a fake qsub/qstat script instead
(tests/test_job_runner.py), the same way local_runner.py was
validated before real g09 access was available. Script *generation*
is checked against the site's actual template line for line.

run_pbs() is the queue-based counterpart to local_runner.run_local() --
same contract (blocks until the job is done, returns the log path,
raises GaussianRunError/GaussianCancelledError on failure/cancellation)
so pipeline.py/ts_search.py/irc.py/verification.py can call either one
interchangeably through executor.run_com()'s "local"/"pbs" switch,
wired up end to end via webapp.py's/cli.py's session-wide PBS toggle.
Still can't be validated against a real PBS system from here -- same
fake-qsub/qstat testing as everything else in this module, plus code
review for the parts that genuinely can't be faked (a real cluster's
actual queue behavior).
"""

from __future__ import annotations

import re
import subprocess
import threading
import time
from pathlib import Path
from typing import Optional

from .local_runner import GaussianRunError, GaussianCancelledError

# These four defaults match one real site's PBS setup exactly (Intel
# setvars path, g09root, ppn, walltime) -- confirmed against the site's
# own submission script. They will NOT be right for a different
# cluster: edit these constants (or pass the equivalent keyword to
# run_pbs()/build_pbs_script()) for yours before relying on PBS mode.
DEFAULT_G09ROOT = "/apps/scratch/compile"
DEFAULT_INTEL_SETVARS = "/apps/compilers/intel/oneapi/setvars.sh"
DEFAULT_WALLTIME = "9999:00:00"
DEFAULT_PPN = 30


class PBSError(RuntimeError):
    """Raised when qsub/qstat can't be found or exit with an error."""


def build_pbs_script(
    com_path: str,
    job_name: Optional[str] = None,
    nodes: int = 1,
    ppn: int = DEFAULT_PPN,
    walltime: str = DEFAULT_WALLTIME,
    g09root: str = DEFAULT_G09ROOT,
    intel_setvars: str = DEFAULT_INTEL_SETVARS,
) -> str:
    """
    Build the PBS script text for `com_path`. job_name defaults to the
    .com file's stem, matching the site script's `$input` convention:
    the PBS job name (#PBS -N) and the file g09 is pointed at
    (<job_name>.com) both come from it, so com_path should already be
    named <job_name>.com.

    One deliberate difference from the original bash version: it
    builds ". $g09root/g09/bsd/g09.profile" with a *double*-quoted
    echo, so $g09root there gets expanded immediately by the shell
    generating the script (the one running qsub), not the shell that
    later executes it on the compute node -- it only produces the
    right path if g09root happens to already be set in that outer
    shell's environment. Here g09root is a single explicit parameter
    used for both the export line and the source line, so the two are
    guaranteed to match regardless of what's in the environment.
    """
    stem = Path(com_path).stem
    job_name = job_name or stem

    lines = [
        "#!/bin/bash -l",
        f"#PBS -N {job_name}",
        f"#PBS -l nodes={nodes}:ppn={ppn},walltime={walltime}",
        "#PBS -j oe",
        "#" + "*" * 40,
        "cd $PBS_O_WORKDIR",
        f" source {intel_setvars}",
        f"export g09root={g09root}",
        f". {g09root}/g09/bsd/g09.profile",
        "export GAUSS_SCRDIR=/scratch/$PBS_JOBID",
        "mkdir -p $GAUSS_SCRDIR",
        "#" + "*" * 16 + "  Submit the GAUSSIAN job " + "*" * 34,
        f"g09 {stem}.com",
        "rm -rf $GAUSS_SCRDIR",
        "#" + "*" * 40,
    ]
    return "\n".join(lines) + "\n"


def write_pbs_script(com_path: str, **kwargs) -> str:
    """Write the PBS script alongside com_path (as <stem>.sh) and return its path."""
    sh_path = str(Path(com_path).with_suffix(".sh"))
    with open(sh_path, "w") as f:
        f.write(build_pbs_script(com_path, **kwargs))
    return sh_path


PBS_TEMPLATE_PLACEHOLDER = "__JOB__"


def default_pbs_template() -> str:
    """The starting point for a user-editable PBS script template --
    build_pbs_script()'s usual output, but with the job name/`.com`
    filename left as the PBS_TEMPLATE_PLACEHOLDER (`__JOB__`) instead
    of a real one, so the same edited text works for every job
    submitted this session (see render_pbs_script())."""
    return build_pbs_script(f"{PBS_TEMPLATE_PLACEHOLDER}.com")


def render_pbs_script(template: str, job_name: str) -> str:
    """Substitute PBS_TEMPLATE_PLACEHOLDER in a user-edited PBS script
    template with one job's actual name/`.com`-file stem -- it appears
    twice in the default template (`#PBS -N __JOB__` and
    `g09 __JOB__.com`), and every occurrence is substituted, not just
    the first, in case the user's edit added more."""
    return template.replace(PBS_TEMPLATE_PLACEHOLDER, job_name)


def write_pbs_script_from_template(com_path: str, template: str) -> str:
    """Like write_pbs_script(), but rendering a user-edited template
    (render_pbs_script()) instead of build_pbs_script()'s generated
    text -- for a session-wide custom PBS script the user edited
    themselves (the GUI's text box / the CLI's $EDITOR prompt) rather
    than the module's DEFAULT_* constants."""
    stem = Path(com_path).stem
    sh_path = str(Path(com_path).with_suffix(".sh"))
    with open(sh_path, "w") as f:
        f.write(render_pbs_script(template, stem))
    return sh_path


def submit_job(sh_path: str) -> str:
    """
    Submit a PBS script with `qsub`, run from the script's own
    directory (matching the site script's `cd $PBS_O_WORKDIR` -- qsub
    needs to be invoked from where the .com file actually lives).
    Returns the job id qsub prints back (e.g. "12345.headnode").
    """
    result = subprocess.run(
        ["qsub", Path(sh_path).name],
        cwd=str(Path(sh_path).parent) or ".",
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise PBSError(f"qsub failed on {sh_path}: {result.stderr.strip() or result.stdout.strip()}")
    job_id = result.stdout.strip()
    if not job_id:
        raise PBSError(f"qsub for {sh_path} produced no job id on stdout")
    return job_id


_QSTAT_STATE_RE = re.compile(r"job_state\s*=\s*(\S+)")

# PBS/Torque single-letter job states -> our three-way status.
_STATE_MAP = {
    "Q": "queued",
    "H": "queued",  # held
    "R": "running",
    "E": "running",  # exiting
    "C": "done",
}


def job_status(job_id: str) -> str:
    """
    Poll qstat for a job's status: "queued" | "running" | "done" |
    "unknown". "done" doesn't distinguish success from failure -- as
    with local_runner.py, check the .log file itself for that
    (check_normal_termination / parser.parse_log), since a job can
    leave the queue after crashing just as easily as after finishing.
    Most PBS/Torque setups drop a job from qstat entirely once it's
    fully completed and cleared, which is also treated as "done".
    """
    result = subprocess.run(["qstat", "-f", job_id], capture_output=True, text=True)
    if result.returncode != 0:
        return "done"
    m = _QSTAT_STATE_RE.search(result.stdout)
    if not m:
        return "unknown"
    return _STATE_MAP.get(m.group(1), "unknown")


def wait_for_job(job_id: str, poll_seconds: int = 30) -> str:
    """Block, polling job_status(), until the job is no longer queued/running."""
    while True:
        status = job_status(job_id)
        if status != "queued" and status != "running":
            return status
        time.sleep(poll_seconds)


def cancel_pbs_job(job_id: str) -> None:
    """`qdel` a queued/running job -- the PBS equivalent of run_local()'s
    SIGTERM/SIGKILL, used by run_pbs() when `cancel_event` fires. Best-
    effort: a job that's already finished/left the queue by the time
    this runs just no-ops (qdel on an unknown job id fails harmlessly)."""
    subprocess.run(["qdel", job_id], capture_output=True, text=True)


def run_pbs(
    com_path: str,
    cancel_event: Optional[threading.Event] = None,
    poll_seconds: int = 30,
    **pbs_kwargs,
) -> str:
    """
    Run a .com file by submitting it to PBS instead of running g09
    directly -- the queue-based counterpart to local_runner.run_local(),
    with the same contract (blocks until done, returns the log path,
    raises GaussianRunError on failure) so callers can use either one
    interchangeably.

    If `<the .com file's directory>/pbs_template.sh` exists, it's used
    as a user-edited script template (render_pbs_script()) instead of
    build_pbs_script()'s generated text -- that's how the GUI's text
    box / the CLI's $EDITOR prompt actually take effect, written once
    per job/study rather than threaded through every function that
    might run a Gaussian job. Falls back to write_pbs_script()
    (**pbs_kwargs passed straight to build_pbs_script() --
    nodes/ppn/walltime/g09root/intel_setvars) when no such file exists.

    Qsubs whichever script was written, then polls job_status() every
    `poll_seconds` until the job leaves the queue.

    If `cancel_event` gets set while the job is still queued/running,
    it's qdel'd and GaussianCancelledError is raised -- log_path will
    have whatever partial output Gaussian had written before that,
    same as run_local()'s cancellation.

    Like local_runner.check_normal_termination(), "the job left the
    queue" doesn't by itself mean it succeeded -- callers already check
    the log for that; this only guarantees a log file exists once it
    returns.
    """
    log_path = str(Path(com_path).with_suffix(".log"))
    template_path = Path(com_path).parent / "pbs_template.sh"

    try:
        if template_path.exists():
            sh_path = write_pbs_script_from_template(com_path, template_path.read_text())
        else:
            sh_path = write_pbs_script(com_path, **pbs_kwargs)
        job_id = submit_job(sh_path)
    except PBSError as e:
        raise GaussianRunError(f"PBS submission failed for {com_path}: {e}") from e

    while True:
        status = job_status(job_id)
        if status not in ("queued", "running"):
            break
        if cancel_event is not None and cancel_event.is_set():
            cancel_pbs_job(job_id)
            raise GaussianCancelledError(
                f"Cancelled -- PBS job {job_id} was qdel'd (see {log_path} for whatever it had written)"
            )
        time.sleep(poll_seconds)

    if not Path(log_path).exists():
        raise GaussianRunError(f"PBS job {job_id} left the queue but no {log_path} was ever produced")
    return log_path
