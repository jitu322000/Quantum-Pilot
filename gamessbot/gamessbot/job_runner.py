"""
job_runner.py

Generates the PBS submission script for a GAMESS .inp file -- matching
this user's own real PBS/GAMESS submission script exactly (module-loaded
GCC/OpenMPI toolchain, a per-job SCR scratch dir under /scratch/$PBS_JOBID,
cleaned up after) -- and submits it with qsub, as the queue-based
alternative to running rungms directly on this machine (local_runner.py).

Mirrors gaussbot/job_runner.py's design (same PBS_TEMPLATE_PLACEHOLDER /
render_pbs_script() / run_pbs() contract) so the GUI/CLI PBS toggle works
the same way here as it does there; only the script content itself
differs (this project's own real GAMESS site template, not gaussbot's
Gaussian one).

No real PBS system is reachable from here to test submission/polling
against, so this is validated with a fake qsub/qstat script instead
(tests/test_job_runner.py), the same way gaussbot's own job_runner.py
was validated. Script *generation* is checked against the user's actual
template line for line.
"""

from __future__ import annotations

import re
import subprocess
import threading
import time
from pathlib import Path
from typing import Optional, Sequence

from .local_runner import GamessCancelledError, GamessRunError

# These defaults match this user's real site setup exactly (module
# names, rungms-dev path, ppn/ncpus split) -- confirmed against their
# own submission script. They will NOT be right for a different
# cluster: edit these constants (or pass the equivalent keyword to
# run_pbs()/build_pbs_script()) for yours before relying on PBS mode.
DEFAULT_PPN = 60
DEFAULT_WALLTIME = "9999:00:00"
DEFAULT_MODULES = ("compilers/gcc-9.3.0", "compilers/openmpi-4.1.1")
DEFAULT_RUNGMS_DEV_PATH = "/apps/scratch/compile/gamess_2024/rungms-dev"
DEFAULT_NCPUS = 10


class PBSError(RuntimeError):
    """Raised when qsub/qstat can't be found or exit with an error."""


def build_pbs_script(
    inp_path: str,
    job_name: Optional[str] = None,
    nodes: int = 1,
    ppn: int = DEFAULT_PPN,
    walltime: str = DEFAULT_WALLTIME,
    rungms_path: str = DEFAULT_RUNGMS_DEV_PATH,
    ncpus: int = DEFAULT_NCPUS,
    modules: Sequence[str] = DEFAULT_MODULES,
) -> str:
    """
    Build the PBS script text for `inp_path`. job_name defaults to the
    .inp file's stem, matching the site script's `$input`/`$Project`
    convention: the PBS job name (#PBS -N), the `Project=` variable, and
    the file rungms-dev is pointed at (<Project>.log) all come from it,
    so inp_path should already be named <job_name>.inp.

    `rungms-dev` is invoked as `<rungms_path> $Project 00 <ncpus>
    <ncpus>`, same repeated-ncpus shape local_runner.run_gamess() uses
    for local runs -- confirmed against the user's own script, which
    calls it with "00 10 10" even though its own #PBS -l line requests
    ppn=60 (a deliberate site convention of requesting more slots than
    actually used, preserved here verbatim rather than "corrected").
    """
    stem = Path(inp_path).stem
    job_name = job_name or stem

    lines = [
        "#!/bin/bash -l",
        f"#PBS -N {job_name}",
        f"#PBS -l nodes={nodes}:ppn={ppn},walltime={walltime}",
        "#PBS -j oe",
        "#" + "*" * 40,
        "cd $PBS_O_WORKDIR",
        f"Project={stem}",
        "rm -rf $Project.dat",
        "export SCR=/scratch/$PBS_JOBID",
        "mkdir -p $SCR",
        "#*********** fix the problem of sysv semaphores ******************",
        "#*** This may occur if the GAMESS job terminated abrouptly ******* ",
        "#***************** Load intel ************************************",
    ]
    for module in modules:
        lines.append(f"module load {module}")
    lines += [
        "#****************  Submit the gamess job **********************************",
        f"{rungms_path} $Project 00 {ncpus} {ncpus}  > $Project.log",
        "#****************  Remove the scratch directory *****************",
        "rm -rf $SCR",
        "#****************************************************************",
    ]
    return "\n".join(lines) + "\n"


def write_pbs_script(inp_path: str, **kwargs) -> str:
    """Write the PBS script alongside inp_path (as <stem>.sh) and return its path."""
    sh_path = str(Path(inp_path).with_suffix(".sh"))
    with open(sh_path, "w") as f:
        f.write(build_pbs_script(inp_path, **kwargs))
    return sh_path


PBS_TEMPLATE_PLACEHOLDER = "__JOB__"


def default_pbs_template() -> str:
    """The starting point for a user-editable PBS script template --
    build_pbs_script()'s usual output, but with the job name/`.inp`
    filename left as the PBS_TEMPLATE_PLACEHOLDER (`__JOB__`) instead
    of a real one, so the same edited text works for every job
    submitted this session (see render_pbs_script())."""
    return build_pbs_script(f"{PBS_TEMPLATE_PLACEHOLDER}.inp")


def render_pbs_script(template: str, job_name: str) -> str:
    """Substitute PBS_TEMPLATE_PLACEHOLDER in a user-edited PBS script
    template with one job's actual name/`.inp`-file stem -- it appears
    twice in the default template (`#PBS -N __JOB__` and
    `Project=__JOB__`), and every occurrence is substituted, not just
    the first, in case the user's edit added more."""
    return template.replace(PBS_TEMPLATE_PLACEHOLDER, job_name)


def write_pbs_script_from_template(inp_path: str, template: str) -> str:
    """Like write_pbs_script(), but rendering a user-edited template
    (render_pbs_script()) instead of build_pbs_script()'s generated
    text -- for a session-wide custom PBS script the user edited
    themselves (the GUI's text box / the CLI's $EDITOR prompt) rather
    than the module's DEFAULT_* constants."""
    stem = Path(inp_path).stem
    sh_path = str(Path(inp_path).with_suffix(".sh"))
    with open(sh_path, "w") as f:
        f.write(render_pbs_script(template, stem))
    return sh_path


def submit_job(sh_path: str) -> str:
    """
    Submit a PBS script with `qsub`, run from the script's own
    directory (matching the site script's `cd $PBS_O_WORKDIR` -- qsub
    needs to be invoked from where the .inp file actually lives).
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
    (check_normal_termination / parser.parse_rhf_log/parse_cis_log/
    parse_casscf_log), since a job can leave the queue after crashing
    just as easily as after finishing. Most PBS/Torque setups drop a
    job from qstat entirely once it's fully completed and cleared,
    which is also treated as "done".
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
    """`qdel` a queued/running job -- the PBS equivalent of run_gamess()'s
    SIGTERM/SIGKILL, used by run_pbs() when `cancel_event` fires. Best-
    effort: a job that's already finished/left the queue by the time
    this runs just no-ops (qdel on an unknown job id fails harmlessly)."""
    subprocess.run(["qdel", job_id], capture_output=True, text=True)


def run_pbs(
    inp_path: str,
    cancel_event: Optional[threading.Event] = None,
    poll_seconds: int = 30,
    **pbs_kwargs,
) -> str:
    """
    Run a .inp file by submitting it to PBS instead of running rungms
    directly -- the queue-based counterpart to local_runner.run_gamess(),
    with the same contract (blocks until done, returns the log path,
    raises GamessRunError on failure) so callers can use either one
    interchangeably.

    If `<the .inp file's directory>/pbs_template.sh` exists, it's used
    as a user-edited script template (render_pbs_script()) instead of
    build_pbs_script()'s generated text -- that's how the GUI's text
    box / the CLI's $EDITOR prompt actually take effect, written once
    per job/study rather than threaded through every function that
    might run a GAMESS job. Falls back to write_pbs_script()
    (**pbs_kwargs passed straight to build_pbs_script() --
    nodes/ppn/walltime/rungms_path/ncpus/modules) when no such file
    exists.

    Qsubs whichever script was written, then polls job_status() every
    `poll_seconds` until the job leaves the queue.

    If `cancel_event` gets set while the job is still queued/running,
    it's qdel'd and GamessCancelledError is raised -- log_path will
    have whatever partial output GAMESS had written before that, same
    as run_gamess()'s cancellation.

    Unlike local_runner.run_gamess(), no .dat move step is needed here:
    the site's rungms-dev script writes GAMESS's output directly to
    $PBS_O_WORKDIR (the .inp's own directory) -- confirmed by the
    user's own template, which only ever *deletes* a stale
    `$Project.dat` there beforehand, never copies one back from `$SCR`
    afterward.

    Like local_runner.check_normal_termination(), "the job left the
    queue" doesn't by itself mean it succeeded -- callers already check
    the log for that; this only guarantees a log file exists once it
    returns.
    """
    log_path = str(Path(inp_path).with_suffix(".log"))
    template_path = Path(inp_path).parent / "pbs_template.sh"

    try:
        if template_path.exists():
            sh_path = write_pbs_script_from_template(inp_path, template_path.read_text())
        else:
            sh_path = write_pbs_script(inp_path, **pbs_kwargs)
        job_id = submit_job(sh_path)
    except PBSError as e:
        raise GamessRunError(f"PBS submission failed for {inp_path}: {e}") from e

    while True:
        status = job_status(job_id)
        if status not in ("queued", "running"):
            break
        if cancel_event is not None and cancel_event.is_set():
            cancel_pbs_job(job_id)
            raise GamessCancelledError(
                f"Cancelled -- PBS job {job_id} was qdel'd (see {log_path} for whatever it had written)"
            )
        time.sleep(poll_seconds)

    if not Path(log_path).exists():
        raise GamessRunError(f"PBS job {job_id} left the queue but no {log_path} was ever produced")
    return log_path
