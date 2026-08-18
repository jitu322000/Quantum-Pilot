import os
import stat
import threading
from pathlib import Path

import pytest

from gaussbot.job_runner import (
    PBSError,
    build_pbs_script,
    default_pbs_template,
    job_status,
    render_pbs_script,
    run_pbs,
    submit_job,
    wait_for_job,
    write_pbs_script,
)
from gaussbot.local_runner import GaussianCancelledError, GaussianRunError

# Fake qsub/qstat stand-ins, same approach as test_local_runner.py's
# fake g09: exercise the real subprocess wiring without a real PBS
# system to submit to.
FAKE_QSUB_OK = """#!/bin/sh
echo "12345.headnode"
"""

FAKE_QSUB_FAIL = """#!/bin/sh
echo "qsub: submit error" >&2
exit 1
"""

FAKE_QSTAT_STATE_TEMPLATE = """#!/bin/sh
echo "job_state = {state}"
"""

FAKE_QSTAT_MISSING = """#!/bin/sh
echo "qstat: Unknown Job Id" >&2
exit 153
"""

# Ticks through R, R, C on successive calls -- for wait_for_job.
FAKE_QSTAT_PROGRESSING = """#!/bin/sh
n=$(cat "$FAKE_QSTAT_COUNTER" 2>/dev/null || echo 0)
n=$((n + 1))
echo $n > "$FAKE_QSTAT_COUNTER"
if [ "$n" -lt 3 ]; then
    echo "job_state = R"
else
    echo "job_state = C"
fi
"""

FAKE_QSTAT_ALWAYS_RUNNING = """#!/bin/sh
echo "job_state = R"
"""

FAKE_QDEL = """#!/bin/sh
echo "$1" >> "$FAKE_QDEL_LOG"
"""

SAMPLE_COM = "%chk=job.chk\n#p PM6 Opt\n\ntitle\n\n0 1\nH 0.0 0.0 0.0\n\n"


def _install_fake(tmp_path, monkeypatch, name, script_text):
    script = tmp_path / name
    script.write_text(script_text)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    return script


def test_build_pbs_script_matches_site_template():
    text = build_pbs_script("jobs/myrxn/myrxn_ts.com", g09root="/apps/scratch/compile")
    expected = "\n".join(
        [
            "#!/bin/bash -l",
            "#PBS -N myrxn_ts",
            "#PBS -l nodes=1:ppn=30,walltime=9999:00:00",
            "#PBS -j oe",
            "#" + "*" * 40,
            "cd $PBS_O_WORKDIR",
            " source /apps/compilers/intel/oneapi/setvars.sh",
            "export g09root=/apps/scratch/compile",
            ". /apps/scratch/compile/g09/bsd/g09.profile",
            "export GAUSS_SCRDIR=/scratch/$PBS_JOBID",
            "mkdir -p $GAUSS_SCRDIR",
            "#" + "*" * 16 + "  Submit the GAUSSIAN job " + "*" * 34,
            "g09 myrxn_ts.com",
            "rm -rf $GAUSS_SCRDIR",
            "#" + "*" * 40,
        ]
    ) + "\n"
    assert text == expected


def test_build_pbs_script_custom_ppn_and_walltime():
    text = build_pbs_script("job.com", ppn=8, walltime="24:00:00")
    assert "#PBS -l nodes=1:ppn=8,walltime=24:00:00" in text


def test_write_pbs_script_writes_file_next_to_com(tmp_path):
    com = tmp_path / "myrxn_ts.com"
    com.write_text(SAMPLE_COM)

    sh_path = write_pbs_script(str(com))

    assert sh_path == str(tmp_path / "myrxn_ts.sh")
    assert Path(sh_path).exists()
    assert "#PBS -N myrxn_ts" in Path(sh_path).read_text()


def test_submit_job_returns_job_id(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    sh = tmp_path / "job.sh"
    sh.write_text("#!/bin/bash -l\necho hi\n")

    job_id = submit_job(str(sh))

    assert job_id == "12345.headnode"


def test_submit_job_nonzero_exit_raises(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_FAIL)
    sh = tmp_path / "job.sh"
    sh.write_text("#!/bin/bash -l\necho hi\n")

    with pytest.raises(PBSError):
        submit_job(str(sh))


@pytest.mark.parametrize(
    "pbs_state,expected",
    [("Q", "queued"), ("H", "queued"), ("R", "running"), ("E", "running"), ("C", "done")],
)
def test_job_status_maps_pbs_states(tmp_path, monkeypatch, pbs_state, expected):
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_STATE_TEMPLATE.format(state=pbs_state))
    assert job_status("12345.headnode") == expected


def test_job_status_missing_job_is_done(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_MISSING)
    assert job_status("12345.headnode") == "done"


def test_wait_for_job_polls_until_done(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_PROGRESSING)
    counter_file = tmp_path / "counter"
    monkeypatch.setenv("FAKE_QSTAT_COUNTER", str(counter_file))
    monkeypatch.setattr("gaussbot.job_runner.time.sleep", lambda _: None)  # no real waiting in tests

    status = wait_for_job("12345.headnode", poll_seconds=0)

    assert status == "done"
    assert counter_file.read_text().strip() == "3"  # took 3 polls to flip R -> C


def test_run_pbs_success(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_PROGRESSING)
    monkeypatch.setenv("FAKE_QSTAT_COUNTER", str(tmp_path / "counter"))
    monkeypatch.setattr("gaussbot.job_runner.time.sleep", lambda _: None)

    com = tmp_path / "job.com"
    com.write_text(SAMPLE_COM)
    # The fake qsub/qstat never actually runs g09, so the log has to be
    # pre-created here to stand in for what a real PBS job would have
    # produced by the time it leaves the queue.
    (tmp_path / "job.log").write_text("Normal termination of Gaussian\n")

    log_path = run_pbs(str(com), poll_seconds=0, g09root="/apps/scratch/compile")

    assert log_path == str(tmp_path / "job.log")


def test_run_pbs_submission_failure_raises_gaussian_run_error(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_FAIL)
    com = tmp_path / "job.com"
    com.write_text(SAMPLE_COM)

    with pytest.raises(GaussianRunError):
        run_pbs(str(com))


def test_run_pbs_missing_log_after_queue_exit_raises(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_STATE_TEMPLATE.format(state="C"))
    com = tmp_path / "job.com"
    com.write_text(SAMPLE_COM)
    # Deliberately no job.log written -- simulates a job that left the
    # queue without ever producing output (e.g. crashed before g09 ran).

    with pytest.raises(GaussianRunError):
        run_pbs(str(com))


def test_run_pbs_cancel_event_qdels_and_raises(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_ALWAYS_RUNNING)
    qdel_log = tmp_path / "qdel.log"
    _install_fake(tmp_path, monkeypatch, "qdel", FAKE_QDEL)
    monkeypatch.setenv("FAKE_QDEL_LOG", str(qdel_log))
    monkeypatch.setattr("gaussbot.job_runner.time.sleep", lambda _: None)

    com = tmp_path / "job.com"
    com.write_text(SAMPLE_COM)

    cancel_event = threading.Event()
    cancel_event.set()  # already set before the first poll -- cancel on the first check

    with pytest.raises(GaussianCancelledError):
        run_pbs(str(com), cancel_event=cancel_event, poll_seconds=0)

    assert qdel_log.read_text().strip() == "12345.headnode"


def test_render_pbs_script_substitutes_every_placeholder():
    rendered = render_pbs_script(default_pbs_template(), "myrxn_ts")
    assert "__JOB__" not in rendered
    assert "#PBS -N myrxn_ts" in rendered
    assert "g09 myrxn_ts.com" in rendered


def test_run_pbs_uses_template_file_next_to_com(tmp_path, monkeypatch):
    """A pbs_template.sh sitting next to the .com file should be used
    (with __JOB__ substituted) instead of the module-default script --
    this is how the GUI text box / CLI $EDITOR edit actually takes
    effect, without threading a template through every caller."""
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_STATE_TEMPLATE.format(state="C"))

    custom_template = "#!/bin/bash -l\n#PBS -N __JOB__\necho custom for __JOB__\ng09 __JOB__.com\n"
    (tmp_path / "pbs_template.sh").write_text(custom_template)

    com = tmp_path / "myjob.com"
    com.write_text(SAMPLE_COM)
    (tmp_path / "myjob.log").write_text("Normal termination of Gaussian\n")

    log_path = run_pbs(str(com), poll_seconds=0)

    assert log_path == str(tmp_path / "myjob.log")
    written_sh = (tmp_path / "myjob.sh").read_text()
    assert "echo custom for myjob" in written_sh
    assert "__JOB__" not in written_sh
