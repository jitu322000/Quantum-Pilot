import os
import stat
import threading
from pathlib import Path

import pytest

from gamessbot.job_runner import (
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
from gamessbot.local_runner import GamessCancelledError, GamessRunError

# Fake qsub/qstat stand-ins, same approach gaussbot's own test_job_runner.py
# uses: exercise the real subprocess wiring without a real PBS system to
# submit to.
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

SAMPLE_INP = " $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG=0 MULT=1 $END\n $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END\n"


def _install_fake(tmp_path, monkeypatch, name, script_text):
    script = tmp_path / name
    script.write_text(script_text)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    return script


def test_build_pbs_script_matches_site_template():
    text = build_pbs_script("jobs/water/water_rhf.inp")
    expected = "\n".join(
        [
            "#!/bin/bash -l",
            "#PBS -N water_rhf",
            "#PBS -l nodes=1:ppn=60,walltime=9999:00:00",
            "#PBS -j oe",
            "#" + "*" * 40,
            "cd $PBS_O_WORKDIR",
            "Project=water_rhf",
            "rm -rf $Project.dat",
            "export SCR=/scratch/$PBS_JOBID",
            "mkdir -p $SCR",
            "#*********** fix the problem of sysv semaphores ******************",
            "#*** This may occur if the GAMESS job terminated abrouptly ******* ",
            "#***************** Load intel ************************************",
            "module load compilers/gcc-9.3.0",
            "module load compilers/openmpi-4.1.1",
            "#****************  Submit the gamess job **********************************",
            "/apps/scratch/compile/gamess_2024/rungms-dev $Project 00 10 10  > $Project.log",
            "#****************  Remove the scratch directory *****************",
            "rm -rf $SCR",
            "#****************************************************************",
        ]
    ) + "\n"
    assert text == expected


def test_build_pbs_script_custom_ppn_and_ncpus():
    text = build_pbs_script("job.inp", ppn=8, ncpus=4)
    assert "#PBS -l nodes=1:ppn=8,walltime=9999:00:00" in text
    assert "rungms-dev $Project 00 4 4  > $Project.log" in text


def test_write_pbs_script_writes_file_next_to_inp(tmp_path):
    inp = tmp_path / "water_rhf.inp"
    inp.write_text(SAMPLE_INP)

    sh_path = write_pbs_script(str(inp))

    assert sh_path == str(tmp_path / "water_rhf.sh")
    assert Path(sh_path).exists()
    assert "#PBS -N water_rhf" in Path(sh_path).read_text()


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
    monkeypatch.setattr("gamessbot.job_runner.time.sleep", lambda _: None)  # no real waiting in tests

    status = wait_for_job("12345.headnode", poll_seconds=0)

    assert status == "done"
    assert counter_file.read_text().strip() == "3"  # took 3 polls to flip R -> C


def test_run_pbs_success(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_PROGRESSING)
    monkeypatch.setenv("FAKE_QSTAT_COUNTER", str(tmp_path / "counter"))
    monkeypatch.setattr("gamessbot.job_runner.time.sleep", lambda _: None)

    inp = tmp_path / "job.inp"
    inp.write_text(SAMPLE_INP)
    # The fake qsub/qstat never actually runs rungms-dev, so the log has to
    # be pre-created here to stand in for what a real PBS job would have
    # produced by the time it leaves the queue.
    (tmp_path / "job.log").write_text("EXECUTION OF GAMESS TERMINATED NORMALLY\n")

    log_path = run_pbs(str(inp), poll_seconds=0)

    assert log_path == str(tmp_path / "job.log")


def test_run_pbs_submission_failure_raises_gamess_run_error(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_FAIL)
    inp = tmp_path / "job.inp"
    inp.write_text(SAMPLE_INP)

    with pytest.raises(GamessRunError):
        run_pbs(str(inp))


def test_run_pbs_missing_log_after_queue_exit_raises(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_STATE_TEMPLATE.format(state="C"))
    inp = tmp_path / "job.inp"
    inp.write_text(SAMPLE_INP)
    # Deliberately no job.log written -- simulates a job that left the
    # queue without ever producing output (e.g. crashed before rungms-dev ran).

    with pytest.raises(GamessRunError):
        run_pbs(str(inp))


def test_run_pbs_cancel_event_qdels_and_raises(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_ALWAYS_RUNNING)
    qdel_log = tmp_path / "qdel.log"
    _install_fake(tmp_path, monkeypatch, "qdel", FAKE_QDEL)
    monkeypatch.setenv("FAKE_QDEL_LOG", str(qdel_log))
    monkeypatch.setattr("gamessbot.job_runner.time.sleep", lambda _: None)

    inp = tmp_path / "job.inp"
    inp.write_text(SAMPLE_INP)

    cancel_event = threading.Event()
    cancel_event.set()  # already set before the first poll -- cancel on the first check

    with pytest.raises(GamessCancelledError):
        run_pbs(str(inp), cancel_event=cancel_event, poll_seconds=0)

    assert qdel_log.read_text().strip() == "12345.headnode"


def test_render_pbs_script_substitutes_every_placeholder():
    rendered = render_pbs_script(default_pbs_template(), "myjob")
    assert "__JOB__" not in rendered
    assert "#PBS -N myjob" in rendered
    assert "Project=myjob" in rendered


def test_run_pbs_uses_template_file_next_to_inp(tmp_path, monkeypatch):
    """A pbs_template.sh sitting next to the .inp file should be used
    (with __JOB__ substituted) instead of the module-default script --
    this is how the GUI text box / CLI $EDITOR edit actually takes
    effect, without threading a template through every caller."""
    _install_fake(tmp_path, monkeypatch, "qsub", FAKE_QSUB_OK)
    _install_fake(tmp_path, monkeypatch, "qstat", FAKE_QSTAT_STATE_TEMPLATE.format(state="C"))

    custom_template = "#!/bin/bash -l\n#PBS -N __JOB__\nProject=__JOB__\necho custom for __JOB__\n"
    (tmp_path / "pbs_template.sh").write_text(custom_template)

    inp = tmp_path / "myjob.inp"
    inp.write_text(SAMPLE_INP)
    (tmp_path / "myjob.log").write_text("EXECUTION OF GAMESS TERMINATED NORMALLY\n")

    log_path = run_pbs(str(inp), poll_seconds=0)

    assert log_path == str(tmp_path / "myjob.log")
    written_sh = (tmp_path / "myjob.sh").read_text()
    assert "echo custom for myjob" in written_sh
    assert "__JOB__" not in written_sh
