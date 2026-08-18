"""
level_select.py

Interactive picks for the GAMESS-specific questions: basis set, charge/
multiplicity, SOSCF vs. DIIS, how many CIS states, and where to find
rungms/scratch on this machine -- same questionary-based pattern as
gaussbot's own level_select.py.
"""

from __future__ import annotations

import os
from typing import Optional, Tuple

import questionary

# GBASIS keywords confirmed against the $BASIS section of the GAMESS
# manual (GBASIS=STO/N21/N31/N311/DZV/TZV/CCn...). Curated to a short
# list of common choices, mirroring gaussbot's BASIS_SETS + Custom...
# pattern rather than exposing every keyword in the manual.
GBASIS_CHOICES = {
    "STO-3G": "GBASIS=STO NGAUSS=3",
    "3-21G": "GBASIS=N21 NGAUSS=3",
    "6-31G": "GBASIS=N31 NGAUSS=6",
    "6-31G(d)": "GBASIS=N31 NGAUSS=6 NDFUNC=1",
    "6-31G(d,p)": "GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1",
    "6-311G(d,p)": "GBASIS=N311 NGAUSS=6 NDFUNC=1 NPFUNC=1",
    "DZV (double zeta valence)": "GBASIS=DZV",
    "TZV (triple zeta valence)": "GBASIS=TZV",
    "cc-pVDZ": "GBASIS=CCD",
    "cc-pVTZ": "GBASIS=CCT",
}
_CUSTOM = "Custom..."

# Overridable via env vars (set by install.sh's generated run_gui.sh/
# run_cli.sh wrappers, or manually) so nothing here is tied to one
# person's machine -- "rungms" alone (bare command name) assumes it's
# on PATH, and the scratch default is under the current user's home,
# not a hardcoded absolute path. The GUI/CLI still let you override
# either per-job regardless of these defaults.
DEFAULT_RUNGMS_PATH = os.environ.get("GAMESSBOT_RUNGMS_PATH", "rungms")
DEFAULT_SCRATCH_DIR = os.environ.get("GAMESSBOT_SCRATCH_DIR", os.path.expanduser("~/gamess_scratch"))
DEFAULT_MEM_MWORDS = 1  # RHF/CIS -- cheap, single-reference
DEFAULT_CASSCF_MEM_MWORDS = 20  # CASSCF/XMCQDPT -- multireference, can be far more memory-hungry


def prompt_for_gbasis() -> str:
    """Returns the exact text to go between "$BASIS" and "$END",
    e.g. "GBASIS=STO NGAUSS=3"."""
    choice = questionary.select("Basis set:", choices=list(GBASIS_CHOICES.keys()) + [_CUSTOM]).ask()
    if choice is None:
        raise KeyboardInterrupt
    if choice == _CUSTOM:
        return (
            questionary.text(
                "GBASIS line (the text between \"$BASIS\" and \"$END\", "
                "e.g. \"GBASIS=N31 NGAUSS=6 NDFUNC=1\"):"
            ).ask()
            or ""
        ).strip()
    return GBASIS_CHOICES[choice]


def prompt_for_charge_mult() -> Tuple[int, int]:
    charge = int((questionary.text("Charge:", default="0").ask() or "0").strip())
    mult = int((questionary.text("Multiplicity:", default="1").ask() or "1").strip())
    return charge, mult


def prompt_for_soscf() -> bool:
    """SOSCF is the real GAMESS default for RHF -- offer it as the
    default choice, with DIIS as the explicit alternative, per your own
    framing ("the user can be asked to turn on either with soscf being
    default"). If SOSCF fails to converge, run_rhf_staged() falls back
    to DIIS automatically regardless of this choice."""
    return bool(
        questionary.confirm("Use SOSCF for the RHF convergence (DIIS otherwise)?", default=True).ask()
    )


def prompt_for_cis() -> Tuple[bool, int]:
    """Whether to run CIS on top of RHF, and if so, how many states
    (NSTATE)."""
    run_cis = bool(questionary.confirm("Run CIS on top of RHF?", default=False).ask())
    if not run_cis:
        return False, 0
    nstate = int((questionary.text("Number of CIS states (NSTATE):", default="5").ask() or "5").strip())
    return True, nstate


def prompt_for_casscf() -> bool:
    """Whether to go on to CASSCF once RHF/CIS are done -- the CIS
    excited states are used to suggest an active space (NMCC/NDOC/NVAL),
    which you'll get to confirm or adjust before it actually runs."""
    return bool(
        questionary.confirm(
            "Run CASSCF on top of RHF/CIS (uses the CIS excitations to suggest an active space)?",
            default=False,
        ).ask()
    )


def prompt_for_casscf_nstate(default: int) -> int:
    return int(
        (questionary.text("Number of CASSCF states (state-averaged):", default=str(default)).ask() or str(default)).strip()
        or default
    )


def prompt_for_xmcqdpt() -> bool:
    """Whether to go on to XMCQDPT (dynamic correlation on top of the
    converged CASSCF reference) once CASSCF is done -- reuses the same
    active space and state weights, no further confirmation needed."""
    return bool(questionary.confirm("Run XMCQDPT on top of the converged CASSCF?", default=False).ask())


def prompt_for_transitn() -> bool:
    """Whether to go on to RUNTYP=TRANSITN (radiative transition
    moments / oscillator strengths, relative to S0) once CASSCF is done
    -- restarts from the same converged CASSCF orbitals XMCQDPT does,
    independent of whether XMCQDPT was also requested."""
    return bool(
        questionary.confirm(
            "Run an oscillator-strength (TRANSITN) calculation on top of the converged CASSCF?", default=False,
        ).ask()
    )


def prompt_for_mo_source() -> str:
    """Where a later active-space/state combination should start its
    CASSCF orbital optimization from -- the previous combination's own
    converged (and possibly MO-reordered) orbitals, or fresh
    closed-shell orbitals, per your request that this be a choice
    rather than always reusing the previous combo's (possibly
    reordered) orbitals. Returns "previous" or "rhf"."""
    choice = questionary.select(
        "Start this combination's CASSCF from:",
        choices=[
            questionary.Choice("The previous combination's optimized orbitals", value="previous"),
            questionary.Choice("Fresh closed-shell (RHF) orbitals", value="rhf"),
        ],
    ).ask()
    return choice or "rhf"


def prompt_for_another_active_space_combo() -> bool:
    """Whether to try another active-space/state combination on the
    same CIS excitations, per your request to compare multiple
    active-space/nstate results side by side -- mirrors gaussbot's own
    "Run another study?" loop pattern."""
    return bool(
        questionary.confirm(
            "Try another active-space/state combination to compare against this one?", default=False,
        ).ask()
    )


def prompt_for_casscf_mem_mwords(default: int = DEFAULT_CASSCF_MEM_MWORDS) -> int:
    """
    Memory for CASSCF/XMCQDPT, asked separately from RHF/CIS's -- a
    multireference calculation over a real active space is far more
    memory-hungry than a single-reference RHF/CIS run, so reusing the
    same MWORDS default for both undersizes the expensive stage (a real
    CAS run can genuinely need MWORDS in the hundreds; this default is
    a conservative starting point for a local workstation, not a hint
    at what a serious active space needs).
    """
    return int(
        (questionary.text("Memory (MWORDS) for CASSCF/XMCQDPT:", default=str(default)).ask() or str(default)).strip()
        or default
    )


def prompt_for_gamess_settings() -> Tuple[str, str, int, int]:
    """
    rungms path / scratch dir / ncpus / mem(mwords) for RHF/CIS --
    configurable rather than hardcoded, pre-filled with your own values
    as defaults. CASSCF/XMCQDPT get their own, separately-asked memory
    setting (prompt_for_casscf_mem_mwords()) since they're typically far
    more memory-hungry than a single-reference RHF/CIS run.
    """
    rungms_path = (
        questionary.text("Path to rungms:", default=DEFAULT_RUNGMS_PATH).ask() or DEFAULT_RUNGMS_PATH
    ).strip()
    scratch_dir = (
        questionary.text("GAMESS scratch directory:", default=DEFAULT_SCRATCH_DIR).ask()
        or DEFAULT_SCRATCH_DIR
    ).strip()
    ncpus = int((questionary.text("Number of CPUs:", default="1").ask() or "1").strip())
    mem_mwords = int(
        (questionary.text("Memory (MWORDS) for RHF/CIS:", default=str(DEFAULT_MEM_MWORDS)).ask()
         or str(DEFAULT_MEM_MWORDS)).strip()
    )
    return rungms_path, scratch_dir, ncpus, mem_mwords


def detect_system_resources() -> Tuple[int, int]:
    """Return (cpu_count, total_mem_gb) for this machine."""
    cpus = os.cpu_count() or 1
    try:
        mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        mem_gb = max(1, int(mem_bytes / (1024**3)))
    except (ValueError, OSError, AttributeError):
        mem_gb = 8
    return cpus, mem_gb


def prompt_for_executor() -> str:
    """
    Ask once per (interactive) session whether to run GAMESS jobs
    directly on this machine or submit them to a PBS queue instead
    (job_runner.py -- its DEFAULT_RUNGMS_DEV_PATH/DEFAULT_PPN/
    DEFAULT_WALLTIME/DEFAULT_MODULES constants match one real site's
    setup; edit them for your own cluster's paths/settings before
    relying on this). A session-wide choice, not per-stage -- every
    stage run afterward (RHF, CIS, CASSCF, XMCQDPT) until gamessbot is
    restarted uses whichever is picked here. Returns "local" or "pbs".
    """
    use_pbs = bool(
        questionary.confirm(
            "Submit GAMESS jobs to a PBS queue instead of running them directly on this "
            "machine? (needs qsub/qstat/qdel on PATH -- you'll get a chance to edit the PBS "
            "script itself, including paths/resources, right after this)",
            default=False,
        ).ask()
    )
    return "pbs" if use_pbs else "local"


def prompt_for_pbs_template() -> Optional[str]:
    """
    Show the default PBS script (job_runner.default_pbs_template() --
    build_pbs_script()'s usual text, with the job name/`.inp` filename
    left as the __JOB__ placeholder) and offer to edit it in $EDITOR
    before it's used for every PBS submission this session -- same
    idea as the GUI's editable text box. Keep __JOB__ intact wherever
    it appears; it's substituted with each job's actual name at
    submission time (job_runner.render_pbs_script()). Returns the
    (possibly edited) template text -- never None, since even
    declining the edit still returns the default text explicitly
    rather than falling through to job_runner.py's own defaults
    silently.
    """
    from .job_runner import default_pbs_template

    default_template = default_pbs_template()
    print("\nDefault PBS script template for this session (edit paths/resources for your cluster):")
    print(default_template)
    edit = questionary.confirm(
        "Edit this template in $EDITOR before it's used for every PBS submission this session?",
        default=True,
    ).ask()
    if not edit:
        return default_template

    import subprocess
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as f:
        f.write(default_template)
        tmp_path = f.name
    editor = os.environ.get("EDITOR", "nano")
    subprocess.call([editor, tmp_path])
    with open(tmp_path) as f:
        edited = f.read()
    os.unlink(tmp_path)
    return edited
