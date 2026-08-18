"""
level_select.py

Interactive picks for the questions that come up around the actual
Gaussian study: what kind of run this is, which method/basis set (or
GenECP mixed basis) to use, how much of this machine to give
Gaussian, and whether/how to run the IRC.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import questionary

from .input_builder import BasisGroups

METHODS = ["HF", "B3LYP", "M06-2X", "wB97XD", "MP2"]
BASIS_SETS = [
    "STO-3G",
    "3-21G",
    "6-31G(d)",
    "6-31+G(d,p)",
    "6-311G(d,p)",
    "6-311+G(d,p)",
    "cc-pVDZ",
    "cc-pVTZ",
    "def2-TZVP",
]
_CUSTOM = "Custom..."
_MIXED_BASIS = "Mixed basis set (GenECP) -- different basis per element"


@dataclass
class LevelChoice:
    """What prompt_for_method_basis() came back with -- either a plain
    basis set, or a mixed (GenECP) one (basis_groups/ecp_groups set,
    basis left blank)."""

    method: str
    basis: str = ""
    basis_groups: Optional[BasisGroups] = None
    ecp_groups: Optional[BasisGroups] = None


def prompt_for_genecp_assignment(elements: List[str]) -> Tuple[BasisGroups, BasisGroups]:
    """
    Mixed (GenECP) basis set: ask which elements need an effective core
    potential (heavy atoms / transition metals -- e.g. LanL2DZ, SDD),
    what name to use for each, and what conventional basis set covers
    everything else. Elements sharing the same ECP name are grouped
    together automatically, matching Gaussian's basis/pseudopotential
    block format. Returns (basis_groups, ecp_groups) ready for
    input_builder.build_input.
    """
    elements = sorted(set(elements))
    ecp_elements = (
        questionary.checkbox(
            "Which elements need an effective core potential (heavy atoms / "
            "transition metals)? Space to select, enter to confirm -- none "
            "selected just uses one basis set for everything.",
            choices=elements,
        ).ask()
        or []
    )

    light_elements = [e for e in elements if e not in ecp_elements]
    basis_groups: BasisGroups = []
    if light_elements:
        light_basis = questionary.select(
            f"Basis set for {', '.join(light_elements)}:", choices=BASIS_SETS + [_CUSTOM]
        ).ask()
        if light_basis == _CUSTOM:
            light_basis = (questionary.text("Basis set:").ask() or "").strip()
        basis_groups.append((light_elements, light_basis))

    ecp_groups: BasisGroups = []
    names: Dict[str, List[str]] = {}
    for el in ecp_elements:
        name = (
            questionary.text(
                f"ECP/basis name for {el} (e.g. LanL2DZ, SDD, def2-TZVP):", default="LanL2DZ"
            ).ask()
            or "LanL2DZ"
        ).strip()
        names.setdefault(name, []).append(el)
    for name, els in names.items():
        basis_groups.append((els, name))
        ecp_groups.append((els, name))

    return basis_groups, ecp_groups


def prompt_for_method_basis(elements: Optional[List[str]] = None) -> LevelChoice:
    """
    GaussView-style pick: a method off a short list (or type your own),
    then a basis set the same way -- or, if `elements` is given (the
    distinct elements actually in the structures being studied), the
    option to assign a mixed GenECP basis per element instead.
    """
    method = questionary.select("Method for the final optimization:", choices=METHODS + [_CUSTOM]).ask()
    if method is None:
        raise KeyboardInterrupt
    if method == _CUSTOM:
        method = (questionary.text("Method (e.g. B3LYP, M062X, MP2):").ask() or "").strip()

    basis_choices = BASIS_SETS + [_CUSTOM]
    if elements:
        basis_choices = basis_choices + [_MIXED_BASIS]

    basis = questionary.select("Basis set:", choices=basis_choices).ask()
    if basis is None:
        raise KeyboardInterrupt

    if basis == _MIXED_BASIS:
        basis_groups, ecp_groups = prompt_for_genecp_assignment(elements)
        return LevelChoice(method=method, basis_groups=basis_groups, ecp_groups=ecp_groups)

    if basis == _CUSTOM:
        basis = (questionary.text("Basis set (e.g. 6-31G(d)):").ask() or "").strip()

    return LevelChoice(method=method, basis=basis)


def detect_system_resources() -> Tuple[int, int]:
    """Return (cpu_count, total_mem_gb) for this machine."""
    cpus = os.cpu_count() or 1
    try:
        mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        mem_gb = max(1, int(mem_bytes / (1024**3)))
    except (ValueError, OSError, AttributeError):
        mem_gb = 8  # couldn't detect -- a conservative guess
    return cpus, mem_gb


def prompt_for_resources(default_nprocs: int = 2, default_mem_gb: int = 2) -> Tuple[int, int]:
    """
    Default %nprocshared/%mem to a conservative (2, 2) for now, with
    the option to ask for more. This used to default to 1/4 of the
    machine's resources, but that's more than needed while still
    validating the pipeline against real reactions -- (2, 2) keeps
    test runs cheap and out of each other's way; bump the defaults
    once you trust it enough to actually want more.
    """
    total_cpus, total_mem_gb = detect_system_resources()

    use_default = questionary.confirm(
        f"Use default resources -- %nprocshared={default_nprocs}, %mem={default_mem_gb}GB "
        f"(this machine has {total_cpus} cores / {total_mem_gb}GB total)?",
        default=True,
    ).ask()
    if use_default or use_default is None:
        return default_nprocs, default_mem_gb

    nprocs = int(
        questionary.text(
            f"nprocshared (machine has {total_cpus} cores):", default=str(default_nprocs)
        ).ask()
        or default_nprocs
    )
    mem_gb = int(
        questionary.text(
            f"mem in GB (machine has {total_mem_gb}GB):", default=str(default_mem_gb)
        ).ask()
        or default_mem_gb
    )
    return nprocs, mem_gb


def prompt_for_skip_pm6() -> bool:
    """
    PM6 pre-optimization is the default -- it's cheap and gets a
    structure into a reasonable starting shape before the expensive
    final-level optimization. Ask whether to skip it and start directly
    at the final level instead, for whenever the input geometry is
    already trustworthy enough not to need it.
    """
    return bool(
        questionary.confirm(
            "Skip PM6 pre-optimization and start directly at the final level you'll choose next?",
            default=False,
        ).ask()
    )


def prompt_for_xyz_export() -> bool:
    """Ask whether to convert the final reactant/product/TS (or single
    structure) logs to .xyz via Open Babel once the study finishes --
    handy for dropping straight into a paper's SI."""
    return bool(
        questionary.confirm(
            "Convert the final structure(s) to .xyz files once the study finishes (needs Open Babel)?",
            default=False,
        ).ask()
    )


def prompt_for_cdxml_export() -> bool:
    """Ask whether to convert the final reactant/product/TS (or single
    structure) logs to .cdxml (ChemDraw) via Open Babel once the study
    finishes -- same idea as prompt_for_xyz_export(), a different
    format for dropping into a paper or a ChemDraw figure."""
    return bool(
        questionary.confirm(
            "Convert the final structure(s) to .cdxml (ChemDraw) files once the study finishes (needs Open Babel)?",
            default=False,
        ).ask()
    )


def prompt_for_gcorr() -> Optional[Tuple[float, str]]:
    """
    Ask whether to compute G-corr (thermal correction to Gibbs free
    energy, via Gaussian's freqchk) for each structure independently --
    no reference subtraction, just the raw correction per molecule.
    Returns (temperature, unit) if wanted, None to skip.
    """
    want = questionary.confirm(
        "Compute G-corr (thermal correction to Gibbs free energy) for each structure?",
        default=False,
    ).ask()
    if not want:
        return None
    unit = questionary.select("Temperature unit:", choices=["K", "C", "F"]).ask() or "K"
    default_temp = "298.15" if unit == "K" else "25" if unit == "C" else "77"
    temp = questionary.text(f"Temperature ({unit}):", default=default_temp).ask()
    try:
        return float(temp), unit
    except (TypeError, ValueError):
        return 298.15, "K"


def prompt_for_irc() -> Optional[int]:
    """
    Ask whether to run the IRC verification step, and if so how many
    points per direction (default 25). Returns None to skip it.
    """
    want_irc = questionary.confirm(
        "Run the IRC to verify the TS actually connects the reactant and product?",
        default=True,
    ).ask()
    if not want_irc:
        return None
    points = questionary.text("IRC points per direction:", default="25").ask()
    try:
        return int(points)
    except (TypeError, ValueError):
        return 25


def prompt_for_max_repair_attempts(default: int = 4) -> int:
    """
    How many times the repair loop (pipeline.repair_and_optimize) may
    retry a non-converged or imaginary-frequency-flagged optimization
    before giving up. Higher isn't free -- each attempt is a full Opt
    Freq -- but a geometry that's genuinely close to converging can
    need more than the default to get there.
    """
    answer = questionary.text(
        "Max repair attempts per optimization (non-convergence retries, imaginary-frequency kicks):",
        default=str(default),
    ).ask()
    try:
        return max(1, int(answer))
    except (TypeError, ValueError):
        return default


def prompt_for_endpoint_recovery() -> bool:
    """
    Ask whether to keep going if the reactant or product can't be
    optimized after every fallback, instead of aborting the study --
    the missing side is then recovered from the TS (once found) via
    TS-mode distortion and/or IRC-endpoint reopt, if either is also
    enabled. Off by default, matching this study's other opt-in
    cross-checks.
    """
    return bool(
        questionary.confirm(
            "If the reactant or product fails to optimize, keep going and try to recover it "
            "from the TS afterward (via TS-mode distortion / IRC-endpoint reopt) instead of "
            "aborting the study?",
            default=False,
        ).ask()
    )


def prompt_for_ts_distortion_check() -> bool:
    """
    Ask whether to also cross-check the reactant/product by displacing
    the converged TS along its own imaginary mode (+ and -) and
    reoptimizing each side -- an independent regeneration of the
    reactant/product straight from the TS, compared against the
    guess-based ones by energy. Off by default: adds two extra Opt
    Freq jobs at the final level.
    """
    return bool(
        questionary.confirm(
            "Also cross-check the reactant/product by distorting the TS along its "
            "imaginary mode and reoptimizing both sides? (adds 2 more Opt Freq jobs)",
            default=False,
        ).ask()
    )


def prompt_for_irc_endpoint_reopt() -> bool:
    """
    Ask whether to reoptimize the IRC's forward/reverse endpoints to
    proper minima (the raw endpoint is just the last point of a
    fixed-step path) -- a third candidate reactant/product, compared
    the same way. Only meaningful once the IRC has actually run. Off
    by default: adds two extra Opt Freq jobs at the final level.
    """
    return bool(
        questionary.confirm(
            "Also reoptimize the IRC's forward/reverse endpoints to proper minima? "
            "(adds 2 more Opt Freq jobs)",
            default=False,
        ).ask()
    )


def prompt_for_executor() -> str:
    """
    Ask once per (interactive) session whether to run Gaussian jobs
    directly on this machine or submit them to a PBS queue instead
    (job_runner.py -- its DEFAULT_G09ROOT/DEFAULT_INTEL_SETVARS/
    DEFAULT_PPN/DEFAULT_WALLTIME constants match one real site's setup;
    edit them for your own cluster's paths/settings before relying on
    this). A session-wide choice, not per-study -- every study run
    afterward (until gaussbot is restarted, or --executor overrides it)
    uses whichever is picked here. Returns "local" or "pbs".
    """
    use_pbs = bool(
        questionary.confirm(
            "Submit Gaussian jobs to a PBS queue instead of running them directly on this "
            "machine? (needs qsub/qstat/qdel on PATH -- you'll get a chance to edit the PBS "
            "script itself, including paths/resources, right after this)",
            default=False,
        ).ask()
    )
    return "pbs" if use_pbs else "local"


def prompt_for_pbs_template() -> Optional[str]:
    """
    Show the default PBS script (job_runner.default_pbs_template() --
    build_pbs_script()'s usual text, with the job name/`.com` filename
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


_MECHANISM = "Reaction mechanism study (reactant + product + TS)"
_SINGLE = "Single geometry optimization"
_TS_SEARCH = "Transition state search (TS only, no forced IRC/energetics)"


def prompt_for_study_type() -> str:
    """
    Ask what this run is for: a full reaction-mechanism study
    (reactant + product + TS + IRC + energetics), a single geometry
    optimization -- e.g. to hand off to a separate single-point
    energy, NMR, or scan calculation later -- or a dedicated TS search
    (just the saddle point itself, via search_ts_staged() in
    ts_search.py, without committing to a full mechanism study).
    Returns "single", "mechanism", or "ts_search".
    """
    choice = questionary.select("What do you want to do?", choices=[_MECHANISM, _SINGLE, _TS_SEARCH]).ask()
    if choice == _SINGLE:
        return "single"
    if choice == _TS_SEARCH:
        return "ts_search"
    return "mechanism"


def prompt_for_ts_search_skip_pm6() -> bool:
    """Same idea as prompt_for_skip_pm6(), for the dedicated TS Search
    section's own PM6 stage (search_ts_staged() in ts_search.py)."""
    return bool(
        questionary.confirm(
            "Skip the PM6 pre-stage and search directly at the final level you'll choose next?",
            default=False,
        ).ask()
    )


def prompt_for_ts_calcall_fallback() -> bool:
    """Ask whether search_ts_staged() should fall back to the CalcAll
    (full-Hessian-every-step) route when the usual CalcFC route can't
    locate a saddle point -- off by default, since CalcAll is much
    more expensive and not every user wants to spend the resources on
    it automatically."""
    return bool(
        questionary.confirm(
            "If the usual (CalcFC) search can't find the TS, fall back to the more "
            "expensive CalcAll (full-Hessian) brute-force route?",
            default=False,
        ).ask()
    )


def prompt_for_ts_recovery() -> bool:
    """Ask whether to try recovering the reactant/product straight from
    the found TS (verification.recover_endpoints_from_ts()) once the
    search succeeds."""
    return bool(
        questionary.confirm(
            "Once the TS is found, also try to recover the reactant/product from it "
            "(displace along the imaginary mode, reoptimize both sides)?",
            default=False,
        ).ask()
    )
