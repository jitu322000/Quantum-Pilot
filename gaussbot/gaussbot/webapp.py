"""
webapp.py

A local web GUI wrapped around the existing pipeline -- cli.py is the
only other place that talks to a human, and it does that through
questionary/input(); this does it through a browser instead. Every
actual computation still goes through the same functions cli.py
calls (pipeline.py, ts_search.py, irc.py, energetics.py) completely
unchanged -- this file only orchestrates them: accepts a job over
HTTP, runs it on a background thread (so the request returns
immediately), and lets the frontend poll for progress.

Run: gaussbot-gui   (opens a browser tab at http://127.0.0.1:8765)
"""

from __future__ import annotations

import os
import threading
import time
import uuid
import webbrowser
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import uvicorn
from fastapi import FastAPI, HTTPException, UploadFile
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .intake import load_structure, check_reaction_match
from .level_select import LevelChoice
from .local_runner import GaussianCancelledError
from .pipeline import OptOutcome, preopt_with_escalation, repair_and_optimize
from .ts_search import run_ts_search, search_ts_staged
from .irc import run_irc
from .verification import (
    run_ts_distortion_check, run_irc_endpoint_reopt, select_best_candidate,
    format_candidate_comparison, recover_endpoints_from_ts, DistortionOutcome,
)
from .energetics import build_energy_report, format_energy_report, HARTREE_TO
from .geometry import Structure, from_smiles, from_pubchem
from .thermochem import run_freqchk, to_kelvin
from .xyz_export import convert_log_to_xyz, convert_log_to_cdxml
from . import llm_assist

WEB_DIR = Path(__file__).parent / "web"
UPLOAD_DIR = Path(os.environ.get("GAUSSBOT_UPLOAD_DIR", "/tmp/gaussbot_uploads"))
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="gaussbot")
app.mount("/static", StaticFiles(directory=str(WEB_DIR)), name="static")


@app.middleware("http")
async def _no_cache_static(request, call_next):
    """Without this, browsers can serve /static/style.css or app.js
    from their own heuristic cache on a plain reload (no Cache-Control
    header means the browser is free to guess a freshness lifetime).
    Forcing revalidation on every request means a plain reload always
    picks up the latest file (StaticFiles's own ETag/Last-Modified
    still make that cheap -- a 304 when nothing changed)."""
    response = await call_next(request)
    if request.url.path.startswith("/static/"):
        response.headers["Cache-Control"] = "no-cache"
    return response

# token -> (Structure, saved file path) for uploaded files -- lets the
# frontend ask "what elements are in this file?" (for the GenECP
# picker) before a job is actually submitted.
_UPLOADS: Dict[str, "_Upload"] = {}


@dataclass
class _Upload:
    structure: Structure
    path: str


@dataclass
class JobState:
    id: str
    status: str = "running"  # running | done | error | cancelled
    stage: str = "starting"
    log: List[str] = field(default_factory=list)
    result: Optional[dict] = None
    error: Optional[str] = None
    cancel_event: threading.Event = field(default_factory=threading.Event)
    req: Optional["JobRequest"] = None  # the request that started this job -- replayed on restart
    # stage_name -> the OptOutcome/TSOutcome that stage produced, so a restart after
    # cancellation can skip straight past whatever already finished cleanly rather
    # than re-running it from scratch. Only the expensive, deterministic stages
    # (pre-opt, final opt, TS search) are checkpointed -- the opt-in extras (IRC,
    # TS-distortion, IRC-endpoint reopt) are cheap enough by comparison to just
    # re-run on restart rather than also threading resume logic through them.
    checkpoints: dict = field(default_factory=dict)


_JOBS: Dict[str, JobState] = {}
_JOBS_LOCK = threading.Lock()


def _log(job: JobState, message: str) -> None:
    with _JOBS_LOCK:
        job.log.append(message)


def _set_stage(job: JobState, stage: str) -> None:
    with _JOBS_LOCK:
        job.stage = stage
        job.log.append(f"--- {stage} ---")


def _fail(job: JobState, message: str) -> None:
    with _JOBS_LOCK:
        job.status = "error"
        job.error = message
        job.log.append(message)


# ---------------------------------------------------------------- API models

class GenECPChoice(BaseModel):
    ecp_assignments: Dict[str, str] = {}  # element -> ECP/basis name
    light_basis: str = ""  # basis for every element not in ecp_assignments


class LevelRequest(BaseModel):
    method: str
    basis: str = ""
    genecp: Optional[GenECPChoice] = None


class JobRequest(BaseModel):
    study_type: str  # "single" | "mechanism" | "ts_search"
    job_name: str
    reactant_token: Optional[str] = None
    product_token: Optional[str] = None
    ts_guess_token: Optional[str] = None
    reactant_id: Optional[str] = None
    product_id: Optional[str] = None
    nprocs: int = 2
    mem_gb: int = 2
    max_repair_attempts: int = 4
    level: LevelRequest
    run_irc: bool = True
    irc_points: int = 25
    energy_units: List[str] = ["kcal/mol"]
    run_ts_distortion: bool = False
    distortion_factor: float = 1.0
    run_irc_endpoint_reopt: bool = False
    recover_missing_endpoint: bool = False
    skip_pm6_preopt: bool = False
    ts_recover_endpoints: bool = False
    ts_use_calcall: bool = False
    executor: str = "local"  # "local" (run g09 directly) or "pbs" (submit via qsub, see job_runner.py)
    pbs_script_template: Optional[str] = None  # user-edited PBS script (see job_runner.render_pbs_script) -- only used when executor == "pbs"
    export_xyz: bool = False
    export_cdxml: bool = False
    compute_gcorr: bool = False
    gcorr_temperature: float = 298.15
    gcorr_temperature_unit: str = "K"


def _level_choice(level: LevelRequest, elements: List[str]) -> LevelChoice:
    if level.genecp is None:
        return LevelChoice(method=level.method, basis=level.basis)

    ecp_elements = list(level.genecp.ecp_assignments.keys())
    light_elements = [e for e in elements if e not in ecp_elements]
    basis_groups = []
    if light_elements:
        basis_groups.append((light_elements, level.genecp.light_basis))

    ecp_groups = []
    by_name: Dict[str, List[str]] = {}
    for element, name in level.genecp.ecp_assignments.items():
        by_name.setdefault(name, []).append(element)
    for name, els in by_name.items():
        basis_groups.append((els, name))
        ecp_groups.append((els, name))

    return LevelChoice(method=level.method, basis_groups=basis_groups, ecp_groups=ecp_groups)


# ------------------------------------------------------------- job runners

def _checkpoint(resume_from: Optional[JobState], name: str):
    """If the previous (cancelled) run already finished this stage
    cleanly, return its outcome so the caller can skip re-running it;
    None if there's nothing to reuse."""
    if resume_from is None:
        return None
    outcome = resume_from.checkpoints.get(name)
    return outcome if outcome is not None and outcome.success else None


def _maybe_write_pbs_template(out_dir: str, req: JobRequest) -> None:
    """If PBS mode is on and the user edited a script template (GUI
    text box), write it once as <out_dir>/pbs_template.sh --
    job_runner.run_pbs() picks it up automatically from there for
    every job this study submits, so this is the only place the
    template needs to be threaded through, not every pipeline call."""
    if req.executor == "pbs" and req.pbs_script_template:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "pbs_template.sh"), "w") as f:
            f.write(req.pbs_script_template)


def _run_single_job(job: JobState, req: JobRequest, structure: Structure, resume_from: Optional[JobState] = None) -> None:
    out_dir = os.path.join("jobs", req.job_name)
    _maybe_write_pbs_template(out_dir, req)

    pre = None
    if req.skip_pm6_preopt:
        _log(job, "Skipping PM6 pre-optimization -- starting directly at the final level.")
    else:
        pre = _checkpoint(resume_from, "preopt")
        if pre is not None:
            _log(job, f"Reusing already-completed pre-optimization from before the restart ({pre.log_path}).")
        else:
            _set_stage(job, f"Pre-optimizing {req.job_name} (PM6)")
            pre = preopt_with_escalation(
                structure, out_dir, "structure", pubchem_query=req.reactant_id, nprocs=req.nprocs, mem_gb=req.mem_gb,
                max_repair_attempts=req.max_repair_attempts, cancel_event=job.cancel_event, executor=req.executor,
            )
            for line in pre.log:
                _log(job, line)
        job.checkpoints["preopt"] = pre
        if not pre.success:
            _fail(job, f"Couldn't get {req.job_name} to a clean PM6 minimum after every fallback.")
            return

    starting_structure = pre.structure if pre else structure
    elements = sorted({a[0] for a in starting_structure.atoms})
    level = _level_choice(req.level, elements)

    final = _checkpoint(resume_from, "final")
    if final is not None:
        _log(job, f"Reusing already-completed final-level optimization from before the restart ({final.log_path}).")
    else:
        _set_stage(job, "Optimizing at the final level")
        final = repair_and_optimize(
            starting_structure, os.path.join(out_dir, "structure_final.com"),
            method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
            nprocs=req.nprocs, mem_gb=req.mem_gb, max_repair_attempts=req.max_repair_attempts,
            cancel_event=job.cancel_event, executor=req.executor,
        )
        for line in final.log:
            _log(job, line)
    job.checkpoints["final"] = final
    if not final.success:
        _fail(job, f"{req.job_name} didn't reach a clean minimum at the final level.")
        return

    xyz_files = {}
    if req.export_xyz:
        xyz_path = convert_log_to_xyz(final.log_path)
        if xyz_path:
            xyz_files["final"] = xyz_path
        else:
            _log(job, f"Couldn't convert {final.log_path} to .xyz (is obabel installed?).")

    cdxml_files = {}
    if req.export_cdxml:
        cdxml_path = convert_log_to_cdxml(final.log_path)
        if cdxml_path:
            cdxml_files["final"] = cdxml_path
        else:
            _log(job, f"Couldn't convert {final.log_path} to .cdxml (is obabel installed?).")

    g_corr_hartree = None
    if req.compute_gcorr:
        temp_k = to_kelvin(req.gcorr_temperature, req.gcorr_temperature_unit)
        g_corr_hartree = run_freqchk(final.log_path.replace(".log", ".chk"), temp_k)
        if g_corr_hartree is None:
            _log(job, f"Couldn't compute G-corr from {final.log_path.replace('.log', '.chk')} (is freqchk installed/working?).")

    with _JOBS_LOCK:
        job.result = {
            "type": "single",
            "job_name": req.job_name,
            "scf_energy": final.result.scf_energy,
            "zpe_energy": final.result.zpe_energy,
            "g_corr": g_corr_hartree,
            "out_dir": out_dir,
            "log_files": {"final": final.log_path},
            "xyz_files": xyz_files,
            "cdxml_files": cdxml_files,
        }
        job.status = "done"
    _set_stage(job, "Done")


def _optimize_reactant_or_product(
    job: JobState, req: JobRequest, structure: Structure, label: str, level: LevelChoice,
    resume_from: Optional[JobState],
) -> Optional[OptOutcome]:
    """Runs the same PM6-then-final optimization pipeline
    _run_mechanism_job() uses, but for a single side and used by
    _run_ts_search_job() -- so a directly-supplied reactant/product
    also ends up with a real converged geometry/energy of its own,
    not just a raw reference used for TS-guess interpolation, per your
    request to be able to extract reaction-coordinate data from a
    TS-search-style study too. Returns None (never raises) on any
    failure -- this is a best-effort addition on top of the TS search,
    never something that should block it."""
    out_dir = os.path.join("jobs", req.job_name)
    pre = None
    if req.skip_pm6_preopt:
        _log(job, f"Skipping PM6 pre-optimization for the {label} -- starting directly at the final level.")
    else:
        pre = _checkpoint(resume_from, f"{label}_preopt")
        if pre is not None:
            _log(job, f"Reusing already-completed {label} pre-optimization from before the restart ({pre.log_path}).")
        else:
            _set_stage(job, f"Pre-optimizing {label} (PM6)")
            pre = preopt_with_escalation(
                structure, out_dir, label, pubchem_query=(req.reactant_id if label == "reactant" else req.product_id),
                nprocs=req.nprocs, mem_gb=req.mem_gb, max_repair_attempts=req.max_repair_attempts,
                cancel_event=job.cancel_event, executor=req.executor,
            )
            for line in pre.log:
                _log(job, line)
        job.checkpoints[f"{label}_preopt"] = pre
        if not pre.success:
            _log(job, f"{label.capitalize()} didn't reach a clean PM6 minimum -- skipping the reaction-coordinate "
                       "report for this side (the TS search itself is unaffected).")
            return None

    final = _checkpoint(resume_from, f"{label}_final")
    if final is not None:
        _log(job, f"Reusing already-completed {label} final-level optimization from before the restart ({final.log_path}).")
    else:
        _set_stage(job, f"Reoptimizing {label} at the final level")
        final = repair_and_optimize(
            pre.structure if pre else structure, os.path.join(out_dir, f"{label}_final.com"),
            method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
            nprocs=req.nprocs, mem_gb=req.mem_gb, max_repair_attempts=req.max_repair_attempts,
            cancel_event=job.cancel_event, executor=req.executor,
        )
        for line in final.log:
            _log(job, line)
    job.checkpoints[f"{label}_final"] = final
    if not final.success:
        _log(job, f"{label.capitalize()} didn't reach a clean minimum at the final level -- skipping the "
                   "reaction-coordinate report for this side (the TS search itself is unaffected).")
        return None
    return final


def _run_ts_search_job(
    job: JobState, req: JobRequest,
    reactant: Optional[Structure], product: Optional[Structure], ts_guess: Optional[Structure],
    resume_from: Optional[JobState] = None,
) -> None:
    out_dir = os.path.join("jobs", req.job_name)
    _maybe_write_pbs_template(out_dir, req)
    elements = sorted(
        {a[0] for a in (reactant.atoms if reactant else [])}
        | {a[0] for a in (product.atoms if product else [])}
        | {a[0] for a in (ts_guess.atoms if ts_guess else [])}
    )
    level = _level_choice(req.level, elements)

    # Independently optimize whichever of reactant/product was actually
    # given -- see _optimize_reactant_or_product(). The TS search itself
    # below is untouched, it still uses the original raw reactant/product
    # for guess interpolation and overlap validation, same as before.
    r_final = _optimize_reactant_or_product(job, req, reactant, "reactant", level, resume_from) if reactant else None
    p_final = _optimize_reactant_or_product(job, req, product, "product", level, resume_from) if product else None

    ts_outcome = _checkpoint(resume_from, "ts_search")
    if ts_outcome is not None:
        _log(job, f"Reusing already-completed TS search from before the restart ({ts_outcome.log_path}).")
    else:
        _set_stage(job, "Searching for the TS")
        ts_outcome = search_ts_staged(
            out_dir, method=level.method, basis=level.basis,
            basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
            nprocs=req.nprocs, mem_gb=req.mem_gb,
            reactant=reactant, product=product, ts_guess=ts_guess,
            skip_pm6=req.skip_pm6_preopt, use_calcall_fallback=req.ts_use_calcall, cancel_event=job.cancel_event,
            executor=req.executor,
        )
        for line in ts_outcome.log:
            _log(job, line)
    job.checkpoints["ts_search"] = ts_outcome
    if not ts_outcome.success:
        _fail(job, "Couldn't find a validated TS. See the log for what was tried.")
        return

    xyz_files = {}
    if req.export_xyz:
        xyz_path = convert_log_to_xyz(ts_outcome.log_path)
        if xyz_path:
            xyz_files["ts"] = xyz_path
        else:
            _log(job, f"Couldn't convert {ts_outcome.log_path} to .xyz (is obabel installed?).")

    cdxml_files = {}
    if req.export_cdxml:
        cdxml_path = convert_log_to_cdxml(ts_outcome.log_path)
        if cdxml_path:
            cdxml_files["ts"] = cdxml_path
        else:
            _log(job, f"Couldn't convert {ts_outcome.log_path} to .cdxml (is obabel installed?).")

    recovery_result = None
    reactant_candidates: Dict[str, Optional[OptOutcome]] = {"guess": r_final}
    product_candidates: Dict[str, Optional[OptOutcome]] = {"guess": p_final}
    if req.ts_recover_endpoints:
        _set_stage(job, "Recovering reactant/product from the TS")
        recovery = recover_endpoints_from_ts(
            ts_outcome, out_dir, level, nprocs=req.nprocs, mem_gb=req.mem_gb,
            reactant_ref=reactant, product_ref=product, cancel_event=job.cancel_event, executor=req.executor,
        )
        for line in recovery.log:
            _log(job, line)
        if isinstance(recovery, DistortionOutcome):
            recovery_result = {
                "classified": True,
                "reactant_log": recovery.reactant_candidate.outcome.log_path if recovery.reactant_candidate else None,
                "product_log": recovery.product_candidate.outcome.log_path if recovery.product_candidate else None,
            }
            reactant_candidates["ts_distortion"] = recovery.reactant_candidate.outcome if recovery.reactant_candidate else None
            product_candidates["ts_distortion"] = recovery.product_candidate.outcome if recovery.product_candidate else None
        else:
            recovery_result = {
                "classified": False,
                "side_a_log": recovery.side_a.log_path if recovery.side_a else None,
                "side_b_log": recovery.side_b.log_path if recovery.side_b else None,
            }

    # A reaction-coordinate report (energies/G-corr/xyz/cdxml for
    # reactant+TS+product, same shape as the mechanism study's) is only
    # built once there's a usable reactant AND product -- either given
    # directly and optimized above, or recovered from the TS just now.
    reaction_coordinate: dict = {}
    try:
        reactant_source, r_winner = select_best_candidate(reactant_candidates)
        product_source, p_winner = select_best_candidate(product_candidates)
    except ValueError:
        r_winner = p_winner = None
        if reactant is not None or product is not None or req.ts_recover_endpoints:
            _log(job, "Couldn't get a usable reactant and product for a reaction-coordinate report "
                       "(need both sides converged, either given directly or recovered from the TS).")

    log_files = {"ts": ts_outcome.log_path}
    if r_winner is not None and p_winner is not None:
        _log(job, format_candidate_comparison("reactant", reactant_candidates, reactant_source))
        _log(job, format_candidate_comparison("product", product_candidates, product_source))
        _set_stage(job, "Computing energetics")
        units = req.energy_units or ["kcal/mol"]
        energies = {}
        for unit in units:
            report = build_energy_report(r_winner.result, ts_outcome.result, p_winner.result, unit)
            _log(job, format_energy_report(report))
            energies[unit] = {"ts": report.ts, "product": report.product}

        log_files["reactant"] = r_winner.log_path
        log_files["product"] = p_winner.log_path

        if req.export_xyz:
            for label, log_path in (("reactant", r_winner.log_path), ("product", p_winner.log_path)):
                xyz_path = convert_log_to_xyz(log_path)
                if xyz_path:
                    xyz_files[label] = xyz_path
                else:
                    _log(job, f"Couldn't convert {log_path} to .xyz (is obabel installed?).")

        if req.export_cdxml:
            for label, log_path in (("reactant", r_winner.log_path), ("product", p_winner.log_path)):
                cdxml_path = convert_log_to_cdxml(log_path)
                if cdxml_path:
                    cdxml_files[label] = cdxml_path
                else:
                    _log(job, f"Couldn't convert {log_path} to .cdxml (is obabel installed?).")

        g_corr = {}
        if req.compute_gcorr:
            _set_stage(job, "Computing G-corr")
            temp_k = to_kelvin(req.gcorr_temperature, req.gcorr_temperature_unit)
            for label, log_path in (("reactant", r_winner.log_path), ("ts", ts_outcome.log_path), ("product", p_winner.log_path)):
                g_corr_hartree = run_freqchk(log_path.replace(".log", ".chk"), temp_k)
                if g_corr_hartree is None:
                    _log(job, f"Couldn't compute G-corr from {log_path.replace('.log', '.chk')} (is freqchk installed/working?).")
                    continue
                g_corr[label] = {unit: g_corr_hartree * HARTREE_TO[unit] for unit in units}

        reaction_coordinate = {
            "energy_units": units,
            "energies": energies,
            "reactant_source": reactant_source,
            "product_source": product_source,
            "g_corr": g_corr,
        }

    with _JOBS_LOCK:
        job.result = {
            "type": "ts_search",
            "job_name": req.job_name,
            "imaginary_freq": ts_outcome.result.imaginary_freqs[0],
            "match_overlap": ts_outcome.match_overlap,
            "log_files": log_files,
            "xyz_files": xyz_files,
            "cdxml_files": cdxml_files,
            "recovery": recovery_result,
            "out_dir": out_dir,
            **reaction_coordinate,
        }
        job.status = "done"
    _set_stage(job, "Done")


def _run_mechanism_job(
    job: JobState, req: JobRequest, reactant: Structure, product: Structure, ts_guess: Optional[Structure],
    resume_from: Optional[JobState] = None,
) -> None:
    check_reaction_match(reactant, product)
    r_elems = [a[0] for a in reactant.atoms]
    p_elems = [a[0] for a in product.atoms]
    if r_elems != p_elems:
        _fail(job, "Reactant/product atom count or order didn't match -- a TS search needs the same atoms in the same order on both ends.")
        return

    out_dir = os.path.join("jobs", req.job_name)
    _maybe_write_pbs_template(out_dir, req)
    recover = req.recover_missing_endpoint

    r_pre = None
    if req.skip_pm6_preopt:
        _log(job, "Skipping PM6 pre-optimization for the reactant -- starting directly at the final level.")
    else:
        r_pre = _checkpoint(resume_from, "reactant_preopt")
        if r_pre is not None:
            _log(job, f"Reusing already-completed reactant pre-optimization from before the restart ({r_pre.log_path}).")
        else:
            _set_stage(job, "Pre-optimizing reactant (PM6)")
            r_pre = preopt_with_escalation(
                reactant, out_dir, "reactant", pubchem_query=req.reactant_id, nprocs=req.nprocs, mem_gb=req.mem_gb,
                max_repair_attempts=req.max_repair_attempts, cancel_event=job.cancel_event, executor=req.executor,
            )
            for line in r_pre.log:
                _log(job, line)
        job.checkpoints["reactant_preopt"] = r_pre
        if not r_pre.success:
            if not recover:
                _fail(job, "Couldn't get the reactant to a clean PM6 minimum after every fallback.")
                return
            _log(job, "Reactant didn't reach a clean PM6 minimum after every fallback -- recovery is enabled, "
                       "skipping the final-level reoptimization for it and continuing with the product. TS search "
                       "will use the reactant's best (non-converged) attempt as a rough starting point; enabling "
                       "TS-mode distortion and/or IRC-endpoint reopt may recover a proper reactant afterward.")

    p_pre = None
    if req.skip_pm6_preopt:
        _log(job, "Skipping PM6 pre-optimization for the product -- starting directly at the final level.")
    else:
        p_pre = _checkpoint(resume_from, "product_preopt")
        if p_pre is not None:
            _log(job, f"Reusing already-completed product pre-optimization from before the restart ({p_pre.log_path}).")
        else:
            _set_stage(job, "Pre-optimizing product (PM6)")
            p_pre = preopt_with_escalation(
                product, out_dir, "product", pubchem_query=req.product_id, nprocs=req.nprocs, mem_gb=req.mem_gb,
                max_repair_attempts=req.max_repair_attempts, cancel_event=job.cancel_event, executor=req.executor,
            )
            for line in p_pre.log:
                _log(job, line)
        job.checkpoints["product_preopt"] = p_pre
        if not p_pre.success:
            if not recover:
                _fail(job, "Couldn't get the product to a clean PM6 minimum after every fallback.")
                return
            _log(job, "Product didn't reach a clean PM6 minimum after every fallback -- recovery is enabled, "
                       "skipping the final-level reoptimization for it and continuing with the reactant.")

    # A side is ready for the final-level stage if PM6 pre-opt was
    # skipped entirely (req.skip_pm6_preopt -- r_pre/p_pre stay None,
    # same sentinel endpoint recovery already uses for "nothing to
    # start final-opt from yet") or actually succeeded.
    r_preopt_ok = r_pre is None or r_pre.success
    p_preopt_ok = p_pre is None or p_pre.success
    if not r_preopt_ok and not p_preopt_ok:
        _fail(job, "Neither the reactant nor the product reached a clean PM6 minimum -- nothing usable to search a TS from.")
        return

    elements = sorted(
        {a[0] for a in (r_pre.structure if r_pre else reactant).atoms}
        | {a[0] for a in (p_pre.structure if p_pre else product).atoms}
    )
    level = _level_choice(req.level, elements)

    r_final = None
    if r_preopt_ok:
        r_final = _checkpoint(resume_from, "reactant_final")
        if r_final is not None:
            _log(job, f"Reusing already-completed reactant final-level optimization from before the restart ({r_final.log_path}).")
        else:
            _set_stage(job, "Reoptimizing reactant at the final level")
            r_final = repair_and_optimize(
                r_pre.structure if r_pre else reactant, os.path.join(out_dir, "reactant_final.com"),
                method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
                nprocs=req.nprocs, mem_gb=req.mem_gb, max_repair_attempts=req.max_repair_attempts,
                cancel_event=job.cancel_event, executor=req.executor,
            )
            for line in r_final.log:
                _log(job, line)
        job.checkpoints["reactant_final"] = r_final
        if not r_final.success:
            if not recover:
                _fail(job, "Reactant didn't reach a clean minimum at the final level.")
                return
            _log(job, "Reactant didn't reach a clean minimum at the final level either -- recovery is enabled, "
                       "continuing without a verified reactant.")

    p_final = None
    if p_preopt_ok:
        p_final = _checkpoint(resume_from, "product_final")
        if p_final is not None:
            _log(job, f"Reusing already-completed product final-level optimization from before the restart ({p_final.log_path}).")
        else:
            _set_stage(job, "Reoptimizing product at the final level")
            p_final = repair_and_optimize(
                p_pre.structure if p_pre else product, os.path.join(out_dir, "product_final.com"),
                method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
                nprocs=req.nprocs, mem_gb=req.mem_gb, max_repair_attempts=req.max_repair_attempts,
                cancel_event=job.cancel_event, executor=req.executor,
            )
            for line in p_final.log:
                _log(job, line)
        job.checkpoints["product_final"] = p_final
        if not p_final.success:
            if not recover:
                _fail(job, "Product didn't reach a clean minimum at the final level.")
                return
            _log(job, "Product didn't reach a clean minimum at the final level either -- recovery is enabled, "
                       "continuing without a verified product.")

    reactant_missing = r_final is None or not r_final.success
    product_missing = p_final is None or not p_final.success
    if reactant_missing and product_missing:
        _fail(job, "Neither the reactant nor the product reached a clean minimum at the final level -- nothing usable to search a TS from.")
        return

    reactant_for_ts = (r_final.structure if r_final else None) or (r_pre.structure if r_pre else None) or reactant
    product_for_ts = (p_final.structure if p_final else None) or (p_pre.structure if p_pre else None) or product
    if reactant_missing or product_missing:
        _log(job, f"Using {'the reactant' if reactant_missing else 'the product'}'s best available (non-converged) "
                   "geometry as the TS-search reference for that side -- not independently verified as a minimum.")

    ts_outcome = _checkpoint(resume_from, "ts")
    if ts_outcome is not None:
        _log(job, f"Reusing already-completed TS search from before the restart ({ts_outcome.log_path}).")
    else:
        _set_stage(job, "Searching for the TS")
        ts_com = os.path.join(out_dir, "ts.com")
        ts_outcome = run_ts_search(
            reactant_for_ts, product_for_ts, ts_com,
            method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
            nprocs=req.nprocs, mem_gb=req.mem_gb, ts_guess=ts_guess, cancel_event=job.cancel_event,
            executor=req.executor,
        )
        for line in ts_outcome.log:
            _log(job, line)
    job.checkpoints["ts"] = ts_outcome
    if not ts_outcome.success:
        _fail(job, "Couldn't find a TS whose imaginary mode matches the reaction.")
        return

    reactant_trusted, product_trusted = not reactant_missing, not product_missing

    distortion_outcome = None
    if req.run_ts_distortion:
        _set_stage(job, "Cross-checking via TS-mode distortion")
        distortion_outcome = run_ts_distortion_check(
            ts_outcome, reactant_for_ts, product_for_ts, out_dir, level,
            nprocs=req.nprocs, mem_gb=req.mem_gb, factor=req.distortion_factor, cancel_event=job.cancel_event,
            reactant_trusted=reactant_trusted, product_trusted=product_trusted, executor=req.executor,
        )
        for line in distortion_outcome.log:
            _log(job, line)

    irc_outcome = None
    if req.run_irc:
        _set_stage(job, f"Running the IRC ({req.irc_points} points/direction)")
        irc_com = os.path.join(out_dir, "irc.com")
        irc_outcome = run_irc(
            ts_outcome.structure, reactant_for_ts, product_for_ts, irc_com,
            method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
            nprocs=req.nprocs, mem_gb=req.mem_gb, maxpoints=req.irc_points, cancel_event=job.cancel_event,
            reactant_trusted=reactant_trusted, product_trusted=product_trusted, executor=req.executor,
        )
        for line in irc_outcome.log:
            _log(job, line)
    else:
        _log(job, "Skipping IRC verification.")

    irc_reopt_outcome = None
    if (
        req.run_irc_endpoint_reopt
        and irc_outcome is not None
        and irc_outcome.forward_structure is not None
        and irc_outcome.reverse_structure is not None
    ):
        _set_stage(job, "Reoptimizing the IRC endpoints")
        irc_reopt_outcome = run_irc_endpoint_reopt(
            irc_outcome, reactant_for_ts, product_for_ts, out_dir, level,
            nprocs=req.nprocs, mem_gb=req.mem_gb, cancel_event=job.cancel_event,
            reactant_trusted=reactant_trusted, product_trusted=product_trusted, executor=req.executor,
        )
        for line in irc_reopt_outcome.log:
            _log(job, line)

    reactant_candidates = {
        "guess": r_final if (r_final is not None and r_final.success) else None,
        "ts_distortion": distortion_outcome.reactant_candidate.outcome if distortion_outcome and distortion_outcome.reactant_candidate else None,
        "irc_reopt": irc_reopt_outcome.reactant_candidate.outcome if irc_reopt_outcome and irc_reopt_outcome.reactant_candidate else None,
    }
    product_candidates = {
        "guess": p_final if (p_final is not None and p_final.success) else None,
        "ts_distortion": distortion_outcome.product_candidate.outcome if distortion_outcome and distortion_outcome.product_candidate else None,
        "irc_reopt": irc_reopt_outcome.product_candidate.outcome if irc_reopt_outcome and irc_reopt_outcome.product_candidate else None,
    }

    unrecovered = []
    reactant_source = product_source = None
    r_winner = p_winner = None
    try:
        reactant_source, r_winner = select_best_candidate(reactant_candidates)
    except ValueError:
        unrecovered.append("reactant")
    try:
        product_source, p_winner = select_best_candidate(product_candidates)
    except ValueError:
        unrecovered.append("product")
    if unrecovered:
        _fail(job, (
            f"Couldn't get a usable {' or '.join(unrecovered)} for the energy report even with recovery "
            f"enabled -- the TS was found (imaginary freq {ts_outcome.result.imaginary_freqs[0]:.1f} cm^-1) "
            f"but neither the original optimization nor TS-mode distortion/IRC-endpoint reopt produced a "
            f"usable {' or '.join(unrecovered)}. Try enabling TS-mode distortion and/or IRC-endpoint reopt, "
            "or supply a better starting geometry."
        ))
        return

    _log(job, format_candidate_comparison("reactant", reactant_candidates, reactant_source))
    _log(job, format_candidate_comparison("product", product_candidates, product_source))
    reactant_result, product_result = r_winner.result, p_winner.result
    reactant_log_path, product_log_path = r_winner.log_path, p_winner.log_path

    _set_stage(job, "Computing energetics")
    units = req.energy_units or ["kcal/mol"]
    energies = {}
    for unit in units:
        report = build_energy_report(reactant_result, ts_outcome.result, product_result, unit)
        _log(job, format_energy_report(report))
        energies[unit] = {"ts": report.ts, "product": report.product}

    log_files = {
        "reactant": reactant_log_path,
        "product": product_log_path,
        "ts": ts_outcome.log_path,
    }
    if irc_outcome is not None:
        log_files["irc"] = irc_outcome.log_path

    xyz_files = {}
    if req.export_xyz:
        for label, log_path in (("reactant", reactant_log_path), ("product", product_log_path), ("ts", ts_outcome.log_path)):
            xyz_path = convert_log_to_xyz(log_path)
            if xyz_path:
                xyz_files[label] = xyz_path
            else:
                _log(job, f"Couldn't convert {log_path} to .xyz (is obabel installed?).")

    cdxml_files = {}
    if req.export_cdxml:
        for label, log_path in (("reactant", reactant_log_path), ("product", product_log_path), ("ts", ts_outcome.log_path)):
            cdxml_path = convert_log_to_cdxml(log_path)
            if cdxml_path:
                cdxml_files[label] = cdxml_path
            else:
                _log(job, f"Couldn't convert {log_path} to .cdxml (is obabel installed?).")

    g_corr = {}
    if req.compute_gcorr:
        _set_stage(job, "Computing G-corr")
        temp_k = to_kelvin(req.gcorr_temperature, req.gcorr_temperature_unit)
        for label, log_path in (("reactant", reactant_log_path), ("ts", ts_outcome.log_path), ("product", product_log_path)):
            g_corr_hartree = run_freqchk(log_path.replace(".log", ".chk"), temp_k)
            if g_corr_hartree is None:
                _log(job, f"Couldn't compute G-corr from {log_path.replace('.log', '.chk')} (is freqchk installed/working?).")
                continue
            g_corr[label] = {unit: g_corr_hartree * HARTREE_TO[unit] for unit in units}

    with _JOBS_LOCK:
        job.result = {
            "type": "mechanism",
            "job_name": req.job_name,
            "energy_units": units,
            "energies": energies,
            "ts_imaginary_freq": ts_outcome.result.imaginary_freqs[0],
            "ts_match_overlap": ts_outcome.match_overlap,
            "irc_ran": req.run_irc,
            "irc_connected": irc_outcome.success if irc_outcome else None,
            "reactant_source": reactant_source,
            "product_source": product_source,
            "reactant_recovered": reactant_missing,
            "product_recovered": product_missing,
            "log_files": log_files,
            "xyz_files": xyz_files,
            "cdxml_files": cdxml_files,
            "g_corr": g_corr,
            "gcorr_temperature_k": to_kelvin(req.gcorr_temperature, req.gcorr_temperature_unit) if req.compute_gcorr else None,
            "out_dir": out_dir,
        }
        job.status = "done"
    _set_stage(job, "Done")


def _run_job_thread(job: JobState, req: JobRequest, resume_from: Optional[JobState] = None) -> None:
    job.req = req
    try:
        if req.study_type == "single":
            upload = _UPLOADS.get(req.reactant_token or "")
            if upload is None:
                _fail(job, "Structure file not found -- try uploading it again.")
                return
            _run_single_job(job, req, upload.structure, resume_from=resume_from)
        elif req.study_type == "ts_search":
            r_upload = _UPLOADS.get(req.reactant_token or "")
            p_upload = _UPLOADS.get(req.product_token or "")
            ts_upload = _UPLOADS.get(req.ts_guess_token or "")
            reactant = r_upload.structure if r_upload else None
            product = p_upload.structure if p_upload else None
            ts_guess = ts_upload.structure if ts_upload else None
            if ts_guess is None and (reactant is None or product is None):
                _fail(job, "Need a TS guess, or both reactant and product -- try uploading again.")
                return
            _run_ts_search_job(job, req, reactant, product, ts_guess, resume_from=resume_from)
        else:
            r_upload = _UPLOADS.get(req.reactant_token or "")
            p_upload = _UPLOADS.get(req.product_token or "")
            if r_upload is None or p_upload is None:
                _fail(job, "Reactant/product file not found -- try uploading them again.")
                return
            ts_guess = None
            if req.ts_guess_token:
                ts_upload = _UPLOADS.get(req.ts_guess_token)
                ts_guess = ts_upload.structure if ts_upload else None
            _run_mechanism_job(job, req, r_upload.structure, p_upload.structure, ts_guess, resume_from=resume_from)
    except GaussianCancelledError:
        with _JOBS_LOCK:
            job.status = "cancelled"
            job.log.append("Cancelled -- Gaussian process terminated. Whatever finished before this point is kept below.")
    except Exception as e:  # noqa: BLE001 -- surface any unexpected failure to the UI instead of a dead job
        _fail(job, f"Unexpected error: {e}")


# --------------------------------------------------------------- endpoints

@app.get("/")
def index() -> FileResponse:
    return FileResponse(str(WEB_DIR / "index.html"))


@app.post("/api/analyze")
async def analyze(file: UploadFile):
    token = uuid.uuid4().hex
    suffix = Path(file.filename or "structure.xyz").suffix or ".xyz"
    dest = UPLOAD_DIR / f"{token}{suffix}"
    dest.write_bytes(await file.read())

    try:
        structure = load_structure(str(dest), label=Path(file.filename or "structure").stem)
    except Exception as e:
        raise HTTPException(400, f"Couldn't parse {file.filename}: {e}")

    _UPLOADS[token] = _Upload(structure=structure, path=str(dest))
    elements = sorted({a[0] for a in structure.atoms})
    return {"token": token, "elements": elements, "atom_count": len(structure.atoms)}


class SmilesRequest(BaseModel):
    smiles: str
    label: str = "structure"


@app.post("/api/analyze_smiles")
def analyze_smiles(req: SmilesRequest):
    """Same shape/purpose as /api/analyze, for a SMILES string instead
    of an uploaded file -- the GUI's equivalent of the CLI's SMILES
    intake option (embeds with RDKit via geometry.from_smiles)."""
    smiles = req.smiles.strip()
    if not smiles:
        raise HTTPException(400, "Enter a SMILES string first.")
    try:
        structure = from_smiles(smiles, label=req.label or "structure")
    except Exception as e:
        raise HTTPException(400, f"Couldn't embed SMILES {smiles!r}: {e}")

    token = uuid.uuid4().hex
    dest = UPLOAD_DIR / f"{token}.xyz"
    structure.to_xyz_file(str(dest))
    _UPLOADS[token] = _Upload(structure=structure, path=str(dest))
    elements = sorted({a[0] for a in structure.atoms})
    return {"token": token, "elements": elements, "atom_count": len(structure.atoms)}


def _register_upload(structure: Structure) -> str:
    token = uuid.uuid4().hex
    dest = UPLOAD_DIR / f"{token}.xyz"
    structure.to_xyz_file(str(dest))
    _UPLOADS[token] = _Upload(structure=structure, path=str(dest))
    return token


class RefineRequest(BaseModel):
    token: str
    description: str


@app.get("/api/refine/available")
def refine_available():
    return {"available": llm_assist.is_available()}


@app.post("/api/refine")
def refine(req: RefineRequest):
    """Ask Claude to read a free-text description alongside an already
    loaded structure and suggest a better starting geometry -- see
    llm_assist.py for what this is (and isn't) allowed to do. Any
    alternate structure offered back has already had its formula
    checked against the original before being registered as a new
    upload token."""
    upload = _UPLOADS.get(req.token)
    if upload is None:
        raise HTTPException(400, "Unknown structure -- upload or load it again.")
    if not req.description.strip():
        raise HTTPException(400, "Enter a description first.")

    original_formula = upload.structure.formula()
    try:
        hint = llm_assist.refine_structure(
            description=req.description,
            smiles=None,
            formula=original_formula,
            label=upload.structure.label,
        )
    except llm_assist.LLMUnavailableError as e:
        raise HTTPException(503, str(e))

    result = {
        "matches_description": hint.matches_description,
        "notes": hint.notes,
        "suggested_charge": hint.suggested_charge,
        "suggested_multiplicity": hint.suggested_multiplicity,
        "pubchem_token": None,
        "pubchem_elements": None,
        "pubchem_atom_count": None,
        "smiles_token": None,
        "suggested_smiles": None,
        "smiles_elements": None,
        "smiles_atom_count": None,
    }

    if hint.suggested_pubchem_query:
        try:
            pubchem_structure = from_pubchem(hint.suggested_pubchem_query, label=upload.structure.label)
            pubchem_formula = pubchem_structure.formula()
            if pubchem_formula == original_formula:
                result["pubchem_token"] = _register_upload(pubchem_structure)
                result["pubchem_elements"] = sorted({a[0] for a in pubchem_structure.atoms})
                result["pubchem_atom_count"] = len(pubchem_structure.atoms)
                result["notes"] += (
                    f" (PubChem match for '{hint.suggested_pubchem_query}', "
                    f"same formula {original_formula} -- available to use.)"
                )
            else:
                result["notes"] += (
                    f" (Tried PubChem for '{hint.suggested_pubchem_query}' but its formula "
                    f"{pubchem_formula} doesn't match yours ({original_formula}) -- not offering it.)"
                )
        except Exception as e:
            result["notes"] += f" (PubChem lookup for '{hint.suggested_pubchem_query}' failed: {e})"

    if hint.suggested_smiles:
        try:
            smiles_structure = from_smiles(hint.suggested_smiles, label=upload.structure.label)
            result["smiles_token"] = _register_upload(smiles_structure)
            result["suggested_smiles"] = hint.suggested_smiles
            result["smiles_elements"] = sorted({a[0] for a in smiles_structure.atoms})
            result["smiles_atom_count"] = len(smiles_structure.atoms)
        except Exception as e:
            result["notes"] += f" (Couldn't embed the suggested SMILES {hint.suggested_smiles!r}: {e})"

    return result


@app.get("/api/units")
def units():
    return {"units": list(HARTREE_TO)}


@app.get("/api/pbs_template")
def pbs_template():
    """The starting text for the GUI's editable PBS script box --
    job_runner.default_pbs_template()'s output, with __JOB__ standing
    in for the job name/`.com` filename (substituted per-job at
    submission time, see job_runner.render_pbs_script())."""
    from .job_runner import default_pbs_template

    return {"template": default_pbs_template()}


@app.post("/api/jobs")
def create_job(req: JobRequest):
    if req.study_type not in ("single", "mechanism", "ts_search"):
        raise HTTPException(400, "study_type must be 'single', 'mechanism', or 'ts_search'")

    if req.study_type == "single":
        if not req.reactant_token or req.reactant_token not in _UPLOADS:
            raise HTTPException(400, "Unknown reactant/structure file -- upload it again.")
    elif req.study_type == "mechanism":
        if not req.reactant_token or req.reactant_token not in _UPLOADS:
            raise HTTPException(400, "Unknown reactant/structure file -- upload it again.")
        if not req.product_token or req.product_token not in _UPLOADS:
            raise HTTPException(400, "Unknown product file -- upload it again.")
        if not req.energy_units:
            raise HTTPException(400, "Pick at least one energy unit.")
    else:  # ts_search
        has_guess = bool(req.ts_guess_token) and req.ts_guess_token in _UPLOADS
        has_pair = (
            bool(req.reactant_token) and req.reactant_token in _UPLOADS
            and bool(req.product_token) and req.product_token in _UPLOADS
        )
        if not has_guess and not has_pair:
            raise HTTPException(400, "Provide a TS guess, or both a reactant and a product.")

    job_id = uuid.uuid4().hex
    job = JobState(id=job_id, req=req)
    with _JOBS_LOCK:
        _JOBS[job_id] = job

    thread = threading.Thread(target=_run_job_thread, args=(job, req), daemon=True)
    thread.start()
    return {"job_id": job_id}


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "Unknown job id")
        return {
            "status": job.status,
            "stage": job.stage,
            "log": job.log,
            "result": job.result,
            "error": job.error,
        }


@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str):
    """Terminate whatever Gaussian process this job currently has
    running (see local_runner.run_local's cancel_event polling) and
    stop it from starting another. Whatever finished before this point
    stays on disk and in the job's log -- nothing is deleted."""
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "Unknown job id")
        if job.status != "running":
            raise HTTPException(400, f"Job is {job.status}, not running -- nothing to cancel.")
    job.cancel_event.set()
    return {"status": "cancelling"}


@app.post("/api/jobs/{job_id}/restart")
def restart_job(job_id: str):
    """Re-run a cancelled job's original request, skipping straight
    past whatever stage(s) already finished cleanly before it was
    cancelled (see JobState.checkpoints / _checkpoint()) -- resumes
    from the last fully completed stage, not from scratch and not
    mid-optimization-cycle. Returns a new job id; the old one's log
    and checkpoints are left as they were for reference."""
    with _JOBS_LOCK:
        old_job = _JOBS.get(job_id)
        if old_job is None:
            raise HTTPException(404, "Unknown job id")
        if old_job.status != "cancelled":
            raise HTTPException(400, f"Job is {old_job.status}, not cancelled -- nothing to restart.")
        if old_job.req is None:
            raise HTTPException(400, "No original request stored for this job -- can't restart it.")

        new_job_id = uuid.uuid4().hex
        new_job = JobState(id=new_job_id)
        _JOBS[new_job_id] = new_job

    thread = threading.Thread(target=_run_job_thread, args=(new_job, old_job.req, old_job), daemon=True)
    thread.start()
    return {"job_id": new_job_id}


@app.delete("/api/jobs/{job_id}")
def delete_job(job_id: str):
    """Remove a study's entry from the Studies list -- just the
    in-memory JobState the GUI tracks (log, stage, result), never the
    .com/.log files a study actually produced under jobs/<name>/,
    which are untouched. Refuses to delete a job that's still running
    -- stop it first, so there's no orphaned Gaussian process left
    with no card to track or cancel it."""
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "Unknown job id")
        if job.status == "running":
            raise HTTPException(400, "Stop the study before removing it.")
        del _JOBS[job_id]
    return {"status": "deleted"}


def main() -> None:
    port = int(os.environ.get("GAUSSBOT_GUI_PORT", "8765"))
    url = f"http://127.0.0.1:{port}"
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    print(f"gaussbot GUI running at {url} (Ctrl+C to stop)")
    uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning")


if __name__ == "__main__":
    main()
