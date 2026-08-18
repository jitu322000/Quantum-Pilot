"""
webapp.py

A local web GUI wrapped around the RHF/CIS pipeline, same shape as
gaussbot's own webapp.py: accepts a job over HTTP, runs it on a
background thread (so the request returns immediately), and lets the
frontend poll for progress. Every actual computation goes through the
same functions cli.py calls (intake.py, rhf.py, cis.py) unchanged.

Run: gamessbot-gui   (opens a browser tab at http://127.0.0.1:8767)
"""

from __future__ import annotations

import os
import threading
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

from gaussbot.geometry import Structure, from_file, from_smiles
from gaussbot.level_select import BASIS_SETS as GAUSSIAN_BASIS_SETS
from gaussbot.level_select import METHODS as GAUSSIAN_METHODS

from .active_space import ActiveSpaceSuggestion, build_active_space
from .casscf import (
    CASSCFOutcome,
    mo_source_from_casscf_outcome,
    run_casscf_staged,
    run_casscf_with_smaller_active_space_recovery,
    suggest_active_space_from_cis,
)
from .cis import CISOutcome, run_cis
from .energetics import build_energy_table, format_latex_table
from .gamess_input import data_block_from_gamess_inp, data_block_from_gaussian_log
from .intake import optimize_and_convert
from .level_select import DEFAULT_CASSCF_MEM_MWORDS, DEFAULT_RUNGMS_PATH, DEFAULT_SCRATCH_DIR, GBASIS_CHOICES
from .local_runner import GamessCancelledError
from .orbital_character import OrbitalCharacter, parse_orbital_character
from .rhf import RHFOutcome, run_rhf_staged
from .transitn import TransitnOutcome, run_transitn
from .xmcqdpt import XMCQDPTOutcome, run_xmcqdpt

WEB_DIR = Path(__file__).parent / "web"
UPLOAD_DIR = Path(os.environ.get("GAMESSBOT_UPLOAD_DIR", "/tmp/gamessbot_uploads"))
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="gamessbot")
app.mount("/static", StaticFiles(directory=str(WEB_DIR)), name="static")


@app.middleware("http")
async def _no_cache_static(request, call_next):
    """Without this, browsers can serve /static/style.css or app.js
    from their own heuristic cache on a plain reload (no Cache-Control
    header means the browser is free to guess a freshness lifetime) --
    confirmed this bit you after a CSS fix that curl/the server showed
    as correct but a normal reload kept showing stale. Forcing
    revalidation on every request means a plain reload always picks up
    the latest file (StaticFiles's own ETag/Last-Modified still make
    that cheap -- a 304 when nothing changed)."""
    response = await call_next(request)
    if request.url.path.startswith("/static/"):
        response.headers["Cache-Control"] = "no-cache"
    return response


@dataclass
class _Upload:
    kind: str  # "gaussian_log" | "gamess_inp" | "guess"
    path: str
    structure: Optional[Structure] = None  # set for kind == "guess"


_UPLOADS: Dict[str, _Upload] = {}


@dataclass
class JobState:
    id: str
    status: str = "running"  # running | awaiting_active_space | done | error | cancelled
    stage: str = "starting"
    log: List[str] = field(default_factory=list)
    result: Optional[dict] = None
    error: Optional[str] = None
    cancel_event: threading.Event = field(default_factory=threading.Event)
    # Stashed here (not in `result`, which is the JSON-serializable view the
    # frontend polls) so the /casscf continuation endpoint below can pick up
    # exactly where the RHF/CIS thread left off, without rerunning either.
    rhf_outcome: Optional[RHFOutcome] = None
    active_space_suggestion: Optional[ActiveSpaceSuggestion] = None
    out_dir: Optional[str] = None
    rungms_path: Optional[str] = None
    scratch_dir: Optional[str] = None
    ncpus: int = 1
    mem_mwords: int = 1
    executor: str = "local"
    # The most recently *converged* CASSCF outcome across every combo run
    # so far (this batch or an earlier one, via "try another combination")
    #, lets a later combo start its orbital optimization from these
    # optimized orbitals instead of the original RHF/closed-shell ones,
    # per your request.
    last_casscf_outcome: Optional[CASSCFOutcome] = None
    # Strictly increasing across the whole job (every batch, including
    # ones added later via "try another combination") -- baked into
    # each combo's file prefix so two combos that happen to share the
    # same active-space/nstate never overwrite each other's .inp/.log/
    # .dat files. See _run_casscf_batch_thread/_run_one_combo.
    next_combo_index: int = 1


_JOBS: Dict[str, JobState] = {}
_JOBS_LOCK = threading.Lock()


def _log(job: JobState, message: str) -> None:
    with _JOBS_LOCK:
        job.log.append(message)


def _set_stage(job: JobState, stage: str) -> None:
    with _JOBS_LOCK:
        job.stage = stage
        job.log.append(f"--- {stage},-")


def _fail(job: JobState, message: str) -> None:
    with _JOBS_LOCK:
        job.status = "error"
        job.error = message
        job.log.append(message)


#,-------------------------------------------------------------- API models

class JobRequest(BaseModel):
    job_name: str
    geometry_source: str  # "gaussian_log" | "gamess_inp" | "guess"
    upload_token: Optional[str] = None
    nprocs: int = 2  # Gaussian resources, only used for geometry_source == "guess"
    mem_gb: int = 2
    gaussian_method: str = "B3LYP"  # final-level Gaussian optimization, only used for geometry_source == "guess"
    gaussian_basis: str = "6-31G(d)"
    charge: int = 0
    mult: int = 1
    gbasis_line: str = "GBASIS=STO NGAUSS=3"
    use_soscf: bool = True
    run_cis: bool = False
    nstate: int = 5
    run_casscf: bool = False  # needs run_cis, its excitations drive the active-space suggestion
    casscf_threshold: float = 0.20
    rungms_path: str = DEFAULT_RUNGMS_PATH
    scratch_dir: str = DEFAULT_SCRATCH_DIR
    ncpus: int = 1
    mem_mwords: int = 1
    executor: str = "local"  # "local" (run rungms directly) or "pbs" (submit via qsub, see job_runner.py)
    pbs_script_template: Optional[str] = None  # user-edited PBS script (see job_runner.render_pbs_script), only used when executor == "pbs"


def _transition_dict(t) -> dict:
    return {"from_mo": t.from_mo, "to_mo": t.to_mo, "coefficient": t.coefficient}


def _state_dict(s) -> dict:
    return {
        "index": s.index, "energy": s.energy, "spin": s.spin, "space_sym": s.space_sym,
        "transitions": [_transition_dict(t) for t in s.transitions],
    }


def _active_space_dict(s: ActiveSpaceSuggestion, orbital_character: Optional[Dict[int, "OrbitalCharacter"]] = None) -> dict:
    return {
        "nmcc": s.nmcc, "ndoc": s.ndoc, "nval": s.nval, "n_occ": s.n_occ, "norb": s.norb,
        "occ_selected": s.occ_selected, "virt_selected": s.virt_selected,
        "occ_dropped": s.occ_dropped, "virt_dropped": s.virt_dropped,
        "scores": {str(k): v for k, v in s.scores.items()}, "capped": s.capped,
        "orbital_character": (
            {str(mo): {"label": c.label, "coefficient": c.coefficient} for mo, c in orbital_character.items()}
            if orbital_character else {}
        ),
    }


def _casscf_result_dict(outcome: CASSCFOutcome) -> dict:
    result = outcome.result
    return {
        "success": outcome.success, "log_path": outcome.log_path,
        "converged": result.converged if result else False,
        "normal_termination": result.normal_termination if result else False,
        "final_energy": result.final_energy if result else None,
        "state_energies": result.state_energies if result else [],
        "active_space": _active_space_dict(outcome.active_space),
    }


def _xmcqdpt_result_dict(outcome: XMCQDPTOutcome) -> dict:
    result = outcome.result
    return {
        "success": outcome.success, "log_path": outcome.log_path,
        "normal_termination": result.normal_termination if result else False,
        "casscf_converged": result.casscf_converged if result else False,
        "mcqdpt_state_energies": result.mcqdpt_state_energies if result else [],
    }


def _transitn_result_dict(outcome: TransitnOutcome) -> dict:
    return {
        "success": outcome.success, "log_path": outcome.log_path,
        "oscillator_strengths": {str(k): v for k, v in outcome.oscillator_strengths.items()},
    }


def _energy_table_dict(
    casscf_outcome: CASSCFOutcome, xmcqdpt_outcome: Optional[XMCQDPTOutcome], nstate: int,
    oscillator_strengths: Optional[Dict[int, float]] = None,
) -> dict:
    xmcqdpt_result = xmcqdpt_outcome.result if (xmcqdpt_outcome and xmcqdpt_outcome.success) else None
    rows = build_energy_table(casscf_outcome.result, xmcqdpt_result, oscillator_strengths)
    latex = format_latex_table(rows, nstate=nstate, active_space=casscf_outcome.active_space)
    return {
        "rows": [
            {
                "state_index": r.state_index, "label": r.label, "casscf_ev": r.casscf_ev,
                "xmcqdpt_ev": r.xmcqdpt_ev, "oscillator_strength": r.oscillator_strength,
            }
            for r in rows
        ],
        "latex": latex,
    }


def _maybe_write_pbs_template(out_dir: str, req: JobRequest) -> None:
    """If PBS mode is on and the user edited a script template (GUI
    text box), write it once as <out_dir>/pbs_template.sh,
    job_runner.run_pbs() picks it up automatically from there for
    every stage in this job."""
    if req.executor == "pbs" and req.pbs_script_template:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "pbs_template.sh"), "w") as f:
            f.write(req.pbs_script_template)


def _run_job_thread(job: JobState, req: JobRequest) -> None:
    try:
        upload = _UPLOADS.get(req.upload_token or "")

        if req.geometry_source == "gaussian_log":
            if upload is None:
                _fail(job, "Gaussian log not found, try uploading it again.")
                return
            _set_stage(job, "Converting the Gaussian log to a GAMESS $DATA block")
            data_block = data_block_from_gaussian_log(upload.path, title=req.job_name)
        elif req.geometry_source == "gamess_inp":
            if upload is None:
                _fail(job, "GAMESS input/punch file not found, try uploading it again.")
                return
            _set_stage(job, "Extracting the $DATA block")
            data_block = data_block_from_gamess_inp(upload.path, title=req.job_name)
        else:  # guess
            if upload is None or upload.structure is None:
                _fail(job, "Guess geometry not found, try uploading/entering it again.")
                return
            out_dir = os.path.join("jobs", req.job_name)
            _set_stage(job, "Pre-optimizing the guess geometry with Gaussian (PM6)")
            try:
                data_block, pre, final = optimize_and_convert(
                    upload.structure, out_dir, req.job_name, method=req.gaussian_method, basis=req.gaussian_basis,
                    nprocs=req.nprocs, mem_gb=req.mem_gb, cancel_event=job.cancel_event,
                )
            except ValueError as e:
                _fail(job, str(e))
                return
            for line in pre.log:
                _log(job, line)
            for line in final.log:
                _log(job, line)

        out_dir = os.path.join("jobs", req.job_name)
        os.makedirs(out_dir, exist_ok=True)
        _maybe_write_pbs_template(out_dir, req)

        _set_stage(job, "Running RHF")
        rhf_outcome: RHFOutcome = run_rhf_staged(
            data_block, out_dir, charge=req.charge, mult=req.mult, gbasis_line=req.gbasis_line,
            rungms_path=req.rungms_path, scratch_dir=req.scratch_dir, ncpus=req.ncpus,
            mem_mwords=req.mem_mwords, use_soscf=req.use_soscf, cancel_event=job.cancel_event,
            executor=req.executor,
        )
        for line in rhf_outcome.trail:
            _log(job, line)
        if not rhf_outcome.success:
            _fail(job, f"RHF did not converge. See {rhf_outcome.log_path}.")
            return

        result = {
            "job_name": req.job_name,
            "rhf": {
                "energy": rhf_outcome.energy, "norb": rhf_outcome.norb, "log_path": rhf_outcome.log_path,
            },
            "cis": None,
            "out_dir": out_dir,
        }

        if req.run_cis:
            _set_stage(job, f"Running CIS ({req.nstate} states)")
            cis_outcome: CISOutcome = run_cis(
                rhf_outcome, out_dir, nstate=req.nstate, rungms_path=req.rungms_path,
                scratch_dir=req.scratch_dir, ncpus=req.ncpus, mem_mwords=req.mem_mwords,
                cancel_event=job.cancel_event, executor=req.executor,
            )
            if not cis_outcome.success:
                _fail(job, f"CIS did not terminate normally. See {cis_outcome.log_path}.")
                return
            result["cis"] = {
                "log_path": cis_outcome.log_path,
                "states": [_state_dict(s) for s in cis_outcome.states],
            }

            if req.run_casscf:
                _set_stage(job, "Suggesting a CASSCF active space from the CIS excitations")
                try:
                    suggestion = suggest_active_space_from_cis(
                        cis_outcome, rhf_outcome, threshold=req.casscf_threshold,
                    )
                except ValueError as e:
                    _fail(job, f"Couldn't suggest an active space: {e}")
                    return
                orbital_character = parse_orbital_character(
                    rhf_outcome.log_path, suggestion.occ_selected + suggestion.virt_selected,
                )
                result["active_space_suggestion"] = _active_space_dict(suggestion, orbital_character)
                with _JOBS_LOCK:
                    job.rhf_outcome = rhf_outcome
                    job.active_space_suggestion = suggestion
                    job.out_dir = out_dir
                    job.rungms_path = req.rungms_path
                    job.scratch_dir = req.scratch_dir
                    job.ncpus = req.ncpus
                    job.mem_mwords = req.mem_mwords
                    job.executor = req.executor
                    job.result = result
                    job.status = "awaiting_active_space"
                _log(job, "Suggested an active space, confirm it (or adjust it) to continue to CASSCF.")
                return

        with _JOBS_LOCK:
            job.result = result
            job.status = "done"
        _set_stage(job, "Done")

    except GamessCancelledError:
        with _JOBS_LOCK:
            job.status = "cancelled"
            job.log.append("Cancelled, the GAMESS process was terminated. Whatever finished before this point is kept above.")
    except Exception as e:  # noqa: BLE001, surface any unexpected failure to the UI instead of a dead job
        _fail(job, f"Unexpected error: {e}")


#,------------------------------------------------------------- endpoints

@app.get("/")
def index() -> FileResponse:
    return FileResponse(str(WEB_DIR / "index.html"))


@app.post("/api/upload")
async def upload(file: UploadFile, kind: str):
    if kind not in ("gaussian_log", "gamess_inp", "guess"):
        raise HTTPException(400, "kind must be 'gaussian_log', 'gamess_inp', or 'guess'")

    token = uuid.uuid4().hex
    suffix = Path(file.filename or "structure").suffix or ".txt"
    dest = UPLOAD_DIR / f"{token}{suffix}"
    dest.write_bytes(await file.read())

    structure = None
    if kind == "guess":
        try:
            structure = from_file(str(dest), label=Path(file.filename or "structure").stem)
        except Exception as e:
            raise HTTPException(400, f"Couldn't parse {file.filename}: {e}")

    _UPLOADS[token] = _Upload(kind=kind, path=str(dest), structure=structure)
    elements = sorted({a[0] for a in structure.atoms}) if structure else None
    return {"token": token, "elements": elements}


class SmilesRequest(BaseModel):
    smiles: str
    label: str = "structure"


@app.post("/api/analyze_smiles")
def analyze_smiles(req: SmilesRequest):
    """The GUI's equivalent of the CLI's SMILES intake option for the
    guess-geometry source, embeds with RDKit via geometry.from_smiles."""
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
    _UPLOADS[token] = _Upload(kind="guess", path=str(dest), structure=structure)
    elements = sorted({a[0] for a in structure.atoms})
    return {"token": token, "elements": elements}


@app.get("/api/gbasis_choices")
def gbasis_choices():
    return {"choices": GBASIS_CHOICES}


@app.get("/api/gaussian_levels")
def gaussian_levels():
    """Method/basis choices for the guess-geometry source's final-level
    Gaussian optimization, the same lists gaussbot's own GUI offers."""
    return {"methods": GAUSSIAN_METHODS, "basis_sets": GAUSSIAN_BASIS_SETS}


@app.get("/api/gamess_defaults")
def gamess_defaults():
    return {"rungms_path": DEFAULT_RUNGMS_PATH, "scratch_dir": DEFAULT_SCRATCH_DIR}


@app.get("/api/pbs_template")
def pbs_template():
    """The starting text for the GUI's editable PBS script box,
    job_runner.default_pbs_template()'s output, with __JOB__ standing
    in for the actual job name/`.inp` stem (substituted per-job at
    submission time, see job_runner.render_pbs_script())."""
    from .job_runner import default_pbs_template

    return {"template": default_pbs_template()}


@app.post("/api/jobs")
def create_job(req: JobRequest):
    if req.geometry_source not in ("gaussian_log", "gamess_inp", "guess"):
        raise HTTPException(400, "geometry_source must be 'gaussian_log', 'gamess_inp', or 'guess'")
    if not req.upload_token or req.upload_token not in _UPLOADS:
        raise HTTPException(400, "Unknown geometry source, upload/enter it again.")

    job_id = uuid.uuid4().hex
    job = JobState(id=job_id)
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
            "status": job.status, "stage": job.stage, "log": job.log,
            "result": job.result, "error": job.error,
        }


@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "Unknown job id")
        if job.status != "running":
            raise HTTPException(400, f"Job is {job.status}, not running, nothing to cancel.")
    job.cancel_event.set()
    return {"status": "cancelling"}


class ComboSpec(BaseModel):
    nstate: int = 3
    occ_selected: Optional[List[int]] = None  # omit both to accept the suggestion as-is
    virt_selected: Optional[List[int]] = None
    mem_mwords: int = DEFAULT_CASSCF_MEM_MWORDS  # separate from RHF/CIS's, multireference runs are far more memory-hungry
    run_xmcqdpt: bool = False  # dynamic correlation on top, once CASSCF converges, same active space, no further confirmation
    run_transitn: bool = False  # oscillator strengths (RUNTYP=TRANSITN) on top, once CASSCF converges, independent of run_xmcqdpt
    allow_smaller_active_space_recovery: bool = True  # per your request: if retrying with more iterations still fails, fall back to a smaller active space first, then regrow
    smaller_max_electrons: int = 4
    smaller_max_orbitals: int = 4
    mo_source: str = "rhf"  # "rhf" (start fresh from the closed-shell orbitals) or "previous" (start from the last converged combo's optimized orbitals)


class CASSCFBatchRequest(BaseModel):
    combos: List[ComboSpec]


def _run_one_combo(
    job: JobState, req: ComboSpec, active_space: ActiveSpaceSuggestion, mo_source: RHFOutcome, combo_index: int,
) -> Optional[CASSCFOutcome]:
    """Runs one active-space/state combination end to end (CASSCF ->
    optional XMCQDPT -> optional TRANSITN -> energy table), named
    cas-{e}{o}-sa{n}-c{combo_index} (and xpt-.../optical-... for its
    own follow-ups) so multiple combos can coexist in the same job's
    out_dir, per your request to compare several active-space/state
    results side by side. `combo_index` (see JobState.next_combo_index)
    is what actually guarantees the filenames never collide -- two
    combos sharing the same active space/nstate would otherwise
    produce the identical prefix and silently overwrite each other's
    files, which happened for real with "try another combination"
    reusing an earlier combo's active space. `mo_source` supplies the
    CASSCF starting orbitals, either job.rhf_outcome (closed-shell) or
    an RHFOutcome-shaped adapter built from a previous combo's
    optimized orbitals (see _run_casscf_batch_thread), per your
    request to optionally continue from a previous combination's
    optimized MOs instead of always restarting from the closed-shell
    orbitals. Appends its result to job.result["combos"] and returns
    the CASSCFOutcome (so the caller can offer it as the next combo's
    mo_source), None is never returned. Never raises for an expected
    GAMESS failure (non-convergence etc.), those are recorded on the
    combo so the batch can move on to the next one; only cancellation/
    unexpected errors propagate to the caller."""
    n_electrons = 2 * active_space.ndoc
    n_orbitals = active_space.ndoc + active_space.nval
    cas_prefix = f"cas-{n_electrons}{n_orbitals}-sa{req.nstate}-c{combo_index}"
    _set_stage(job, f"Running CASSCF ({req.nstate} states, NMCC={active_space.nmcc} "
                     f"NDOC={active_space.ndoc} NVAL={active_space.nval})")
    staged = run_casscf_staged(
        mo_source, active_space, job.out_dir, nstate=req.nstate,
        rungms_path=job.rungms_path, scratch_dir=job.scratch_dir, ncpus=job.ncpus,
        mem_mwords=req.mem_mwords, cancel_event=job.cancel_event, executor=job.executor,
        name_prefix=cas_prefix,
    )
    for line in staged.trail:
        _log(job, line)
    casscf_outcome = staged.outcome

    if not casscf_outcome.success and staged.exhausted and req.allow_smaller_active_space_recovery:
        _set_stage(job, f"CASSCF still hasn't converged, trying a smaller active space "
                         f"(max {req.smaller_max_electrons} electrons, {req.smaller_max_orbitals} orbitals)")
        casscf_outcome, active_space, recovery_trail = run_casscf_with_smaller_active_space_recovery(
            mo_source, active_space, job.out_dir, nstate=req.nstate,
            rungms_path=job.rungms_path, scratch_dir=job.scratch_dir,
            smaller_max_electrons=req.smaller_max_electrons, smaller_max_orbitals=req.smaller_max_orbitals,
            ncpus=job.ncpus, mem_mwords=req.mem_mwords, cancel_event=job.cancel_event, executor=job.executor,
            name_prefix=cas_prefix,
        )
        for line in recovery_trail:
            _log(job, line)
        n_electrons = 2 * active_space.ndoc
        n_orbitals = active_space.ndoc + active_space.nval

    combo: dict = {"casscf": _casscf_result_dict(casscf_outcome)}

    if not casscf_outcome.success:
        note = (
            "ran to completion but the orbital optimization did not converge, try a "
            "different active space or a different starting geometry"
            if casscf_outcome.result and casscf_outcome.result.normal_termination
            else "GAMESS exited abnormally, check the log, e.g. for an unviable active space"
        )
        combo["error"] = f"CASSCF did not succeed ({note}). See {casscf_outcome.log_path}."
        _log(job, combo["error"])
        with _JOBS_LOCK:
            job.result.setdefault("combos", []).append(combo)
        return casscf_outcome

    combo_tag = f"{n_electrons}{n_orbitals}-sa{req.nstate}-c{combo_index}"

    xmcqdpt_outcome = None
    if req.run_xmcqdpt:
        _set_stage(job, f"Running XMCQDPT ({req.nstate} states)")
        xmcqdpt_outcome = run_xmcqdpt(
            casscf_outcome, job.out_dir, nstate=req.nstate,
            rungms_path=job.rungms_path, scratch_dir=job.scratch_dir, ncpus=job.ncpus,
            mem_mwords=req.mem_mwords, cancel_event=job.cancel_event, executor=job.executor,
            name=f"xpt-{combo_tag}",
        )
        combo["xmcqdpt"] = _xmcqdpt_result_dict(xmcqdpt_outcome)
        if not xmcqdpt_outcome.success:
            _log(job, f"XMCQDPT did not succeed, see {xmcqdpt_outcome.log_path}. "
                      "Continuing with CASSCF-only energies in the table below.")
            xmcqdpt_outcome = None

    oscillator_strengths = None
    if req.run_transitn:
        _set_stage(job, f"Running TRANSITN (oscillator strengths, {req.nstate} states)")
        transitn_outcome = run_transitn(
            casscf_outcome, job.out_dir, nstate=req.nstate,
            rungms_path=job.rungms_path, scratch_dir=job.scratch_dir, ncpus=job.ncpus,
            mem_mwords=req.mem_mwords, cancel_event=job.cancel_event, executor=job.executor,
            name=f"optical-{combo_tag}",
        )
        combo["transitn"] = _transitn_result_dict(transitn_outcome)
        if transitn_outcome.success:
            oscillator_strengths = transitn_outcome.oscillator_strengths
        else:
            _log(job, f"TRANSITN did not succeed, see {transitn_outcome.log_path}.")

    combo["energy_table"] = _energy_table_dict(casscf_outcome, xmcqdpt_outcome, req.nstate, oscillator_strengths)

    with _JOBS_LOCK:
        job.result.setdefault("combos", []).append(combo)
    return casscf_outcome


def _run_casscf_batch_thread(job: JobState, combos: List[ComboSpec]) -> None:
    """Runs every requested active-space/state combination in order,
    one after another, without pausing for further confirmation,
    per your request that all combinations be specified upfront and
    then run straight through (CASSCF -> XMCQDPT -> TRANSITN as
    configured for each). Each combo's active space is looked up
    against the suggestion fresh, so one combo's smaller-active-space
    recovery never affects the next combo's active-space selection,
    but its *starting orbitals* can optionally carry forward from the
    previous combo's converged CASSCF orbitals (req.mo_source ==
    "previous"), per your request, tracked via job.last_casscf_outcome
    so this also works across separate "try another combination"
    batches, not just within one."""
    try:
        suggestion = job.active_space_suggestion
        for i, req in enumerate(combos):
            if job.cancel_event.is_set():
                raise GamessCancelledError()
            _log(job, f"--- Combination {i + 1} of {len(combos)},-")
            if req.occ_selected is not None or req.virt_selected is not None:
                active_space = build_active_space(
                    suggestion.n_occ, suggestion.norb,
                    req.occ_selected if req.occ_selected is not None else suggestion.occ_selected,
                    req.virt_selected if req.virt_selected is not None else suggestion.virt_selected,
                )
            else:
                active_space = suggestion

            if req.mo_source == "previous" and job.last_casscf_outcome is not None and job.last_casscf_outcome.success:
                mo_source = mo_source_from_casscf_outcome(job.last_casscf_outcome)
                _log(job, "Starting this combination's CASSCF from the previous combination's optimized orbitals.")
            else:
                if req.mo_source == "previous":
                    _log(job, "No converged previous combination to start from, using the closed-shell orbitals instead.")
                mo_source = job.rhf_outcome

            with _JOBS_LOCK:
                combo_index = job.next_combo_index
                job.next_combo_index += 1

            casscf_outcome = _run_one_combo(job, req, active_space, mo_source, combo_index)
            if casscf_outcome is not None and casscf_outcome.success:
                with _JOBS_LOCK:
                    job.last_casscf_outcome = casscf_outcome

        with _JOBS_LOCK:
            job.status = "done"
        _set_stage(job, "Done")
    except GamessCancelledError:
        with _JOBS_LOCK:
            job.status = "cancelled"
            job.log.append("Cancelled, the GAMESS process was terminated. Whatever finished before this point is kept above.")
    except Exception as e:  # noqa: BLE001
        _fail(job, f"Unexpected error: {e}")


@app.post("/api/jobs/{job_id}/casscf")
def confirm_casscf(job_id: str, req: CASSCFBatchRequest):
    """Continues a job that's either at 'awaiting_active_space' (the
    first batch) or 'done' (a "try another combination" batch added
    afterwards, see the GUI's try-another button), runs every
    requested active-space/state combination (each either the
    suggested active space as-is, or the user's own edited MO lists)
    one after another in a single background thread, per your request
    that all combinations in a batch be specified upfront rather than
    one at a time with a re-prompt in between."""
    if not req.combos:
        raise HTTPException(400, "Add at least one active-space/state combination.")
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "Unknown job id")
        if job.status not in ("awaiting_active_space", "done"):
            raise HTTPException(400, f"Job is {job.status}, not ready to run more combinations.")
        job.status = "running"
        job.cancel_event = threading.Event()

    thread = threading.Thread(target=_run_casscf_batch_thread, args=(job, req.combos), daemon=True)
    thread.start()
    return {"status": "running"}


@app.delete("/api/jobs/{job_id}")
def delete_job(job_id: str):
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            raise HTTPException(404, "Unknown job id")
        if job.status == "running":
            raise HTTPException(400, "Stop the study before removing it.")
        del _JOBS[job_id]
    return {"status": "deleted"}


def main() -> None:
    port = int(os.environ.get("GAMESSBOT_GUI_PORT", "8766"))
    url = f"http://127.0.0.1:{port}"
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    print(f"gamessbot GUI running at {url} (Ctrl+C to stop)")
    uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning")


if __name__ == "__main__":
    main()
