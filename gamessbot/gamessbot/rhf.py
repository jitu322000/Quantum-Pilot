"""
rhf.py

Runs the RHF stage with the convergence ladder you described: SOSCF is
tried first (the real GAMESS default for RHF). If it fails to converge,
fall back to DIIS from a fresh guess. If DIIS converges, take its
orbitals and retry SOSCF once more (SOSCF usually gives a better-behaved
final wavefunction than DIIS, so it's worth a second attempt once there's
a decent starting guess) -- if that also fails, the DIIS result stands
as the final answer, since it's already a valid converged RHF.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import List, Optional

from .executor import run_gms
from .gamess_input import build_rhf_input, extract_vec_annotation, extract_vec_block
from .local_runner import GamessCancelledError, GamessRunError, check_normal_termination
from .parser import parse_rhf_log


@dataclass
class RHFOutcome:
    """What run_rhf_staged() came back with."""

    success: bool
    log_path: Optional[str]
    energy: Optional[float]
    norb: Optional[int]
    vec_block: Optional[str]  # the converged $VEC...$END block, for CIS's GUESS=MOREAD
    data_block: str
    charge: int
    mult: int
    gbasis_line: str
    trail: List[str] = field(default_factory=list)  # human-readable account of what was tried
    vec_annotation: str = ""  # the "--- CLOSED SHELL ORBITALS --- ..." heading above vec_block, if any


def _attempt(
    inp_name: str,
    out_dir: str,
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    mem_mwords: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int,
    use_soscf: bool,
    guess: str,
    vec_block: Optional[str],
    norb: Optional[int],
    cancel_event: Optional[threading.Event],
    executor: str = "local",
    vec_annotation: str = "",
) -> tuple:
    """Builds, runs, and checks one RHF attempt. Returns (success, log_path)."""
    inp_text = build_rhf_input(
        data_block, charge=charge, mult=mult, gbasis_line=gbasis_line,
        mem_mwords=mem_mwords, use_soscf=use_soscf, guess=guess,
        vec_block=vec_block, norb=norb, vec_annotation=vec_annotation,
    )
    inp_path = os.path.join(out_dir, inp_name)
    with open(inp_path, "w") as f:
        f.write(inp_text)

    log_path = run_gms(inp_path, rungms_path, scratch_dir, ncpus=ncpus, executor=executor, cancel_event=cancel_event)
    return check_normal_termination(log_path), log_path


def run_rhf_staged(
    data_block: str,
    out_dir: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    mem_mwords: int = 1,
    use_soscf: bool = True,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> RHFOutcome:
    """
    `use_soscf` is the user's own choice of which convergence method to
    start from (SOSCF is the default). If they chose DIIS outright,
    that's just run directly -- the automatic SOSCF->DIIS->SOSCF(MOREAD)
    ladder only kicks in when starting from SOSCF and it fails to
    converge on its own, per your description ("if it doesn't [converge]
    at soscf ... it should automatically switch to DIIS").

    Raises GamessRunError/GamessCancelledError straight through if
    rungms itself can't run (missing executable, cancelled) -- those
    aren't convergence failures, they're setup/execution failures the
    caller needs to see immediately.
    """
    trail: List[str] = []

    if not use_soscf:
        ok, log_path = _attempt(
            "rhf_diis.inp", out_dir, data_block, charge, mult, gbasis_line, mem_mwords,
            rungms_path, scratch_dir, ncpus, use_soscf=False, guess="HUCKEL",
            vec_block=None, norb=None, cancel_event=cancel_event, executor=executor,
        )
        if not ok:
            trail.append("DIIS (HUCKEL guess, user-selected): did not converge")
            return RHFOutcome(
                success=False, log_path=log_path, energy=None, norb=None, vec_block=None,
                data_block=data_block, charge=charge, mult=mult, gbasis_line=gbasis_line, trail=trail,
            )
        trail.append("DIIS (HUCKEL guess, user-selected): converged")
        result = parse_rhf_log(log_path)
        with open(os.path.join(out_dir, "rhf_diis.dat")) as f:
            dat_text = f.read()
        vec_block = extract_vec_block(dat_text)
        return RHFOutcome(
            success=True, log_path=log_path, energy=result.energy, norb=result.norb,
            vec_block=vec_block, data_block=data_block, charge=charge, mult=mult,
            gbasis_line=gbasis_line, trail=trail, vec_annotation=extract_vec_annotation(dat_text),
        )

    ok, log_path = _attempt(
        "rhf_soscf.inp", out_dir, data_block, charge, mult, gbasis_line, mem_mwords,
        rungms_path, scratch_dir, ncpus, use_soscf=True, guess="HUCKEL",
        vec_block=None, norb=None, cancel_event=cancel_event, executor=executor,
    )
    if ok:
        trail.append("SOSCF (HUCKEL guess): converged")
        result = parse_rhf_log(log_path)
        with open(os.path.join(out_dir, "rhf_soscf.dat")) as f:
            dat_text = f.read()
        vec_block = extract_vec_block(dat_text)
        return RHFOutcome(
            success=True, log_path=log_path, energy=result.energy, norb=result.norb,
            vec_block=vec_block, data_block=data_block, charge=charge, mult=mult,
            gbasis_line=gbasis_line, trail=trail, vec_annotation=extract_vec_annotation(dat_text),
        )
    trail.append("SOSCF (HUCKEL guess): did not converge, falling back to DIIS")

    ok, diis_log_path = _attempt(
        "rhf_diis.inp", out_dir, data_block, charge, mult, gbasis_line, mem_mwords,
        rungms_path, scratch_dir, ncpus, use_soscf=False, guess="HUCKEL",
        vec_block=None, norb=None, cancel_event=cancel_event, executor=executor,
    )
    if not ok:
        trail.append("DIIS (HUCKEL guess): did not converge either -- giving up")
        return RHFOutcome(
            success=False, log_path=diis_log_path, energy=None, norb=None, vec_block=None,
            data_block=data_block, charge=charge, mult=mult, gbasis_line=gbasis_line, trail=trail,
        )
    trail.append("DIIS (HUCKEL guess): converged")

    diis_result = parse_rhf_log(diis_log_path)
    with open(os.path.join(out_dir, "rhf_diis.dat")) as f:
        diis_dat_text = f.read()
    diis_vec_block = extract_vec_block(diis_dat_text)
    diis_vec_annotation = extract_vec_annotation(diis_dat_text)

    ok, moread_log_path = _attempt(
        "rhf_soscf_moread.inp", out_dir, data_block, charge, mult, gbasis_line, mem_mwords,
        rungms_path, scratch_dir, ncpus, use_soscf=True, guess="MOREAD",
        vec_block=diis_vec_block, norb=diis_result.norb, cancel_event=cancel_event, executor=executor,
        vec_annotation=diis_vec_annotation,
    )
    if ok:
        trail.append("SOSCF (MOREAD from DIIS orbitals): converged")
        result = parse_rhf_log(moread_log_path)
        with open(os.path.join(out_dir, "rhf_soscf_moread.dat")) as f:
            dat_text = f.read()
        vec_block = extract_vec_block(dat_text)
        return RHFOutcome(
            success=True, log_path=moread_log_path, energy=result.energy, norb=result.norb,
            vec_block=vec_block, data_block=data_block, charge=charge, mult=mult,
            gbasis_line=gbasis_line, trail=trail, vec_annotation=extract_vec_annotation(dat_text),
        )

    trail.append("SOSCF (MOREAD from DIIS orbitals): did not converge -- keeping the DIIS result as final")
    return RHFOutcome(
        success=True, log_path=diis_log_path, energy=diis_result.energy, norb=diis_result.norb,
        vec_block=diis_vec_block, data_block=data_block, charge=charge, mult=mult,
        gbasis_line=gbasis_line, trail=trail, vec_annotation=diis_vec_annotation,
    )
