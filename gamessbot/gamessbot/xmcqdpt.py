"""
xmcqdpt.py

Runs the XMCQDPT (MCQDPT2) stage on top of a converged CASSCF: reuses
the same active space and the CASSCF punch's own "OPTIMIZED MCSCF MO-S"
orbitals (gamess_input.extract_optimized_mcscf_vec_block()), through
gamess_input.build_xmcqdpt_input() and the same local_runner.run_gamess()
everything else uses.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass
from typing import List, Optional

from .casscf import CASSCFOutcome
from .executor import run_gms
from .gamess_input import build_xmcqdpt_input
from .local_runner import check_normal_termination
from .parser import XMCQDPTResult, parse_xmcqdpt_log


@dataclass
class XMCQDPTOutcome:
    """What run_xmcqdpt() came back with. `success` requires the
    underlying CASSCF reference to have converged AND the perturbative
    "*** MCQDPT2 ENERGIES ***" section to have actually produced
    `nstate` state energies -- a run can terminate normally without
    reaching that section (e.g. if the CASSCF reference redone inside
    this job fails to converge, same as the standalone CASSCF stage)."""

    success: bool
    log_path: str
    result: Optional[XMCQDPTResult]


def run_xmcqdpt(
    casscf_outcome: CASSCFOutcome,
    out_dir: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    mem_mwords: int = 1,
    wstate: Optional[List[float]] = None,
    edshft: float = 0.04,
    xzero: bool = True,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
    name: str = "xmcqdpt",
) -> XMCQDPTOutcome:
    """
    Requires a successful/converged `casscf_outcome` with its
    optimized_vec_block set -- that's what this restarts from, using
    the same active space (casscf_outcome.active_space). `name`
    (default "xmcqdpt", giving "xmcqdpt.inp") lets multiple active-
    space/state combinations coexist in the same out_dir without
    overwriting each other's files -- e.g. "xpt-44-sa3".
    """
    if not casscf_outcome.success or casscf_outcome.optimized_vec_block is None:
        raise ValueError(
            "run_xmcqdpt() needs a successful, converged CASSCF outcome with optimized orbitals to restart from"
        )

    inp_text = build_xmcqdpt_input(
        casscf_outcome.data_block, charge=casscf_outcome.charge, mult=casscf_outcome.mult,
        gbasis_line=casscf_outcome.gbasis_line, vec_block=casscf_outcome.optimized_vec_block,
        norb=casscf_outcome.norb, active_space=casscf_outcome.active_space, nstate=nstate,
        wstate=wstate, mem_mwords=mem_mwords, edshft=edshft, xzero=xzero,
        vec_annotation=casscf_outcome.optimized_vec_annotation,
    )
    inp_path = os.path.join(out_dir, f"{name}.inp")
    with open(inp_path, "w") as f:
        f.write(inp_text)

    log_path = run_gms(inp_path, rungms_path, scratch_dir, ncpus=ncpus, executor=executor, cancel_event=cancel_event)

    ran_cleanly = check_normal_termination(log_path)
    result = parse_xmcqdpt_log(log_path, nstate=nstate) if ran_cleanly else None
    success = (
        ran_cleanly and result is not None and result.casscf_converged
        and len(result.mcqdpt_state_energies) == nstate
    )
    return XMCQDPTOutcome(success=success, log_path=log_path, result=result)
