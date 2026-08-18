"""
transitn.py

Runs the RUNTYP=TRANSITN (oscillator strength) stage on top of a
converged CASSCF: restarts from the same CASSCF orbitals xmcqdpt.py
uses (the punch's own "OPTIMIZED MCSCF MO-S" block), relabeled to
$VEC1 (gamess_input.relabel_vec_group()) as GAMESS's own $TRANST group
requires, through gamess_input.build_transitn_input() and the same
executor.run_gms() everything else uses.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .casscf import CASSCFOutcome
from .executor import run_gms
from .gamess_input import build_transitn_input, relabel_vec_group
from .local_runner import check_normal_termination
from .parser import parse_transitn_log


@dataclass
class TransitnOutcome:
    """What run_transitn() came back with. `success` only means GAMESS
    ran to completion (TERMINATED NORMALLY) -- oscillator_strengths is
    keyed by state_index (1 = S0, 2 = S1, ...), matching the same
    1-based convention CASSCF/XMCQDPT use elsewhere; S0 itself is never
    a key (see parser.parse_transitn_log())."""

    success: bool
    log_path: str
    oscillator_strengths: Dict[int, float] = field(default_factory=dict)


def run_transitn(
    casscf_outcome: CASSCFOutcome,
    out_dir: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    mem_mwords: int = 1,
    wstate: Optional[List[float]] = None,
    itermx: int = 120,
    name: str = "transitn",
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> TransitnOutcome:
    """
    Requires a successful/converged `casscf_outcome` with its
    optimized_vec_block set -- that's what this restarts from (the
    same orbitals xmcqdpt.run_xmcqdpt() uses), reusing the same active
    space (casscf_outcome.active_space) and state count, relabeled to
    $VEC1 as GAMESS's own $TRANST group requires (see
    gamess_input.build_transitn_input()'s docstring).
    """
    if not casscf_outcome.success or casscf_outcome.optimized_vec_block is None:
        raise ValueError(
            "run_transitn() needs a successful, converged CASSCF outcome with optimized orbitals to restart from"
        )

    vec1_block = relabel_vec_group(casscf_outcome.optimized_vec_block, "VEC1")

    inp_text = build_transitn_input(
        casscf_outcome.data_block, charge=casscf_outcome.charge, mult=casscf_outcome.mult,
        gbasis_line=casscf_outcome.gbasis_line, vec_block=vec1_block,
        norb=casscf_outcome.norb, active_space=casscf_outcome.active_space, nstate=nstate,
        wstate=wstate, mem_mwords=mem_mwords, itermx=itermx,
        vec_annotation=casscf_outcome.optimized_vec_annotation,
    )
    inp_path = os.path.join(out_dir, f"{name}.inp")
    with open(inp_path, "w") as f:
        f.write(inp_text)

    log_path = run_gms(inp_path, rungms_path, scratch_dir, ncpus=ncpus, executor=executor, cancel_event=cancel_event)

    success = check_normal_termination(log_path)
    oscillator_strengths = parse_transitn_log(log_path, nstate=nstate) if success else {}
    return TransitnOutcome(success=success, log_path=log_path, oscillator_strengths=oscillator_strengths)
