"""
cis.py

Runs CIS on top of an already-converged RHF: reads the RHF's own $VEC as
GUESS=MOREAD (so this converges fast, per your description) and requests
`nstate` singlet excited states.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import List, Optional

from .executor import run_gms
from .gamess_input import build_cis_input
from .local_runner import check_normal_termination
from .parser import ExcitedState, parse_cis_log
from .rhf import RHFOutcome


@dataclass
class CISOutcome:
    """What run_cis() came back with."""

    success: bool
    log_path: str
    states: List[ExcitedState] = field(default_factory=list)


def run_cis(
    rhf_outcome: RHFOutcome,
    out_dir: str,
    nstate: int,
    rungms_path: str,
    scratch_dir: str,
    ncpus: int = 1,
    mem_mwords: int = 1,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> CISOutcome:
    """
    Requires a successful `rhf_outcome` (its data_block/charge/mult/
    gbasis_line/norb/vec_block feed the CIS input directly). Raises
    ValueError if `rhf_outcome.success` is False -- there's nothing to
    build CIS's GUESS=MOREAD from otherwise.
    """
    if not rhf_outcome.success:
        raise ValueError("run_cis() needs a successful RHF outcome to read orbitals from")

    inp_text = build_cis_input(
        rhf_outcome.data_block, charge=rhf_outcome.charge, mult=rhf_outcome.mult,
        gbasis_line=rhf_outcome.gbasis_line, nstate=nstate, norb=rhf_outcome.norb,
        vec_block=rhf_outcome.vec_block, mem_mwords=mem_mwords,
        vec_annotation=rhf_outcome.vec_annotation,
    )
    inp_path = os.path.join(out_dir, "cis.inp")
    with open(inp_path, "w") as f:
        f.write(inp_text)

    log_path = run_gms(inp_path, rungms_path, scratch_dir, ncpus=ncpus, executor=executor, cancel_event=cancel_event)

    success = check_normal_termination(log_path)
    states = parse_cis_log(log_path) if success else []
    return CISOutcome(success=success, log_path=log_path, states=states)
