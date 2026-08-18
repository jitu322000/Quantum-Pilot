"""
parser.py

Reads finished GAMESS .log files and pulls out what the pipeline needs:
did RHF terminate normally and converge (and what NORB is, for the CIS
stage's GUESS=MOREAD), and what excited states CIS found.

Hand-rolled regex parsing against fixed-format lines GAMESS prints,
confirmed against real logs on this machine (a GAMESS example log and a
real CIS run) -- same practice as gaussbot's own parser.py.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

_NORB_RE = re.compile(r"NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS\s*=\s*(\d+)")
_VARIATION_SPACE_NORB_RE = re.compile(r"TOTAL NUMBER OF MOS IN VARIATION SPACE\s*=\s*(\d+)")
_NOCC_RE = re.compile(r"NUMBER OF OCCUPIED ORBITALS \(ALPHA\)\s*=\s*(\d+)")
_RHF_ENERGY_RE = re.compile(r"FINAL RHF ENERGY IS\s+(-?\d+\.\d+)")

_EXCITED_STATE_RE = re.compile(
    r"^\s*EXCITED STATE\s+(\d+)\s+ENERGY=\s+(-?\d+\.\d+)\s+S\s*=\s*(-?\d+\.\d+)\s+SPACE SYM\s*=\s*(\S+)"
)
_TRANSITION_ROW_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s*$")

_FINAL_MCSCF_ENERGY_RE = re.compile(r"FINAL MCSCF ENERGY IS\s+(-?\d+\.\d+)\s+AFTER\s+(\d+)\s+ITERATIONS")
_CASSCF_STATE_RE = re.compile(r"^\s*STATE #\s*(\d+)\s+ENERGY\s*=\s*(-?\d+\.\d+)")

_MCQDPT_ENERGIES_MARKER = "*** MCQDPT2 ENERGIES ***"

_CI_STATE_PAIR_RE = re.compile(r"^\s*CI STATE NUMBER=\s*(\d+)\s+(\d+)\s+STATE MULTIPLICITY=")
_OSCILLATOR_STRENGTH_RE = re.compile(r"OSCILLATOR STRENGTH\s*=\s*(-?\d+\.\d+)")


@dataclass
class RHFResult:
    """What parse_rhf_log() pulled out of an RHF .log file."""

    normal_termination: bool
    energy: Optional[float] = None
    norb: Optional[int] = None
    n_occ: Optional[int] = None  # number of doubly-occupied orbitals -- the NMCC+NDOC boundary for CASSCF


@dataclass
class Transition:
    """One FROM MO -> TO MO contribution to a CIS excited state."""

    from_mo: int
    to_mo: int
    coefficient: float


@dataclass
class ExcitedState:
    """One CIS excited state block."""

    index: int
    energy: float
    spin: float
    space_sym: str
    transitions: List[Transition] = field(default_factory=list)


def parse_rhf_log(log_path: str) -> RHFResult:
    """
    normal_termination is the ground truth for whether the RHF run
    succeeded (confirmed exact string: "TERMINATED NORMALLY" -- see
    local_runner.check_normal_termination()). energy/norb are only
    meaningful when normal_termination is True.
    """
    with open(log_path) as f:
        text = f.read()

    normal_termination = "TERMINATED NORMALLY" in text

    energy = None
    m = _RHF_ENERGY_RE.search(text)
    if m is not None:
        energy = float(m.group(1))

    # With ISPHER=1 (gamess_input._needs_ispher() -- correlation-consistent
    # and other spherical-only basis families), GAMESS drops "spherical
    # contaminant" AOs (e.g. cartesian d has 6 components, spherical d has
    # 5) and prints the real MOREAD-relevant count as "TOTAL NUMBER OF MOS
    # IN VARIATION SPACE=" -- confirmed against a real ethene/cc-pVDZ log
    # where this was 38, not the 40 cartesian basis functions, and using
    # the cartesian count for $GUESS ... NORB=/$VEC caused a real rungms
    # crash ("PREMATURE END OF ORBITAL INPUT ... LOOKING FOR ORBITAL 39")
    # since the punched $VEC block only ever has 38 MOs' worth of rows.
    # Basis families that don't trigger ISPHER (no d/f shells, or ISPHER
    # not needed) never print this line, so the cartesian count -- which
    # already equals the variation-space size there -- is used instead.
    norb = None
    m = _VARIATION_SPACE_NORB_RE.search(text)
    if m is not None:
        norb = int(m.group(1))
    else:
        m = _NORB_RE.search(text)
        if m is not None:
            norb = int(m.group(1))

    n_occ = None
    m = _NOCC_RE.search(text)
    if m is not None:
        n_occ = int(m.group(1))

    return RHFResult(normal_termination=normal_termination, energy=energy, norb=norb, n_occ=n_occ)


def parse_cis_log(log_path: str) -> List[ExcitedState]:
    """
    Parses every "EXCITED STATE n  ENERGY=...  S = ...  SPACE SYM = ..."
    block and the FROM MO / TO MO / SAP COEFFICIENT table under it, in
    order. Returns an empty list if the run didn't reach the CIS section
    at all (e.g. it failed) -- callers should check
    local_runner.check_normal_termination() first.
    """
    with open(log_path) as f:
        lines = f.readlines()

    states: List[ExcitedState] = []
    current: Optional[ExcitedState] = None
    in_table = False

    for line in lines:
        m = _EXCITED_STATE_RE.match(line)
        if m is not None:
            current = ExcitedState(
                index=int(m.group(1)),
                energy=float(m.group(2)),
                spin=float(m.group(3)),
                space_sym=m.group(4),
            )
            states.append(current)
            in_table = False
            continue

        if current is None:
            continue

        if "FROM MO" in line and "TO MO" in line:
            in_table = True
            continue

        if in_table:
            m = _TRANSITION_ROW_RE.match(line)
            if m is not None:
                current.transitions.append(
                    Transition(
                        from_mo=int(m.group(1)),
                        to_mo=int(m.group(2)),
                        coefficient=float(m.group(3)),
                    )
                )
            elif line.strip() and "---" not in line:
                in_table = False

    return states


@dataclass
class CASSCFResult:
    """What parse_casscf_log() pulled out of a CASSCF .log file."""

    normal_termination: bool
    converged: bool  # did the orbital optimization itself converge, not just "ran without crashing"
    final_energy: Optional[float] = None
    iterations: Optional[int] = None
    state_energies: List[Tuple[int, float]] = field(default_factory=list)  # (state index, energy)


def parse_casscf_log(log_path: str, nstate: Optional[int] = None) -> CASSCFResult:
    """
    normal_termination is only "GAMESS didn't crash" (same "TERMINATED
    NORMALLY" check as RHF/CIS) -- it is NOT evidence the MCSCF orbital
    optimization actually converged. Confirmed against a real run that
    exhausted MAXIT without converging: GAMESS still terminates
    normally, but prints "MCSCF IS NOT CONVERGED!" and a placeholder
    "FINAL MCSCF ENERGY IS 0.0000000000" (not a real energy -- ignore it
    when `converged` is False). The real convergence marker is
    "LAGRANGIAN CONVERGED".

    state_energies are the FINAL per-state energies -- confirmed against
    a real CASSCF log, GAMESS prints "STATE #  n  ENERGY = ..." blocks
    twice (once mid-iteration, once after "FINAL MCSCF ENERGY IS ...
    AFTER n ITERATIONS" as the "-MCCI- BASED ON OPTIMIZED ORBITALS"
    analysis) -- only the ones after that line are used here, though
    when `converged` is False these are the best the optimizer reached,
    not truly converged values either.

    `nstate`, when given, stops collecting state_energies once that many
    have been found. Pass it whenever the caller knows NSTATE -- without
    it, a log with more content after the state block (e.g. a combined
    CASSCF+XMCQDPT run, which reprints "STATE # n ENERGY=" blocks
    several times further down for its own restart/analysis steps,
    confirmed against a real such log) would otherwise sweep up those
    unrelated later blocks too.
    """
    with open(log_path) as f:
        text = f.read()

    normal_termination = "TERMINATED NORMALLY" in text
    converged = "LAGRANGIAN CONVERGED" in text

    final_energy = None
    iterations = None
    m = _FINAL_MCSCF_ENERGY_RE.search(text)
    if m is not None:
        iterations = int(m.group(2))
        if converged:
            final_energy = float(m.group(1))

    state_energies: List[Tuple[int, float]] = []
    tail = text[m.end():] if m is not None else text
    for line in tail.splitlines():
        sm = _CASSCF_STATE_RE.match(line)
        if sm is not None:
            state_energies.append((int(sm.group(1)), float(sm.group(2))))
            if nstate is not None and len(state_energies) >= nstate:
                break

    return CASSCFResult(
        normal_termination=normal_termination, converged=converged, final_energy=final_energy,
        iterations=iterations, state_energies=state_energies,
    )


@dataclass
class XMCQDPTResult:
    """What parse_xmcqdpt_log() pulled out of an XMCQDPT (MCQDPT2) .log
    file -- a GAMESS XMCQDPT run redoes the CASSCF reference calculation
    in the same job before the perturbative correction, so both sets of
    state energies are available here."""

    normal_termination: bool
    casscf_converged: bool
    casscf_state_energies: List[Tuple[int, float]] = field(default_factory=list)
    mcqdpt_state_energies: List[Tuple[int, float]] = field(default_factory=list)


def parse_xmcqdpt_log(log_path: str, nstate: int) -> XMCQDPTResult:
    """
    Confirmed against a real XMCQDPT log: the CASSCF reference section
    (same "FINAL MCSCF ENERGY IS" + "STATE # n ENERGY=" shape parse_casscf_log
    already handles) comes first, followed later by a
    "*** MCQDPT2 ENERGIES ***" marker ("(FROM DIAGONALIZATION OF 2ND
    ORDER EFFECTIVE HAMILTONIAN)") and then the final perturbatively-
    corrected "STATE # n ENERGY=" block -- the numbers actually meant
    for downstream use ("These are the energies that should be used for
    subsequent analysis"). `nstate` bounds both scans the same way
    parse_casscf_log() does, since GAMESS reprints state-energy-shaped
    blocks multiple times further down in a combined log like this one.
    """
    casscf = parse_casscf_log(log_path, nstate=nstate)

    with open(log_path) as f:
        text = f.read()
    normal_termination = "TERMINATED NORMALLY" in text

    mcqdpt_state_energies: List[Tuple[int, float]] = []
    idx = text.find(_MCQDPT_ENERGIES_MARKER)
    if idx != -1:
        for line in text[idx:].splitlines():
            sm = _CASSCF_STATE_RE.match(line)
            if sm is not None:
                mcqdpt_state_energies.append((int(sm.group(1)), float(sm.group(2))))
                if len(mcqdpt_state_energies) >= nstate:
                    break

    return XMCQDPTResult(
        normal_termination=normal_termination, casscf_converged=casscf.converged,
        casscf_state_energies=casscf.state_energies, mcqdpt_state_energies=mcqdpt_state_energies,
    )


def parse_transitn_log(log_path: str, nstate: int) -> Dict[int, float]:
    """
    Parses oscillator strengths from a GAMESS RUNTYP=TRANSITN log
    (OPERAT=DM). GAMESS prints one "CI STATE NUMBER=  i  j STATE
    MULTIPLICITY= ..." block per (bra, ket) state pair -- this looks
    specifically for the "1 j" pairs (transitions FROM state 1, i.e.
    S0), matching the 1-based state_index convention CASSCF/XMCQDPT
    already use elsewhere in this package (state_index=1 is S0,
    state_index=j is S(j-1)), and reads the "OSCILLATOR STRENGTH ="
    line that follows within that same block.

    State 1 vs itself has no such line -- GAMESS instead prints "THE
    NEXT PAIR ARE THE SAME STATE, SO THIS IS AN EXPECTATION VALUE,
    RATHER THAN A TRANSITION MOMENT" -- so S0 is simply absent from the
    returned dict; callers should treat a missing state_index as "no
    oscillator strength for this state" (true for S0 itself), not a
    parsing failure.

    Confirmed against a real production TRANSITN log from this user's
    own prior work.
    """
    with open(log_path) as f:
        lines = f.readlines()

    result: Dict[int, float] = {}
    current_ket: Optional[int] = None
    for line in lines:
        m = _CI_STATE_PAIR_RE.match(line)
        if m is not None:
            bra, ket = int(m.group(1)), int(m.group(2))
            current_ket = ket if bra == 1 else None
            continue
        if current_ket is not None:
            fm = _OSCILLATOR_STRENGTH_RE.search(line)
            if fm is not None:
                result[current_ket] = float(fm.group(1))
                current_ket = None  # only one OSCILLATOR STRENGTH line per block
        if len(result) >= nstate - 1:
            break
    return result
