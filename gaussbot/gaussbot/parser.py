"""
parser.py

Reads a finished Gaussian .log file and pulls out what the pipeline
needs to decide its next move: did it terminate normally, did the
optimization converge, and -- the real test of whether an Opt Freq
result is a genuine minimum -- are all the vibrational frequencies
positive, or is at least one imaginary (negative)?

Hand-rolled regex parsing rather than cclib: the only things the
pipeline actually branches on (imaginary frequencies, final geometry,
SCF energy) are a handful of fixed-format lines Gaussian has printed
unchanged since G03, so a full log-format library is more than this
needs.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from rdkit import Chem

from .geometry import Atom

_PERIODIC_TABLE = Chem.GetPeriodicTable()

_SCF_RE = re.compile(r"SCF Done:\s+E\([^)]*\)\s*=\s*(-?\d+\.\d+)")
_ZPE_SUM_RE = re.compile(r"Sum of electronic and zero-point Energies=\s*(-?\d+\.\d+)")
_FREQ_RE = re.compile(r"^\s*Frequencies\s+--\s+(.+)$")
_ORIENTATION_HEADER_RE = re.compile(r"(Standard|Input) orientation:")
_ORIENTATION_ROW_RE = re.compile(
    r"^\s*\d+\s+(\d+)\s+-?\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)
_ATOM_DISP_HEADER_RE = re.compile(r"^\s*Atom\s+AN\s+X\s+Y\s+Z")
_IRC_PATH2_RE = re.compile(r"Path Number:\s*2\b")
_CONVERGENCE_ROW_RE = re.compile(
    r"^\s*(Maximum Force|RMS\s+Force|Maximum Displacement|RMS\s+Displacement)"
    r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(YES|NO)\s*$"
)
_CHARGE_MULT_RE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")


@dataclass
class LogResult:
    """What parse_log() pulled out of a Gaussian .log file."""

    normal_termination: bool
    stationary_point_found: bool
    scf_energy: Optional[float]
    frequencies: List[float] = field(default_factory=list)
    final_geometry: Optional[List[Atom]] = None
    zpe_energy: Optional[float] = None  # "Sum of electronic and zero-point Energies"
    charge: Optional[int] = None
    multiplicity: Optional[int] = None

    @property
    def imaginary_freqs(self) -> List[float]:
        return [f for f in self.frequencies if f < 0]

    @property
    def is_minimum(self) -> bool:
        """True if this looks like a genuine minimum: converged, and every
        vibrational mode is real (no imaginary frequencies)."""
        return (
            self.normal_termination
            and self.stationary_point_found
            and bool(self.frequencies)
            and not self.imaginary_freqs
        )

    @property
    def energy(self) -> Optional[float]:
        """The energy to use for barrier/reaction-energy reporting: ZPE-
        corrected (the standard 0 K reference) when available, the raw
        SCF electronic energy otherwise."""
        return self.zpe_energy if self.zpe_energy is not None else self.scf_energy


@dataclass
class ConvergenceItem:
    """One row of Gaussian's optimization convergence table (e.g.
    'Maximum Force  0.043770  0.000450  NO')."""

    value: float
    threshold: float
    converged: bool


@dataclass
class ConvergenceStatus:
    """The four criteria from Gaussian's convergence table, from the
    most recent 'Item Value Threshold Converged?' block in a log --
    whichever ones the log actually printed (a crashed-before-any-step
    log leaves all four None)."""

    max_force: Optional[ConvergenceItem] = None
    rms_force: Optional[ConvergenceItem] = None
    max_displacement: Optional[ConvergenceItem] = None
    rms_displacement: Optional[ConvergenceItem] = None

    @property
    def items(self) -> List[Tuple[str, ConvergenceItem]]:
        return [
            (name, item)
            for name, item in [
                ("Maximum Force", self.max_force),
                ("RMS Force", self.rms_force),
                ("Maximum Displacement", self.max_displacement),
                ("RMS Displacement", self.rms_displacement),
            ]
            if item is not None
        ]

    @property
    def n_converged(self) -> int:
        return sum(1 for _, item in self.items if item.converged)

    @property
    def n_total(self) -> int:
        return len(self.items)


@dataclass
class IRCResult:
    """What parse_irc_log() pulled out of a Gaussian IRC .log file."""

    normal_termination: bool
    completed: bool  # "Reaction path calculation complete" found
    forward_endpoint: Optional[List[Atom]]
    reverse_endpoint: Optional[List[Atom]]


def _parse_final_geometry(lines: List[str]) -> Optional[List[Atom]]:
    """Find the last orientation block in the log and read off its atoms.

    For an Opt Freq job this is the geometry the frequency analysis
    (and thus the "is this a minimum" check) was actually run on.
    """
    header_idx = None
    for i, line in enumerate(lines):
        if _ORIENTATION_HEADER_RE.search(line):
            header_idx = i

    if header_idx is None:
        return None

    # header, dashes, 2 column-label lines, dashes -- then atom rows.
    row_start = header_idx + 5
    atoms: List[Atom] = []
    for line in lines[row_start:]:
        m = _ORIENTATION_ROW_RE.match(line)
        if not m:
            break
        atomic_num, x, y, z = m.groups()
        symbol = _PERIODIC_TABLE.GetElementSymbol(int(atomic_num))
        atoms.append((symbol, float(x), float(y), float(z)))

    return atoms or None


def last_geometry(log_path: str) -> Optional[List[Atom]]:
    """
    Read the last orientation block in a log, whatever state the job
    ended up in -- a converged minimum, a TS, or an optimization that
    ran out of cycles without converging. Used to continue a
    non-converged optimization from wherever it actually got to,
    rather than restarting from the original guess every time.
    """
    with open(log_path) as f:
        lines = f.readlines()
    return _parse_final_geometry(lines)


_CONVERGENCE_FIELD_BY_LABEL = {
    "Maximum Force": "max_force",
    "RMS Force": "rms_force",
    "Maximum Displacement": "max_displacement",
    "RMS Displacement": "rms_displacement",
}


def last_convergence_status(log_path: str) -> Optional[ConvergenceStatus]:
    """
    The LAST 'Item Value Threshold Converged?' block in the log -- the
    most recent convergence check Gaussian actually ran, whether or not
    the job went on to terminate normally. Same "read whatever's there,
    converged or not" spirit as last_geometry() -- used to decide *how*
    to repair a non-converged optimization (pipeline.py) instead of
    always reacting the same way regardless of how close it actually
    got. Returns None if the log has no convergence table at all (e.g.
    it crashed before any optimization step ran).
    """
    with open(log_path) as f:
        lines = f.readlines()

    last_block_start = None
    for i, line in enumerate(lines):
        if "Converged?" in line and "Threshold" in line:
            last_block_start = i
    if last_block_start is None:
        return None

    fields: dict = {}
    for line in lines[last_block_start + 1 : last_block_start + 5]:
        m = _CONVERGENCE_ROW_RE.match(line)
        if not m:
            continue
        label, value, threshold, converged = m.groups()
        field_name = _CONVERGENCE_FIELD_BY_LABEL[re.sub(r"\s+", " ", label)]
        fields[field_name] = ConvergenceItem(float(value), float(threshold), converged == "YES")

    return ConvergenceStatus(**fields) if fields else None


def _iter_freq_blocks(lines: List[str]):
    """
    Yield (frequencies, atom_displacement_rows) for each block in a
    Gaussian frequency printout. `frequencies` is the list of (up to 3)
    mode frequencies in that block, in left-to-right column order.
    `atom_displacement_rows[i]` is the flat list of 3*len(frequencies)
    displacement floats for atom i, in the same column order.
    """
    n = len(lines)
    i = 0
    while i < n:
        m = _FREQ_RE.match(lines[i])
        if not m:
            i += 1
            continue
        freqs = [float(v) for v in m.group(1).split()]
        n_modes = len(freqs)

        header_idx = None
        for j in range(i + 1, min(i + 10, n)):
            if _ATOM_DISP_HEADER_RE.match(lines[j]):
                header_idx = j
                break
        if header_idx is None:
            i += 1
            continue

        atom_rows: List[List[float]] = []
        k = header_idx + 1
        while k < n:
            parts = lines[k].split()
            if len(parts) != 2 + 3 * n_modes:
                break
            try:
                atom_rows.append([float(v) for v in parts[2:]])
            except ValueError:
                break
            k += 1

        yield freqs, atom_rows
        i = k


def imaginary_mode_displacement(log_path: str, which: int = 0) -> Optional[List[Tuple[float, float, float]]]:
    """
    Return the per-atom Cartesian displacement vector for the `which`-th
    imaginary vibrational mode in the log (0 = the first/most negative
    one, which is the one that actually matters almost always), or None
    if the log has no imaginary frequency.
    """
    with open(log_path) as f:
        lines = f.readlines()

    seen = 0
    for freqs, atom_rows in _iter_freq_blocks(lines):
        for col, freq in enumerate(freqs):
            if freq < 0:
                if seen == which:
                    return [
                        (row[col * 3], row[col * 3 + 1], row[col * 3 + 2])
                        for row in atom_rows
                    ]
                seen += 1
    return None


def parse_log(log_path: str) -> LogResult:
    """Parse a finished Gaussian .log file."""
    with open(log_path) as f:
        text = f.read()
    lines = text.splitlines()

    normal_termination = "Normal termination of Gaussian" in text
    stationary_point_found = "Stationary point found" in text

    scf_matches = _SCF_RE.findall(text)
    scf_energy = float(scf_matches[-1]) if scf_matches else None

    zpe_matches = _ZPE_SUM_RE.findall(text)
    zpe_energy = float(zpe_matches[-1]) if zpe_matches else None

    frequencies: List[float] = []
    for line in lines:
        m = _FREQ_RE.match(line)
        if m:
            frequencies.extend(float(v) for v in m.group(1).split())

    final_geometry = _parse_final_geometry(lines)

    # Last match, not first -- a restarted/reread job reprints this
    # block, so the last one is the one that actually produced
    # final_geometry above.
    charge_mult_matches = _CHARGE_MULT_RE.findall(text)
    charge, multiplicity = (
        (int(charge_mult_matches[-1][0]), int(charge_mult_matches[-1][1]))
        if charge_mult_matches
        else (None, None)
    )

    return LogResult(
        normal_termination=normal_termination,
        stationary_point_found=stationary_point_found,
        scf_energy=scf_energy,
        frequencies=frequencies,
        final_geometry=final_geometry,
        zpe_energy=zpe_energy,
        charge=charge,
        multiplicity=multiplicity,
    )


def parse_irc_log(log_path: str) -> IRCResult:
    """
    Parse a finished IRC .log. Gaussian runs the forward direction
    first (Path Number 1) then the reverse direction (Path Number 2);
    the endpoint of each is the last orientation block printed before
    (forward) / anywhere after (reverse) the first "Path Number: 2"
    marker -- simpler and just as reliable as tracking every
    intermediate corrector-step block, since only the final geometry
    of each direction matters here.
    """
    with open(log_path) as f:
        lines = f.readlines()
    text = "".join(lines)

    normal_termination = "Normal termination of Gaussian" in text
    completed = "Reaction path calculation complete" in text

    path2_start = None
    for i, line in enumerate(lines):
        if _IRC_PATH2_RE.search(line):
            path2_start = i
            break

    if path2_start is not None:
        forward_endpoint = _parse_final_geometry(lines[:path2_start])
        reverse_endpoint = _parse_final_geometry(lines[path2_start:])
    else:
        # Only one direction ran (or the path markers weren't found) --
        # whatever's there is the one endpoint we have.
        forward_endpoint = _parse_final_geometry(lines)
        reverse_endpoint = None

    return IRCResult(
        normal_termination=normal_termination,
        completed=completed,
        forward_endpoint=forward_endpoint,
        reverse_endpoint=reverse_endpoint,
    )
