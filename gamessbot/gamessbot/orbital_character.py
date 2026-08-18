"""
orbital_character.py

Reads which atom/basis-function each MO in a suggested active space is
mostly made of, by scanning the RHF log's own EIGENVECTORS printout for
the largest-magnitude coefficient in each MO's column -- lets you see
at a glance what an active orbital chemically IS (e.g. "O1 2p"), not
just its index number, per your request. Meant to be shown collapsed
(a dropdown), not directly on the results screen.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List

_HEADER_RE = re.compile(r"^(?:\s+\d+)+\s*$")
_ROW_RE = re.compile(r"^\s*\d+\s+([A-Za-z]+)\s+(\d+)\s+(\S+)\s+(.+)$")


@dataclass
class OrbitalCharacter:
    atom_symbol: str
    atom_number: int
    shell: str  # e.g. "S", "X", "Y", "Z", "XX", ... (GAMESS's own basis-function shell labels)
    coefficient: float

    @property
    def label(self) -> str:
        return f"{self.atom_symbol}{self.atom_number} {self.shell}"


def parse_orbital_character(log_path: str, mo_indices: List[int]) -> Dict[int, OrbitalCharacter]:
    """
    For each MO in `mo_indices`, finds the basis function with the
    largest |coefficient| in that MO's own printed column of the log's
    "EIGENVECTORS" table -- confirmed format against a real RHF log:

                      1          2          3          4          5
                  -20.2363    -1.2687    -0.6270    -0.4465    -0.3893
                     A          A          A          A          A
        1  O  1  S    0.994102  -0.232724 ...
        ...

    Returns only the MOs actually found (missing ones -- e.g. requested
    indices beyond the printed range -- are simply absent).
    """
    with open(log_path) as f:
        lines = f.read().splitlines()

    start = next((i for i, line in enumerate(lines) if line.strip() == "EIGENVECTORS"), None)
    if start is None:
        return {}

    remaining = set(mo_indices)
    result: Dict[int, OrbitalCharacter] = {}
    i = start
    while i < len(lines) and remaining:
        line = lines[i]
        if _HEADER_RE.match(line):
            mo_numbers = [int(x) for x in line.split()]
            # next line: orbital energies; then: symmetry labels; then basis-function rows until blank
            j = i + 3
            rows = []
            while j < len(lines) and lines[j].strip():
                m = _ROW_RE.match(lines[j])
                if m:
                    coeffs = [float(x) for x in m.group(4).split()]
                    rows.append((m.group(1), int(m.group(2)), m.group(3), coeffs))
                j += 1
            for col, mo in enumerate(mo_numbers):
                if mo in remaining and rows and all(col < len(r[3]) for r in rows):
                    best = max(rows, key=lambda r: abs(r[3][col]))
                    result[mo] = OrbitalCharacter(
                        atom_symbol=best[0], atom_number=best[1], shell=best[2], coefficient=best[3][col],
                    )
                    remaining.discard(mo)
            i = j
        else:
            i += 1
    return result
