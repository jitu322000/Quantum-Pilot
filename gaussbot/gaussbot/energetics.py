"""
energetics.py

Reports the barrier and reaction energy for a completed reactant/TS/
product study, relative to the reactant (= 0), in whatever unit the
user wants -- the last step once a TS has been found and verified.

Uses the ZPE-corrected energy (LogResult.energy: "Sum of electronic
and zero-point Energies" when the log has it, the standard 0 K
reference for a barrier height) rather than the bare electronic SCF
energy.
"""

from __future__ import annotations

from dataclasses import dataclass

import questionary

from .parser import LogResult

# Hartree -> unit.
HARTREE_TO = {
    "kcal/mol": 627.5094740631,
    "kJ/mol": 2625.499639,
    "eV": 27.211386245988,
    "Hartree": 1.0,
}


def prompt_for_energy_unit() -> str:
    return questionary.select("Report energies in:", choices=list(HARTREE_TO)).ask() or "kcal/mol"


@dataclass
class EnergyReport:
    unit: str
    reactant: float  # always 0.0, by definition -- the reference point
    ts: float
    product: float


def build_energy_report(
    reactant_result: LogResult, ts_result: LogResult, product_result: LogResult, unit: str
) -> EnergyReport:
    if unit not in HARTREE_TO:
        raise ValueError(f"Unknown unit {unit!r}. Choose from: {sorted(HARTREE_TO)}")
    factor = HARTREE_TO[unit]
    e_r = reactant_result.energy
    e_ts = ts_result.energy
    e_p = product_result.energy
    if e_r is None or e_ts is None or e_p is None:
        raise ValueError("One of reactant/TS/product has no usable energy in its log -- can't build a report")
    return EnergyReport(unit=unit, reactant=0.0, ts=(e_ts - e_r) * factor, product=(e_p - e_r) * factor)


def format_energy_report(report: EnergyReport) -> str:
    lines = [
        f"Relative energies ({report.unit}), reactant = 0:",
        f"  Reactant:        {report.reactant:>10.2f}",
        f"  TS:              {report.ts:>10.2f}   (forward barrier)",
        f"  Product:         {report.product:>10.2f}   "
        f"({'exothermic' if report.product < 0 else 'endothermic'})",
        f"  Reverse barrier: {report.ts - report.product:>10.2f}   (TS - product)",
    ]
    return "\n".join(lines)
