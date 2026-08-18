"""
energetics.py

Turns CASSCF/XMCQDPT state energies (Hartree, absolute) into the table
form actually used for analysis: vertical excitation energies in eV,
relative to the ground state (S0), matching standard practice in
multireference spectroscopy work (e.g. S. Rajagopala Reddy's own
CASSCF/XMCQDPT publications report "VEE" this way, not absolute
energies) -- and a LaTeX rendering of that table for papers/reports.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .active_space import ActiveSpaceSuggestion
from .parser import CASSCFResult, XMCQDPTResult

HARTREE_TO_EV = 27.211386245988  # CODATA 2018


@dataclass
class EnergyTableRow:
    state_index: int  # 1-based, as GAMESS numbers it (1 = S0)
    label: str  # "S0", "S1", ...
    casscf_ev: Optional[float]  # vertical excitation energy relative to S0, eV
    xmcqdpt_ev: Optional[float] = None
    oscillator_strength: Optional[float] = None  # relative to S0 (parser.parse_transitn_log()) -- S0's own row is always None


def build_energy_table(
    casscf_result: CASSCFResult,
    xmcqdpt_result: Optional[XMCQDPTResult] = None,
    oscillator_strengths: Optional[Dict[int, float]] = None,
) -> List[EnergyTableRow]:
    """
    Builds one row per CASSCF state, each energy expressed as a vertical
    excitation energy (eV) relative to that calculation's own S0 (state
    index 1) -- so S0 always reads 0.00. XMCQDPT columns are filled in
    only for states present in `xmcqdpt_result.mcqdpt_state_energies`
    (None otherwise, e.g. if XMCQDPT wasn't run). `oscillator_strengths`
    (transitn.run_transitn()'s own output, state_index -> f, relative
    to S0) fills in the oscillator-strength column the same way -- pass
    None (the default) to omit that column entirely, e.g. when the
    RUNTYP=TRANSITN stage wasn't run for this active space/state combo.
    """
    casscf_by_index = dict(casscf_result.state_energies)
    if not casscf_by_index:
        return []
    casscf_reference = casscf_by_index[min(casscf_by_index)]

    xmcqdpt_by_index = dict(xmcqdpt_result.mcqdpt_state_energies) if xmcqdpt_result else {}
    xmcqdpt_reference = xmcqdpt_by_index.get(min(xmcqdpt_by_index)) if xmcqdpt_by_index else None

    rows = []
    for index in sorted(casscf_by_index):
        xmcqdpt_ev = None
        if index in xmcqdpt_by_index and xmcqdpt_reference is not None:
            xmcqdpt_ev = (xmcqdpt_by_index[index] - xmcqdpt_reference) * HARTREE_TO_EV
        rows.append(EnergyTableRow(
            state_index=index, label=f"S{index - 1}",
            casscf_ev=(casscf_by_index[index] - casscf_reference) * HARTREE_TO_EV,
            xmcqdpt_ev=xmcqdpt_ev,
            oscillator_strength=(oscillator_strengths or {}).get(index),
        ))
    return rows


def format_latex_table(
    rows: List[EnergyTableRow],
    nstate: int,
    active_space: ActiveSpaceSuggestion,
    caption: Optional[str] = None,
    label: Optional[str] = None,
) -> str:
    """
    Renders `rows` as a LaTeX table with a caption in the "SA5 (8,8)"
    style (state-averaged over `nstate` states, active space expressed
    as (electrons, orbitals)) -- ready to drop into a paper or report.
    """
    n_electrons = 2 * active_space.ndoc
    n_orbitals = active_space.ndoc + active_space.nval
    if caption is None:
        caption = (
            f"SA{nstate} ({n_electrons},{n_orbitals}) vertical excitation energies "
            "at the CASSCF and XMCQDPT levels of theory."
        )
    if label is None:
        label = f"tab:sa{nstate}-{n_electrons}-{n_orbitals}"

    # The oscillator-strength column is only added when at least one row
    # actually has one (i.e. RUNTYP=TRANSITN was run for this active
    # space/state combo) -- omitted entirely otherwise, rather than a
    # column of nothing but "--".
    has_f = any(row.oscillator_strength is not None for row in rows)

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption{{{caption}}}",
        f"\\label{{{label}}}",
        r"\begin{tabular}{lccc}" if has_f else r"\begin{tabular}{lcc}",
        r"\hline",
        (r"State & CASSCF (eV) & XMCQDPT (eV) & $f$ \\" if has_f else r"State & CASSCF (eV) & XMCQDPT (eV) \\"),
        r"\hline",
    ]
    for row in rows:
        casscf_str = f"{row.casscf_ev:.2f}" if row.casscf_ev is not None else "-"
        xmcqdpt_str = f"{row.xmcqdpt_ev:.2f}" if row.xmcqdpt_ev is not None else "-"
        state_tex = f"S$_{{{row.state_index - 1}}}$"
        row_tex = f"{state_tex} & {casscf_str} & {xmcqdpt_str}"
        if has_f:
            f_str = f"{row.oscillator_strength:.4f}" if row.oscillator_strength is not None else "-"
            row_tex += f" & {f_str}"
        lines.append(row_tex + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)
