"""
thermochem.py

Optional, best-effort extra: the thermal correction to Gibbs free
energy (G-corr) for a single molecule, at a user-chosen temperature,
straight from Gaussian's own `freqchk` utility run against that
molecule's checkpoint file -- not something parser.py's regex approach
can get, since it depends on the temperature/pressure/isotope choices
freqchk makes at run time, not just reading a fixed line out of a log.

`freqchk` is interactive, not flag-driven -- there's no command-line
way to hand it a temperature. The prompt sequence below (checkpoint
file, "write Hyperchem files?", temperature, pressure, frequency scale
factor, "use principal isotope masses?") was verified by hand against
a real .chk from this project; skipping the last prompt makes freqchk
loop asking for every atom's isotope mass one at a time, so it's
answered explicitly ("y", the default) rather than left out.
"""

from __future__ import annotations

import re
import subprocess
from typing import Optional

_GCORR_RE = re.compile(r"Thermal correction to Gibbs Free Energy=\s*(-?\d+\.\d+)")

_TO_KELVIN = {
    "K": lambda t: t,
    "C": lambda t: t + 273.15,
    "F": lambda t: (t - 32) * 5.0 / 9.0 + 273.15,
}


def to_kelvin(value: float, unit: str) -> float:
    """`unit` is "K", "C", or "F" -- whatever the temperature field's
    unit selector sent. Raises ValueError for anything else, same as a
    bad energy unit would in energetics.py."""
    try:
        return _TO_KELVIN[unit.upper()](value)
    except KeyError:
        raise ValueError(f"Unknown temperature unit {unit!r}. Choose from: K, C, F")


def run_freqchk(chk_path: str, temperature_kelvin: float, timeout: int = 60) -> Optional[float]:
    """
    Thermal correction to Gibbs Free Energy (Hartree/particle) for the
    structure in `chk_path`, at `temperature_kelvin`, unscaled (no
    empirical frequency scaling -- consistent with the rest of gaussbot
    never applying one) and using default isotope masses.

    This is optional and best-effort: returns None rather than raising
    if `freqchk` isn't on PATH, the checkpoint has no frequency data,
    it times out, or the expected output line just isn't there --
    callers should skip this value entirely rather than fail the whole
    study over it.
    """
    stdin = f"{chk_path}\nn\n{temperature_kelvin}\n1\n1\ny\n"
    try:
        result = subprocess.run(
            ["freqchk"], input=stdin, capture_output=True, text=True, timeout=timeout,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    match = _GCORR_RE.search(result.stdout)
    return float(match.group(1)) if match else None
