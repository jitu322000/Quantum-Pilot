"""
gamessbot — an orchestration layer around GAMESS for multireference
electronic structure studies (guess geometry -> RHF -> CIS -> CASSCF ->
XMCQDPT). This round covers RHF and CIS; CASSCF/XMCQDPT and UHF/ROHF
are future work.

This package does NOT do any quantum chemistry itself. It builds .inp
files, hands them to GAMESS (via rungms), and parses the logs that come
back. For the "guess geometry" input path, it reuses gaussbot's own
Gaussian-based geometry optimization pipeline rather than duplicating it.
"""

from .gamess_input import build_rhf_input, build_cis_input

__all__ = [
    "build_rhf_input",
    "build_cis_input",
]
