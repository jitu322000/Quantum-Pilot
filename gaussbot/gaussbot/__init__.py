"""
gaussbot — a small orchestration layer around Gaussian for reaction
mechanism studies (guess geometry -> PM6 pre-opt -> QST2/QST3 TS search
-> higher-level refinement -> IRC -> barrier extraction).

This package does NOT do any quantum chemistry itself. It builds .com
files, hands them to Gaussian via PBS, and parses the logs that come
back. Gaussian does the actual work.
"""

from .geometry import Structure, from_smiles, from_file, from_pubchem, from_gaussian_log
from .input_builder import build_input, write_com

__all__ = [
    "Structure",
    "from_smiles",
    "from_file",
    "from_pubchem",
    "from_gaussian_log",
    "build_input",
    "write_com",
]
