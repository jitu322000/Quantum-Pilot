"""
xyz_export.py

Optional, best-effort extra: convert a finished Gaussian .log to a
plain .xyz file via Open Babel, for dropping straight into a paper's
supporting information -- geometry.Structure.to_xyz_file() writes an
.xyz from a Structure gaussbot already has in memory, but this is for
handing someone the log's *own* final geometry as a standalone file
without needing to load it back into Python first.
"""

from __future__ import annotations

import subprocess
from typing import Optional


def convert_log_to_xyz(log_path: str, xyz_path: Optional[str] = None, timeout: int = 30) -> Optional[str]:
    """
    `obabel -ig09 <log_path> -oxyz -O <xyz_path>`. `xyz_path` defaults
    to `log_path` with its extension swapped for `.xyz`.

    Optional and best-effort: returns None rather than raising if
    `obabel` isn't on PATH, the conversion fails, or it times out --
    callers should skip this file entirely rather than fail the whole
    study over it. Returns `xyz_path` on success.
    """
    xyz_path = xyz_path or f"{log_path.rsplit('.', 1)[0]}.xyz"
    try:
        result = subprocess.run(
            ["obabel", "-ig09", log_path, "-oxyz", "-O", xyz_path],
            capture_output=True, text=True, timeout=timeout,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    return xyz_path if result.returncode == 0 else None


def convert_log_to_cdxml(log_path: str, cdxml_path: Optional[str] = None, timeout: int = 30) -> Optional[str]:
    """
    `obabel -ig09 <log_path> -ocdxml -O <cdxml_path>`. `cdxml_path`
    defaults to `log_path` with its extension swapped for `.cdxml`.

    Same best-effort contract as convert_log_to_xyz(): returns None
    rather than raising if `obabel` isn't on PATH, the conversion
    fails, or it times out. Returns `cdxml_path` on success.
    """
    cdxml_path = cdxml_path or f"{log_path.rsplit('.', 1)[0]}.cdxml"
    try:
        result = subprocess.run(
            ["obabel", "-ig09", log_path, "-ocdxml", "-O", cdxml_path],
            capture_output=True, text=True, timeout=timeout,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    return cdxml_path if result.returncode == 0 else None
