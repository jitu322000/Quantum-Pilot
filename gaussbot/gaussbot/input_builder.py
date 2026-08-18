"""
input_builder.py

Turns one or more Structure objects plus a job type into a Gaussian
.com file. This is the only place that knows Gaussian input syntax —
everything else in the pipeline just passes Structures around.
"""

from __future__ import annotations

import os
import re
import warnings
from typing import List, Optional, Tuple

from .geometry import Structure

# A "basis group" is (element_symbols, name) -- e.g. (["H","C","O"], "6-31G(d,p)").
BasisGroups = List[Tuple[List[str], str]]

# {mb} is filled in with "method" or "method/basis" (PM6 has no basis).
# NoSymm on every route: (1) an optimization constrained to a detected
# point group can get trapped at a symmetric saddle point instead of
# reaching a lower-symmetry minimum -- exactly the failure mode the
# imaginary-frequency repair loop exists to fix, so let the optimizer
# roam free instead of fighting symmetry constraints; (2) it makes
# Gaussian report geometry as "Input orientation" throughout instead
# of reorienting to a symmetry-adapted frame, which keeps every
# geometry/displacement-vector we parse back in the same Cartesian
# frame we sent in.
JOB_ROUTES = {
    "pm6_opt": "{mb} Opt NoSymm",
    "pm6_opt_freq": "{mb} Opt Freq NoSymm",
    "qst2": "{mb} Opt=QST2 Freq NoSymm",
    "qst3": "{mb} Opt=QST3 Freq NoSymm",
    "ts_highlevel": "{mb} Opt=(TS,CalcFC,NoEigenTest) Freq NoSymm",
    # Brute-force TS fallback: CalcAll recomputes the full Hessian at
    # every step instead of just the first (CalcFC) -- far more likely
    # to actually locate a genuine saddle point when CalcFC can't, at
    # a much higher cost per step.
    "ts_calcall": "{mb} Opt=(TS,CalcAll,NoEigenTest) Freq NoSymm",
    "opt_highlevel": "{mb} Opt Freq NoSymm",
    # No "Both": not valid IRC syntax on this G09 revision (RevC.01) --
    # syntax-errors the whole route. CalcFC alone already runs forward
    # and reverse by default (confirmed against a real run).
    "irc": "{mb} IRC=(CalcFC) NoSymm",
}

STRUCTURE_COUNTS = {
    "pm6_opt": 1,
    "pm6_opt_freq": 1,
    "qst2": 2,   # reactant, product
    "qst3": 3,   # reactant, product, guessed TS
    "ts_highlevel": 1,
    "ts_calcall": 1,
    "opt_highlevel": 1,
    "irc": 1,
}


def _inject_keyword_option(route: str, keyword: str, option: str) -> str:
    """
    Fold an extra option (e.g. "MaxCycles=100" into Opt, "MaxPoints=25"
    into IRC) into a route keyword, whatever form it's already in --
    bare "Opt", "Opt=QST2", "Opt=(...)" -- so retries don't end up
    with two conflicting copies of the same keyword on the route line.
    """
    paren_re = re.compile(rf"\b{keyword}=\(([^)]*)\)")
    value_re = re.compile(rf"\b{keyword}=([A-Za-z0-9]+)\b")
    bare_re = re.compile(rf"\b{keyword}\b(?!=)")

    m = paren_re.search(route)
    if m:
        return paren_re.sub(f"{keyword}=({m.group(1)},{option})", route, count=1)
    m = value_re.search(route)
    if m:
        return value_re.sub(f"{keyword}=({m.group(1)},{option})", route, count=1)
    if bare_re.search(route):
        return bare_re.sub(f"{keyword}=({option})", route, count=1)
    raise ValueError(f"Route {route!r} has no {keyword} keyword to add {option!r} to")


def build_input(
    job_type: str,
    structures: List[Structure],
    method: str = "PM6",
    basis: str = "",
    basis_groups: Optional[BasisGroups] = None,
    ecp_groups: Optional[BasisGroups] = None,
    nprocs: int = 2,
    mem_gb: int = 2,
    chk: Optional[str] = None,
    extra_route: str = "",
    opt_maxcycles: Optional[int] = None,
    irc_maxpoints: Optional[int] = None,
) -> str:
    """
    Build the full text of a Gaussian .com file.

    job_type: one of JOB_ROUTES, e.g. "pm6_opt", "qst2", "ts_highlevel", "irc"
    structures: the geometries the job type expects, in order
                (qst2 -> [reactant, product]; qst3 -> [reactant, product, guess_ts])
    method/basis: e.g. method="PM6" (no basis needed), or
                  method="B3LYP", basis="6-31G(d)"
    basis_groups/ecp_groups: a mixed (GenECP) basis instead of a single
                  `basis` string -- e.g. light atoms at 6-31G(d,p),
                  a transition metal at LanL2DZ with its ECP. Each is a
                  list of (element_symbols, name) pairs. `basis` is
                  ignored when `basis_groups` is given; ecp_groups may
                  be empty (a per-element basis with no ECP -> "gen"
                  instead of "genecp", still needs Pseudo=Read since
                  the basis itself is still being read from the input).
    extra_route: appended verbatim, e.g. "Geom=AllCheck Guess=Read" for a
                 restart, or an empirical dispersion keyword
    """
    if job_type not in JOB_ROUTES:
        raise ValueError(f"Unknown job_type {job_type!r}. Choose from: {sorted(JOB_ROUTES)}")

    expected = STRUCTURE_COUNTS[job_type]
    if len(structures) != expected:
        raise ValueError(f"{job_type} needs {expected} structure(s), got {len(structures)}")

    charges = {s.charge for s in structures}
    if len(structures) > 1 and len(charges) > 1:
        warnings.warn(
            f"Structures for a {job_type} job have different formal charges {charges} — "
            "double check the reactant/product/TS are really the same overall species."
        )

    if basis_groups is not None:
        method_basis = f"{method}/genecp" if ecp_groups else f"{method}/gen"
    else:
        method_basis = f"{method}/{basis}" if basis else method
    route = JOB_ROUTES[job_type].format(mb=method_basis)
    if opt_maxcycles is not None:
        route = _inject_keyword_option(route, "Opt", f"MaxCycles={opt_maxcycles}")
    if irc_maxpoints is not None:
        route = _inject_keyword_option(route, "IRC", f"MaxPoints={irc_maxpoints}")
    if basis_groups is not None:
        route = f"{route} Pseudo=Read"
    if extra_route:
        route = f"{route} {extra_route}"

    chk_name = chk or job_type

    lines = [
        f"%chk={chk_name}.chk",
        f"%mem={mem_gb}GB",
        f"%nprocshared={nprocs}",
        f"#p {route}",
        "",
    ]
    for struct in structures:
        lines.append(struct.label)
        lines.append("")
        lines.append(f"{struct.charge} {struct.multiplicity}")
        lines.append(struct.to_xyz_block())
        lines.append("")

    if basis_groups is not None:
        for elements, name in basis_groups:
            lines.append(" ".join(elements) + " 0")
            lines.append(name)
            lines.append("****")
        if ecp_groups:
            lines.append("")
            for elements, name in ecp_groups:
                lines.append(" ".join(elements) + " 0")
                lines.append(name)
        lines.append("")
        lines.append("")

    return "\n".join(lines) + "\n"


def write_com(content: str, path: str) -> str:
    """Write built input text to disk. Returns the path, for chaining."""
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w") as f:
        f.write(content)
    return path
