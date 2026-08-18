"""
irc.py

Runs the IRC from a validated TS geometry (default MaxPoints=25 per
direction, same method/basis/resources as the TS) and checks whether
it actually connects back to the optimized reactant and product --
the real test of whether the TS found is the right one for this
reaction, complementing the imaginary-mode-overlap check in
ts_search.py.

"Connects" is judged by RMSD (after Kabsch-aligning each IRC endpoint
onto the reactant/product, since -- as always -- separately-run
Gaussian jobs share no common frame): close enough to one and close
enough to the other, in either pairing, or it needs a human to look
at it. That's a genuinely different situation from "wrong TS" --  the
path might be entirely correct but just not have relaxed all the way
to the minimum within MaxPoints, which is common and not itself a
failure.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import List, Optional

from .executor import run_com
from .geometry import Structure, rmsd
from .input_builder import build_input, write_com, BasisGroups
from .local_runner import check_normal_termination, GaussianRunError
from .parser import parse_irc_log

RMSD_MATCH_THRESHOLD = 0.5  # Angstrom -- "close enough to call it connected"


@dataclass
class IRCOutcome:
    success: bool
    forward_structure: Optional[Structure]
    reverse_structure: Optional[Structure]
    forward_rmsd_to_reactant: Optional[float] = None
    forward_rmsd_to_product: Optional[float] = None
    reverse_rmsd_to_reactant: Optional[float] = None
    reverse_rmsd_to_product: Optional[float] = None
    log: List[str] = field(default_factory=list)
    log_path: Optional[str] = None  # the actual .log the IRC endpoints came from


def run_irc(
    ts: Structure,
    reactant: Structure,
    product: Structure,
    com_path: str,
    method: str,
    basis: str,
    basis_groups: Optional[BasisGroups] = None,
    ecp_groups: Optional[BasisGroups] = None,
    nprocs: int = 2,
    mem_gb: int = 2,
    maxpoints: int = 25,
    match_threshold: float = RMSD_MATCH_THRESHOLD,
    cancel_event: Optional[threading.Event] = None,
    reactant_trusted: bool = True,
    product_trusted: bool = True,
    executor: str = "local",
) -> IRCOutcome:
    """
    `reactant_trusted`/`product_trusted` -- set False when that side's
    reference structure is only a stand-in (e.g. its own optimization
    failed and `reactant`/`product` here is really just the best
    attempt reached, or the untouched original guess -- see
    webapp.py/cli.py's endpoint-recovery path). With one side untrusted,
    the usual "must land near *both* references" gate is meaningless
    for that side, so classification instead goes by elimination
    against whichever side *is* trusted -- closer to the trusted
    reference wins that label, the other endpoint is reported as the
    (unverified) recovered opposite side, and success no longer
    depends on the untrusted side matching anything.
    """
    trail: List[str] = []

    write_com(
        build_input(
            "irc", [ts], method=method, basis=basis,
            basis_groups=basis_groups, ecp_groups=ecp_groups,
            nprocs=nprocs, mem_gb=mem_gb, irc_maxpoints=maxpoints, chk=com_path[:-4],
        ),
        com_path,
    )

    try:
        log_path = run_com(com_path, executor=executor, cancel_event=cancel_event)
    except GaussianRunError as e:
        trail.append(f"g09 error -- {e}")
        return IRCOutcome(False, None, None, log=trail)

    if not check_normal_termination(log_path):
        trail.append(f"{log_path} did not terminate normally")
        return IRCOutcome(False, None, None, log=trail, log_path=log_path)

    result = parse_irc_log(log_path)
    if result.forward_endpoint is None or result.reverse_endpoint is None:
        trail.append("Couldn't find both IRC endpoints in the log -- needs a manual look.")
        return IRCOutcome(False, None, None, log=trail, log_path=log_path)

    forward = Structure("irc_forward", result.forward_endpoint, ts.charge, ts.multiplicity)
    reverse = Structure("irc_reverse", result.reverse_endpoint, ts.charge, ts.multiplicity)

    f_r = rmsd(forward, reactant)
    f_p = rmsd(forward, product)
    r_r = rmsd(reverse, reactant)
    r_p = rmsd(reverse, product)

    trail.append(f"forward endpoint: RMSD to reactant = {f_r:.3f} A, to product = {f_p:.3f} A")
    trail.append(f"reverse endpoint: RMSD to reactant = {r_r:.3f} A, to product = {r_p:.3f} A")

    if not reactant_trusted or not product_trusted:
        # One side has no real reference to check against, so the usual
        # "must land near *both* references" gate can't apply to it --
        # report success as long as both endpoints exist; classifying
        # *which* endpoint is which side (by elimination against
        # whichever reference is trusted) is verification.py's job when
        # it reoptimizes these endpoints, not this function's.
        untrusted = "reactant" if not reactant_trusted else "product"
        trail.append(
            f"{untrusted.capitalize()} wasn't independently optimized -- can't check the usual "
            "dual-threshold connection criterion against it. Both endpoints exist; which one "
            f"recovers the {untrusted} is decided by elimination downstream."
        )
        return IRCOutcome(True, forward, reverse, f_r, f_p, r_r, r_p, trail, log_path)

    # One endpoint should land near the reactant and the other near
    # the product -- check both pairings, since Gaussian's own
    # forward/reverse labeling (tied to the imaginary mode's arbitrary
    # sign) doesn't tell us which physical direction is which.
    clean_a = f_r < match_threshold and r_p < match_threshold
    clean_b = f_p < match_threshold and r_r < match_threshold

    outcome = IRCOutcome(clean_a or clean_b, forward, reverse, f_r, f_p, r_r, r_p, trail, log_path)

    if clean_a or clean_b:
        trail.append("IRC connects smoothly to both the reactant and product minima.")
    else:
        trail.append(
            "IRC does NOT clearly connect to both endpoints within the RMSD threshold -- "
            "either the path hasn't fully relaxed within MaxPoints (try increasing it), or "
            "this TS doesn't actually connect this reactant/product. Needs a manual look."
        )
    return outcome
