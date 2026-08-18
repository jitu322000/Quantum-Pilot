"""
intake.py

The three ways a study can start, all funneled down to a single $DATA
block:

  1. An already-optimized Gaussian .log -> data_block_from_optimized_log()
  2. An existing GAMESS input/punch file -> data_block_from_gamess_file()
  3. A guess geometry (SMILES / uploaded file) -> optimize_and_convert(),
     which reuses gaussbot's own two-stage Gaussian pipeline -- PM6
     pre-optimization, THEN a real optimization at a user-chosen final
     level (repair_and_optimize) -- and converts the FINAL-level log,
     not the PM6 log. A PM6 geometry is only a cheap starting shape;
     GAMESS needs the actual final-level minimum, same as gaussbot's own
     reactant/product/TS stages never stop at PM6 either.
"""

from __future__ import annotations

import os
import threading
from typing import Optional, Tuple

from gaussbot.geometry import Structure, from_file, from_smiles
from gaussbot.level_select import prompt_for_method_basis
from gaussbot.pipeline import OptOutcome, preopt_with_escalation, repair_and_optimize

from .gamess_input import data_block_from_gamess_inp, data_block_from_gaussian_log


def data_block_from_optimized_log(log_path: str, title: Optional[str] = None) -> str:
    """Source 1: the user already has an optimized Gaussian .log."""
    return data_block_from_gaussian_log(log_path, title=title)


def data_block_from_gamess_file(inp_path: str, title: Optional[str] = None) -> str:
    """Source 2: the user already has a GAMESS input/punch file for
    this molecule."""
    return data_block_from_gamess_inp(inp_path, title=title)


def optimize_and_convert(
    structure: Structure,
    out_dir: str,
    tag: str,
    method: str,
    basis: str,
    nprocs: int = 2,
    mem_gb: int = 2,
    cancel_event: Optional[threading.Event] = None,
    executor: str = "local",
) -> Tuple[str, OptOutcome, OptOutcome]:
    """
    Source 3: a guess geometry. Runs gaussbot's full two-stage pipeline:
    PM6 pre-optimization (cheap starting shape), then a real optimization
    at `method`/`basis` (repair_and_optimize) -- and converts the
    FINAL-level log to a $DATA block via data_block_from_optimized_log().

    Returns (data_block, preopt_outcome, final_outcome). Raises
    ValueError if either stage never reached a usable stationary point --
    `outcome.log` on the failing stage has the human-readable trail of
    what was tried.
    """
    pre = preopt_with_escalation(
        structure, out_dir, tag, nprocs=nprocs, mem_gb=mem_gb,
        cancel_event=cancel_event, executor=executor,
    )
    if not pre.success:
        raise ValueError(
            f"Could not reach a PM6 pre-optimized geometry: {'; '.join(pre.log)}"
        )

    final = repair_and_optimize(
        pre.structure, os.path.join(out_dir, f"{tag}_final.com"), method=method, basis=basis,
        nprocs=nprocs, mem_gb=mem_gb, cancel_event=cancel_event, executor=executor,
    )
    if not final.success:
        raise ValueError(
            f"Could not reach a final-level ({method}/{basis or 'no basis'}) minimum "
            f"for GAMESS input: {'; '.join(final.log)}"
        )

    data_block = data_block_from_optimized_log(final.log_path, title=tag)
    return data_block, pre, final


def prompt_for_geometry_source(
    out_dir: str,
    tag: str,
    nprocs: int = 2,
    mem_gb: int = 2,
    executor: str = "local",
) -> str:
    """
    Ask which of the three sources to use and load it -- the
    interactive counterpart to calling data_block_from_optimized_log()/
    data_block_from_gamess_file()/optimize_and_convert() directly for a
    scripted invocation. Returns the $DATA block either way.
    """
    print("\nHow do you want to provide the geometry?")
    print("  1) Already-optimized Gaussian .log file")
    print("  2) Existing GAMESS input/punch file (has its own $DATA)")
    print("  3) Guess geometry (SMILES or file) -- optimized with Gaussian first")
    choice = input("Choice [1-3]: ").strip()

    if choice == "1":
        path = input("Path to the optimized Gaussian .log: ").strip()
        return data_block_from_optimized_log(path, title=tag)
    if choice == "2":
        path = input("Path to the GAMESS input/punch file: ").strip()
        return data_block_from_gamess_file(path, title=tag)
    if choice == "3":
        print("  a) SMILES string")
        print("  b) File (.xyz/.mol/.sdf/.pdb)")
        sub = input("  Choice [a/b]: ").strip().lower()
        if sub == "a":
            smiles = input("SMILES: ").strip()
            structure = from_smiles(smiles, label=tag)
        else:
            path = input("Path to the geometry file: ").strip()
            structure = from_file(path, label=tag)
        print("\nPick the level for the final-level Gaussian optimization (after a PM6 pre-opt):")
        level = prompt_for_method_basis()
        print("\nOptimizing the guess geometry with Gaussian before handing off to GAMESS:")
        data_block, pre, final = optimize_and_convert(
            structure, out_dir, tag, method=level.method, basis=level.basis,
            nprocs=nprocs, mem_gb=mem_gb, executor=executor,
        )
        for line in pre.log:
            print(f"  {line}")
        for line in final.log:
            print(f"  {line}")
        return data_block
    raise ValueError(f"Unrecognized choice: {choice!r}")
