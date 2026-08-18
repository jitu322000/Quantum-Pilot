"""
intake.py

Turns a reactant and a product -- each an existing Gaussian .com
file, an .xyz/.mol/.sdf/.pdb file, or a plain SMILES string -- into
Structures ready for build_input(). Two front doors: load_structure()
for a known file path (what a non-interactive CLI invocation uses),
and prompt_for_structure()/prompt_for_reaction() for walking a person
through the choice interactively.
"""

from __future__ import annotations

import os
import re
from typing import Optional, Tuple

from .geometry import Structure, from_smiles, from_file


def from_com_file(path: str, label: Optional[str] = None) -> Structure:
    """
    Parse the first geometry block out of an existing Gaussian .com
    file -- useful when the user already has a hand-built or
    previously-run input they want to reuse as a starting geometry.

    Handles the standard layout this package itself writes: link0
    (%...) lines, a route (#...) line, blank line, title, blank line,
    charge/multiplicity, Cartesian coordinates, blank line. Does not
    handle Z-matrices.
    """
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]

    i = 0
    seen_route = False
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped.startswith("%"):
            i += 1
            continue
        if stripped.startswith("#"):
            seen_route = True
            i += 1
            continue
        if seen_route and stripped == "":
            i += 1
            break
        i += 1
    if not seen_route:
        raise ValueError(f"{path} doesn't look like a Gaussian input file (no route line found)")

    while i < len(lines) and lines[i].strip() == "":
        i += 1
    title_line = lines[i].strip()
    i += 1

    while i < len(lines) and lines[i].strip() == "":
        i += 1
    charge_mult = lines[i].split()
    if len(charge_mult) < 2:
        raise ValueError(f"Couldn't find a charge/multiplicity line in {path}")
    charge, multiplicity = int(charge_mult[0]), int(charge_mult[1])
    i += 1

    atom_re = re.compile(r"^\s*([A-Za-z]{1,2})\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)")
    atoms = []
    while i < len(lines) and lines[i].strip() != "":
        m = atom_re.match(lines[i])
        if not m:
            break
        sym, x, y, z = m.groups()
        atoms.append((sym, float(x), float(y), float(z)))
        i += 1

    if not atoms:
        raise ValueError(f"No atom coordinates found in {path}")

    return Structure(label=label or title_line or path, atoms=atoms, charge=charge, multiplicity=multiplicity)


def load_structure(path: str, label: Optional[str] = None) -> Structure:
    """
    Load a Structure straight from a file path, picking the loader by
    extension -- the non-interactive counterpart to
    prompt_for_structure(), for a CLI invocation that already knows
    what file it wants.
    """
    label = label or os.path.splitext(os.path.basename(path))[0]
    ext = os.path.splitext(path)[1].lower()
    if ext == ".com":
        return from_com_file(path, label=label)
    return from_file(path, label=label)


def check_reaction_match(reactant: Structure, product: Structure) -> bool:
    """
    QST2/QST3 need identical atom count and order between reactant
    and product. Prints a warning (or confirmation) and returns
    whether they match.
    """
    r_elems = [a[0] for a in reactant.atoms]
    p_elems = [a[0] for a in product.atoms]
    if r_elems != p_elems:
        print(
            "\nWarning: reactant and product don't have matching atom "
            "count/order -- QST2/QST3 needs the same atoms in the same "
            f"order.\n  reactant: {r_elems}\n  product:  {p_elems}"
        )
        return False
    print(f"\nOK: {len(r_elems)} atoms, matching order in both structures.")
    return True


def prompt_for_structure(role: str) -> Structure:
    """
    Ask the user how they want to supply one geometry (reactant or
    product) and load it. role is used in the prompts, e.g. "reactant".
    """
    print(f"\nHow do you want to provide the {role}?")
    print("  1) Gaussian .com file")
    print("  2) .xyz file")
    print("  3) .mol / .sdf / .pdb file")
    print("  4) SMILES string")
    choice = input("Choice [1-4]: ").strip()

    label = input(f"Label for this structure (default: {role}): ").strip() or role

    if choice == "1":
        path = input(f"Path to {role} .com file: ").strip()
        return from_com_file(path, label=label)
    if choice == "2":
        path = input(f"Path to {role} .xyz file: ").strip()
        return from_file(path, label=label)
    if choice == "3":
        path = input(f"Path to {role} file: ").strip()
        return from_file(path, label=label)
    if choice == "4":
        smiles = input(f"SMILES for {role}: ").strip()
        return from_smiles(smiles, label=label)
    raise ValueError(f"Unrecognized choice: {choice!r}")


def prompt_for_reaction() -> Tuple[Structure, Structure]:
    """
    Prompt for both reactant and product, then flag (not block) a
    mismatch -- QST2/QST3 need identical atom count and order, so
    this is the first thing worth checking before building an input.
    """
    reactant = prompt_for_structure("reactant")
    product = prompt_for_structure("product")
    check_reaction_match(reactant, product)
    return reactant, product


def prompt_for_ts_guess() -> Optional[Structure]:
    """
    Ask if the user has a transition-state guess to provide. Blank ->
    None, meaning the pipeline should generate its own guess by
    interpolating between the optimized reactant and product.
    """
    path = input(
        "\nTransition-state guess file (.com/.xyz/.mol/.sdf/.pdb), "
        "or leave blank to generate one automatically: "
    ).strip()
    if not path:
        return None
    return load_structure(path, label="ts_guess")


def prompt_for_ts_search_inputs() -> Tuple[Optional[Structure], Optional[Structure], Optional[Structure]]:
    """
    Ask what to provide for a dedicated TS search: a TS guess, a
    reactant+product pair, or both together. Returns
    (reactant, product, ts_guess) -- any may be None, but at least a
    ts_guess or both reactant and product is required.

    Unlike the full mechanism study, reactant/product here are used
    exactly as given -- no PM6 pre-opt, no final-level reopt -- so a
    non-.log upload (not already a known-optimized Gaussian output)
    gets a warning that it should already be a good geometry.
    """
    print("\nWhat do you want to provide?")
    print("  1) TS guess only")
    print("  2) Reactant + product only (used as-is -- ideally already-optimized Gaussian .log files)")
    print("  3) Both -- reactant + product (to validate the TS) and a TS guess (as the starting geometry)")
    choice = input("Choice [1-3]: ").strip()

    reactant = product = ts_guess = None
    if choice in ("2", "3"):
        r_path = input("Path to reactant file (.log strongly preferred): ").strip()
        reactant = load_structure(r_path, label="reactant")
        if not r_path.lower().endswith(".log"):
            print(
                "  Warning: not a .log file -- this section doesn't PM6-preopt or reoptimize "
                "the reactant/product, so make sure it's already a good/optimized geometry, or "
                "the TS search may not find the transition state you're looking for."
            )
        p_path = input("Path to product file (.log strongly preferred): ").strip()
        product = load_structure(p_path, label="product")
        if not p_path.lower().endswith(".log"):
            print(
                "  Warning: not a .log file -- this section doesn't PM6-preopt or reoptimize "
                "the reactant/product, so make sure it's already a good/optimized geometry, or "
                "the TS search may not find the transition state you're looking for."
            )
    if choice in ("1", "3"):
        ts_path = input("Path to TS guess file: ").strip()
        ts_guess = load_structure(ts_path, label="ts_guess")

    if ts_guess is None and (reactant is None or product is None):
        raise ValueError("Need a TS guess, or both reactant and product, to run a TS search.")

    return reactant, product, ts_guess
