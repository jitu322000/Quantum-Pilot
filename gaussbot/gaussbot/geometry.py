"""
geometry.py

Turns a guess geometry — a SMILES string, an uploaded structure file
(.mol/.sdf, .xyz, .pdb, a ChemDraw .cdx/.cdxml sketch, or an already-
optimized Gaussian .log), or a PubChem lookup — into a Structure: a
plain list of (element, x, y, z) atoms plus charge and spin
multiplicity, ready to drop into a Gaussian .com file.

RDKit does the embedding/optimization. Nothing here talks to Gaussian.
"""

from __future__ import annotations

import dataclasses
import math
import os
import random
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


Atom = Tuple[str, float, float, float]


@dataclass
class Structure:
    """A single molecular geometry ready for a Gaussian input block."""

    label: str
    atoms: List[Atom]
    charge: int = 0
    multiplicity: int = 1

    def to_xyz_block(self) -> str:
        """Cartesian block in the format Gaussian expects (no header)."""
        lines = [f"{sym:<2} {x:>14.8f} {y:>14.8f} {z:>14.8f}" for sym, x, y, z in self.atoms]
        return "\n".join(lines)

    def to_xyz_file(self, path: str) -> None:
        """Write a standalone .xyz file (useful for visual sanity checks)."""
        with open(path, "w") as f:
            f.write(f"{len(self.atoms)}\n{self.label}\n")
            for sym, x, y, z in self.atoms:
                f.write(f"{sym} {x:.8f} {y:.8f} {z:.8f}\n")

    def formula(self) -> str:
        """Molecular formula in Hill notation (C first, H second, then
        the rest alphabetically) -- e.g. "C2H4O". Used to sanity-check
        that an alternate structure (a PubChem hit, an LLM-suggested
        SMILES) is actually the same molecule before offering it."""
        counts: dict = {}
        for sym, _, _, _ in self.atoms:
            counts[sym] = counts.get(sym, 0) + 1
        parts = []
        if "C" in counts:
            parts.append(("C", counts.pop("C")))
            if "H" in counts:
                parts.append(("H", counts.pop("H")))
        for sym in sorted(counts):
            parts.append((sym, counts[sym]))
        return "".join(f"{sym}{n if n > 1 else ''}" for sym, n in parts)


def displace_structure(
    structure: Structure,
    vectors: List[Tuple[float, float, float]],
    scale: float = 0.3,
) -> Structure:
    """
    Return a copy of `structure` with each atom nudged along the given
    per-atom displacement vectors -- typically an imaginary vibrational
    mode read back from a Gaussian log -- scaled by `scale` (mode
    vectors are ~unit-normalized, so `scale` is roughly the step in
    Angstroms). Used to kick a geometry that optimized to a saddle
    point off of it before reoptimizing.
    """
    if len(vectors) != len(structure.atoms):
        raise ValueError(
            f"Got {len(vectors)} displacement vectors for {len(structure.atoms)} atoms"
        )
    new_atoms = [
        (sym, x + dx * scale, y + dy * scale, z + dz * scale)
        for (sym, x, y, z), (dx, dy, dz) in zip(structure.atoms, vectors)
    ]
    return dataclasses.replace(structure, atoms=new_atoms)


def jitter_structure(
    structure: Structure,
    scale: float = 0.02,
    rng: Optional[random.Random] = None,
) -> Structure:
    """
    Return a copy of `structure` with each atom nudged by a small
    random displacement (uniform in [-scale, scale] Angstrom per
    axis). Used to kick a stalled optimization retry off an exact
    repeat of the same geometry -- e.g. continuing a non-converged
    optimization from where it left off, where the raw last geometry
    alone can land right back in the same oscillation.
    """
    rng = rng or random
    new_atoms = [
        (sym, x + rng.uniform(-scale, scale), y + rng.uniform(-scale, scale), z + rng.uniform(-scale, scale))
        for sym, x, y, z in structure.atoms
    ]
    return dataclasses.replace(structure, atoms=new_atoms)


def align_to(mobile: Structure, reference: Structure) -> Structure:
    """
    Rigidly rotate+translate `mobile` (Kabsch, least-squares fit) onto
    `reference` -- same atom order required. Every Gaussian job we run
    is its own independent calculation and (even with NoSymm) has no
    idea what frame any *other* job's geometry is sitting in, so
    comparing or interpolating between two separately-optimized
    structures (e.g. a reactant and a product) isn't meaningful until
    they're brought into a common frame first.
    """
    if len(mobile.atoms) != len(reference.atoms):
        raise ValueError(
            f"Got {len(mobile.atoms)} atoms for mobile, {len(reference.atoms)} for reference"
        )
    mobile_elems = [a[0] for a in mobile.atoms]
    reference_elems = [a[0] for a in reference.atoms]
    if mobile_elems != reference_elems:
        raise ValueError(f"atom order mismatch:\n  mobile:    {mobile_elems}\n  reference: {reference_elems}")

    m = np.array([[x, y, z] for _, x, y, z in mobile.atoms])
    r = np.array([[x, y, z] for _, x, y, z in reference.atoms])
    m_centroid = m.mean(axis=0)
    r_centroid = r.mean(axis=0)
    p = m - m_centroid
    q = r - r_centroid

    u, _, vt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    aligned = (rotation @ p.T).T + r_centroid

    new_atoms = [(sym, *coord) for sym, coord in zip(mobile_elems, aligned.tolist())]
    return dataclasses.replace(mobile, atoms=new_atoms)


def rmsd(mobile: Structure, reference: Structure) -> float:
    """
    Root-mean-square deviation between `mobile` and `reference`, after
    Kabsch-aligning `mobile` onto `reference` (see align_to) -- the
    standard "how close are these two independently-obtained geometries,
    really" check used throughout this pipeline (matching an IRC
    endpoint to the optimized reactant/product, classifying which side
    of a TS-mode distortion a reoptimized structure landed on, etc.).
    """
    aligned = align_to(mobile, reference=reference)
    n = len(reference.atoms)
    total = sum(
        (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
        for (_, ax, ay, az), (_, bx, by, bz) in zip(aligned.atoms, reference.atoms)
    )
    return math.sqrt(total / n)


def interpolate_structures(
    reactant: Structure,
    product: Structure,
    fraction: float = 0.5,
    label: str = "ts_guess",
) -> Structure:
    """
    A crude "linear synchronous transit" TS guess: align the product
    onto the reactant (see align_to) and interpolate each atom's
    Cartesian position `fraction` of the way from reactant to the
    aligned product. Reactant and product must have the same atoms in
    the same order (same requirement as QST2/QST3).
    """
    aligned_product = align_to(product, reference=reactant)
    r_coords = np.array([[x, y, z] for _, x, y, z in reactant.atoms])
    p_coords = np.array([[x, y, z] for _, x, y, z in aligned_product.atoms])
    guess_coords = r_coords + (p_coords - r_coords) * fraction

    elems = [a[0] for a in reactant.atoms]
    atoms = [(sym, *coord) for sym, coord in zip(elems, guess_coords.tolist())]
    return Structure(label=label, atoms=atoms, charge=reactant.charge, multiplicity=reactant.multiplicity)


def _mol_to_structure(mol: Chem.Mol, label: str, multiplicity: int, conf_id: int = 0) -> Structure:
    conf = mol.GetConformer(conf_id)
    atoms: List[Atom] = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append((atom.GetSymbol(), pos.x, pos.y, pos.z))
    charge = Chem.GetFormalCharge(mol)
    return Structure(label=label, atoms=atoms, charge=charge, multiplicity=multiplicity)


def _embed_3d(mol: Chem.Mol, n_confs: int = 20, seed: int = 42, context: str = "") -> Tuple[Chem.Mol, int]:
    """
    Generate several conformers with ETKDGv3, do a quick MMFF94 cleanup
    on each, and return (mol, best_conf_id) for the lowest-energy one.
    This is a force-field guess, not a QM-quality geometry — that's
    what the downstream PM6 optimization is for. Shared by from_smiles()
    and from_chemdraw() -- anywhere the only geometry information
    available going in is a 2D/topological one (a SMILES string, a flat
    ChemDraw sketch) rather than real 3D coordinates. `context` is just
    for the error message (e.g. the SMILES string or file path), so a
    failure says what it was trying to embed.
    """
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if not conf_ids:
        raise ValueError(f"RDKit failed to embed any conformer{f' for {context!r}' if context else ''}")

    energies = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    # energies: list of (not_converged_flag, energy) per conformer
    best_conf_id = min(conf_ids, key=lambda cid: energies[cid][1])
    return mol, best_conf_id


def from_smiles(
    smiles: str,
    label: str = "molecule",
    multiplicity: int = 1,
    n_confs: int = 20,
    seed: int = 42,
) -> Structure:
    """
    Embed a 3D geometry from a SMILES string -- see _embed_3d() for how.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    mol, best_conf_id = _embed_3d(mol, n_confs=n_confs, seed=seed, context=smiles)
    return _mol_to_structure(mol, label=label, multiplicity=multiplicity, conf_id=best_conf_id)


def from_file(
    path: str,
    label: Optional[str] = None,
    multiplicity: int = 1,
    charge_override: Optional[int] = None,
) -> Structure:
    """
    Load a geometry from a file — .mol/.sdf (e.g. a ChemDraw export),
    .pdb, or .xyz.

    .xyz files carry no bond/charge information, so RDKit's bond
    perception (rdDetermineBonds) is used to infer formal charge;
    pass charge_override if you already know it (e.g. a charged TS).
    """
    label = label or os.path.splitext(os.path.basename(path))[0]
    ext = os.path.splitext(path)[1].lower()

    if ext in (".mol", ".sdf"):
        mol = Chem.MolFromMolFile(path, removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit could not parse {path}")
        return _mol_to_structure(mol, label=label, multiplicity=multiplicity)

    if ext == ".pdb":
        mol = Chem.MolFromPDBFile(path, removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit could not parse {path}")
        return _mol_to_structure(mol, label=label, multiplicity=multiplicity)

    if ext == ".xyz":
        raw = Chem.MolFromXYZFile(path)
        if raw is None:
            raise ValueError(f"RDKit could not parse {path}")
        mol = Chem.Mol(raw)
        try:
            from rdkit.Chem import rdDetermineBonds

            rdDetermineBonds.DetermineBonds(mol, charge=charge_override or 0)
        except Exception:
            pass  # fall back to no connectivity; coordinates are still fine
        struct = _mol_to_structure(mol, label=label, multiplicity=multiplicity)
        if charge_override is not None:
            struct.charge = charge_override
        return struct

    if ext in (".cdx", ".cdxml"):
        return from_chemdraw(path, label=label, multiplicity=multiplicity)

    if ext == ".log":
        return from_gaussian_log(path, label=label)

    raise ValueError(f"Unsupported file type: {ext} (expected .mol, .sdf, .pdb, .xyz, .cdx, .cdxml, or .log)")


def from_chemdraw(
    path: str,
    label: Optional[str] = None,
    multiplicity: int = 1,
    n_confs: int = 20,
    seed: int = 42,
) -> Structure:
    """
    Load a ChemDraw file (.cdx binary or .cdxml XML) and embed it in
    3D. ChemDraw files encode a flat 2D drawing, not a real starting
    geometry -- confirmed by hand (export a real structure to .cdxml,
    read it back: the coordinates come back as flat drawing-canvas
    positions, z=0 on every atom, non-physical scale) -- so this reads
    only the atom/bond graph and discards the 2D coordinates, then
    re-embeds with the same RDKit pipeline from_smiles() uses
    (_embed_3d(): ETKDGv3 + MMFF cleanup, lowest-energy conformer kept)
    rather than trusting ChemDraw's 2D layout as if it were 3D.

    Needs Open Babel (`obabel` on PATH) -- RDKit doesn't read ChemDraw's
    own formats directly, so this first converts to an intermediate
    .mol both Open Babel and RDKit understand.
    """
    import subprocess
    import tempfile

    label = label or os.path.splitext(os.path.basename(path))[0]
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    if ext not in ("cdx", "cdxml"):
        raise ValueError(f"Not a ChemDraw file: {path}")

    with tempfile.NamedTemporaryFile(suffix=".mol", delete=False) as tmp:
        tmp_mol_path = tmp.name
    try:
        try:
            result = subprocess.run(
                ["obabel", f"-i{ext}", path, "-omol", "-O", tmp_mol_path],
                capture_output=True, text=True, timeout=30,
            )
        except FileNotFoundError as e:
            raise ValueError(
                "Open Babel ('obabel') isn't on PATH -- needed to read ChemDraw files."
            ) from e
        if result.returncode != 0:
            raise ValueError(f"Open Babel couldn't convert {path}: {result.stderr.strip()}")

        mol = Chem.MolFromMolFile(tmp_mol_path, removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit couldn't parse the structure converted from {path}")
    finally:
        if os.path.exists(tmp_mol_path):
            os.unlink(tmp_mol_path)

    mol = Chem.AddHs(mol)
    mol, best_conf_id = _embed_3d(mol, n_confs=n_confs, seed=seed, context=path)
    return _mol_to_structure(mol, label=label, multiplicity=multiplicity, conf_id=best_conf_id)


def from_gaussian_log(path: str, label: Optional[str] = None) -> Structure:
    """
    Load a Structure straight from an already-finished Gaussian .log --
    for trusting a geometry someone else already optimized (e.g. as a
    reactant/product for a TS search) without re-running PM6 on it.

    Requires the log to actually be a converged stationary point
    (normal termination + "Stationary point found") and to have a
    parseable charge/multiplicity line -- raises clearly otherwise,
    since the entire point of this loader is that the caller can trust
    the geometry it hands back without independently reoptimizing it.
    """
    from .parser import parse_log

    label = label or os.path.splitext(os.path.basename(path))[0]
    result = parse_log(path)

    if not result.normal_termination or not result.stationary_point_found:
        raise ValueError(
            f"{path} isn't a converged Gaussian optimization (normal_termination="
            f"{result.normal_termination}, stationary_point_found={result.stationary_point_found}) "
            "-- from_gaussian_log() only accepts an already-optimized log, since it's meant to be "
            "trusted as-is rather than reoptimized."
        )
    if result.final_geometry is None:
        raise ValueError(f"Couldn't find a final geometry in {path}")
    if result.charge is None or result.multiplicity is None:
        raise ValueError(f"Couldn't find a 'Charge = ... Multiplicity = ...' line in {path}")

    return Structure(label=label, atoms=result.final_geometry, charge=result.charge, multiplicity=result.multiplicity)


def from_pubchem(
    query: str,
    namespace: str = "name",
    label: Optional[str] = None,
    multiplicity: int = 1,
) -> Structure:
    """
    Fetch a geometry from PubChem instead of building one locally —
    useful when a decent starting structure already exists online.

    Tries a 3D SDF record first; falls back to embedding the 2D
    structure with RDKit (see from_smiles) if PubChem has no 3D
    conformer for the compound.

    Requires network access to pubchem.ncbi.nlm.nih.gov.
    """
    import pubchempy as pcp

    label = label or query

    try:
        compounds = pcp.get_compounds(query, namespace, record_type="3d")
        if compounds:
            sdf = compounds[0].to_dict(properties=["record"])
            # pubchempy doesn't expose 3D SDF directly in all versions;
            # simplest robust path is to re-fetch as an SDF file via CID.
            cid = compounds[0].cid
            sdf_text = pcp.get_sdf(cid, "cid", record_type="3d")
            mol = Chem.MolFromMolBlock(sdf_text, removeHs=False)
            if mol is not None:
                return _mol_to_structure(mol, label=label, multiplicity=multiplicity)
    except Exception:
        pass  # fall through to 2D -> embed

    compounds = pcp.get_compounds(query, namespace)
    if not compounds:
        raise ValueError(f"PubChem returned no compounds for {query!r} ({namespace})")
    smiles = compounds[0].isomeric_smiles or compounds[0].canonical_smiles
    return from_smiles(smiles, label=label, multiplicity=multiplicity)
