"""
Toy run of the geometry + input-builder pieces on a small, well-known
isomerization (HCN <-> HNC) — this is just to exercise the code path
and show you what the .com files look like. The geometries are
force-field guesses, not vetted starting structures for a real study.

Run from the repo root:  python examples/build_pm6_and_qst2.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gaussbot import from_smiles, build_input, write_com

reactant = from_smiles("C#N", label="HCN reactant")
product = from_smiles("[C-]#[NH+]", label="HNC product")  # standard isocyanide depiction

# QST2 needs the SAME atoms, in the SAME order, in both structures --
# it's one interpolation between two geometries of one molecule, not
# two different molecules. Always check this before submitting.
assert [a[0] for a in reactant.atoms] == [a[0] for a in product.atoms], (
    "reactant/product element order mismatch -- QST2 will reject this"
)

print(f"reactant: {len(reactant.atoms)} atoms, charge {reactant.charge}, mult {reactant.multiplicity}")
print(f"product:  {len(product.atoms)} atoms, charge {product.charge}, mult {product.multiplicity}")

# Step 1: PM6 pre-optimization inputs — run these first, separately.
write_com(build_input("pm6_opt", [reactant]), "jobs/hcn_reactant_pm6opt.com")
write_com(build_input("pm6_opt", [product]), "jobs/hcn_product_pm6opt.com")

# Step 2: PM6 QST2 TS search. In a real run you'd swap in the
# PM6-optimized geometries from step 1 rather than the raw guesses.
qst2_text = build_input("qst2", [reactant, product], method="PM6")
write_com(qst2_text, "jobs/hcn_ts_qst2.com")

print("\n--- jobs/hcn_ts_qst2.com ---\n")
print(qst2_text)
