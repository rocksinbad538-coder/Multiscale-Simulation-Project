#!/usr/bin/env python3
from pathlib import Path
import csv, math, sys
from collections import Counter

OUT = Path("structures/phase1A/chromophores")
OUT.mkdir(parents=True, exist_ok=True)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception as e:
    print("ERROR: RDKit is required for this script.")
    print("Install with: conda install -c conda-forge rdkit")
    raise SystemExit(1)

# Pyrene, C16H10
SMILES = "c1cc2ccc3cccc4ccc(c1)c2c34"

mol = Chem.MolFromSmiles(SMILES)
if mol is None:
    raise RuntimeError("Could not parse pyrene SMILES.")

mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=12345)
AllChem.UFFOptimizeMolecule(mol, maxIters=10000)

conf = mol.GetConformer()

atoms = []
for i, atom in enumerate(mol.GetAtoms()):
    p = conf.GetAtomPosition(i)
    atoms.append({
        "idx": i + 1,
        "elem": atom.GetSymbol(),
        "x": p.x,
        "y": p.y,
        "z": p.z
    })

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

# Bonds from RDKit topology
bonds = []
for b in mol.GetBonds():
    i = b.GetBeginAtomIdx()
    j = b.GetEndAtomIdx()
    d = dist(atoms[i], atoms[j])
    bonds.append((i+1, j+1, atoms[i]["elem"], atoms[j]["elem"], d, str(b.GetBondType())))

# Nonbonded overlaps
bond_pairs = {tuple(sorted((i,j))) for i,j,_,_,_,_ in bonds}
overlaps = []
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        if tuple(sorted((i+1,j+1))) in bond_pairs:
            continue
        d = dist(atoms[i], atoms[j])
        if d < 1.00:
            overlaps.append((i+1,j+1,atoms[i]["elem"],atoms[j]["elem"],d))

# Planarity via RMS distance from best-fit plane
# simple PCA using numpy
try:
    import numpy as np
except Exception:
    print("ERROR: numpy required.")
    raise SystemExit(1)

coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms if a["elem"] == "C"])
center = coords.mean(axis=0)
X = coords - center
_, _, vh = np.linalg.svd(X, full_matrices=False)
normal = vh[-1]
dist_plane = X @ normal
plane_rms = float(np.sqrt(np.mean(dist_plane**2)))
plane_max = float(np.max(np.abs(dist_plane)))

counts = Counter(a["elem"] for a in atoms)

# XYZ
with (OUT/"pyrene.xyz").open("w") as f:
    f.write(f"{len(atoms)}\n")
    f.write("pyrene C16H10; RDKit UFF optimized; Phase 1A.2\n")
    for a in atoms:
        f.write(f'{a["elem"]:2s} {a["x"]: .8f} {a["y"]: .8f} {a["z"]: .8f}\n')

# PDB
with (OUT/"pyrene.pdb").open("w") as f:
    f.write("REMARK pyrene C16H10; RDKit UFF optimized; Phase 1A.2\n")
    for a in atoms:
        f.write(
            f"HETATM{a['idx']:5d} {a['elem']:>2s}   PYR A   1    "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            f"  1.00  0.00          {a['elem']:>2s}\n"
        )
    f.write("END\n")

with (OUT/"audit_pyrene_bonds.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["atom_i","atom_j","elem_i","elem_j","distance_A","bond_type"])
    for row in bonds:
        w.writerow(row)

with (OUT/"audit_pyrene_overlaps.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["atom_i","atom_j","elem_i","elem_j","distance_A"])
    for row in overlaps:
        w.writerow(row)

with (OUT/"audit_pyrene_summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["formula_expected","C16H10"])
    w.writerow(["C_atoms",counts["C"]])
    w.writerow(["H_atoms",counts["H"]])
    w.writerow(["total_atoms",len(atoms)])
    w.writerow(["bond_count",len(bonds)])
    w.writerow(["nonbonded_overlap_count",len(overlaps)])
    w.writerow(["carbon_plane_rms_A",f"{plane_rms:.6f}"])
    w.writerow(["carbon_plane_max_abs_A",f"{plane_max:.6f}"])

print("Wrote pyrene files to", OUT)
print("C atoms:", counts["C"])
print("H atoms:", counts["H"])
print("total atoms:", len(atoms))
print("bonds:", len(bonds))
print("nonbonded overlaps:", len(overlaps))
print("carbon plane RMS A:", f"{plane_rms:.6f}")
print("carbon plane max abs A:", f"{plane_max:.6f}")
