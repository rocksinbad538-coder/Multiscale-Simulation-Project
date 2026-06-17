#!/usr/bin/env python3
from pathlib import Path
import math, csv
import numpy as np

HBN_XYZ = Path("structures/phase1A/hbn_nt_initial.xyz")
PYR_XYZ = Path("structures/phase1A/chromophores/pyrene.xyz")
OUT = Path("structures/phase1A/hybrid_dry")
OUT.mkdir(parents=True, exist_ok=True)

# Placement parameters
radial_gap = 4.50       # approximate pi-surface to BN-surface separation, Angstrom
axial_margin = 10.0     # avoid open tube edges
n_pyrene = 4

def read_xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    atoms = []
    for k, line in enumerate(lines[2:2+n], start=1):
        parts = line.split()
        atoms.append({
            "elem": parts[0],
            "x": float(parts[1]),
            "y": float(parts[2]),
            "z": float(parts[3]),
            "source_idx": k
        })
    return atoms

def write_xyz(path, atoms, comment):
    with path.open("w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for a in atoms:
            f.write(f'{a["elem"]:2s} {a["x"]: .8f} {a["y"]: .8f} {a["z"]: .8f}\n')

def write_pdb(path, atoms):
    with path.open("w") as f:
        f.write("REMARK h-BN nanotube + 4 pyrene dry hybrid, Phase 1A.3\n")
        for i,a in enumerate(atoms,1):
            res = a.get("resname","HBN")
            chain = a.get("chain","A")
            resid = a.get("resid",1)
            f.write(
                f"HETATM{i:5d} {a['elem']:>2s}   {res:>3s} {chain}{resid:4d}    "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                f"  1.00  0.00          {a['elem']:>2s}\n"
            )
        f.write("END\n")

def center_atoms(atoms):
    coords = np.array([[a["x"],a["y"],a["z"]] for a in atoms])
    c = coords.mean(axis=0)
    out = []
    for a in atoms:
        b = dict(a)
        b["x"] -= c[0]
        b["y"] -= c[1]
        b["z"] -= c[2]
        out.append(b)
    return out

def rotation_matrix_from_axes(e1, e2, e3):
    # columns are target basis vectors
    return np.column_stack([e1, e2, e3])

def transform_pyrene(pyrene_atoms, origin, radial, tangent, axial):
    # pyrene local x -> tangent, y -> axial, z -> radial
    Rmat = rotation_matrix_from_axes(tangent, axial, radial)
    transformed = []
    for a in pyrene_atoms:
        v = np.array([a["x"], a["y"], a["z"]])
        p = origin + Rmat @ v
        b = dict(a)
        b["x"], b["y"], b["z"] = float(p[0]), float(p[1]), float(p[2])
        transformed.append(b)
    return transformed

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

hbn = read_xyz(HBN_XYZ)
pyr = center_atoms(read_xyz(PYR_XYZ))

# annotate hBN
for a in hbn:
    a["resname"] = "HBN"
    a["chain"] = "A"
    a["resid"] = 1
    a["chrom_id"] = 0

coords_hbn = np.array([[a["x"],a["y"],a["z"]] for a in hbn])
radius = float(np.mean(np.sqrt(coords_hbn[:,0]**2 + coords_hbn[:,1]**2)))
zmin = float(coords_hbn[:,2].min())
zmax = float(coords_hbn[:,2].max())
z_positions = np.linspace(zmin + axial_margin, zmax - axial_margin, n_pyrene)

pyrenes = []
mapping = []

for cid in range(1, n_pyrene+1):
    theta = 2.0 * math.pi * (cid-1) / n_pyrene
    radial = np.array([math.cos(theta), math.sin(theta), 0.0])
    tangent = np.array([-math.sin(theta), math.cos(theta), 0.0])
    axial = np.array([0.0, 0.0, 1.0])

    origin = radial * (radius + radial_gap)
    origin[2] = z_positions[cid-1]

    placed = transform_pyrene(pyr, origin, radial, tangent, axial)

    for local_idx, a in enumerate(placed, start=1):
        a["resname"] = "PYR"
        a["chain"] = "P"
        a["resid"] = cid
        a["chrom_id"] = cid
        pyrenes.append(a)
        mapping.append({
            "chromophore_id": cid,
            "local_atom_index": local_idx,
            "element": a["elem"],
            "global_atom_index": len(hbn) + len(pyrenes)
        })

atoms = hbn + pyrenes

# Audits
overlaps = []
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        ai, aj = atoms[i], atoms[j]
        d = dist(ai, aj)

        # Only audit inter-component nonbonded overlaps here.
        # Intramolecular pyrene bonds and hBN bonds were audited separately.
        same_hbn = ai["chrom_id"] == 0 and aj["chrom_id"] == 0
        same_pyr = ai["chrom_id"] != 0 and ai["chrom_id"] == aj["chrom_id"]
        if same_hbn or same_pyr:
            continue

        if d < 1.60:
            overlaps.append((i+1,j+1,ai["elem"],aj["elem"],ai["chrom_id"],aj["chrom_id"],d))

# Distance summaries
rows = []

# pyrene center positions
centers = {}
for cid in range(1, n_pyrene+1):
    subset = [a for a in atoms if a["chrom_id"] == cid]
    c = np.array([[a["x"],a["y"],a["z"]] for a in subset]).mean(axis=0)
    centers[cid] = c

for i in range(1, n_pyrene+1):
    for j in range(i+1, n_pyrene+1):
        dcen = float(np.linalg.norm(centers[i]-centers[j]))
        min_atom = min(
            dist(a,b)
            for a in atoms if a["chrom_id"] == i
            for b in atoms if b["chrom_id"] == j
        )
        rows.append(["pyrene_pyrene",i,j,dcen,min_atom])

for cid in range(1, n_pyrene+1):
    pyr_atoms = [a for a in atoms if a["chrom_id"] == cid]
    hbn_atoms = [a for a in atoms if a["chrom_id"] == 0]
    min_bn = min(dist(a,b) for a in pyr_atoms for b in hbn_atoms)
    c = centers[cid]
    radial_c = math.sqrt(c[0]**2 + c[1]**2)
    rows.append(["pyrene_hBN",cid,0,radial_c-radius,min_bn])

# Write files
write_xyz(
    OUT/"hbn_pyrene_4_dry.xyz",
    atoms,
    "h-BN nanotube + 4 pyrenes dry hybrid; Phase 1A.3"
)
write_pdb(OUT/"hbn_pyrene_4_dry.pdb", atoms)

with (OUT/"chromophore_mapping.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["chromophore_id","local_atom_index","element","global_atom_index"])
    w.writeheader()
    for r in mapping:
        w.writerow(r)

with (OUT/"audit_chromophore_distances.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["type","id_i","id_j","center_or_radial_distance_A","min_atom_atom_distance_A"])
    for r in rows:
        w.writerow(r)

with (OUT/"audit_overlaps.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["atom_i","atom_j","elem_i","elem_j","chrom_i","chrom_j","distance_A"])
    for r in overlaps:
        w.writerow(r)

with (OUT/"audit_summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["hBN_atoms",len(hbn)])
    w.writerow(["pyrene_count",n_pyrene])
    w.writerow(["pyrene_atoms_total",len(pyrenes)])
    w.writerow(["total_atoms",len(atoms)])
    w.writerow(["hBN_radius_A",f"{radius:.6f}"])
    w.writerow(["radial_gap_requested_A",f"{radial_gap:.6f}"])
    w.writerow(["zmin_A",f"{zmin:.6f}"])
    w.writerow(["zmax_A",f"{zmax:.6f}"])
    w.writerow(["intercomponent_overlap_count",len(overlaps)])

print("Wrote dry hBN + 4 pyrene hybrid to", OUT)
print("hBN atoms:", len(hbn))
print("pyrene atoms:", len(pyrenes))
print("total atoms:", len(atoms))
print("intercomponent overlaps:", len(overlaps))
