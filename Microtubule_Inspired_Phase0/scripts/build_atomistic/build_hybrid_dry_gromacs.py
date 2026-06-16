#!/usr/bin/env python3
from pathlib import Path
import csv
import math

XYZ = Path("structures/phase1A/accepted/hybrid_dry/hbn_pyrene_4_dry.xyz")
PYR_ITP_SRC = Path("parameters/phase1A/accepted/pyrene_gromacs/pyrene.itp")
OUT = Path("parameters/phase1A/hybrid_dry_gromacs")
OUT.mkdir(parents=True, exist_ok=True)

gro = OUT / "hbn_pyrene_4_dry.gro"
top = OUT / "hbn_pyrene_4_dry.top"
hbn_itp = OUT / "hbn_fixed_dummy.itp"
pyrene_itp = OUT / "pyrene.itp"
audit = OUT / "audit_hybrid_dry_gromacs.csv"

lines = XYZ.read_text().splitlines()
n = int(lines[0])
atoms = []
for i, line in enumerate(lines[2:2+n], start=1):
    p = line.split()
    atoms.append({
        "idx": i,
        "elem": p[0],
        "x": float(p[1]) / 10.0,  # A -> nm
        "y": float(p[2]) / 10.0,
        "z": float(p[3]) / 10.0,
    })

hbn = atoms[:1680]
pyrs = atoms[1680:]

if len(hbn) != 1680:
    raise SystemExit("Unexpected hBN atom count.")
if len(pyrs) != 104:
    raise SystemExit("Unexpected pyrene atom count.")

# Copy pyrene itp locally
pyrene_itp.write_text(PYR_ITP_SRC.read_text())

# Build hBN dummy/fixed ITP
# Purpose: topology-compatible, non-interacting scaffold placeholder.
# Charges and LJ epsilons are zero; this is NOT physical h-BN parameterization.
with hbn_itp.open("w") as f:
    f.write("; h-BN fixed dummy scaffold for topology assembly only\n")
    f.write("; NOT a physical h-BN force field\n\n")
    f.write("[ moleculetype ]\n")
    f.write("; name  nrexcl\n")
    f.write("HBN     1\n\n")
    f.write("[ atoms ]\n")
    f.write("; nr type resnr residue atom cgnr charge mass\n")
    for a in hbn:
        typ = "B0" if a["elem"] == "B" else "N0"
        mass = 10.81 if a["elem"] == "B" else 14.007
        f.write(f"{a['idx']:6d} {typ:4s} 1 HBN {a['elem']}{a['idx']:04d} {a['idx']:6d} 0.000000 {mass:.6f}\n")

# Determine box
xmin, xmax = min(a["x"] for a in atoms), max(a["x"] for a in atoms)
ymin, ymax = min(a["y"] for a in atoms), max(a["y"] for a in atoms)
zmin, zmax = min(a["z"] for a in atoms), max(a["z"] for a in atoms)
pad = 2.0
lx = (xmax - xmin) + 2*pad
ly = (ymax - ymin) + 2*pad
lz = (zmax - zmin) + 2*pad

# shift to positive box
for a in atoms:
    a["x"] = a["x"] - xmin + pad
    a["y"] = a["y"] - ymin + pad
    a["z"] = a["z"] - zmin + pad

# GRO
with gro.open("w") as f:
    f.write("h-BN fixed dummy scaffold + 4 pyrenes dry Phase 1A.7\n")
    f.write(f"{len(atoms):5d}\n")

    # HBN residue
    for local_i, a in enumerate(hbn, start=1):
        f.write(f"{1%100000:5d}{'HBN':<5s}{a['elem'][:5]:>5s}{a['idx']%100000:5d}"
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}\n")

    # Pyrenes as four residues, atom naming must be compatible enough with topology order
    for k, a in enumerate(pyrs, start=1):
        pyr_id = (k-1)//26 + 1
        local = (k-1)%26 + 1
        resnr = 1 + pyr_id
        atom_name = "C" if local == 1 else (f"C{local-1}" if local <= 16 else ("H" if local == 17 else f"H{local-17}"))
        f.write(f"{resnr%100000:5d}{'PYR':<5s}{atom_name[:5]:>5s}{(1680+k)%100000:5d}"
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}\n")

    f.write(f"{lx:10.5f}{ly:10.5f}{lz:10.5f}\n")

# TOP
with top.open("w") as f:
    f.write("; Dry h-BN + 4 pyrene GROMACS topology, Phase 1A.7\n")
    f.write("; h-BN currently dummy/fixed/non-interacting: not physical h-BN FF\n\n")
    f.write("[ defaults ]\n")
    f.write("; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n")
    f.write("1 2 yes 0.5 0.83333333\n\n")

    f.write("[ atomtypes ]\n")
    f.write("; name at.num mass charge ptype sigma epsilon\n")
    f.write("B0 5 10.810000 0.000000 A 0.000000 0.000000\n")
    f.write("N0 7 14.007000 0.000000 A 0.000000 0.000000\n")
    f.write("; Pyrene atomtypes are defined in pyrene.itp/topology conversion if needed\n\n")

    # pyrene.itp already contains [ moleculetype ], but not [ atomtypes ] if extracted only from moleculetype.
    # The standalone pyrene.top contains atomtypes; include them manually for ca/ha.
    f.write("; GAFF2 atomtypes used by pyrene\n")
    f.write("ca 6 12.010000 0.000000 A 0.33152123 0.4133792\n")
    f.write("ha 1 1.008000 0.000000 A 0.26254785 0.0673624\n\n")

    f.write('#include "hbn_fixed_dummy.itp"\n')
    f.write('#include "pyrene.itp"\n\n')
    f.write("[ system ]\n")
    f.write("Dry h-BN fixed dummy scaffold + 4 pyrenes\n\n")
    f.write("[ molecules ]\n")
    f.write("; molname count\n")
    f.write("HBN 1\n")
    f.write("PYR 4\n")

# Audit
with audit.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["gro_exists", gro.exists()])
    w.writerow(["top_exists", top.exists()])
    w.writerow(["hbn_itp_exists", hbn_itp.exists()])
    w.writerow(["pyrene_itp_exists", pyrene_itp.exists()])
    w.writerow(["hbn_atoms", len(hbn)])
    w.writerow(["pyrene_atoms_total", len(pyrs)])
    w.writerow(["pyrene_count", 4])
    w.writerow(["total_atoms", len(atoms)])
    w.writerow(["box_x_nm", f"{lx:.5f}"])
    w.writerow(["box_y_nm", f"{ly:.5f}"])
    w.writerow(["box_z_nm", f"{lz:.5f}"])

print(audit.read_text())
