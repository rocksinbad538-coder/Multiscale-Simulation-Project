from pathlib import Path
from collections import defaultdict

gro = Path("runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.gro")
outdir = Path("runs/phase1A/day016_md_bath_extraction")
outdir.mkdir(parents=True, exist_ok=True)

lines = gro.read_text(errors="ignore").splitlines()
atoms = lines[2:-1]

records = []
for i, line in enumerate(atoms, start=1):  # GROMACS atom index, 1-based
    if len(line) < 44:
        continue
    resid = int(line[0:5])
    resname = line[5:10].strip()
    atomname = line[10:15].strip()
    atomnr = int(line[15:20])
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    records.append((i, resid, resname, atomname, atomnr, x, y, z))

pyr = [r for r in records if r[2] == "PYR"]

by_resid = defaultdict(list)
for r in pyr:
    by_resid[r[1]].append(r)

summary = outdir / "pyr_residue_blocks_summary.csv"
with summary.open("w") as f:
    f.write("pyr_label,resid,n_atoms,first_index,last_index,center_x_nm,center_y_nm,center_z_nm\n")
    for k, resid in enumerate(sorted(by_resid), start=1):
        block = by_resid[resid]
        xs = [r[5] for r in block]
        ys = [r[6] for r in block]
        zs = [r[7] for r in block]
        f.write(
            f"PYR{k},{resid},{len(block)},{block[0][0]},{block[-1][0]},"
            f"{sum(xs)/len(xs):.6f},{sum(ys)/len(ys):.6f},{sum(zs)/len(zs):.6f}\n"
        )

detail = outdir / "pyr_atom_index_detail.csv"
with detail.open("w") as f:
    f.write("pyr_label,resid,atom_index,atom_name,x_nm,y_nm,z_nm\n")
    for k, resid in enumerate(sorted(by_resid), start=1):
        for r in by_resid[resid]:
            f.write(f"PYR{k},{resid},{r[0]},{r[3]},{r[5]:.6f},{r[6]:.6f},{r[7]:.6f}\n")

print(summary)
print(summary.read_text())
print(detail)
