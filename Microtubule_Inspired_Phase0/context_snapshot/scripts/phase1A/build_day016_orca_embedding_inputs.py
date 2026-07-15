from pathlib import Path
import pandas as pd

BASE = Path("runs/phase1A/day016_md_bath_extraction")
MANIFEST = BASE / "orca_pilot_manifest/day016_orca_pilot_manifest.csv"
OUT = BASE / "orca_embedding_pilot_inputs"
OUT.mkdir(parents=True, exist_ok=True)

CHARGES = {
    "OW": 0.0,
    "HW1": 0.5564,
    "HW2": 0.5564,
    "MW": -1.1128,
}

def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    atoms = []
    for line in lines[2:2+n]:
        p = line.split()
        atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return atoms

def write_orca_input(stem, qm_atoms, pc_file, outdir):
    inp = outdir / f"{stem}.inp"
    xyz = outdir / f"{stem}_qm.xyz"

    with xyz.open("w") as f:
        f.write(f"{len(qm_atoms)}\n{stem} QM pyrene only\n")
        for e,x,y,z in qm_atoms:
            f.write(f"{e:2s} {x:14.6f} {y:14.6f} {z:14.6f}\n")

    with inp.open("w") as f:
        f.write("! wB97X-D3 def2-SVP def2/J RIJCOSX TightSCF\n")
        f.write("%maxcore 4096\n")
        f.write("%tddft\n")
        f.write("  nroots 10\n")
        f.write("  tda true\n")
        f.write("end\n")
        f.write(f"%pointcharges \"{pc_file.name}\"\n")
        f.write("* xyz 0 1\n")
        for e,x,y,z in qm_atoms:
            f.write(f"{e:2s} {x:14.6f} {y:14.6f} {z:14.6f}\n")
        f.write("*\n")

    return inp, xyz

manifest = pd.read_csv(MANIFEST)
pilot = manifest[manifest["type"] == "monomer"].copy()

rows = []
for _, row in pilot.iterrows():
    atoms = read_xyz(row["xyz_file"])

    # First 26 atoms are PYR by construction.
    qm_atoms = atoms[:26]
    env_atoms = atoms[26:]

    stem = f"frame{int(row.frame):03d}_{row.cluster}_embedding"
    outdir = OUT / stem
    outdir.mkdir(parents=True, exist_ok=True)

    pc_path = outdir / f"{stem}.pc"
    n_pc = 0
    q_sum = 0.0

    with pc_path.open("w") as f:
        f.write(f"{len(env_atoms)}\n")
        for e,x,y,z in env_atoms:
            # The XYZ element was converted to O/H/X; infer charges by water order.
            # local_qm_clusters writes waters as OW, HW1, HW2, MW converted to O,H,H,X.
            # Use 4-site repeating order after the 26 PYR atoms.
            mod = n_pc % 4
            if mod == 0:
                q = CHARGES["OW"]
            elif mod == 1:
                q = CHARGES["HW1"]
            elif mod == 2:
                q = CHARGES["HW2"]
            else:
                q = CHARGES["MW"]
            f.write(f"{q: .8f} {x:14.6f} {y:14.6f} {z:14.6f}\n")
            q_sum += q
            n_pc += 1

    inp, xyz = write_orca_input(stem, qm_atoms, pc_path, outdir)

    rows.append({
        "frame": row.frame,
        "cluster": row.cluster,
        "role": "site_energy_embedding_pilot",
        "n_qm_atoms": len(qm_atoms),
        "n_point_charges": n_pc,
        "point_charge_total": q_sum,
        "orca_input": str(inp),
        "qm_xyz": str(xyz),
        "point_charge_file": str(pc_path),
    })

df = pd.DataFrame(rows)
df.to_csv(OUT / "orca_embedding_pilot_manifest.csv", index=False)

print(df.to_string(index=False))
print("Wrote:", OUT)
