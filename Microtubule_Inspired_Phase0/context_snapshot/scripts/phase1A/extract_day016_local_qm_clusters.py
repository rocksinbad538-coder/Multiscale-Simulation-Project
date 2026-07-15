from pathlib import Path
import numpy as np
import pandas as pd

FRAMEDIR = Path("runs/phase1A/day016_md_bath_extraction/frames_0to100ps_every5ps")
OUTDIR = Path("runs/phase1A/day016_md_bath_extraction/local_qm_clusters")
OUTDIR.mkdir(parents=True, exist_ok=True)

MAPPING = [
    ("PYR2", 1681, 1706),
    ("PYR3", 1707, 1732),
    ("PYR4", 1733, 1758),
    ("PYR5", 1759, 1784),
]

PAIRS = [
    ("PYR2_PYR3", ["PYR2", "PYR3"]),
    ("PYR3_PYR4", ["PYR3", "PYR4"]),
    ("PYR4_PYR5", ["PYR4", "PYR5"]),
]

CUTOFF_NM = 0.50  # 5 Å first hydration shell proxy

def parse_gro(path):
    lines = path.read_text(errors="ignore").splitlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    box = lines[-1]
    atoms = []
    for idx, line in enumerate(lines[2:-1], start=1):
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20])
        xyz = np.array([
            float(line[20:28]),
            float(line[28:36]),
            float(line[36:44]),
        ])
        atoms.append({
            "index": idx,
            "resid": resid,
            "resname": resname,
            "atomname": atomname,
            "atomnr": atomnr,
            "xyz": xyz,
            "line": line,
        })
    return title, natoms, atoms, box

def write_xyz(path, selected_atoms):
    elem_map = {
        "C": "C", "H": "H",
        "OW": "O", "HW1": "H", "HW2": "H", "MW": "X",
    }
    with path.open("w") as f:
        f.write(f"{len(selected_atoms)}\n")
        f.write(f"{path.stem}\n")
        for a in selected_atoms:
            name = a["atomname"]
            elem = elem_map.get(name, name[0])
            x, y, z = a["xyz"] * 10.0  # nm -> Å
            f.write(f"{elem:2s} {x:14.6f} {y:14.6f} {z:14.6f}\n")

def min_dist_to_ref(sol_atoms, ref_xyz):
    xyz = np.array([a["xyz"] for a in sol_atoms])
    d = np.linalg.norm(xyz[:, None, :] - ref_xyz[None, :, :], axis=2)
    return float(d.min())

summary = []

for frame in sorted(FRAMEDIR.glob("frame*.gro"), key=lambda p: int(p.stem.replace("frame",""))):
    frame_id = int(frame.stem.replace("frame",""))
    _, _, atoms, _ = parse_gro(frame)

    atom_by_idx = {a["index"]: a for a in atoms}

    pyr_atoms = {}
    for label, start, end in MAPPING:
        block = [atom_by_idx[i] for i in range(start, end + 1)]
        pyr_atoms[label] = block

    waters = {}
    for a in atoms:
        if a["resname"] == "SOL":
            waters.setdefault(a["resid"], []).append(a)

    # Monomer clusters
    for label, _, _ in MAPPING:
        ref = np.array([a["xyz"] for a in pyr_atoms[label]])
        selected_waters = []
        for resid, wat_atoms in waters.items():
            if min_dist_to_ref(wat_atoms, ref) <= CUTOFF_NM:
                selected_waters.extend(wat_atoms)

        selected = pyr_atoms[label] + selected_waters
        out = OUTDIR / f"frame{frame_id:03d}_{label}_water5A.xyz"
        write_xyz(out, selected)

        summary.append({
            "frame": frame_id,
            "cluster": label,
            "type": "monomer",
            "n_pyr_atoms": len(pyr_atoms[label]),
            "n_water_atoms": len(selected_waters),
            "n_water_molecules": len(selected_waters) // 4,
            "n_total_atoms": len(selected),
            "xyz_file": str(out),
        })

    # Pair clusters
    for pair_label, labels in PAIRS:
        ref_atoms = []
        for lab in labels:
            ref_atoms.extend(pyr_atoms[lab])
        ref = np.array([a["xyz"] for a in ref_atoms])

        selected_waters = []
        for resid, wat_atoms in waters.items():
            if min_dist_to_ref(wat_atoms, ref) <= CUTOFF_NM:
                selected_waters.extend(wat_atoms)

        selected = ref_atoms + selected_waters
        out = OUTDIR / f"frame{frame_id:03d}_{pair_label}_water5A.xyz"
        write_xyz(out, selected)

        summary.append({
            "frame": frame_id,
            "cluster": pair_label,
            "type": "pair",
            "n_pyr_atoms": len(ref_atoms),
            "n_water_atoms": len(selected_waters),
            "n_water_molecules": len(selected_waters) // 4,
            "n_total_atoms": len(selected),
            "xyz_file": str(out),
        })

df = pd.DataFrame(summary)
df.to_csv(OUTDIR / "local_qm_cluster_summary.csv", index=False)

print(df.groupby(["type", "cluster"])[["n_water_molecules", "n_total_atoms"]].describe().to_string())
print("Wrote:", OUTDIR / "local_qm_cluster_summary.csv")
