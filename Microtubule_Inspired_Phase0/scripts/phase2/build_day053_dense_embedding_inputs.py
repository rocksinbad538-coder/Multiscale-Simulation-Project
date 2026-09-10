from pathlib import Path
import numpy as np
import pandas as pd

FRAMEDIR = Path(
    "runs/phase2/day053_dense_electronic_pilot/"
    "frames_100to100p5ps_every50fs"
)

BASE = Path("runs/phase2/day053_dense_electronic_pilot")
CLUSTER_DIR = BASE / "local_qm_clusters"
ORCA_DIR = BASE / "orca_inputs"

CLUSTER_DIR.mkdir(parents=True, exist_ok=True)
ORCA_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF_NM = 0.50

MAPPING = {
    "PYR3": (1707, 1732),
    "PYR4": (1733, 1758),
}

CHARGES = {
    "OW": 0.0,
    "HW1": 0.5564,
    "HW2": 0.5564,
    "MW": -1.1128,
}


def parse_gro(path):
    lines = path.read_text(errors="ignore").splitlines()
    natoms = int(lines[1].strip())

    atoms = []
    for idx, line in enumerate(lines[2:-1], start=1):
        atoms.append({
            "index": idx,
            "resid": int(line[0:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "xyz": np.array([
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            ]),
        })

    box = np.array([float(x) for x in lines[-1].split()[:3]])
    return natoms, atoms, box


def min_dist(sol_atoms, ref_xyz):
    xyz = np.array([a["xyz"] for a in sol_atoms])
    d = np.linalg.norm(
        xyz[:, None, :] - ref_xyz[None, :, :],
        axis=2,
    )
    return float(d.min())


def write_qm_xyz(path, atoms, stem):
    with path.open("w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{stem} QM pyrene only\n")
        for a in atoms:
            x, y, z = a["xyz"] * 10.0
            element = "C" if a["atomname"].startswith("C") else "H"
            f.write(
                f"{element:2s} "
                f"{x:14.6f} {y:14.6f} {z:14.6f}\n"
            )


def write_pc(path, waters):
    ordered = []
    for resid in sorted(waters):
        ordered.extend(waters[resid])

    with path.open("w") as f:
        f.write(f"{len(ordered)}\n")

        for i, a in enumerate(ordered):
            mod = i % 4

            if mod == 0:
                q = CHARGES["OW"]
            elif mod == 1:
                q = CHARGES["HW1"]
            elif mod == 2:
                q = CHARGES["HW2"]
            else:
                q = CHARGES["MW"]

            x, y, z = a["xyz"] * 10.0
            f.write(
                f"{q: .8f} "
                f"{x:14.6f} {y:14.6f} {z:14.6f}\n"
            )

    return len(ordered)


def write_inp(path, qm_atoms, pc_name):
    with path.open("w") as f:
        f.write("! wB97X-D3 def2-SVP def2/J RIJCOSX TightSCF\n")
        f.write("%maxcore 4096\n")
        f.write("%tddft\n")
        f.write("  nroots 10\n")
        f.write("  tda true\n")
        f.write("end\n")
        f.write(f'%pointcharges "{pc_name}"\n')
        f.write("* xyz 0 1\n")

        for a in qm_atoms:
            x, y, z = a["xyz"] * 10.0
            element = "C" if a["atomname"].startswith("C") else "H"
            f.write(
                f"{element:2s} "
                f"{x:14.6f} {y:14.6f} {z:14.6f}\n"
            )

        f.write("*\n")


rows = []

frames = sorted(
    FRAMEDIR.glob("frame*.gro"),
    key=lambda p: int(p.stem.replace("frame", "")),
)

for gro in frames:
    frame_id = int(gro.stem.replace("frame", ""))
    time_ps = 100.0 + 0.05 * frame_id

    natoms, atoms, box = parse_gro(gro)

    atom_by_idx = {a["index"]: a for a in atoms}

    waters = {}
    for a in atoms:
        if a["resname"] == "SOL":
            waters.setdefault(a["resid"], []).append(a)

    for site, (start, end) in MAPPING.items():
        qm_atoms = [
            atom_by_idx[i]
            for i in range(start, end + 1)
        ]

        ref_xyz = np.array([a["xyz"] for a in qm_atoms])

        selected_waters = {}

        for resid, wat_atoms in waters.items():
            if min_dist(wat_atoms, ref_xyz) <= CUTOFF_NM:
                selected_waters[resid] = wat_atoms

        stem = f"dense{frame_id:03d}_{site}"

        xyz_path = CLUSTER_DIR / f"{stem}_water5A.xyz"

        cluster_atoms = list(qm_atoms)
        for resid in sorted(selected_waters):
            cluster_atoms.extend(selected_waters[resid])

        with xyz_path.open("w") as f:
            f.write(f"{len(cluster_atoms)}\n")
            f.write(f"{stem} 5A hydration-shell cluster\n")

            for a in cluster_atoms:
                name = a["atomname"]
                elem = {
                    "C": "C",
                    "H": "H",
                    "OW": "O",
                    "HW1": "H",
                    "HW2": "H",
                    "MW": "X",
                }.get(name, name[0])

                x, y, z = a["xyz"] * 10.0
                f.write(
                    f"{elem:2s} "
                    f"{x:14.6f} {y:14.6f} {z:14.6f}\n"
                )

        outdir = ORCA_DIR / stem
        outdir.mkdir(parents=True, exist_ok=True)

        qm_path = outdir / f"{stem}_qm.xyz"
        pc_path = outdir / f"{stem}.pc"
        inp_path = outdir / f"{stem}.inp"

        write_qm_xyz(qm_path, qm_atoms, stem)
        n_pc = write_pc(pc_path, selected_waters)
        write_inp(inp_path, qm_atoms, pc_path.name)

        rows.append({
            "frame_id": frame_id,
            "time_ps": time_ps,
            "site": site,
            "source_gro": str(gro),
            "cutoff_nm": CUTOFF_NM,
            "n_qm_atoms": len(qm_atoms),
            "n_water_molecules": len(selected_waters),
            "n_point_charges": n_pc,
            "point_charge_total": 0.0,
            "orca_input": str(inp_path),
            "qm_xyz": str(qm_path),
            "point_charge_file": str(pc_path),
        })


df = pd.DataFrame(rows)
manifest = BASE / "DAY053_DENSE_ELECTRONIC_PILOT_MANIFEST.csv"
df.to_csv(manifest, index=False)

print("=" * 78)
print("DAY053 DENSE ELECTRONIC PILOT INPUT BUILD")
print("=" * 78)
print("N_FRAMES =", df["frame_id"].nunique())
print("N_SITES =", df["site"].nunique())
print("N_JOBS =", len(df))
print("TIME_MIN_PS =", df["time_ps"].min())
print("TIME_MAX_PS =", df["time_ps"].max())
print("UNIQUE_DT_PS =", sorted(df["time_ps"].drop_duplicates().diff().dropna().unique()))
print("POINT_CHARGE_MOD4_OK =", bool((df["n_point_charges"] % 4 == 0).all()))
print("MIN_WATERS =", int(df["n_water_molecules"].min()))
print("MAX_WATERS =", int(df["n_water_molecules"].max()))
print("MANIFEST =", manifest)
