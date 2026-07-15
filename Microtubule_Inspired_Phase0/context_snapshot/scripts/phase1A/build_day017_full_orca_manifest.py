from pathlib import Path
import re
import pandas as pd

BASE = Path("runs/phase1A/day016_md_bath_extraction")
CLUSTER_DIR = BASE / "local_qm_clusters"
OUTDIR = BASE / "orca_pilot_manifest"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT = OUTDIR / "day016_orca_pilot_manifest.csv"
BACKUP = OUTDIR / "day016_orca_pilot_manifest_pilot_backup.csv"

chromophores = ["PYR2", "PYR3", "PYR4", "PYR5"]
pat = re.compile(r"frame(\d{3})_(PYR[2-5])_water5A\.xyz$")

rows = []

for xyz in sorted(CLUSTER_DIR.glob("frame*_PYR*_water5A.xyz")):
    m = pat.match(xyz.name)
    if not m:
        continue

    frame = int(m.group(1))
    cluster = m.group(2)

    lines = xyz.read_text(errors="ignore").splitlines()
    try:
        n_total_atoms = int(lines[0].strip())
    except Exception:
        raise SystemExit(f"Could not read atom count from {xyz}")

    # Pyrene monomer has 26 QM atoms; remaining atoms are water-shell atoms.
    n_pyr_atoms = 26
    n_water_atoms = n_total_atoms - n_pyr_atoms
    if n_water_atoms < 0 or n_water_atoms % 4 != 0:
        raise SystemExit(
            f"Unexpected atom count for {xyz}: total={n_total_atoms}, "
            f"water_atoms={n_water_atoms}"
        )

    rows.append({
        "frame": frame,
        "cluster": cluster,
        "type": "monomer",
        "recommended_role": "site_energy",
        "n_pyr_atoms": n_pyr_atoms,
        "n_water_molecules": n_water_atoms // 4,
        "n_total_atoms": n_total_atoms,
        "calculation_level": "production_screening",
        "xyz_file": str(xyz),
        "note": (
            "Full-production monomer site-energy input reconstructed from "
            "local QM cluster XYZ files; water shell will be converted to "
            "electrostatic point-charge embedding."
        ),
    })

df = pd.DataFrame(rows).sort_values(["frame", "cluster"]).reset_index(drop=True)

expected_frames = list(range(21))
found_frames = sorted(df["frame"].unique().tolist())
if found_frames != expected_frames:
    raise SystemExit(f"Unexpected frames. Found {found_frames}, expected {expected_frames}")

counts = df.groupby("frame")["cluster"].nunique()
bad = counts[counts != 4]
if len(bad):
    raise SystemExit("Incomplete chromophore coverage:\n" + bad.to_string())

if OUT.exists() and not BACKUP.exists():
    BACKUP.write_text(OUT.read_text())

df.to_csv(OUT, index=False)

print("Wrote:", OUT)
print("Backup:", BACKUP if BACKUP.exists() else "not created")
print("n_frames:", df["frame"].nunique())
print("n_monomer_jobs:", len(df))
print()
print(df.groupby("frame")["cluster"].apply(lambda x: ",".join(x)).to_string())
