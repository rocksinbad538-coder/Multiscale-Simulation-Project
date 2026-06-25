from pathlib import Path
import json
import numpy as np
import pandas as pd

IN = Path("runs/phase1A/day016_md_bath_extraction/orca_embedding_analysis/embedding_pilot_summary.csv")
OUT = Path("runs/phase1A/day016_md_bath_extraction/site_energy_trajectory")
OUT.mkdir(parents=True, exist_ok=True)

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]

df = pd.read_csv(IN)

required = ["frame", "cluster", "S1_eV", "terminated_normally", "tddft_finished"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing required columns: {missing}")

bad = df[~(df["terminated_normally"] & df["tddft_finished"])]
if len(bad):
    raise SystemExit("Some ORCA jobs did not finish correctly:\n" + bad.to_string(index=False))

site = df.pivot(index="frame", columns="cluster", values="S1_eV").reset_index()
missing_chrom = [c for c in CHROMOPHORES if c not in site.columns]
if missing_chrom:
    raise SystemExit(f"Missing chromophores in site-energy table: {missing_chrom}")

site = site[["frame"] + CHROMOPHORES].sort_values("frame")
site.to_csv(OUT / "site_energy_trajectory.csv", index=False)

values = site[CHROMOPHORES]
centered = site.copy()
centered[CHROMOPHORES] = values.sub(values.mean(axis=1), axis=0)
centered.to_csv(OUT / "site_energy_trajectory_centered.csv", index=False)

global_centered = site.copy()
global_centered[CHROMOPHORES] = values - values.values.mean()
global_centered.to_csv(OUT / "site_energy_trajectory_global_centered.csv", index=False)

np.save(OUT / "site_energy_trajectory_eV.npy", site[CHROMOPHORES].to_numpy())
np.save(OUT / "site_energy_trajectory_centered_eV.npy", centered[CHROMOPHORES].to_numpy())

hamdir = OUT / "hamiltonian_diagonal_frames"
hamdir.mkdir(exist_ok=True)

for _, row in centered.iterrows():
    frame = int(row["frame"])
    H = np.diag([row[c] for c in CHROMOPHORES])
    pd.DataFrame(H, index=CHROMOPHORES, columns=CHROMOPHORES).to_csv(
        hamdir / f"H_diagonal_centered_frame{frame:03d}.csv"
    )

metadata = {
    "day": "016",
    "source": str(IN),
    "quantity": "TDDFT/TDA S1 site energies",
    "units": "eV",
    "centering": {
        "site_energy_trajectory_centered.csv": "per-frame mean removed",
        "site_energy_trajectory_global_centered.csv": "global mean over all frames and chromophores removed"
    },
    "frames": [int(x) for x in site["frame"].tolist()],
    "chromophores": CHROMOPHORES,
    "embedding": "electrostatic point-charge embedding from MD water shell",
    "software": "ORCA",
    "hamiltonian_note": "Only diagonal site-energy terms are included here. Off-diagonal excitonic couplings J_ij(t) are not included yet.",
}
(OUT / "metadata.json").write_text(json.dumps(metadata, indent=2))

summary = OUT / "SITE_ENERGY_TRAJECTORY_DAY016.md"
with summary.open("w") as f:
    f.write("# Day016 site-energy trajectory\n\n")
    f.write("## Source\n\n")
    f.write(f"- Input: `{IN}`\n")
    f.write("- Quantity: embedded TDDFT/TDA S1 excitation energy\n")
    f.write("- Units: eV\n")
    f.write("- Chromophore order: PYR2, PYR3, PYR4, PYR5\n\n")
    f.write("## Site energies\n\n")
    f.write(site.to_string(index=False))
    f.write("\n\n## Per-frame centered diagonal terms\n\n")
    f.write(centered.to_string(index=False))
    f.write("\n\n## Statistics\n\n")
    stats = site[CHROMOPHORES].describe().T
    f.write(stats.to_string())
    f.write("\n\n## Interpretation\n\n")
    f.write(
        "This table is the first MD-derived embedded-TDDFT site-energy trajectory. "
        "It provides the diagonal part of the time-dependent excitonic Hamiltonian. "
        "At this stage, off-diagonal couplings are not included.\n"
    )

print("\nSite-energy trajectory:")
print(site.to_string(index=False))
print("\nPer-frame centered:")
print(centered.to_string(index=False))
print("\nGlobal-centered:")
print(global_centered.to_string(index=False))
print("\nWrote:", OUT)
