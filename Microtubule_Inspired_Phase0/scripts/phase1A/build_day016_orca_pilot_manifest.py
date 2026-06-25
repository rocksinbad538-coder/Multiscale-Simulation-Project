from pathlib import Path
import pandas as pd

BASE = Path("runs/phase1A/day016_md_bath_extraction")
CLUSTERS = BASE / "local_qm_clusters"
OUT = BASE / "orca_pilot_manifest"
OUT.mkdir(parents=True, exist_ok=True)

summary = pd.read_csv(CLUSTERS / "local_qm_cluster_summary.csv")

pilot_frames = [0, 10, 20]
pilot_clusters = ["PYR2", "PYR3", "PYR4", "PYR5", "PYR2_PYR3", "PYR3_PYR4", "PYR4_PYR5"]

pilot = summary[
    summary["frame"].isin(pilot_frames)
    & summary["cluster"].isin(pilot_clusters)
].copy()

pilot["recommended_role"] = pilot["type"].map({
    "monomer": "site_energy",
    "pair": "coupling"
})

pilot["calculation_level"] = "pilot_screening"
pilot["note"] = (
    "Geometry extracted from 100 ps accepted MD trajectory; "
    "water shell is explicit in XYZ but should be reduced or converted to embedding before systematic TDDFT."
)

pilot = pilot[
    [
        "frame",
        "cluster",
        "type",
        "recommended_role",
        "n_pyr_atoms",
        "n_water_molecules",
        "n_total_atoms",
        "calculation_level",
        "xyz_file",
        "note",
    ]
]

pilot.to_csv(OUT / "day016_orca_pilot_manifest.csv", index=False)

md = """# Day016 ORCA Pilot Manifest

## Purpose

This manifest defines a small pilot subset for testing the MD-to-electronic-extraction workflow before launching systematic calculations.

## Pilot frames

- frame000: beginning of accepted 100 ps trajectory
- frame010: middle of accepted 100 ps trajectory
- frame020: end of accepted 100 ps trajectory

## Pilot clusters

Monomers:

- PYR2
- PYR3
- PYR4
- PYR5

Pairs:

- PYR2_PYR3
- PYR3_PYR4
- PYR4_PYR5

## Important caveat

The current extracted clusters contain explicit nearby water molecules. These are useful for geometry/electrostatic screening but are too large for routine TDDFT across all frames and all clusters.

The production route should either:

1. Convert nearby waters into point-charge embedding, or
2. Use a smaller explicit water shell, or
3. Run only a limited number of explicit-water TDDFT pilots to calibrate a cheaper embedding approximation.

## Current scientific role

This pilot does not yet estimate sub-100 fs bath correlation times because the accepted trajectory is sampled every 0.5 ps. It is used to validate the MD-frame-to-electronic-input pipeline.
"""

(OUT / "DAY016_ORCA_PILOT_MANIFEST.md").write_text(md)

print(pilot.to_string(index=False))
print("Wrote:", OUT)
