from pathlib import Path
import numpy as np
import pandas as pd

gro = Path("runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.gro")
mapfile = Path("runs/phase1A/day016_md_bath_extraction/md_to_model_pyr_mapping.csv")
outdir = Path("runs/phase1A/day016_md_bath_extraction")

mapping = pd.read_csv(mapfile)

lines = gro.read_text(errors="ignore").splitlines()
atoms = lines[2:-1]

coords = {}
for line in atoms:
    idx = int(line[15:20])
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    coords[idx] = np.array([x, y, z], float)

centers = []
for _, row in mapping.iterrows():
    inds = list(range(int(row.atom_start), int(row.atom_end) + 1))
    xyz = np.array([coords[i] for i in inds])
    c = xyz.mean(axis=0)
    centers.append({
        "md_label": row.md_label,
        "model_label": row.model_label,
        "md_resid": row.md_resid,
        "center_x_nm": c[0],
        "center_y_nm": c[1],
        "center_z_nm": c[2],
    })

centers = pd.DataFrame(centers)
centers.to_csv(outdir / "pyr_model_centers_final_gro.csv", index=False)

pairs = []
for i in range(len(centers)):
    for j in range(i+1, len(centers)):
        ci = centers.loc[i, ["center_x_nm","center_y_nm","center_z_nm"]].to_numpy(float)
        cj = centers.loc[j, ["center_x_nm","center_y_nm","center_z_nm"]].to_numpy(float)
        d = np.linalg.norm(ci-cj)
        pairs.append({
            "pair": f"{centers.loc[i,'model_label']}-{centers.loc[j,'model_label']}",
            "md_pair": f"{centers.loc[i,'md_label']}-{centers.loc[j,'md_label']}",
            "distance_nm": d,
            "distance_A": 10*d,
        })

pairs = pd.DataFrame(pairs)
pairs.to_csv(outdir / "pyr_pair_distances_final_gro.csv", index=False)

print("Centers:")
print(centers.to_string(index=False))
print("\nPair distances:")
print(pairs.to_string(index=False))
