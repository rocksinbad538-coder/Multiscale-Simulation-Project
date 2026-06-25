from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path("runs/phase1A/day016_md_bath_extraction")

traj = pd.read_csv(
    ROOT /
    "site_energy_trajectory" /
    "site_energy_trajectory.csv"
)

outdir = ROOT / "hamiltonian_diagonals"

outdir.mkdir(parents=True, exist_ok=True)

chromophores = ["PYR2","PYR3","PYR4","PYR5"]

metadata = {
    "description":
        "Diagonal excitonic Hamiltonians extracted from embedded TDDFT calculations.",
    "units":"eV",
    "method":"TDA-TDDFT electrostatic embedding",
    "hamiltonian_type":"diagonal",
    "chromophores":chromophores,
    "source":"Day016 production pilot"
}

with open(outdir/"metadata.json","w") as f:
    json.dump(metadata,f,indent=2)

summary=[]

for _,row in traj.iterrows():

    frame=int(row.frame)

    H=np.diag(
        row[chromophores].to_numpy(dtype=float)
    )

    np.save(
        outdir/f"Hdiag_frame{frame:03d}.npy",
        H
    )

    pd.DataFrame(
        H,
        index=chromophores,
        columns=chromophores
    ).to_csv(
        outdir/f"Hdiag_frame{frame:03d}.csv"
    )

    summary.append({
        "frame":frame,
        "trace_eV":float(np.trace(H)),
        "mean_diag_eV":float(np.mean(np.diag(H)))
    })

pd.DataFrame(summary).to_csv(
    outdir/"hamiltonian_summary.csv",
    index=False
)

print("Hamiltonians exported to")
print(outdir)
