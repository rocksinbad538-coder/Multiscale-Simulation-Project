#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

def participation_ratio(v):
    w = np.abs(v) ** 2
    return 1.0 / np.sum(w ** 2)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--e2", type=float, required=True)
    p.add_argument("--e3", type=float, required=True)
    p.add_argument("--e4", type=float, required=True)
    p.add_argument("--e5", type=float, required=True)
    p.add_argument("--j23", type=float, default=0.0, help="meV")
    p.add_argument("--j24", type=float, default=0.0, help="meV")
    p.add_argument("--j25", type=float, default=0.0, help="meV")
    p.add_argument("--j34", type=float, required=True, help="meV")
    p.add_argument("--j35", type=float, default=0.0, help="meV")
    p.add_argument("--j45", type=float, required=True, help="meV")
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sites = ["PYR2", "PYR3", "PYR4", "PYR5"]
    e = np.array([args.e2, args.e3, args.e4, args.e5], dtype=float)

    H = np.diag(e)
    couplings = {
        (0, 1): args.j23,
        (0, 2): args.j24,
        (0, 3): args.j25,
        (1, 2): args.j34,
        (1, 3): args.j35,
        (2, 3): args.j45,
    }

    for (i, j), mev in couplings.items():
        H[i, j] = H[j, i] = mev / 1000.0

    vals, vecs = np.linalg.eigh(H)

    h_df = pd.DataFrame(H, index=sites, columns=sites)
    h_df.to_csv(outdir / f"{args.label}_hamiltonian_eV.csv")

    rows = []
    for k, val in enumerate(vals):
        v = vecs[:, k]
        w = np.abs(v) ** 2
        dom = sites[int(np.argmax(w))]
        rows.append({
            "case": args.label,
            "exciton_state": f"X{k+1}",
            "energy_eV": val,
            "shift_meV_from_mean_site": (val - np.mean(e)) * 1000.0,
            "participation_ratio": participation_ratio(v),
            "dominant_site": dom,
            **{f"weight_{s}": w[i] for i, s in enumerate(sites)}
        })

    pd.DataFrame(rows).to_csv(outdir / f"{args.label}_exciton_eigenstates.csv", index=False)

    ev_df = pd.DataFrame(vecs, index=sites, columns=[f"X{i+1}" for i in range(4)])
    ev_df.to_csv(outdir / f"{args.label}_exciton_eigenvectors.csv")

    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\nWrote: {outdir}")

if __name__ == "__main__":
    main()
