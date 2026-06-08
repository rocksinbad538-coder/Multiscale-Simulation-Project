#!/usr/bin/env python3

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_lammps_dump(path):
    path = Path(path)
    with path.open("r", errors="ignore") as f:
        while True:
            line = f.readline()
            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError("Expected NUMBER OF ATOMS block")
            n_atoms = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("Expected BOX BOUNDS block")
            bounds = [f.readline().split() for _ in range(3)]

            line = f.readline()
            if not line.startswith("ITEM: ATOMS"):
                raise RuntimeError("Expected ATOMS block")

            cols = line.split()[2:]
            idx = {c: i for i, c in enumerate(cols)}

            required = ["id", "mol", "type", "x", "y", "z"]
            missing = [c for c in required if c not in idx]
            if missing:
                raise RuntimeError(f"Missing columns in dump: {missing}. Available: {cols}")

            data = np.empty((n_atoms, 6), dtype=float)

            for i in range(n_atoms):
                p = f.readline().split()
                data[i, 0] = int(p[idx["id"]])
                data[i, 1] = int(p[idx["mol"]])
                data[i, 2] = int(p[idx["type"]])
                data[i, 3] = float(p[idx["x"]])
                data[i, 4] = float(p[idx["y"]])
                data[i, 5] = float(p[idx["z"]])

            yield timestep, data, bounds


def compute_dipole_unit_vectors(frame, oxygen_type, hydrogen_type):
    """
    Reconstruct water dipole directions from O and two H atoms by molecule ID.

    Dipole direction used here:
        u = normalize( (r_H1 + r_H2)/2 - r_O )

    This gives the molecular orientation from oxygen toward the hydrogen midpoint.
    """
    oxy = frame[frame[:, 2] == oxygen_type]
    hyd = frame[frame[:, 2] == hydrogen_type]

    if len(oxy) == 0 or len(hyd) == 0:
        raise RuntimeError("No oxygen or hydrogen atoms found with the requested types.")

    oxy_by_mol = {}
    for row in oxy:
        mol = int(row[1])
        oxy_by_mol[mol] = row[3:6]

    h_by_mol = {}
    for row in hyd:
        mol = int(row[1])
        h_by_mol.setdefault(mol, []).append(row[3:6])

    mols = []
    uvecs = []

    for mol, ro in oxy_by_mol.items():
        hs = h_by_mol.get(mol, [])
        if len(hs) != 2:
            continue

        rh = 0.5 * (np.asarray(hs[0]) + np.asarray(hs[1]))
        u = rh - ro
        norm = np.linalg.norm(u)

        if norm <= 0.0:
            continue

        mols.append(mol)
        uvecs.append(u / norm)

    if not uvecs:
        raise RuntimeError("No valid water dipoles reconstructed.")

    return np.asarray(mols, dtype=int), np.asarray(uvecs, dtype=float)


def multi_origin_autocorrelation(times, mol_ids_per_frame, uvecs_per_frame):
    """
    Multi-time-origin orientational autocorrelation.

    For each lag k:
        C(k) = average over all origins i and molecules common to frames i and i+k
               of u_i(mol) dot u_{i+k}(mol)
    """
    n = len(times)
    records = []

    # Convert each frame into mol -> index mapping for fast intersection.
    maps = []
    for mols in mol_ids_per_frame:
        maps.append({int(m): j for j, m in enumerate(mols)})

    for lag in range(n):
        values = []
        n_pairs_total = 0

        for i0 in range(0, n - lag):
            i1 = i0 + lag

            mols0 = maps[i0]
            mols1 = maps[i1]
            common = sorted(set(mols0).intersection(mols1))

            if not common:
                continue

            idx0 = np.asarray([mols0[m] for m in common], dtype=int)
            idx1 = np.asarray([mols1[m] for m in common], dtype=int)

            u0 = uvecs_per_frame[i0][idx0]
            u1 = uvecs_per_frame[i1][idx1]

            dots = np.sum(u0 * u1, axis=1)
            values.append(dots)
            n_pairs_total += len(dots)

        if values:
            all_values = np.concatenate(values)
            c1 = float(np.mean(all_values))
            c1_std = float(np.std(all_values))
        else:
            c1 = np.nan
            c1_std = np.nan

        records.append({
            "lag_index": lag,
            "lag_steps": int(times[lag] - times[0]) if lag < len(times) else lag,
            "n_time_origins": n - lag,
            "n_molecule_pairs": n_pairs_total,
            "C1_mu": c1,
            "C1_mu_std": c1_std,
        })

    return pd.DataFrame(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajectory", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--oxygen-type", type=int, default=3)
    ap.add_argument("--hydrogen-type", type=int, default=4)
    ap.add_argument("--timestep-fs", type=float, default=0.10)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    times = []
    mol_ids_per_frame = []
    uvecs_per_frame = []

    for iframe, (step, frame, _bounds) in enumerate(read_lammps_dump(args.trajectory), start=1):
        print(f"[{args.label}] processing frame {iframe}, step {step}")

        mols, uvecs = compute_dipole_unit_vectors(
            frame,
            oxygen_type=args.oxygen_type,
            hydrogen_type=args.hydrogen_type,
        )

        times.append(step)
        mol_ids_per_frame.append(mols)
        uvecs_per_frame.append(uvecs)

    times = np.asarray(times, dtype=int)

    acf = multi_origin_autocorrelation(times, mol_ids_per_frame, uvecs_per_frame)

    acf["lag_time_ps"] = acf["lag_steps"] * args.timestep_fs / 1000.0

    csv = outdir / f"{args.label}_water_dipole_autocorrelation.csv"
    acf.to_csv(csv, index=False)

    plt.figure(figsize=(8, 5))
    plt.plot(acf["lag_time_ps"], acf["C1_mu"], marker="o", label="C_mu(t)")
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("Lag time, ps")
    plt.ylabel("Dipole orientational autocorrelation")
    plt.title(f"Water dipole autocorrelation: {args.label}")
    plt.legend()
    plt.tight_layout()

    png = outdir / f"{args.label}_water_dipole_autocorrelation.png"
    plt.savefig(png, dpi=200)
    plt.close()

    print("Water dipole autocorrelation analysis complete.")
    print(f"CSV: {csv}")
    print(f"PNG: {png}")


if __name__ == "__main__":
    main()
