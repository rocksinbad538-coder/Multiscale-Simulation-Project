#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors


if len(sys.argv) != 4:
    raise SystemExit(
        "usage: run_F46D.py topology.tpr trajectory.xtc output_dir"
    )

topology = Path(sys.argv[1])
trajectory = Path(sys.argv[2])
outdir = Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Load trajectory
# ---------------------------------------------------------------------

u = mda.Universe(str(topology), str(trajectory))

times = np.asarray(
    [ts.time for ts in u.trajectory],
    dtype=float,
)

if len(times) < 20:
    raise RuntimeError(
        f"Insufficient trajectory frames: {len(times)}"
    )

dt = np.diff(times)

if not np.allclose(
    dt,
    dt[0],
    rtol=1e-6,
    atol=1e-8,
):
    raise RuntimeError(
        "Trajectory frame times are not uniformly sampled."
    )

dt_ps = float(dt[0])

print(f"N_FRAMES={len(times)}")
print(f"FIRST_TIME_PS={times[0]:.9f}")
print(f"LAST_TIME_PS={times[-1]:.9f}")
print(f"FRAME_SPACING_PS={dt_ps:.9f}")


# ---------------------------------------------------------------------
# Water topology detection
#
# Canonical expected TIP4P/2005:
#   resname SOL
#   OW HW1 HW2 MW
#
# We intentionally build orientation from OW -> midpoint(HW1,HW2).
# The virtual M site is not needed for normalized orientational dynamics.
# ---------------------------------------------------------------------

water = u.select_atoms("resname SOL")

if len(water) == 0:
    raise RuntimeError(
        "No atoms with resname SOL were found."
    )

unique_names = sorted(set(water.names.tolist()))

print("WATER_ATOM_NAMES=" + ",".join(unique_names))
print(f"N_WATER_ATOMS={len(water)}")
print(f"N_WATER_RESIDUES={len(water.residues)}")

oxygen = u.select_atoms("resname SOL and name OW")
h1 = u.select_atoms("resname SOL and name HW1")
h2 = u.select_atoms("resname SOL and name HW2")

if not (
    len(oxygen) > 0
    and len(oxygen) == len(h1) == len(h2)
):
    print(
        "Canonical OW/HW1/HW2 selection failed; "
        "attempting residue-based detection."
    )

    o_idx = []
    h1_idx = []
    h2_idx = []

    for res in water.residues:
        O = [
            a for a in res.atoms
            if a.name.upper().startswith("O")
        ]

        H = [
            a for a in res.atoms
            if a.name.upper().startswith("H")
        ]

        if len(O) == 1 and len(H) == 2:
            o_idx.append(O[0].index)
            h1_idx.append(H[0].index)
            h2_idx.append(H[1].index)

    if not o_idx:
        raise RuntimeError(
            "Could not identify water O/H atoms."
        )

    oxygen = u.atoms[np.asarray(o_idx)]
    h1 = u.atoms[np.asarray(h1_idx)]
    h2 = u.atoms[np.asarray(h2_idx)]


nwater = len(oxygen)

if nwater != len(water.residues):
    raise RuntimeError(
        f"Water mapping incomplete: "
        f"{nwater} mapped vs "
        f"{len(water.residues)} residues"
    )

print(f"N_WATER={nwater}")


# ---------------------------------------------------------------------
# Build molecular dipole unit vectors.
#
# PBC treatment:
# minimum-image vectors H-O are reconstructed independently for every
# frame using the instantaneous simulation cell.
# ---------------------------------------------------------------------

vectors = np.empty(
    (len(times), nwater, 3),
    dtype=np.float32,
)

collective = np.empty(
    (len(times), 3),
    dtype=np.float64,
)

for iframe, ts in enumerate(u.trajectory):

    box = ts.dimensions

    if box is None or len(box) < 6:
        raise RuntimeError(
            "Missing periodic-box information."
        )

    rO = oxygen.positions.astype(np.float64)
    rH1 = h1.positions.astype(np.float64)
    rH2 = h2.positions.astype(np.float64)

    OH1 = minimize_vectors(
        rH1 - rO,
        box,
    )

    OH2 = minimize_vectors(
        rH2 - rO,
        box,
    )

    mu = 0.5 * (OH1 + OH2)

    norm = np.linalg.norm(
        mu,
        axis=1,
    )

    if np.any(norm <= 0.0):
        raise RuntimeError(
            f"Zero dipole vector at frame {iframe}"
        )

    uvec = mu / norm[:, None]

    vectors[iframe] = uvec.astype(np.float32)

    collective[iframe] = np.mean(
        uvec,
        axis=0,
    )

    if (
        iframe == 0
        or iframe == len(times) - 1
        or iframe % 25 == 0
    ):
        print(
            f"FRAME={iframe:4d}/"
            f"{len(times)-1:4d} "
            f"TIME_PS={ts.time:.3f}"
        )


# ---------------------------------------------------------------------
# Correlations with all available time origins.
# ---------------------------------------------------------------------

def molecular_correlations(v: np.ndarray):
    nframe = v.shape[0]

    C1 = np.empty(nframe)
    C2 = np.empty(nframe)
    origins = np.empty(
        nframe,
        dtype=int,
    )

    for lag in range(nframe):

        left = v[: nframe-lag]
        right = v[lag:]

        dots = np.einsum(
            "twi,twi->tw",
            left,
            right,
            optimize=True,
        ).astype(np.float64)

        C1[lag] = np.mean(dots)

        C2[lag] = np.mean(
            0.5 * (
                3.0 * dots*dots - 1.0
            )
        )

        origins[lag] = nframe - lag

    return C1, C2, origins


def vector_autocorrelation(
    x: np.ndarray,
    subtract_mean: bool = True,
):
    y = np.asarray(
        x,
        dtype=np.float64,
    )

    if subtract_mean:
        y = y - np.mean(
            y,
            axis=0,
            keepdims=True,
        )

    denominator = np.mean(
        np.sum(y*y, axis=1)
    )

    if denominator <= 0.0:
        raise RuntimeError(
            "Collective correlation denominator <= 0"
        )

    nframe = len(y)
    out = np.empty(nframe)

    for lag in range(nframe):

        left = y[: nframe-lag]
        right = y[lag:]

        out[lag] = (
            np.mean(
                np.einsum(
                    "ti,ti->t",
                    left,
                    right,
                )
            )
            / denominator
        )

    return out


C1, C2, origins = molecular_correlations(
    vectors
)

CP = vector_autocorrelation(
    collective,
    subtract_mean=True,
)


# ---------------------------------------------------------------------
# Block analysis
#
# Four equal contiguous blocks provide a first reproducibility estimate.
# This is not treated as a final uncertainty model.
# ---------------------------------------------------------------------

nblocks = 4
block_edges = np.linspace(
    0,
    len(times),
    nblocks + 1,
    dtype=int,
)

block_rows = []
block_curves = []

for ib in range(nblocks):

    a = int(block_edges[ib])
    b = int(block_edges[ib+1])

    vb = vectors[a:b]
    pb = collective[a:b]

    if len(vb) < 10:
        continue

    c1b, c2b, ob = molecular_correlations(
        vb
    )

    cpb = vector_autocorrelation(
        pb,
        subtract_mean=True,
    )

    tb = np.arange(
        len(vb),
        dtype=float,
    ) * dt_ps

    for j in range(len(tb)):
        block_curves.append(
            {
                "block": ib + 1,
                "lag_ps": tb[j],
                "C1": c1b[j],
                "C2": c2b[j],
                "C_collective": cpb[j],
                "n_time_origins": int(ob[j]),
            }
        )

    block_rows.append(
        {
            "block": ib + 1,
            "start_time_ps": float(times[a]),
            "end_time_ps": float(times[b-1]),
            "n_frames": int(b-a),
        }
    )


# ---------------------------------------------------------------------
# Conservative numerical integral.
#
# Integral only until first zero crossing OR 20 ps, whichever is earlier.
# It is explicitly labelled finite-window tau_int, not a converged
# microscopic relaxation lifetime.
# ---------------------------------------------------------------------

def finite_integral(
    time_ps,
    corr,
    cap_ps=20.0,
):
    mask = time_ps <= cap_ps

    t = time_ps[mask]
    c = corr[mask]

    if len(t) < 2:
        return np.nan, np.nan

    nonpos = np.where(
        c[1:] <= 0.0
    )[0]

    if len(nonpos):
        stop = int(
            nonpos[0] + 2
        )
        t = t[:stop]
        c = c[:stop]

    value = float(
        np.trapezoid(
            c,
            t,
        )
    )

    return value, float(t[-1])


lags = np.arange(
    len(times),
    dtype=float,
) * dt_ps

tau_C1, win_C1 = finite_integral(
    lags,
    C1,
)

tau_C2, win_C2 = finite_integral(
    lags,
    C2,
)

tau_CP, win_CP = finite_integral(
    lags,
    CP,
)


# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

df = pd.DataFrame(
    {
        "lag_ps": lags,
        "n_time_origins": origins,
        "C1_molecular": C1,
        "C2_molecular": C2,
        "C_collective_water": CP,
    }
)

df.to_csv(
    outdir / "F46D_water_dipole_correlations.csv",
    index=False,
)

pd.DataFrame(
    block_curves
).to_csv(
    outdir / "F46D_block_correlations.csv",
    index=False,
)

pd.DataFrame(
    block_rows
).to_csv(
    outdir / "F46D_block_manifest.csv",
    index=False,
)

collective_df = pd.DataFrame(
    {
        "time_ps": times,
        "Pux": collective[:,0],
        "Puy": collective[:,1],
        "Puz": collective[:,2],
        "P_magnitude": np.linalg.norm(
            collective,
            axis=1,
        ),
    }
)

collective_df.to_csv(
    outdir / "F46D_collective_orientation_timeseries.csv",
    index=False,
)

summary = {
    "gate": "F46D",
    "temperature_K": 300.0,
    "trajectory": str(trajectory),
    "topology": str(topology),

    "n_frames": int(len(times)),
    "n_water": int(nwater),

    "first_time_ps": float(times[0]),
    "last_time_ps": float(times[-1]),
    "frame_spacing_ps": dt_ps,

    "sampling_nyquist_ps_inv": float(
        1.0 / (2.0 * dt_ps)
    ),

    "sampling_nyquist_THz": float(
        1.0 / (2.0 * dt_ps)
    ),

    "C1_final": float(C1[-1]),
    "C2_final": float(C2[-1]),
    "C_collective_final": float(CP[-1]),

    "finite_window_tau_C1_ps": tau_C1,
    "finite_window_tau_C1_limit_ps": win_C1,

    "finite_window_tau_C2_ps": tau_C2,
    "finite_window_tau_C2_limit_ps": win_C2,

    "finite_window_tau_collective_ps": tau_CP,
    "finite_window_tau_collective_limit_ps": win_CP,

    "mean_collective_orientation_x":
        float(np.mean(collective[:,0])),

    "mean_collective_orientation_y":
        float(np.mean(collective[:,1])),

    "mean_collective_orientation_z":
        float(np.mean(collective[:,2])),

    "mean_collective_orientation_magnitude":
        float(
            np.mean(
                np.linalg.norm(
                    collective,
                    axis=1,
                )
            )
        ),

    "interpretation_boundary":
        (
            "Normalized water orientational and collective "
            "relaxation only. This is not yet epsilon(omega), "
            "a dielectric constant, or a microscopic excitonic "
            "bath spectral density."
        ),
}

with (
    outdir
    / "F46D_summary.json"
).open("w") as fh:
    json.dump(
        summary,
        fh,
        indent=2,
    )

print()
print("================================================================")
print("F46D ANALYSIS COMPLETE")
print("================================================================")

for k, v in summary.items():
    print(f"{k}={v}")

print()
print(
    "OUTPUT="
    + str(
        outdir
        / "F46D_water_dipole_correlations.csv"
    )
)

print(
    "SUMMARY="
    + str(
        outdir
        / "F46D_summary.json"
    )
)

print(
    "NEW_MD_SIMULATIONS=0"
)

print(
    "NEW_QM_CALCULATIONS=0"
)

print(
    "NEW_EXCITON_DYNAMICS=0"
)

print(
    "ACTUAL_M3_DATA_ANALYSIS=YES"
)
