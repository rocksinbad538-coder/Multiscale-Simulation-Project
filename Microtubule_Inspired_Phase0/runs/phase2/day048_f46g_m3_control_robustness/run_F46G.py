#!/usr/bin/env python3

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import minimize_vectors
from scipy.optimize import curve_fit
from scipy.special import gamma


ENG_TPR = Path(sys.argv[1])
ENG_XTC = Path(sys.argv[2])
CTRL_TPR = Path(sys.argv[3])
CTRL_XTC = Path(sys.argv[4])
OUT = Path(sys.argv[5])

OUT.mkdir(parents=True, exist_ok=True)


# =====================================================================
# Models
# =====================================================================

def kww(t, tau, beta):
    return np.exp(
        -np.power(t / tau, beta)
    )


def kww_integrated_tau(tau, beta):
    return float(
        tau
        / beta
        * gamma(
            1.0 / beta
        )
    )


# =====================================================================
# System inspection
# =====================================================================

def topology_summary(tpr, label):

    u = mda.Universe(str(tpr))

    residue_counts = {}

    for resname in sorted(
        set(u.residues.resnames.tolist())
    ):
        residue_counts[resname] = int(
            np.sum(
                u.residues.resnames == resname
            )
        )

    atom_counts = {}

    for name in sorted(
        set(u.atoms.names.tolist())
    ):
        atom_counts[name] = int(
            np.sum(
                u.atoms.names == name
            )
        )

    summary = {
        "label": label,
        "n_atoms": int(len(u.atoms)),
        "n_residues": int(len(u.residues)),
        "residue_counts": residue_counts,
        "atom_name_counts": atom_counts,
    }

    try:
        charges = np.asarray(
            u.atoms.charges,
            dtype=float,
        )

        summary["total_charge_e"] = float(
            np.sum(charges)
        )

        summary["sum_abs_charge_e"] = float(
            np.sum(
                np.abs(charges)
            )
        )

    except Exception:
        summary["total_charge_e"] = None
        summary["sum_abs_charge_e"] = None

    return summary


# =====================================================================
# Trajectory loading
# =====================================================================

def load_water_vectors(tpr, xtc, label):

    u = mda.Universe(
        str(tpr),
        str(xtc),
    )

    water = u.select_atoms(
        "resname SOL"
    )

    O = u.select_atoms(
        "resname SOL and name OW"
    )

    H1 = u.select_atoms(
        "resname SOL and name HW1"
    )

    H2 = u.select_atoms(
        "resname SOL and name HW2"
    )

    if not (
        len(O) > 0
        and len(O) == len(H1)
        and len(O) == len(H2)
        and len(O) == len(water.residues)
    ):
        raise RuntimeError(
            f"{label}: TIP4P/2005 mapping failed"
        )

    times = np.asarray(
        [ts.time for ts in u.trajectory],
        dtype=float,
    )

    dts = np.diff(times)

    if len(dts) == 0:
        raise RuntimeError(
            f"{label}: insufficient frames"
        )

    if not np.allclose(
        dts,
        dts[0],
        rtol=1e-7,
        atol=1e-8,
    ):
        raise RuntimeError(
            f"{label}: nonuniform frame spacing"
        )

    vectors = np.empty(
        (
            len(times),
            len(O),
            3,
        ),
        dtype=np.float32,
    )

    for iframe, ts in enumerate(u.trajectory):

        box = ts.dimensions

        rO = O.positions.astype(
            np.float64
        )

        OH1 = minimize_vectors(
            H1.positions.astype(
                np.float64
            ) - rO,
            box,
        )

        OH2 = minimize_vectors(
            H2.positions.astype(
                np.float64
            ) - rO,
            box,
        )

        mu = 0.5 * (
            OH1 + OH2
        )

        norm = np.linalg.norm(
            mu,
            axis=1,
        )

        if np.any(norm <= 0.0):
            raise RuntimeError(
                f"{label}: zero-length dipole"
            )

        vectors[iframe] = (
            mu / norm[:, None]
        ).astype(np.float32)

    return {
        "label": label,
        "times": times,
        "dt_ps": float(dts[0]),
        "vectors": vectors,
        "n_water": int(len(O)),
    }


# =====================================================================
# Correlations
# =====================================================================

def molecular_corr(v):

    n = len(v)

    C1 = np.empty(
        n,
        dtype=float,
    )

    C2 = np.empty(
        n,
        dtype=float,
    )

    for lag in range(n):

        dots = np.einsum(
            "twi,twi->tw",
            v[:n-lag],
            v[lag:],
            optimize=True,
        ).astype(
            np.float64
        )

        C1[lag] = float(
            np.mean(dots)
        )

        C2[lag] = float(
            np.mean(
                0.5 * (
                    3.0 * dots*dots - 1.0
                )
            )
        )

    return C1, C2


def fit_kww(t, c, fit_max_ps=10.0):

    mask = (
        (t >= 0.0)
        & (t <= fit_max_ps)
        & np.isfinite(c)
        & (c > 0.0)
    )

    x = t[mask]
    y = c[mask]

    if len(x) < 8:
        raise RuntimeError(
            "Insufficient points for block KWW fit"
        )

    p, _ = curve_fit(
        kww,
        x,
        y,
        p0=[
            5.0,
            0.8,
        ],
        bounds=(
            [
                0.01,
                0.10,
            ],
            [
                100.0,
                2.0,
            ],
        ),
        maxfev=100000,
    )

    tau = float(p[0])
    beta = float(p[1])

    return {
        "tau_kww_ps":
            tau,

        "beta":
            beta,

        "tau_integrated_ps":
            kww_integrated_tau(
                tau,
                beta,
            ),
    }


def integrate_to(
    t,
    c,
    cutoff_ps,
):

    mask = t <= cutoff_ps

    return float(
        np.trapezoid(
            c[mask],
            t[mask],
        )
    )


# =====================================================================
# Matched 25 ps blocks
# =====================================================================

def matched_blocks(system):

    dt = system["dt_ps"]

    target_block_ps = 25.0

    frames_per_block = int(
        round(
            target_block_ps / dt
        )
    ) + 1

    # Explicit matched windows:
    # [0,25], [25,50], [50,75], [75,100]
    starts_ps = [
        0.0,
        25.0,
        50.0,
        75.0,
    ]

    rows = []

    for ib, start_ps in enumerate(
        starts_ps,
        start=1,
    ):

        end_ps = (
            start_ps
            + target_block_ps
        )

        times = system["times"]

        mask = (
            (times >= times[0] + start_ps - 1e-8)
            &
            (times <= times[0] + end_ps + 1e-8)
        )

        idx = np.where(
            mask
        )[0]

        if len(idx) != frames_per_block:
            raise RuntimeError(
                f"{system['label']} block {ib}: "
                f"expected {frames_per_block} frames, "
                f"found {len(idx)}"
            )

        v = system[
            "vectors"
        ][idx]

        C1, C2 = molecular_corr(
            v
        )

        lag = (
            np.arange(
                len(v),
                dtype=float,
            )
            * dt
        )

        f1 = fit_kww(
            lag,
            C1,
            fit_max_ps=10.0,
        )

        f2 = fit_kww(
            lag,
            C2,
            fit_max_ps=10.0,
        )

        rows.append(
            {
                "system":
                    system["label"],

                "block":
                    ib,

                "start_ps":
                    start_ps,

                "end_ps":
                    end_ps,

                "C1_tau_integrated_ps":
                    f1[
                        "tau_integrated_ps"
                    ],

                "C1_tau_kww_ps":
                    f1[
                        "tau_kww_ps"
                    ],

                "C1_beta":
                    f1[
                        "beta"
                    ],

                "C1_integral_10ps":
                    integrate_to(
                        lag,
                        C1,
                        10.0,
                    ),

                "C2_tau_integrated_ps":
                    f2[
                        "tau_integrated_ps"
                    ],

                "C2_tau_kww_ps":
                    f2[
                        "tau_kww_ps"
                    ],

                "C2_beta":
                    f2[
                        "beta"
                    ],

                "C2_integral_10ps":
                    integrate_to(
                        lag,
                        C2,
                        10.0,
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# =====================================================================
# Statistics
# =====================================================================

def stats(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    return {
        "mean":
            float(
                np.mean(x)
            ),

        "sd":
            float(
                np.std(
                    x,
                    ddof=1,
                )
            ),

        "sem":
            float(
                np.std(
                    x,
                    ddof=1,
                )
                / np.sqrt(
                    len(x)
                )
            ),

        "min":
            float(
                np.min(x)
            ),

        "max":
            float(
                np.max(x)
            ),
    }


def hedges_g(x, y):

    x = np.asarray(
        x,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )

    nx = len(x)
    ny = len(y)

    sx2 = np.var(
        x,
        ddof=1,
    )

    sy2 = np.var(
        y,
        ddof=1,
    )

    sp = np.sqrt(
        (
            (nx-1)*sx2
            +
            (ny-1)*sy2
        )
        /
        (
            nx+ny-2
        )
    )

    if sp == 0:
        return np.nan

    d = (
        np.mean(x)
        -
        np.mean(y)
    ) / sp

    correction = (
        1.0
        -
        3.0
        /
        (
            4.0*(nx+ny)
            - 9.0
        )
    )

    return float(
        correction * d
    )


def exact_permutation_p(
    x,
    y,
):

    x = np.asarray(
        x,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )

    pooled = np.concatenate(
        [
            x,
            y,
        ]
    )

    nx = len(x)

    observed = abs(
        np.mean(x)
        -
        np.mean(y)
    )

    diffs = []

    for comb in itertools.combinations(
        range(len(pooled)),
        nx,
    ):

        idx_x = np.asarray(
            comb,
            dtype=int,
        )

        mask = np.ones(
            len(pooled),
            dtype=bool,
        )

        mask[idx_x] = False

        a = pooled[idx_x]
        b = pooled[mask]

        diffs.append(
            abs(
                np.mean(a)
                -
                np.mean(b)
            )
        )

    diffs = np.asarray(
        diffs,
        dtype=float,
    )

    return float(
        np.mean(
            diffs
            >=
            observed
            - 1e-12
        )
    )


def compare_metric(
    eng_df,
    ctrl_df,
    metric,
):

    x = eng_df[
        metric
    ].to_numpy(
        dtype=float
    )

    y = ctrl_df[
        metric
    ].to_numpy(
        dtype=float
    )

    eng_mean = float(
        np.mean(x)
    )

    ctrl_mean = float(
        np.mean(y)
    )

    delta = (
        eng_mean
        -
        ctrl_mean
    )

    ratio = (
        eng_mean
        /
        ctrl_mean
    )

    return {
        "engineered":
            stats(x),

        "control":
            stats(y),

        "difference_ps":
            float(delta),

        "ratio_engineered_over_control":
            float(ratio),

        "percent_change":
            float(
                100.0
                * (
                    ratio - 1.0
                )
            ),

        "hedges_g":
            hedges_g(
                x,
                y,
            ),

        "exact_two_sided_permutation_p":
            exact_permutation_p(
                x,
                y,
            ),
    }


# =====================================================================
# Execute
# =====================================================================

print(
    "Inspecting topology composition..."
)

eng_top = topology_summary(
    ENG_TPR,
    "engineered",
)

ctrl_top = topology_summary(
    CTRL_TPR,
    "control",
)


print()
print("ENGINEERED_TOPOLOGY")
print(
    json.dumps(
        eng_top,
        indent=2,
    )
)

print()
print("CONTROL_TOPOLOGY")
print(
    json.dumps(
        ctrl_top,
        indent=2,
    )
)


print()
print(
    "Loading trajectories and reconstructing water dipoles..."
)

eng = load_water_vectors(
    ENG_TPR,
    ENG_XTC,
    "engineered",
)

ctrl = load_water_vectors(
    CTRL_TPR,
    CTRL_XTC,
    "control",
)


print(
    f"ENGINEERED_FRAMES={len(eng['times'])}"
)

print(
    f"CONTROL_FRAMES={len(ctrl['times'])}"
)

print(
    f"ENGINEERED_DURATION_PS="
    f"{eng['times'][-1]-eng['times'][0]}"
)

print(
    f"CONTROL_DURATION_PS="
    f"{ctrl['times'][-1]-ctrl['times'][0]}"
)

print(
    f"FRAME_SPACING_ENGINEERED_PS="
    f"{eng['dt_ps']}"
)

print(
    f"FRAME_SPACING_CONTROL_PS="
    f"{ctrl['dt_ps']}"
)


eng_blocks = matched_blocks(
    eng
)

ctrl_blocks = matched_blocks(
    ctrl
)

all_blocks = pd.concat(
    [
        eng_blocks,
        ctrl_blocks,
    ],
    ignore_index=True,
)

all_blocks.to_csv(
    OUT
    / "F46G_matched_25ps_block_results.csv",
    index=False,
)


metrics = [
    "C1_tau_integrated_ps",
    "C1_integral_10ps",
    "C2_tau_integrated_ps",
    "C2_integral_10ps",
]

comparisons = {
    metric:
        compare_metric(
            eng_blocks,
            ctrl_blocks,
            metric,
        )
    for metric in metrics
}


summary = {
    "gate":
        "F46G",

    "topology_engineered":
        eng_top,

    "topology_control":
        ctrl_top,

    "matched_block_duration_ps":
        25.0,

    "n_blocks_per_system":
        4,

    "comparisons":
        comparisons,

    "interpretation_boundary":
        (
            "Permutation tests use only four contiguous "
            "25-ps blocks per trajectory and therefore "
            "provide a robustness diagnostic rather than "
            "a definitive independent-replica significance "
            "test. Causal interpretation additionally "
            "requires the matched-control provenance to "
            "establish exactly which physical components "
            "differ between systems."
        ),
}


with (
    OUT
    / "F46G_summary.json"
).open("w") as fh:

    json.dump(
        summary,
        fh,
        indent=2,
    )


print()
print(
    "============================================================"
)

print(
    "MATCHED 25-PS BLOCK RESULTS"
)

print(
    "============================================================"
)

print(
    all_blocks.to_string(
        index=False
    )
)


print()
print(
    "============================================================"
)

print(
    "ENGINEERED vs CONTROL"
)

print(
    "============================================================"
)

for metric, result in comparisons.items():

    print()
    print(
        f"METRIC={metric}"
    )

    print(
        f"ENGINEERED_MEAN="
        f"{result['engineered']['mean']}"
    )

    print(
        f"ENGINEERED_SD="
        f"{result['engineered']['sd']}"
    )

    print(
        f"CONTROL_MEAN="
        f"{result['control']['mean']}"
    )

    print(
        f"CONTROL_SD="
        f"{result['control']['sd']}"
    )

    print(
        f"DIFFERENCE="
        f"{result['difference_ps']}"
    )

    print(
        f"PERCENT_CHANGE="
        f"{result['percent_change']}"
    )

    print(
        f"HEDGES_G="
        f"{result['hedges_g']}"
    )

    print(
        f"EXACT_PERMUTATION_P="
        f"{result['exact_two_sided_permutation_p']}"
    )


print()
print(
    "SUMMARY_FILE="
    + str(
        OUT
        / "F46G_summary.json"
    )
)

print(
    "F46G_COMPLETE=YES"
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
