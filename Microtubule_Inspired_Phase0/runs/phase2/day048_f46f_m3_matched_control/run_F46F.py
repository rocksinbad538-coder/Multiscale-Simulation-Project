#!/usr/bin/env python3

from pathlib import Path
import json
import sys

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


def kww(t, tau, beta):
    return np.exp(-(t/tau)**beta)


def load_water_vectors(tpr, xtc, label):
    u = mda.Universe(str(tpr), str(xtc))

    water = u.select_atoms("resname SOL")

    if len(water) == 0:
        raise RuntimeError(
            f"{label}: no SOL residues"
        )

    O = u.select_atoms("resname SOL and name OW")
    H1 = u.select_atoms("resname SOL and name HW1")
    H2 = u.select_atoms("resname SOL and name HW2")

    if not (
        len(O) > 0
        and len(O) == len(H1) == len(H2)
        and len(O) == len(water.residues)
    ):
        raise RuntimeError(
            f"{label}: canonical TIP4P/2005 "
            "OW/HW1/HW2 mapping failed"
        )

    times = np.asarray(
        [ts.time for ts in u.trajectory],
        dtype=float,
    )

    if len(times) < 40:
        raise RuntimeError(
            f"{label}: insufficient frames={len(times)}"
        )

    dt = np.diff(times)

    if not np.allclose(
        dt,
        dt[0],
        rtol=1e-6,
        atol=1e-8,
    ):
        raise RuntimeError(
            f"{label}: nonuniform frame spacing"
        )

    vectors = np.empty(
        (len(times), len(O), 3),
        dtype=np.float32,
    )

    for iframe, ts in enumerate(u.trajectory):
        box = ts.dimensions

        if box is None:
            raise RuntimeError(
                f"{label}: missing box"
            )

        rO = O.positions.astype(float)

        OH1 = minimize_vectors(
            H1.positions.astype(float) - rO,
            box,
        )

        OH2 = minimize_vectors(
            H2.positions.astype(float) - rO,
            box,
        )

        mu = 0.5 * (OH1 + OH2)

        norms = np.linalg.norm(
            mu,
            axis=1,
        )

        vectors[iframe] = (
            mu / norms[:, None]
        ).astype(np.float32)

    return {
        "label": label,
        "times": times,
        "dt_ps": float(dt[0]),
        "n_water": len(O),
        "vectors": vectors,
        "atom_names": sorted(
            set(
                water.names.tolist()
            )
        ),
    }


def molecular_corr(v):
    n = len(v)

    c1 = np.empty(n)
    c2 = np.empty(n)
    origins = np.arange(
        n,
        0,
        -1,
        dtype=int,
    )

    for lag in range(n):
        dots = np.einsum(
            "twi,twi->tw",
            v[:n-lag],
            v[lag:],
            optimize=True,
        ).astype(float)

        c1[lag] = np.mean(dots)

        c2[lag] = np.mean(
            0.5 * (
                3.0*dots*dots - 1.0
            )
        )

    return c1, c2, origins


def fit_kww(t, c, fit_max=20.0):
    mask = (
        (t <= fit_max)
        & (c > 0.0)
        & np.isfinite(c)
    )

    x = t[mask]
    y = c[mask]

    if len(x) < 8:
        raise RuntimeError(
            "Insufficient positive points for KWW fit"
        )

    p, _ = curve_fit(
        kww,
        x,
        y,
        p0=[5.0, 0.8],
        bounds=(
            [0.01, 0.10],
            [100.0, 2.0],
        ),
        maxfev=100000,
    )

    tau = float(p[0])
    beta = float(p[1])

    integrated = float(
        tau
        / beta
        * gamma(
            1.0 / beta
        )
    )

    pred = kww(
        x,
        tau,
        beta,
    )

    rss = float(
        np.sum(
            (y-pred)**2
        )
    )

    return {
        "tau_ps": tau,
        "beta": beta,
        "integrated_tau_ps": integrated,
        "rss": rss,
        "n_fit": len(x),
    }


def integral(t, c, cutoff):
    m = t <= cutoff

    return float(
        np.trapezoid(
            c[m],
            t[m],
        )
    )


def analyze(system):
    t = (
        np.arange(
            len(system["times"]),
            dtype=float,
        )
        * system["dt_ps"]
    )

    c1, c2, origins = molecular_corr(
        system["vectors"]
    )

    f1 = fit_kww(
        t,
        c1,
        20.0,
    )

    f2 = fit_kww(
        t,
        c2,
        20.0,
    )

    df = pd.DataFrame({
        "lag_ps": t,
        "n_time_origins": origins,
        "C1": c1,
        "C2": c2,
    })

    df.to_csv(
        OUT
        / f"{system['label']}_correlations.csv",
        index=False,
    )

    return {
        "label": system["label"],
        "n_frames": len(t),
        "duration_ps": float(t[-1]),
        "dt_ps": system["dt_ps"],
        "n_water": system["n_water"],
        "atom_names": system["atom_names"],
        "C1_fit": f1,
        "C2_fit": f2,
        "C1_integral_20ps": integral(
            t,
            c1,
            20.0,
        ),
        "C1_integral_40ps": integral(
            t,
            c1,
            min(
                40.0,
                t[-1],
            ),
        ),
        "C2_integral_20ps": integral(
            t,
            c2,
            20.0,
        ),
        "C2_integral_40ps": integral(
            t,
            c2,
            min(
                40.0,
                t[-1],
            ),
        ),
    }


print("Loading engineered system...")
eng = load_water_vectors(
    ENG_TPR,
    ENG_XTC,
    "engineered",
)

print("Loading matched control...")
ctrl = load_water_vectors(
    CTRL_TPR,
    CTRL_XTC,
    "control",
)

print()
print("ENGINEERED")
print(
    f"n_water={eng['n_water']} "
    f"frames={len(eng['times'])} "
    f"dt_ps={eng['dt_ps']} "
    f"duration_ps={eng['times'][-1]-eng['times'][0]}"
)

print("CONTROL")
print(
    f"n_water={ctrl['n_water']} "
    f"frames={len(ctrl['times'])} "
    f"dt_ps={ctrl['dt_ps']} "
    f"duration_ps={ctrl['times'][-1]-ctrl['times'][0]}"
)


# ------------------------------------------------------------
# Compatibility gate
# ------------------------------------------------------------

compatible_water_model = (
    eng["atom_names"]
    ==
    ctrl["atom_names"]
)

water_count_ratio = (
    ctrl["n_water"]
    /
    eng["n_water"]
)

if not compatible_water_model:
    print("CONTROL_COMPATIBILITY=FAIL")
    print(
        "REASON=WATER_ATOM_MODEL_DIFFERENT"
    )
    raise SystemExit(3)

if not (
    0.95
    <= water_count_ratio
    <= 1.05
):
    print("CONTROL_COMPATIBILITY=FAIL")
    print(
        "REASON=WATER_COUNT_DIFFERS_BY_MORE_THAN_5_PERCENT"
    )
    print(
        f"WATER_COUNT_RATIO={water_count_ratio}"
    )
    raise SystemExit(4)

print("CONTROL_WATER_MODEL_COMPATIBILITY=PASS")
print(
    f"WATER_COUNT_RATIO={water_count_ratio}"
)


# ------------------------------------------------------------
# Quantitative analysis
# ------------------------------------------------------------

eng_r = analyze(eng)
ctrl_r = analyze(ctrl)


ratio_C1 = (
    eng_r["C1_fit"]["integrated_tau_ps"]
    /
    ctrl_r["C1_fit"]["integrated_tau_ps"]
)

ratio_C2 = (
    eng_r["C2_fit"]["integrated_tau_ps"]
    /
    ctrl_r["C2_fit"]["integrated_tau_ps"]
)

enh_C1 = 100.0 * (
    ratio_C1 - 1.0
)

enh_C2 = 100.0 * (
    ratio_C2 - 1.0
)


summary = {
    "gate": "F46F",

    "engineered": eng_r,
    "control": ctrl_r,

    "ratio_tau1_engineered_over_control":
        ratio_C1,

    "percent_change_tau1":
        enh_C1,

    "ratio_tau2_engineered_over_control":
        ratio_C2,

    "percent_change_tau2":
        enh_C2,

    "interpretation_boundary": (
        "This comparison quantifies water orientational "
        "relaxation only. A positive change means slower "
        "water orientational relaxation in the engineered "
        "system; it does not by itself establish enhanced "
        "quantum coherence."
    ),
}


with (
    OUT / "F46F_summary.json"
).open("w") as fh:
    json.dump(
        summary,
        fh,
        indent=2,
    )


print()
print("============================================================")
print("F46F RESULTS")
print("============================================================")

print(
    "ENGINEERED_TAU1_PS="
    f"{eng_r['C1_fit']['integrated_tau_ps']}"
)

print(
    "CONTROL_TAU1_PS="
    f"{ctrl_r['C1_fit']['integrated_tau_ps']}"
)

print(
    "TAU1_RATIO_ENGINEERED_OVER_CONTROL="
    f"{ratio_C1}"
)

print(
    "TAU1_PERCENT_CHANGE="
    f"{enh_C1}"
)

print()

print(
    "ENGINEERED_TAU2_PS="
    f"{eng_r['C2_fit']['integrated_tau_ps']}"
)

print(
    "CONTROL_TAU2_PS="
    f"{ctrl_r['C2_fit']['integrated_tau_ps']}"
)

print(
    "TAU2_RATIO_ENGINEERED_OVER_CONTROL="
    f"{ratio_C2}"
)

print(
    "TAU2_PERCENT_CHANGE="
    f"{enh_C2}"
)

print()

print("F46F_COMPLETE=YES")
print("NEW_MD_SIMULATIONS=0")
print("NEW_QM_CALCULATIONS=0")
print("NEW_EXCITON_DYNAMICS=0")
