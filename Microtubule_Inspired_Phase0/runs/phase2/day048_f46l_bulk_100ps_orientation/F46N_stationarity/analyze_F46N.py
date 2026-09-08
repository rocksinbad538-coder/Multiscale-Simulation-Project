from pathlib import Path
import sys

import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import minimize_vectors
from scipy.optimize import curve_fit
from scipy.special import gamma


OUT = Path(sys.argv[1])

systems = {
    "engineered": (
        Path(sys.argv[2]),
        Path(sys.argv[3]),
    ),
    "bulk": (
        Path(sys.argv[4]),
        Path(sys.argv[5]),
    ),
}


def kww(t, tau, beta):
    return np.exp(-(t/tau)**beta)


def load_vectors(tpr, xtc):

    u = mda.Universe(
        str(tpr),
        str(xtc),
    )

    waters = u.residues[
        u.residues.resnames == "SOL"
    ]

    ow = []
    h1 = []
    h2 = []

    for res in waters:

        amap = {
            a.name: a.index
            for a in res.atoms
        }

        ow.append(amap["OW"])
        h1.append(amap["HW1"])
        h2.append(amap["HW2"])

    ow = np.asarray(ow)
    h1 = np.asarray(h1)
    h2 = np.asarray(h2)

    times = []
    vectors = []

    for ts in u.trajectory:

        p = u.atoms.positions
        box = ts.dimensions

        rO = p[ow]

        v1 = minimize_vectors(
            p[h1]-rO,
            box,
        )

        v2 = minimize_vectors(
            p[h2]-rO,
            box,
        )

        mu = 0.5*(v1+v2)

        mu /= np.linalg.norm(
            mu,
            axis=1,
        )[:,None]

        times.append(
            float(ts.time)
        )

        vectors.append(
            mu.astype(np.float32)
        )

    t = np.asarray(times)
    v = np.asarray(vectors)

    t -= t[0]

    return t, v


def correlations(v):

    nf = len(v)

    c1 = np.empty(nf)
    c2 = np.empty(nf)

    for lag in range(nf):

        dots = np.sum(
            v[:nf-lag]
            * v[lag:],
            axis=2,
            dtype=np.float64,
        )

        c1[lag] = np.mean(dots)

        c2[lag] = np.mean(
            0.5*(3*dots*dots-1)
        )

    return c1, c2


def fit_corr(t, c, maxfit):

    mask = (
        (t <= maxfit)
        & (c > 0)
        & np.isfinite(c)
    )

    x = t[mask]
    y = c[mask]

    p, _ = curve_fit(
        kww,
        x,
        y,
        p0=(3.0,0.8),
        bounds=(
            (1e-5,0.05),
            (1000,2),
        ),
        maxfev=50000,
    )

    tau = float(p[0])
    beta = float(p[1])

    tint = float(
        tau/beta
        * gamma(1/beta)
    )

    return tau, beta, tint


def integral(t,c,cutoff):

    m = t <= cutoff

    return float(
        np.trapezoid(
            c[m],
            t[m],
        )
    )


windows = [
    (0,100),
    (25,100),
    (50,100),
    (75,100),
]

rows = []


for system, (
    tpr,
    xtc,
) in systems.items():

    print(
        f"LOADING_{system.upper()}=YES"
    )

    time, vectors = load_vectors(
        tpr,
        xtc,
    )

    for start, end in windows:

        mask = (
            (time >= start)
            & (time <= end)
        )

        t = (
            time[mask]
            - time[mask][0]
        )

        v = vectors[mask]

        c1,c2 = correlations(v)

        duration = (
            end-start
        )

        # Fit only a conservative fraction
        # of each available window.
        maxfit = min(
            20.0,
            max(
                10.0,
                0.4*duration
            )
        )

        t1,b1,i1 = fit_corr(
            t,
            c1,
            maxfit,
        )

        t2,b2,i2 = fit_corr(
            t,
            c2,
            maxfit,
        )

        cutoff = min(
            10.0,
            duration
        )

        rows.append({
            "system":system,
            "start_ps":start,
            "end_ps":end,
            "duration_ps":duration,
            "fit_limit_ps":maxfit,

            "C1_tau_kww_ps":t1,
            "C1_beta":b1,
            "C1_tau_integrated_ps":i1,
            "C1_integral_10ps":
                integral(
                    t,c1,cutoff
                ),

            "C2_tau_kww_ps":t2,
            "C2_beta":b2,
            "C2_tau_integrated_ps":i2,
            "C2_integral_10ps":
                integral(
                    t,c2,cutoff
                ),
        })


df = pd.DataFrame(rows)

df.to_csv(
    OUT/"F46N_window_results.csv",
    index=False,
)


print()
print("="*72)
print("F46N WINDOW RESULTS")
print("="*72)

print(
    df.to_string(index=False)
)


print()
print("="*72)
print("ENGINEERED / BULK RATIOS BY WINDOW")
print("="*72)

metrics = [
    "C1_tau_integrated_ps",
    "C1_integral_10ps",
    "C2_tau_integrated_ps",
    "C2_integral_10ps",
]

comparison_rows=[]


for start,end in windows:

    e = df[
        (df.system=="engineered")
        & (df.start_ps==start)
        & (df.end_ps==end)
    ].iloc[0]

    b = df[
        (df.system=="bulk")
        & (df.start_ps==start)
        & (df.end_ps==end)
    ].iloc[0]

    print()
    print(
        f"WINDOW_PS={start}-{end}"
    )

    for metric in metrics:

        ratio = (
            e[metric]
            / b[metric]
        )

        pct = (
            100*(ratio-1)
        )

        print(
            f"{metric}_RATIO="
            f"{ratio:.12f}"
        )

        print(
            f"{metric}_PERCENT="
            f"{pct:.6f}"
        )

        comparison_rows.append({
            "start_ps":start,
            "end_ps":end,
            "metric":metric,
            "engineered":
                float(e[metric]),
            "bulk":
                float(b[metric]),
            "ratio":
                float(ratio),
            "percent_change":
                float(pct),
        })


pd.DataFrame(
    comparison_rows
).to_csv(
    OUT/"F46N_window_comparison.csv",
    index=False,
)


print()
print("="*72)
print("STATIONARITY DIAGNOSTIC")
print("="*72)

for metric in metrics:

    sub = pd.DataFrame(
        comparison_rows
    )

    sub = sub[
        sub.metric==metric
    ]

    vals = (
        sub.sort_values(
            "start_ps"
        )
        ["percent_change"]
        .to_numpy()
    )

    print(
        f"{metric}_PERCENT_SEQUENCE="
        + ",".join(
            f"{x:.6f}"
            for x in vals
        )
    )

    monotonic_decline = (
        np.all(
            np.diff(vals) <= 1e-10
        )
    )

    print(
        f"{metric}_MONOTONIC_DECLINE="
        f"{'YES' if monotonic_decline else 'NO'}"
    )


print()
print(
    "F46N_ANALYSIS_COMPLETE=YES"
)
