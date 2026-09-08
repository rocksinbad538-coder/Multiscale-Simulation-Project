#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.optimize import curve_fit
from scipy.special import gamma


IN = Path(
    "runs/phase2/"
    "day048_f46d_m3_atomistic_dipole_relaxation"
)

OUT = Path(
    "runs/phase2/"
    "day048_f46e_m3_relaxation_quantification"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ================================================================
# Utilities
# ================================================================

def monoexp(t, tau):
    return np.exp(-t / tau)


def kww(t, tau, beta):
    return np.exp(
        -np.power(t / tau, beta)
    )


def aic(n, rss, k):
    rss = max(
        float(rss),
        np.finfo(float).tiny,
    )

    return (
        n * np.log(rss / n)
        + 2 * k
    )


def integrate_to(t, c, cutoff):
    mask = t <= cutoff

    if np.count_nonzero(mask) < 2:
        return np.nan

    return float(
        np.trapezoid(
            c[mask],
            t[mask],
        )
    )


def first_zero_integral(t, c):
    idx = np.where(c[1:] <= 0.0)[0]

    if len(idx):
        stop = int(idx[0] + 2)
    else:
        stop = len(c)

    tt = t[:stop]
    cc = c[:stop]

    return (
        float(
            np.trapezoid(
                cc,
                tt,
            )
        ),
        float(tt[-1]),
        int(stop),
    )


def fit_models(t, c, fit_max=20.0):

    mask = (
        (t >= 0.0)
        & (t <= fit_max)
        & np.isfinite(c)
        & (c > 0.0)
    )

    x = t[mask]
    y = c[mask]

    if len(x) < 8:
        raise RuntimeError(
            "Too few points for fitting."
        )

    # ------------------------------------------------------------
    # Monoexponential
    # ------------------------------------------------------------

    p_mono, _ = curve_fit(
        monoexp,
        x,
        y,
        p0=[5.0],
        bounds=(
            [0.01],
            [100.0],
        ),
        maxfev=100000,
    )

    pred_mono = monoexp(
        x,
        *p_mono,
    )

    rss_mono = float(
        np.sum(
            (y - pred_mono)**2
        )
    )

    # ------------------------------------------------------------
    # Kohlrausch-Williams-Watts
    # ------------------------------------------------------------

    p_kww, _ = curve_fit(
        kww,
        x,
        y,
        p0=[5.0, 0.8],
        bounds=(
            [0.01, 0.10],
            [100.0, 2.00],
        ),
        maxfev=100000,
    )

    pred_kww = kww(
        x,
        *p_kww,
    )

    rss_kww = float(
        np.sum(
            (y - pred_kww)**2
        )
    )

    tau_kww = float(
        p_kww[0]
    )

    beta_kww = float(
        p_kww[1]
    )

    integral_kww = float(
        tau_kww
        / beta_kww
        * gamma(
            1.0 / beta_kww
        )
    )

    return {
        "n_fit": int(len(x)),

        "mono_tau_ps":
            float(p_mono[0]),

        "mono_rss":
            rss_mono,

        "mono_aic":
            float(
                aic(
                    len(x),
                    rss_mono,
                    1,
                )
            ),

        "kww_tau_ps":
            tau_kww,

        "kww_beta":
            beta_kww,

        "kww_integrated_tau_ps":
            integral_kww,

        "kww_rss":
            rss_kww,

        "kww_aic":
            float(
                aic(
                    len(x),
                    rss_kww,
                    2,
                )
            ),

        "delta_AIC_mono_minus_kww":
            float(
                aic(
                    len(x),
                    rss_mono,
                    1,
                )
                -
                aic(
                    len(x),
                    rss_kww,
                    2,
                )
            ),
    }


# ================================================================
# Full trajectory
# ================================================================

df = pd.read_csv(
    IN
    / "F46D_water_dipole_correlations.csv"
)

t = df["lag_ps"].to_numpy(
    dtype=float
)

C1 = df["C1_molecular"].to_numpy(
    dtype=float
)

C2 = df["C2_molecular"].to_numpy(
    dtype=float
)

CP = df["C_collective_water"].to_numpy(
    dtype=float
)

origins = df["n_time_origins"].to_numpy(
    dtype=int
)


fit_C1 = fit_models(
    t,
    C1,
    fit_max=20.0,
)

fit_C2 = fit_models(
    t,
    C2,
    fit_max=20.0,
)


# ================================================================
# Integral convergence
# ================================================================

cutoffs = [
    2.0,
    5.0,
    10.0,
    15.0,
    20.0,
    30.0,
    40.0,
    50.0,
]

conv_rows = []

for cutoff in cutoffs:
    conv_rows.append(
        {
            "cutoff_ps": cutoff,

            "integral_C1_ps":
                integrate_to(
                    t,
                    C1,
                    cutoff,
                ),

            "integral_C2_ps":
                integrate_to(
                    t,
                    C2,
                    cutoff,
                ),

            "C1_at_cutoff":
                float(
                    C1[
                        np.argmin(
                            np.abs(
                                t-cutoff
                            )
                        )
                    ]
                ),

            "C2_at_cutoff":
                float(
                    C2[
                        np.argmin(
                            np.abs(
                                t-cutoff
                            )
                        )
                    ]
                ),

            "n_origins_at_cutoff":
                int(
                    origins[
                        np.argmin(
                            np.abs(
                                t-cutoff
                            )
                        )
                    ]
                ),
        }
    )

conv = pd.DataFrame(
    conv_rows
)

conv.to_csv(
    OUT
    / "F46E_integral_convergence.csv",
    index=False,
)


# ================================================================
# Collective correlation
#
# Do not fit the long negative tail.
# Quantify first-zero integral and statistical support only.
# ================================================================

cp_integral, cp_zero_window, cp_stop = (
    first_zero_integral(
        t,
        CP,
    )
)

cp_zero_origins = int(
    origins[
        cp_stop - 1
    ]
)


# ================================================================
# Block reproducibility
# ================================================================

blocks = pd.read_csv(
    IN
    / "F46D_block_correlations.csv"
)

block_results = []

for block_id, g in blocks.groupby(
    "block",
    sort=True,
):

    g = g.sort_values(
        "lag_ps"
    )

    tb = g["lag_ps"].to_numpy(
        dtype=float
    )

    c1b = g["C1"].to_numpy(
        dtype=float
    )

    c2b = g["C2"].to_numpy(
        dtype=float
    )

    cpb = g[
        "C_collective"
    ].to_numpy(
        dtype=float
    )

    ob = g[
        "n_time_origins"
    ].to_numpy(
        dtype=int
    )

    # Use max 10 ps for ~25 ps blocks.
    # This retains a reasonable number of origins.
    fit1 = fit_models(
        tb,
        c1b,
        fit_max=10.0,
    )

    fit2 = fit_models(
        tb,
        c2b,
        fit_max=10.0,
    )

    cp_int, cp_win, cp_n = (
        first_zero_integral(
            tb,
            cpb,
        )
    )

    block_results.append(
        {
            "block":
                int(block_id),

            "C1_integral_10ps":
                integrate_to(
                    tb,
                    c1b,
                    10.0,
                ),

            "C2_integral_10ps":
                integrate_to(
                    tb,
                    c2b,
                    10.0,
                ),

            "C1_KWW_tau_ps":
                fit1[
                    "kww_tau_ps"
                ],

            "C1_KWW_beta":
                fit1[
                    "kww_beta"
                ],

            "C1_KWW_integrated_tau_ps":
                fit1[
                    "kww_integrated_tau_ps"
                ],

            "C2_KWW_tau_ps":
                fit2[
                    "kww_tau_ps"
                ],

            "C2_KWW_beta":
                fit2[
                    "kww_beta"
                ],

            "C2_KWW_integrated_tau_ps":
                fit2[
                    "kww_integrated_tau_ps"
                ],

            "collective_first_zero_integral_ps":
                cp_int,

            "collective_first_zero_window_ps":
                cp_win,

            "collective_origins_at_zero_window":
                int(
                    ob[
                        min(
                            cp_n-1,
                            len(ob)-1,
                        )
                    ]
                ),
        }
    )


block_df = pd.DataFrame(
    block_results
)

block_df.to_csv(
    OUT
    / "F46E_block_relaxation_summary.csv",
    index=False,
)


# ================================================================
# Block statistics
# ================================================================

def stat(series):
    x = np.asarray(
        series,
        dtype=float
    )

    return {
        "mean":
            float(
                np.mean(x)
            ),

        "std":
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


block_stats = {
    "C1_KWW_integrated_tau_ps":
        stat(
            block_df[
                "C1_KWW_integrated_tau_ps"
            ]
        ),

    "C2_KWW_integrated_tau_ps":
        stat(
            block_df[
                "C2_KWW_integrated_tau_ps"
            ]
        ),

    "C1_integral_10ps":
        stat(
            block_df[
                "C1_integral_10ps"
            ]
        ),

    "C2_integral_10ps":
        stat(
            block_df[
                "C2_integral_10ps"
            ]
        ),

    "collective_first_zero_integral_ps":
        stat(
            block_df[
                "collective_first_zero_integral_ps"
            ]
        ),
}


# ================================================================
# Summary
# ================================================================

summary = {
    "gate": "F46E",

    "temperature_K":
        300.0,

    "n_frames":
        int(len(df)),

    "frame_spacing_ps":
        float(
            t[1] - t[0]
        ),

    "fit_window_ps":
        20.0,

    "C1_full_fit":
        fit_C1,

    "C2_full_fit":
        fit_C2,

    "C1_numeric_integral_20ps":
        integrate_to(
            t,
            C1,
            20.0,
        ),

    "C1_numeric_integral_40ps":
        integrate_to(
            t,
            C1,
            40.0,
        ),

    "C2_numeric_integral_20ps":
        integrate_to(
            t,
            C2,
            20.0,
        ),

    "C2_numeric_integral_40ps":
        integrate_to(
            t,
            C2,
            40.0,
        ),

    "collective_first_zero_integral_ps":
        cp_integral,

    "collective_first_zero_window_ps":
        cp_zero_window,

    "collective_n_time_origins_at_zero_window":
        cp_zero_origins,

    "block_statistics":
        block_stats,

    "interpretation_boundary":
        (
            "These are effective orientational "
            "relaxation times of the atomistic "
            "water environment at 300 K. They are "
            "not excitonic T2/T2*, microscopic "
            "dephasing rates, or epsilon(omega)."
        ),
}


with (
    OUT
    / "F46E_summary.json"
).open("w") as fh:

    json.dump(
        summary,
        fh,
        indent=2,
    )


# ================================================================
# Human-readable report
# ================================================================

print(
    "============================================================"
)

print(
    "F46E RELAXATION QUANTIFICATION COMPLETE"
)

print(
    "============================================================"
)

print()
print("FULL C1 FIT")
for k, v in fit_C1.items():
    print(f"{k}={v}")

print()
print("FULL C2 FIT")
for k, v in fit_C2.items():
    print(f"{k}={v}")

print()
print("NUMERICAL INTEGRAL CONVERGENCE")
print(
    conv.to_string(
        index=False
    )
)

print()
print("BLOCK RESULTS")
print(
    block_df.to_string(
        index=False
    )
)

print()
print("BLOCK STATISTICS")
print(
    json.dumps(
        block_stats,
        indent=2,
    )
)

print()
print("COLLECTIVE CORRELATION")
print(
    f"first_zero_integral_ps="
    f"{cp_integral}"
)

print(
    f"first_zero_window_ps="
    f"{cp_zero_window}"
)

print(
    f"n_time_origins_at_zero_window="
    f"{cp_zero_origins}"
)

print()
print(
    "SUMMARY_FILE="
    + str(
        OUT
        / "F46E_summary.json"
    )
)

print(
    "ACTUAL_M3_ANALYSIS=YES"
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
