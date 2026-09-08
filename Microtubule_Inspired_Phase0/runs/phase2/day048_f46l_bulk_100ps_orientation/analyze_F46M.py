from pathlib import Path
from itertools import combinations
import json
import math
import sys

import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import minimize_vectors
from scipy.optimize import curve_fit
from scipy.special import gamma


OUT = Path(sys.argv[1])

ENG_TPR = Path(sys.argv[2])
ENG_XTC = Path(sys.argv[3])

BULK_TPR = Path(sys.argv[4])
BULK_XTC = Path(sys.argv[5])
BULK_GRO = Path(sys.argv[6])


# ============================================================
# BASIC HELPERS
# ============================================================

def mono(t, tau):
    return np.exp(-t / tau)


def kww(t, tau, beta):
    return np.exp(-(t / tau) ** beta)


def read_box_from_gro(path):

    lines = path.read_text().splitlines()
    parts = lines[-1].split()

    if len(parts) < 3:
        raise RuntimeError(
            "Could not read orthorhombic box from GRO"
        )

    lx, ly, lz = map(float, parts[:3])

    return lx, ly, lz


def density_from_water_count(nwater, gro):

    lx, ly, lz = read_box_from_gro(gro)

    volume_nm3 = lx * ly * lz

    mw_g_mol = 18.01528
    NA = 6.02214076e23

    mass_g = (
        nwater
        * mw_g_mol
        / NA
    )

    volume_cm3 = (
        volume_nm3
        * 1e-21
    )

    density_g_cm3 = (
        mass_g
        / volume_cm3
    )

    return {
        "box_nm":
            [lx, ly, lz],
        "volume_nm3":
            volume_nm3,
        "density_g_cm3":
            density_g_cm3,
        "density_kg_m3":
            1000.0 * density_g_cm3,
    }


# ============================================================
# WATER DIPOLE VECTOR EXTRACTION
# ============================================================

def load_vectors(tpr, xtc, label):

    u = mda.Universe(
        str(tpr),
        str(xtc),
    )

    waters = u.residues[
        u.residues.resnames == "SOL"
    ]

    nwater = len(waters)

    print(
        f"{label}_N_WATER={nwater}"
    )

    ow = []
    h1 = []
    h2 = []

    for res in waters:

        amap = {
            a.name: a.index
            for a in res.atoms
        }

        needed = {
            "OW",
            "HW1",
            "HW2",
        }

        if not needed.issubset(
            amap.keys()
        ):
            raise RuntimeError(
                f"{label}: missing expected "
                f"TIP4P/2005 atoms in residue "
                f"{res.resid}"
            )

        ow.append(amap["OW"])
        h1.append(amap["HW1"])
        h2.append(amap["HW2"])

    ow = np.asarray(ow)
    h1 = np.asarray(h1)
    h2 = np.asarray(h2)

    vectors = []
    times = []

    for ts in u.trajectory:

        box = ts.dimensions

        if box is None:
            raise RuntimeError(
                f"{label}: missing box"
            )

        pos = u.atoms.positions

        rO = pos[ow]

        OH1 = minimize_vectors(
            pos[h1] - rO,
            box,
        )

        OH2 = minimize_vectors(
            pos[h2] - rO,
            box,
        )

        # Molecular dipole direction:
        # O -> midpoint(H1,H2)
        mu = 0.5 * (
            OH1 + OH2
        )

        norm = np.linalg.norm(
            mu,
            axis=1,
        )

        if np.any(norm <= 0):
            raise RuntimeError(
                f"{label}: zero molecular vector"
            )

        mu = (
            mu
            / norm[:, None]
        )

        vectors.append(
            mu.astype(np.float32)
        )

        times.append(
            float(ts.time)
        )

    times = np.asarray(times)
    vectors = np.asarray(vectors)

    if len(times) != 201:
        raise RuntimeError(
            f"{label}: expected 201 frames; "
            f"got {len(times)}"
        )

    if not np.allclose(
        np.diff(times),
        0.5,
        atol=1e-5,
        rtol=0,
    ):
        raise RuntimeError(
            f"{label}: sampling not exactly 0.5 ps"
        )

    print(
        f"{label}_N_FRAMES={len(times)}"
    )

    print(
        f"{label}_DURATION_PS="
        f"{times[-1]-times[0]}"
    )

    return (
        times-times[0],
        vectors,
    )


# ============================================================
# CORRELATIONS
# ============================================================

def correlations(vectors):

    nf = vectors.shape[0]

    c1 = np.zeros(nf)
    c2 = np.zeros(nf)
    origins = np.zeros(
        nf,
        dtype=int,
    )

    for lag in range(nf):

        a = vectors[:nf-lag]
        b = vectors[lag:]

        dots = np.sum(
            a*b,
            axis=2,
            dtype=np.float64,
        )

        c1[lag] = np.mean(dots)

        c2[lag] = np.mean(
            0.5
            * (
                3.0*dots*dots
                - 1.0
            )
        )

        origins[lag] = (
            nf-lag
        )

    return (
        c1,
        c2,
        origins,
    )


# ============================================================
# FITS
# ============================================================

def fit_relaxation(
    time,
    corr,
    fit_limit=20.0,
):

    mask = (
        (time <= fit_limit)
        & np.isfinite(corr)
        & (corr > 0)
    )

    t = time[mask]
    y = corr[mask]

    if len(t) < 8:
        raise RuntimeError(
            "Insufficient positive correlation "
            "points for fit"
        )

    pm, _ = curve_fit(
        mono,
        t,
        y,
        p0=(3.0,),
        bounds=(
            (1e-6,),
            (1000.0,),
        ),
        maxfev=50000,
    )

    tau_mono = float(pm[0])

    mono_pred = mono(
        t,
        tau_mono,
    )

    rss_mono = float(
        np.sum(
            (y-mono_pred)**2
        )
    )

    pk, _ = curve_fit(
        kww,
        t,
        y,
        p0=(3.0, 0.8),
        bounds=(
            (1e-6, 0.05),
            (1000.0, 2.0),
        ),
        maxfev=50000,
    )

    tau = float(pk[0])
    beta = float(pk[1])

    kww_pred = kww(
        t,
        tau,
        beta,
    )

    rss_kww = float(
        np.sum(
            (y-kww_pred)**2
        )
    )

    n = len(t)

    def AIC(rss, k):

        rss = max(
            rss,
            np.finfo(float).tiny,
        )

        return (
            n*np.log(rss/n)
            + 2*k
        )

    aic_mono = float(
        AIC(
            rss_mono,
            1,
        )
    )

    aic_kww = float(
        AIC(
            rss_kww,
            2,
        )
    )

    integrated = float(
        (
            tau
            / beta
        )
        * gamma(
            1.0/beta
        )
    )

    return {
        "mono_tau_ps":
            tau_mono,
        "mono_rss":
            rss_mono,
        "mono_aic":
            aic_mono,
        "kww_tau_ps":
            tau,
        "kww_beta":
            beta,
        "kww_integrated_tau_ps":
            integrated,
        "kww_rss":
            rss_kww,
        "kww_aic":
            aic_kww,
        "delta_aic_mono_minus_kww":
            (
                aic_mono
                - aic_kww
            ),
    }


def integrate_to(
    time,
    corr,
    cutoff,
):

    mask = (
        time <= cutoff
    )

    return float(
        np.trapezoid(
            corr[mask],
            time[mask],
        )
    )


# ============================================================
# SYSTEM ANALYSIS
# ============================================================

def analyze(
    label,
    tpr,
    xtc,
):

    time, vectors = (
        load_vectors(
            tpr,
            xtc,
            label,
        )
    )

    c1, c2, origins = (
        correlations(
            vectors
        )
    )

    f1 = fit_relaxation(
        time,
        c1,
        20.0,
    )

    f2 = fit_relaxation(
        time,
        c2,
        20.0,
    )

    integrals = {}

    for cutoff in [
        2,
        5,
        10,
        20,
        30,
        40,
        50,
    ]:

        integrals[
            f"C1_{cutoff}"
        ] = integrate_to(
            time,
            c1,
            cutoff,
        )

        integrals[
            f"C2_{cutoff}"
        ] = integrate_to(
            time,
            c2,
            cutoff,
        )

    pd.DataFrame({
        "time_ps":
            time,
        "C1":
            c1,
        "C2":
            c2,
        "n_origins":
            origins,
    }).to_csv(
        OUT
        / f"{label}_C1_C2.csv",
        index=False,
    )

    # --------------------------------------------------------
    # 25 ps blocks
    # --------------------------------------------------------

    block_rows = []

    windows = [
        (0, 25),
        (25, 50),
        (50, 75),
        (75, 100),
    ]

    for ib, (
        start,
        stop,
    ) in enumerate(
        windows,
        1,
    ):

        mask = (
            (time >= start)
            & (time <= stop)
        )

        bt = (
            time[mask]
            - time[mask][0]
        )

        bv = vectors[mask]

        bc1, bc2, _ = (
            correlations(
                bv
            )
        )

        bf1 = fit_relaxation(
            bt,
            bc1,
            10.0,
        )

        bf2 = fit_relaxation(
            bt,
            bc2,
            10.0,
        )

        block_rows.append({
            "system":
                label,
            "block":
                ib,
            "start_ps":
                start,
            "end_ps":
                stop,

            "C1_tau_integrated_ps":
                bf1[
                    "kww_integrated_tau_ps"
                ],

            "C1_tau_kww_ps":
                bf1[
                    "kww_tau_ps"
                ],

            "C1_beta":
                bf1[
                    "kww_beta"
                ],

            "C1_integral_10ps":
                integrate_to(
                    bt,
                    bc1,
                    10.0,
                ),

            "C2_tau_integrated_ps":
                bf2[
                    "kww_integrated_tau_ps"
                ],

            "C2_tau_kww_ps":
                bf2[
                    "kww_tau_ps"
                ],

            "C2_beta":
                bf2[
                    "kww_beta"
                ],

            "C2_integral_10ps":
                integrate_to(
                    bt,
                    bc2,
                    10.0,
                ),
        })

    result = {
        "C1":
            f1,
        "C2":
            f2,
        "integrals_ps":
            integrals,
        "n_water":
            int(
                vectors.shape[1]
            ),
    }

    return (
        result,
        pd.DataFrame(
            block_rows
        ),
    )


# ============================================================
# EFFECT-SIZE HELPERS
# ============================================================

def hedges_g(a, b):

    a = np.asarray(a)
    b = np.asarray(b)

    n1 = len(a)
    n2 = len(b)

    s1 = np.var(
        a,
        ddof=1,
    )

    s2 = np.var(
        b,
        ddof=1,
    )

    pooled = math.sqrt(
        (
            (n1-1)*s1
            + (n2-1)*s2
        )
        / (
            n1+n2-2
        )
    )

    if pooled == 0:
        return float("nan")

    d = (
        np.mean(a)
        - np.mean(b)
    ) / pooled

    df = (
        n1+n2-2
    )

    J = (
        1
        - 3
        / (
            4*df
            - 1
        )
    )

    return float(
        J*d
    )


def permutation_p(a, b):

    a = np.asarray(a)
    b = np.asarray(b)

    allv = np.concatenate(
        [a, b]
    )

    obs = abs(
        np.mean(a)
        - np.mean(b)
    )

    na = len(a)

    count = 0
    total = 0

    for comb in combinations(
        range(len(allv)),
        na,
    ):

        mask = np.zeros(
            len(allv),
            dtype=bool,
        )

        mask[list(comb)] = True

        x = allv[mask]
        y = allv[~mask]

        diff = abs(
            np.mean(x)
            - np.mean(y)
        )

        total += 1

        if diff >= (
            obs - 1e-15
        ):
            count += 1

    return float(
        count/total
    )


# ============================================================
# RUN
# ============================================================

print()
print("ANALYZING_ENGINEERED=YES")

eng, eng_blocks = analyze(
    "engineered",
    ENG_TPR,
    ENG_XTC,
)

print()
print("ANALYZING_BULK=YES")

bulk, bulk_blocks = analyze(
    "bulk",
    BULK_TPR,
    BULK_XTC,
)


density = density_from_water_count(
    bulk["n_water"],
    BULK_GRO,
)


blocks = pd.concat(
    [
        eng_blocks,
        bulk_blocks,
    ],
    ignore_index=True,
)

blocks.to_csv(
    OUT
    / "F46M_matched_25ps_blocks.csv",
    index=False,
)


metrics = [
    "C1_tau_integrated_ps",
    "C1_integral_10ps",
    "C2_tau_integrated_ps",
    "C2_integral_10ps",
]


block_comparison = {}

for metric in metrics:

    e = eng_blocks[
        metric
    ].to_numpy()

    b = bulk_blocks[
        metric
    ].to_numpy()

    em = float(
        np.mean(e)
    )

    bm = float(
        np.mean(b)
    )

    block_comparison[
        metric
    ] = {
        "engineered_mean":
            em,

        "engineered_sd":
            float(
                np.std(
                    e,
                    ddof=1,
                )
            ),

        "bulk_mean":
            bm,

        "bulk_sd":
            float(
                np.std(
                    b,
                    ddof=1,
                )
            ),

        "difference":
            em-bm,

        "ratio":
            em/bm,

        "percent_change":
            100.0
            * (
                em/bm
                - 1.0
            ),

        "hedges_g":
            hedges_g(
                e,
                b,
            ),

        "exact_permutation_p":
            permutation_p(
                e,
                b,
            ),
    }


C1_ratio = (
    eng["C1"][
        "kww_integrated_tau_ps"
    ]
    / bulk["C1"][
        "kww_integrated_tau_ps"
    ]
)

C2_ratio = (
    eng["C2"][
        "kww_integrated_tau_ps"
    ]
    / bulk["C2"][
        "kww_integrated_tau_ps"
    ]
)


summary = {
    "bulk_density":
        density,

    "engineered":
        eng,

    "bulk":
        bulk,

    "full_comparison": {
        "C1_ratio_engineered_over_bulk":
            C1_ratio,

        "C1_percent_change":
            100.0*(C1_ratio-1.0),

        "C2_ratio_engineered_over_bulk":
            C2_ratio,

        "C2_percent_change":
            100.0*(C2_ratio-1.0),
    },

    "block_comparison":
        block_comparison,

    "interpretation_boundary": [
        (
            "Comparison is confined HBN-PYR "
            "environment versus equilibrated "
            "bulk TIP4P/2005 water."
        ),

        (
            "The difference represents the combined "
            "effect of confinement, interfaces, "
            "electrostatics, and pyrene chromophores."
        ),

        (
            "C1 and C2 relaxation times are "
            "environmental molecular orientational "
            "times and are not quantum T2 or T2*."
        ),

        (
            "Four contiguous 25 ps blocks provide "
            "a robustness diagnostic but are not "
            "independent replicas."
        ),
    ],
}


with (
    OUT
    / "F46M_summary.json"
).open("w") as fh:

    json.dump(
        summary,
        fh,
        indent=2,
    )


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("="*72)
print("F46M BULK NVT DENSITY")
print("="*72)

print(
    "BULK_BOX_NM="
    + ",".join(
        f"{x:.8f}"
        for x in density[
            "box_nm"
        ]
    )
)

print(
    f"BULK_VOLUME_NM3="
    f"{density['volume_nm3']:.9f}"
)

print(
    f"BULK_DENSITY_KG_M3="
    f"{density['density_kg_m3']:.6f}"
)


print()
print("="*72)
print("F46M FULL-TRAJECTORY RESULTS")
print("="*72)

for label, r in [
    ("ENGINEERED", eng),
    ("BULK", bulk),
]:

    print()
    print(label)

    print(
        "C1_KWW_TAU_PS="
        f"{r['C1']['kww_tau_ps']:.12f}"
    )

    print(
        "C1_BETA="
        f"{r['C1']['kww_beta']:.12f}"
    )

    print(
        "C1_INTEGRATED_KWW_PS="
        f"{r['C1']['kww_integrated_tau_ps']:.12f}"
    )

    print(
        "C1_DELTA_AIC_MONO_MINUS_KWW="
        f"{r['C1']['delta_aic_mono_minus_kww']:.12f}"
    )

    print(
        "C1_NUMERIC_INT_10PS="
        f"{r['integrals_ps']['C1_10']:.12f}"
    )

    print(
        "C1_NUMERIC_INT_20PS="
        f"{r['integrals_ps']['C1_20']:.12f}"
    )

    print(
        "C1_NUMERIC_INT_40PS="
        f"{r['integrals_ps']['C1_40']:.12f}"
    )

    print(
        "C2_KWW_TAU_PS="
        f"{r['C2']['kww_tau_ps']:.12f}"
    )

    print(
        "C2_BETA="
        f"{r['C2']['kww_beta']:.12f}"
    )

    print(
        "C2_INTEGRATED_KWW_PS="
        f"{r['C2']['kww_integrated_tau_ps']:.12f}"
    )

    print(
        "C2_DELTA_AIC_MONO_MINUS_KWW="
        f"{r['C2']['delta_aic_mono_minus_kww']:.12f}"
    )

    print(
        "C2_NUMERIC_INT_10PS="
        f"{r['integrals_ps']['C2_10']:.12f}"
    )

    print(
        "C2_NUMERIC_INT_20PS="
        f"{r['integrals_ps']['C2_20']:.12f}"
    )

    print(
        "C2_NUMERIC_INT_40PS="
        f"{r['integrals_ps']['C2_40']:.12f}"
    )


print()
print("="*72)
print("ENGINEERED vs BULK — FULL TRAJECTORY")
print("="*72)

print(
    f"C1_RATIO_ENGINEERED_OVER_BULK="
    f"{C1_ratio:.12f}"
)

print(
    f"C1_PERCENT_CHANGE="
    f"{100*(C1_ratio-1):.6f}"
)

print(
    f"C2_RATIO_ENGINEERED_OVER_BULK="
    f"{C2_ratio:.12f}"
)

print(
    f"C2_PERCENT_CHANGE="
    f"{100*(C2_ratio-1):.6f}"
)


print()
print("="*72)
print("MATCHED 25-PS BLOCK RESULTS")
print("="*72)

print(
    blocks.to_string(
        index=False
    )
)


print()
print("="*72)
print("BLOCK ROBUSTNESS COMPARISON")
print("="*72)

for metric, r in (
    block_comparison.items()
):

    print()
    print(
        f"METRIC={metric}"
    )

    print(
        f"ENGINEERED_MEAN="
        f"{r['engineered_mean']}"
    )

    print(
        f"ENGINEERED_SD="
        f"{r['engineered_sd']}"
    )

    print(
        f"BULK_MEAN="
        f"{r['bulk_mean']}"
    )

    print(
        f"BULK_SD="
        f"{r['bulk_sd']}"
    )

    print(
        f"PERCENT_CHANGE="
        f"{r['percent_change']}"
    )

    print(
        f"HEDGES_G="
        f"{r['hedges_g']}"
    )

    print(
        f"EXACT_PERMUTATION_P="
        f"{r['exact_permutation_p']}"
    )


print()
print(
    "SUMMARY_FILE="
    + str(
        OUT
        / "F46M_summary.json"
    )
)

print(
    "F46M_ANALYSIS_COMPLETE=YES"
)
