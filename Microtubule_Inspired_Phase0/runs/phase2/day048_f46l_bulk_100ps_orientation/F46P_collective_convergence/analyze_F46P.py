from pathlib import Path
import sys
import numpy as np
import pandas as pd

OUT = Path(sys.argv[1])
SRC = Path(sys.argv[2])

EPS0 = 8.8541878128e-12
KB = 1.380649e-23
T = 300.0
DEBYE_TO_CM = 3.33564e-30

FILES = {
    "engineered": SRC / "engineered_total_water_dipole_D.csv",
    "bulk": SRC / "bulk_total_water_dipole_D.csv",
}

VOLUMES_NM3 = {
    "engineered": 547.747619629,
    "bulk": 501.540008545,
}


def fluctuation_tensor(M_D, volume_nm3):

    M = M_D * DEBYE_TO_CM

    dM = (
        M
        - np.mean(M, axis=0, keepdims=True)
    )

    cov = (
        dM.T @ dM
        / len(dM)
    )

    V = volume_nm3 * 1e-27

    tensor = (
        cov
        / (
            EPS0
            * V
            * KB
            * T
        )
    )

    return tensor


def isotropy_metrics(tensor):

    diag = np.diag(tensor)

    mean_diag = np.mean(diag)
    sd_diag = np.std(diag, ddof=1)

    cv = (
        sd_diag / mean_diag
        if mean_diag != 0
        else np.nan
    )

    max_min = (
        np.max(diag) / np.min(diag)
        if np.min(diag) > 0
        else np.nan
    )

    return {
        "diag_mean": float(mean_diag),
        "diag_sd": float(sd_diag),
        "diag_cv": float(cv),
        "diag_max_min_ratio": float(max_min),
    }


rows = []

block_sets = {
    "25ps_block1": (0,25),
    "25ps_block2": (25,50),
    "25ps_block3": (50,75),
    "25ps_block4": (75,100),

    "50ps_block1": (0,50),
    "50ps_block2": (50,100),

    "cum_0_25": (0,25),
    "cum_0_50": (0,50),
    "cum_0_75": (0,75),
    "cum_0_100": (0,100),
}


for system, path in FILES.items():

    df = pd.read_csv(path)

    time = df["time_ps"].to_numpy()

    M = df[
        ["Mx_D","My_D","Mz_D"]
    ].to_numpy()

    for name, (start,end) in block_sets.items():

        mask = (
            (time >= start)
            & (time <= end)
        )

        X = M[mask]

        tensor = fluctuation_tensor(
            X,
            VOLUMES_NM3[system],
        )

        iso = isotropy_metrics(tensor)

        mean_M = np.mean(X,axis=0)

        # Linear drift in each dipole component.
        tx = time[mask]

        slopes = []

        for j in range(3):
            slopes.append(
                float(
                    np.polyfit(
                        tx,
                        X[:,j],
                        1,
                    )[0]
                )
            )

        rows.append({
            "system":system,
            "window":name,
            "start_ps":start,
            "end_ps":end,
            "n_frames":int(np.sum(mask)),

            "XX":float(tensor[0,0]),
            "YY":float(tensor[1,1]),
            "ZZ":float(tensor[2,2]),

            "XY":float(tensor[0,1]),
            "XZ":float(tensor[0,2]),
            "YZ":float(tensor[1,2]),

            "trace":float(np.trace(tensor)),

            "diag_cv":
                iso["diag_cv"],

            "diag_max_min_ratio":
                iso["diag_max_min_ratio"],

            "mean_Mx_D":float(mean_M[0]),
            "mean_My_D":float(mean_M[1]),
            "mean_Mz_D":float(mean_M[2]),

            "drift_Mx_D_per_ps":
                slopes[0],

            "drift_My_D_per_ps":
                slopes[1],

            "drift_Mz_D_per_ps":
                slopes[2],
        })


res = pd.DataFrame(rows)

res.to_csv(
    OUT/"F46P_tensor_windows.csv",
    index=False
)


print()
print("="*76)
print("F46P — 25 ps INDEPENDENT BLOCKS")
print("="*76)

sel = res[
    res.window.str.startswith(
        "25ps_block"
    )
]

print(
    sel[
        [
            "system",
            "window",
            "XX","YY","ZZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
            "drift_Mx_D_per_ps",
            "drift_My_D_per_ps",
            "drift_Mz_D_per_ps",
        ]
    ].to_string(index=False)
)


print()
print("="*76)
print("F46P — 50 ps INDEPENDENT BLOCKS")
print("="*76)

sel50 = res[
    res.window.str.startswith(
        "50ps_block"
    )
]

print(
    sel50[
        [
            "system",
            "window",
            "XX","YY","ZZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    ].to_string(index=False)
)


print()
print("="*76)
print("F46P — CUMULATIVE CONVERGENCE")
print("="*76)

cum = res[
    res.window.str.startswith(
        "cum_"
    )
]

print(
    cum[
        [
            "system",
            "window",
            "XX","YY","ZZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    ].to_string(index=False)
)


print()
print("="*76)
print("F46P — ENGINEERED/BULK RATIOS AT MATCHED WINDOWS")
print("="*76)

ratio_rows=[]

for window in [
    "25ps_block1",
    "25ps_block2",
    "25ps_block3",
    "25ps_block4",
    "50ps_block1",
    "50ps_block2",
    "cum_0_25",
    "cum_0_50",
    "cum_0_75",
    "cum_0_100",
]:

    e = res[
        (res.system=="engineered")
        & (res.window==window)
    ].iloc[0]

    b = res[
        (res.system=="bulk")
        & (res.window==window)
    ].iloc[0]

    row = {
        "window":window
    }

    for metric in [
        "XX","YY","ZZ","trace"
    ]:

        row[
            metric+"_ratio"
        ] = (
            e[metric]
            / b[metric]
        )

    ratio_rows.append(row)


rat = pd.DataFrame(
    ratio_rows
)

rat.to_csv(
    OUT/"F46P_engineered_bulk_ratios.csv",
    index=False
)

print(
    rat.to_string(index=False)
)


print()
print("="*76)
print("F46P — BULK ISOTROPY CONVERGENCE")
print("="*76)

bulk_cum = cum[
    cum.system=="bulk"
].sort_values("end_ps")

for _,r in bulk_cum.iterrows():

    print(
        f"{r['window']}: "
        f"CV={r['diag_cv']:.6f} "
        f"MAX_MIN={r['diag_max_min_ratio']:.6f} "
        f"TRACE={r['trace']:.6f}"
    )


full_bulk = res[
    (res.system=="bulk")
    & (res.window=="cum_0_100")
].iloc[0]

full_eng = res[
    (res.system=="engineered")
    & (res.window=="cum_0_100")
].iloc[0]


# Conservative diagnostic only.
bulk_iso_pass = (
    full_bulk[
        "diag_max_min_ratio"
    ] <= 1.5
)

bulk_cv_pass = (
    full_bulk[
        "diag_cv"
    ] <= 0.20
)

print()
print(
    "BULK_ISOTROPY_MAXMIN_GATE="
    + (
        "PASS"
        if bulk_iso_pass
        else "FAIL"
    )
)

print(
    "BULK_ISOTROPY_CV_GATE="
    + (
        "PASS"
        if bulk_cv_pass
        else "FAIL"
    )
)

print(
    f"FULL_TRACE_RATIO_ENGINEERED_OVER_BULK="
    f"{full_eng['trace']/full_bulk['trace']:.12f}"
)


if (
    bulk_iso_pass
    and bulk_cv_pass
):
    decision = (
        "COLLECTIVE_FLUCTUATIONS_PLAUSIBLY_CONVERGED"
    )
else:
    decision = (
        "LONGER_COLLECTIVE_SAMPLING_REQUIRED"
    )

print(
    f"F46P_DECISION={decision}"
)

print(
    "F46P_ANALYSIS_COMPLETE=YES"
)
