from pathlib import Path
import sys
import json

import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import minimize_vectors


OUT = Path(sys.argv[1])

SYSTEMS = {
    "engineered": (
        Path(sys.argv[2]),
        Path(sys.argv[3]),
    ),
    "bulk": (
        Path(sys.argv[4]),
        Path(sys.argv[5]),
    ),
}


# ============================================================
# CONSTANTS
# ============================================================

QH = +0.5564
QM = -1.1128

E_CHARGE = 1.602176634e-19
NM_TO_M = 1e-9
DEBYE_C_M = 3.33564e-30

EPS0 = 8.8541878128e-12
KB = 1.380649e-23
TEMP_K = 300.0


# ============================================================
# LOAD TOTAL WATER DIPOLE
# ============================================================

def load_total_dipole(tpr, xtc, label):

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
    mw = []

    for res in waters:

        amap = {
            a.name: a.index
            for a in res.atoms
        }

        required = {
            "OW",
            "HW1",
            "HW2",
            "MW",
        }

        if not required.issubset(
            amap.keys()
        ):
            raise RuntimeError(
                f"{label}: unexpected water atom names "
                f"in residue {res.resid}"
            )

        ow.append(amap["OW"])
        h1.append(amap["HW1"])
        h2.append(amap["HW2"])
        mw.append(amap["MW"])

    ow = np.asarray(ow)
    h1 = np.asarray(h1)
    h2 = np.asarray(h2)
    mw = np.asarray(mw)

    times = []
    M_e_nm = []
    boxes = []

    for ts in u.trajectory:

        box = ts.dimensions

        if box is None:
            raise RuntimeError(
                f"{label}: no box at frame {ts.frame}"
            )

        p = u.atoms.positions * 0.1
        # MDAnalysis coordinates are Å; convert to nm.

        rO = p[ow]

        # minimize_vectors expects a box in same length units.
        box_nm = box.copy()
        box_nm[:3] *= 0.1

        OH1 = minimize_vectors(
            p[h1] - rO,
            box_nm,
        )

        OH2 = minimize_vectors(
            p[h2] - rO,
            box_nm,
        )

        OM = minimize_vectors(
            p[mw] - rO,
            box_nm,
        )

        mu_mol = (
            QH * OH1
            + QH * OH2
            + QM * OM
        )

        M = np.sum(
            mu_mol,
            axis=0,
            dtype=np.float64,
        )

        times.append(
            float(ts.time)
        )

        M_e_nm.append(M)

        boxes.append(
            box_nm[:3].copy()
        )

    times = np.asarray(times)
    times -= times[0]

    M_e_nm = np.asarray(M_e_nm)
    boxes = np.asarray(boxes)

    if len(times) != 201:
        raise RuntimeError(
            f"{label}: expected 201 frames, "
            f"got {len(times)}"
        )

    if not np.allclose(
        np.diff(times),
        0.5,
        atol=1e-5,
        rtol=0,
    ):
        raise RuntimeError(
            f"{label}: invalid sampling"
        )

    return {
        "time_ps": times,
        "M_e_nm": M_e_nm,
        "M_D": (
            M_e_nm
            * E_CHARGE
            * NM_TO_M
            / DEBYE_C_M
        ),
        "boxes_nm": boxes,
        "n_water": len(waters),
    }


# ============================================================
# CONNECTED VECTOR AUTOCORRELATION
# ============================================================

def connected_vector_corr(M):

    X = M - np.mean(
        M,
        axis=0,
        keepdims=True,
    )

    nf = len(X)

    C = np.empty(nf)
    origins = np.empty(
        nf,
        dtype=int,
    )

    denom = np.mean(
        np.sum(
            X*X,
            axis=1,
        )
    )

    for lag in range(nf):

        dot = np.sum(
            X[:nf-lag]
            * X[lag:],
            axis=1,
        )

        C[lag] = (
            np.mean(dot)
            / denom
        )

        origins[lag] = nf-lag

    return (
        C,
        origins,
    )


def integral_until(
    t,
    c,
    cutoff,
):

    mask = t <= cutoff

    return float(
        np.trapezoid(
            c[mask],
            t[mask],
        )
    )


def first_zero_integral(t,c):

    idx = np.where(
        c <= 0
    )[0]

    if len(idx) == 0:
        stop = len(c)
        zero_t = None
    else:
        stop = int(idx[0]) + 1
        zero_t = float(
            t[idx[0]]
        )

    value = float(
        np.trapezoid(
            c[:stop],
            t[:stop],
        )
    )

    return (
        value,
        zero_t,
    )


# ============================================================
# FLUCTUATION TENSOR
# ============================================================

def fluctuation_tensor(M_e_nm, boxes_nm):

    # Convert e*nm -> C*m.
    M = (
        M_e_nm
        * E_CHARGE
        * NM_TO_M
    )

    dM = (
        M
        - np.mean(
            M,
            axis=0,
            keepdims=True,
        )
    )

    cov = (
        dM.T @ dM
        / len(dM)
    )

    volume_nm3 = float(
        np.mean(
            np.prod(
                boxes_nm,
                axis=1,
            )
        )
    )

    V_m3 = (
        volume_nm3
        * 1e-27
    )

    # Dimensionless polarization fluctuation tensor:
    # beta / (eps0 V) <dMi dMj>
    #
    # For homogeneous isotropic bulk under the appropriate
    # electrostatic ensemble this is related to dielectric
    # susceptibility. For the confined heterogeneous system
    # we explicitly do NOT call this epsilon.
    tensor = (
        cov
        / (
            EPS0
            * V_m3
            * KB
            * TEMP_K
        )
    )

    return (
        cov,
        tensor,
        volume_nm3,
    )


# ============================================================
# WINDOW ANALYSIS
# ============================================================

def window_metrics(
    time,
    M,
    start,
    stop,
):

    mask = (
        (time >= start)
        & (time <= stop)
    )

    t = (
        time[mask]
        - time[mask][0]
    )

    X = M[mask]

    c, origins = (
        connected_vector_corr(X)
    )

    iz, zt = first_zero_integral(
        t,
        c,
    )

    return {
        "start_ps": start,
        "end_ps": stop,
        "duration_ps": stop-start,

        "C_final":
            float(c[-1]),

        "C_2ps":
            float(
                c[
                    np.argmin(
                        np.abs(t-2)
                    )
                ]
            ),

        "C_5ps":
            float(
                c[
                    np.argmin(
                        np.abs(t-5)
                    )
                ]
            ),

        "C_10ps":
            float(
                c[
                    np.argmin(
                        np.abs(t-10)
                    )
                ]
            ),

        "integral_5ps":
            integral_until(
                t,c,5
            ),

        "integral_10ps":
            integral_until(
                t,c,10
            ),

        "first_zero_integral_ps":
            iz,

        "first_zero_time_ps":
            zt,

        "n_origins_final":
            int(origins[-1]),
    }


# ============================================================
# MAIN
# ============================================================

results = {}
window_rows = []

windows = [
    (0,100),
    (25,100),
    (50,100),
    (75,100),
]


for label, (
    tpr,
    xtc,
) in SYSTEMS.items():

    print()
    print(
        f"LOADING_{label.upper()}=YES"
    )

    dat = load_total_dipole(
        tpr,
        xtc,
        label,
    )

    time = dat["time_ps"]
    M_e_nm = dat["M_e_nm"]
    M_D = dat["M_D"]

    C, origins = (
        connected_vector_corr(
            M_e_nm
        )
    )

    cov, tensor, volume = (
        fluctuation_tensor(
            M_e_nm,
            dat["boxes_nm"],
        )
    )

    iz, zt = first_zero_integral(
        time,
        C,
    )

    mean_M_D = np.mean(
        M_D,
        axis=0,
    )

    sd_M_D = np.std(
        M_D,
        axis=0,
        ddof=1,
    )

    results[label] = {
        "n_water":
            dat["n_water"],

        "mean_M_D":
            mean_M_D.tolist(),

        "sd_M_D":
            sd_M_D.tolist(),

        "mean_M_magnitude_D":
            float(
                np.mean(
                    np.linalg.norm(
                        M_D,
                        axis=1,
                    )
                )
            ),

        "volume_nm3":
            volume,

        "fluctuation_covariance_C2m2":
            cov.tolist(),

        "dimensionless_fluctuation_tensor":
            tensor.tolist(),

        "trace_fluctuation_tensor":
            float(
                np.trace(tensor)
            ),

        "collective_C_2ps":
            float(C[4]),

        "collective_C_5ps":
            float(C[10]),

        "collective_C_10ps":
            float(C[20]),

        "collective_C_20ps":
            float(C[40]),

        "collective_first_zero_integral_ps":
            iz,

        "collective_first_zero_time_ps":
            zt,
    }

    pd.DataFrame({
        "time_ps":
            time,

        "C_M_connected":
            C,

        "n_origins":
            origins,
    }).to_csv(
        OUT
        / f"{label}_collective_C_M.csv",
        index=False,
    )

    pd.DataFrame({
        "time_ps":
            time,

        "Mx_D":
            M_D[:,0],

        "My_D":
            M_D[:,1],

        "Mz_D":
            M_D[:,2],
    }).to_csv(
        OUT
        / f"{label}_total_water_dipole_D.csv",
        index=False,
    )

    for start, stop in windows:

        row = window_metrics(
            time,
            M_e_nm,
            start,
            stop,
        )

        row["system"] = label

        window_rows.append(
            row
        )


window_df = pd.DataFrame(
    window_rows
)

window_df.to_csv(
    OUT
    / "F46O_window_collective_metrics.csv",
    index=False,
)


# ============================================================
# PRINT
# ============================================================

print()
print("="*72)
print("F46O FULL TRAJECTORY COLLECTIVE POLARIZATION")
print("="*72)

for label in [
    "engineered",
    "bulk",
]:

    r = results[label]

    print()
    print(label.upper())

    print(
        f"N_WATER={r['n_water']}"
    )

    print(
        f"VOLUME_NM3="
        f"{r['volume_nm3']:.9f}"
    )

    print(
        "MEAN_M_D="
        + ",".join(
            f"{x:.9f}"
            for x in r[
                "mean_M_D"
            ]
        )
    )

    print(
        "SD_M_D="
        + ",".join(
            f"{x:.9f}"
            for x in r[
                "sd_M_D"
            ]
        )
    )

    print(
        f"TRACE_FLUCTUATION_TENSOR="
        f"{r['trace_fluctuation_tensor']:.12f}"
    )

    tensor = np.asarray(
        r[
            "dimensionless_fluctuation_tensor"
        ]
    )

    for i, axis in enumerate(
        ["X","Y","Z"]
    ):

        print(
            f"FLUCTUATION_{axis}{axis}="
            f"{tensor[i,i]:.12f}"
        )

    print(
        f"C_M_2PS="
        f"{r['collective_C_2ps']:.12f}"
    )

    print(
        f"C_M_5PS="
        f"{r['collective_C_5ps']:.12f}"
    )

    print(
        f"C_M_10PS="
        f"{r['collective_C_10ps']:.12f}"
    )

    print(
        f"C_M_20PS="
        f"{r['collective_C_20ps']:.12f}"
    )

    print(
        f"FIRST_ZERO_INTEGRAL_PS="
        f"{r['collective_first_zero_integral_ps']:.12f}"
    )

    print(
        f"FIRST_ZERO_TIME_PS="
        f"{r['collective_first_zero_time_ps']}"
    )


print()
print("="*72)
print("F46O WINDOWED COLLECTIVE STATIONARITY")
print("="*72)

print(
    window_df[
        [
            "system",
            "start_ps",
            "end_ps",
            "C_2ps",
            "C_5ps",
            "C_10ps",
            "integral_5ps",
            "integral_10ps",
            "first_zero_integral_ps",
            "first_zero_time_ps",
        ]
    ].to_string(
        index=False
    )
)


print()
print("="*72)
print("F46O ENGINEERED/BULK FLUCTUATION RATIOS")
print("="*72)

E = results[
    "engineered"
]

B = results[
    "bulk"
]

print(
    "TRACE_RATIO_ENGINEERED_OVER_BULK="
    f"{E['trace_fluctuation_tensor'] / B['trace_fluctuation_tensor']:.12f}"
)

Et = np.asarray(
    E[
        "dimensionless_fluctuation_tensor"
    ]
)

Bt = np.asarray(
    B[
        "dimensionless_fluctuation_tensor"
    ]
)

for i, axis in enumerate(
    ["X","Y","Z"]
):

    print(
        f"{axis}{axis}_RATIO_ENGINEERED_OVER_BULK="
        f"{Et[i,i]/Bt[i,i]:.12f}"
    )


summary = {
    "gate":
        "F46O",

    "engineered":
        results["engineered"],

    "bulk":
        results["bulk"],

    "window_metrics":
        window_df.to_dict(
            orient="records"
        ),

    "interpretation_boundary": [
        (
            "The fluctuation tensor is a charge-weighted "
            "water polarization fluctuation observable."
        ),
        (
            "For bulk homogeneous water it is related to "
            "dielectric susceptibility under the appropriate "
            "electrostatic ensemble."
        ),
        (
            "For the confined heterogeneous HBN-PYR system "
            "it is not interpreted directly as the dielectric "
            "constant of the whole device."
        ),
        (
            "The frozen solute contributes no dynamical "
            "polarization in this trajectory."
        ),
        (
            "The 0.5 ps sampling limits the Nyquist frequency "
            "to 1 THz."
        ),
    ],
}

with (
    OUT
    / "F46O_summary.json"
).open("w") as fh:

    json.dump(
        summary,
        fh,
        indent=2,
    )


print()
print(
    "SUMMARY_FILE="
    + str(
        OUT
        / "F46O_summary.json"
    )
)

print(
    "F46O_ANALYSIS_COMPLETE=YES"
)
