#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd

ROOT = Path(".").resolve()

SRC = ROOT / (
    "runs/phase2/day049_f46t_bulk_500ps_convergence/"
    "F46T_bulk_collective_series_0_500.npz"
)

F46R = ROOT / (
    "runs/phase2/day048_f46q_collective_200ps_extension/"
    "F46R_collective_200ps/analyze_F46R.py"
)

OUT = ROOT / "runs/phase2/day049_f46u_bulk_sampling_strategy"
OUT.mkdir(parents=True, exist_ok=True)

# ================================================================
# Reuse validated tensor implementation
# ================================================================

spec = importlib.util.spec_from_file_location(
    "f46r_module",
    F46R,
)

f46r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f46r)

dat = np.load(SRC)

time = np.asarray(
    dat["time_ps"],
    dtype=float,
)

M = np.asarray(
    dat["M"],
    dtype=float,
)

boxes = np.asarray(
    dat["boxes"],
    dtype=float,
)

dt = float(
    np.median(np.diff(time))
)

print("="*84)
print("F46U — BULK COLLECTIVE SAMPLING-STRATEGY DIAGNOSTIC")
print("="*84)

print(f"N_FRAMES={len(time)}")
print(f"FIRST_PS={time[0]}")
print(f"LAST_PS={time[-1]}")
print(f"DT_PS={dt}")

gate = (
    len(time) == 1001
    and abs(time[0]) < 1e-6
    and abs(time[-1]-500.0) < 1e-6
    and np.allclose(
        np.diff(time),
        0.5,
        atol=1e-6,
        rtol=0,
    )
)

print(
    "INPUT_TRAJECTORY_GATE="
    + ("PASS" if gate else "FAIL")
)

if not gate:
    raise RuntimeError(
        "F46U input trajectory invalid."
    )

# ================================================================
# Helpers
# ================================================================

def tensor_window(start, end):

    mask = (
        (time >= start)
        & (time <= end)
    )

    T,V = f46r.fluctuation_tensor(
        M[mask],
        boxes[mask],
    )

    met = f46r.tensor_metrics(T)

    return {
        "start_ps": float(start),
        "end_ps": float(end),
        "duration_ps": float(end-start),
        "n_frames": int(mask.sum()),
        "XX": float(met["XX"]),
        "YY": float(met["YY"]),
        "ZZ": float(met["ZZ"]),
        "trace": float(met["trace"]),
        "diag_cv": float(
            met["diag_cv"]
        ),
        "diag_max_min_ratio":
            float(
                met[
                    "diag_max_min_ratio"
                ]
            ),
    }


def acf_fft(x):

    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)

    n = len(x)

    nfft = 1 << (
        2*n - 1
    ).bit_length()

    fx = np.fft.rfft(
        x,
        n=nfft,
    )

    ac = np.fft.irfft(
        fx*np.conjugate(fx),
        n=nfft,
    )[:n]

    norm = np.arange(
        n,0,-1,
        dtype=float,
    )

    ac = ac / norm

    if ac[0] <= 0:
        return np.ones(n)

    return ac/ac[0]


def tau_initial_positive(acf, dt):

    """
    Integrated autocorrelation time using
    the initial-positive truncation.

    tau_int = dt * (1/2 + sum_{k>=1} rho_k)
    """

    s = 0.5

    cutoff = 0

    for k in range(
        1,len(acf)
    ):
        if not np.isfinite(acf[k]):
            break

        if acf[k] <= 0:
            break

        s += acf[k]
        cutoff = k

    return (
        float(dt*s),
        int(cutoff),
    )


def acf_metrics(name,x):

    ac = acf_fft(x)

    tau,cut = (
        tau_initial_positive(
            ac,
            dt,
        )
    )

    total_duration = (
        time[-1]-time[0]
    )

    # Approximate number of statistically
    # independent observations.
    neff = (
        total_duration
        / (2.0*tau)
        if tau > 0
        else np.nan
    )

    def at(tps):
        i = int(
            round(tps/dt)
        )
        i = min(
            i,
            len(ac)-1,
        )
        return float(ac[i])

    return {
        "observable": name,
        "tau_int_ps": tau,
        "positive_cutoff_ps":
            cut*dt,
        "effective_samples":
            float(neff),
        "acf_2ps": at(2),
        "acf_5ps": at(5),
        "acf_10ps": at(10),
        "acf_20ps": at(20),
        "acf_50ps": at(50),
    }, ac


# ================================================================
# Rolling tensor windows
# ================================================================

rolling100=[]

for start in np.arange(
    0,400.0001,25.0
):
    rolling100.append(
        tensor_window(
            start,
            start+100,
        )
    )

rolling200=[]

for start in np.arange(
    0,300.0001,25.0
):
    rolling200.append(
        tensor_window(
            start,
            start+200,
        )
    )

r100 = pd.DataFrame(
    rolling100
)

r200 = pd.DataFrame(
    rolling200
)

r100.to_csv(
    OUT/"F46U_rolling_100ps_tensor.csv",
    index=False,
)

r200.to_csv(
    OUT/"F46U_rolling_200ps_tensor.csv",
    index=False,
)

# ================================================================
# Mean-polarization windows
# ================================================================

mean_rows=[]

for start in range(
    0,500,50
):

    end=start+50

    mask=(
        (time>=start)
        & (time<=end)
    )

    mm=np.mean(
        M[mask],
        axis=0,
    )

    mean_rows.append({
        "start_ps":start,
        "end_ps":end,
        "Mx_mean":float(mm[0]),
        "My_mean":float(mm[1]),
        "Mz_mean":float(mm[2]),
        "Mnorm_mean":
            float(np.linalg.norm(mm)),
    })

means=pd.DataFrame(mean_rows)

means.to_csv(
    OUT/"F46U_mean_polarization_50ps_windows.csv",
    index=False,
)

# ================================================================
# Autocorrelation of M and tensor-building observables
# ================================================================

dM = (
    M - np.mean(M,axis=0)
)

observables = {
    "Mx": dM[:,0],
    "My": dM[:,1],
    "Mz": dM[:,2],

    "Mx2": dM[:,0]**2,
    "My2": dM[:,1]**2,
    "Mz2": dM[:,2]**2,

    "MxMy": dM[:,0]*dM[:,1],
    "MxMz": dM[:,0]*dM[:,2],
    "MyMz": dM[:,1]*dM[:,2],
}

acf_rows=[]
acf_series={}

for name,x in observables.items():

    row,ac = acf_metrics(
        name,x
    )

    acf_rows.append(row)

    acf_series[name]=ac

acf_df = pd.DataFrame(
    acf_rows
)

acf_df.to_csv(
    OUT/"F46U_autocorrelation_times.csv",
    index=False,
)

maxlag = min(
    len(time),
    int(100/dt)+1,
)

acf_out = pd.DataFrame({
    "lag_ps":
        np.arange(maxlag)*dt
})

for name,ac in acf_series.items():
    acf_out[name] = (
        ac[:maxlag]
    )

acf_out.to_csv(
    OUT/"F46U_tensor_observable_acf.csv",
    index=False,
)

# ================================================================
# Late-window stability
# ================================================================

late100 = r100[
    r100.start_ps >= 200
].copy()

late200 = r200[
    r200.start_ps >= 100
].copy()

def relative_cv(series):

    x=np.asarray(
        series,
        dtype=float,
    )

    return float(
        np.std(x,ddof=1)
        / abs(np.mean(x))
    )


rolling100_cv = {
    q:relative_cv(
        late100[q]
    )
    for q in [
        "XX","YY","ZZ","trace"
    ]
}

rolling200_cv = {
    q:relative_cv(
        late200[q]
    )
    for q in [
        "XX","YY","ZZ","trace"
    ]
}

# Tensor-building effective samples
tensor_obs = acf_df[
    acf_df.observable.isin(
        ["Mx2","My2","Mz2"]
    )
]

min_tensor_neff = float(
    tensor_obs[
        "effective_samples"
    ].min()
)

max_tensor_tau = float(
    tensor_obs[
        "tau_int_ps"
    ].max()
)

# ================================================================
# Decision logic
# ================================================================

#
# This is an operational recommendation,
# not a physical theorem.
#
# If diagonal variance observables have very
# long autocorrelation / very low Neff, one
# longer trajectory remains valuable.
#
# If Neff is not extremely small but rolling
# estimates remain highly variable, independent
# replicas are preferable for uncertainty.
#

if (
    min_tensor_neff < 15
    or max_tensor_tau > 20
):
    recommendation = (
        "LONGER_TRAJECTORY_OR_LONG_REPLICAS_REQUIRED"
    )

elif (
    max(
        rolling200_cv.values()
    ) > 0.20
):
    recommendation = (
        "INDEPENDENT_REPLICAS_PREFERRED"
    )

else:
    recommendation = (
        "SAMPLING_CLOSE_TO_CONVERGENCE"
    )

# ================================================================
# Print
# ================================================================

print()
print("="*84)
print("F46U ROLLING 100 ps TENSORS")
print("="*84)

print(
    r100[
        [
            "start_ps",
            "end_ps",
            "XX","YY","ZZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    ].to_string(
        index=False
    )
)

print()
print("="*84)
print("F46U ROLLING 200 ps TENSORS")
print("="*84)

print(
    r200[
        [
            "start_ps",
            "end_ps",
            "XX","YY","ZZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    ].to_string(
        index=False
    )
)

print()
print("="*84)
print("F46U POLARIZATION MEAN — 50 ps WINDOWS")
print("="*84)

print(
    means.to_string(
        index=False
    )
)

print()
print("="*84)
print("F46U AUTOCORRELATION / EFFECTIVE SAMPLE SIZE")
print("="*84)

print(
    acf_df.to_string(
        index=False
    )
)

print()
print("="*84)
print("F46U LATE-WINDOW STABILITY")
print("="*84)

for k,v in rolling100_cv.items():
    print(
        f"LATE_ROLLING100_{k}_CV={v}"
    )

for k,v in rolling200_cv.items():
    print(
        f"LATE_ROLLING200_{k}_CV={v}"
    )

print(
    f"MAX_TENSOR_OBSERVABLE_TAU_INT_PS="
    f"{max_tensor_tau}"
)

print(
    f"MIN_TENSOR_OBSERVABLE_EFFECTIVE_SAMPLES="
    f"{min_tensor_neff}"
)

print()
print(
    "F46U_SAMPLING_RECOMMENDATION="
    + recommendation
)

summary = {
    "max_tensor_observable_tau_int_ps":
        max_tensor_tau,
    "min_tensor_observable_effective_samples":
        min_tensor_neff,
    "late_rolling100_relative_cv":
        rolling100_cv,
    "late_rolling200_relative_cv":
        rolling200_cv,
    "sampling_recommendation":
        recommendation,
}

(OUT/"F46U_summary.json").write_text(
    json.dumps(
        summary,
        indent=2,
    )
)

print("F46U_ANALYSIS_COMPLETE=YES")
