#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
import json
import MDAnalysis as mda

ROOT = Path(".").resolve()

F46R = ROOT / (
    "runs/phase2/day048_f46q_collective_200ps_extension/"
    "F46R_collective_200ps/analyze_F46R.py"
)

R1_NPZ = ROOT / (
    "runs/phase2/day049_f46t_bulk_500ps_convergence/"
    "F46T_bulk_collective_series_0_500.npz"
)

R2_TPR = ROOT / (
    "runs/phase2/day049_f46w_bulk_independent_replicas/"
    "R2/production_300ps.tpr"
)

R2_XTC = ROOT / (
    "runs/phase2/day049_f46w_bulk_independent_replicas/"
    "R2/production_300ps.xtc"
)

OUT = ROOT / "runs/phase2/day049_f46x_preliminary_R1_R2"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Reuse validated F46R physics
# ------------------------------------------------------------------

spec = importlib.util.spec_from_file_location(
    "f46r_module",
    F46R,
)
f46r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f46r)

print("="*84)
print("F46X PRELIMINARY — R1/R2 BULK REPLICA COMPARISON")
print("Matched 300 ps production duration")
print("="*84)

# ------------------------------------------------------------------
# R1: use first 300 ps of validated 0–500 ps trajectory
# ------------------------------------------------------------------

d1 = np.load(R1_NPZ)

t1 = np.asarray(d1["time_ps"], dtype=float)
M1 = np.asarray(d1["M"], dtype=float)
b1 = np.asarray(d1["boxes"], dtype=float)

mask1 = (
    (t1 >= 0.0)
    & (t1 <= 300.0)
)

t1 = t1[mask1]
M1 = M1[mask1]
b1 = b1[mask1]

print()
print("R1")
print(f"N_FRAMES={len(t1)}")
print(f"FIRST_PS={t1[0]}")
print(f"LAST_PS={t1[-1]}")
print(f"DURATION_PS={t1[-1]-t1[0]}")
print(f"DT_PS={np.median(np.diff(t1))}")

r1_gate = (
    len(t1) == 601
    and abs((t1[-1]-t1[0])-300.0) < 1e-6
    and np.allclose(np.diff(t1),0.5,atol=1e-6,rtol=0)
)

print(
    "R1_MATCHED_300PS_GATE="
    + ("PASS" if r1_gate else "FAIL")
)

# ------------------------------------------------------------------
# R2: load completed independent production
# ------------------------------------------------------------------

t2, frames2, nw2 = f46r.load_segment(
    R2_TPR,
    R2_XTC,
)

t2 = np.asarray(t2,dtype=float)

M2 = np.asarray(
    [x["M_e_nm"] for x in frames2],
    dtype=float,
)

b2 = np.asarray(
    [x["box_nm"] for x in frames2],
    dtype=float,
)

print()
print("R2")
print(f"N_WATER={nw2}")
print(f"N_FRAMES={len(t2)}")
print(f"FIRST_PS={t2[0]}")
print(f"LAST_PS={t2[-1]}")
print(f"DURATION_PS={t2[-1]-t2[0]}")
print(f"DT_PS={np.median(np.diff(t2))}")

r2_gate = (
    len(t2) == 601
    and abs((t2[-1]-t2[0])-300.0) < 1e-6
    and np.allclose(np.diff(t2),0.5,atol=1e-6,rtol=0)
)

print(
    "R2_300PS_GATE="
    + ("PASS" if r2_gate else "FAIL")
)

if not (r1_gate and r2_gate):
    raise RuntimeError(
        "Matched-duration trajectory gate failed."
    )

# ------------------------------------------------------------------
# Tensor helper
# ------------------------------------------------------------------

def tensor_record(replica,M,boxes):

    T,V = f46r.fluctuation_tensor(
        M,
        boxes,
    )

    met = f46r.tensor_metrics(T)

    return {
        "replica": replica,
        "XX": float(met["XX"]),
        "YY": float(met["YY"]),
        "ZZ": float(met["ZZ"]),
        "XY": float(met["XY"]),
        "XZ": float(met["XZ"]),
        "YZ": float(met["YZ"]),
        "trace": float(met["trace"]),
        "diag_cv": float(met["diag_cv"]),
        "diag_max_min_ratio":
            float(met["diag_max_min_ratio"]),
        "volume_nm3": float(V),
    }, np.asarray(T,dtype=float)

r1, T1 = tensor_record(
    "R1_300ps",
    M1,
    b1,
)

r2, T2 = tensor_record(
    "R2_300ps",
    M2,
    b2,
)

# ------------------------------------------------------------------
# Preliminary equal-weight replica estimator
# ------------------------------------------------------------------

Tmean = 0.5 * (T1 + T2)
mean_met = f46r.tensor_metrics(Tmean)

repdf = pd.DataFrame([r1,r2])

components = ["XX","YY","ZZ","trace"]

stats = []

for c in components:

    vals = repdf[c].to_numpy(dtype=float)

    stats.append({
        "observable": c,
        "mean_R1_R2": float(np.mean(vals)),
        "sd_between_replicas":
            float(np.std(vals,ddof=1)),
        "sem_between_replicas":
            float(np.std(vals,ddof=1)/np.sqrt(2)),
        "relative_difference_percent":
            float(
                abs(vals[1]-vals[0])
                / abs(np.mean(vals))
                * 100.0
            ),
    })

stats = pd.DataFrame(stats)

# ------------------------------------------------------------------
# Print
# ------------------------------------------------------------------

print()
print("="*84)
print("INDIVIDUAL 300 ps REPLICA TENSORS")
print("="*84)

print(
    repdf[
        [
            "replica",
            "XX","YY","ZZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    ].to_string(index=False)
)

print()
print("="*84)
print("PRELIMINARY R1/R2 EQUAL-WEIGHT TENSOR")
print("="*84)

for k in [
    "XX","YY","ZZ",
    "XY","XZ","YZ",
    "trace",
    "diag_cv",
    "diag_max_min_ratio",
]:
    print(
        f"R1_R2_MEAN_{k.upper()}="
        f"{mean_met[k]}"
    )

prelim_gate = (
    mean_met["diag_cv"] <= 0.20
    and
    mean_met["diag_max_min_ratio"] <= 1.5
)

print(
    "R1_R2_PRELIM_ISOTROPY_GATE="
    + ("PASS" if prelim_gate else "FAIL")
)

print()
print("="*84)
print("BETWEEN-REPLICA DIAGNOSTICS")
print("="*84)

print(
    stats.to_string(index=False)
)

repdf.to_csv(
    OUT/"F46X_preliminary_replica_tensors.csv",
    index=False,
)

stats.to_csv(
    OUT/"F46X_preliminary_between_replica_stats.csv",
    index=False,
)

summary = {
    "matched_duration_ps": 300,
    "R1": r1,
    "R2": r2,
    "R1_R2_mean_tensor": {
        k: float(mean_met[k])
        for k in [
            "XX","YY","ZZ",
            "XY","XZ","YZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    },
    "preliminary_isotropy_gate":
        bool(prelim_gate),
    "status":
        "PRELIMINARY_R3_PENDING",
}

(
    OUT/"F46X_preliminary_summary.json"
).write_text(
    json.dumps(summary,indent=2)
)

print()
print("F46X_PRELIMINARY_STATUS=R1_R2_COMPLETE_R3_PENDING")
print("F46X_PRELIMINARY_ANALYSIS_COMPLETE=YES")
