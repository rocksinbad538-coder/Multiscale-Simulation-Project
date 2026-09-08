#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
import json

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

R3_TPR = ROOT / (
    "runs/phase2/day049_f46w_bulk_independent_replicas/"
    "R3/production_300ps.tpr"
)
R3_XTC = ROOT / (
    "runs/phase2/day049_f46w_bulk_independent_replicas/"
    "R3/production_300ps.xtc"
)

OUT = ROOT / "runs/phase2/day049_f46x_preliminary_R1_R2_R3_200ps"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# validated F46R functions
# ------------------------------------------------------------------

spec = importlib.util.spec_from_file_location(
    "f46r_module",
    F46R,
)

f46r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f46r)

print("="*88)
print("F46X PRELIMINARY — R1/R2/R3 MATCHED 200 ps BULK COMPARISON")
print("="*88)

def load_from_npz_200():

    d = np.load(R1_NPZ)

    t = np.asarray(d["time_ps"], dtype=float)
    M = np.asarray(d["M"], dtype=float)
    boxes = np.asarray(d["boxes"], dtype=float)

    mask = (
        (t >= 0.0)
        & (t <= 200.0)
    )

    return (
        t[mask],
        M[mask],
        boxes[mask],
    )

def load_xtc_200(tpr, xtc):

    t, frames, nw = f46r.load_segment(
        tpr,
        xtc,
    )

    t = np.asarray(t, dtype=float)

    M = np.asarray(
        [x["M_e_nm"] for x in frames],
        dtype=float,
    )

    boxes = np.asarray(
        [x["box_nm"] for x in frames],
        dtype=float,
    )

    # R3 may still be actively writing.
    mask = (
        (t >= 0.0)
        & (t <= 200.0)
    )

    return (
        t[mask],
        M[mask],
        boxes[mask],
        nw,
    )

def gate_200(label,t):

    ok = (
        len(t) == 401
        and abs(t[0]) < 1e-6
        and abs(t[-1]-200.0) < 1e-6
        and np.allclose(
            np.diff(t),
            0.5,
            atol=1e-6,
            rtol=0,
        )
    )

    print()
    print(label)
    print(f"N_FRAMES={len(t)}")
    print(f"FIRST_PS={t[0] if len(t) else 'NA'}")
    print(f"LAST_PS={t[-1] if len(t) else 'NA'}")

    if len(t) > 1:
        print(
            f"DT_PS={np.median(np.diff(t))}"
        )

    print(
        f"{label}_200PS_GATE="
        + ("PASS" if ok else "FAIL")
    )

    return ok

t1,M1,b1 = load_from_npz_200()
t2,M2,b2,nw2 = load_xtc_200(
    R2_TPR,R2_XTC
)
t3,M3,b3,nw3 = load_xtc_200(
    R3_TPR,R3_XTC
)

g1 = gate_200("R1",t1)
g2 = gate_200("R2",t2)
g3 = gate_200("R3",t3)

if not (g1 and g2 and g3):
    raise RuntimeError(
        "Matched 200 ps gate failed."
    )

if nw2 != nw3:
    raise RuntimeError(
        "R2/R3 water counts differ."
    )

def tensor_rec(name,M,b):

    T,V = f46r.fluctuation_tensor(
        M,b
    )

    met = f46r.tensor_metrics(T)

    return (
        {
            "replica":name,
            "XX":float(met["XX"]),
            "YY":float(met["YY"]),
            "ZZ":float(met["ZZ"]),
            "XY":float(met["XY"]),
            "XZ":float(met["XZ"]),
            "YZ":float(met["YZ"]),
            "trace":float(met["trace"]),
            "diag_cv":
                float(met["diag_cv"]),
            "diag_max_min_ratio":
                float(
                    met["diag_max_min_ratio"]
                ),
        },
        np.asarray(T,dtype=float)
    )

r1,T1 = tensor_rec("R1_200ps",M1,b1)
r2,T2 = tensor_rec("R2_200ps",M2,b2)
r3,T3 = tensor_rec("R3_200ps",M3,b3)

rep = pd.DataFrame(
    [r1,r2,r3]
)

T12 = (T1+T2)/2.0
T123 = (T1+T2+T3)/3.0

m12 = f46r.tensor_metrics(T12)
m123 = f46r.tensor_metrics(T123)

def gate(m):
    return (
        m["diag_cv"] <= 0.20
        and
        m["diag_max_min_ratio"] <= 1.5
    )

print()
print("="*88)
print("INDIVIDUAL 200 ps TENSORS")
print("="*88)

print(
    rep[
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
print("="*88)
print("R1+R2 MATCHED 200 ps")
print("="*88)

for k in [
    "XX","YY","ZZ",
    "trace",
    "diag_cv",
    "diag_max_min_ratio",
]:
    print(
        f"R1_R2_200_{k.upper()}="
        f"{m12[k]}"
    )

print(
    "R1_R2_200_ISOTROPY_GATE="
    + ("PASS" if gate(m12) else "FAIL")
)

print()
print("="*88)
print("R1+R2+R3 MATCHED 200 ps")
print("="*88)

for k in [
    "XX","YY","ZZ",
    "XY","XZ","YZ",
    "trace",
    "diag_cv",
    "diag_max_min_ratio",
]:
    print(
        f"R1_R2_R3_200_{k.upper()}="
        f"{m123[k]}"
    )

final_gate = gate(m123)

print(
    "R1_R2_R3_200_ISOTROPY_GATE="
    + ("PASS" if final_gate else "FAIL")
)

# Between-replica stats
rows=[]

for c in [
    "XX","YY","ZZ","trace"
]:
    vals=rep[c].to_numpy()

    rows.append({
        "observable":c,
        "mean":float(np.mean(vals)),
        "sd":float(np.std(vals,ddof=1)),
        "sem":float(
            np.std(vals,ddof=1)
            / np.sqrt(3)
        ),
        "cv_between_replicas":
            float(
                np.std(vals,ddof=1)
                / abs(np.mean(vals))
            ),
    })

stats=pd.DataFrame(rows)

print()
print("="*88)
print("BETWEEN-REPLICA STATISTICS")
print("="*88)

print(
    stats.to_string(index=False)
)

rep.to_csv(
    OUT/"F46X_prelim_200ps_replica_tensors.csv",
    index=False,
)

stats.to_csv(
    OUT/"F46X_prelim_200ps_replica_stats.csv",
    index=False,
)

summary={
    "matched_duration_ps":200,
    "ensemble_3rep":{
        k:float(m123[k])
        for k in [
            "XX","YY","ZZ",
            "XY","XZ","YZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    },
    "ensemble_isotropy_gate":
        bool(final_gate),
    "status":
        "PRELIMINARY_R3_PRODUCTION_STILL_RUNNING",
}

(
    OUT/"F46X_prelim_200ps_summary.json"
).write_text(
    json.dumps(
        summary,
        indent=2,
    )
)

print()
print(
    "F46X_PRELIM_200PS_STATUS="
    "R1_R2_R3_COMPLETE_AT_MATCHED_200PS"
)
print(
    "FINAL_300PS_ANALYSIS_PENDING=YES"
)
