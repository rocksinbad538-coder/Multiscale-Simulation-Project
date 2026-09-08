#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd

ROOT = Path(".").resolve()

F46R = ROOT / (
    "runs/phase2/day048_f46q_collective_200ps_extension/"
    "F46R_collective_200ps/analyze_F46R.py"
)

R1_NPZ = ROOT / (
    "runs/phase2/day049_f46t_bulk_500ps_convergence/"
    "F46T_bulk_collective_series_0_500.npz"
)

REP_ROOT = ROOT / (
    "runs/phase2/day049_f46w_bulk_independent_replicas"
)

OUT = ROOT / (
    "runs/phase2/day049_f46x_bulk_replica_convergence_final"
)
OUT.mkdir(parents=True, exist_ok=True)

# ======================================================================
# VALIDATED FUNCTIONS
# ======================================================================

spec = importlib.util.spec_from_file_location(
    "f46r_module",
    F46R,
)

f46r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f46r)

print("="*92)
print("F46X FINAL — BULK THREE-REPLICA COLLECTIVE-POLARIZATION CONVERGENCE")
print("Matched duration: 300 ps per replica")
print("="*92)

# ======================================================================
# LOADERS
# ======================================================================

def load_r1():

    d = np.load(R1_NPZ)

    t = np.asarray(
        d["time_ps"],
        dtype=float,
    )

    M = np.asarray(
        d["M"],
        dtype=float,
    )

    boxes = np.asarray(
        d["boxes"],
        dtype=float,
    )

    mask = (
        (t >= 0.0)
        & (t <= 300.0)
    )

    return (
        t[mask],
        M[mask],
        boxes[mask],
        16634,
    )


def load_replica(name):

    d = REP_ROOT/name

    t,frames,nwater = (
        f46r.load_segment(
            d/"production_300ps.tpr",
            d/"production_300ps.xtc",
        )
    )

    t=np.asarray(
        t,
        dtype=float,
    )

    M=np.asarray(
        [
            f["M_e_nm"]
            for f in frames
        ],
        dtype=float,
    )

    boxes=np.asarray(
        [
            f["box_nm"]
            for f in frames
        ],
        dtype=float,
    )

    return (
        t,M,boxes,nwater
    )


systems = {
    "R1": load_r1(),
    "R2": load_replica("R2"),
    "R3": load_replica("R3"),
}

# ======================================================================
# TRAJECTORY GATES
# ======================================================================

overall=True

for name,(t,M,b,nw) in systems.items():

    dt=np.diff(t)

    gate=(
        len(t)==601
        and abs(t[0])<1e-6
        and abs(t[-1]-300.0)<1e-6
        and np.allclose(
            dt,
            0.5,
            atol=1e-6,
            rtol=0,
        )
        and nw==16634
    )

    print()
    print(name)
    print(f"N_WATER={nw}")
    print(f"N_FRAMES={len(t)}")
    print(f"FIRST_PS={t[0]}")
    print(f"LAST_PS={t[-1]}")
    print(
        f"DT_PS="
        f"{np.median(dt)}"
    )
    print(
        f"{name}_300PS_GATE="
        + ("PASS" if gate else "FAIL")
    )

    overall &= gate

print()
print(
    "F46X_MATCHED_TRAJECTORY_GATE="
    + ("PASS" if overall else "FAIL")
)

if not overall:
    raise RuntimeError(
        "Matched 300 ps trajectory validation failed."
    )

# ======================================================================
# TENSOR PER REPLICA
# ======================================================================

records=[]
tensors={}

for name,(t,M,b,nw) in systems.items():

    T,V = f46r.fluctuation_tensor(
        M,b
    )

    met = f46r.tensor_metrics(T)

    tensors[name] = np.asarray(
        T,
        dtype=float,
    )

    records.append({
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
                met[
                    "diag_max_min_ratio"
                ]
            ),
        "volume_nm3":
            float(V),
    })

rep = pd.DataFrame(records)

# ======================================================================
# EQUAL-WEIGHT ENSEMBLE TENSOR
# ======================================================================

Tmean = (
    tensors["R1"]
    + tensors["R2"]
    + tensors["R3"]
)/3.0

mean_met = f46r.tensor_metrics(
    Tmean
)

iso_gate = (
    mean_met["diag_cv"] <= 0.20
    and
    mean_met[
        "diag_max_min_ratio"
    ] <= 1.5
)

# ======================================================================
# BETWEEN-REPLICA STATISTICS
# ======================================================================

stats=[]

for c in [
    "XX","YY","ZZ",
    "XY","XZ","YZ",
    "trace",
]:

    vals=rep[c].to_numpy(
        dtype=float
    )

    mean=float(
        np.mean(vals)
    )

    sd=float(
        np.std(
            vals,
            ddof=1,
        )
    )

    sem=float(
        sd/np.sqrt(3)
    )

    stats.append({
        "observable":c,
        "mean":mean,
        "sd_between_replicas":sd,
        "sem_between_replicas":sem,
        "cv_between_replicas":
            (
                float(sd/abs(mean))
                if mean != 0
                else np.nan
            ),
    })

stats=pd.DataFrame(stats)

# ======================================================================
# 200 ps vs 300 ps ENSEMBLE STABILITY
# ======================================================================

PRE200 = ROOT / (
    "runs/phase2/day049_f46x_preliminary_R1_R2_R3_200ps/"
    "F46X_prelim_200ps_summary.json"
)

if PRE200.exists():

    s200=json.loads(
        PRE200.read_text()
    )

    e200=s200[
        "ensemble_3rep"
    ]

    stability=[]

    for c in [
        "XX","YY","ZZ","trace"
    ]:

        v200=float(
            e200[c]
        )

        v300=float(
            mean_met[c]
        )

        rel=(
            abs(v300-v200)
            / abs(v200)
            *100.0
        )

        stability.append({
            "observable":c,
            "value_200ps":v200,
            "value_300ps":v300,
            "relative_change_percent":
                rel,
        })

    stab=pd.DataFrame(
        stability
    )

else:

    stab=pd.DataFrame()

# ======================================================================
# OFF-DIAGONAL SCALE
# ======================================================================

diag_mean = np.mean([
    abs(mean_met["XX"]),
    abs(mean_met["YY"]),
    abs(mean_met["ZZ"]),
])

offdiag_max = max(
    abs(mean_met["XY"]),
    abs(mean_met["XZ"]),
    abs(mean_met["YZ"]),
)

offdiag_fraction = (
    offdiag_max/diag_mean
)

# ======================================================================
# FINAL DECISION
# ======================================================================

trace_stat = stats[
    stats.observable=="trace"
].iloc[0]

trace_cv = float(
    trace_stat.cv_between_replicas
)

# Primary hard gate remains isotropy.
# Trace CV and 200->300 changes are supporting convergence diagnostics.

bulk_converged = bool(
    iso_gate
)

# ======================================================================
# OUTPUTS
# ======================================================================

rep.to_csv(
    OUT/
    "F46X_final_replica_tensors.csv",
    index=False,
)

stats.to_csv(
    OUT/
    "F46X_final_between_replica_statistics.csv",
    index=False,
)

if not stab.empty:
    stab.to_csv(
        OUT/
        "F46X_final_200_vs_300_stability.csv",
        index=False,
    )

summary={
    "matched_duration_ps":300,
    "replicas":3,

    "ensemble_tensor":{
        k:float(mean_met[k])
        for k in [
            "XX","YY","ZZ",
            "XY","XZ","YZ",
            "trace",
            "diag_cv",
            "diag_max_min_ratio",
        ]
    },

    "ensemble_isotropy_gate":
        bool(iso_gate),

    "trace_between_replica_cv":
        trace_cv,

    "max_offdiag_over_mean_diag":
        float(offdiag_fraction),

    "bulk_collective_converged":
        bulk_converged,

    "decision":
        (
            "BULK_REFERENCE_CONVERGED"
            if bulk_converged
            else
            "BULK_REFERENCE_NOT_CONVERGED"
        ),
}

(
    OUT/"F46X_FINAL_SUMMARY.json"
).write_text(
    json.dumps(
        summary,
        indent=2,
    )
)

# ======================================================================
# PRINT
# ======================================================================

print()
print("="*92)
print("INDIVIDUAL REPLICA TENSORS — MATCHED 300 ps")
print("="*92)

print(
    rep[
        [
            "replica",
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
print("="*92)
print("THREE-REPLICA EQUAL-WEIGHT ENSEMBLE TENSOR")
print("="*92)

for k in [
    "XX","YY","ZZ",
    "XY","XZ","YZ",
    "trace",
    "diag_cv",
    "diag_max_min_ratio",
]:
    print(
        f"ENSEMBLE_{k.upper()}="
        f"{mean_met[k]}"
    )

print(
    "ENSEMBLE_ISOTROPY_GATE="
    + (
        "PASS"
        if iso_gate
        else "FAIL"
    )
)

print()
print("="*92)
print("BETWEEN-REPLICA STATISTICS")
print("="*92)

print(
    stats.to_string(
        index=False
    )
)

print()
print(
    "TRACE_BETWEEN_REPLICA_CV="
    f"{trace_cv}"
)

print(
    "MAX_OFFDIAG_OVER_MEAN_DIAG="
    f"{offdiag_fraction}"
)

if not stab.empty:

    print()
    print("="*92)
    print("ENSEMBLE STABILITY — 200 ps vs 300 ps")
    print("="*92)

    print(
        stab.to_string(
            index=False
        )
    )

    for _,r in stab.iterrows():

        print(
            f"CHANGE_200_TO_300_"
            f"{r.observable}_PERCENT="
            f"{r.relative_change_percent}"
        )

print()
print("="*92)
print("F46X FINAL DECISION")
print("="*92)

print(
    "F46X_BULK_COLLECTIVE_CONVERGENCE="
    + (
        "PASS"
        if bulk_converged
        else "FAIL"
    )
)

print(
    "F46X_DECISION="
    + summary["decision"]
)

print(
    "F46X_FINAL_ANALYSIS_COMPLETE=YES"
)
