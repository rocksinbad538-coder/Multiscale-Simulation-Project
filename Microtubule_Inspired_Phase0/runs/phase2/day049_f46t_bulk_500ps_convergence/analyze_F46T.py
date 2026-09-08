#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd

# ======================================================================
# PATHS
# ======================================================================

ROOT = Path(".").resolve()

F46R_PATH = ROOT / (
    "runs/phase2/day048_f46q_collective_200ps_extension/"
    "F46R_collective_200ps/analyze_F46R.py"
)

OUT = ROOT / "runs/phase2/day049_f46t_bulk_500ps_convergence"
OUT.mkdir(parents=True, exist_ok=True)

TPR = ROOT / (
    "runs/phase2/day048_f46l_bulk_100ps_orientation/"
    "nvt_bulk_100ps.tpr"
)

SEGMENTS = [
    (
        "0_100",
        ROOT / (
            "runs/phase2/day048_f46l_bulk_100ps_orientation/"
            "nvt_bulk_100ps.xtc"
        ),
    ),
    (
        "100_200",
        ROOT / (
            "runs/phase2/day048_f46q_collective_200ps_extension/"
            "bulk/bulk_extend_100ps.part0002.xtc"
        ),
    ),
    (
        "200_500",
        ROOT / (
            "runs/phase2/day048_f46s_bulk_500ps/"
            "bulk_500ps.part0003.xtc"
        ),
    ),
]

# ======================================================================
# LOAD VALIDATED F46R FUNCTIONS
# ======================================================================

spec = importlib.util.spec_from_file_location(
    "f46r_module",
    F46R_PATH,
)

f46r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f46r)

required = [
    "load_segment",
    "fluctuation_tensor",
    "tensor_metrics",
    "connected_corr",
    "first_zero_integral",
    "integral_to",
]

for name in required:
    if not hasattr(f46r, name):
        raise RuntimeError(
            f"Required validated F46R function missing: {name}"
        )

print("="*84)
print("F46T — BULK 0–500 ps COLLECTIVE-POLARIZATION CONVERGENCE")
print("="*84)
print("F46R_FUNCTION_REUSE_GATE=PASS")

# ======================================================================
# LOAD THE THREE SEGMENTS
# ======================================================================

segments = []

for label, xtc in SEGMENTS:

    print()
    print(f"LOADING_SEGMENT={label}")

    result = f46r.load_segment(
        TPR,
        xtc,
    )

    if not isinstance(result, tuple):
        raise RuntimeError(
            f"{label}: unexpected load_segment return type "
            f"{type(result)}"
        )

    if len(result) != 3:
        raise RuntimeError(
            f"{label}: expected load_segment to return 3 values, "
            f"got {len(result)}"
        )

    time, frames, nwater = result

    time = np.asarray(
        time,
        dtype=float,
    )

    M = np.asarray(
        [
            frame["M_e_nm"]
            for frame in frames
        ],
        dtype=float,
    )

    boxes = np.asarray(
        [
            frame["box_nm"]
            for frame in frames
        ],
        dtype=float,
    )

    if not (
        len(time)
        == len(M)
        == len(boxes)
        == len(frames)
    ):
        raise RuntimeError(
            f"{label}: inconsistent segment lengths "
            f"time={len(time)} "
            f"M={len(M)} "
            f"boxes={len(boxes)} "
            f"frames={len(frames)}"
        )

    print(f"N_WATER={nwater}")
    print(f"N_FRAMES={len(time)}")
    print(f"FIRST_PS={time[0]}")
    print(f"LAST_PS={time[-1]}")

    if len(time) > 1:
        print(
            f"DT_PS="
            f"{np.median(np.diff(time))}"
        )

    segments.append(
        {
            "label": label,
            "time": time,
            "M": M,
            "boxes": boxes,
            "nwater": nwater,
        }
    )

# ======================================================================
# CONSISTENCY
# ======================================================================

nw = {
    s["nwater"]
    for s in segments
}

if len(nw) != 1:
    raise RuntimeError(
        f"Water count inconsistent: {nw}"
    )

nwater = list(nw)[0]

# Concatenate.
time = np.concatenate(
    [s["time"] for s in segments]
)

M = np.concatenate(
    [s["M"] for s in segments],
    axis=0,
)

boxes = np.concatenate(
    [s["boxes"] for s in segments],
    axis=0,
)

# Deduplicate the shared 100 ps and 200 ps frames.
order = np.argsort(time)

time = time[order]
M = M[order]
boxes = boxes[order]

rounded = np.round(
    time,
    decimals=6,
)

_, unique_idx = np.unique(
    rounded,
    return_index=True,
)

unique_idx = np.sort(unique_idx)

time = time[unique_idx]
M = M[unique_idx]
boxes = boxes[unique_idx]

print()
print("="*84)
print("F46T COMBINED TRAJECTORY")
print("="*84)

print(f"BULK_N_WATER={nwater}")
print(f"BULK_N_FRAMES={len(time)}")
print(f"BULK_FIRST_PS={time[0]}")
print(f"BULK_LAST_PS={time[-1]}")

dt = np.diff(time)

print(
    f"BULK_DT_PS="
    f"{np.median(dt)}"
)

complete = (
    len(time) == 1001
    and abs(time[0] - 0.0) < 1e-4
    and abs(time[-1] - 500.0) < 1e-4
    and np.allclose(
        dt,
        0.5,
        atol=1e-4,
        rtol=0,
    )
)

print(
    "BULK_0_500_TRAJECTORY_GATE="
    + ("PASS" if complete else "FAIL")
)

if not complete:
    raise RuntimeError(
        "Combined bulk 0–500 ps trajectory failed validation."
    )

# ======================================================================
# HELPERS
# ======================================================================

def tensor_row(
    window_type,
    start,
    end,
):

    mask = (
        (time >= start)
        & (time <= end)
    )

    t = time[mask]
    x = M[mask]
    b = boxes[mask]

    T, V = f46r.fluctuation_tensor(
        x,
        b,
    )

    met = f46r.tensor_metrics(T)

    return {
        "window_type": window_type,
        "start_ps": start,
        "end_ps": end,
        "duration_ps": end-start,
        "n_frames": len(t),
        "XX": met["XX"],
        "YY": met["YY"],
        "ZZ": met["ZZ"],
        "XY": met["XY"],
        "XZ": met["XZ"],
        "YZ": met["YZ"],
        "trace": met["trace"],
        "diag_cv": met["diag_cv"],
        "diag_max_min_ratio":
            met["diag_max_min_ratio"],
        "volume_nm3": V,
    }


def corr_row(
    start,
    end,
):

    mask = (
        (time >= start)
        & (time <= end)
    )

    t_abs = time[mask]
    x = M[mask]

    t = (
        t_abs
        - t_abs[0]
    )

    C, origins = f46r.connected_corr(
        x
    )

    iz, zt = f46r.first_zero_integral(
        t,
        C,
    )

    def at(xps):
        i = np.argmin(
            np.abs(t - xps)
        )
        return float(C[i])

    return {
        "start_ps": start,
        "end_ps": end,
        "duration_ps": end-start,
        "n_frames": len(t),
        "C_2ps": at(2),
        "C_5ps": at(5),
        "C_10ps": at(10),
        "integral_5ps":
            f46r.integral_to(
                t,C,5
            ),
        "integral_10ps":
            f46r.integral_to(
                t,C,10
            ),
        "first_zero_integral_ps":
            iz,
        "first_zero_time_ps":
            zt,
    }

# ======================================================================
# CUMULATIVE TENSORS
# ======================================================================

cumulative_ends = [
    100,
    200,
    300,
    400,
    500,
]

cum_rows = [
    tensor_row(
        "cumulative",
        0,
        end,
    )
    for end in cumulative_ends
]

cum = pd.DataFrame(cum_rows)

# ======================================================================
# INDEPENDENT 100 ps BLOCKS
# ======================================================================

blocks = [
    (0,100),
    (100,200),
    (200,300),
    (300,400),
    (400,500),
]

block_rows = [
    tensor_row(
        "block100",
        start,
        end,
    )
    for start,end in blocks
]

block = pd.DataFrame(
    block_rows
)

# ======================================================================
# LATE 200 ps AND 300 ps WINDOWS
# ======================================================================

late_windows = [
    (200,500),
    (300,500),
    (350,500),
    (400,500),
]

late_rows = [
    tensor_row(
        "late_window",
        start,
        end,
    )
    for start,end in late_windows
]

late = pd.DataFrame(
    late_rows
)

# ======================================================================
# COLLECTIVE CORRELATION WINDOWS
# ======================================================================

corr_windows = [
    (0,500),
    (100,500),
    (200,500),
    (300,500),
    (400,500),
]

corr = pd.DataFrame(
    [
        corr_row(start,end)
        for start,end
        in corr_windows
    ]
)

# ======================================================================
# 500 ps GATE
# ======================================================================

final = cum[
    cum.end_ps == 500
].iloc[0]

cv_gate = (
    final.diag_cv
    <= 0.20
)

ratio_gate = (
    final.diag_max_min_ratio
    <= 1.5
)

isotropy_gate = (
    cv_gate
    and ratio_gate
)

# ----------------------------------------------------------------------
# Secondary convergence diagnostics.
# These are diagnostics, not substitutes for the frozen isotropy gate.
# ----------------------------------------------------------------------

last3 = block[
    block.start_ps >= 200
].copy()

trace_mean = float(
    last3["trace"].mean()
)

trace_sd = float(
    last3["trace"].std(
        ddof=1
    )
)

trace_cv_late = (
    trace_sd
    / abs(trace_mean)
    if trace_mean != 0
    else np.nan
)

# Component variability across 200–300, 300–400, 400–500 blocks.
component_block_cv = {}

for comp in [
    "XX","YY","ZZ"
]:
    vals = last3[comp].to_numpy()

    component_block_cv[comp] = float(
        np.std(
            vals,
            ddof=1,
        )
        / abs(
            np.mean(vals)
        )
    )

# ======================================================================
# OUTPUT
# ======================================================================

cum.to_csv(
    OUT/"F46T_bulk_cumulative_tensor.csv",
    index=False,
)

block.to_csv(
    OUT/"F46T_bulk_100ps_blocks.csv",
    index=False,
)

late.to_csv(
    OUT/"F46T_bulk_late_windows.csv",
    index=False,
)

corr.to_csv(
    OUT/"F46T_bulk_collective_correlation_windows.csv",
    index=False,
)

np.savez_compressed(
    OUT/"F46T_bulk_collective_series_0_500.npz",
    time_ps=time,
    M=M,
    boxes=boxes,
)

summary = {
    "trajectory": {
        "nwater": int(nwater),
        "frames": int(len(time)),
        "first_ps": float(time[0]),
        "last_ps": float(time[-1]),
        "dt_ps": float(
            np.median(dt)
        ),
    },
    "bulk_500ps": {
        "XX": float(final.XX),
        "YY": float(final.YY),
        "ZZ": float(final.ZZ),
        "trace": float(
            final["trace"]
        ),
        "diag_cv": float(
            final.diag_cv
        ),
        "diag_max_min_ratio":
            float(
                final.diag_max_min_ratio
            ),
    },
    "gates": {
        "diag_cv_le_0p20":
            bool(cv_gate),
        "max_min_le_1p5":
            bool(ratio_gate),
        "bulk_500ps_isotropy_gate":
            bool(isotropy_gate),
    },
    "late_block_diagnostics": {
        "trace_cv_200_500":
            float(trace_cv_late),
        "component_block_cv":
            component_block_cv,
    },
}

(OUT/"F46T_summary.json").write_text(
    json.dumps(
        summary,
        indent=2,
    )
)

# ======================================================================
# PRINT
# ======================================================================

print()
print("="*84)
print("F46T BULK CUMULATIVE ISOTROPY")
print("="*84)

print(
    cum[
        [
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
print("F46T INDEPENDENT 100 ps BLOCKS")
print("="*84)

print(
    block[
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
print("F46T LATE WINDOWS")
print("="*84)

print(
    late[
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
print("F46T COLLECTIVE CORRELATION WINDOWS")
print("="*84)

print(
    corr.to_string(
        index=False
    )
)

print()
print("="*84)
print("F46T FINAL BULK 500 ps GATE")
print("="*84)

print(
    f"BULK_500PS_XX="
    f"{final.XX}"
)

print(
    f"BULK_500PS_YY="
    f"{final.YY}"
)

print(
    f"BULK_500PS_ZZ="
    f"{final.ZZ}"
)

print(
    f"BULK_500PS_TRACE="
    f"{final['trace']}"
)

print(
    f"BULK_500PS_DIAG_CV="
    f"{final.diag_cv}"
)

print(
    "BULK_500PS_DIAG_MAX_MIN_RATIO="
    f"{final.diag_max_min_ratio}"
)

print(
    "BULK_500PS_CV_GATE="
    + (
        "PASS"
        if cv_gate
        else "FAIL"
    )
)

print(
    "BULK_500PS_MAXMIN_GATE="
    + (
        "PASS"
        if ratio_gate
        else "FAIL"
    )
)

print(
    "BULK_500PS_ISOTROPY_GATE="
    + (
        "PASS"
        if isotropy_gate
        else "FAIL"
    )
)

print()
print(
    "LATE_200_500_TRACE_CV="
    f"{trace_cv_late}"
)

for comp,val in (
    component_block_cv.items()
):
    print(
        f"LATE_100PS_BLOCK_{comp}_CV="
        f"{val}"
    )

print()
print(
    "F46T_DECISION="
    + (
        "BULK_CONVERGED_PROCEED_TO_ENGINEERED_MATCHED_DURATION"
        if isotropy_gate
        else
        "BULK_NOT_YET_CONVERGED_ADJUDICATE_SAMPLING_STRATEGY"
    )
)

print("F46T_ANALYSIS_COMPLETE=YES")
