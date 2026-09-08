#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import sys

import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import minimize_vectors


# ======================================================================
# CONSTANTS
# ======================================================================

QH = +0.5564
QM = -1.1128

E_CHARGE = 1.602176634e-19
NM_TO_M = 1.0e-9
DEBYE_CM = 3.33564e-30

EPS0 = 8.8541878128e-12
KB = 1.380649e-23
TEMP_K = 300.0


# ======================================================================
# UTILITIES
# ======================================================================

def discover_extension(directory, prefix):

    d = Path(directory)

    candidates = sorted(
        list(d.glob(f"{prefix}.xtc"))
        + list(d.glob(f"{prefix}.part*.xtc"))
    )

    if not candidates:
        return None

    # When -noappend was used, the relevant continuation is normally
    # the part-numbered XTC. Choose the largest candidate defensively.
    return max(
        candidates,
        key=lambda p: p.stat().st_size,
    )


def load_segment(tpr, xtc):

    u = mda.Universe(
        str(tpr),
        str(xtc),
    )

    times = []
    frames = []

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

        if not required.issubset(amap):
            raise RuntimeError(
                f"Unexpected TIP4P/2005 water definition "
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

    for ts in u.trajectory:

        p_nm = (
            u.atoms.positions
            * 0.1
        )

        box_nm = ts.dimensions.copy()
        box_nm[:3] *= 0.1

        rO = p_nm[ow]

        OH1 = minimize_vectors(
            p_nm[h1] - rO,
            box_nm,
        )

        OH2 = minimize_vectors(
            p_nm[h2] - rO,
            box_nm,
        )

        OM = minimize_vectors(
            p_nm[mw] - rO,
            box_nm,
        )

        mu = (
            QH * OH1
            + QH * OH2
            + QM * OM
        )

        M = np.sum(
            mu,
            axis=0,
            dtype=np.float64,
        )

        times.append(
            float(ts.time)
        )

        frames.append({
            "M_e_nm": M,
            "box_nm": box_nm[:3].copy(),
        })

    return (
        np.asarray(times),
        frames,
        len(waters),
    )


def combine_segments(
    original_tpr,
    original_xtc,
    extension_tpr,
    extension_xtc,
):

    t0, f0, nw0 = load_segment(
        original_tpr,
        original_xtc,
    )

    t1, f1, nw1 = load_segment(
        extension_tpr,
        extension_xtc,
    )

    if nw0 != nw1:
        raise RuntimeError(
            "Water count differs between segments."
        )

    all_t = np.concatenate(
        [t0, t1]
    )

    all_f = f0 + f1

    # Deduplicate by absolute trajectory time.
    # The checkpoint continuation may repeat 100.0 ps.
    keep = []

    seen = set()

    for i,t in enumerate(all_t):

        key = round(
            float(t),
            6,
        )

        if key not in seen:
            seen.add(key)
            keep.append(i)

    time = all_t[keep]

    M = np.asarray([
        all_f[i]["M_e_nm"]
        for i in keep
    ])

    boxes = np.asarray([
        all_f[i]["box_nm"]
        for i in keep
    ])

    order = np.argsort(time)

    time = time[order]
    M = M[order]
    boxes = boxes[order]

    return (
        time,
        M,
        boxes,
        nw0,
    )


# ======================================================================
# FLUCTUATION TENSOR
# ======================================================================

def fluctuation_tensor(
    M_e_nm,
    boxes_nm,
):

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

    V = (
        volume_nm3
        * 1e-27
    )

    tensor = (
        cov
        / (
            EPS0
            * V
            * KB
            * TEMP_K
        )
    )

    return (
        tensor,
        volume_nm3,
    )


def tensor_metrics(T):

    diag = np.diag(T)

    mean = float(
        np.mean(diag)
    )

    sd = float(
        np.std(
            diag,
            ddof=1,
        )
    )

    cv = (
        sd / mean
        if mean != 0
        else np.nan
    )

    maxmin = (
        float(
            np.max(diag)
            / np.min(diag)
        )
        if np.min(diag) > 0
        else np.nan
    )

    return {
        "XX": float(T[0,0]),
        "YY": float(T[1,1]),
        "ZZ": float(T[2,2]),
        "XY": float(T[0,1]),
        "XZ": float(T[0,2]),
        "YZ": float(T[1,2]),
        "trace": float(np.trace(T)),
        "diag_cv": cv,
        "diag_max_min_ratio": maxmin,
    }


# ======================================================================
# COLLECTIVE CORRELATION
# ======================================================================

def connected_corr(M):

    X = (
        M
        - np.mean(
            M,
            axis=0,
            keepdims=True,
        )
    )

    denom = np.mean(
        np.sum(
            X*X,
            axis=1,
        )
    )

    nf = len(X)

    C = np.empty(nf)
    origins = np.empty(
        nf,
        dtype=int,
    )

    for lag in range(nf):

        dots = np.sum(
            X[:nf-lag]
            * X[lag:],
            axis=1,
        )

        C[lag] = (
            np.mean(dots)
            / denom
        )

        origins[lag] = (
            nf-lag
        )

    return C, origins


def first_zero_integral(
    time,
    corr,
):

    idx = np.where(
        corr <= 0
    )[0]

    if len(idx) == 0:
        stop = len(corr)
        zero = None
    else:
        stop = int(idx[0]) + 1
        zero = float(
            time[idx[0]]
        )

    val = float(
        np.trapezoid(
            corr[:stop],
            time[:stop],
        )
    )

    return val, zero


def integral_to(
    time,
    corr,
    cutoff,
):

    m = (
        time <= cutoff
    )

    return float(
        np.trapezoid(
            corr[m],
            time[m],
        )
    )


# ======================================================================
# WINDOW ANALYSIS
# ======================================================================

def analyze_windows(
    system,
    time,
    M,
    boxes,
):

    rows = []

    # Cumulative convergence
    for end in range(
        25,
        201,
        25,
    ):

        mask = (
            (time >= 0)
            & (time <= end)
        )

        T, V = fluctuation_tensor(
            M[mask],
            boxes[mask],
        )

        met = tensor_metrics(T)

        rows.append({
            "system": system,
            "window_type": "cumulative",
            "start_ps": 0,
            "end_ps": end,
            "n_frames": int(
                np.sum(mask)
            ),
            "volume_nm3": V,
            **met,
        })

    # Independent 50 ps blocks
    for start in [
        0,50,100,150
    ]:

        end = start + 50

        mask = (
            (time >= start)
            & (time <= end)
        )

        T,V = fluctuation_tensor(
            M[mask],
            boxes[mask],
        )

        met = tensor_metrics(T)

        rows.append({
            "system": system,
            "window_type": "block50",
            "start_ps": start,
            "end_ps": end,
            "n_frames": int(
                np.sum(mask)
            ),
            "volume_nm3": V,
            **met,
        })

    # Independent 100 ps blocks
    for start in [
        0,100
    ]:

        end = start + 100

        mask = (
            (time >= start)
            & (time <= end)
        )

        T,V = fluctuation_tensor(
            M[mask],
            boxes[mask],
        )

        met = tensor_metrics(T)

        rows.append({
            "system": system,
            "window_type": "block100",
            "start_ps": start,
            "end_ps": end,
            "n_frames": int(
                np.sum(mask)
            ),
            "volume_nm3": V,
            **met,
        })

    return pd.DataFrame(rows)


def correlation_windows(
    system,
    time,
    M,
):

    rows = []

    windows = [
        (0,200),
        (50,200),
        (100,200),
        (150,200),
    ]

    for start,end in windows:

        mask = (
            (time >= start)
            & (time <= end)
        )

        t = (
            time[mask]
            - time[mask][0]
        )

        X = M[mask]

        C, origins = (
            connected_corr(X)
        )

        iz,zt = (
            first_zero_integral(
                t,C
            )
        )

        def at(x):
            return float(
                C[
                    np.argmin(
                        np.abs(t-x)
                    )
                ]
            )

        rows.append({
            "system": system,
            "start_ps": start,
            "end_ps": end,
            "duration_ps": end-start,
            "C_2ps": at(2),
            "C_5ps": at(5),
            "C_10ps": at(10),
            "integral_5ps":
                integral_to(
                    t,C,5
                ),
            "integral_10ps":
                integral_to(
                    t,C,10
                ),
            "first_zero_integral_ps":
                iz,
            "first_zero_time_ps":
                zt,
        })

    return pd.DataFrame(rows)


# ======================================================================
# MAIN
# ======================================================================

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        required=True,
    )

    ap.add_argument(
        "--preflight-only",
        action="store_true",
    )

    args = ap.parse_args()

    root = Path(args.root)
    out = root / "F46R_collective_200ps"

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    phase_root = Path(
        "runs/phase1A/accepted/"
        "hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute"
    )

    bulk_root = Path(
        "runs/phase2/day048_f46l_bulk_100ps_orientation"
    )

    eng_ext_dir = (
        root/"engineered"
    )

    bulk_ext_dir = (
        root/"bulk"
    )

    eng_ext = discover_extension(
        eng_ext_dir,
        "engineered_extend_100ps",
    )

    bulk_ext = discover_extension(
        bulk_ext_dir,
        "bulk_extend_100ps",
    )

    print("="*76)
    print("F46R — 0-200 ps COLLECTIVE POLARIZATION ANALYSIS")
    print("="*76)

    print(
        f"ENGINEERED_EXTENSION_XTC={eng_ext}"
    )

    print(
        f"BULK_EXTENSION_XTC={bulk_ext}"
    )

    if args.preflight_only:

        print(
            "ENGINEERED_EXTENSION_AVAILABLE="
            + (
                "YES"
                if eng_ext is not None
                else "NO"
            )
        )

        print(
            "BULK_EXTENSION_AVAILABLE="
            + (
                "YES"
                if bulk_ext is not None
                else "NO"
            )
        )

        print(
            "F46R_PREFLIGHT_COMPLETE=YES"
        )

        return

    if eng_ext is None:
        raise RuntimeError(
            "Engineered extension XTC missing."
        )

    if bulk_ext is None:
        raise RuntimeError(
            "Bulk extension XTC missing."
        )

    systems = {
        "engineered": {
            "orig_tpr":
                phase_root/
                "nvt_100ps_frozenSolute.tpr",

            "orig_xtc":
                phase_root/
                "nvt_100ps_frozenSolute.xtc",

            "ext_tpr":
                eng_ext_dir/
                "engineered_200ps.tpr",

            "ext_xtc":
                eng_ext,
        },

        "bulk": {
            "orig_tpr":
                bulk_root/
                "nvt_bulk_100ps.tpr",

            "orig_xtc":
                bulk_root/
                "nvt_bulk_100ps.xtc",

            "ext_tpr":
                bulk_ext_dir/
                "bulk_200ps.tpr",

            "ext_xtc":
                bulk_ext,
        },
    }

    tensors = []
    correlations = []

    full = {}

    for label,dat in systems.items():

        print()
        print(
            f"LOADING_{label.upper()}=YES"
        )

        time,M,boxes,nwater = (
            combine_segments(
                dat["orig_tpr"],
                dat["orig_xtc"],
                dat["ext_tpr"],
                dat["ext_xtc"],
            )
        )

        print(
            f"{label.upper()}_N_WATER="
            f"{nwater}"
        )

        print(
            f"{label.upper()}_N_FRAMES="
            f"{len(time)}"
        )

        print(
            f"{label.upper()}_FIRST_PS="
            f"{time[0]}"
        )

        print(
            f"{label.upper()}_LAST_PS="
            f"{time[-1]}"
        )

        dt = np.diff(time)

        complete = (
            len(time) == 401
            and abs(time[0]) < 1e-4
            and abs(time[-1]-200) < 1e-4
            and np.allclose(
                dt,
                0.5,
                atol=1e-4,
                rtol=0,
            )
        )

        print(
            f"{label.upper()}_TRAJECTORY_GATE="
            + (
                "PASS"
                if complete
                else "FAIL"
            )
        )

        if not complete:
            raise RuntimeError(
                f"{label}: incomplete 0-200 ps trajectory"
            )

        T,V = fluctuation_tensor(
            M,
            boxes,
        )

        met = tensor_metrics(T)

        C,origins = connected_corr(M)

        iz,zt = first_zero_integral(
            time,
            C,
        )

        full[label] = {
            "tensor": met,
            "volume_nm3": V,
            "first_zero_integral_ps":
                iz,
            "first_zero_time_ps":
                zt,
            "C_2ps":
                float(C[4]),
            "C_5ps":
                float(C[10]),
            "C_10ps":
                float(C[20]),
        }

        tw = analyze_windows(
            label,
            time,
            M,
            boxes,
        )

        cw = correlation_windows(
            label,
            time,
            M,
        )

        tensors.append(tw)
        correlations.append(cw)

        pd.DataFrame({
            "time_ps": time,
            "C_M_connected": C,
            "n_origins": origins,
        }).to_csv(
            out/
            f"{label}_collective_C_M_200ps.csv",
            index=False,
        )

    tdf = pd.concat(
        tensors,
        ignore_index=True,
    )

    cdf = pd.concat(
        correlations,
        ignore_index=True,
    )

    tdf.to_csv(
        out/"F46R_tensor_convergence.csv",
        index=False,
    )

    cdf.to_csv(
        out/"F46R_collective_correlation_windows.csv",
        index=False,
    )

    print()
    print("="*76)
    print("F46R FULL 0-200 ps RESULTS")
    print("="*76)

    for label in [
        "engineered",
        "bulk",
    ]:

        r = full[label]

        print()
        print(label.upper())

        for k,v in (
            r["tensor"].items()
        ):
            print(
                f"{k.upper()}={v}"
            )

        print(
            "FIRST_ZERO_INTEGRAL_PS="
            f"{r['first_zero_integral_ps']}"
        )

        print(
            "FIRST_ZERO_TIME_PS="
            f"{r['first_zero_time_ps']}"
        )

        print(
            f"C_M_2PS={r['C_2ps']}"
        )

        print(
            f"C_M_5PS={r['C_5ps']}"
        )

        print(
            f"C_M_10PS={r['C_10ps']}"
        )

    print()
    print("="*76)
    print("F46R BULK CUMULATIVE ISOTROPY")
    print("="*76)

    bulk_cum = tdf[
        (tdf.system=="bulk")
        & (
            tdf.window_type
            =="cumulative"
        )
    ].sort_values(
        "end_ps"
    )

    print(
        bulk_cum[
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
    print("="*76)
    print("F46R INDEPENDENT 100 ps BLOCKS")
    print("="*76)

    print(
        tdf[
            tdf.window_type
            =="block100"
        ][
            [
                "system",
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
    print("="*76)
    print("F46R COLLECTIVE CORRELATION STATIONARITY")
    print("="*76)

    print(
        cdf.to_string(
            index=False
        )
    )

    E = full["engineered"]["tensor"]
    B = full["bulk"]["tensor"]

    trace_ratio = (
        E["trace"]
        / B["trace"]
    )

    print()
    print("="*76)
    print("F46R ENGINEERED / BULK")
    print("="*76)

    print(
        "TRACE_RATIO_ENGINEERED_OVER_BULK="
        f"{trace_ratio:.12f}"
    )

    print(
        "XX_RATIO_ENGINEERED_OVER_BULK="
        f"{E['XX']/B['XX']:.12f}"
    )

    print(
        "YY_RATIO_ENGINEERED_OVER_BULK="
        f"{E['YY']/B['YY']:.12f}"
    )

    print(
        "ZZ_RATIO_ENGINEERED_OVER_BULK="
        f"{E['ZZ']/B['ZZ']:.12f}"
    )

    bulk_iso = (
        B["diag_max_min_ratio"] <= 1.5
        and B["diag_cv"] <= 0.20
    )

    print(
        "BULK_200PS_ISOTROPY_GATE="
        + (
            "PASS"
            if bulk_iso
            else "FAIL"
        )
    )

    if bulk_iso:
        decision = (
            "M3P5_LOW_FREQUENCY_RESPONSE_"
            "CAN_PROCEED_TO_NEXT_VALIDATION"
        )
    else:
        decision = (
            "M3P5_COLLECTIVE_SAMPLING_"
            "STILL_NOT_CONVERGED"
        )

    print(
        f"F46R_DECISION={decision}"
    )

    summary = {
        "full_results": full,
        "bulk_isotropy_gate":
            bulk_iso,
        "trace_ratio_engineered_over_bulk":
            trace_ratio,
        "decision":
            decision,
        "interpretation_boundary": [
            (
                "No quantum coherence lifetime "
                "is inferred from this analysis."
            ),
            (
                "The engineered frozen-solute trajectory "
                "contains a known large static solute "
                "LJ/virial contribution."
            ),
            (
                "Water-solute and water-water geometry "
                "showed no pathological short contacts."
            ),
            (
                "A dielectric-response interpretation "
                "requires collective convergence."
            ),
        ],
    }

    with (
        out/"F46R_summary.json"
    ).open("w") as fh:

        json.dump(
            summary,
            fh,
            indent=2,
        )

    print()
    print(
        "F46R_ANALYSIS_COMPLETE=YES"
    )


if __name__ == "__main__":
    main()
