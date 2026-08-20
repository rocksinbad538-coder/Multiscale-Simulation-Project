#!/usr/bin/env python3

from pathlib import Path
import json
import numpy as np
import pandas as pd

from md_analysis.statistics import (
    integrated_autocorrelation_time,
    effective_sample_size,
    statistical_inefficiency,
)

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
)

TEMPERATURES = [150, 200, 250, 300, 350]

DT_FS = 0.25
OUTPUT_EVERY_STEPS = 100
FRAME_INTERVAL_PS = DT_FS * OUTPUT_EVERY_STEPS / 1000.0

MIN_REMAINING_FRACTION = 0.25
N_T0_CANDIDATES = 200
MIN_NEFF = 50.0
N_BLOCKS = 5

MAX_BLOCK_DRIFT = {
    "Rg": 0.02,
    "AlignedRMSD": 0.05,
    "PotentialEnergy": 0.02,
}


def relative_change(a, b):
    scale = max(abs(a), abs(b), 1.0e-12)
    return abs(a-b) / scale


def block_means(x, nblocks=N_BLOCKS):
    x = np.asarray(x, dtype=float)

    edges = np.linspace(
        0,
        len(x),
        nblocks + 1,
        dtype=int,
    )

    means = []

    for i in range(nblocks):
        part = x[edges[i]:edges[i+1]]

        if len(part):
            means.append(float(np.mean(part)))

    return means


def max_adjacent_block_drift(means):
    if len(means) < 2:
        return float("inf")

    return max(
        relative_change(means[i], means[i+1])
        for i in range(len(means)-1)
    )


def detect_equilibration(values):
    """
    Scan candidate equilibration indices t0.

    For each suffix x[t0:], estimate statistical inefficiency g and

        N_eff = N_remaining / g.

    Select t0 maximizing N_eff.

    Search is restricted so at least MIN_REMAINING_FRACTION of the
    trajectory remains.
    """

    x = np.asarray(values, dtype=float)
    n = len(x)

    if n < 100:
        raise RuntimeError(
            f"Insufficient samples for equilibration detection: {n}"
        )

    max_t0 = int(
        n * (1.0 - MIN_REMAINING_FRACTION)
    )

    candidates = np.unique(
        np.linspace(
            0,
            max_t0,
            min(N_T0_CANDIDATES, max_t0 + 1),
            dtype=int,
        )
    )

    best = None

    for t0 in candidates:

        tail = x[t0:]

        if len(tail) < 50:
            continue

        g = statistical_inefficiency(tail)
        neff = len(tail) / g

        record = {
            "t0_index": int(t0),
            "n_remaining": int(len(tail)),
            "g": float(g),
            "N_eff": float(neff),
        }

        if best is None or record["N_eff"] > best["N_eff"]:
            best = record

    if best is None:
        raise RuntimeError(
            "No valid equilibration candidate."
        )

    return best


def analyze_equilibrated_series(name, values, start):
    x = np.asarray(values, dtype=float)[start:]

    tau = integrated_autocorrelation_time(x)
    g = statistical_inefficiency(x)
    neff = effective_sample_size(x)

    means = block_means(x)
    drift = max_adjacent_block_drift(means)

    spacing_frames = int(
        np.ceil(max(1.0, g))
    )

    stable = (
        neff >= MIN_NEFF
        and drift <= MAX_BLOCK_DRIFT[name]
    )

    return {
        "observable": name,
        "n_samples": int(len(x)),
        "tau_int_samples": float(tau),
        "tau_int_ps": float(tau * FRAME_INTERVAL_PS),
        "statistical_inefficiency": float(g),
        "N_eff": float(neff),
        "block_means": means,
        "max_adjacent_block_mean_drift": float(drift),
        "maximum_allowed_block_drift":
            float(MAX_BLOCK_DRIFT[name]),
        "recommended_spacing_frames":
            int(spacing_frames),
        "recommended_spacing_ps":
            float(spacing_frames * FRAME_INTERVAL_PS),
        "stable": bool(stable),
    }


summary = []

for T in TEMPERATURES:

    analysis = CAMPAIGN / f"{T}K" / "analysis"

    thermo_file = analysis / "thermodynamics.csv"
    traj_file = analysis / "trajectory_summary.csv"
    rmsd_file = analysis / "aligned_rmsd.csv"

    required = [thermo_file, traj_file, rmsd_file]

    missing = [
        p.name for p in required if not p.exists()
    ]

    if missing:
        print(f"{T}K: SKIP — missing {missing}")
        continue

    thermo = pd.read_csv(thermo_file)
    traj = pd.read_csv(traj_file)
    rmsd = pd.read_csv(rmsd_file)

    n = min(
        len(thermo),
        len(traj),
        len(rmsd),
    )

    if n < 100:
        print(f"{T}K: SKIP — only {n} samples")
        continue

    step_traj = traj["timestep"].iloc[:n].to_numpy(dtype=int)
    step_rmsd = rmsd["timestep"].iloc[:n].to_numpy(dtype=int)
    step_thermo = thermo["Step"].iloc[:n].to_numpy(dtype=int)

    if not (
        np.array_equal(step_traj, step_rmsd)
        and np.array_equal(step_traj, step_thermo)
    ):
        raise RuntimeError(
            f"{T}K: timestep grids are not identical "
            "between trajectory, RMSD and thermodynamics."
        )

    rg = traj["Rg"].iloc[:n].to_numpy(dtype=float)
    ar = rmsd["AlignedRMSD"].iloc[:n].to_numpy(dtype=float)
    pe = thermo["PotEng"].iloc[:n].to_numpy(dtype=float)

    equilibration = {
        "Rg": detect_equilibration(rg),
        "AlignedRMSD": detect_equilibration(ar),
        "PotentialEnergy": detect_equilibration(pe),
    }

    # Conservative common equilibration point.
    global_t0 = max(
        item["t0_index"]
        for item in equilibration.values()
    )

    global_t0_step = int(
        step_traj[global_t0]
    )

    global_t0_ns = (
        global_t0_step
        * DT_FS
        / 1.0e6
    )

    observables = [
        analyze_equilibrated_series(
            "Rg", rg, global_t0
        ),
        analyze_equilibrated_series(
            "AlignedRMSD", ar, global_t0
        ),
        analyze_equilibrated_series(
            "PotentialEnergy", pe, global_t0
        ),
    ]

    spacing_frames = max(
        x["recommended_spacing_frames"]
        for x in observables
    )

    spacing_steps = (
        spacing_frames
        * OUTPUT_EVERY_STEPS
    )

    spacing_ps = (
        spacing_steps
        * DT_FS
        / 1000.0
    )

    remaining_frames = n-global_t0

    max_independent_snapshots = (
        1 + (remaining_frames-1)//spacing_frames
        if remaining_frames > 0
        else 0
    )

    min_neff = min(
        x["N_eff"]
        for x in observables
    )

    stable = all(
        x["stable"]
        for x in observables
    )

    record = {
        "temperature_K": T,
        "n_total_samples": int(n),
        "equilibration_by_observable":
            equilibration,
        "global_t0_index":
            int(global_t0),
        "global_t0_step":
            int(global_t0_step),
        "global_t0_ns":
            float(global_t0_ns),
        "equilibrated_frame_count":
            int(remaining_frames),
        "minimum_N_eff":
            float(min_neff),
        "recommended_spacing_frames":
            int(spacing_frames),
        "recommended_spacing_steps":
            int(spacing_steps),
        "recommended_spacing_ps":
            float(spacing_ps),
        "maximum_independent_snapshots":
            int(max_independent_snapshots),
        "observables":
            observables,
        "simulation_stable":
            bool(stable),
    }

    summary.append(record)

    print()
    print("="*72)
    print(f"{T} K")
    print("="*72)

    for name, item in equilibration.items():
        print(
            f"{name:20s} "
            f"t0_index={item['t0_index']:8d} "
            f"Neff_at_optimum={item['N_eff']:.2f}"
        )

    print(
        f"GLOBAL t0 = step {global_t0_step} "
        f"({global_t0_ns:.6f} ns)"
    )

    for obs in observables:
        print(
            f"{obs['observable']:20s} "
            f"tau={obs['tau_int_samples']:.3f} frames  "
            f"g={obs['statistical_inefficiency']:.3f}  "
            f"Neff={obs['N_eff']:.2f}  "
            f"drift={obs['max_adjacent_block_mean_drift']:.5f}  "
            f"{'PASS' if obs['stable'] else 'FAIL'}"
        )

    print(
        "Independent snapshot spacing = "
        f"{spacing_frames} frames / "
        f"{spacing_steps} steps / "
        f"{spacing_ps:.6f} ps"
    )

    print(
        "Maximum independent snapshots = "
        f"{max_independent_snapshots}"
    )

    print(
        "TEMPERATURE_GATE="
        f"{'PASS' if stable else 'FAIL'}"
    )


OUT = CAMPAIGN / "ensemble_convergence"
OUT.mkdir(parents=True, exist_ok=True)

json_file = OUT / "ensemble_convergence_gate.json"
csv_file = OUT / "ensemble_convergence_gate.csv"
status_file = OUT / "ENSEMBLE_CONVERGENCE_STATUS.txt"

json_file.write_text(
    json.dumps(summary, indent=2) + "\n"
)

rows = []

for record in summary:

    for obs in record["observables"]:

        rows.append({
            "Temperature_K":
                record["temperature_K"],
            "Global_t0_index":
                record["global_t0_index"],
            "Global_t0_step":
                record["global_t0_step"],
            "Global_t0_ns":
                record["global_t0_ns"],
            "Observable":
                obs["observable"],
            "TauInt_samples":
                obs["tau_int_samples"],
            "TauInt_ps":
                obs["tau_int_ps"],
            "StatisticalInefficiency":
                obs["statistical_inefficiency"],
            "N_eff":
                obs["N_eff"],
            "BlockDrift":
                obs[
                    "max_adjacent_block_mean_drift"
                ],
            "Stable":
                obs["stable"],
            "Spacing_frames":
                record["recommended_spacing_frames"],
            "Spacing_steps":
                record["recommended_spacing_steps"],
            "Spacing_ps":
                record["recommended_spacing_ps"],
            "MaxIndependentSnapshots":
                record["maximum_independent_snapshots"],
        })

pd.DataFrame(rows).to_csv(
    csv_file,
    index=False
)

all_stable = (
    len(summary) == len(TEMPERATURES)
    and all(
        r["simulation_stable"]
        for r in summary
    )
)

status = (
    "PHASE5_ENSEMBLE_CONVERGENCE_PASS"
    if all_stable
    else
    "PHASE5_ENSEMBLE_CONVERGENCE_NOT_YET_PASS"
)

status_file.write_text(status + "\n")

print()
print("="*90)
print("PHASE5-D45")
print("AUTOMATED EQUILIBRATION + CONVERGENCE GATE")
print("="*90)
print(f"RESULTS={json_file}")
print(f"STATUS={status}")
print("="*90)
