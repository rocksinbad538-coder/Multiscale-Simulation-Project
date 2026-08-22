#!/usr/bin/env python3

from pathlib import Path
import json
import time

import numpy as np
import pandas as pd

from md_analysis.statistics import (
    integrated_autocorrelation_time,
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

FRAME_INTERVAL_PS = (
    DT_FS
    * OUTPUT_EVERY_STEPS
    / 1000.0
)

MIN_REMAINING_FRACTION = 0.25

# Full-resolution two-stage t0 search.
N_COARSE_CANDIDATES = 41
N_REFINE_CANDIDATES = 21

MIN_NEFF = 50.0
N_BLOCKS = 5

# Operational stability criteria.
# These are explicitly reported as analysis thresholds,
# not treated as universal physical constants.
MAX_BLOCK_DRIFT = {
    "Rg": 0.02,
    "AlignedRMSD": 0.05,
    "PotentialEnergy": 0.02,
}


def correlation_statistics(x):
    """
    Compute autocorrelation-derived quantities once.

    tau_int is obtained with FFT ACF + Geyer IPS.
    """
    x = np.asarray(
        x,
        dtype=float,
    )

    tau = (
        integrated_autocorrelation_time(
            x
        )
    )

    g = max(
        1.0,
        2.0 * tau,
    )

    neff = (
        len(x) / g
    )

    return {
        "tau_int_samples":
            float(tau),
        "statistical_inefficiency":
            float(g),
        "N_eff":
            float(neff),
    }


def evaluate_t0(x, t0):

    tail = x[t0:]

    if len(tail) < 50:
        return None

    stats = correlation_statistics(
        tail
    )

    return {
        "t0_index":
            int(t0),
        "n_remaining":
            int(len(tail)),
        **stats,
    }


def detect_equilibration(values):
    """
    Two-stage, full-resolution search for t0.

    Stage 1:
        coarse scan from frame 0 to the point where at least
        MIN_REMAINING_FRACTION of the trajectory remains.

    Stage 2:
        refine around the best coarse candidate.

    Objective:
        maximize N_eff(t0).

    No trajectory decimation is used.
    """

    x = np.asarray(
        values,
        dtype=float,
    )

    n = len(x)

    if n < 100:
        raise RuntimeError(
            f"Insufficient samples: {n}"
        )

    max_t0 = int(
        n
        * (
            1.0
            - MIN_REMAINING_FRACTION
        )
    )

    coarse = np.unique(
        np.linspace(
            0,
            max_t0,
            N_COARSE_CANDIDATES,
            dtype=int,
        )
    )

    coarse_results = []

    for t0 in coarse:

        result = evaluate_t0(
            x,
            int(t0),
        )

        if result is not None:
            coarse_results.append(
                result
            )

    if not coarse_results:
        raise RuntimeError(
            "No valid coarse t0 candidates."
        )

    best_coarse = max(
        coarse_results,
        key=lambda r: r["N_eff"],
    )

    coarse_step = max(
        1,
        int(
            np.ceil(
                max_t0
                / (
                    N_COARSE_CANDIDATES
                    - 1
                )
            )
        ),
    )

    refine_lo = max(
        0,
        best_coarse["t0_index"]
        - coarse_step,
    )

    refine_hi = min(
        max_t0,
        best_coarse["t0_index"]
        + coarse_step,
    )

    refine = np.unique(
        np.linspace(
            refine_lo,
            refine_hi,
            N_REFINE_CANDIDATES,
            dtype=int,
        )
    )

    refine_results = []

    for t0 in refine:

        result = evaluate_t0(
            x,
            int(t0),
        )

        if result is not None:
            refine_results.append(
                result
            )

    candidates = (
        coarse_results
        + refine_results
    )

    best = max(
        candidates,
        key=lambda r: r["N_eff"],
    )

    best = dict(best)

    best[
        "search_method"
    ] = (
        "TWO_STAGE_FULL_RESOLUTION_GRID"
    )

    best[
        "coarse_candidates"
    ] = len(coarse)

    best[
        "refine_candidates"
    ] = len(refine)

    return best


def relative_change(a, b):

    scale = max(
        abs(a),
        abs(b),
        1.0e-12,
    )

    return (
        abs(a-b)
        / scale
    )


def block_means(
    x,
    nblocks=N_BLOCKS,
):

    x = np.asarray(
        x,
        dtype=float,
    )

    edges = np.linspace(
        0,
        len(x),
        nblocks + 1,
        dtype=int,
    )

    means = []

    for i in range(nblocks):

        part = x[
            edges[i]:
            edges[i+1]
        ]

        if len(part):

            means.append(
                float(
                    np.mean(part)
                )
            )

    return means


def max_adjacent_block_drift(
    means
):

    if len(means) < 2:
        return float("inf")

    return max(
        relative_change(
            means[i],
            means[i+1],
        )
        for i in range(
            len(means)-1
        )
    )


def analyze_final_series(
    name,
    values,
    start,
):

    x = np.asarray(
        values,
        dtype=float,
    )[start:]

    stats = correlation_statistics(
        x
    )

    means = block_means(
        x
    )

    drift = (
        max_adjacent_block_drift(
            means
        )
    )

    g = stats[
        "statistical_inefficiency"
    ]

    spacing_frames = int(
        np.ceil(g)
    )

    stable = (
        stats["N_eff"]
        >= MIN_NEFF
        and
        drift
        <= MAX_BLOCK_DRIFT[name]
    )

    return {
        "observable":
            name,
        "n_samples":
            int(len(x)),
        "tau_int_samples":
            stats[
                "tau_int_samples"
            ],
        "tau_int_ps":
            float(
                stats[
                    "tau_int_samples"
                ]
                * FRAME_INTERVAL_PS
            ),
        "statistical_inefficiency":
            float(g),
        "N_eff":
            stats["N_eff"],
        "block_means":
            means,
        "max_adjacent_block_mean_drift":
            float(drift),
        "maximum_allowed_block_drift":
            float(
                MAX_BLOCK_DRIFT[name]
            ),
        "recommended_spacing_frames":
            int(
                spacing_frames
            ),
        "recommended_spacing_ps":
            float(
                spacing_frames
                * FRAME_INTERVAL_PS
            ),
        "stable":
            bool(stable),
    }


summary = []

campaign_start = time.time()


for T in TEMPERATURES:

    temperature_start = (
        time.time()
    )

    analysis = (
        CAMPAIGN
        / f"{T}K"
        / "analysis"
    )

    thermo_file = (
        analysis
        / "thermodynamics.csv"
    )

    traj_file = (
        analysis
        / "trajectory_summary.csv"
    )

    rmsd_file = (
        analysis
        / "aligned_rmsd.csv"
    )

    required = [
        thermo_file,
        traj_file,
        rmsd_file,
    ]

    missing = [
        p.name
        for p in required
        if not p.exists()
    ]

    if missing:
        raise RuntimeError(
            f"{T}K missing files: "
            f"{missing}"
        )

    thermo = pd.read_csv(
        thermo_file
    )

    traj = pd.read_csv(
        traj_file
    )

    rmsd = pd.read_csv(
        rmsd_file
    )

    n = min(
        len(thermo),
        len(traj),
        len(rmsd),
    )

    if n != 800001:
        raise RuntimeError(
            f"{T}K expected 800001 "
            f"samples, found {n}."
        )

    step_traj = (
        traj["timestep"]
        .iloc[:n]
        .to_numpy(dtype=int)
    )

    step_rmsd = (
        rmsd["timestep"]
        .iloc[:n]
        .to_numpy(dtype=int)
    )

    step_thermo = (
        thermo["Step"]
        .iloc[:n]
        .to_numpy(dtype=int)
    )

    if not (
        np.array_equal(
            step_traj,
            step_rmsd,
        )
        and
        np.array_equal(
            step_traj,
            step_thermo,
        )
    ):

        raise RuntimeError(
            f"{T}K timestep grids "
            "are inconsistent."
        )

    rg = (
        traj["Rg"]
        .iloc[:n]
        .to_numpy(dtype=float)
    )

    aligned = (
        rmsd["AlignedRMSD"]
        .iloc[:n]
        .to_numpy(dtype=float)
    )

    pe = (
        thermo["PotEng"]
        .iloc[:n]
        .to_numpy(dtype=float)
    )

    if not (
        np.all(np.isfinite(rg))
        and
        np.all(np.isfinite(aligned))
        and
        np.all(np.isfinite(pe))
    ):

        raise RuntimeError(
            f"{T}K contains "
            "non-finite observables."
        )

    print()
    print(
        "="*72
    )
    print(
        f"{T} K — EQUILIBRATION SEARCH"
    )
    print(
        "="*72
    )

    equilibration = {}

    for name, series in (
        ("Rg", rg),
        (
            "AlignedRMSD",
            aligned,
        ),
        (
            "PotentialEnergy",
            pe,
        ),
    ):

        start = time.time()

        result = (
            detect_equilibration(
                series
            )
        )

        elapsed = (
            time.time()-start
        )

        result[
            "search_elapsed_seconds"
        ] = float(elapsed)

        equilibration[name] = (
            result
        )

        print(
            f"{name:20s} "
            f"t0_index="
            f"{result['t0_index']:7d}  "
            f"Neff="
            f"{result['N_eff']:10.2f}  "
            f"time="
            f"{elapsed:7.2f}s"
        )

    # Conservative common equilibration point:
    # all observables must be past their individual t0.
    global_t0 = max(
        item["t0_index"]
        for item
        in equilibration.values()
    )

    global_t0_step = int(
        step_traj[
            global_t0
        ]
    )

    global_t0_ns = (
        global_t0_step
        * DT_FS
        / 1.0e6
    )

    observables = [
        analyze_final_series(
            "Rg",
            rg,
            global_t0,
        ),
        analyze_final_series(
            "AlignedRMSD",
            aligned,
            global_t0,
        ),
        analyze_final_series(
            "PotentialEnergy",
            pe,
            global_t0,
        ),
    ]

    spacing_frames = max(
        item[
            "recommended_spacing_frames"
        ]
        for item
        in observables
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

    remaining_frames = (
        n-global_t0
    )

    max_independent_snapshots = (
        1
        + (
            remaining_frames-1
        )
        // spacing_frames
    )

    minimum_neff = min(
        item["N_eff"]
        for item
        in observables
    )

    stable = all(
        item["stable"]
        for item
        in observables
    )

    elapsed_temperature = (
        time.time()
        - temperature_start
    )

    record = {
        "temperature_K":
            int(T),
        "n_total_samples":
            int(n),
        "equilibration_search_method":
            "TWO_STAGE_FULL_RESOLUTION_GRID",
        "equilibration_by_observable":
            equilibration,
        "global_t0_index":
            int(global_t0),
        "global_t0_step":
            int(global_t0_step),
        "global_t0_ns":
            float(global_t0_ns),
        "equilibrated_frame_count":
            int(
                remaining_frames
            ),
        "minimum_N_eff":
            float(
                minimum_neff
            ),
        "recommended_spacing_frames":
            int(
                spacing_frames
            ),
        "recommended_spacing_steps":
            int(
                spacing_steps
            ),
        "recommended_spacing_ps":
            float(
                spacing_ps
            ),
        "maximum_independent_snapshots":
            int(
                max_independent_snapshots
            ),
        "observables":
            observables,
        "simulation_stable":
            bool(stable),
        "elapsed_seconds":
            float(
                elapsed_temperature
            ),
    }

    summary.append(
        record
    )

    print()
    print(
        f"GLOBAL_t0_STEP="
        f"{global_t0_step}"
    )

    print(
        f"GLOBAL_t0_NS="
        f"{global_t0_ns:.6f}"
    )

    for obs in observables:

        print(
            f"{obs['observable']:20s} "
            f"tau="
            f"{obs['tau_int_samples']:9.3f} "
            f"frames  "
            f"g="
            f"{obs['statistical_inefficiency']:9.3f}  "
            f"Neff="
            f"{obs['N_eff']:10.2f}  "
            f"drift="
            f"{obs['max_adjacent_block_mean_drift']:.6f}  "
            f"{'PASS' if obs['stable'] else 'FAIL'}"
        )

    print(
        "INDEPENDENT_SPACING="
        f"{spacing_frames} frames / "
        f"{spacing_steps} steps / "
        f"{spacing_ps:.6f} ps"
    )

    print(
        "MAX_INDEPENDENT_SNAPSHOTS="
        f"{max_independent_snapshots}"
    )

    print(
        "TEMPERATURE_GATE="
        f"{'PASS' if stable else 'FAIL'}"
    )

    print(
        f"TEMPERATURE_ELAPSED_S="
        f"{elapsed_temperature:.2f}"
    )


OUT = (
    CAMPAIGN
    / "ensemble_convergence"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

json_file = (
    OUT
    / "ensemble_convergence_gate.json"
)

csv_file = (
    OUT
    / "ensemble_convergence_gate.csv"
)

status_file = (
    OUT
    / "ENSEMBLE_CONVERGENCE_STATUS.txt"
)


json_file.write_text(
    json.dumps(
        summary,
        indent=2,
    )
    + "\n"
)


rows = []

for record in summary:

    for obs in record[
        "observables"
    ]:

        rows.append({
            "Temperature_K":
                record[
                    "temperature_K"
                ],
            "Global_t0_index":
                record[
                    "global_t0_index"
                ],
            "Global_t0_step":
                record[
                    "global_t0_step"
                ],
            "Global_t0_ns":
                record[
                    "global_t0_ns"
                ],
            "Observable":
                obs[
                    "observable"
                ],
            "TauInt_samples":
                obs[
                    "tau_int_samples"
                ],
            "TauInt_ps":
                obs[
                    "tau_int_ps"
                ],
            "StatisticalInefficiency":
                obs[
                    "statistical_inefficiency"
                ],
            "N_eff":
                obs[
                    "N_eff"
                ],
            "BlockDrift":
                obs[
                    "max_adjacent_block_mean_drift"
                ],
            "Stable":
                obs[
                    "stable"
                ],
            "Spacing_frames":
                record[
                    "recommended_spacing_frames"
                ],
            "Spacing_steps":
                record[
                    "recommended_spacing_steps"
                ],
            "Spacing_ps":
                record[
                    "recommended_spacing_ps"
                ],
            "MaxIndependentSnapshots":
                record[
                    "maximum_independent_snapshots"
                ],
        })


pd.DataFrame(
    rows
).to_csv(
    csv_file,
    index=False,
)


all_stable = (
    len(summary)
    == len(TEMPERATURES)
    and
    all(
        record[
            "simulation_stable"
        ]
        for record
        in summary
    )
)


status = (
    "PHASE5_ENSEMBLE_CONVERGENCE_PASS"
    if all_stable
    else
    "PHASE5_ENSEMBLE_CONVERGENCE_NOT_YET_PASS"
)


status_file.write_text(
    status + "\n"
)


elapsed_campaign = (
    time.time()
    - campaign_start
)


print()
print(
    "="*90
)

print(
    "PHASE5-E15"
)

print(
    "FULL-RESOLUTION STATISTICAL CONVERGENCE GATE"
)

print(
    "="*90
)

print(
    f"RESULTS={json_file}"
)

print(
    f"STATUS={status}"
)

print(
    f"TOTAL_ELAPSED_SECONDS="
    f"{elapsed_campaign:.2f}"
)

print(
    "="*90
)
