#!/usr/bin/env python3

from pathlib import Path
import json
import shutil
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
)

TEMPERATURES = [150, 200, 250, 300, 350]

TARGET_SNAPSHOTS = 50

GATE_FILE = (
    CAMPAIGN
    / "ensemble_convergence"
    / "ensemble_convergence_gate.json"
)

FINAL_GATE_FILE = (
    CAMPAIGN
    / "ensemble_convergence"
    / "PHASE5_FINAL_ENSEMBLE_CONVERGENCE.json"
)

STATUS_FILE = (
    CAMPAIGN
    / "ensemble_convergence"
    / "ENSEMBLE_CONVERGENCE_STATUS.txt"
)

OUTROOT = (
    CAMPAIGN
    / "representative_ensemble"
)


if not GATE_FILE.exists():
    raise RuntimeError(
        "Primary convergence gate does not exist."
    )

if not FINAL_GATE_FILE.exists():
    raise RuntimeError(
        "Final convergence decision does not exist. "
        "Run PHASE5-E19 first."
    )

if not STATUS_FILE.exists():
    raise RuntimeError(
        "Convergence status does not exist."
    )

status = STATUS_FILE.read_text().strip()

if status != "PHASE5_ENSEMBLE_CONVERGENCE_PASS":
    raise RuntimeError(
        f"Ensemble convergence gate has not passed: {status}"
    )


gate = json.loads(
    GATE_FILE.read_text()
)

final_gate = json.loads(
    FINAL_GATE_FILE.read_text()
)

if not final_gate.get("final_pass", False):
    raise RuntimeError(
        "Final ensemble convergence decision is not PASS."
    )

gate_by_temperature = {
    int(record["temperature_K"]): record
    for record in gate
}

final_by_temperature = {
    int(record["temperature_K"]): record
    for record in final_gate["temperatures"]
}


if OUTROOT.exists():
    shutil.rmtree(OUTROOT)

OUTROOT.mkdir(
    parents=True,
    exist_ok=True
)

manifest = []


def read_frames(path):

    frames = []

    with path.open() as f:

        while True:

            line = f.readline()

            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            step = int(f.readline())

            header = [
                line,
                str(step) + "\n",
            ]

            # NUMBER OF ATOMS header/value,
            # BOX header + 3 bounds,
            # ATOMS header = 7 additional lines.
            for _ in range(7):
                header.append(
                    f.readline()
                )

            natoms = int(header[3])

            atom_header = header[-1]

            if not (
                atom_header.startswith("ITEM: ATOMS")
                and "xu" in atom_header
                and "yu" in atom_header
                and "zu" in atom_header
            ):
                raise RuntimeError(
                    f"{path}: expected unwrapped xu/yu/zu coordinates."
                )

            atoms = [
                f.readline()
                for _ in range(natoms)
            ]

            frames.append({
                "step": step,
                "header": header,
                "atoms": atoms,
            })

    return frames


for T in TEMPERATURES:

    if T not in gate_by_temperature:
        raise RuntimeError(
            f"No convergence record for {T} K."
        )

    record = gate_by_temperature[T]

    if T not in final_by_temperature:
        raise RuntimeError(
            f"No final convergence decision for {T} K."
        )

    final_record = final_by_temperature[T]

    if not final_record["final_temperature_pass"]:
        raise RuntimeError(
            f"{T} K did not pass final convergence decision."
        )

    burnin_step = int(
        record["global_t0_step"]
    )

    spacing_steps = int(
        record["recommended_spacing_steps"]
    )

    dump = (
        CAMPAIGN
        / f"{T}K"
        / "production.xyz"
    )

    if not dump.exists():
        raise RuntimeError(
            f"Missing trajectory: {dump}"
        )

    frames = read_frames(dump)

    eligible = [
        frame
        for frame in frames
        if frame["step"] >= burnin_step
    ]

    if not eligible:
        raise RuntimeError(
            f"{T} K: no frames after burn-in."
        )

    # --------------------------------------------------------
    # Uniform temporal coverage of the complete equilibrated
    # region while preserving the minimum statistical spacing.
    # --------------------------------------------------------

    if len(eligible) < TARGET_SNAPSHOTS:
        raise RuntimeError(
            f"{T} K: only {len(eligible)} eligible frames "
            f"for {TARGET_SNAPSHOTS} requested snapshots."
        )

    # Select approximately uniformly across the full equilibrated
    # trajectory, including both temporal endpoints.
    target_indices = np.linspace(
        0,
        len(eligible)-1,
        TARGET_SNAPSHOTS,
        dtype=int,
    )

    selected = [
        eligible[i]
        for i in target_indices
    ]

    selected_steps = [
        int(frame["step"])
        for frame in selected
    ]

    if len(set(selected_steps)) != TARGET_SNAPSHOTS:
        raise RuntimeError(
            f"{T} K: duplicate selected timesteps."
        )

    deltas = [
        b-a
        for a,b in zip(
            selected_steps[:-1],
            selected_steps[1:],
        )
    ]

    minimum_actual_spacing = min(deltas)

    if minimum_actual_spacing < spacing_steps:
        raise RuntimeError(
            f"{T} K: temporal coverage selection violates "
            f"statistical spacing: actual="
            f"{minimum_actual_spacing}, required={spacing_steps}."
        )

    selected_span_steps = (
        selected_steps[-1]
        - selected_steps[0]
    )

    eligible_span_steps = (
        int(eligible[-1]["step"])
        - int(eligible[0]["step"])
    )

    temporal_coverage_fraction = (
        selected_span_steps
        / eligible_span_steps
        if eligible_span_steps > 0
        else 1.0
    )

    outdir = OUTROOT / f"{T}K"
    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    snapshot_records = []

    for n, frame in enumerate(
        selected,
        start=1
    ):

        outfile = (
            outdir
            / f"snapshot_{n:03d}.dump"
        )

        with outfile.open("w") as g:
            g.writelines(frame["header"])
            g.writelines(frame["atoms"])

        snapshot_records.append({
            "snapshot": n,
            "step": int(frame["step"]),
            "file": outfile.name,
        })

    manifest.append({
        "temperature_K": T,
        "burn_in_step": burnin_step,
        "burn_in_ns":
            float(record["global_t0_ns"]),
        "spacing_steps":
            spacing_steps,
        "spacing_ps":
            float(
                record["recommended_spacing_ps"]
            ),
        "target_snapshots":
            TARGET_SNAPSHOTS,
        "available_independent_snapshots":
            int(
                record[
                    "maximum_independent_snapshots"
                ]
            ),
        "selected_snapshots":
            len(snapshot_records),
        "minimum_required_spacing_steps":
            spacing_steps,
        "minimum_actual_spacing_steps":
            int(minimum_actual_spacing),
        "selected_first_step":
            int(selected_steps[0]),
        "selected_last_step":
            int(selected_steps[-1]),
        "temporal_coverage_fraction":
            float(temporal_coverage_fraction),
        "selection_strategy":
            "UNIFORM_FULL_EQUILIBRATED_WINDOW_WITH_MINIMUM_STATISTICAL_SPACING",
        "snapshots":
            snapshot_records,
    })

    print(
        f"{T} K: "
        f"burn-in={burnin_step} steps, "
        f"spacing={spacing_steps} steps, "
        f"selected={len(snapshot_records)}"
    )


manifest_file = (
    OUTROOT
    / "ensemble_manifest.json"
)

manifest_file.write_text(
    json.dumps(
        manifest,
        indent=2
    )
    + "\n"
)

print()
print("="*90)
print("PHASE5-D46")
print("STATISTICALLY-INDEPENDENT REPRESENTATIVE ENSEMBLE")
print("="*90)
print(manifest_file)
print("REPRESENTATIVE_ENSEMBLE_READY")
