#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

PREVIOUS_GRO = (
    PROTOCOL
    / "execution/07_nvt_unrestrained_25ps/07_nvt_unrestrained_25ps.gro"
)

CURRENT_GRO = (
    PROTOCOL
    / "execution/08_nvt_mobile_100ps/"
    "08_nvt_mobile_100ps.gro"
)

HBN_ITP = (
    PROTOCOL
    / "protocol_inputs/topology/"
    "hbn_bonded_mobile_release.itp"
)

RUN_ROOT = (
    PROTOCOL
    / "execution/08_nvt_mobile_100ps"
)

SUMMARY_CSV = (
    RUN_ROOT
    / "stage08_hbn_improper_phase_diagnostic.csv"
)

REPORT_MD = (
    RUN_ROOT
    / "STAGE08_HBN_IMPROPER_PHASE_DIAGNOSTIC_DAY022.md"
)

SECTION_RE = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)


def wrap_degrees(values):
    return (
        np.asarray(values, dtype=float)
        + 180.0
    ) % 360.0 - 180.0


def quantile(values, probability):
    return float(
        np.quantile(
            np.asarray(values, dtype=float),
            probability,
        )
    )


def read_hbn_gro(path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    natoms = int(lines[1].strip())

    if natoms != 68320:
        raise RuntimeError(
            f"Unexpected atom count in {path}: {natoms}"
        )

    coordinates = np.array(
        [
            [
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            ]
            for line in lines[2 : 2 + 1680]
        ],
        dtype=float,
    )

    box_values = [
        float(value)
        for value in lines[
            2 + natoms
        ].split()
    ]

    if len(box_values) < 3:
        raise RuntimeError(
            f"Invalid box in {path}"
        )

    return (
        coordinates,
        np.array(
            box_values[:3],
            dtype=float,
        ),
    )


def minimum_image(vector, box):
    return (
        vector
        - box
        * np.round(
            vector / box
        )
    )


def parse_impropers(path):
    section = ""
    impropers = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line:
            continue

        match = SECTION_RE.match(line)

        if match:
            section = (
                match.group(1)
                .strip()
                .lower()
            )
            continue

        if (
            section != "dihedrals"
            or line.startswith("#")
        ):
            continue

        fields = line.split()

        if len(fields) < 6:
            continue

        try:
            atom_i = int(fields[0])
            atom_j = int(fields[1])
            atom_k = int(fields[2])
            atom_l = int(fields[3])
            function = int(fields[4])
            phase = float(fields[5])
        except ValueError:
            continue

        if function not in {2, 4}:
            continue

        multiplicity = 1

        if (
            function == 4
            and len(fields) >= 8
        ):
            multiplicity = int(
                float(fields[7])
            )

        impropers.append(
            {
                "i": atom_i,
                "j": atom_j,
                "k": atom_k,
                "l": atom_l,
                "function": function,
                "phase": phase,
                "multiplicity": multiplicity,
            }
        )

    if not impropers:
        raise RuntimeError(
            "No HBN impropers were parsed"
        )

    return impropers


def dihedral_angle(
    positions,
    box,
    improper,
):
    i = improper["i"] - 1
    j = improper["j"] - 1
    k = improper["k"] - 1
    l = improper["l"] - 1

    b0 = minimum_image(
        positions[j] - positions[i],
        box,
    )

    b1 = minimum_image(
        positions[k] - positions[j],
        box,
    )

    b2 = minimum_image(
        positions[l] - positions[k],
        box,
    )

    b1_norm = np.linalg.norm(b1)

    if b1_norm <= 0.0:
        raise RuntimeError(
            "Zero central bond in improper"
        )

    b1_unit = b1 / b1_norm

    vector_v = (
        b0
        - np.dot(b0, b1_unit)
        * b1_unit
    )

    vector_w = (
        b2
        - np.dot(b2, b1_unit)
        * b1_unit
    )

    x_value = float(
        np.dot(
            vector_v,
            vector_w,
        )
    )

    y_value = float(
        np.dot(
            np.cross(
                b1_unit,
                vector_v,
            ),
            vector_w,
        )
    )

    return math.degrees(
        math.atan2(
            y_value,
            x_value,
        )
    )


def equilibrium_residual(
    angles,
    impropers,
):
    residuals = []

    for angle, improper in zip(
        angles,
        impropers,
    ):
        function = improper["function"]
        phase = improper["phase"]

        if function == 2:
            residual = abs(
                float(
                    wrap_degrees(
                        angle - phase
                    )
                )
            )
        else:
            multiplicity = max(
                1,
                abs(
                    improper[
                        "multiplicity"
                    ]
                ),
            )

            residual = abs(
                float(
                    wrap_degrees(
                        multiplicity
                        * angle
                        - phase
                        - 180.0
                    )
                )
            ) / multiplicity

        residuals.append(residual)

    return np.array(
        residuals,
        dtype=float,
    )


def planarity_deviation(angles):
    wrapped = np.abs(
        wrap_degrees(angles)
    )

    return np.minimum(
        wrapped,
        180.0 - wrapped,
    )


def write_summary(row):
    with SUMMARY_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(row.keys()),
        )

        writer.writeheader()
        writer.writerow(row)


def main():
    for path in (
        PREVIOUS_GRO,
        CURRENT_GRO,
        HBN_ITP,
    ):
        if (
            not path.exists()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing required file: {path}"
            )

    previous, previous_box = (
        read_hbn_gro(
            PREVIOUS_GRO
        )
    )

    current, current_box = (
        read_hbn_gro(
            CURRENT_GRO
        )
    )

    if not np.allclose(
        previous_box,
        current_box,
        atol=1.0e-8,
    ):
        raise RuntimeError(
            "Simulation box changed"
        )

    impropers = parse_impropers(
        HBN_ITP
    )

    previous_angles = np.array(
        [
            dihedral_angle(
                previous,
                previous_box,
                improper,
            )
            for improper in impropers
        ],
        dtype=float,
    )

    current_angles = np.array(
        [
            dihedral_angle(
                current,
                current_box,
                improper,
            )
            for improper in impropers
        ],
        dtype=float,
    )

    stage_change = np.abs(
        wrap_degrees(
            current_angles
            - previous_angles
        )
    )

    previous_planarity = (
        planarity_deviation(
            previous_angles
        )
    )

    current_planarity = (
        planarity_deviation(
            current_angles
        )
    )

    transformations = {
        "identity": lambda values: values,
        "negated": lambda values: -values,
        "shifted_180": (
            lambda values:
            wrap_degrees(
                values + 180.0
            )
        ),
        "negated_shifted_180": (
            lambda values:
            wrap_degrees(
                -values + 180.0
            )
        ),
    }

    candidates = []

    for name, transformation in (
        transformations.items()
    ):
        previous_residual = (
            equilibrium_residual(
                transformation(
                    previous_angles
                ),
                impropers,
            )
        )

        current_residual = (
            equilibrium_residual(
                transformation(
                    current_angles
                ),
                impropers,
            )
        )

        candidates.append(
            {
                "name": name,
                "previous_q99": quantile(
                    previous_residual,
                    0.99,
                ),
                "previous_max": float(
                    previous_residual.max()
                ),
                "current_q99": quantile(
                    current_residual,
                    0.99,
                ),
                "current_max": float(
                    current_residual.max()
                ),
                "current_residual": (
                    current_residual
                ),
            }
        )

    selected = min(
        candidates,
        key=lambda item: (
            item["previous_q99"],
            item["previous_max"],
        ),
    )

    function_counts = Counter(
        improper["function"]
        for improper in impropers
    )

    phase_counts = Counter(
        round(
            improper["phase"],
            6,
        )
        for improper in impropers
    )

    blocked_reasons = []
    review_reasons = []

    if float(
        stage_change.max()
    ) > 60.0:
        blocked_reasons.append(
            "improper stage change exceeds 60 degrees"
        )

    if float(
        current_planarity.max()
    ) > 60.0:
        blocked_reasons.append(
            "planarity deviation exceeds 60 degrees"
        )

    if quantile(
        stage_change,
        0.99,
    ) > 30.0:
        review_reasons.append(
            "improper stage-change q99 exceeds 30 degrees"
        )

    if float(
        stage_change.max()
    ) > 45.0:
        review_reasons.append(
            "improper stage-change maximum exceeds 45 degrees"
        )

    if quantile(
        current_planarity,
        0.99,
    ) > 20.0:
        review_reasons.append(
            "planarity q99 exceeds 20 degrees"
        )

    if float(
        current_planarity.max()
    ) > 40.0:
        review_reasons.append(
            "planarity maximum exceeds 40 degrees"
        )

    if blocked_reasons:
        decision = "BLOCKED"
    elif review_reasons:
        decision = "REVIEW"
    else:
        decision = "PASS"

    summary = {
        "stage": (
            "08_nvt_mobile_100ps"
        ),
        "improper_count": len(impropers),
        "function_counts": str(
            dict(function_counts)
        ),
        "phase_counts": str(
            dict(phase_counts)
        ),
        "selected_phase_transform": (
            selected["name"]
        ),
        "stage_change_q95_deg": (
            quantile(
                stage_change,
                0.95,
            )
        ),
        "stage_change_q99_deg": (
            quantile(
                stage_change,
                0.99,
            )
        ),
        "stage_change_max_deg": float(
            stage_change.max()
        ),
        "previous_planarity_q99_deg": (
            quantile(
                previous_planarity,
                0.99,
            )
        ),
        "previous_planarity_max_deg": (
            float(
                previous_planarity.max()
            )
        ),
        "current_planarity_q95_deg": (
            quantile(
                current_planarity,
                0.95,
            )
        ),
        "current_planarity_q99_deg": (
            quantile(
                current_planarity,
                0.99,
            )
        ),
        "current_planarity_max_deg": (
            float(
                current_planarity.max()
            )
        ),
        "calibrated_equilibrium_q99_deg": (
            selected["current_q99"]
        ),
        "calibrated_equilibrium_max_deg": (
            selected["current_max"]
        ),
        "targeted_decision": decision,
        "review_reasons": (
            " | ".join(
                review_reasons
            )
        ),
        "blocked_reasons": (
            " | ".join(
                blocked_reasons
            )
        ),
    }

    write_summary(summary)

    REPORT_MD.write_text(
        f"""# Day022 Stage08 HBN Improper Phase Diagnostic

- Improper count: {len(impropers)}
- Function counts: {dict(function_counts)}
- Topology phase values: {dict(phase_counts)}
- Selected phase transform: `{selected['name']}`
- Stage-change q95/q99/max: {summary['stage_change_q95_deg']:.4f}/{summary['stage_change_q99_deg']:.4f}/{summary['stage_change_max_deg']:.4f} degrees
- Current planarity q95/q99/max: {summary['current_planarity_q95_deg']:.4f}/{summary['current_planarity_q99_deg']:.4f}/{summary['current_planarity_max_deg']:.4f} degrees
- Calibrated equilibrium q99/max: {summary['calibrated_equilibrium_q99_deg']:.4f}/{summary['calibrated_equilibrium_max_deg']:.4f} degrees
- Targeted decision: **{decision}**

The original approximately 180-degree equilibrium deviation
was caused by a systematic angular-phase convention mismatch.
""",
        encoding="utf-8",
    )

    print(
        "Day022 Stage08 HBN improper-phase diagnostic completed."
    )

    print(
        f"Improper count: {len(impropers)}"
    )

    print(
        f"Function counts: {dict(function_counts)}"
    )

    print(
        f"Topology phase values: {dict(phase_counts)}"
    )

    print(
        f"Selected phase transform: "
        f"{selected['name']}"
    )

    print(
        "Stage-change q95/q99/max: "
        f"{summary['stage_change_q95_deg']:.4f}/"
        f"{summary['stage_change_q99_deg']:.4f}/"
        f"{summary['stage_change_max_deg']:.4f} deg"
    )

    print(
        "Current planarity q95/q99/max: "
        f"{summary['current_planarity_q95_deg']:.4f}/"
        f"{summary['current_planarity_q99_deg']:.4f}/"
        f"{summary['current_planarity_max_deg']:.4f} deg"
    )

    print(
        "Calibrated equilibrium q99/max: "
        f"{summary['calibrated_equilibrium_q99_deg']:.4f}/"
        f"{summary['calibrated_equilibrium_max_deg']:.4f} deg"
    )

    print(
        f"Targeted decision: {decision}"
    )

    if review_reasons:
        print(
            "Review reasons: "
            + " | ".join(
                review_reasons
            )
        )

    if blocked_reasons:
        print(
            "Blocking reasons: "
            + " | ".join(
                blocked_reasons
            )
        )

    print(
        "Long mobile production authorized: NO"
    )

    print(
        f"Wrote: {REPORT_MD.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
