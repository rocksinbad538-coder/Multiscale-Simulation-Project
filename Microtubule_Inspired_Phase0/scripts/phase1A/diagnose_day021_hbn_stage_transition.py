#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_ITP = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "protocol_inputs/topology/hbn_bonded_mobile_release.itp"
)

HBN_FIRST = 1
HBN_LAST = 1680

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)

NONFINITE_PATTERN = re.compile(
    r"(?<![A-Za-z])"
    r"(?:nan|[-+]?inf(?:inity)?)"
    r"(?![A-Za-z])",
    re.IGNORECASE,
)

HARMLESS_EPSILON_RF = re.compile(
    r"^\s*epsilon-rf\s*=\s*"
    r"(?:inf|infinity)\s*$",
    re.IGNORECASE,
)

SERIOUS_PATTERN = re.compile(
    r"LINCS WARNING|"
    r"Fatal error|"
    r"SETTLE.*(?:error|failed|cannot)|"
    r"SHAKE.*(?:failed|did not converge)|"
    r"constraint.*(?:error|failed)",
    re.IGNORECASE,
)


def resolve_path(value: str) -> Path:
    path = Path(value)

    if path.is_absolute():
        return path

    return PROJECT_ROOT / path


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(PROJECT_ROOT)
        )
    except ValueError:
        return str(path.resolve())


def read_gro(
    path: Path,
) -> tuple[np.ndarray, np.ndarray]:
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
            for line in lines[
                2 : 2 + natoms
            ]
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


def minimum_image(
    vectors: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        vectors
        - box
        * np.round(
            vectors / box
        )
    )


def kabsch_residuals(
    reference: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
    reference_centered = (
        reference
        - reference.mean(axis=0)
    )

    target_centered = (
        target
        - target.mean(axis=0)
    )

    covariance = (
        target_centered.T
        @ reference_centered
    )

    u_matrix, _, vt_matrix = (
        np.linalg.svd(covariance)
    )

    correction = np.eye(3)

    correction[2, 2] = np.sign(
        np.linalg.det(
            u_matrix @ vt_matrix
        )
    )

    rotation = (
        u_matrix
        @ correction
        @ vt_matrix
    )

    aligned = (
        target_centered
        @ rotation
    )

    return np.linalg.norm(
        aligned - reference_centered,
        axis=1,
    )


def parse_itp(
    path: Path,
) -> tuple[
    list[tuple[int, int, float]],
    list[tuple[int, int, int, float]],
]:
    section = ""

    bonds: list[
        tuple[int, int, float]
    ] = []

    angles: list[
        tuple[int, int, int, float]
    ] = []

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

        match = SECTION_PATTERN.match(line)

        if match:
            section = (
                match.group(1)
                .strip()
                .lower()
            )
            continue

        if line.startswith("#"):
            continue

        fields = line.split()

        if section == "bonds" and len(fields) >= 5:
            try:
                bonds.append(
                    (
                        int(fields[0]),
                        int(fields[1]),
                        float(fields[3]),
                    )
                )
            except ValueError:
                pass

        elif section == "angles" and len(fields) >= 6:
            try:
                angles.append(
                    (
                        int(fields[0]),
                        int(fields[1]),
                        int(fields[2]),
                        float(fields[4]),
                    )
                )
            except ValueError:
                pass

    if not bonds or not angles:
        raise RuntimeError(
            "Could not parse HBN bonds and angles"
        )

    return bonds, angles


def bond_lengths(
    positions: np.ndarray,
    box: np.ndarray,
    bonds: list[tuple[int, int, float]],
) -> np.ndarray:
    vectors = np.array(
        [
            positions[j - 1]
            - positions[i - 1]
            for i, j, _ in bonds
        ],
        dtype=float,
    )

    vectors = minimum_image(
        vectors,
        box,
    )

    return np.linalg.norm(
        vectors,
        axis=1,
    )


def angle_values(
    positions: np.ndarray,
    box: np.ndarray,
    angles: list[
        tuple[int, int, int, float]
    ],
) -> np.ndarray:
    values = []

    for i, j, k, _ in angles:
        vector_ji = minimum_image(
            (
                positions[i - 1]
                - positions[j - 1]
            )[None, :],
            box,
        )[0]

        vector_jk = minimum_image(
            (
                positions[k - 1]
                - positions[j - 1]
            )[None, :],
            box,
        )[0]

        denominator = (
            np.linalg.norm(vector_ji)
            * np.linalg.norm(vector_jk)
        )

        if denominator <= 0.0:
            raise RuntimeError(
                "Zero-length angle vector"
            )

        cosine = float(
            np.dot(
                vector_ji,
                vector_jk,
            )
            / denominator
        )

        cosine = max(
            -1.0,
            min(
                1.0,
                cosine,
            ),
        )

        values.append(
            math.degrees(
                math.acos(cosine)
            )
        )

    return np.array(
        values,
        dtype=float,
    )


def audit_logs(
    paths: tuple[Path, ...],
) -> tuple[int, int, int]:
    harmless = 0
    harmful = 0
    serious = 0

    for path in paths:
        for line in path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines():
            if SERIOUS_PATTERN.search(line):
                serious += 1

            if not NONFINITE_PATTERN.search(line):
                continue

            if HARMLESS_EPSILON_RF.match(line):
                harmless += 1
            else:
                harmful += 1

    return harmless, harmful, serious


def audit_xvg(
    run_root: Path,
) -> tuple[int, int]:
    rows = 0
    nonfinite = 0

    for path in sorted(
        run_root.glob("*.xvg")
    ):
        for raw_line in path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines():
            line = raw_line.strip()

            if (
                not line
                or line.startswith("#")
                or line.startswith("@")
            ):
                continue

            rows += 1

            for token in line.split():
                try:
                    value = float(token)
                except ValueError:
                    continue

                if not math.isfinite(value):
                    nonfinite += 1

    return rows, nonfinite


def quantile(
    values: np.ndarray,
    probability: float,
) -> float:
    return float(
        np.quantile(
            values,
            probability,
        )
    )


def write_summary(
    path: Path,
    row: dict[str, object],
) -> None:
    with path.open(
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


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage-name",
        required=True,
    )

    parser.add_argument(
        "--previous-gro",
        required=True,
    )

    parser.add_argument(
        "--current-gro",
        required=True,
    )

    parser.add_argument(
        "--mdrun-log",
        required=True,
    )

    parser.add_argument(
        "--console-log",
        required=True,
    )

    parser.add_argument(
        "--run-root",
        required=True,
    )

    parser.add_argument(
        "--hbn-itp",
        default=str(DEFAULT_ITP),
    )

    args = parser.parse_args()

    previous_gro = resolve_path(
        args.previous_gro
    )

    current_gro = resolve_path(
        args.current_gro
    )

    mdrun_log = resolve_path(
        args.mdrun_log
    )

    console_log = resolve_path(
        args.console_log
    )

    run_root = resolve_path(
        args.run_root
    )

    hbn_itp = resolve_path(
        args.hbn_itp
    )

    required = (
        previous_gro,
        current_gro,
        mdrun_log,
        console_log,
        hbn_itp,
    )

    missing = [
        path
        for path in required
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing required inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    previous, previous_box = read_gro(
        previous_gro
    )

    current, current_box = read_gro(
        current_gro
    )

    if not np.allclose(
        previous_box,
        current_box,
        atol=1.0e-8,
    ):
        raise RuntimeError(
            "NVT box changed"
        )

    previous_hbn = previous[
        HBN_FIRST - 1 : HBN_LAST
    ]

    current_hbn = current[
        HBN_FIRST - 1 : HBN_LAST
    ]

    displacement_vectors = minimum_image(
        current_hbn - previous_hbn,
        current_box,
    )

    displacements = np.linalg.norm(
        displacement_vectors,
        axis=1,
    )

    unwrapped_current = (
        previous_hbn
        + displacement_vectors
    )

    aligned_residuals = kabsch_residuals(
        previous_hbn,
        unwrapped_current,
    )

    bonds, angles = parse_itp(
        hbn_itp
    )

    bond_equilibrium = np.array(
        [
            equilibrium
            for _, _, equilibrium in bonds
        ],
        dtype=float,
    )

    angle_equilibrium = np.array(
        [
            equilibrium
            for _, _, _, equilibrium
            in angles
        ],
        dtype=float,
    )

    previous_bonds = bond_lengths(
        previous_hbn,
        previous_box,
        bonds,
    )

    current_bonds = bond_lengths(
        current_hbn,
        current_box,
        bonds,
    )

    previous_angles = angle_values(
        previous_hbn,
        previous_box,
        angles,
    )

    current_angles = angle_values(
        current_hbn,
        current_box,
        angles,
    )

    bond_deviation = np.abs(
        current_bonds
        - bond_equilibrium
    )

    angle_deviation = np.abs(
        current_angles
        - angle_equilibrium
    )

    bond_change = np.abs(
        current_bonds
        - previous_bonds
    )

    angle_change = np.abs(
        current_angles
        - previous_angles
    )

    harmless, harmful, serious = audit_logs(
        (
            mdrun_log,
            console_log,
        )
    )

    xvg_rows, xvg_nonfinite = audit_xvg(
        run_root
    )

    blocked_reasons = []
    review_flags = []

    if harmful:
        blocked_reasons.append(
            "harmful non-finite log value"
        )

    if serious:
        blocked_reasons.append(
            "serious instability signature"
        )

    if xvg_nonfinite:
        blocked_reasons.append(
            "non-finite XVG value"
        )

    if float(current_bonds.max()) > 0.20:
        blocked_reasons.append(
            "HBN bond exceeds 0.20 nm"
        )

    if float(
        np.sqrt(
            np.mean(
                aligned_residuals ** 2
            )
        )
    ) > 0.08:
        review_flags.append(
            "aligned HBN RMS exceeds 0.08 nm"
        )

    if quantile(
        bond_deviation,
        0.99,
    ) > 0.015:
        review_flags.append(
            "q99 bond deviation exceeds 0.015 nm"
        )

    if float(
        bond_deviation.max()
    ) > 0.03:
        review_flags.append(
            "maximum bond deviation exceeds 0.03 nm"
        )

    if quantile(
        angle_deviation,
        0.95,
    ) > 20.0:
        review_flags.append(
            "q95 angle deviation exceeds 20 degrees"
        )

    if float(
        angle_deviation.max()
    ) > 60.0:
        review_flags.append(
            "maximum angle deviation exceeds 60 degrees"
        )

    if blocked_reasons:
        screen = "BLOCKED"
    elif review_flags:
        screen = "REVIEW"
    else:
        screen = "STABLE_CANDIDATE"

    summary = {
        "stage": args.stage_name,
        "structural_screen": screen,
        "harmless_epsilon_rf_matches": harmless,
        "harmful_nonfinite_matches": harmful,
        "serious_instability_matches": serious,
        "numeric_xvg_rows": xvg_rows,
        "numeric_xvg_nonfinite_values": xvg_nonfinite,
        "HBN_incremental_rms_nm": float(
            np.sqrt(
                np.mean(
                    displacements ** 2
                )
            )
        ),
        "HBN_incremental_q95_nm": quantile(
            displacements,
            0.95,
        ),
        "HBN_incremental_q99_nm": quantile(
            displacements,
            0.99,
        ),
        "HBN_incremental_max_nm": float(
            displacements.max()
        ),
        "HBN_aligned_rms_nm": float(
            np.sqrt(
                np.mean(
                    aligned_residuals ** 2
                )
            )
        ),
        "HBN_aligned_max_nm": float(
            aligned_residuals.max()
        ),
        "HBN_atoms_above_0p10_nm": int(
            np.count_nonzero(
                displacements > 0.10
            )
        ),
        "HBN_atoms_above_0p15_nm": int(
            np.count_nonzero(
                displacements > 0.15
            )
        ),
        "HBN_bond_count": len(bonds),
        "HBN_bond_length_mean_nm": float(
            current_bonds.mean()
        ),
        "HBN_bond_length_min_nm": float(
            current_bonds.min()
        ),
        "HBN_bond_length_max_nm": float(
            current_bonds.max()
        ),
        "HBN_bond_deviation_q95_nm": quantile(
            bond_deviation,
            0.95,
        ),
        "HBN_bond_deviation_q99_nm": quantile(
            bond_deviation,
            0.99,
        ),
        "HBN_bond_deviation_max_nm": float(
            bond_deviation.max()
        ),
        "HBN_stage_bond_change_q95_nm": quantile(
            bond_change,
            0.95,
        ),
        "HBN_stage_bond_change_max_nm": float(
            bond_change.max()
        ),
        "HBN_angle_count": len(angles),
        "HBN_angle_deviation_q95_deg": quantile(
            angle_deviation,
            0.95,
        ),
        "HBN_angle_deviation_q99_deg": quantile(
            angle_deviation,
            0.99,
        ),
        "HBN_angle_deviation_max_deg": float(
            angle_deviation.max()
        ),
        "HBN_stage_angle_change_q95_deg": quantile(
            angle_change,
            0.95,
        ),
        "HBN_stage_angle_change_max_deg": float(
            angle_change.max()
        ),
        "review_flags": (
            " | ".join(review_flags)
        ),
        "blocked_reasons": (
            " | ".join(blocked_reasons)
        ),
    }

    summary_path = (
        run_root
        / f"{args.stage_name}_structural_summary.csv"
    )

    report_path = (
        run_root
        / f"{args.stage_name.upper()}_STRUCTURAL_DIAGNOSTIC_DAY021.md"
    )

    write_summary(
        summary_path,
        summary,
    )

    report_path.write_text(
        f"""# Day021 {args.stage_name} Structural Diagnostic

- Structural screen: **{screen}**
- HBN incremental RMS/max: {summary['HBN_incremental_rms_nm']:.8f}/{summary['HBN_incremental_max_nm']:.8f} nm
- HBN aligned RMS/max: {summary['HBN_aligned_rms_nm']:.8f}/{summary['HBN_aligned_max_nm']:.8f} nm
- Bond length mean/min/max: {summary['HBN_bond_length_mean_nm']:.8f}/{summary['HBN_bond_length_min_nm']:.8f}/{summary['HBN_bond_length_max_nm']:.8f} nm
- Bond deviation q95/q99/max: {summary['HBN_bond_deviation_q95_nm']:.8f}/{summary['HBN_bond_deviation_q99_nm']:.8f}/{summary['HBN_bond_deviation_max_nm']:.8f} nm
- Angle deviation q95/q99/max: {summary['HBN_angle_deviation_q95_deg']:.4f}/{summary['HBN_angle_deviation_q99_deg']:.4f}/{summary['HBN_angle_deviation_max_deg']:.4f} degrees
- Serious instability signatures: {serious}
- Harmful non-finite values: {harmful}
- Non-finite XVG values: {xvg_nonfinite}

The unrestrained stage remains unauthorized until this diagnostic is interpreted.
""",
        encoding="utf-8",
    )

    print(
        "Day021 HBN stage-transition diagnostic completed."
    )

    print(
        f"Stage: {args.stage_name}"
    )

    print(
        f"Serious instability signatures: {serious}"
    )

    print(
        f"Harmful non-finite values: {harmful}"
    )

    print(
        f"Non-finite XVG values: {xvg_nonfinite}"
    )

    print(
        "HBN incremental RMS/q95/q99/max: "
        f"{summary['HBN_incremental_rms_nm']:.8f}/"
        f"{summary['HBN_incremental_q95_nm']:.8f}/"
        f"{summary['HBN_incremental_q99_nm']:.8f}/"
        f"{summary['HBN_incremental_max_nm']:.8f} nm"
    )

    print(
        "HBN aligned RMS/max: "
        f"{summary['HBN_aligned_rms_nm']:.8f}/"
        f"{summary['HBN_aligned_max_nm']:.8f} nm"
    )

    print(
        "HBN bond length mean/min/max: "
        f"{summary['HBN_bond_length_mean_nm']:.8f}/"
        f"{summary['HBN_bond_length_min_nm']:.8f}/"
        f"{summary['HBN_bond_length_max_nm']:.8f} nm"
    )

    print(
        "HBN bond deviation q95/q99/max: "
        f"{summary['HBN_bond_deviation_q95_nm']:.8f}/"
        f"{summary['HBN_bond_deviation_q99_nm']:.8f}/"
        f"{summary['HBN_bond_deviation_max_nm']:.8f} nm"
    )

    print(
        "HBN angle deviation q95/q99/max: "
        f"{summary['HBN_angle_deviation_q95_deg']:.4f}/"
        f"{summary['HBN_angle_deviation_q99_deg']:.4f}/"
        f"{summary['HBN_angle_deviation_max_deg']:.4f} deg"
    )

    print(
        f"Structural screen: {screen}"
    )

    if review_flags:
        print(
            "Review flags: "
            + " | ".join(review_flags)
        )

    if blocked_reasons:
        print(
            "Blocking reasons: "
            + " | ".join(blocked_reasons)
        )

    print(
        "Stage05 authorized: NO"
    )

    print(
        f"Wrote: {relative(report_path)}"
    )


if __name__ == "__main__":
    main()
