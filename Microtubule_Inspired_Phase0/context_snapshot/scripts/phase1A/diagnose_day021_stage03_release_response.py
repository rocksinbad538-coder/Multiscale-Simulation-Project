#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/03_nvt_k1000_2ps"
)

START_GRO = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/02_nvt_k10000_1ps/"
    "02_nvt_k10000_1ps.gro"
)

REFERENCE_GRO = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute/"
    "nvt_100ps_frozenSolute.gro"
)

FINAL_GRO = (
    RUN_ROOT
    / "03_nvt_k1000_2ps.gro"
)

MDRUN_LOG = (
    RUN_ROOT
    / "03_nvt_k1000_2ps.log"
)

CONSOLE_LOG = (
    RUN_ROOT
    / "03_nvt_k1000_2ps_mdrun_console.log"
)

HBN_ITP = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "protocol_inputs/topology/"
    "hbn_bonded_mobile_release.itp"
)

SUMMARY_CSV = (
    RUN_ROOT
    / "stage03_release_response_diagnostic.csv"
)

TOP_DISPLACEMENTS_CSV = (
    RUN_ROOT
    / "stage03_top_hbn_displacements.csv"
)

BOND_DIAGNOSTIC_CSV = (
    RUN_ROOT
    / "stage03_hbn_bond_diagnostic.csv"
)

REPORT_MD = (
    RUN_ROOT
    / "STAGE03_RELEASE_RESPONSE_DIAGNOSTIC_DAY021.md"
)

HBN_FIRST = 1
HBN_LAST = 1680

PYR_FIRST = 1681
PYR_LAST = 1784

INCREMENTAL_THRESHOLD_NM = 0.08
CUMULATIVE_THRESHOLD_NM = 0.10

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


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(
                PROJECT_ROOT
            )
        )
    except ValueError:
        return str(path.resolve())


def require_inputs() -> None:
    required = (
        START_GRO,
        REFERENCE_GRO,
        FINAL_GRO,
        MDRUN_LOG,
        CONSOLE_LOG,
        HBN_ITP,
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
            "Missing or empty required inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def read_gro(
    path: Path,
) -> tuple[
    list[dict[str, object]],
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    natoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]

    if len(atom_lines) != natoms:
        raise RuntimeError(
            f"Atom-count mismatch in {path}"
        )

    atoms: list[dict[str, object]] = []

    for global_index, line in enumerate(
        atom_lines,
        start=1,
    ):
        atoms.append(
            {
                "global_index": global_index,
                "residue_number": int(
                    line[0:5]
                ),
                "residue_name": (
                    line[5:10].strip()
                ),
                "atom_name": (
                    line[10:15].strip()
                ),
                "gro_atom_number": int(
                    line[15:20]
                ),
                "position": np.array(
                    [
                        float(line[20:28]),
                        float(line[28:36]),
                        float(line[36:44]),
                    ],
                    dtype=float,
                ),
            }
        )

    box_fields = [
        float(value)
        for value in lines[
            2 + natoms
        ].split()
    ]

    if len(box_fields) < 3:
        raise RuntimeError(
            f"Invalid GRO box in {path}"
        )

    box = np.array(
        box_fields[:3],
        dtype=float,
    )

    return atoms, box


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


def group_positions(
    atoms: list[dict[str, object]],
    first_atom: int,
    last_atom: int,
) -> np.ndarray:
    return np.array(
        [
            atoms[index - 1]["position"]
            for index in range(
                first_atom,
                last_atom + 1,
            )
        ],
        dtype=float,
    )


def kabsch_residuals(
    reference: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
    reference_center = (
        reference.mean(axis=0)
    )

    target_center = (
        target.mean(axis=0)
    )

    reference_centered = (
        reference
        - reference_center
    )

    target_centered = (
        target
        - target_center
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

    differences = (
        aligned
        - reference_centered
    )

    return np.linalg.norm(
        differences,
        axis=1,
    )


def displacement_analysis(
    source_positions: np.ndarray,
    target_positions: np.ndarray,
    box: np.ndarray,
) -> dict[str, object]:
    displacement_vectors = (
        minimum_image(
            target_positions
            - source_positions,
            box,
        )
    )

    magnitudes = np.linalg.norm(
        displacement_vectors,
        axis=1,
    )

    mean_vector = (
        displacement_vectors.mean(
            axis=0
        )
    )

    translation_removed_vectors = (
        displacement_vectors
        - mean_vector
    )

    translation_removed = np.linalg.norm(
        translation_removed_vectors,
        axis=1,
    )

    aligned_residuals = (
        kabsch_residuals(
            source_positions,
            target_positions,
        )
    )

    return {
        "vectors": displacement_vectors,
        "magnitudes": magnitudes,
        "mean_vector": mean_vector,
        "translation_removed": (
            translation_removed
        ),
        "aligned_residuals": (
            aligned_residuals
        ),
        "rms": float(
            np.sqrt(
                np.mean(
                    magnitudes ** 2
                )
            )
        ),
        "maximum": float(
            magnitudes.max()
        ),
        "maximum_local_index": int(
            magnitudes.argmax()
        ),
        "median": float(
            np.median(magnitudes)
        ),
        "q95": float(
            np.quantile(
                magnitudes,
                0.95,
            )
        ),
        "q99": float(
            np.quantile(
                magnitudes,
                0.99,
            )
        ),
        "translation_nm": float(
            np.linalg.norm(
                mean_vector
            )
        ),
        "translation_removed_rms": float(
            np.sqrt(
                np.mean(
                    translation_removed ** 2
                )
            )
        ),
        "translation_removed_max": float(
            translation_removed.max()
        ),
        "aligned_rms": float(
            np.sqrt(
                np.mean(
                    aligned_residuals ** 2
                )
            )
        ),
        "aligned_max": float(
            aligned_residuals.max()
        ),
        "count_above_0p05": int(
            np.count_nonzero(
                magnitudes > 0.05
            )
        ),
        "count_above_0p08": int(
            np.count_nonzero(
                magnitudes > 0.08
            )
        ),
        "count_above_0p10": int(
            np.count_nonzero(
                magnitudes > 0.10
            )
        ),
    }


def parse_hbn_bonds(
    path: Path,
) -> list[tuple[int, int]]:
    current_section = ""
    bonds: list[tuple[int, int]] = []

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

        section_match = (
            SECTION_PATTERN.match(line)
        )

        if section_match:
            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            continue

        if current_section != "bonds":
            continue

        if line.startswith("#"):
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        try:
            atom_i = int(fields[0])
            atom_j = int(fields[1])
        except ValueError:
            continue

        bonds.append(
            (
                atom_i,
                atom_j,
            )
        )

    if not bonds:
        raise RuntimeError(
            "No HBN bonds were parsed"
        )

    return bonds


def bond_lengths(
    positions: np.ndarray,
    bonds: list[tuple[int, int]],
    box: np.ndarray,
) -> np.ndarray:
    vectors = np.array(
        [
            positions[j - 1]
            - positions[i - 1]
            for i, j in bonds
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


def audit_logs() -> tuple[
    list[tuple[str, int, str]],
    list[tuple[str, int, str]],
    list[tuple[str, int, str]],
]:
    harmless_hits = []
    harmful_nonfinite_hits = []
    serious_hits = []

    for path in (
        MDRUN_LOG,
        CONSOLE_LOG,
    ):
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            if SERIOUS_PATTERN.search(line):
                serious_hits.append(
                    (
                        path.name,
                        line_number,
                        line.strip(),
                    )
                )

            if not NONFINITE_PATTERN.search(
                line
            ):
                continue

            hit = (
                path.name,
                line_number,
                line.strip(),
            )

            if HARMLESS_EPSILON_RF.match(
                line
            ):
                harmless_hits.append(hit)
            else:
                harmful_nonfinite_hits.append(
                    hit
                )

    return (
        harmless_hits,
        harmful_nonfinite_hits,
        serious_hits,
    )


def audit_xvg() -> tuple[int, int]:
    numeric_rows = 0
    nonfinite_values = 0

    for path in sorted(
        RUN_ROOT.glob("*.xvg")
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

            numeric_rows += 1

            for token in line.split():
                try:
                    value = float(token)
                except ValueError:
                    continue

                if not math.isfinite(value):
                    nonfinite_values += 1

    return (
        numeric_rows,
        nonfinite_values,
    )


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No data available for {path}"
        )

    fieldnames: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for key in row:
            if key in seen:
                continue

            seen.add(key)
            fieldnames.append(key)

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    key: row.get(key, "")
                    for key in fieldnames
                }
            )


def main() -> None:
    require_inputs()

    start_atoms, start_box = read_gro(
        START_GRO
    )

    reference_atoms, reference_box = (
        read_gro(
            REFERENCE_GRO
        )
    )

    final_atoms, final_box = read_gro(
        FINAL_GRO
    )

    if not (
        len(start_atoms)
        == len(reference_atoms)
        == len(final_atoms)
        == 68320
    ):
        raise RuntimeError(
            "Unexpected atom count"
        )

    if not np.allclose(
        start_box,
        final_box,
        atol=1.0e-8,
    ):
        raise RuntimeError(
            "NVT box changed unexpectedly"
        )

    hbn_start = group_positions(
        start_atoms,
        HBN_FIRST,
        HBN_LAST,
    )

    hbn_reference = group_positions(
        reference_atoms,
        HBN_FIRST,
        HBN_LAST,
    )

    hbn_final = group_positions(
        final_atoms,
        HBN_FIRST,
        HBN_LAST,
    )

    pyr_start = group_positions(
        start_atoms,
        PYR_FIRST,
        PYR_LAST,
    )

    pyr_reference = group_positions(
        reference_atoms,
        PYR_FIRST,
        PYR_LAST,
    )

    pyr_final = group_positions(
        final_atoms,
        PYR_FIRST,
        PYR_LAST,
    )

    hbn_incremental = (
        displacement_analysis(
            hbn_start,
            hbn_final,
            final_box,
        )
    )

    hbn_cumulative = (
        displacement_analysis(
            hbn_reference,
            hbn_final,
            final_box,
        )
    )

    pyr_incremental = (
        displacement_analysis(
            pyr_start,
            pyr_final,
            final_box,
        )
    )

    pyr_cumulative = (
        displacement_analysis(
            pyr_reference,
            pyr_final,
            final_box,
        )
    )

    bonds = parse_hbn_bonds(
        HBN_ITP
    )

    start_bonds = bond_lengths(
        hbn_start,
        bonds,
        start_box,
    )

    final_bonds = bond_lengths(
        hbn_final,
        bonds,
        final_box,
    )

    bond_changes = (
        final_bonds
        - start_bonds
    )

    absolute_bond_changes = np.abs(
        bond_changes
    )

    top_indices = np.argsort(
        hbn_incremental["magnitudes"]
    )[::-1][:20]

    top_rows: list[
        dict[str, object]
    ] = []

    for rank, local_index in enumerate(
        top_indices,
        start=1,
    ):
        global_index = (
            HBN_FIRST
            + int(local_index)
        )

        atom = final_atoms[
            global_index - 1
        ]

        top_rows.append(
            {
                "rank": rank,
                "global_atom_index": (
                    global_index
                ),
                "atom_name": (
                    atom["atom_name"]
                ),
                "residue_name": (
                    atom["residue_name"]
                ),
                "incremental_displacement_nm": float(
                    hbn_incremental[
                        "magnitudes"
                    ][local_index]
                ),
                "cumulative_displacement_nm": float(
                    hbn_cumulative[
                        "magnitudes"
                    ][local_index]
                ),
                "translation_removed_incremental_nm": float(
                    hbn_incremental[
                        "translation_removed"
                    ][local_index]
                ),
                "aligned_incremental_residual_nm": float(
                    hbn_incremental[
                        "aligned_residuals"
                    ][local_index]
                ),
            }
        )

    bond_order = np.argsort(
        absolute_bond_changes
    )[::-1][:20]

    bond_rows: list[
        dict[str, object]
    ] = []

    for rank, bond_index in enumerate(
        bond_order,
        start=1,
    ):
        atom_i, atom_j = bonds[
            int(bond_index)
        ]

        bond_rows.append(
            {
                "rank": rank,
                "atom_i": atom_i,
                "atom_j": atom_j,
                "start_length_nm": float(
                    start_bonds[bond_index]
                ),
                "final_length_nm": float(
                    final_bonds[bond_index]
                ),
                "change_nm": float(
                    bond_changes[bond_index]
                ),
                "absolute_change_nm": float(
                    absolute_bond_changes[
                        bond_index
                    ]
                ),
            }
        )

    (
        harmless_hits,
        harmful_nonfinite_hits,
        serious_hits,
    ) = audit_logs()

    xvg_rows, xvg_nonfinite = (
        audit_xvg()
    )

    blocked_reasons = []
    revise_reasons = []

    if harmful_nonfinite_hits:
        blocked_reasons.append(
            "harmful non-finite log values"
        )

    if serious_hits:
        blocked_reasons.append(
            "serious instability signatures"
        )

    if xvg_nonfinite:
        blocked_reasons.append(
            "non-finite XVG values"
        )

    if (
        hbn_incremental["maximum"]
        > INCREMENTAL_THRESHOLD_NM
    ):
        revise_reasons.append(
            "HBN incremental maximum "
            "displacement exceeds 0.08 nm"
        )

    if (
        hbn_cumulative["maximum"]
        > CUMULATIVE_THRESHOLD_NM
    ):
        revise_reasons.append(
            "HBN cumulative maximum "
            "displacement exceeds 0.10 nm"
        )

    if blocked_reasons:
        decision = "BLOCKED"
    elif revise_reasons:
        decision = "REVISE"
    else:
        decision = "PASS"

    summary = {
        "stage": "03_nvt_k1000_2ps",
        "decision": decision,
        "harmless_epsilon_rf_matches": (
            len(harmless_hits)
        ),
        "harmful_nonfinite_matches": (
            len(harmful_nonfinite_hits)
        ),
        "serious_instability_matches": (
            len(serious_hits)
        ),
        "numeric_xvg_rows": xvg_rows,
        "numeric_xvg_nonfinite_values": (
            xvg_nonfinite
        ),
        "HBN_incremental_rms_nm": (
            hbn_incremental["rms"]
        ),
        "HBN_incremental_median_nm": (
            hbn_incremental["median"]
        ),
        "HBN_incremental_q95_nm": (
            hbn_incremental["q95"]
        ),
        "HBN_incremental_q99_nm": (
            hbn_incremental["q99"]
        ),
        "HBN_incremental_max_nm": (
            hbn_incremental["maximum"]
        ),
        "HBN_incremental_max_atom": (
            HBN_FIRST
            + hbn_incremental[
                "maximum_local_index"
            ]
        ),
        "HBN_atoms_above_0p05_nm": (
            hbn_incremental[
                "count_above_0p05"
            ]
        ),
        "HBN_atoms_above_0p08_nm": (
            hbn_incremental[
                "count_above_0p08"
            ]
        ),
        "HBN_atoms_above_0p10_nm": (
            hbn_incremental[
                "count_above_0p10"
            ]
        ),
        "HBN_translation_nm": (
            hbn_incremental[
                "translation_nm"
            ]
        ),
        "HBN_translation_removed_rms_nm": (
            hbn_incremental[
                "translation_removed_rms"
            ]
        ),
        "HBN_translation_removed_max_nm": (
            hbn_incremental[
                "translation_removed_max"
            ]
        ),
        "HBN_aligned_rms_nm": (
            hbn_incremental[
                "aligned_rms"
            ]
        ),
        "HBN_aligned_max_nm": (
            hbn_incremental[
                "aligned_max"
            ]
        ),
        "HBN_cumulative_max_nm": (
            hbn_cumulative["maximum"]
        ),
        "PYR_incremental_rms_nm": (
            pyr_incremental["rms"]
        ),
        "PYR_incremental_max_nm": (
            pyr_incremental["maximum"]
        ),
        "PYR_aligned_rms_nm": (
            pyr_incremental[
                "aligned_rms"
            ]
        ),
        "PYR_aligned_max_nm": (
            pyr_incremental[
                "aligned_max"
            ]
        ),
        "HBN_bond_count": len(bonds),
        "HBN_start_bond_mean_nm": float(
            start_bonds.mean()
        ),
        "HBN_final_bond_mean_nm": float(
            final_bonds.mean()
        ),
        "HBN_max_absolute_bond_change_nm": float(
            absolute_bond_changes.max()
        ),
        "HBN_bond_change_q95_nm": float(
            np.quantile(
                absolute_bond_changes,
                0.95,
            )
        ),
        "revise_reasons": (
            " | ".join(revise_reasons)
        ),
        "blocked_reasons": (
            " | ".join(blocked_reasons)
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    write_csv(
        TOP_DISPLACEMENTS_CSV,
        top_rows,
    )

    write_csv(
        BOND_DIAGNOSTIC_CSV,
        bond_rows,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Stage03 Release-Response Diagnostic\n\n"
        )

        handle.write(
            f"- Decision: **{decision}**\n"
        )

        handle.write(
            f"- Harmless `epsilon-rf = inf` matches: "
            f"{len(harmless_hits)}\n"
        )

        handle.write(
            f"- Harmful non-finite matches: "
            f"{len(harmful_nonfinite_hits)}\n"
        )

        handle.write(
            f"- Serious instability matches: "
            f"{len(serious_hits)}\n"
        )

        handle.write(
            f"- Non-finite XVG values: "
            f"{xvg_nonfinite}\n\n"
        )

        handle.write(
            "## HBN displacement response\n\n"
        )

        handle.write(
            f"- Incremental RMS/max: "
            f"{hbn_incremental['rms']:.8f}/"
            f"{hbn_incremental['maximum']:.8f} nm\n"
        )

        handle.write(
            f"- Incremental q95/q99: "
            f"{hbn_incremental['q95']:.8f}/"
            f"{hbn_incremental['q99']:.8f} nm\n"
        )

        handle.write(
            f"- Atoms above 0.05/0.08/0.10 nm: "
            f"{hbn_incremental['count_above_0p05']}/"
            f"{hbn_incremental['count_above_0p08']}/"
            f"{hbn_incremental['count_above_0p10']}\n"
        )

        handle.write(
            f"- Rigid translation magnitude: "
            f"{hbn_incremental['translation_nm']:.8f} nm\n"
        )

        handle.write(
            f"- Translation-removed RMS/max: "
            f"{hbn_incremental['translation_removed_rms']:.8f}/"
            f"{hbn_incremental['translation_removed_max']:.8f} nm\n"
        )

        handle.write(
            f"- Kabsch-aligned RMS/max: "
            f"{hbn_incremental['aligned_rms']:.8f}/"
            f"{hbn_incremental['aligned_max']:.8f} nm\n\n"
        )

        handle.write(
            "## HBN bonded-network response\n\n"
        )

        handle.write(
            f"- Bonds analyzed: {len(bonds)}\n"
        )

        handle.write(
            f"- Mean bond length start/final: "
            f"{start_bonds.mean():.8f}/"
            f"{final_bonds.mean():.8f} nm\n"
        )

        handle.write(
            f"- Maximum absolute bond change: "
            f"{absolute_bond_changes.max():.8f} nm\n"
        )

        handle.write(
            f"- q95 absolute bond change: "
            f"{np.quantile(absolute_bond_changes, 0.95):.8f} nm\n\n"
        )

        handle.write(
            "Stage04 remains unauthorized pending interpretation.\n"
        )

    print(
        "Day021 Stage03 release-response diagnostic completed."
    )

    print(
        f"Harmless epsilon-rf matches: "
        f"{len(harmless_hits)}"
    )

    print(
        f"Harmful non-finite matches: "
        f"{len(harmful_nonfinite_hits)}"
    )

    print(
        f"Serious instability signatures: "
        f"{len(serious_hits)}"
    )

    print(
        f"Non-finite XVG values: "
        f"{xvg_nonfinite}"
    )

    print(
        "HBN incremental RMS/median/q95/q99/max: "
        f"{hbn_incremental['rms']:.8f}/"
        f"{hbn_incremental['median']:.8f}/"
        f"{hbn_incremental['q95']:.8f}/"
        f"{hbn_incremental['q99']:.8f}/"
        f"{hbn_incremental['maximum']:.8f} nm"
    )

    print(
        "HBN atoms above 0.05/0.08/0.10 nm: "
        f"{hbn_incremental['count_above_0p05']}/"
        f"{hbn_incremental['count_above_0p08']}/"
        f"{hbn_incremental['count_above_0p10']}"
    )

    print(
        "HBN rigid translation / aligned RMS/max: "
        f"{hbn_incremental['translation_nm']:.8f}/"
        f"{hbn_incremental['aligned_rms']:.8f}/"
        f"{hbn_incremental['aligned_max']:.8f} nm"
    )

    print(
        "HBN bond mean start/final: "
        f"{start_bonds.mean():.8f}/"
        f"{final_bonds.mean():.8f} nm"
    )

    print(
        "HBN bond-change q95/max: "
        f"{np.quantile(absolute_bond_changes, 0.95):.8f}/"
        f"{absolute_bond_changes.max():.8f} nm"
    )

    print(
        f"Stage03 diagnostic decision: "
        f"{decision}"
    )

    print(
        "Stage04 authorized: NO"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
