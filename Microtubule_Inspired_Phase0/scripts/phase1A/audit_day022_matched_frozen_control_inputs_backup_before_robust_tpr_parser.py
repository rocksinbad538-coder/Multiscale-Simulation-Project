#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

EXECUTION = PROTOCOL / "execution"

STAGE02 = "02_nvt_k10000_1ps"

STAGE02_RUN = (
    EXECUTION
    / STAGE02
)

STAGE02_TPR = (
    STAGE02_RUN
    / f"{STAGE02}.tpr"
)

STAGE02_MDP = (
    PROTOCOL
    / "protocol_inputs/mdp"
    / f"{STAGE02}.mdp"
)

INDEX_FILE = (
    PROTOCOL
    / "protocol_inputs"
    / "mobile_release_index.ndx"
)

TOPOLOGY_FILE = (
    PROTOCOL
    / "protocol_inputs/topology"
    / "hbn_pyrene_mobile_release.top"
)

ACCEPTED_FROZEN_ENDPOINT = (
    ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute/"
    "nvt_100ps_frozenSolute.gro"
)

OUTPUT_ROOT = (
    PROTOCOL
    / "matched_frozen_control_144ps"
    / "input_audit"
)

TPR_DUMP = (
    OUTPUT_ROOT
    / "stage02_initial_state_tpr_dump.txt"
)

COORDINATE_MATCH_CSV = (
    OUTPUT_ROOT
    / "stage02_initial_coordinate_candidate_matches.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "matched_frozen_control_input_audit_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "MATCHED_FROZEN_CONTROL_INPUT_AUDIT_DAY022.md"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26
SOLUTE_COUNT = HBN_COUNT + PYR_COUNT * PYR_ATOMS

TOTAL_ATOMS = 68320
WATER_ATOM_COUNT = TOTAL_ATOMS - SOLUTE_COUNT

FLOAT_PATTERN = (
    r"[-+]?"
    r"(?:"
    r"(?:\d+(?:\.\d*)?)"
    r"|"
    r"(?:\.\d+)"
    r")"
    r"(?:[eE][-+]?\d+)?"
)

VECTOR_PATTERN = re.compile(
    rf"^\s*(x|v)\[\s*(\d+)\s*\]\s*=\s*"
    rf"\{{\s*({FLOAT_PATTERN})\s*,\s*"
    rf"({FLOAT_PATTERN})\s*,\s*"
    rf"({FLOAT_PATTERN})\s*\}}"
)

BOX_PATTERN = re.compile(
    rf"^\s*box\[\s*(\d+)\s*\]\s*=\s*"
    rf"\{{\s*({FLOAT_PATTERN})\s*,\s*"
    rf"({FLOAT_PATTERN})\s*,\s*"
    rf"({FLOAT_PATTERN})\s*\}}"
)


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def locate_gmx() -> str:
    executable = shutil.which("gmx")

    if executable:
        return executable

    fallback = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if fallback.exists():
        return str(fallback)

    raise RuntimeError(
        "Could not locate GROMACS"
    )


def canonical_key(value: str) -> str:
    return (
        value.strip()
        .lower()
        .replace("_", "-")
    )


def parse_mdp(path: Path) -> dict[str, str]:
    settings: dict[str, str] = {}

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line or "=" not in line:
            continue

        key, value = line.split(
            "=",
            1,
        )

        settings[
            canonical_key(key)
        ] = value.strip()

    return settings


def parse_index_groups(path: Path) -> list[str]:
    groups = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.strip()

        if (
            line.startswith("[")
            and line.endswith("]")
        ):
            groups.append(
                line[1:-1].strip()
            )

    return groups


def dump_tpr(
    gmx: str,
) -> str:
    completed = subprocess.run(
        [
            gmx,
            "dump",
            "-s",
            str(STAGE02_TPR),
        ],
        cwd=STAGE02_RUN,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    TPR_DUMP.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "gmx dump failed. "
            f"See {TPR_DUMP}"
        )

    return completed.stdout


def parse_tpr_vectors(
    text: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    positions: dict[int, tuple[float, float, float]] = {}
    velocities: dict[int, tuple[float, float, float]] = {}
    box_rows: dict[int, tuple[float, float, float]] = {}

    for line in text.splitlines():
        vector_match = VECTOR_PATTERN.match(
            line
        )

        if vector_match is not None:
            kind = vector_match.group(1)
            index = int(
                vector_match.group(2)
            )

            vector = (
                float(
                    vector_match.group(3)
                ),
                float(
                    vector_match.group(4)
                ),
                float(
                    vector_match.group(5)
                ),
            )

            if kind == "x":
                positions[index] = vector
            else:
                velocities[index] = vector

            continue

        box_match = BOX_PATTERN.match(
            line
        )

        if box_match is not None:
            index = int(
                box_match.group(1)
            )

            box_rows[index] = (
                float(
                    box_match.group(2)
                ),
                float(
                    box_match.group(3)
                ),
                float(
                    box_match.group(4)
                ),
            )

    if len(positions) != TOTAL_ATOMS:
        raise RuntimeError(
            "Could not parse all TPR positions. "
            f"Expected {TOTAL_ATOMS}, found "
            f"{len(positions)}"
        )

    if len(velocities) != TOTAL_ATOMS:
        raise RuntimeError(
            "Could not parse all TPR velocities. "
            f"Expected {TOTAL_ATOMS}, found "
            f"{len(velocities)}"
        )

    if len(box_rows) < 3:
        raise RuntimeError(
            "Could not parse the TPR box matrix"
        )

    position_array = np.array(
        [
            positions[index]
            for index in range(
                TOTAL_ATOMS
            )
        ],
        dtype=np.float64,
    )

    velocity_array = np.array(
        [
            velocities[index]
            for index in range(
                TOTAL_ATOMS
            )
        ],
        dtype=np.float64,
    )

    box_matrix = np.array(
        [
            box_rows[index]
            for index in range(3)
        ],
        dtype=np.float64,
    )

    return (
        position_array,
        velocity_array,
        box_matrix,
    )


def read_gro_positions(
    path: Path,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    natoms = int(
        lines[1].strip()
    )

    if natoms != TOTAL_ATOMS:
        raise RuntimeError(
            f"Expected {TOTAL_ATOMS} atoms in "
            f"{path}; found {natoms}"
        )

    positions = np.empty(
        (natoms, 3),
        dtype=np.float64,
    )

    for atom_index, line in enumerate(
        lines[2 : 2 + natoms]
    ):
        positions[
            atom_index
        ] = (
            float(line[20:28]),
            float(line[28:36]),
            float(line[36:44]),
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

    box = np.array(
        box_values[:3],
        dtype=np.float64,
    )

    return positions, box


def minimum_image(
    displacement: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        displacement
        - box
        * np.round(
            displacement / box
        )
    )


def candidate_gro_files() -> list[Path]:
    candidates = []

    explicit = [
        ACCEPTED_FROZEN_ENDPOINT,
        (
            EXECUTION
            / "00_em_k100000"
            / "00_em_k100000.gro"
        ),
        (
            EXECUTION
            / "01_em_k10000"
            / "01_em_k10000.gro"
        ),
        (
            EXECUTION
            / STAGE02
            / f"{STAGE02}.gro"
        ),
    ]

    for path in explicit:
        if path.exists():
            candidates.append(path)

    for stage_pattern in (
        "00_*",
        "01_*",
        "02_*",
    ):
        for path in EXECUTION.glob(
            f"{stage_pattern}/*.gro"
        ):
            if path not in candidates:
                candidates.append(path)

    return sorted(
        candidates,
        key=lambda path: str(path),
    )


def compare_coordinate_candidates(
    tpr_positions: np.ndarray,
) -> list[dict[str, object]]:
    rows = []

    for path in candidate_gro_files():
        try:
            gro_positions, box = (
                read_gro_positions(path)
            )

            displacement = minimum_image(
                gro_positions - tpr_positions,
                box,
            )

            atom_distance = np.linalg.norm(
                displacement,
                axis=1,
            )

            rows.append(
                {
                    "candidate_gro": relative(path),
                    "status": "ANALYZED",
                    "all_atom_rms_nm": float(
                        np.sqrt(
                            np.mean(
                                atom_distance ** 2
                            )
                        )
                    ),
                    "all_atom_max_nm": float(
                        atom_distance.max()
                    ),
                    "water_atom_rms_nm": float(
                        np.sqrt(
                            np.mean(
                                atom_distance[
                                    SOLUTE_COUNT:
                                ]
                                ** 2
                            )
                        )
                    ),
                    "water_atom_max_nm": float(
                        atom_distance[
                            SOLUTE_COUNT:
                        ].max()
                    ),
                }
            )

        except Exception as error:
            rows.append(
                {
                    "candidate_gro": relative(path),
                    "status": "ERROR",
                    "error": str(error),
                }
            )

    analyzed = [
        row
        for row in rows
        if row["status"] == "ANALYZED"
    ]

    analyzed.sort(
        key=lambda row: float(
            row["all_atom_rms_nm"]
        )
    )

    errors = [
        row
        for row in rows
        if row["status"] != "ANALYZED"
    ]

    return analyzed + errors


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )


def array_hash(
    array: np.ndarray,
) -> str:
    canonical = np.asarray(
        array,
        dtype="<f8",
    )

    return hashlib.sha256(
        canonical.tobytes(
            order="C"
        )
    ).hexdigest()


def velocity_metrics(
    velocities: np.ndarray,
    start: int,
    stop: int,
) -> dict[str, float]:
    subset = velocities[
        start:stop
    ]

    speed = np.linalg.norm(
        subset,
        axis=1,
    )

    nonzero = np.any(
        np.abs(subset) > 0.0,
        axis=1,
    )

    return {
        "atom_count": int(
            len(subset)
        ),
        "nonzero_count": int(
            np.count_nonzero(
                nonzero
            )
        ),
        "nonzero_fraction": float(
            np.mean(nonzero)
        ),
        "speed_mean_nm_per_ps": float(
            speed.mean()
        ),
        "speed_std_nm_per_ps": float(
            speed.std()
        ),
        "speed_min_nm_per_ps": float(
            speed.min()
        ),
        "speed_max_nm_per_ps": float(
            speed.max()
        ),
    }


def main() -> None:
    for path in (
        STAGE02_TPR,
        STAGE02_MDP,
        INDEX_FILE,
        TOPOLOGY_FILE,
    ):
        require_file(path)

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    gmx = locate_gmx()

    mdp = parse_mdp(
        STAGE02_MDP
    )

    index_groups = parse_index_groups(
        INDEX_FILE
    )

    dump_text = dump_tpr(
        gmx
    )

    (
        positions,
        velocities,
        box_matrix,
    ) = parse_tpr_vectors(
        dump_text
    )

    coordinate_rows = (
        compare_coordinate_candidates(
            positions
        )
    )

    write_csv(
        COORDINATE_MATCH_CSV,
        coordinate_rows,
    )

    best_match = next(
        (
            row
            for row in coordinate_rows
            if row["status"] == "ANALYZED"
        ),
        None,
    )

    if best_match is None:
        raise RuntimeError(
            "No usable GRO coordinate candidate "
            "was identified"
        )

    hbn_metrics = velocity_metrics(
        velocities,
        0,
        HBN_COUNT,
    )

    pyr_metrics = velocity_metrics(
        velocities,
        HBN_COUNT,
        SOLUTE_COUNT,
    )

    water_metrics = velocity_metrics(
        velocities,
        SOLUTE_COUNT,
        TOTAL_ATOMS,
    )

    selected_mdp_fields = [
        "integrator",
        "dt",
        "nsteps",
        "continuation",
        "gen-vel",
        "gen-seed",
        "tcoupl",
        "tc-grps",
        "tau-t",
        "ref-t",
        "freezegrps",
        "freezedim",
        "comm-mode",
        "comm-grps",
        "nstcomm",
        "constraints",
        "constraint-algorithm",
        "pbc",
        "nstxout-compressed",
    ]

    summary: dict[str, object] = {
        "stage": STAGE02,
        "tpr": relative(
            STAGE02_TPR
        ),
        "source_mdp": relative(
            STAGE02_MDP
        ),
        "topology": relative(
            TOPOLOGY_FILE
        ),
        "index_file": relative(
            INDEX_FILE
        ),
        "atom_count": len(
            positions
        ),
        "solute_atom_count": (
            SOLUTE_COUNT
        ),
        "water_atom_count": (
            WATER_ATOM_COUNT
        ),
        "initial_coordinate_sha256": (
            array_hash(
                positions
            )
        ),
        "initial_water_velocity_sha256": (
            array_hash(
                velocities[
                    SOLUTE_COUNT:
                ]
            )
        ),
        "initial_all_velocity_sha256": (
            array_hash(
                velocities
            )
        ),
        "HBN_nonzero_velocity_fraction": (
            hbn_metrics[
                "nonzero_fraction"
            ]
        ),
        "PYR_nonzero_velocity_fraction": (
            pyr_metrics[
                "nonzero_fraction"
            ]
        ),
        "water_nonzero_velocity_fraction": (
            water_metrics[
                "nonzero_fraction"
            ]
        ),
        "water_speed_mean_nm_per_ps": (
            water_metrics[
                "speed_mean_nm_per_ps"
            ]
        ),
        "water_speed_std_nm_per_ps": (
            water_metrics[
                "speed_std_nm_per_ps"
            ]
        ),
        "water_speed_min_nm_per_ps": (
            water_metrics[
                "speed_min_nm_per_ps"
            ]
        ),
        "water_speed_max_nm_per_ps": (
            water_metrics[
                "speed_max_nm_per_ps"
            ]
        ),
        "best_coordinate_match": (
            best_match[
                "candidate_gro"
            ]
        ),
        "best_coordinate_match_rms_nm": (
            best_match[
                "all_atom_rms_nm"
            ]
        ),
        "best_coordinate_match_max_nm": (
            best_match[
                "all_atom_max_nm"
            ]
        ),
        "index_group_names": (
            " | ".join(
                index_groups
            )
        ),
        "matched_control_execution_authorized": (
            False
        ),
        "required_next_step": (
            "BUILD_CONTROL_TPR_AND_VERIFY_"
            "INITIAL_WATER_VELOCITY_HASH"
        ),
    }

    for field in selected_mdp_fields:
        summary[
            f"mdp_{field}"
        ] = mdp.get(
            field,
            "MISSING",
        )

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    report_lines = [
        "# Matched Frozen-Control Input Audit",
        "",
        "## Stage02 authoritative initial state",
        "",
        f"- TPR: `{relative(STAGE02_TPR)}`",
        f"- Atoms: **{len(positions)}**",
        (
            "- Initial coordinate SHA256: "
            f"`{summary['initial_coordinate_sha256']}`"
        ),
        (
            "- Initial water-velocity SHA256: "
            f"`{summary['initial_water_velocity_sha256']}`"
        ),
        "",
        "## Initial velocity populations",
        "",
        (
            "- HBN nonzero fraction: "
            f"{hbn_metrics['nonzero_fraction']:.6f}"
        ),
        (
            "- PYR nonzero fraction: "
            f"{pyr_metrics['nonzero_fraction']:.6f}"
        ),
        (
            "- Water nonzero fraction: "
            f"{water_metrics['nonzero_fraction']:.6f}"
        ),
        (
            "- Water speed mean/std: "
            f"{water_metrics['speed_mean_nm_per_ps']:.8f}/"
            f"{water_metrics['speed_std_nm_per_ps']:.8f} "
            "nm ps^-1"
        ),
        "",
        "## Closest GRO representation",
        "",
        (
            f"- File: `{best_match['candidate_gro']}`"
        ),
        (
            "- All-atom RMS/max difference from TPR coordinates: "
            f"{float(best_match['all_atom_rms_nm']):.8f}/"
            f"{float(best_match['all_atom_max_nm']):.8f} nm"
        ),
        "",
        "## Stage02 MDP settings",
        "",
    ]

    for field in selected_mdp_fields:
        report_lines.append(
            f"- `{field}`: `{mdp.get(field, 'MISSING')}`"
        )

    report_lines.extend(
        [
            "",
            "## Index groups",
            "",
            (
                "- "
                + ", ".join(
                    f"`{group}`"
                    for group in index_groups
                )
            ),
            "",
            "## Decision",
            "",
            (
                "- Matched frozen-control execution authorized: "
                "**NO**"
            ),
            (
                "- Next requirement: construct the frozen-control "
                "TPR from the authoritative Stage02 initial state "
                "and verify that its water-velocity SHA256 is "
                "identical before running MD."
            ),
        ]
    )

    REPORT_MD.write_text(
        "\n".join(
            report_lines
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        "Day022 matched frozen-control input audit completed."
    )

    print(
        "Stage02 initial atoms: "
        f"{len(positions)}/{TOTAL_ATOMS}"
    )

    print(
        "Stage02 MDP gen-vel / gen-seed / continuation: "
        f"{mdp.get('gen-vel', 'MISSING')} / "
        f"{mdp.get('gen-seed', 'MISSING')} / "
        f"{mdp.get('continuation', 'MISSING')}"
    )

    print(
        "Stage02 tc-grps: "
        f"{mdp.get('tc-grps', 'MISSING')}"
    )

    print(
        "Stage02 freezegrps / freezedim: "
        f"{mdp.get('freezegrps', 'MISSING')} / "
        f"{mdp.get('freezedim', 'MISSING')}"
    )

    print(
        "Initial HBN/PYR/water nonzero velocity fractions: "
        f"{hbn_metrics['nonzero_fraction']:.6f}/"
        f"{pyr_metrics['nonzero_fraction']:.6f}/"
        f"{water_metrics['nonzero_fraction']:.6f}"
    )

    print(
        "Initial water speed mean/std/min/max: "
        f"{water_metrics['speed_mean_nm_per_ps']:.8f}/"
        f"{water_metrics['speed_std_nm_per_ps']:.8f}/"
        f"{water_metrics['speed_min_nm_per_ps']:.8f}/"
        f"{water_metrics['speed_max_nm_per_ps']:.8f} nm/ps"
    )

    print(
        "Initial water velocity SHA256: "
        f"{summary['initial_water_velocity_sha256']}"
    )

    print(
        "Best GRO coordinate match: "
        f"{best_match['candidate_gro']}"
    )

    print(
        "Best GRO RMS/max difference: "
        f"{float(best_match['all_atom_rms_nm']):.8f}/"
        f"{float(best_match['all_atom_max_nm']):.8f} nm"
    )

    print(
        "Available index groups: "
        + " | ".join(
            index_groups
        )
    )

    print(
        "Matched frozen-control execution authorized: NO"
    )

    print(
        "Required next step: "
        "BUILD_CONTROL_TPR_AND_VERIFY_INITIAL_WATER_VELOCITY_HASH"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
