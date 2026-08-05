#!/usr/bin/env python3
"""
DAY040 / D040-A3

Local 37-atom mapping feasibility audit.

Purpose
-------
Evaluate whether the adopted Phase 1A-F 37-atom working charge model
can be mapped into the selected accepted GROMACS HBN topology.

Critical compatibility issue
----------------------------
The adopted charge model contains B17 N14 H6, while the selected HBN
component currently contains only B0 and N0 atom types. Therefore this
block:

1. tests direct parent-index anchors encoded as P:<index>;
2. determines the most plausible topology-index convention;
3. performs rigid alignment from direct anchors;
4. constructs a one-to-one B/N mapping candidate by element and spatial
   proximity;
5. audits whether the six adopted H atoms already exist in the target
   HBN topology;
6. keeps all topology modification blocked.

No topology or coordinate file is modified.
No charge is written.
No GROMACS calculation is executed.
"""

from __future__ import annotations

import csv
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

A2_REPORT = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_target_selection"
    / "QM_F06_UPPER_V7A_R1_PHASE1A_G_TARGET_SELECTION.json"
)

EXPECTED_A2_DECISION = (
    "D040_A2_FORCE_FIELD_TARGET_SELECTION_PASS_"
    "LOCAL_37_ATOM_MAPPING_DESIGN_AUTHORIZED"
)

CLOSURE_POINTER = (
    ROOT
    / "runs/phase1A"
    / "LATEST_PHASE1A_F_CHARGE_MODEL_CLOSURE.txt"
)

TARGET_DIR = (
    ROOT
    / "parameters/phase1A/accepted"
    / "hybrid_hbnBonded_kang2000_improperGeo100_validated"
)

TARGET_HBN_ITP = (
    TARGET_DIR
    / "hbn_bonded_candidate_kang2000_improperGeo100.itp"
)

TARGET_COORDINATES = (
    ROOT
    / "parameters/phase1A/accepted"
    / "hybrid_hydrated_gap45_pyr5shift_clean032"
    / "hbn_pyrene_4_tip4p2005_solvated_gap45_pyr5shift_clean032.gro"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_local_mapping_feasibility"
)

MAPPING_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_LOCAL_37_ATOM_MAPPING_CANDIDATE.csv"
)

ANCHOR_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_PARENT_INDEX_ANCHOR_AUDIT.csv"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_LOCAL_37_ATOM_MAPPING_FEASIBILITY.json"
)

EXPECTED_CHARGE_COUNT = 37
EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 6,
}

EXPECTED_HBN_TOPOLOGY_ATOMS = 1680

DIRECT_ANCHOR_PATTERN = re.compile(
    r"^P:(\d+)$"
)

MAX_ANCHOR_RMSD_A_FOR_PLAUSIBILITY = 0.25
MAX_MAPPING_DISTANCE_A_FOR_REVIEW = 0.35
MAX_MAPPING_DISTANCE_A_FOR_AUTHORIZATION = 0.20


def parse_itp_atoms(path: Path) -> list[dict]:
    require_file(path)

    active_section = None
    atoms = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        content = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not content:
            continue

        section_match = re.match(
            r"^\[\s*([^\]]+?)\s*\]$",
            content,
        )

        if section_match:
            active_section = (
                section_match.group(1)
                .strip()
                .lower()
            )
            continue

        if content.startswith("#"):
            continue

        if active_section != "atoms":
            continue

        tokens = content.split()

        if len(tokens) < 7:
            continue

        try:
            atom_index_1based = int(
                tokens[0]
            )
            charge_e = float(
                tokens[6]
            )
            mass = (
                float(tokens[7])
                if len(tokens) >= 8
                else None
            )
        except ValueError:
            continue

        atom_type = tokens[1]

        if atom_type.upper().startswith("B"):
            element = "B"
        elif atom_type.upper().startswith("N"):
            element = "N"
        elif atom_type.upper().startswith("H"):
            element = "H"
        else:
            element = "UNKNOWN"

        atoms.append(
            {
                "topology_index_1based": (
                    atom_index_1based
                ),
                "topology_index_0based": (
                    atom_index_1based - 1
                ),
                "atom_type": atom_type,
                "residue_number": int(
                    tokens[2]
                ),
                "residue_name": tokens[3],
                "atom_name": tokens[4],
                "charge_group": int(
                    tokens[5]
                ),
                "charge_e": charge_e,
                "mass": mass,
                "element": element,
            }
        )

    return atoms


def parse_gro(path: Path) -> dict:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            "GRO file is too short"
        )

    try:
        atom_count = int(
            lines[1].strip()
        )
    except ValueError as error:
        raise RuntimeError(
            "Could not parse GRO atom count"
        ) from error

    atom_lines = lines[
        2:2 + atom_count
    ]

    if len(atom_lines) != atom_count:
        raise RuntimeError(
            "GRO atom-line count mismatch"
        )

    atoms = []

    for sequence_index, line in enumerate(
        atom_lines
    ):
        if len(line) < 44:
            raise RuntimeError(
                "Malformed GRO atom record at "
                f"sequence {sequence_index}"
            )

        residue_number = int(
            line[0:5]
        )
        residue_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        file_atom_number = int(
            line[15:20]
        )

        x_nm = float(line[20:28])
        y_nm = float(line[28:36])
        z_nm = float(line[36:44])

        atoms.append(
            {
                "gro_sequence_index_0based": (
                    sequence_index
                ),
                "gro_sequence_index_1based": (
                    sequence_index + 1
                ),
                "gro_file_atom_number": (
                    file_atom_number
                ),
                "residue_number": (
                    residue_number
                ),
                "residue_name": (
                    residue_name
                ),
                "atom_name": atom_name,
                "x_A": 10.0 * x_nm,
                "y_A": 10.0 * y_nm,
                "z_A": 10.0 * z_nm,
            }
        )

    box_tokens = (
        lines[2 + atom_count].split()
        if len(lines) > 2 + atom_count
        else []
    )

    return {
        "title": lines[0],
        "atom_count": atom_count,
        "atoms": atoms,
        "box_tokens": box_tokens,
    }


def load_adopted_charges(
    path: Path,
) -> list[dict]:
    require_file(path)

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    rows.sort(
        key=lambda row: int(
            row[
                "real_atom_sequence_index"
            ]
        )
    )

    parsed = []

    for row in rows:
        parsed.append(
            {
                "real_atom_sequence_index": int(
                    row[
                        "real_atom_sequence_index"
                    ]
                ),
                "original_atom_index_0based": int(
                    row[
                        "original_atom_index_0based"
                    ]
                ),
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row["atom_role"],
                "node_type": row["node_type"],
                "x_A": float(row["x_A"]),
                "y_A": float(row["y_A"]),
                "z_A": float(row["z_A"]),
                "adopted_working_charge_e": float(
                    row[
                        "adopted_working_charge_e"
                    ]
                ),
            }
        )

    return parsed


def kabsch_transform(
    source_xyz: np.ndarray,
    target_xyz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    if source_xyz.shape != target_xyz.shape:
        raise RuntimeError(
            "Kabsch source/target shape mismatch"
        )

    if len(source_xyz) < 3:
        raise RuntimeError(
            "At least three anchors are required"
        )

    source_centroid = np.mean(
        source_xyz,
        axis=0,
    )

    target_centroid = np.mean(
        target_xyz,
        axis=0,
    )

    source_centered = (
        source_xyz - source_centroid
    )

    target_centered = (
        target_xyz - target_centroid
    )

    covariance = (
        source_centered.T
        @ target_centered
    )

    u_matrix, _, vt_matrix = (
        np.linalg.svd(covariance)
    )

    determinant_sign = np.sign(
        np.linalg.det(
            u_matrix @ vt_matrix
        )
    )

    correction = np.eye(3)
    correction[2, 2] = (
        determinant_sign
        if determinant_sign != 0.0
        else 1.0
    )

    rotation = (
        u_matrix
        @ correction
        @ vt_matrix
    )

    translation = (
        target_centroid
        - source_centroid
        @ rotation
    )

    transformed = (
        source_xyz @ rotation
        + translation
    )

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (
                        transformed
                        - target_xyz
                    )
                    ** 2,
                    axis=1,
                )
            )
        )
    )

    return (
        rotation,
        translation,
        rmsd,
    )


def greedy_unique_assignment(
    query_rows: list[dict],
    transformed_xyz: np.ndarray,
    topology_rows: list[dict],
    topology_xyz: np.ndarray,
    locked_topology_indices: set[int],
) -> list[dict]:
    candidate_pairs = []

    for query_position, query_row in enumerate(
        query_rows
    ):
        query_element = query_row[
            "element"
        ]

        for topology_position, topology_row in enumerate(
            topology_rows
        ):
            if (
                topology_row["element"]
                != query_element
            ):
                continue

            topology_index_0based = (
                topology_row[
                    "topology_index_0based"
                ]
            )

            if topology_index_0based in locked_topology_indices:
                continue

            distance_A = float(
                np.linalg.norm(
                    transformed_xyz[
                        query_position
                    ]
                    - topology_xyz[
                        topology_position
                    ]
                )
            )

            candidate_pairs.append(
                (
                    distance_A,
                    query_position,
                    topology_position,
                )
            )

    candidate_pairs.sort(
        key=lambda item: item[0]
    )

    assigned_queries = set()
    assigned_topology_positions = set()
    assignments = []

    for (
        distance_A,
        query_position,
        topology_position,
    ) in candidate_pairs:
        if query_position in assigned_queries:
            continue

        if topology_position in assigned_topology_positions:
            continue

        assigned_queries.add(
            query_position
        )
        assigned_topology_positions.add(
            topology_position
        )

        assignments.append(
            {
                "query_position": (
                    query_position
                ),
                "topology_position": (
                    topology_position
                ),
                "distance_A": (
                    distance_A
                ),
            }
        )

    return assignments


print("=" * 100)
print("DAY040 / D040-A3 — LOCAL 37-ATOM MAPPING FEASIBILITY")
print("=" * 100)


print("\n[1] UPSTREAM AUTHORIZATION")

require_file(A2_REPORT)
require_file(CLOSURE_POINTER)

a2_report = load_json(
    A2_REPORT
)

if (
    a2_report.get("decision")
    != EXPECTED_A2_DECISION
):
    raise RuntimeError(
        "Unexpected D040-A2 decision.\n"
        f"Observed: {a2_report.get('decision')}"
    )

if (
    a2_report.get(
        "authorizations",
        {},
    ).get(
        "local_37_atom_mapping_design_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Local mapping design is not authorized"
    )

closure_dir = (
    ROOT
    / CLOSURE_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

charges_csv = (
    closure_dir
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.csv"
)

require_file(charges_csv)
require_file(TARGET_HBN_ITP)
require_file(TARGET_COORDINATES)

print("D040_A2_decision_gate = PASS")
print("local_mapping_design_authorization_gate = PASS")
print("topology_modification_blocked_gate = PASS")
print("charge_mapping_execution_blocked_gate = PASS")


print("\n[2] LOAD ADOPTED 37-ATOM CONTRACT")

charge_rows = load_adopted_charges(
    charges_csv
)

if len(charge_rows) != EXPECTED_CHARGE_COUNT:
    raise RuntimeError(
        "Expected 37 adopted charge rows"
    )

composition = {}

for row in charge_rows:
    composition[
        row["element"]
    ] = (
        composition.get(
            row["element"],
            0,
        )
        + 1
    )

if composition != EXPECTED_COMPOSITION:
    raise RuntimeError(
        "Unexpected adopted composition.\n"
        f"Observed: {composition}"
    )

charge_xyz = np.asarray(
    [
        [
            row["x_A"],
            row["y_A"],
            row["z_A"],
        ]
        for row in charge_rows
    ],
    dtype=float,
)

print(
    f"adopted_atom_count = "
    f"{len(charge_rows)}"
)
print(
    f"adopted_composition = "
    f"{composition}"
)
print(
    f"adopted_charge_sum_e = "
    f"{sum(row['adopted_working_charge_e'] for row in charge_rows):.16g}"
)


print("\n[3] LOAD SELECTED HBN TOPOLOGY AND COORDINATES")

topology_atoms = parse_itp_atoms(
    TARGET_HBN_ITP
)

gro = parse_gro(
    TARGET_COORDINATES
)

if (
    len(topology_atoms)
    != EXPECTED_HBN_TOPOLOGY_ATOMS
):
    raise RuntimeError(
        "Unexpected HBN topology atom count.\n"
        f"Observed: {len(topology_atoms)}"
    )

if gro["atom_count"] < len(
    topology_atoms
):
    raise RuntimeError(
        "Coordinate file has fewer atoms than HBN topology"
    )

hbn_gro_atoms = gro["atoms"][
    :len(topology_atoms)
]

for position, (
    topology_atom,
    gro_atom,
) in enumerate(
    zip(
        topology_atoms,
        hbn_gro_atoms,
    )
):
    topology_atom[
        "gro_sequence_index_0based"
    ] = gro_atom[
        "gro_sequence_index_0based"
    ]

    topology_atom[
        "gro_sequence_index_1based"
    ] = gro_atom[
        "gro_sequence_index_1based"
    ]

    topology_atom[
        "gro_atom_name"
    ] = gro_atom[
        "atom_name"
    ]

    topology_atom[
        "gro_residue_name"
    ] = gro_atom[
        "residue_name"
    ]

    topology_atom["x_A"] = (
        gro_atom["x_A"]
    )

    topology_atom["y_A"] = (
        gro_atom["y_A"]
    )

    topology_atom["z_A"] = (
        gro_atom["z_A"]
    )

topology_composition = {}

for row in topology_atoms:
    topology_composition[
        row["element"]
    ] = (
        topology_composition.get(
            row["element"],
            0,
        )
        + 1
    )

topology_xyz = np.asarray(
    [
        [
            row["x_A"],
            row["y_A"],
            row["z_A"],
        ]
        for row in topology_atoms
    ],
    dtype=float,
)

print(
    f"hbn_topology_atom_count = "
    f"{len(topology_atoms)}"
)
print(
    f"coordinate_file_total_atom_count = "
    f"{gro['atom_count']}"
)
print(
    f"hbn_topology_composition = "
    f"{topology_composition}"
)
print(
    f"hbn_topology_hydrogen_count = "
    f"{topology_composition.get('H', 0)}"
)
print(
    f"coordinate_source = "
    f"{TARGET_COORDINATES.relative_to(ROOT)}"
)


print("\n[4] DIRECT PARENT-INDEX ANCHOR INVENTORY")

anchor_candidates = []

for charge_position, row in enumerate(
    charge_rows
):
    match = DIRECT_ANCHOR_PATTERN.match(
        row["atom_id"]
    )

    if not match:
        continue

    encoded_index = int(
        match.group(1)
    )

    anchor_candidates.append(
        {
            "charge_position": (
                charge_position
            ),
            "atom_id": row["atom_id"],
            "element": row["element"],
            "encoded_index": (
                encoded_index
            ),
        }
    )

print(
    f"direct_parent_anchor_count = "
    f"{len(anchor_candidates)}"
)

for row in anchor_candidates:
    print(
        f"atom_id={row['atom_id']} "
        f"element={row['element']} "
        f"encoded_index={row['encoded_index']}"
    )


print("\n[5] TEST INDEX CONVENTIONS")

convention_results = []

conventions = (
    (
        "P_INDEX_IS_TOPOLOGY_0BASED",
        lambda encoded: encoded,
    ),
    (
        "P_INDEX_IS_TOPOLOGY_1BASED",
        lambda encoded: encoded - 1,
    ),
)

for (
    convention_name,
    converter,
) in conventions:
    source_points = []
    target_points = []
    valid_anchor_rows = []
    invalid_count = 0
    element_mismatch_count = 0

    for anchor in anchor_candidates:
        topology_index_0based = (
            converter(
                anchor[
                    "encoded_index"
                ]
            )
        )

        if not (
            0
            <= topology_index_0based
            < len(topology_atoms)
        ):
            invalid_count += 1
            continue

        topology_row = topology_atoms[
            topology_index_0based
        ]

        if (
            topology_row["element"]
            != anchor["element"]
        ):
            element_mismatch_count += 1
            continue

        source_points.append(
            charge_xyz[
                anchor[
                    "charge_position"
                ]
            ]
        )

        target_points.append(
            topology_xyz[
                topology_index_0based
            ]
        )

        valid_anchor_rows.append(
            {
                **anchor,
                "topology_index_0based": (
                    topology_index_0based
                ),
                "topology_index_1based": (
                    topology_index_0based
                    + 1
                ),
                "topology_atom_type": (
                    topology_row[
                        "atom_type"
                    ]
                ),
                "topology_atom_name": (
                    topology_row[
                        "atom_name"
                    ]
                ),
            }
        )

    result = {
        "convention": convention_name,
        "valid_anchor_count": len(
            valid_anchor_rows
        ),
        "invalid_index_count": (
            invalid_count
        ),
        "element_mismatch_count": (
            element_mismatch_count
        ),
        "anchor_rows": valid_anchor_rows,
        "rotation": None,
        "translation": None,
        "anchor_RMSD_A": None,
    }

    if len(valid_anchor_rows) >= 3:
        (
            rotation,
            translation,
            anchor_rmsd_A,
        ) = kabsch_transform(
            np.asarray(
                source_points,
                dtype=float,
            ),
            np.asarray(
                target_points,
                dtype=float,
            ),
        )

        result[
            "rotation"
        ] = rotation

        result[
            "translation"
        ] = translation

        result[
            "anchor_RMSD_A"
        ] = anchor_rmsd_A

    convention_results.append(
        result
    )

    print(
        f"convention={convention_name} "
        f"valid_anchors={result['valid_anchor_count']} "
        f"invalid_indices={invalid_count} "
        f"element_mismatches={element_mismatch_count} "
        f"anchor_RMSD_A={result['anchor_RMSD_A']}"
    )

plausible_conventions = [
    result
    for result in convention_results
    if (
        result[
            "anchor_RMSD_A"
        ]
        is not None
    )
]

if not plausible_conventions:
    best_convention = None
else:
    best_convention = min(
        plausible_conventions,
        key=lambda result: (
            result[
                "anchor_RMSD_A"
            ]
        ),
    )

if best_convention is None:
    print(
        "best_index_convention = NONE"
    )
else:
    print(
        "best_index_convention = "
        f"{best_convention['convention']}"
    )
    print(
        "best_anchor_RMSD_A = "
        f"{best_convention['anchor_RMSD_A']:.16g}"
    )


print("\n[6] B/N ONE-TO-ONE MAPPING CANDIDATE")

mapping_records = []
anchor_records_output = []

if (
    best_convention is None
    or best_convention[
        "anchor_RMSD_A"
    ]
    > MAX_ANCHOR_RMSD_A_FOR_PLAUSIBILITY
):
    alignment_authorized = False
    transformed_charge_xyz = None

    print(
        "alignment_authorized = False"
    )
    print(
        "mapping_candidate_status = "
        "BLOCKED_BY_ANCHOR_ALIGNMENT"
    )

else:
    alignment_authorized = True

    rotation = best_convention[
        "rotation"
    ]

    translation = best_convention[
        "translation"
    ]

    transformed_charge_xyz = (
        charge_xyz @ rotation
        + translation
    )

    locked_topology_indices = set()
    directly_mapped_charge_positions = set()

    for anchor in best_convention[
        "anchor_rows"
    ]:
        charge_position = anchor[
            "charge_position"
        ]

        topology_index_0based = anchor[
            "topology_index_0based"
        ]

        locked_topology_indices.add(
            topology_index_0based
        )

        directly_mapped_charge_positions.add(
            charge_position
        )

        distance_A = float(
            np.linalg.norm(
                transformed_charge_xyz[
                    charge_position
                ]
                - topology_xyz[
                    topology_index_0based
                ]
            )
        )

        mapping_records.append(
            {
                "real_atom_sequence_index": (
                    charge_rows[
                        charge_position
                    ][
                        "real_atom_sequence_index"
                    ]
                ),
                "atom_id": (
                    charge_rows[
                        charge_position
                    ]["atom_id"]
                ),
                "element": (
                    charge_rows[
                        charge_position
                    ]["element"]
                ),
                "atom_role": (
                    charge_rows[
                        charge_position
                    ]["atom_role"]
                ),
                "adopted_working_charge_e": (
                    charge_rows[
                        charge_position
                    ][
                        "adopted_working_charge_e"
                    ]
                ),
                "mapping_status": (
                    "DIRECT_PARENT_INDEX_ANCHOR"
                ),
                "topology_index_0based": (
                    topology_index_0based
                ),
                "topology_index_1based": (
                    topology_index_0based
                    + 1
                ),
                "topology_atom_type": (
                    topology_atoms[
                        topology_index_0based
                    ]["atom_type"]
                ),
                "topology_atom_name": (
                    topology_atoms[
                        topology_index_0based
                    ]["atom_name"]
                ),
                "mapping_distance_A": (
                    distance_A
                ),
            }
        )

        anchor_records_output.append(
            {
                "convention": (
                    best_convention[
                        "convention"
                    ]
                ),
                "atom_id": (
                    charge_rows[
                        charge_position
                    ]["atom_id"]
                ),
                "element": (
                    charge_rows[
                        charge_position
                    ]["element"]
                ),
                "encoded_index": (
                    int(
                        DIRECT_ANCHOR_PATTERN.match(
                            charge_rows[
                                charge_position
                            ]["atom_id"]
                        ).group(1)
                    )
                ),
                "topology_index_0based": (
                    topology_index_0based
                ),
                "topology_index_1based": (
                    topology_index_0based
                    + 1
                ),
                "mapping_distance_A": (
                    distance_A
                ),
            }
        )

    remaining_charge_positions = [
        index
        for index, row in enumerate(
            charge_rows
        )
        if (
            row["element"]
            in {"B", "N"}
            and index
            not in directly_mapped_charge_positions
        )
    ]

    remaining_query_rows = [
        charge_rows[index]
        for index in remaining_charge_positions
    ]

    remaining_transformed_xyz = (
        transformed_charge_xyz[
            remaining_charge_positions
        ]
    )

    assignments = (
        greedy_unique_assignment(
            remaining_query_rows,
            remaining_transformed_xyz,
            topology_atoms,
            topology_xyz,
            locked_topology_indices,
        )
    )

    for assignment in assignments:
        local_query_position = (
            assignment[
                "query_position"
            ]
        )

        charge_position = (
            remaining_charge_positions[
                local_query_position
            ]
        )

        topology_position = (
            assignment[
                "topology_position"
            ]
        )

        charge_row = charge_rows[
            charge_position
        ]

        topology_row = topology_atoms[
            topology_position
        ]

        mapping_records.append(
            {
                "real_atom_sequence_index": (
                    charge_row[
                        "real_atom_sequence_index"
                    ]
                ),
                "atom_id": (
                    charge_row["atom_id"]
                ),
                "element": (
                    charge_row["element"]
                ),
                "atom_role": (
                    charge_row["atom_role"]
                ),
                "adopted_working_charge_e": (
                    charge_row[
                        "adopted_working_charge_e"
                    ]
                ),
                "mapping_status": (
                    "ALIGNED_NEAREST_SAME_ELEMENT"
                ),
                "topology_index_0based": (
                    topology_row[
                        "topology_index_0based"
                    ]
                ),
                "topology_index_1based": (
                    topology_row[
                        "topology_index_1based"
                    ]
                ),
                "topology_atom_type": (
                    topology_row[
                        "atom_type"
                    ]
                ),
                "topology_atom_name": (
                    topology_row[
                        "atom_name"
                    ]
                ),
                "mapping_distance_A": (
                    assignment[
                        "distance_A"
                    ]
                ),
            }
        )

    mapped_charge_positions = {
        int(
            record[
                "real_atom_sequence_index"
            ]
        )
        for record in mapping_records
    }

    for charge_row in charge_rows:
        if (
            charge_row[
                "real_atom_sequence_index"
            ]
            in mapped_charge_positions
        ):
            continue

        if charge_row["element"] == "H":
            status = (
                "UNMAPPED_H_TOPOLOGY_ATOM_ABSENT"
            )
        else:
            status = (
                "UNMAPPED_NO_UNIQUE_CANDIDATE"
            )

        mapping_records.append(
            {
                "real_atom_sequence_index": (
                    charge_row[
                        "real_atom_sequence_index"
                    ]
                ),
                "atom_id": (
                    charge_row["atom_id"]
                ),
                "element": (
                    charge_row["element"]
                ),
                "atom_role": (
                    charge_row["atom_role"]
                ),
                "adopted_working_charge_e": (
                    charge_row[
                        "adopted_working_charge_e"
                    ]
                ),
                "mapping_status": status,
                "topology_index_0based": "",
                "topology_index_1based": "",
                "topology_atom_type": "",
                "topology_atom_name": "",
                "mapping_distance_A": "",
            }
        )

    print(
        "alignment_authorized = True"
    )
    print(
        f"mapping_record_count = "
        f"{len(mapping_records)}"
    )


print("\n[7] MAPPING DIAGNOSTICS")

mapping_records.sort(
    key=lambda row: int(
        row[
            "real_atom_sequence_index"
        ]
    )
)

mapped_BN_records = [
    row
    for row in mapping_records
    if (
        row["element"]
        in {"B", "N"}
        and row[
            "mapping_status"
        ]
        in {
            "DIRECT_PARENT_INDEX_ANCHOR",
            "ALIGNED_NEAREST_SAME_ELEMENT",
        }
    )
]

unmapped_H_records = [
    row
    for row in mapping_records
    if (
        row["element"] == "H"
        and row[
            "mapping_status"
        ]
        == "UNMAPPED_H_TOPOLOGY_ATOM_ABSENT"
    )
]

BN_distances = np.asarray(
    [
        float(
            row[
                "mapping_distance_A"
            ]
        )
        for row in mapped_BN_records
    ],
    dtype=float,
) if mapped_BN_records else np.asarray(
    [],
    dtype=float,
)

topology_indices_mapped = [
    int(
        row[
            "topology_index_0based"
        ]
    )
    for row in mapped_BN_records
]

unique_topology_mapping_gate = (
    len(
        topology_indices_mapped
    )
    == len(
        set(
            topology_indices_mapped
        )
    )
)

mapped_BN_count = len(
    mapped_BN_records
)

unmapped_H_count = len(
    unmapped_H_records
)

mapping_distance_max_A = (
    float(
        np.max(BN_distances)
    )
    if len(BN_distances)
    else None
)

mapping_distance_mean_A = (
    float(
        np.mean(BN_distances)
    )
    if len(BN_distances)
    else None
)

mapping_distance_RMS_A = (
    float(
        np.sqrt(
            np.mean(
                BN_distances ** 2
            )
        )
    )
    if len(BN_distances)
    else None
)

print(
    f"mapped_BN_count = "
    f"{mapped_BN_count}/31"
)
print(
    f"unmapped_H_count = "
    f"{unmapped_H_count}/6"
)
print(
    f"unique_topology_mapping = "
    f"{unique_topology_mapping_gate}"
)
print(
    f"mapping_distance_mean_A = "
    f"{mapping_distance_mean_A}"
)
print(
    f"mapping_distance_RMS_A = "
    f"{mapping_distance_RMS_A}"
)
print(
    f"mapping_distance_max_A = "
    f"{mapping_distance_max_A}"
)

for row in mapping_records:
    print(
        f"real_index={row['real_atom_sequence_index']:2d} "
        f"atom_id={row['atom_id']:<28s} "
        f"element={row['element']} "
        f"status={row['mapping_status']:<36s} "
        f"topology_index_1based={row['topology_index_1based']} "
        f"distance_A={row['mapping_distance_A']}"
    )


print("\n[8] SCIENTIFIC INTERPRETATION")

hydrogen_topology_augmentation_required = (
    topology_composition.get(
        "H",
        0,
    )
    == 0
    and composition.get(
        "H",
        0,
    )
    == 6
)

BN_mapping_review_authorized = (
    alignment_authorized
    and mapped_BN_count == 31
    and unique_topology_mapping_gate
    and mapping_distance_max_A
    is not None
    and mapping_distance_max_A
    <= MAX_MAPPING_DISTANCE_A_FOR_REVIEW
)

direct_mapping_execution_authorized = (
    BN_mapping_review_authorized
    and mapping_distance_max_A
    <= MAX_MAPPING_DISTANCE_A_FOR_AUTHORIZATION
    and not hydrogen_topology_augmentation_required
)

print(
    "hydrogen_topology_augmentation_required = "
    f"{hydrogen_topology_augmentation_required}"
)
print(
    "BN_mapping_review_authorized = "
    f"{BN_mapping_review_authorized}"
)
print(
    "direct_mapping_execution_authorized = "
    f"{direct_mapping_execution_authorized}"
)

if hydrogen_topology_augmentation_required:
    print(
        "interpretation = EXISTING_HBN_TOPOLOGY_CANNOT_"
        "REPRESENT_COMPLETE_B17_N14_H6_CHARGE_MODEL"
    )
    print(
        "required_next_action = DESIGN_LOCAL_HYDROGEN_"
        "TOPOLOGY_AUGMENTATION_AND_BONDED_PARAMETERS"
    )
else:
    print(
        "interpretation = COMPLETE_37_ATOM_TARGET_"
        "REPRESENTATION_EXISTS"
    )


print("\n[9] WRITE OUTPUTS")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

mapping_fieldnames = [
    "real_atom_sequence_index",
    "atom_id",
    "element",
    "atom_role",
    "adopted_working_charge_e",
    "mapping_status",
    "topology_index_0based",
    "topology_index_1based",
    "topology_atom_type",
    "topology_atom_name",
    "mapping_distance_A",
]

with MAPPING_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=mapping_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        mapping_records
    )

anchor_fieldnames = [
    "convention",
    "atom_id",
    "element",
    "encoded_index",
    "topology_index_0based",
    "topology_index_1based",
    "mapping_distance_A",
]

with ANCHOR_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=anchor_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        anchor_records_output
    )


print("\n[10] GATES")

gates = {
    "D040_A2_decision_gate": True,
    "adopted_charge_count_37_gate": (
        len(charge_rows) == 37
    ),
    "adopted_composition_gate": (
        composition
        == EXPECTED_COMPOSITION
    ),
    "selected_HBN_atom_count_gate": (
        len(topology_atoms)
        == EXPECTED_HBN_TOPOLOGY_ATOMS
    ),
    "selected_coordinate_coverage_gate": (
        gro["atom_count"]
        >= len(topology_atoms)
    ),
    "selected_HBN_BN_only_gate": (
        topology_composition.get(
            "H",
            0,
        )
        == 0
        and topology_composition.get(
            "B",
            0,
        )
        + topology_composition.get(
            "N",
            0,
        )
        == EXPECTED_HBN_TOPOLOGY_ATOMS
    ),
    "direct_anchor_inventory_gate": (
        len(anchor_candidates) >= 3
    ),
    "mapping_output_created_gate": (
        MAPPING_CSV.is_file()
        and MAPPING_CSV.stat().st_size > 0
    ),
    "anchor_output_created_gate": (
        ANCHOR_CSV.is_file()
        and ANCHOR_CSV.stat().st_size > 0
    ),
    "no_topology_modified_gate": True,
    "no_charge_written_gate": True,
    "force_field_adoption_blocked_gate": True,
    "MD_execution_blocked_gate": True,
}

for gate_name, value in gates.items():
    print(
        f"{gate_name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[11] WRITE REPORT")

decision = (
    "D040_A3_LOCAL_37_ATOM_MAPPING_FEASIBILITY_PASS_"
    "HYDROGEN_TOPOLOGY_AUGMENTATION_REVIEW_AUTHORIZED"
    if (
        all_gates_pass
        and hydrogen_topology_augmentation_required
    )
    else (
        "D040_A3_LOCAL_37_ATOM_MAPPING_FEASIBILITY_PASS_"
        "COMPLETE_MAPPING_REVIEW_AUTHORIZED"
        if all_gates_pass
        else
        "D040_A3_LOCAL_37_ATOM_MAPPING_FEASIBILITY_"
        "REVIEW_REQUIRED"
    )
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_identity": {
        "A2_report": {
            "path": str(
                A2_REPORT.relative_to(ROOT)
            ),
            "sha256": sha256(
                A2_REPORT
            ),
        },
        "adopted_charges_csv": {
            "path": str(
                charges_csv.relative_to(ROOT)
            ),
            "sha256": sha256(
                charges_csv
            ),
        },
        "target_HBN_itp": {
            "path": str(
                TARGET_HBN_ITP.relative_to(ROOT)
            ),
            "sha256": sha256(
                TARGET_HBN_ITP
            ),
        },
        "target_coordinates": {
            "path": str(
                TARGET_COORDINATES.relative_to(ROOT)
            ),
            "sha256": sha256(
                TARGET_COORDINATES
            ),
        },
    },
    "adopted_charge_contract": {
        "atom_count": len(
            charge_rows
        ),
        "composition": (
            composition
        ),
    },
    "selected_HBN_target": {
        "atom_count": len(
            topology_atoms
        ),
        "composition": (
            topology_composition
        ),
        "contains_hydrogen": (
            topology_composition.get(
                "H",
                0,
            )
            > 0
        ),
    },
    "anchor_analysis": {
        "anchor_candidate_count": (
            len(anchor_candidates)
        ),
        "convention_results": [
            {
                key: value
                for key, value in result.items()
                if key
                not in {
                    "rotation",
                    "translation",
                }
            }
            for result in convention_results
        ],
        "selected_convention": (
            best_convention[
                "convention"
            ]
            if best_convention
            is not None
            else None
        ),
        "selected_anchor_RMSD_A": (
            best_convention[
                "anchor_RMSD_A"
            ]
            if best_convention
            is not None
            else None
        ),
    },
    "mapping_diagnostics": {
        "mapped_BN_count": (
            mapped_BN_count
        ),
        "expected_BN_count": 31,
        "unmapped_H_count": (
            unmapped_H_count
        ),
        "expected_H_count": 6,
        "unique_topology_mapping": (
            unique_topology_mapping_gate
        ),
        "distance_mean_A": (
            mapping_distance_mean_A
        ),
        "distance_RMS_A": (
            mapping_distance_RMS_A
        ),
        "distance_max_A": (
            mapping_distance_max_A
        ),
    },
    "scientific_interpretation": {
        "hydrogen_topology_augmentation_required": (
            hydrogen_topology_augmentation_required
        ),
        "BN_mapping_review_authorized": (
            BN_mapping_review_authorized
        ),
        "direct_mapping_execution_authorized": (
            direct_mapping_execution_authorized
        ),
        "complete_existing_topology_representation": (
            not hydrogen_topology_augmentation_required
        ),
    },
    "mapping_records": (
        mapping_records
    ),
    "gates": gates,
    "authorizations": {
        "BN_mapping_candidate_review_authorized": (
            BN_mapping_review_authorized
        ),
        "hydrogen_topology_augmentation_review_authorized": (
            all_gates_pass
            and hydrogen_topology_augmentation_required
        ),
        "charge_to_topology_mapping_execution_authorized": (
            direct_mapping_execution_authorized
        ),
        "topology_modification_authorized": False,
        "new_atom_type_definition_authorized": False,
        "bonded_parameter_modification_authorized": False,
        "force_field_adoption_authorized": False,
        "energy_execution_authorized": False,
        "minimization_execution_authorized": False,
        "validation_MD_execution_authorized": False,
        "production_MD_authorized": False,
    },
    "next_required_block": {
        "name": (
            "D040_A4_LOCAL_HYDROGEN_TOPOLOGY_"
            "AUGMENTATION_DESIGN"
            if hydrogen_topology_augmentation_required
            else
            "D040_A4_MAPPING_CANDIDATE_SCIENTIFIC_AUDIT"
        ),
        "required_actions": (
            [
                (
                    "Determine whether the six retained H atoms "
                    "are physical passivants required in the final "
                    "classical model."
                ),
                (
                    "Identify their B/N parent atoms and required "
                    "B-H or N-H bonds."
                ),
                (
                    "Define atom types, nonbonded parameters, masses, "
                    "bond lengths, angles and any impropers without "
                    "editing the accepted topology."
                ),
                (
                    "Reconcile atom counts and coordinate insertion "
                    "before charge mapping."
                ),
            ]
            if hydrogen_topology_augmentation_required
            else
            [
                (
                    "Audit connectivity and local-environment "
                    "consistency for every proposed mapping."
                ),
                (
                    "Authorize mapping execution only after all "
                    "one-to-one and topology diagnostics pass."
                ),
            ]
        ),
    },
}

REPORT_JSON.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"mapping_csv = {MAPPING_CSV}")
print(
    f"mapping_csv_sha256 = "
    f"{sha256(MAPPING_CSV)}"
)
print(f"anchor_csv = {ANCHOR_CSV}")
print(
    f"anchor_csv_sha256 = "
    f"{sha256(ANCHOR_CSV)}"
)
print(f"report_json = {REPORT_JSON}")
print(
    f"report_json_sha256 = "
    f"{sha256(REPORT_JSON)}"
)


print("\n[12] DECISION")

print(f"decision={decision}")
print(
    "BN_mapping_candidate_review_authorized="
    f"{BN_mapping_review_authorized}"
)
print(
    "hydrogen_topology_augmentation_review_authorized="
    f"{all_gates_pass and hydrogen_topology_augmentation_required}"
)
print(
    "charge_to_topology_mapping_execution_authorized="
    f"{direct_mapping_execution_authorized}"
)
print(
    "topology_modification_authorized=False"
)
print(
    "force_field_adoption_authorized=False"
)
print(
    "validation_MD_execution_authorized=False"
)
print(
    "production_MD_authorized=False"
)
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
