#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil


ROOT = Path(__file__).resolve().parents[2]

POST_QM_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_post_qm"
)

POST_QM_REPORT = (
    POST_QM_DIR
    / "QM_F06_UPPER_V7A_R1_POST_QM_AUDIT.json"
)

POST_QM_FINAL_XYZ = (
    POST_QM_DIR
    / "QM_F06_UPPER_V7A_R1_FINAL.xyz"
)

COORDINATE_CONSISTENCY_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_coordinate_consistency"
)

COORDINATE_CONSISTENCY_REPORT = (
    COORDINATE_CONSISTENCY_DIR
    / "QM_F06_UPPER_V7A_R1_COORDINATE_CONSISTENCY.json"
)

EXEC_PARENT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_EXECUTION_FILE = (
    EXEC_PARENT
    / "LATEST_V7A_R1_EXECUTION.txt"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_coordinate_adoption"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_XYZ = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_ADOPTED_FINAL.xyz"
)

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_ADOPTED_atom_role_provenance_map.csv"
)

OUTPUT_EDGES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_ADOPTED_nominal_edges.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_COORDINATE_ADOPTION.json"
)


EXPECTED_POST_QM_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "POST_QM_GATE_PASS_"
    "RESP_INPUT_PREPARATION_AUTHORIZED"
)

EXPECTED_COORDINATE_CONSISTENCY_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "COORDINATE_CONSISTENCY_PASS_"
    "COORDINATE_ADOPTION_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52

EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_EDGE_COUNT = 57

COORDINATE_IDENTITY_TOLERANCE_A = 5.0e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="strict",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"XYZ file is incomplete: {path}"
        )

    try:
        atom_count = int(
            lines[0].strip()
        )
    except ValueError as error:
        raise RuntimeError(
            f"Invalid XYZ atom count: {path}"
        ) from error

    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ atom record "
                f"at index {index}: {path}"
            )

        try:
            coordinates = tuple(
                float(value)
                for value in fields[1:4]
            )
        except ValueError as error:
            raise RuntimeError(
                f"Invalid XYZ coordinates "
                f"at index {index}: {path}"
            ) from error

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": coordinates,
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"XYZ declared {atom_count} atoms "
            f"but only {len(atoms)} were parsed: {path}"
        )

    return atoms


def read_csv(path: Path) -> list[dict]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def coordinate_difference(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def main() -> None:
    if not POST_QM_REPORT.is_file():
        raise RuntimeError(
            "Final post-QM report is not available: "
            f"{POST_QM_REPORT}"
        )

    post_qm = json.loads(
        POST_QM_REPORT.read_text(
            encoding="utf-8"
        )
    )

    post_qm_decision = post_qm.get(
        "decision"
    )

    post_qm_authorized = (
        post_qm.get(
            "authorizations",
            {},
        ).get(
            "RESP_input_preparation_authorized"
        )
        is True
    )

    if (
        post_qm_decision
        != EXPECTED_POST_QM_DECISION
        or not post_qm_authorized
    ):
        raise RuntimeError(
            "Final coordinate adoption is not authorized "
            "by the post-QM gate: "
            f"decision={post_qm_decision}; "
            f"RESP_input_preparation_authorized="
            f"{post_qm_authorized}"
        )

    if (
        not POST_QM_FINAL_XYZ.is_file()
        or POST_QM_FINAL_XYZ.stat().st_size == 0
    ):
        raise RuntimeError(
            "Post-QM audited final XYZ is missing: "
            f"{POST_QM_FINAL_XYZ}"
        )

    if not COORDINATE_CONSISTENCY_REPORT.is_file():
        raise RuntimeError(
            "Coordinate-consistency report is not available: "
            f"{COORDINATE_CONSISTENCY_REPORT}"
        )

    coordinate_consistency = json.loads(
        COORDINATE_CONSISTENCY_REPORT.read_text(
            encoding="utf-8"
        )
    )

    coordinate_consistency_decision = (
        coordinate_consistency.get(
            "decision"
        )
    )

    coordinate_adoption_authorized = (
        coordinate_consistency.get(
            "authorizations",
            {},
        ).get(
            "coordinate_adoption_authorized"
        )
        is True
    )

    if (
        coordinate_consistency_decision
        != EXPECTED_COORDINATE_CONSISTENCY_DECISION
        or not coordinate_adoption_authorized
    ):
        raise RuntimeError(
            "Final coordinate adoption is not authorized "
            "by the coordinate-consistency gate: "
            f"decision={coordinate_consistency_decision}; "
            f"coordinate_adoption_authorized="
            f"{coordinate_adoption_authorized}"
        )


    if not LATEST_EXECUTION_FILE.is_file():
        raise RuntimeError(
            "Latest R1 execution pointer is missing: "
            f"{LATEST_EXECUTION_FILE}"
        )

    execution_relative = (
        LATEST_EXECUTION_FILE
        .read_text(encoding="utf-8")
        .strip()
    )

    execution_dir = (
        ROOT
        / execution_relative
    )

    source_map = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_constraint_map.csv"
    )

    source_edges = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_nominal_edges.csv"
    )

    source_orca_xyz = (
        execution_dir
        / "v7a_r1.xyz"
    )

    for path in (
        source_map,
        source_edges,
        source_orca_xyz,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing coordinate-adoption source: {path}"
            )

    audited_atoms = read_xyz(
        POST_QM_FINAL_XYZ
    )

    source_atoms = read_xyz(
        source_orca_xyz
    )

    map_rows = read_csv(
        source_map
    )

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    edge_rows = read_csv(
        source_edges
    )

    gates = {}

    gates["post_QM_decision_authorizes_adoption"] = (
        post_qm_decision
        == EXPECTED_POST_QM_DECISION
        and post_qm_authorized
    )

    gates[
        "coordinate_consistency_gate_authorizes_adoption"
    ] = (
        coordinate_consistency_decision
        == EXPECTED_COORDINATE_CONSISTENCY_DECISION
        and coordinate_adoption_authorized
    )

    gates["audited_XYZ_atom_count_52"] = (
        len(audited_atoms)
        == EXPECTED_ATOM_COUNT
    )

    gates["source_ORCA_XYZ_atom_count_52"] = (
        len(source_atoms)
        == EXPECTED_ATOM_COUNT
    )

    gates["map_atom_count_52"] = (
        len(map_rows)
        == EXPECTED_ATOM_COUNT
    )

    gates["nominal_edge_count_57"] = (
        len(edge_rows)
        == EXPECTED_EDGE_COUNT
    )

    audited_composition = Counter(
        atom["element"]
        for atom in audited_atoms
    )

    source_composition = Counter(
        atom["element"]
        for atom in source_atoms
    )

    gates["audited_composition_B17_N14_H21"] = (
        dict(audited_composition)
        == EXPECTED_COMPOSITION
    )

    gates["source_composition_B17_N14_H21"] = (
        dict(source_composition)
        == EXPECTED_COMPOSITION
    )

    identity_failures = []
    coordinate_differences = []

    comparable_count = min(
        len(audited_atoms),
        len(source_atoms),
        len(map_rows),
    )

    for index in range(comparable_count):
        audited_atom = audited_atoms[index]
        source_atom = source_atoms[index]
        map_row = map_rows[index]

        if not (
            audited_atom["element"]
            == source_atom["element"]
            == map_row["element"]
        ):
            identity_failures.append({
                "index_0based": index,
                "atom_id": map_row["atom_id"],
                "audited_element": (
                    audited_atom["element"]
                ),
                "source_element": (
                    source_atom["element"]
                ),
                "map_element": map_row["element"],
            })

        coordinate_differences.append(
            coordinate_difference(
                audited_atom["xyz_A"],
                source_atom["xyz_A"],
            )
        )

    maximum_coordinate_difference = (
        max(coordinate_differences)
        if coordinate_differences
        else float("inf")
    )

    gates["atom_identity_and_order"] = (
        comparable_count
        == EXPECTED_ATOM_COUNT
        and len(identity_failures) == 0
    )

    gates["audited_coordinates_match_ORCA_final"] = (
        maximum_coordinate_difference
        <= COORDINATE_IDENTITY_TOLERANCE_A
    )

    all_atom_ids = [
        row["atom_id"]
        for row in map_rows
    ]

    gates["atom_ids_unique"] = (
        len(all_atom_ids)
        == len(set(all_atom_ids))
        == EXPECTED_ATOM_COUNT
    )

    edge_pairs = [
        tuple(sorted((
            row["first_atom"],
            row["second_atom"],
        )))
        for row in edge_rows
    ]

    gates["nominal_edges_unique"] = (
        len(edge_pairs)
        == len(set(edge_pairs))
        == EXPECTED_EDGE_COUNT
    )

    atom_id_set = set(all_atom_ids)

    gates["all_nominal_edge_atoms_present"] = (
        all(
            first in atom_id_set
            and second in atom_id_set
            for first, second in edge_pairs
        )
    )

    passed = all(gates.values())

    if not passed:
        failed_gates = [
            name
            for name, value in gates.items()
            if not value
        ]

        raise RuntimeError(
            "Final coordinate-adoption gate failed: "
            + " | ".join(failed_gates)
        )

    shutil.copy2(
        POST_QM_FINAL_XYZ,
        OUTPUT_XYZ,
    )

    shutil.copy2(
        source_map,
        OUTPUT_MAP,
    )

    shutil.copy2(
        source_edges,
        OUTPUT_EDGES,
    )

    adopted_atoms = read_xyz(
        OUTPUT_XYZ
    )

    copied_map_rows = read_csv(
        OUTPUT_MAP
    )

    copied_edge_rows = read_csv(
        OUTPUT_EDGES
    )

    copied_coordinate_differences = [
        coordinate_difference(
            adopted_atoms[index]["xyz_A"],
            audited_atoms[index]["xyz_A"],
        )
        for index in range(
            EXPECTED_ATOM_COUNT
        )
    ]

    maximum_copy_coordinate_difference = max(
        copied_coordinate_differences
    )

    post_copy_gates = {
        "adopted_XYZ_hash_matches_source": (
            sha256(OUTPUT_XYZ)
            == sha256(POST_QM_FINAL_XYZ)
        ),
        "adopted_map_hash_matches_source": (
            sha256(OUTPUT_MAP)
            == sha256(source_map)
        ),
        "adopted_edges_hash_matches_source": (
            sha256(OUTPUT_EDGES)
            == sha256(source_edges)
        ),
        "adopted_XYZ_atom_count_52": (
            len(adopted_atoms)
            == EXPECTED_ATOM_COUNT
        ),
        "adopted_map_atom_count_52": (
            len(copied_map_rows)
            == EXPECTED_ATOM_COUNT
        ),
        "adopted_edge_count_57": (
            len(copied_edge_rows)
            == EXPECTED_EDGE_COUNT
        ),
        "adopted_coordinate_copy_exact": (
            maximum_copy_coordinate_difference
            <= COORDINATE_IDENTITY_TOLERANCE_A
        ),
    }

    if not all(post_copy_gates.values()):
        failed_post_copy_gates = [
            name
            for name, value
            in post_copy_gates.items()
            if not value
        ]

        raise RuntimeError(
            "Coordinate copy-integrity gate failed: "
            + " | ".join(
                failed_post_copy_gates
            )
        )

    decision = (
        "QM_F06_UPPER_V7A_"
        "FINAL_COORDINATES_ADOPTED_"
        "RESP_INPUT_DESIGN_AUTHORIZED"
    )

    report = {
        "model": "QM_F06_UPPER_V7A",
        "source_model": "QM_F06_UPPER_V7A_R1",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "decision": decision,
        "source_execution_directory": str(
            execution_dir
        ),
        "post_QM_report": str(
            POST_QM_REPORT
        ),
        "post_QM_decision": (
            post_qm_decision
        ),
        "coordinate_source": str(
            POST_QM_FINAL_XYZ
        ),
        "coordinate_adoption_policy": (
            "BYTE_IDENTICAL_COPY_OF_INDEPENDENTLY_"
            "AUDITED_POST_QM_FINAL_GEOMETRY"
        ),
        "gates": gates,
        "post_copy_gates": (
            post_copy_gates
        ),
        "summary": {
            "atom_count": len(
                adopted_atoms
            ),
            "composition": dict(
                sorted(
                    Counter(
                        atom["element"]
                        for atom in adopted_atoms
                    ).items()
                )
            ),
            "nominal_edge_count": len(
                copied_edge_rows
            ),
            "maximum_audited_vs_ORCA_coordinate_difference_A": (
                maximum_coordinate_difference
            ),
            "maximum_source_vs_adopted_coordinate_difference_A": (
                maximum_copy_coordinate_difference
            ),
            "coordinate_identity_tolerance_A": (
                COORDINATE_IDENTITY_TOLERANCE_A
            ),
        },
        "outputs": {
            "adopted_final_XYZ": str(
                OUTPUT_XYZ
            ),
            "adopted_atom_role_provenance_map": str(
                OUTPUT_MAP
            ),
            "adopted_nominal_edges": str(
                OUTPUT_EDGES
            ),
        },
        "sha256": {
            "post_QM_report": sha256(
                POST_QM_REPORT
            ),
            "source_post_QM_final_XYZ": sha256(
                POST_QM_FINAL_XYZ
            ),
            "source_ORCA_final_XYZ": sha256(
                source_orca_xyz
            ),
            "source_map": sha256(
                source_map
            ),
            "source_edges": sha256(
                source_edges
            ),
            "adopted_final_XYZ": sha256(
                OUTPUT_XYZ
            ),
            "adopted_map": sha256(
                OUTPUT_MAP
            ),
            "adopted_edges": sha256(
                OUTPUT_EDGES
            ),
        },
        "authorizations": {
            "RESP_input_design_authorized": True,
            "RESP_input_preparation_authorized": (
                False
            ),
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": (
                False
            ),
            "MD_authorized": False,
        },

        "upstream_coordinate_consistency_gate": {
            "report": str(
                COORDINATE_CONSISTENCY_REPORT
            ),
            "decision": (
                coordinate_consistency_decision
            ),
            "coordinate_adoption_authorized": (
                coordinate_adoption_authorized
            ),
            "expected_decision": (
                EXPECTED_COORDINATE_CONSISTENCY_DECISION
            ),
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A "
        "FINAL COORDINATE ADOPTION"
    )
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:62s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Post-copy integrity gates:")

    for name, value in post_copy_gates.items():
        print(
            f"{name:62s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(
        "Atoms:",
        len(adopted_atoms),
    )
    print(
        "Composition:",
        report["summary"]["composition"],
    )
    print(
        "Nominal edges:",
        len(copied_edge_rows),
    )
    print(
        "Maximum audited/ORCA difference A:",
        maximum_coordinate_difference,
    )
    print(
        "Maximum source/adopted difference A:",
        maximum_copy_coordinate_difference,
    )

    print()
    print("Decision:", decision)
    print("Adopted XYZ:", OUTPUT_XYZ)
    print("Adopted map:", OUTPUT_MAP)
    print("Adopted edges:", OUTPUT_EDGES)
    print("Report:", OUTPUT_REPORT)

    print()
    print(
        "RESP input design authorized: True"
    )
    print(
        "RESP input preparation authorized: False"
    )
    print("RESP execution authorized: False")
    print(
        "Force-field adoption authorized: False"
    )
    print("MD authorized: False")


if __name__ == "__main__":
    main()
