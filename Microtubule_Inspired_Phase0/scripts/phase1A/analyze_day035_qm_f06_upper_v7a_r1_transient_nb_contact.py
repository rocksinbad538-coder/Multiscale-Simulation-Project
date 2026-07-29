#!/usr/bin/env python3

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import csv
import importlib.util
import json
import sys


ROOT = Path(__file__).resolve().parents[2]

AUDITOR_PATH = (
    ROOT
    / "scripts"
    / "phase1A"
    / "audit_day035_qm_f06_upper_v7a_r1_live_trajectory.py"
)

EXEC_PARENT = (
    ROOT
    / "runs"
    / "phase1A"
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_FILE = (
    EXEC_PARENT
    / "LATEST_V7A_R1_EXECUTION.txt"
)

OUTPUT_DIR = (
    ROOT
    / "runs"
    / "phase1A"
    / "day035_qm_f06_upper_v7a_r1_transient_nb_contact"
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_transient_NB_contact.csv"
)

OUTPUT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_TRANSIENT_NB_CONTACT.json"
)

FIRST_ATOM = "A:UPPER:13:1"
SECOND_ATOM = "S:1710"

EXPECTED_FIRST_ELEMENT = "N"
EXPECTED_SECOND_ELEMENT = "B"

GEOMETRIC_BOND_THRESHOLD_A = 1.90


def load_auditor_module():
    if not AUDITOR_PATH.is_file():
        raise RuntimeError(
            f"Missing auditor module: {AUDITOR_PATH}"
        )

    specification = importlib.util.spec_from_file_location(
        "qm_f06_v7a_r1_live_auditor",
        AUDITOR_PATH,
    )

    if specification is None or specification.loader is None:
        raise RuntimeError(
            "Could not create auditor import specification"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    sys.modules[specification.name] = module
    specification.loader.exec_module(module)

    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    auditor = load_auditor_module()

    if not LATEST_FILE.is_file():
        raise RuntimeError(
            f"Missing execution pointer: {LATEST_FILE}"
        )

    execution_relative = (
        LATEST_FILE.read_text(
            encoding="utf-8"
        ).strip()
    )

    if not execution_relative:
        raise RuntimeError(
            "Latest execution pointer is empty"
        )

    execution_dir = ROOT / execution_relative

    trajectory_path = (
        execution_dir
        / "v7a_r1_trj.xyz"
    )

    map_path = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_constraint_map.csv"
    )

    for path in (
        trajectory_path,
        map_path,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing required input: {path}"
            )

    map_rows = read_csv(map_path)

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    atom_ids = [
        row["atom_id"]
        for row in map_rows
    ]

    elements = {
        row["atom_id"]: row["element"]
        for row in map_rows
    }

    missing_atoms = [
        atom_id
        for atom_id in (
            FIRST_ATOM,
            SECOND_ATOM,
        )
        if atom_id not in elements
    ]

    if missing_atoms:
        raise RuntimeError(
            "Missing target atom identities: "
            + ", ".join(missing_atoms)
        )

    if elements[FIRST_ATOM] != EXPECTED_FIRST_ELEMENT:
        raise RuntimeError(
            f"{FIRST_ATOM} has element "
            f"{elements[FIRST_ATOM]}, expected "
            f"{EXPECTED_FIRST_ELEMENT}"
        )

    if elements[SECOND_ATOM] != EXPECTED_SECOND_ELEMENT:
        raise RuntimeError(
            f"{SECOND_ATOM} has element "
            f"{elements[SECOND_ATOM]}, expected "
            f"{EXPECTED_SECOND_ELEMENT}"
        )

    first_index = atom_ids.index(FIRST_ATOM)
    second_index = atom_ids.index(SECOND_ATOM)

    frames = auditor.read_trajectory(
        trajectory_path
    )

    if not frames:
        raise RuntimeError(
            "No complete trajectory frames available"
        )

    records = []

    for frame in frames:
        atoms = frame["atoms"]

        if len(atoms) != len(atom_ids):
            raise RuntimeError(
                "Complete frame atom-count mismatch at "
                f"frame {frame['frame_index_0based']}: "
                f"{len(atoms)} versus {len(atom_ids)}"
            )

        if atoms[first_index]["element"] != EXPECTED_FIRST_ELEMENT:
            raise RuntimeError(
                "Element-order mismatch for first atom at "
                f"frame {frame['frame_index_0based']}"
            )

        if atoms[second_index]["element"] != EXPECTED_SECOND_ELEMENT:
            raise RuntimeError(
                "Element-order mismatch for second atom at "
                f"frame {frame['frame_index_0based']}"
            )

        first_xyz = atoms[first_index]["xyz_A"]
        second_xyz = atoms[second_index]["xyz_A"]

        separation = auditor.distance(
            first_xyz,
            second_xyz,
        )

        margin = (
            separation
            - GEOMETRIC_BOND_THRESHOLD_A
        )

        contact_classification = (
            "GEOMETRIC_EDGE_PRESENT"
            if separation
            <= GEOMETRIC_BOND_THRESHOLD_A
            else "GEOMETRIC_EDGE_ABSENT"
        )

        energy = auditor.extract_energy(
            frame["comment"]
        )

        records.append({
            "frame_index_0based": (
                frame["frame_index_0based"]
            ),
            "energy_Eh": energy,
            "first_atom": FIRST_ATOM,
            "first_element": EXPECTED_FIRST_ELEMENT,
            "first_index_0based": first_index,
            "second_atom": SECOND_ATOM,
            "second_element": EXPECTED_SECOND_ELEMENT,
            "second_index_0based": second_index,
            "distance_A": separation,
            "geometric_bond_threshold_A": (
                GEOMETRIC_BOND_THRESHOLD_A
            ),
            "margin_to_threshold_A": margin,
            "contact_classification": (
                contact_classification
            ),
        })

    present_records = [
        record
        for record in records
        if record["contact_classification"]
        == "GEOMETRIC_EDGE_PRESENT"
    ]

    minimum_record = min(
        records,
        key=lambda record: record["distance_A"],
    )

    first_present_frame = (
        present_records[0]["frame_index_0based"]
        if present_records
        else None
    )

    last_present_frame = (
        present_records[-1]["frame_index_0based"]
        if present_records
        else None
    )

    latest_record = records[-1]

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = list(records[0].keys())

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(records)

    report = {
        "generated_at_utc": (
            datetime.now(timezone.utc)
            .isoformat()
        ),
        "system": "QM_F06_UPPER_V7A_R1",
        "execution_directory": str(
            execution_dir
        ),
        "trajectory": str(
            trajectory_path
        ),
        "contact": {
            "first_atom": FIRST_ATOM,
            "first_element": (
                EXPECTED_FIRST_ELEMENT
            ),
            "first_index_0based": (
                first_index
            ),
            "second_atom": SECOND_ATOM,
            "second_element": (
                EXPECTED_SECOND_ELEMENT
            ),
            "second_index_0based": (
                second_index
            ),
            "pair_class": "B-N",
            "geometric_bond_threshold_A": (
                GEOMETRIC_BOND_THRESHOLD_A
            ),
        },
        "complete_frame_count": len(records),
        "geometric_edge_present_frame_count": (
            len(present_records)
        ),
        "geometric_edge_present_frames": [
            record["frame_index_0based"]
            for record in present_records
        ],
        "first_present_frame": (
            first_present_frame
        ),
        "last_present_frame": (
            last_present_frame
        ),
        "minimum_distance_A": (
            minimum_record["distance_A"]
        ),
        "minimum_distance_frame": (
            minimum_record["frame_index_0based"]
        ),
        "minimum_margin_to_threshold_A": (
            minimum_record[
                "margin_to_threshold_A"
            ]
        ),
        "latest_frame": (
            latest_record["frame_index_0based"]
        ),
        "latest_energy_Eh": (
            latest_record["energy_Eh"]
        ),
        "latest_distance_A": (
            latest_record["distance_A"]
        ),
        "latest_margin_to_threshold_A": (
            latest_record[
                "margin_to_threshold_A"
            ]
        ),
        "latest_contact_classification": (
            latest_record[
                "contact_classification"
            ]
        ),
        "contact_absent_in_latest_frame": (
            latest_record[
                "contact_classification"
            ]
            == "GEOMETRIC_EDGE_ABSENT"
        ),
        "interpretation": (
            "The N-B pair crossed the geometric-edge "
            "threshold transiently during the optimization. "
            "This report characterizes the distance history "
            "but does not authorize post-QM acceptance, "
            "RESP, or MD."
        ),
        "post_QM_acceptance_authorized": False,
        "RESP_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_JSON.write_text(
        json.dumps(
            report,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "TRANSIENT N-B CONTACT ANALYSIS"
    )
    print("=" * 112)

    print(
        f"Contact: {FIRST_ATOM} "
        f"({EXPECTED_FIRST_ELEMENT}) -- "
        f"{SECOND_ATOM} "
        f"({EXPECTED_SECOND_ELEMENT})"
    )
    print(
        "Geometric threshold A: "
        f"{GEOMETRIC_BOND_THRESHOLD_A:.6f}"
    )
    print(
        "Complete frames: "
        f"{len(records)}"
    )
    print()

    for record in records:
        print(
            f"frame={record['frame_index_0based']:3d} | "
            f"E={record['energy_Eh']!s:>20} | "
            f"d_NB={record['distance_A']:.9f} A | "
            f"margin={record['margin_to_threshold_A']:+.9f} A | "
            f"{record['contact_classification']}"
        )

    print()
    print(
        "Geometric-edge-present frames: "
        f"{report['geometric_edge_present_frames']}"
    )
    print(
        "Minimum distance: "
        f"{report['minimum_distance_A']:.9f} A "
        f"at frame "
        f"{report['minimum_distance_frame']}"
    )
    print(
        "Minimum margin: "
        f"{report['minimum_margin_to_threshold_A']:+.9f} A"
    )
    print(
        "Latest frame: "
        f"{report['latest_frame']}"
    )
    print(
        "Latest distance: "
        f"{report['latest_distance_A']:.9f} A"
    )
    print(
        "Latest margin: "
        f"{report['latest_margin_to_threshold_A']:+.9f} A"
    )
    print(
        "Contact absent in latest frame: "
        f"{report['contact_absent_in_latest_frame']}"
    )
    print()
    print(f"CSV: {OUTPUT_CSV}")
    print(f"Report: {OUTPUT_JSON}")
    print()
    print("Post-QM acceptance authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
