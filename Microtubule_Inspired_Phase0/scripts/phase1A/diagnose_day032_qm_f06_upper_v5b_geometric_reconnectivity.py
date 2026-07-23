#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

LATEST_POINTER = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_executions/"
    "LATEST_V5B_EXECUTION.txt"
)

PROVENANCE_MAP = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

PRE_QM_VALENCE = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit/"
    "QM_F06_UPPER_V5B_valence.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5b_geometric_reconnectivity"
)

CSV_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_geometric_reconnectivity.csv"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_GEOMETRIC_RECONNECTIVITY.json"
)


def read_csv(path: Path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0])
    atoms = []

    for line in lines[2:2 + count]:
        fields = line.split()

        atoms.append({
            "element": fields[0],
            "xyz": tuple(map(float, fields[1:4])),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    return atoms


def distance(first, second):
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def geometric_bond(
    first_element: str,
    second_element: str,
    value: float,
) -> bool:
    pair = frozenset((
        first_element,
        second_element,
    ))

    if pair == frozenset(("B", "N")):
        return 1.25 <= value <= 1.85

    if pair == frozenset(("B", "H")):
        return 0.95 <= value <= 1.35

    if pair == frozenset(("N", "H")):
        return 0.85 <= value <= 1.25

    return False


def main():
    execution_relative = (
        LATEST_POINTER
        .read_text(encoding="utf-8")
        .strip()
    )

    execution_dir = ROOT / execution_relative
    current_xyz = execution_dir / "v5b.xyz"

    mapping_rows = read_csv(PROVENANCE_MAP)
    valence_rows = read_csv(PRE_QM_VALENCE)
    atoms = read_xyz(current_xyz)

    if len(mapping_rows) != len(atoms):
        raise RuntimeError(
            "Map/XYZ atom-count mismatch."
        )

    id_by_index = {
        int(row["index_0based"]): row["atom_id"]
        for row in mapping_rows
    }

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping_rows
    }

    nominal_neighbors = {
        row["atom_id"]: set(
            value
            for value in row["neighbors"].split("|")
            if value
        )
        for row in valence_rows
    }

    records = []
    changed = []

    for atom_id, atom_index in index_by_id.items():
        atom = atoms[atom_index]

        geometric_neighbors = set()
        distances = {}

        for other_id, other_index in index_by_id.items():
            if other_id == atom_id:
                continue

            value = distance(
                atom["xyz"],
                atoms[other_index]["xyz"],
            )

            if geometric_bond(
                atom["element"],
                atoms[other_index]["element"],
                value,
            ):
                geometric_neighbors.add(other_id)
                distances[other_id] = value

        nominal = nominal_neighbors.get(
            atom_id,
            set(),
        )

        lost = sorted(
            nominal - geometric_neighbors
        )

        gained = sorted(
            geometric_neighbors - nominal
        )

        changed_flag = bool(
            lost or gained
        )

        record = {
            "atom_id": atom_id,
            "element": atom["element"],
            "nominal_degree": len(nominal),
            "geometric_degree": len(
                geometric_neighbors
            ),
            "nominal_neighbors": "|".join(
                sorted(nominal)
            ),
            "geometric_neighbors": "|".join(
                sorted(geometric_neighbors)
            ),
            "lost_nominal_neighbors": "|".join(
                lost
            ),
            "gained_geometric_neighbors": "|".join(
                gained
            ),
            "geometric_neighbor_distances_A": "|".join(
                f"{neighbor}:{distances[neighbor]:.6f}"
                for neighbor in sorted(distances)
            ),
            "reconnectivity_detected": changed_flag,
        }

        records.append(record)

        if changed_flag:
            changed.append(record)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with CSV_PATH.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0]),
        )
        writer.writeheader()
        writer.writerows(records)

    report = {
        "decision": (
            "QM_F06_UPPER_V5B_GEOMETRIC_"
            "RECONNECTIVITY_DETECTED_REVIEW_REQUIRED"
            if changed
            else
            "QM_F06_UPPER_V5B_NO_GEOMETRIC_"
            "RECONNECTIVITY_DETECTED"
        ),
        "geometry_source": str(
            current_xyz.relative_to(ROOT)
        ),
        "atom_count": len(atoms),
        "changed_atom_count": len(changed),
        "changed_atoms": changed,
        "RESP_authorized": False,
        "MD_authorized": False,
    }

    REPORT_PATH.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 88)
    print("QM_F06 UPPER V5-B GEOMETRIC RECONNECTIVITY")
    print("=" * 88)
    print("Geometry:", current_xyz)
    print("Changed atoms:", len(changed))

    for record in changed:
        print()
        print("Atom:", record["atom_id"])
        print(
            "  Nominal:",
            record["nominal_neighbors"],
        )
        print(
            "  Geometric:",
            record["geometric_neighbors"],
        )
        print(
            "  Lost:",
            record["lost_nominal_neighbors"],
        )
        print(
            "  Gained:",
            record[
                "gained_geometric_neighbors"
            ],
        )

    print()
    print("Decision:", report["decision"])
    print("CSV:", CSV_PATH)
    print("Report:", REPORT_PATH)
    print()
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
