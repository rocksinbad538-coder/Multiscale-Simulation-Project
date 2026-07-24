#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

NODE_FILE = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

EDGE_FILE = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

COORDINATE_FILE = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "23_r2_four_atom_hydrogen_coordinate_embedding/"
    "r2_selected_four_atom_full_coordinates.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_canonical_embedding_conflicts"
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_canonical_close_nonedges.csv"
)

OUTPUT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_CANONICAL_EMBEDDING_CONFLICTS.json"
)

STRICT_THRESHOLD_A = 1.90
DIAGNOSTIC_THRESHOLD_A = 2.10

TARGET_PAIRS = {
    tuple(sorted((
        "S:1739",
        "BR4:UPPER:00:3",
    ))),
    tuple(sorted((
        "BR4:UPPER:14:1",
        "BR4:UPPER:00:4",
    ))),
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def canonical_pair(first: str, second: str) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def distance(
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
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    node_rows = read_csv(NODE_FILE)
    edge_rows = read_csv(EDGE_FILE)
    coordinate_rows = read_csv(COORDINATE_FILE)

    elements = {
        row["node_id"]: row["element"]
        for row in node_rows
    }

    coordinates = {
        row["node_id"]: (
            10.0 * float(row["x_nm"]),
            10.0 * float(row["y_nm"]),
            10.0 * float(row["z_nm"]),
        )
        for row in coordinate_rows
        if row.get("x_nm")
        and row.get("y_nm")
        and row.get("z_nm")
    }

    adjacency: dict[str, set[str]] = defaultdict(set)
    canonical_edges: set[tuple[str, str]] = set()

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        pair = canonical_pair(first, second)
        canonical_edges.add(pair)

        adjacency[first].add(second)
        adjacency[second].add(first)

    atom_ids = sorted(
        atom_id
        for atom_id in coordinates
        if atom_id in elements
    )

    records: list[dict[str, object]] = []

    for first_position, first in enumerate(atom_ids):
        for second in atom_ids[first_position + 1:]:
            first_element = elements[first]
            second_element = elements[second]

            if {first_element, second_element} != {"B", "N"}:
                continue

            pair = canonical_pair(first, second)

            if pair in canonical_edges:
                continue

            value = distance(
                coordinates[first],
                coordinates[second],
            )

            if value > DIAGNOSTIC_THRESHOLD_A:
                continue

            first_h_neighbors = sorted(
                neighbor
                for neighbor in adjacency[first]
                if elements.get(neighbor) == "H"
            )

            second_h_neighbors = sorted(
                neighbor
                for neighbor in adjacency[second]
                if elements.get(neighbor) == "H"
            )

            record = {
                "first_atom": first,
                "first_element": first_element,
                "second_atom": second,
                "second_element": second_element,
                "distance_A": value,
                "strict_conflict": (
                    value <= STRICT_THRESHOLD_A
                ),
                "target_failure_pair": (
                    pair in TARGET_PAIRS
                ),
                "first_degree": len(
                    adjacency[first]
                ),
                "second_degree": len(
                    adjacency[second]
                ),
                "first_neighbors": "|".join(
                    sorted(adjacency[first])
                ),
                "second_neighbors": "|".join(
                    sorted(adjacency[second])
                ),
                "first_H_neighbors": "|".join(
                    first_h_neighbors
                ),
                "second_H_neighbors": "|".join(
                    second_h_neighbors
                ),
            }

            records.append(record)

    records.sort(
        key=lambda row: (
            float(row["distance_A"]),
            str(row["first_atom"]),
            str(row["second_atom"]),
        )
    )

    fieldnames = [
        "rank",
        "first_atom",
        "first_element",
        "second_atom",
        "second_element",
        "distance_A",
        "strict_conflict",
        "target_failure_pair",
        "first_degree",
        "second_degree",
        "first_neighbors",
        "second_neighbors",
        "first_H_neighbors",
        "second_H_neighbors",
    ]

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

        for rank, record in enumerate(records, start=1):
            writer.writerow({
                "rank": rank,
                **record,
            })

    strict_records = [
        record
        for record in records
        if bool(record["strict_conflict"])
    ]

    target_records = [
        record
        for record in records
        if bool(record["target_failure_pair"])
    ]

    report = {
        "decision": (
            "QM_F06_UPPER_CANONICAL_EMBEDDING_"
            "CONFLICTS_DETECTED_TOPOLOGY_OR_"
            "GEOMETRY_REDESIGN_REQUIRED"
            if strict_records
            else
            "QM_F06_UPPER_CANONICAL_EMBEDDING_"
            "NO_STRICT_BN_NONEDGE_CONFLICTS"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "strict_threshold_A": STRICT_THRESHOLD_A,
        "diagnostic_threshold_A": (
            DIAGNOSTIC_THRESHOLD_A
        ),
        "coordinate_count": len(coordinates),
        "canonical_edge_count": len(
            canonical_edges
        ),
        "diagnostic_close_nonedge_count": len(
            records
        ),
        "strict_conflict_count": len(
            strict_records
        ),
        "target_pair_records": target_records,
        "orca_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_JSON.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 110)
    print("QM_F06 UPPER CANONICAL EMBEDDING CONFLICT AUDIT")
    print("=" * 110)
    print("Coordinates:", len(coordinates))
    print("Canonical edges:", len(canonical_edges))
    print(
        "B-N nonedges <= 2.10 A:",
        len(records),
    )
    print(
        "Strict B-N nonedge conflicts <= 1.90 A:",
        len(strict_records),
    )

    print()
    print("Closest B-N nonedges:")

    for rank, record in enumerate(records[:25], start=1):
        marker = (
            "TARGET"
            if record["target_failure_pair"]
            else ""
        )

        print(
            f"{rank:3d} | "
            f"{record['first_atom']:28s} -- "
            f"{record['second_atom']:28s} | "
            f"{float(record['distance_A']):8.5f} A | "
            f"degrees="
            f"{record['first_degree']}/"
            f"{record['second_degree']} "
            f"{marker}"
        )

        if record["target_failure_pair"]:
            print(
                "      first H neighbors:",
                record["first_H_neighbors"] or "none",
            )
            print(
                "      second H neighbors:",
                record["second_H_neighbors"] or "none",
            )

    print()
    print("Decision:", report["decision"])
    print("CSV:", OUTPUT_CSV)
    print("Report:", OUTPUT_JSON)
    print()
    print("ORCA authorized: False")


if __name__ == "__main__":
    main()
