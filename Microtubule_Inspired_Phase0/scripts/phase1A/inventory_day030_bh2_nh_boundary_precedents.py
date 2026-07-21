#!/usr/bin/env python3
"""
Inventory B-H, B-H2, N-H and N-H2 local motifs in Phase 1A XYZ files.

This is an exploratory precedent inventory. It does not classify a
geometry as accepted or rejected and does not authorize V4 construction.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SEARCH_ROOT = ROOT / "runs/phase1A"

OUTPUT = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_final_closure/"
    "QM_F06_BOUNDARY_MOTIF_PRECEDENT_INVENTORY.csv"
)

BH_CUTOFF_A = 1.35
NH_CUTOFF_A = 1.25
BN_CUTOFF_A = 1.90


def read_xyz(path: Path):
    try:
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        atom_count = int(lines[0].strip())

        if len(lines) < atom_count + 2:
            return None

        atoms = []

        for index, line in enumerate(
            lines[2:2 + atom_count]
        ):
            fields = line.split()

            if len(fields) < 4:
                return None

            atoms.append({
                "index": index,
                "element": fields[0],
                "xyz": tuple(
                    float(value)
                    for value in fields[1:4]
                ),
            })

        return atoms

    except Exception:
        return None


def distance(atom_a, atom_b):
    return math.sqrt(sum(
        (a - b) ** 2
        for a, b in zip(
            atom_a["xyz"],
            atom_b["xyz"],
        )
    ))


def classify_file(path: Path):
    text = str(path).lower()

    if "lower" in text:
        fragment_label = "LOWER"
    elif "upper" in text or "v3a" in text:
        fragment_label = "UPPER"
    else:
        fragment_label = "OTHER"

    if any(
        token in text
        for token in (
            "rejected",
            "failed",
            "interim",
            "restart4",
            "v3a2",
        )
    ):
        provisional_status = "REJECTED_OR_DIAGNOSTIC_CONTEXT"
    elif "accepted" in text or "reference" in text:
        provisional_status = "POSSIBLE_ACCEPTED_CONTEXT"
    else:
        provisional_status = "UNCLASSIFIED"

    return fragment_label, provisional_status


def main():
    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    records = []

    for path in sorted(
        SEARCH_ROOT.rglob("*.xyz")
    ):
        atoms = read_xyz(path)

        if not atoms:
            continue

        fragment_label, provisional_status = (
            classify_file(path)
        )

        for center in atoms:
            if center["element"] not in {"B", "N"}:
                continue

            h_cutoff = (
                BH_CUTOFF_A
                if center["element"] == "B"
                else NH_CUTOFF_A
            )

            h_neighbors = sorted(
                (
                    distance(center, atom),
                    atom["index"],
                )
                for atom in atoms
                if atom["element"] == "H"
                and distance(center, atom) <= h_cutoff
            )

            if not h_neighbors:
                continue

            heavy_neighbors = sorted(
                (
                    distance(center, atom),
                    atom["index"],
                    atom["element"],
                )
                for atom in atoms
                if atom["index"] != center["index"]
                and atom["element"] in {"B", "N"}
                and distance(center, atom) <= BN_CUTOFF_A
            )

            motif = (
                f"{center['element']}-H"
                if len(h_neighbors) == 1
                else f"{center['element']}-H{len(h_neighbors)}"
            )

            records.append({
                "path": str(path.relative_to(ROOT)),
                "fragment_label": fragment_label,
                "provisional_context": provisional_status,
                "center_index_0based": center["index"],
                "center_element": center["element"],
                "motif": motif,
                "H_count": len(h_neighbors),
                "heavy_neighbor_count": len(
                    heavy_neighbors
                ),
                "H_distances_A": "|".join(
                    f"{value:.8f}"
                    for value, _ in h_neighbors
                ),
                "heavy_distances_A": "|".join(
                    f"{value:.8f}"
                    for value, _, _ in heavy_neighbors
                ),
                "heavy_elements": "|".join(
                    element
                    for _, _, element in heavy_neighbors
                ),
            })

    fields = [
        "path",
        "fragment_label",
        "provisional_context",
        "center_index_0based",
        "center_element",
        "motif",
        "H_count",
        "heavy_neighbor_count",
        "H_distances_A",
        "heavy_distances_A",
        "heavy_elements",
    ]

    with OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(records)

    lower_records = [
        row for row in records
        if row["fragment_label"] == "LOWER"
    ]

    print("=" * 78)
    print("BOUNDARY MOTIF PRECEDENT INVENTORY")
    print("=" * 78)
    print("Total XYZ motifs:", len(records))
    print("LOWER motifs:", len(lower_records))

    for label, subset in (
        ("ALL", records),
        ("LOWER", lower_records),
    ):
        print()
        print(label)

        for motif in ("B-H", "B-H2", "N-H", "N-H2"):
            count = sum(
                row["motif"] == motif
                for row in subset
            )
            print(f"{motif:8s}: {count}")

    print()
    print("LOWER B-H2 / N-H examples:")

    for row in lower_records:
        if row["motif"] in {"B-H2", "N-H"}:
            print(
                f"{row['motif']:5s} | "
                f"{row['path']} | "
                f"center={row['center_index_0based']} | "
                f"heavy={row['heavy_neighbor_count']} | "
                f"H={row['H_distances_A']}"
            )

    print()
    print("Inventory:", OUTPUT)
    print("Chemical precedent accepted: False")
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
