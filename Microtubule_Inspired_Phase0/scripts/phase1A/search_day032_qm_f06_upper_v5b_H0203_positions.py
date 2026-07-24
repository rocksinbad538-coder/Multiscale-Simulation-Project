#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

START_XYZ = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_start.xyz"
)

MAP_CSV = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5b_H0203_position_search"
)

RANKING_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V5B_H0203_position_ranking.csv"
)

REPORT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V5B_H0203_POSITION_SEARCH.json"
)

TARGET_H = "H4:UPPER:0203:0"
OWNER = "S:1739"

TARGET_BH_A = 1.19

N_THETA = 120
N_PHI = 240

BH_MAX_A = 1.35
NH_MAX_A = 1.25

HH_HARD_A = 1.20
HX_HARD_A = 0.85


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

    return atoms


def distance(first, second):
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def add(first, second):
    return tuple(
        a + b
        for a, b in zip(first, second)
    )


def scale(vector, factor):
    return tuple(
        factor * value
        for value in vector
    )


def normalize(vector):
    magnitude = math.sqrt(
        sum(value * value for value in vector)
    )

    if magnitude == 0.0:
        raise RuntimeError("Zero-length vector.")

    return tuple(
        value / magnitude
        for value in vector
    )


def main():
    mapping = read_csv(MAP_CSV)
    atoms = read_xyz(START_XYZ)

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    id_by_index = {
        int(row["index_0based"]): row["atom_id"]
        for row in mapping
    }

    target_index = index_by_id[TARGET_H]
    owner_index = index_by_id[OWNER]

    owner_xyz = atoms[owner_index]["xyz"]

    candidates = []

    golden_angle = math.pi * (
        3.0 - math.sqrt(5.0)
    )

    sample_count = N_THETA * N_PHI

    # Fibonacci sphere: uniform direction sampling.
    for sample_index in range(sample_count):
        z = (
            1.0
            - 2.0
            * (sample_index + 0.5)
            / sample_count
        )

        radius_xy = math.sqrt(
            max(0.0, 1.0 - z * z)
        )

        phi = sample_index * golden_angle

        direction = (
            radius_xy * math.cos(phi),
            radius_xy * math.sin(phi),
            z,
        )

        candidate_xyz = add(
            owner_xyz,
            scale(
                normalize(direction),
                TARGET_BH_A,
            ),
        )

        extra_bond_owners = []
        hard_contacts = []
        nearest_nonowner = None

        for index, atom in enumerate(atoms):
            if index in {
                target_index,
                owner_index,
            }:
                continue

            value = distance(
                candidate_xyz,
                atom["xyz"],
            )

            if (
                nearest_nonowner is None
                or value < nearest_nonowner["distance_A"]
            ):
                nearest_nonowner = {
                    "atom_id": id_by_index[index],
                    "element": atom["element"],
                    "distance_A": value,
                }

            if atom["element"] == "B":
                if value <= BH_MAX_A:
                    extra_bond_owners.append(
                        id_by_index[index]
                    )

                if value < HX_HARD_A:
                    hard_contacts.append(
                        id_by_index[index]
                    )

            elif atom["element"] == "N":
                if value <= NH_MAX_A:
                    extra_bond_owners.append(
                        id_by_index[index]
                    )

                if value < HX_HARD_A:
                    hard_contacts.append(
                        id_by_index[index]
                    )

            elif atom["element"] == "H":
                if value < HH_HARD_A:
                    hard_contacts.append(
                        id_by_index[index]
                    )

        multi_owner = bool(extra_bond_owners)

        valid = (
            not multi_owner
            and not hard_contacts
            and nearest_nonowner is not None
        )

        clearance = (
            nearest_nonowner["distance_A"]
            if nearest_nonowner
            else 0.0
        )

        candidates.append({
            "sample_index": sample_index,
            "x_A": candidate_xyz[0],
            "y_A": candidate_xyz[1],
            "z_A": candidate_xyz[2],
            "owner_distance_A": distance(
                candidate_xyz,
                owner_xyz,
            ),
            "nearest_nonowner_atom": (
                nearest_nonowner["atom_id"]
                if nearest_nonowner
                else ""
            ),
            "nearest_nonowner_element": (
                nearest_nonowner["element"]
                if nearest_nonowner
                else ""
            ),
            "nearest_nonowner_distance_A": (
                clearance
            ),
            "extra_geometric_owners": "|".join(
                sorted(extra_bond_owners)
            ),
            "hard_contacts": "|".join(
                sorted(hard_contacts)
            ),
            "valid_single_owner_position": valid,
            "score": clearance if valid else -1.0,
        })

    candidates.sort(
        key=lambda row: (
            row["valid_single_owner_position"],
            row["score"],
        ),
        reverse=True,
    )

    for rank, row in enumerate(
        candidates,
        start=1,
    ):
        row["rank"] = rank

    valid_candidates = [
        row
        for row in candidates
        if row["valid_single_owner_position"]
    ]

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "rank",
        "sample_index",
        "x_A",
        "y_A",
        "z_A",
        "owner_distance_A",
        "nearest_nonowner_atom",
        "nearest_nonowner_element",
        "nearest_nonowner_distance_A",
        "extra_geometric_owners",
        "hard_contacts",
        "valid_single_owner_position",
        "score",
    ]

    with RANKING_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(candidates[:500])

    report = {
        "decision": (
            "QM_F06_UPPER_V5B_H0203_"
            "VALID_SINGLE_OWNER_POSITION_FOUND"
            if valid_candidates
            else
            "QM_F06_UPPER_V5B_H0203_"
            "NO_VALID_SINGLE_OWNER_POSITION"
        ),
        "sample_count": sample_count,
        "valid_candidate_count": len(
            valid_candidates
        ),
        "best_candidate": (
            valid_candidates[0]
            if valid_candidates
            else None
        ),
        "orca_authorized": False,
        "RESP_authorized": False,
    }

    REPORT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 96)
    print("QM_F06 UPPER V5-B H0203 POSITION SEARCH")
    print("=" * 96)
    print("Samples:", sample_count)
    print(
        "Valid single-owner positions:",
        len(valid_candidates),
    )

    if valid_candidates:
        best = valid_candidates[0]

        print()
        print("Best candidate:")
        for key, value in best.items():
            print(f"  {key}: {value}")

    print()
    print("Decision:", report["decision"])
    print("Ranking:", RANKING_CSV)
    print("Report:", REPORT_JSON)


if __name__ == "__main__":
    main()
