#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SOURCE_XYZ = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_start.xyz"
)

MAP_CSV = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

VALENCE_CSV = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit/"
    "QM_F06_UPPER_V5B_valence.csv"
)

H_RANKING = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5b_H0203_position_search/"
    "QM_F06_UPPER_V5B_H0203_position_ranking.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_local_geometry_repair"
)

RANKING_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_local_geometry_repair_ranking.csv"
)

BEST_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_LOCAL_GEOMETRY_REPAIR_BEST.xyz"
)

REPORT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_LOCAL_GEOMETRY_REPAIR.json"
)


MOVABLE_HEAVY = [
    "S:1739",
    "BR4:UPPER:14:1",
]

MOVABLE_H = "H4:UPPER:0203:0"

NONCANONICAL_PAIRS = [
    ("S:1739", "BR4:UPPER:00:3"),
    ("BR4:UPPER:14:1", "BR4:UPPER:00:4"),
]

BN_MIN_A = 1.25
BN_MAX_A = 1.85

BH_MIN_A = 0.95
BH_MAX_A = 1.35

NH_MIN_A = 0.85
NH_MAX_A = 1.25

NONCANONICAL_BN_MIN_A = 1.95

HH_HARD_A = 1.20
HX_HARD_A = 0.85
HEAVY_HARD_A = 1.20

SAMPLES = 250000
MAX_HEAVY_SHIFT_A = 0.55
RANDOM_SEED = 20260723


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


def add(first, second):
    return tuple(
        a + b
        for a, b in zip(first, second)
    )


def subtract(first, second):
    return tuple(
        a - b
        for a, b in zip(first, second)
    )


def scale(vector, factor):
    return tuple(
        factor * value
        for value in vector
    )


def norm(vector):
    return math.sqrt(
        sum(value * value for value in vector)
    )


def normalize(vector):
    magnitude = norm(vector)

    if magnitude == 0.0:
        raise RuntimeError("Zero vector.")

    return scale(vector, 1.0 / magnitude)


def canonical_pair(first, second):
    return tuple(sorted((first, second)))


def random_ball_vector(rng, maximum_radius):
    while True:
        vector = (
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
        )

        squared = sum(value * value for value in vector)

        if 0.0 < squared <= 1.0:
            radius_factor = rng.random() ** (1.0 / 3.0)

            return scale(
                normalize(vector),
                maximum_radius * radius_factor,
            )


def bond_range(first_element, second_element):
    pair = frozenset((
        first_element,
        second_element,
    ))

    if pair == frozenset(("B", "N")):
        return BN_MIN_A, BN_MAX_A

    if pair == frozenset(("B", "H")):
        return BH_MIN_A, BH_MAX_A

    if pair == frozenset(("N", "H")):
        return NH_MIN_A, NH_MAX_A

    return None


def main():
    rng = random.Random(RANDOM_SEED)

    mapping = read_csv(MAP_CSV)
    valence = read_csv(VALENCE_CSV)
    atoms = read_xyz(SOURCE_XYZ)
    h_candidates = read_csv(H_RANKING)

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    id_by_index = {
        int(row["index_0based"]): row["atom_id"]
        for row in mapping
    }

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in mapping
    }

    nominal_neighbors = {
        row["atom_id"]: [
            value
            for value in row["neighbors"].split("|")
            if value
        ]
        for row in valence
    }

    nominal_pairs = set()

    for atom_id, neighbors in nominal_neighbors.items():
        for neighbor in neighbors:
            nominal_pairs.add(
                canonical_pair(atom_id, neighbor)
            )

    base_xyz = {
        atom_id: atoms[index]["xyz"]
        for atom_id, index in index_by_id.items()
    }

    valid_h_candidates = [
        row
        for row in h_candidates
        if row[
            "valid_single_owner_position"
        ].strip().lower() == "true"
    ][:300]

    if not valid_h_candidates:
        raise RuntimeError(
            "No valid H0203 positions available."
        )

    records = []

    for sample_index in range(SAMPLES):
        coordinates = dict(base_xyz)

        s_shift = random_ball_vector(
            rng,
            MAX_HEAVY_SHIFT_A,
        )

        br_shift = random_ball_vector(
            rng,
            MAX_HEAVY_SHIFT_A,
        )

        coordinates["S:1739"] = add(
            base_xyz["S:1739"],
            s_shift,
        )

        coordinates["BR4:UPPER:14:1"] = add(
            base_xyz["BR4:UPPER:14:1"],
            br_shift,
        )

        if sample_index < len(valid_h_candidates):
            source_h = valid_h_candidates[sample_index]

            direction = subtract(
                (
                    float(source_h["x_A"]),
                    float(source_h["y_A"]),
                    float(source_h["z_A"]),
                ),
                base_xyz["S:1739"],
            )
        else:
            direction = random_ball_vector(
                rng,
                1.0,
            )

        direction = normalize(direction)

        coordinates[MOVABLE_H] = add(
            coordinates["S:1739"],
            scale(direction, 1.19),
        )

        nominal_bonds_pass = True
        maximum_nominal_deviation = 0.0

        relevant_atoms = set(
            MOVABLE_HEAVY
            + [MOVABLE_H]
        )

        for atom_id in list(relevant_atoms):
            relevant_atoms.update(
                nominal_neighbors.get(atom_id, [])
            )

        for first, second in nominal_pairs:
            if (
                first not in relevant_atoms
                and second not in relevant_atoms
            ):
                continue

            value = distance(
                coordinates[first],
                coordinates[second],
            )

            limits = bond_range(
                element_by_id[first],
                element_by_id[second],
            )

            if limits is None:
                continue

            minimum, maximum = limits

            if not minimum <= value <= maximum:
                nominal_bonds_pass = False
                break

            midpoint = 0.5 * (
                minimum + maximum
            )

            maximum_nominal_deviation = max(
                maximum_nominal_deviation,
                abs(value - midpoint),
            )

        if not nominal_bonds_pass:
            continue

        noncanonical_distances = {}

        for first, second in NONCANONICAL_PAIRS:
            noncanonical_distances[
                canonical_pair(first, second)
            ] = distance(
                coordinates[first],
                coordinates[second],
            )

        noncanonical_clearance_pass = all(
            value >= NONCANONICAL_BN_MIN_A
            for value in noncanonical_distances.values()
        )

        if not noncanonical_clearance_pass:
            continue

        h_extra_owners = []
        hard_contacts = []

        h_xyz = coordinates[MOVABLE_H]

        for atom_id, atom_xyz in coordinates.items():
            if atom_id in {
                MOVABLE_H,
                "S:1739",
            }:
                continue

            value = distance(
                h_xyz,
                atom_xyz,
            )

            element = element_by_id[atom_id]

            if element == "B":
                if value <= BH_MAX_A:
                    h_extra_owners.append(atom_id)

                if value < HX_HARD_A:
                    hard_contacts.append(atom_id)

            elif element == "N":
                if value <= NH_MAX_A:
                    h_extra_owners.append(atom_id)

                if value < HX_HARD_A:
                    hard_contacts.append(atom_id)

            elif element == "H":
                if value < HH_HARD_A:
                    hard_contacts.append(atom_id)

        if h_extra_owners or hard_contacts:
            continue

        heavy_hard_contact = False

        for moved_atom in MOVABLE_HEAVY:
            for atom_id, atom_xyz in coordinates.items():
                if atom_id == moved_atom:
                    continue

                pair = canonical_pair(
                    moved_atom,
                    atom_id,
                )

                if pair in nominal_pairs:
                    continue

                value = distance(
                    coordinates[moved_atom],
                    atom_xyz,
                )

                if (
                    element_by_id[atom_id] != "H"
                    and value < HEAVY_HARD_A
                ):
                    heavy_hard_contact = True
                    break

            if heavy_hard_contact:
                break

        if heavy_hard_contact:
            continue

        maximum_heavy_shift = max(
            norm(s_shift),
            norm(br_shift),
        )

        minimum_noncanonical_distance = min(
            noncanonical_distances.values()
        )

        score = (
            100.0 * minimum_noncanonical_distance
            - 20.0 * maximum_heavy_shift
            - 10.0 * maximum_nominal_deviation
        )

        records.append({
            "sample_index": sample_index,
            "score": score,
            "maximum_heavy_shift_A": (
                maximum_heavy_shift
            ),
            "maximum_nominal_bond_midpoint_deviation_A": (
                maximum_nominal_deviation
            ),
            "S1739_BR4_00_3_distance_A": (
                noncanonical_distances[
                    canonical_pair(
                        "S:1739",
                        "BR4:UPPER:00:3",
                    )
                ]
            ),
            "BR4_14_1_BR4_00_4_distance_A": (
                noncanonical_distances[
                    canonical_pair(
                        "BR4:UPPER:14:1",
                        "BR4:UPPER:00:4",
                    )
                ]
            ),
            "S1739_x_A": coordinates["S:1739"][0],
            "S1739_y_A": coordinates["S:1739"][1],
            "S1739_z_A": coordinates["S:1739"][2],
            "BR4_14_1_x_A": coordinates[
                "BR4:UPPER:14:1"
            ][0],
            "BR4_14_1_y_A": coordinates[
                "BR4:UPPER:14:1"
            ][1],
            "BR4_14_1_z_A": coordinates[
                "BR4:UPPER:14:1"
            ][2],
            "H0203_x_A": coordinates[MOVABLE_H][0],
            "H0203_y_A": coordinates[MOVABLE_H][1],
            "H0203_z_A": coordinates[MOVABLE_H][2],
            "valid": True,
        })

    records.sort(
        key=lambda row: row["score"],
        reverse=True,
    )

    for rank, row in enumerate(records, start=1):
        row["rank"] = rank

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    if records:
        fieldnames = [
            "rank",
            "sample_index",
            "score",
            "maximum_heavy_shift_A",
            "maximum_nominal_bond_midpoint_deviation_A",
            "S1739_BR4_00_3_distance_A",
            "BR4_14_1_BR4_00_4_distance_A",
            "S1739_x_A",
            "S1739_y_A",
            "S1739_z_A",
            "BR4_14_1_x_A",
            "BR4_14_1_y_A",
            "BR4_14_1_z_A",
            "H0203_x_A",
            "H0203_y_A",
            "H0203_z_A",
            "valid",
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
            writer.writerows(records[:1000])

        best = records[0]
        repaired_atoms = [
            dict(atom)
            for atom in atoms
        ]

        replacements = {
            "S:1739": (
                best["S1739_x_A"],
                best["S1739_y_A"],
                best["S1739_z_A"],
            ),
            "BR4:UPPER:14:1": (
                best["BR4_14_1_x_A"],
                best["BR4_14_1_y_A"],
                best["BR4_14_1_z_A"],
            ),
            MOVABLE_H: (
                best["H0203_x_A"],
                best["H0203_y_A"],
                best["H0203_z_A"],
            ),
        }

        for atom_id, xyz in replacements.items():
            repaired_atoms[
                index_by_id[atom_id]
            ]["xyz"] = xyz

        with BEST_XYZ.open(
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(
                f"{len(repaired_atoms)}\n"
            )
            handle.write(
                "QM_F06 UPPER V5-C local geometry "
                "repair candidate; ORCA not authorized\n"
            )

            for atom in repaired_atoms:
                x_value, y_value, z_value = atom["xyz"]

                handle.write(
                    f"{atom['element']:2s} "
                    f"{x_value: .12f} "
                    f"{y_value: .12f} "
                    f"{z_value: .12f}\n"
                )
    else:
        best = None

    report = {
        "decision": (
            "QM_F06_UPPER_V5C_LOCAL_REPAIR_"
            "CANDIDATE_FOUND_PRE_QM_AUDIT_REQUIRED"
            if best
            else
            "QM_F06_UPPER_V5C_LOCAL_REPAIR_"
            "NO_VALID_CANDIDATE"
        ),
        "samples": SAMPLES,
        "valid_candidate_count": len(records),
        "maximum_heavy_shift_A": MAX_HEAVY_SHIFT_A,
        "minimum_nonnominal_BN_clearance_A": (
            NONCANONICAL_BN_MIN_A
        ),
        "best_candidate": best,
        "orca_authorized": False,
        "RESP_authorized": False,
    }

    REPORT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 100)
    print("QM_F06 UPPER V5-C LOCAL GEOMETRY REPAIR SEARCH")
    print("=" * 100)
    print("Samples:", SAMPLES)
    print("Valid candidates:", len(records))

    if best:
        print()
        print("Best candidate:")

        for key, value in best.items():
            print(f"  {key}: {value}")

        print()
        print("Best XYZ:", BEST_XYZ)

    print()
    print("Decision:", report["decision"])
    print("Ranking:", RANKING_CSV)
    print("Report:", REPORT_JSON)
    print()
    print("ORCA authorized: False")


if __name__ == "__main__":
    main()
