#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
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

INPUT_RANKING = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_local_geometry_repair/"
    "QM_F06_UPPER_V5C_local_geometry_repair_ranking.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_robust_rescoring"
)

OUTPUT_RANKING = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_robust_candidate_ranking.csv"
)

BEST_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_ROBUST_START.xyz"
)

REPORT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_ROBUST_RESCORING.json"
)

REPLACED_ATOMS = {
    "S:1739": (
        "S1739_x_A",
        "S1739_y_A",
        "S1739_z_A",
    ),
    "BR4:UPPER:14:1": (
        "BR4_14_1_x_A",
        "BR4_14_1_y_A",
        "BR4_14_1_z_A",
    ),
    "H4:UPPER:0203:0": (
        "H0203_x_A",
        "H0203_y_A",
        "H0203_z_A",
    ),
}

NONCANONICAL_PAIRS = [
    ("S:1739", "BR4:UPPER:00:3"),
    ("BR4:UPPER:14:1", "BR4:UPPER:00:4"),
]

REQUIRED_NONCANONICAL_CLEARANCE_A = 1.95

# Require a real buffer from the raw acceptance boundaries.
MINIMUM_REQUIRED_BOND_MARGIN_A = 0.08


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


def canonical_pair(first, second):
    return tuple(sorted((first, second)))


def bond_limits(first_element, second_element):
    pair = frozenset((
        first_element,
        second_element,
    ))

    if pair == frozenset(("B", "N")):
        return 1.25, 1.85

    if pair == frozenset(("B", "H")):
        return 0.95, 1.35

    if pair == frozenset(("N", "H")):
        return 0.85, 1.25

    return None


def main():
    mapping = read_csv(MAP_CSV)
    valence = read_csv(VALENCE_CSV)
    candidates = read_csv(INPUT_RANKING)
    source_atoms = read_xyz(SOURCE_XYZ)

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in mapping
    }

    base_xyz = {
        atom_id: source_atoms[index]["xyz"]
        for atom_id, index in index_by_id.items()
    }

    nominal_pairs = set()

    for row in valence:
        for neighbor in row["neighbors"].split("|"):
            if neighbor:
                nominal_pairs.add(
                    canonical_pair(
                        row["atom_id"],
                        neighbor,
                    )
                )

    rescored = []

    for candidate in candidates:
        coordinates = dict(base_xyz)

        for atom_id, fields in REPLACED_ATOMS.items():
            coordinates[atom_id] = tuple(
                float(candidate[field])
                for field in fields
            )

        minimum_margin = float("inf")
        limiting_pair = None
        limiting_distance = None
        nominal_pass = True

        for first, second in nominal_pairs:
            limits = bond_limits(
                element_by_id[first],
                element_by_id[second],
            )

            if limits is None:
                continue

            value = distance(
                coordinates[first],
                coordinates[second],
            )

            lower, upper = limits
            margin = min(
                value - lower,
                upper - value,
            )

            if margin < 0.0:
                nominal_pass = False

            if margin < minimum_margin:
                minimum_margin = margin
                limiting_pair = (
                    f"{first}--{second}"
                )
                limiting_distance = value

        noncanonical_distances = {
            f"{first}--{second}": distance(
                coordinates[first],
                coordinates[second],
            )
            for first, second in NONCANONICAL_PAIRS
        }

        minimum_noncanonical_clearance = min(
            noncanonical_distances.values()
        )

        robust_pass = (
            nominal_pass
            and minimum_margin
            >= MINIMUM_REQUIRED_BOND_MARGIN_A
            and minimum_noncanonical_clearance
            >= REQUIRED_NONCANONICAL_CLEARANCE_A
        )

        maximum_shift = float(
            candidate["maximum_heavy_shift_A"]
        )

        # Lexicographic intent encoded numerically:
        # maximize weakest nominal-bond margin,
        # then maximize noncanonical clearance,
        # then minimize heavy-atom displacement.
        robust_score = (
            1000.0 * minimum_margin
            + 100.0 * minimum_noncanonical_clearance
            - 10.0 * maximum_shift
        )

        record = dict(candidate)
        record.update({
            "minimum_nominal_bond_margin_A": (
                minimum_margin
            ),
            "limiting_nominal_pair": limiting_pair,
            "limiting_nominal_distance_A": (
                limiting_distance
            ),
            "minimum_nonnominal_clearance_A": (
                minimum_noncanonical_clearance
            ),
            "robust_gate_pass": robust_pass,
            "robust_score": robust_score,
        })

        rescored.append(record)

    rescored.sort(
        key=lambda row: (
            row["robust_gate_pass"],
            float(
                row[
                    "minimum_nominal_bond_margin_A"
                ]
            ),
            float(
                row[
                    "minimum_nonnominal_clearance_A"
                ]
            ),
            -float(row["maximum_heavy_shift_A"]),
        ),
        reverse=True,
    )

    for rank, row in enumerate(
        rescored,
        start=1,
    ):
        row["robust_rank"] = rank

    robust = [
        row
        for row in rescored
        if row["robust_gate_pass"]
    ]

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "robust_rank",
        *[
            key
            for key in rescored[0]
            if key != "robust_rank"
        ],
    ]

    with OUTPUT_RANKING.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rescored)

    best = robust[0] if robust else None

    if best:
        repaired_atoms = [
            dict(atom)
            for atom in source_atoms
        ]

        for atom_id, fields in REPLACED_ATOMS.items():
            repaired_atoms[
                index_by_id[atom_id]
            ]["xyz"] = tuple(
                float(best[field])
                for field in fields
            )

        with BEST_XYZ.open(
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(
                f"{len(repaired_atoms)}\n"
            )
            handle.write(
                "QM_F06 UPPER V5-C robust local "
                "repair; pending formal pre-QM audit\n"
            )

            for atom in repaired_atoms:
                x_value, y_value, z_value = atom["xyz"]

                handle.write(
                    f"{atom['element']:2s} "
                    f"{x_value: .12f} "
                    f"{y_value: .12f} "
                    f"{z_value: .12f}\n"
                )

    report = {
        "decision": (
            "QM_F06_UPPER_V5C_ROBUST_CANDIDATE_"
            "FOUND_FORMAL_PRE_QM_AUDIT_REQUIRED"
            if best
            else
            "QM_F06_UPPER_V5C_NO_ROBUST_"
            "LOCAL_REPAIR_CANDIDATE"
        ),
        "input_candidate_count": len(candidates),
        "robust_candidate_count": len(robust),
        "minimum_required_bond_margin_A": (
            MINIMUM_REQUIRED_BOND_MARGIN_A
        ),
        "required_nonnominal_clearance_A": (
            REQUIRED_NONCANONICAL_CLEARANCE_A
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
    print("QM_F06 UPPER V5-C ROBUST CANDIDATE RESCORING")
    print("=" * 100)
    print("Input candidates:", len(candidates))
    print("Robust candidates:", len(robust))
    print(
        "Required minimum bond margin A:",
        MINIMUM_REQUIRED_BOND_MARGIN_A,
    )

    if best:
        print()
        print("Best robust candidate:")
        print(
            "  original rank:",
            best["rank"],
        )
        print(
            "  minimum nominal bond margin A:",
            best[
                "minimum_nominal_bond_margin_A"
            ],
        )
        print(
            "  limiting nominal pair:",
            best["limiting_nominal_pair"],
        )
        print(
            "  limiting nominal distance A:",
            best[
                "limiting_nominal_distance_A"
            ],
        )
        print(
            "  minimum nonnominal clearance A:",
            best[
                "minimum_nonnominal_clearance_A"
            ],
        )
        print(
            "  maximum heavy shift A:",
            best["maximum_heavy_shift_A"],
        )
        print("  XYZ:", BEST_XYZ)

    print()
    print("Decision:", report["decision"])
    print("Ranking:", OUTPUT_RANKING)
    print("Report:", REPORT_JSON)
    print("ORCA authorized: False")


if __name__ == "__main__":
    main()
