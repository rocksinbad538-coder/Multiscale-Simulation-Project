#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

BASELINE_XYZ = ROOT / (
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
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_baseline_aware_rescoring"
)

OUTPUT_RANKING = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_baseline_aware_ranking.csv"
)

BEST_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_BASELINE_AWARE_START.xyz"
)

REPORT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_BASELINE_AWARE_RESCORING.json"
)

MODIFIED_ATOMS = {
    "S:1739",
    "BR4:UPPER:14:1",
    "H4:UPPER:0203:0",
}

REPLACEMENT_FIELDS = {
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

NONNOMINAL_PAIRS = [
    ("S:1739", "BR4:UPPER:00:3"),
    ("BR4:UPPER:14:1", "BR4:UPPER:00:4"),
]

MINIMUM_NONNOMINAL_CLEARANCE_A = 1.95

# Applied only to bonds incident on a modified atom.
MINIMUM_LOCAL_BOND_MARGIN_A = 0.04

# Unmodified bonds must not deteriorate relative to the baseline
# by more than numerical/rounding tolerance.
MAXIMUM_UNMODIFIED_MARGIN_LOSS_A = 1.0e-6


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
            f"Incomplete XYZ: expected {count}, found {len(atoms)}"
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


def bond_margin(
    first_element,
    second_element,
    value,
):
    limits = bond_limits(
        first_element,
        second_element,
    )

    if limits is None:
        return None

    lower, upper = limits

    return min(
        value - lower,
        upper - value,
    )


def main():
    mapping = read_csv(MAP_CSV)
    valence = read_csv(VALENCE_CSV)
    candidates = read_csv(INPUT_RANKING)
    baseline_atoms = read_xyz(BASELINE_XYZ)

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in mapping
    }

    baseline_xyz = {
        atom_id: baseline_atoms[index]["xyz"]
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

    baseline_margins = {}

    for first, second in nominal_pairs:
        value = distance(
            baseline_xyz[first],
            baseline_xyz[second],
        )

        margin = bond_margin(
            element_by_id[first],
            element_by_id[second],
            value,
        )

        if margin is not None:
            baseline_margins[
                canonical_pair(first, second)
            ] = margin

    rescored = []

    for candidate in candidates:
        coordinates = dict(baseline_xyz)

        for atom_id, fields in REPLACEMENT_FIELDS.items():
            coordinates[atom_id] = tuple(
                float(candidate[field])
                for field in fields
            )

        all_nominal_in_range = True
        minimum_global_margin = float("inf")
        minimum_local_margin = float("inf")
        local_limiting_pair = ""
        local_limiting_distance = None
        unmodified_margin_loss_max = 0.0

        for first, second in nominal_pairs:
            value = distance(
                coordinates[first],
                coordinates[second],
            )

            margin = bond_margin(
                element_by_id[first],
                element_by_id[second],
                value,
            )

            if margin is None:
                continue

            minimum_global_margin = min(
                minimum_global_margin,
                margin,
            )

            if margin < 0.0:
                all_nominal_in_range = False

            pair = canonical_pair(first, second)

            if (
                first in MODIFIED_ATOMS
                or second in MODIFIED_ATOMS
            ):
                if margin < minimum_local_margin:
                    minimum_local_margin = margin
                    local_limiting_pair = (
                        f"{first}--{second}"
                    )
                    local_limiting_distance = value
            else:
                baseline_margin = baseline_margins[pair]

                margin_loss = (
                    baseline_margin - margin
                )

                unmodified_margin_loss_max = max(
                    unmodified_margin_loss_max,
                    margin_loss,
                )

        nonnominal_distances = {
            f"{first}--{second}": distance(
                coordinates[first],
                coordinates[second],
            )
            for first, second in NONNOMINAL_PAIRS
        }

        minimum_nonnominal_clearance = min(
            nonnominal_distances.values()
        )

        baseline_aware_pass = (
            all_nominal_in_range
            and minimum_local_margin
            >= MINIMUM_LOCAL_BOND_MARGIN_A
            and minimum_nonnominal_clearance
            >= MINIMUM_NONNOMINAL_CLEARANCE_A
            and unmodified_margin_loss_max
            <= MAXIMUM_UNMODIFIED_MARGIN_LOSS_A
        )

        maximum_shift = float(
            candidate["maximum_heavy_shift_A"]
        )

        score = (
            1000.0 * minimum_local_margin
            + 100.0 * minimum_nonnominal_clearance
            - 20.0 * maximum_shift
        )

        record = dict(candidate)
        record.update({
            "minimum_global_nominal_margin_A": (
                minimum_global_margin
            ),
            "minimum_local_nominal_margin_A": (
                minimum_local_margin
            ),
            "local_limiting_pair": (
                local_limiting_pair
            ),
            "local_limiting_distance_A": (
                local_limiting_distance
            ),
            "maximum_unmodified_margin_loss_A": (
                unmodified_margin_loss_max
            ),
            "minimum_nonnominal_clearance_A": (
                minimum_nonnominal_clearance
            ),
            "baseline_aware_gate_pass": (
                baseline_aware_pass
            ),
            "baseline_aware_score": score,
        })

        rescored.append(record)

    rescored.sort(
        key=lambda row: (
            row["baseline_aware_gate_pass"],
            float(
                row[
                    "minimum_local_nominal_margin_A"
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
        row["baseline_aware_rank"] = rank

    passing = [
        row
        for row in rescored
        if row["baseline_aware_gate_pass"]
    ]

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "baseline_aware_rank",
        *[
            key
            for key in rescored[0]
            if key != "baseline_aware_rank"
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

    best = passing[0] if passing else None

    if best:
        output_atoms = [
            dict(atom)
            for atom in baseline_atoms
        ]

        for atom_id, fields in REPLACEMENT_FIELDS.items():
            output_atoms[
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
                f"{len(output_atoms)}\n"
            )
            handle.write(
                "QM_F06 UPPER V5-C baseline-aware "
                "local repair; formal pre-QM audit required\n"
            )

            for atom in output_atoms:
                x_value, y_value, z_value = atom["xyz"]

                handle.write(
                    f"{atom['element']:2s} "
                    f"{x_value: .12f} "
                    f"{y_value: .12f} "
                    f"{z_value: .12f}\n"
                )

    report = {
        "decision": (
            "QM_F06_UPPER_V5C_BASELINE_AWARE_"
            "CANDIDATE_FOUND_FORMAL_AUDIT_REQUIRED"
            if best
            else
            "QM_F06_UPPER_V5C_NO_BASELINE_AWARE_"
            "CANDIDATE"
        ),
        "candidate_count": len(candidates),
        "passing_candidate_count": len(passing),
        "minimum_local_bond_margin_A": (
            MINIMUM_LOCAL_BOND_MARGIN_A
        ),
        "minimum_nonnominal_clearance_A": (
            MINIMUM_NONNOMINAL_CLEARANCE_A
        ),
        "best_candidate": best,
        "orca_authorized": False,
        "RESP_authorized": False,
    }

    REPORT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 104)
    print("QM_F06 UPPER V5-C BASELINE-AWARE RESCORING")
    print("=" * 104)
    print("Candidates analyzed:", len(candidates))
    print("Passing candidates:", len(passing))
    print(
        "Required local bond margin A:",
        MINIMUM_LOCAL_BOND_MARGIN_A,
    )

    if best:
        print()
        print("Best candidate:")
        print(
            "  original rank:",
            best["rank"],
        )
        print(
            "  minimum global margin A:",
            best[
                "minimum_global_nominal_margin_A"
            ],
        )
        print(
            "  minimum local margin A:",
            best[
                "minimum_local_nominal_margin_A"
            ],
        )
        print(
            "  local limiting pair:",
            best["local_limiting_pair"],
        )
        print(
            "  local limiting distance A:",
            best[
                "local_limiting_distance_A"
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
