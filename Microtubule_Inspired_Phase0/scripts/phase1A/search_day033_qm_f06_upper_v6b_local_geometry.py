#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6a_topology_closure"
)

SOURCE_XYZ = INPUT_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_start.xyz"
)

SOURCE_MAP = INPUT_DIR / (
    "QM_F06_UPPER_V6A_atom_role_provenance_map.csv"
)

SOURCE_EDGES = INPUT_DIR / (
    "QM_F06_UPPER_V6A_nominal_edges.csv"
)

SOURCE_REPORT = INPUT_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_REPORT.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_local_geometry_search"
)

BEST_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_BEST.xyz"
)

RANKING_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_local_geometry_ranking.csv"
)

REPORT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_SEARCH.json"
)


SAMPLES = 600_000
RANDOM_SEED = 20260724

S_CENTER = "S:1739"
S_CLOSURE_PARTNER = "BR4:UPPER:00:3"

N_CENTER = "BR4:UPPER:14:1"
N_CLOSURE_PARTNER = "BR4:UPPER:00:4"

N_ATTACHED_CAP = (
    "HCAPV5B:UPPER:BR4_14_1:BR4_14_2"
)

CHANGED_ATOMS = {
    S_CENTER,
    N_CENTER,
    N_ATTACHED_CAP,
}

CLOSURE_PAIRS = {
    tuple(sorted((
        S_CENTER,
        S_CLOSURE_PARTNER,
    ))),
    tuple(sorted((
        N_CENTER,
        N_CLOSURE_PARTNER,
    ))),
}

BN_MIN_A = 1.25
BN_MAX_A = 1.90

BH_MIN_A = 0.90
BH_MAX_A = 1.35

NH_MIN_A = 0.80
NH_MAX_A = 1.25

LOCAL_BOND_MARGIN_A = 0.08
NONNOMINAL_CLEARANCE_MARGIN_A = 0.06

HH_HARD_CONTACT_A = 0.70
HX_HARD_CONTACT_A = 0.75
HEAVY_HEAVY_HARD_CONTACT_A = 1.10

CLOSURE_RADIUS_MIN_A = 1.40
CLOSURE_RADIUS_MAX_A = 1.68

DIRECTION_PERTURBATION_SIGMA = 0.22
MAXIMUM_CENTER_SHIFT_A = 0.95


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path) -> list[dict]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0].strip())
    coordinate_lines = lines[2:2 + count]

    if len(coordinate_lines) != count:
        raise RuntimeError(
            "Incomplete source XYZ."
        )

    atoms = []

    for index, line in enumerate(
        coordinate_lines
    ):
        fields = line.split()

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": tuple(
                map(float, fields[1:4])
            ),
        })

    return atoms


def canonical_pair(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def subtract(a, b):
    return tuple(
        x - y
        for x, y in zip(a, b)
    )


def add(a, b):
    return tuple(
        x + y
        for x, y in zip(a, b)
    )


def scale(vector, factor):
    return tuple(
        factor * value
        for value in vector
    )


def norm(vector) -> float:
    return math.sqrt(
        sum(value * value for value in vector)
    )


def normalize(vector):
    value = norm(vector)

    if value == 0.0:
        raise RuntimeError(
            "Cannot normalize zero vector."
        )

    return tuple(
        component / value
        for component in vector
    )


def distance(a, b) -> float:
    return norm(subtract(a, b))


def bond_limits(
    first_element: str,
    second_element: str,
):
    pair = {
        first_element,
        second_element,
    }

    if pair == {"B", "N"}:
        return BN_MIN_A, BN_MAX_A

    if pair == {"B", "H"}:
        return BH_MIN_A, BH_MAX_A

    if pair == {"N", "H"}:
        return NH_MIN_A, NH_MAX_A

    return None


def hard_contact_threshold(
    first_element: str,
    second_element: str,
) -> float:
    if (
        first_element == "H"
        and second_element == "H"
    ):
        return HH_HARD_CONTACT_A

    if (
        first_element == "H"
        or second_element == "H"
    ):
        return HX_HARD_CONTACT_A

    return HEAVY_HEAVY_HARD_CONTACT_A


def perturbed_radial_position(
    rng: random.Random,
    old_center,
    fixed_partner,
):
    base_direction = normalize(
        subtract(
            old_center,
            fixed_partner,
        )
    )

    trial_direction = normalize((
        base_direction[0]
        + rng.gauss(
            0.0,
            DIRECTION_PERTURBATION_SIGMA,
        ),
        base_direction[1]
        + rng.gauss(
            0.0,
            DIRECTION_PERTURBATION_SIGMA,
        ),
        base_direction[2]
        + rng.gauss(
            0.0,
            DIRECTION_PERTURBATION_SIGMA,
        ),
    ))

    radius = rng.uniform(
        CLOSURE_RADIUS_MIN_A,
        CLOSURE_RADIUS_MAX_A,
    )

    return add(
        fixed_partner,
        scale(trial_direction, radius),
    )


def main() -> None:
    for path in (
        SOURCE_XYZ,
        SOURCE_MAP,
        SOURCE_EDGES,
        SOURCE_REPORT,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    atoms = read_xyz(SOURCE_XYZ)
    map_rows = read_csv(SOURCE_MAP)
    edge_rows = read_csv(SOURCE_EDGES)

    retained_rows = [
        row
        for row in map_rows
        if row["v6a_retained"].strip().lower()
        == "true"
    ]

    retained_rows.sort(
        key=lambda row: int(
            row["v6a_index_0based"]
        )
    )

    if len(atoms) != len(retained_rows):
        raise RuntimeError(
            "XYZ/map retained-atom mismatch."
        )

    id_by_index = []
    index_by_id = {}
    element_by_id = {}
    source_coordinates = {}

    for atom, row in zip(
        atoms,
        retained_rows,
    ):
        atom_id = row["atom_id"]

        if atom["element"] != row["element"]:
            raise RuntimeError(
                f"Element mismatch for {atom_id}."
            )

        id_by_index.append(atom_id)
        index_by_id[atom_id] = atom["index"]
        element_by_id[atom_id] = atom["element"]
        source_coordinates[atom_id] = (
            atom["xyz_A"]
        )

    missing_changed = (
        CHANGED_ATOMS - set(index_by_id)
    )

    if missing_changed:
        raise RuntimeError(
            "Missing changed atoms: "
            + "|".join(sorted(missing_changed))
        )

    nominal_edges = set()
    adjacency = {
        atom_id: set()
        for atom_id in index_by_id
    }

    for row in edge_rows:
        pair = canonical_pair(
            row["first_atom"],
            row["second_atom"],
        )

        nominal_edges.add(pair)
        adjacency[pair[0]].add(pair[1])
        adjacency[pair[1]].add(pair[0])

    for pair in CLOSURE_PAIRS:
        if pair not in nominal_edges:
            raise RuntimeError(
                f"Missing V6-A closure edge: {pair}"
            )

    old_s = source_coordinates[S_CENTER]
    old_n = source_coordinates[N_CENTER]
    old_cap = source_coordinates[N_ATTACHED_CAP]

    fixed_s_partner = source_coordinates[
        S_CLOSURE_PARTNER
    ]

    fixed_n_partner = source_coordinates[
        N_CLOSURE_PARTNER
    ]

    rng = random.Random(RANDOM_SEED)

    valid_records = []

    local_nominal_edges = sorted({
        pair
        for pair in nominal_edges
        if (
            pair[0] in CHANGED_ATOMS
            or pair[1] in CHANGED_ATOMS
        )
    })

    candidate_pairs = []

    atom_ids = list(index_by_id)

    for first_position, first in enumerate(
        atom_ids
    ):
        for second in atom_ids[
            first_position + 1:
        ]:
            pair = canonical_pair(
                first,
                second,
            )

            if not (
                pair[0] in CHANGED_ATOMS
                or pair[1] in CHANGED_ATOMS
            ):
                continue

            candidate_pairs.append(pair)

    for sample_index in range(SAMPLES):
        new_s = perturbed_radial_position(
            rng,
            old_s,
            fixed_s_partner,
        )

        new_n = perturbed_radial_position(
            rng,
            old_n,
            fixed_n_partner,
        )

        s_shift = distance(
            old_s,
            new_s,
        )

        n_shift = distance(
            old_n,
            new_n,
        )

        maximum_center_shift = max(
            s_shift,
            n_shift,
        )

        if (
            maximum_center_shift
            > MAXIMUM_CENTER_SHIFT_A
        ):
            continue

        n_translation = subtract(
            new_n,
            old_n,
        )

        new_cap = add(
            old_cap,
            n_translation,
        )

        coordinates = {
            S_CENTER: new_s,
            N_CENTER: new_n,
            N_ATTACHED_CAP: new_cap,
        }

        def xyz(atom_id):
            return coordinates.get(
                atom_id,
                source_coordinates[atom_id],
            )

        local_bond_records = []
        minimum_local_margin = float("inf")
        valid = True

        for first, second in (
            local_nominal_edges
        ):
            first_element = element_by_id[first]
            second_element = element_by_id[second]

            limits = bond_limits(
                first_element,
                second_element,
            )

            if limits is None:
                valid = False
                break

            minimum, maximum = limits

            value = distance(
                xyz(first),
                xyz(second),
            )

            margin = min(
                value - minimum,
                maximum - value,
            )

            local_bond_records.append((
                first,
                second,
                value,
                margin,
            ))

            minimum_local_margin = min(
                minimum_local_margin,
                margin,
            )

            if margin < LOCAL_BOND_MARGIN_A:
                valid = False
                break

        if not valid:
            continue

        minimum_nonnominal_clearance = (
            float("inf")
        )

        limiting_nonnominal_pair = None

        for first, second in candidate_pairs:
            pair = canonical_pair(
                first,
                second,
            )

            value = distance(
                xyz(first),
                xyz(second),
            )

            first_element = element_by_id[first]
            second_element = element_by_id[second]

            hard_threshold = (
                hard_contact_threshold(
                    first_element,
                    second_element,
                )
            )

            if value < hard_threshold:
                valid = False
                break

            if pair in nominal_edges:
                continue

            limits = bond_limits(
                first_element,
                second_element,
            )

            if limits is None:
                continue

            _, maximum = limits

            clearance = value - maximum

            if (
                clearance
                < minimum_nonnominal_clearance
            ):
                minimum_nonnominal_clearance = (
                    clearance
                )
                limiting_nonnominal_pair = (
                    f"{first}--{second}"
                )

            if (
                clearance
                < NONNOMINAL_CLEARANCE_MARGIN_A
            ):
                valid = False
                break

        if not valid:
            continue

        closure_distances = {
            "--".join(pair): distance(
                xyz(pair[0]),
                xyz(pair[1]),
            )
            for pair in sorted(
                CLOSURE_PAIRS
            )
        }

        maximum_closure_deviation = max(
            abs(value - 1.50)
            for value in closure_distances.values()
        )

        score = (
            500.0 * minimum_local_margin
            + 40.0 * minimum_nonnominal_clearance
            - 15.0 * maximum_center_shift
            - 20.0 * maximum_closure_deviation
        )

        record = {
            "sample_index": sample_index,
            "score": score,
            "minimum_local_bond_margin_A": (
                minimum_local_margin
            ),
            "minimum_nonnominal_clearance_A": (
                minimum_nonnominal_clearance
            ),
            "limiting_nonnominal_pair": (
                limiting_nonnominal_pair or ""
            ),
            "maximum_center_shift_A": (
                maximum_center_shift
            ),
            "S1739_shift_A": s_shift,
            "BR4_14_1_shift_A": n_shift,
            "S1739_BR4_00_3_distance_A": (
                closure_distances[
                    "--".join(
                        canonical_pair(
                            S_CENTER,
                            S_CLOSURE_PARTNER,
                        )
                    )
                ]
            ),
            "BR4_14_1_BR4_00_4_distance_A": (
                closure_distances[
                    "--".join(
                        canonical_pair(
                            N_CENTER,
                            N_CLOSURE_PARTNER,
                        )
                    )
                ]
            ),
            "S1739_x_A": new_s[0],
            "S1739_y_A": new_s[1],
            "S1739_z_A": new_s[2],
            "BR4_14_1_x_A": new_n[0],
            "BR4_14_1_y_A": new_n[1],
            "BR4_14_1_z_A": new_n[2],
            "cap_x_A": new_cap[0],
            "cap_y_A": new_cap[1],
            "cap_z_A": new_cap[2],
        }

        valid_records.append(record)

    valid_records.sort(
        key=lambda row: (
            -row["score"],
            -row[
                "minimum_local_bond_margin_A"
            ],
            -row[
                "minimum_nonnominal_clearance_A"
            ],
            row["maximum_center_shift_A"],
            row["sample_index"],
        )
    )

    for rank, row in enumerate(
        valid_records,
        start=1,
    ):
        row["rank"] = rank

    ranking_fieldnames = [
        "rank",
        "sample_index",
        "score",
        "minimum_local_bond_margin_A",
        "minimum_nonnominal_clearance_A",
        "limiting_nonnominal_pair",
        "maximum_center_shift_A",
        "S1739_shift_A",
        "BR4_14_1_shift_A",
        "S1739_BR4_00_3_distance_A",
        "BR4_14_1_BR4_00_4_distance_A",
        "S1739_x_A",
        "S1739_y_A",
        "S1739_z_A",
        "BR4_14_1_x_A",
        "BR4_14_1_y_A",
        "BR4_14_1_z_A",
        "cap_x_A",
        "cap_y_A",
        "cap_z_A",
    ]

    with RANKING_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ranking_fieldnames,
        )

        writer.writeheader()
        writer.writerows(
            valid_records[:500]
        )

    best = (
        valid_records[0]
        if valid_records
        else None
    )

    if best:
        best_coordinates = dict(
            source_coordinates
        )

        best_coordinates[S_CENTER] = (
            best["S1739_x_A"],
            best["S1739_y_A"],
            best["S1739_z_A"],
        )

        best_coordinates[N_CENTER] = (
            best["BR4_14_1_x_A"],
            best["BR4_14_1_y_A"],
            best["BR4_14_1_z_A"],
        )

        best_coordinates[N_ATTACHED_CAP] = (
            best["cap_x_A"],
            best["cap_y_A"],
            best["cap_z_A"],
        )

        with BEST_XYZ.open(
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(
                f"{len(atom_ids)}\n"
            )
            handle.write(
                "QM_F06 UPPER V6-B local topology-"
                "closure geometry candidate; "
                "ORCA not yet authorized\n"
            )

            for atom_id in atom_ids:
                x_value, y_value, z_value = (
                    best_coordinates[atom_id]
                )

                handle.write(
                    f"{element_by_id[atom_id]:2s} "
                    f"{x_value: .12f} "
                    f"{y_value: .12f} "
                    f"{z_value: .12f}\n"
                )

    report = {
        "decision": (
            "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_"
            "CANDIDATE_FOUND_FORMAL_AUDIT_REQUIRED"
            if best
            else
            "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_"
            "NO_VALID_CANDIDATE"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "sample_count": SAMPLES,
        "random_seed": RANDOM_SEED,
        "valid_candidate_count": len(
            valid_records
        ),
        "best_candidate": best,
        "criteria": {
            "local_bond_margin_A": (
                LOCAL_BOND_MARGIN_A
            ),
            "nonnominal_clearance_margin_A": (
                NONNOMINAL_CLEARANCE_MARGIN_A
            ),
            "closure_radius_min_A": (
                CLOSURE_RADIUS_MIN_A
            ),
            "closure_radius_max_A": (
                CLOSURE_RADIUS_MAX_A
            ),
            "maximum_center_shift_A": (
                MAXIMUM_CENTER_SHIFT_A
            ),
        },
        "files": {
            "best_xyz": (
                str(BEST_XYZ.relative_to(ROOT))
                if best
                else None
            ),
            "ranking": str(
                RANKING_CSV.relative_to(ROOT)
            ),
        },
        "sha256": {
            "best_xyz": (
                sha256(BEST_XYZ)
                if best
                else None
            ),
            "ranking": sha256(
                RANKING_CSV
            ),
        },
        "formal_pre_qm_audit_authorized": (
            best is not None
        ),
        "orca_input_design_authorized": False,
        "orca_execution_authorized": False,
        "RESP_authorized": False,
        "MD_authorized": False,
    }

    REPORT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 104)
    print("QM_F06 UPPER V6-B LOCAL GEOMETRY SEARCH")
    print("=" * 104)
    print("Samples:", SAMPLES)
    print(
        "Valid candidates:",
        len(valid_records),
    )

    if best:
        print()
        print("Best candidate:")

        for key in (
            "rank",
            "sample_index",
            "score",
            "minimum_local_bond_margin_A",
            "minimum_nonnominal_clearance_A",
            "limiting_nonnominal_pair",
            "maximum_center_shift_A",
            "S1739_shift_A",
            "BR4_14_1_shift_A",
            "S1739_BR4_00_3_distance_A",
            "BR4_14_1_BR4_00_4_distance_A",
        ):
            print(f"  {key}: {best[key]}")

        print()
        print("Best XYZ:", BEST_XYZ)

    print()
    print("Decision:", report["decision"])
    print("Ranking:", RANKING_CSV)
    print("Report:", REPORT_JSON)
    print()
    print(
        "Formal pre-QM audit authorized:",
        best is not None,
    )
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
