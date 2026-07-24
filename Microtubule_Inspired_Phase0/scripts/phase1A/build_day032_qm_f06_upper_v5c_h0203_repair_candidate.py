#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
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

POSITION_RANKING = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5b_H0203_position_search/"
    "QM_F06_UPPER_V5B_H0203_position_ranking.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_h0203_repair_candidate"
)

OUTPUT_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_H0203_REPAIR_CANDIDATE.xyz"
)

REPORT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_H0203_REPAIR_CANDIDATE.json"
)

TARGET = "H4:UPPER:0203:0"


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
            "xyz": (
                float(fields[1]),
                float(fields[2]),
                float(fields[3]),
            ),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    return atoms


def main():
    with MAP_CSV.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        mapping = list(csv.DictReader(handle))

    with POSITION_RANKING.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        candidates = list(csv.DictReader(handle))

    valid = [
        row
        for row in candidates
        if row[
            "valid_single_owner_position"
        ].strip().lower() == "true"
    ]

    if not valid:
        raise RuntimeError(
            "No valid H0203 position available."
        )

    best = valid[0]
    atoms = read_xyz(SOURCE_XYZ)

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    target_index = index_by_id[TARGET]

    old_xyz = atoms[target_index]["xyz"]

    new_xyz = (
        float(best["x_A"]),
        float(best["y_A"]),
        float(best["z_A"]),
    )

    atoms[target_index]["xyz"] = new_xyz

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with OUTPUT_XYZ.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(
            "QM_F06 UPPER V5-C candidate; "
            "H0203 single-owner geometric repair only; "
            "ORCA not authorized\n"
        )

        for atom in atoms:
            x_value, y_value, z_value = atom["xyz"]

            handle.write(
                f"{atom['element']:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    report = {
        "decision": (
            "QM_F06_UPPER_V5C_H0203_"
            "REPAIR_CANDIDATE_CONSTRUCTED_"
            "HEAVY_CONTACT_AUDIT_REQUIRED"
        ),
        "atom_count": len(atoms),
        "target_atom": TARGET,
        "old_xyz_A": old_xyz,
        "new_xyz_A": new_xyz,
        "selected_position_rank": int(
            best["rank"]
        ),
        "nearest_nonowner_atom": (
            best["nearest_nonowner_atom"]
        ),
        "nearest_nonowner_distance_A": float(
            best[
                "nearest_nonowner_distance_A"
            ]
        ),
        "orca_authorized": False,
        "RESP_authorized": False,
    }

    REPORT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 96)
    print("QM_F06 UPPER V5-C H0203 REPAIR CANDIDATE")
    print("=" * 96)
    print("Target:", TARGET)
    print("Old xyz A:", old_xyz)
    print("New xyz A:", new_xyz)
    print(
        "Nearest nonowner:",
        best["nearest_nonowner_atom"],
    )
    print(
        "Clearance A:",
        best[
            "nearest_nonowner_distance_A"
        ],
    )
    print()
    print("Candidate XYZ:", OUTPUT_XYZ)
    print("Report:", REPORT_JSON)
    print()
    print("ORCA authorized: False")


if __name__ == "__main__":
    main()
