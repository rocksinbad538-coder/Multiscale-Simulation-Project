#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SOURCE_XYZ = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_baseline_aware_rescoring/"
    "QM_F06_UPPER_V5C_BASELINE_AWARE_START.xyz"
)

SOURCE_MAP = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

SOURCE_CAPS = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_new_artificial_caps.csv"
)

RESCORING_REPORT = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_baseline_aware_rescoring/"
    "QM_F06_UPPER_V5C_BASELINE_AWARE_RESCORING.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_construction"
)

OUTPUT_XYZ = OUTPUT_DIR / "QM_F06_UPPER_V5C_start.xyz"

OUTPUT_MAP = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_atom_role_provenance_map.csv"
)

OUTPUT_CAPS = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_artificial_caps.csv"
)

OUTPUT_REPORT = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_CONSTRUCTION_REPORT.json"
)

REPAIRED_ATOMS = {
    "S:1739",
    "BR4:UPPER:14:1",
    "H4:UPPER:0203:0",
}

EXPECTED_COUNT = 52
EXPECTED_COMPOSITION = {
    "B": 16,
    "N": 13,
    "H": 23,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


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
            f"Incomplete XYZ: expected {count}, "
            f"found {len(atoms)}"
        )

    return atoms


def main() -> None:
    for path in (
        SOURCE_XYZ,
        SOURCE_MAP,
        SOURCE_CAPS,
        RESCORING_REPORT,
    ):
        require_file(path)

    atoms = read_xyz(SOURCE_XYZ)

    composition = {}

    for atom in atoms:
        element = atom["element"]
        composition[element] = (
            composition.get(element, 0) + 1
        )

    if len(atoms) != EXPECTED_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_COUNT} atoms, "
            f"found {len(atoms)}"
        )

    if composition != EXPECTED_COMPOSITION:
        raise RuntimeError(
            f"Composition mismatch: {composition}"
        )

    with SOURCE_MAP.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        mapping = list(csv.DictReader(handle))

    if len(mapping) != len(atoms):
        raise RuntimeError(
            "Provenance-map/XYZ count mismatch."
        )

    fieldnames = list(mapping[0])

    for required in (
        "v5c_coordinate_status",
        "v5c_coordinate_basis",
    ):
        if required not in fieldnames:
            fieldnames.append(required)

    repaired_found = set()

    for row in mapping:
        atom_id = row["atom_id"]

        if atom_id in REPAIRED_ATOMS:
            repaired_found.add(atom_id)
            row["v5c_coordinate_status"] = (
                "LOCAL_GEOMETRY_REPAIRED"
            )
            row["v5c_coordinate_basis"] = (
                "BASELINE_AWARE_CANDIDATE_RANK_1"
            )
        else:
            row["v5c_coordinate_status"] = (
                "RETAINED_FROM_V5B"
            )
            row["v5c_coordinate_basis"] = (
                "UNCHANGED_COORDINATE"
            )

    if repaired_found != REPAIRED_ATOMS:
        raise RuntimeError(
            "Missing repaired atoms in provenance map: "
            f"{sorted(REPAIRED_ATOMS - repaired_found)}"
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copy2(SOURCE_XYZ, OUTPUT_XYZ)
    shutil.copy2(SOURCE_CAPS, OUTPUT_CAPS)

    with OUTPUT_MAP.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(mapping)

    rescoring = json.loads(
        RESCORING_REPORT.read_text(
            encoding="utf-8"
        )
    )

    best = rescoring["best_candidate"]

    report = {
        "decision": (
            "QM_F06_UPPER_V5C_CONSTRUCTED_"
            "FORMAL_PRE_QM_AUDIT_REQUIRED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "atom_count": len(atoms),
        "composition": composition,
        "repaired_atoms": sorted(REPAIRED_ATOMS),
        "selection": {
            "baseline_aware_rank": 1,
            "original_rank": int(best["rank"]),
            "minimum_global_margin_A": float(
                best[
                    "minimum_global_nominal_margin_A"
                ]
            ),
            "minimum_local_margin_A": float(
                best[
                    "minimum_local_nominal_margin_A"
                ]
            ),
            "minimum_nonnominal_clearance_A": float(
                best[
                    "minimum_nonnominal_clearance_A"
                ]
            ),
            "maximum_heavy_shift_A": float(
                best["maximum_heavy_shift_A"]
            ),
        },
        "files": {
            "start_xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "provenance_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
            "artificial_caps": str(
                OUTPUT_CAPS.relative_to(ROOT)
            ),
            "rescoring_report": str(
                RESCORING_REPORT.relative_to(ROOT)
            ),
        },
        "sha256": {
            "start_xyz": sha256(OUTPUT_XYZ),
            "provenance_map": sha256(OUTPUT_MAP),
            "artificial_caps": sha256(
                OUTPUT_CAPS
            ),
            "rescoring_report": sha256(
                RESCORING_REPORT
            ),
        },
        "pre_qm_audit_authorized": True,
        "orca_authorized": False,
        "RESP_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 88)
    print("QM_F06 UPPER V5-C FORMAL CONSTRUCTION")
    print("=" * 88)
    print("Atom count:", len(atoms))
    print("Composition:", composition)
    print("Repaired atoms:")

    for atom_id in sorted(REPAIRED_ATOMS):
        print("  ", atom_id)

    print()
    print("Minimum local margin A:",
          report["selection"]["minimum_local_margin_A"])
    print("Nonnominal clearance A:",
          report["selection"]["minimum_nonnominal_clearance_A"])
    print("Maximum heavy shift A:",
          report["selection"]["maximum_heavy_shift_A"])

    print()
    print("Decision:", report["decision"])
    print("XYZ:", OUTPUT_XYZ)
    print("Map:", OUTPUT_MAP)
    print("Report:", OUTPUT_REPORT)
    print()
    print("Pre-QM audit authorized: True")
    print("ORCA authorized: False")


if __name__ == "__main__":
    main()
