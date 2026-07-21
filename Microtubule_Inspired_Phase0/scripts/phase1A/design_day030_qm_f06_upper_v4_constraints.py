#!/usr/bin/env python3
"""
Design the conservative partial-relaxation constraints for QM_F06 UPPER V4.

Policy:
- retain as fixed only real V3-A2 atoms that were already fixed;
- never fix an artificial boundary cap;
- keep every restored V4 atom mobile;
- preserve the previously mobile V3-A2 boundary atoms as mobile.

This script creates and audits the constraint map only.
It does not generate or authorize an ORCA calculation.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V4_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction"
)

V4_XYZ = V4_DIR / "QM_F06_UPPER_V4_start.xyz"

V4_MAP = (
    V4_DIR
    / "QM_F06_UPPER_V4_atom_role_provenance_map.csv"
)

PRE_QM_REPORT = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_pre_qm_audit/"
    "QM_F06_UPPER_V4_PRE_QM_AUDIT.json"
)

V3A2_MAP = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow/"
    "v3a2_atom_role_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day030_qm_f06_upper_v4_constraint_design"
)

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_constraint_map.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_CONSTRAINT_DESIGN.json"
)

EXPECTED_ATOM_COUNT = 46


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for path in (
        V4_XYZ,
        V4_MAP,
        PRE_QM_REPORT,
        V3A2_MAP,
    ):
        require_file(path)

    pre_qm = json.loads(
        PRE_QM_REPORT.read_text(encoding="utf-8")
    )

    expected_decision = (
        "QM_F06_UPPER_V4_PRE_QM_STRUCTURAL_GATE_PASS_"
        "CONSTRAINT_DESIGN_AUTHORIZED"
    )

    if pre_qm["decision"] != expected_decision:
        raise RuntimeError(
            "Unexpected pre-QM decision: "
            f"{pre_qm['decision']}"
        )

    if not pre_qm["authorization"][
        "constraint_design_authorized"
    ]:
        raise RuntimeError(
            "Constraint design is not authorized."
        )

    v4_rows = read_csv(V4_MAP)
    v3_rows = read_csv(V3A2_MAP)

    if len(v4_rows) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_ATOM_COUNT} V4 atoms; "
            f"found {len(v4_rows)}"
        )

    v3_by_id = {
        row["atom_id"]: row
        for row in v3_rows
    }

    output_rows = []

    for expected_index, row in enumerate(v4_rows):
        index = int(row["index_0based"])

        if index != expected_index:
            raise RuntimeError(
                f"Non-deterministic V4 index at {index}"
            )

        atom_id = row["atom_id"]
        artificial_cap = parse_bool(
            row["artificial_cap"]
        )

        retained_from_v3 = (
            row["coordinate_source"]
            == "RETAINED_VALIDATED_V3A2_START"
        )

        prior = v3_by_id.get(atom_id)

        if artificial_cap:
            fixed = False
            basis = "ALL_ARTIFICIAL_CAPS_MOBILE"

        elif not retained_from_v3:
            fixed = False
            basis = "RESTORED_V4_REAL_ATOM_MOBILE"

        elif prior is None:
            raise RuntimeError(
                "Retained V3 atom missing from V3-A2 map: "
                f"{atom_id}"
            )

        elif parse_bool(prior["v3a2_fixed"]):
            fixed = True
            basis = "RETAINED_VALIDATED_REAL_CORE"

        elif parse_bool(prior["v3a2_mobile"]):
            fixed = False
            basis = "RETAINED_PREVIOUSLY_MOBILE_BOUNDARY"

        else:
            raise RuntimeError(
                "V3-A2 atom has neither fixed nor mobile "
                f"classification: {atom_id}"
            )

        output_rows.append({
            "index_0based": index,
            "atom_id": atom_id,
            "element": row["element"],
            "atom_role": row["atom_role"],
            "node_type": row["node_type"],
            "coordinate_source": row[
                "coordinate_source"
            ],
            "artificial_cap": artificial_cap,
            "retained_from_v3a2": retained_from_v3,
            "v4_fixed": fixed,
            "v4_mobile": not fixed,
            "v4_constraint_basis": basis,
        })

    fixed_rows = [
        row for row in output_rows
        if row["v4_fixed"]
    ]

    mobile_rows = [
        row for row in output_rows
        if row["v4_mobile"]
    ]

    artificial_cap_rows = [
        row for row in output_rows
        if row["artificial_cap"]
    ]

    restored_real_rows = [
        row for row in output_rows
        if row["v4_constraint_basis"]
        == "RESTORED_V4_REAL_ATOM_MOBILE"
    ]

    gates = {
        "atom_count": (
            len(output_rows) == EXPECTED_ATOM_COUNT
        ),
        "exclusive_fixed_mobile": all(
            row["v4_fixed"] != row["v4_mobile"]
            for row in output_rows
        ),
        "all_artificial_caps_mobile": all(
            row["v4_mobile"]
            for row in artificial_cap_rows
        ),
        "all_restored_real_atoms_mobile": all(
            row["v4_mobile"]
            for row in restored_real_rows
        ),
        "no_artificial_cap_fixed": not any(
            row["v4_fixed"]
            for row in artificial_cap_rows
        ),
        "nonempty_fixed_core": len(fixed_rows) >= 3,
        "nonempty_mobile_region": len(mobile_rows) >= 1,
    }

    overall_pass = all(gates.values())

    fieldnames = list(output_rows[0])

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
        writer.writerows(output_rows)

    decision = (
        "QM_F06_UPPER_V4_CONSTRAINT_DESIGN_PASS_"
        "ORCA_INPUT_PREPARATION_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V4_CONSTRAINT_DESIGN_FAIL_"
        "ORCA_INPUT_BLOCKED"
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "policy": (
            "FIX_ONLY_RETAINED_REAL_V3A2_FIXED_CORE_"
            "ALL_CAPS_AND_RESTORED_ATOMS_MOBILE"
        ),
        "atom_count": len(output_rows),
        "fixed_atom_count": len(fixed_rows),
        "mobile_atom_count": len(mobile_rows),
        "fixed_indices": [
            row["index_0based"]
            for row in fixed_rows
        ],
        "mobile_indices": [
            row["index_0based"]
            for row in mobile_rows
        ],
        "artificial_cap_count": len(
            artificial_cap_rows
        ),
        "restored_real_atom_count": len(
            restored_real_rows
        ),
        "constraint_basis_counts": dict(Counter(
            row["v4_constraint_basis"]
            for row in output_rows
        )),
        "gates": gates,
        "overall_pass": overall_pass,
        "files": {
            "constraint_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "v4_xyz": sha256(V4_XYZ),
            "v4_provenance_map": sha256(V4_MAP),
            "pre_qm_report": sha256(
                PRE_QM_REPORT
            ),
            "v3a2_map": sha256(V3A2_MAP),
            "constraint_map": sha256(
                OUTPUT_MAP
            ),
        },
        "authorization": {
            "constraint_design_accepted": overall_pass,
            "orca_input_preparation_authorized": (
                overall_pass
            ),
            "orca_execution_authorized": False,
            "geometric_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 CONSTRAINT DESIGN")
    print("=" * 78)

    for gate, passed in gates.items():
        print(
            f"{gate:38s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Atom count:", len(output_rows))
    print("Fixed atoms:", len(fixed_rows))
    print("Mobile atoms:", len(mobile_rows))
    print(
        "Artificial caps:",
        len(artificial_cap_rows),
    )
    print(
        "Restored real atoms:",
        len(restored_real_rows),
    )
    print(
        "Fixed indices:",
        [
            row["index_0based"]
            for row in fixed_rows
        ],
    )

    print()
    print("Decision:", decision)
    print("Map:", OUTPUT_MAP)
    print("Report:", OUTPUT_REPORT)
    print(
        "ORCA input preparation authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
