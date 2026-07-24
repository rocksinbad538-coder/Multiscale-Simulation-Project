#!/usr/bin/env python3
"""
Formal strategy selection for QM_F06 UPPER V5-B.

V5-B combines:
- the selected real-atom expansion:
    BR4:UPPER:14:1
    P:1637
    S:1737
- removal of the two superseded S:1738 artificial caps;
- replacement passivation at the three new external cuts;
- release of A:UPPER:14:4 to eliminate the fixed/mobile
  asymmetry across canonical edge E:2915;
- reconstruction of H4:UPPER:0203:0 away from fixed N P:1641;
- retention of P:1641 as a fixed core atom.

This script authorizes V5-B construction only. It does not construct
the geometry or authorize ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

MODEL_SELECTION = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_model_selection/"
    "QM_F06_UPPER_V5_MODEL_SELECTION.json"
)

COORDINATE_REPORT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_coordinate_sources/"
    "QM_F06_UPPER_V5_COORDINATE_SOURCE_AUDIT.json"
)

V4_CONSTRAINT_MAP = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_constraint_design/"
    "QM_F06_UPPER_V4_constraint_map.csv"
)

V4_CAPS = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction/"
    "QM_F06_UPPER_V4_artificial_caps.csv"
)

V4_TERMINATION_POINTER = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_executions/"
    "LATEST_V4_EXECUTION.txt"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_strategy"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_STRATEGY.json"
)

SELECTED_REAL_ATOMS = {
    "BR4:UPPER:14:1",
    "P:1637",
    "S:1737",
}

REMOVED_CAPS = {
    "HCAPV4:UPPER:S1738:BR4_14_1",
    "HCAPV4:UPPER:S1738:P1637",
}

NEW_BOUNDARY_CUTS = {
    ("BR4:UPPER:14:1", "BR4:UPPER:14:2"),
    ("P:1637", "P:1636"),
    ("S:1737", "P:1635"),
}

RELEASED_V4_FIXED_ATOMS = {
    "A:UPPER:14:4",
}

REBUILT_PASSIVANT = "H4:UPPER:0203:0"
REBUILT_PASSIVANT_CENTER = "S:1739"
AVOIDED_HEAVY_ATOM = "P:1641"


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


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def main() -> None:
    for path in (
        MODEL_SELECTION,
        COORDINATE_REPORT,
        V4_CONSTRAINT_MAP,
        V4_CAPS,
        V4_TERMINATION_POINTER,
    ):
        require_file(path)

    model = json.loads(
        MODEL_SELECTION.read_text(encoding="utf-8")
    )

    coordinates = json.loads(
        COORDINATE_REPORT.read_text(encoding="utf-8")
    )

    constraint_rows = read_csv(V4_CONSTRAINT_MAP)
    cap_rows = read_csv(V4_CAPS)

    latest_execution = (
        ROOT
        / Path(
            V4_TERMINATION_POINTER.read_text(
                encoding="utf-8"
            ).strip()
        )
    )

    termination_record = (
        latest_execution
        / "QM_F06_UPPER_V4_TERMINATION_RECORD.json"
    )

    require_file(termination_record)

    termination = json.loads(
        termination_record.read_text(encoding="utf-8")
    )

    constraint_by_id = {
        row["atom_id"]: row
        for row in constraint_rows
    }

    cap_ids = {
        row["cap_id"]
        for row in cap_rows
    }

    gates = {
        "v5_model_selected": (
            model["authorization"]["v5_model_selected"]
        ),
        "selected_real_atoms_match": (
            set(model["selected_real_atoms"])
            == SELECTED_REAL_ATOMS
        ),
        "transformed_coordinates_authorized": (
            coordinates["authorization"][
                "transformed_coordinates_authorized"
            ]
        ),
        "all_selected_atoms_have_coordinates": (
            set(coordinates["target_atoms"])
            == SELECTED_REAL_ATOMS
        ),
        "superseded_caps_exist_in_v4": (
            REMOVED_CAPS.issubset(cap_ids)
        ),
        "E2915_mobile_endpoint_confirmed": (
            not parse_bool(
                constraint_by_id[
                    "A:UPPER:13:3"
                ]["v4_fixed"]
            )
        ),
        "E2915_fixed_endpoint_confirmed": (
            parse_bool(
                constraint_by_id[
                    "A:UPPER:14:4"
                ]["v4_fixed"]
            )
        ),
        "P1641_fixed_core_confirmed": (
            parse_bool(
                constraint_by_id[
                    "P:1641"
                ]["v4_fixed"]
            )
        ),
        "S1739_mobile_confirmed": (
            parse_bool(
                constraint_by_id[
                    "S:1739"
                ]["v4_mobile"]
            )
        ),
        "H0203_mobile_confirmed": (
            parse_bool(
                constraint_by_id[
                    REBUILT_PASSIVANT
                ]["v4_mobile"]
            )
        ),
        "v4_failure_record_confirmed": (
            termination["decision"]
            == (
                "QM_F06_UPPER_V4_MANUALLY_TERMINATED_"
                "BOUNDARY_CHEMICAL_REORGANIZATION"
            )
        ),
        "canonical_E2915_rupture_confirmed": (
            termination[
                "canonical_bond_rupture"
            ]["canonical_edge_confirmed"]
        ),
        "hydrogen_transfer_confirmed": (
            termination[
                "hydrogen_transfer"
            ]["hydrogen"]
            == REBUILT_PASSIVANT
            and termination[
                "hydrogen_transfer"
            ]["final_owner"]
            == AVOIDED_HEAVY_ATOM
        ),
        "fixed_core_constraints_respected": (
            termination["fixed_core"][
                "constraints_respected"
            ]
        ),
    }

    overall_pass = all(gates.values())

    decision = (
        "QM_F06_UPPER_V5B_STRATEGY_PASS_"
        "GEOMETRY_CONSTRUCTION_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V5B_STRATEGY_FAIL_"
        "CONSTRUCTION_BLOCKED"
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "strategy": "V5-B",
        "selected_real_atoms": sorted(
            SELECTED_REAL_ATOMS
        ),
        "removed_superseded_caps": sorted(
            REMOVED_CAPS
        ),
        "new_boundary_cuts_to_passivate": [
            list(pair)
            for pair in sorted(
                tuple(sorted(pair))
                for pair in NEW_BOUNDARY_CUTS
            )
        ],
        "constraint_revision": {
            "release_from_v4_fixed_core": sorted(
                RELEASED_V4_FIXED_ATOMS
            ),
            "retain_fixed": [AVOIDED_HEAVY_ATOM],
            "reason": (
                "Eliminate the asymmetric fixed/mobile "
                "partition across canonical edge E:2915."
            ),
        },
        "passivant_repair": {
            "atom_id": REBUILT_PASSIVANT,
            "bonded_center": (
                REBUILT_PASSIVANT_CENTER
            ),
            "avoid_atom": AVOIDED_HEAVY_ATOM,
            "operation": (
                "Reconstruct the B-H direction from the "
                "local heavy-neighbor geometry so that it "
                "points away from P:1641."
            ),
        },
        "gates": gates,
        "overall_pass": overall_pass,
        "sources": {
            "model_selection": str(
                MODEL_SELECTION.relative_to(ROOT)
            ),
            "coordinate_report": str(
                COORDINATE_REPORT.relative_to(ROOT)
            ),
            "v4_constraint_map": str(
                V4_CONSTRAINT_MAP.relative_to(ROOT)
            ),
            "v4_caps": str(
                V4_CAPS.relative_to(ROOT)
            ),
            "v4_termination_record": str(
                termination_record.relative_to(ROOT)
            ),
        },
        "source_sha256": {
            "model_selection": sha256(
                MODEL_SELECTION
            ),
            "coordinate_report": sha256(
                COORDINATE_REPORT
            ),
            "v4_constraint_map": sha256(
                V4_CONSTRAINT_MAP
            ),
            "v4_caps": sha256(V4_CAPS),
            "v4_termination_record": sha256(
                termination_record
            ),
        },
        "authorization": {
            "v5b_strategy_selected": overall_pass,
            "v5b_geometry_construction_authorized": (
                overall_pass
            ),
            "pre_qm_structural_audit_authorized": False,
            "constraint_design_authorized": False,
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 94)
    print("QM_F06 UPPER V5-B STRATEGY")
    print("=" * 94)

    for gate, passed in gates.items():
        print(
            f"{gate:48s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Selected real atoms:", sorted(
        SELECTED_REAL_ATOMS
    ))
    print("Removed caps:", sorted(REMOVED_CAPS))
    print(
        "Released V4 fixed atoms:",
        sorted(RELEASED_V4_FIXED_ATOMS),
    )
    print(
        "Rebuilt passivant:",
        REBUILT_PASSIVANT,
    )
    print("Retained fixed N:", AVOIDED_HEAVY_ATOM)
    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print()
    print(
        "V5-B construction authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
