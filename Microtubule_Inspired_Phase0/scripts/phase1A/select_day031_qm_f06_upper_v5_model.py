#!/usr/bin/env python3
"""
Select the chemically preferred QM_F06 UPPER V5 model.

The selected expansion:
- restores BR4:UPPER:14:1;
- restores P:1637;
- restores S:1737;
- internalizes both problematic S:1738 cuts;
- replaces the P:1637 double-cut boundary by a simple S:1737 cut;
- introduces only simple boundary cuts around the new region.

This script authorizes coordinate-source validation, not construction
or ORCA execution.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SUBSET_REPORT = ROOT / (
    "runs/phase1A/"
    "day031_qm_f06_upper_v5_selective_subsets/"
    "QM_F06_UPPER_V5_SELECTIVE_SUBSET_AUDIT.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day031_qm_f06_upper_v5_model_selection"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_MODEL_SELECTION.json"
)

SELECTED_SUBSET = (
    "BR4:UPPER:14:1+P:1637+S:1737"
)

SELECTED_ATOMS = {
    "BR4:UPPER:14:1",
    "P:1637",
    "S:1737",
}

EXPECTED_METRICS = {
    "added_atom_count": 3,
    "boundary_cut_count": 14,
    "cuts_closed_vs_v4": 2,
    "new_cuts_vs_v4": 3,
    "target_cuts_remaining": 0,
    "multi_cut_boundary_center_count": 1,
}

CLOSED_TARGET_CUTS = [
    ["S:1738", "BR4:UPPER:14:1"],
    ["S:1738", "P:1637"],
]

NEW_BOUNDARY_CUTS = [
    ["BR4:UPPER:14:1", "BR4:UPPER:14:2"],
    ["P:1637", "P:1636"],
    ["S:1737", "P:1635"],
]


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


def main() -> None:
    require_file(SUBSET_REPORT)

    data = json.loads(
        SUBSET_REPORT.read_text(
            encoding="utf-8"
        )
    )

    candidates = data["top_ten"]

    selected = next(
        (
            row
            for row in candidates
            if row["subset_id"] == SELECTED_SUBSET
        ),
        None,
    )

    if selected is None:
        raise RuntimeError(
            "Selected V5 subset was not found "
            "in the audited top-ten models."
        )

    gates = {
        "subset_identity": (
            set(selected["added_atoms"].split("|"))
            == SELECTED_ATOMS
        ),
        "added_atom_count": (
            selected["added_atom_count"]
            == EXPECTED_METRICS["added_atom_count"]
        ),
        "boundary_cut_count": (
            selected["boundary_cut_count"]
            == EXPECTED_METRICS["boundary_cut_count"]
        ),
        "cuts_closed_vs_v4": (
            selected["cuts_closed_vs_v4"]
            == EXPECTED_METRICS["cuts_closed_vs_v4"]
        ),
        "new_cuts_vs_v4": (
            selected["new_cuts_vs_v4"]
            == EXPECTED_METRICS["new_cuts_vs_v4"]
        ),
        "all_target_cuts_closed": (
            selected["target_cuts_remaining"]
            == EXPECTED_METRICS[
                "target_cuts_remaining"
            ]
        ),
        "single_multicut_center_remaining": (
            selected[
                "multi_cut_boundary_center_count"
            ]
            == EXPECTED_METRICS[
                "multi_cut_boundary_center_count"
            ]
        ),
        "remaining_multicut_is_preexisting": (
            selected["multi_cut_boundary_centers"]
            == "A:UPPER:8:4"
        ),
    }

    overall_pass = all(gates.values())

    if not overall_pass:
        decision = (
            "QM_F06_UPPER_V5_MODEL_SELECTION_FAIL"
        )
    else:
        decision = (
            "QM_F06_UPPER_V5_SELECTIVE_MODEL_SELECTED_"
            "COORDINATE_VALIDATION_AUTHORIZED"
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
        "selected_subset": SELECTED_SUBSET,
        "selected_real_atoms": sorted(
            SELECTED_ATOMS
        ),
        "selected_metrics": selected,
        "closed_target_cuts": CLOSED_TARGET_CUTS,
        "new_boundary_cuts": NEW_BOUNDARY_CUTS,
        "selection_basis": {
            "all_problematic_S1738_cuts_internalized": True,
            "P1637_double_cut_removed": True,
            "new_region_multicut_centers": 0,
            "only_remaining_multicut_center": (
                "A:UPPER:8:4"
            ),
            "remaining_multicut_preexists_V5": True,
        },
        "gates": gates,
        "overall_pass": overall_pass,
        "source_report": str(
            SUBSET_REPORT.relative_to(ROOT)
        ),
        "source_report_sha256": sha256(
            SUBSET_REPORT
        ),
        "authorization": {
            "v5_model_selected": overall_pass,
            "coordinate_source_validation_authorized": (
                overall_pass
            ),
            "v5_geometry_construction_authorized": False,
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 90)
    print("QM_F06 UPPER V5 MODEL SELECTION")
    print("=" * 90)

    for gate, passed in gates.items():
        print(
            f"{gate:44s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Selected subset:", SELECTED_SUBSET)
    print(
        "Selected real atoms:",
        sorted(SELECTED_ATOMS),
    )
    print(
        "Boundary cuts:",
        selected["boundary_cut_count"],
    )
    print(
        "Target cuts remaining:",
        selected["target_cuts_remaining"],
    )
    print(
        "New-region multicut centers:",
        0,
    )
    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print()
    print(
        "Coordinate validation authorized:",
        overall_pass,
    )
    print("V5 construction authorized: False")
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
