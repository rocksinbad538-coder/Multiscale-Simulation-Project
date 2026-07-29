#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RUN_ROOT = (
    PROJECT_ROOT
    / "runs"
    / "phase1A"
)

TRAJECTORY_ROOT = (
    RUN_ROOT
    / "day035_qm_f06_upper_v7a_r1_live_trajectory"
)

SEVERITY_CSV = (
    TRAJECTORY_ROOT
    / "QM_F06_UPPER_V7A_R1_live_trajectory_severity.csv"
)

EXECUTION_PARENT = (
    RUN_ROOT
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_V7A_R1_EXECUTION.txt"
)

OUTPUT_ROOT = (
    RUN_ROOT
    / "day035_qm_f06_upper_v7a_r1_execution_continuation"
)

OUTPUT_JSON = (
    OUTPUT_ROOT
    / "QM_F06_UPPER_V7A_R1_EXECUTION_CONTINUATION.json"
)

REQUIRED_CONSECUTIVE_CLEAN_FRAMES = 3


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def parse_int(value: str) -> int:
    return int(float(value or 0))


def main() -> None:
    if not SEVERITY_CSV.is_file():
        raise RuntimeError(
            f"Missing severity CSV: {SEVERITY_CSV}"
        )

    if not LATEST_POINTER.is_file():
        raise RuntimeError(
            f"Missing execution pointer: {LATEST_POINTER}"
        )

    execution_relative = (
        LATEST_POINTER.read_text().strip()
    )

    if not execution_relative:
        raise RuntimeError(
            "Latest execution pointer is empty"
        )

    execution_directory = (
        PROJECT_ROOT / execution_relative
    )

    exit_status_path = (
        execution_directory / "v7a_r1.exit_status"
    )

    execution_running = not (
        exit_status_path.is_file()
        and exit_status_path.stat().st_size > 0
    )

    with SEVERITY_CSV.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            "Severity CSV contains no trajectory frames"
        )

    rows.sort(
        key=lambda row: int(
            row["frame_index_0based"]
        )
    )

    latest = rows[-1]

    clean_tail = []

    for row in reversed(rows):
        if row["severity"] != "STRUCTURALLY_CLEAN":
            break

        clean_tail.append(
            int(row["frame_index_0based"])
        )

    clean_tail.reverse()

    latest_frame_clean = (
        latest["severity"] == "STRUCTURALLY_CLEAN"
    )

    latest_no_lost_edges = (
        parse_int(latest["lost_edge_count"]) == 0
    )

    latest_no_gained_edges = (
        parse_int(latest["gained_edge_count"]) == 0
    )

    latest_no_degree_failures = (
        parse_int(latest["degree_failure_count"]) == 0
    )

    latest_single_component = (
        parse_int(
            latest["connected_component_count"]
        )
        == 1
    )

    latest_no_hard_contacts = (
        parse_int(latest["hard_contact_count"]) == 0
    )

    latest_required_v7a_edges_present = (
        parse_bool(
            latest["required_V7A_edges_present"]
        )
    )

    latest_v6b_forbidden_pair_absent = (
        parse_bool(
            latest["V6B_forbidden_pair_absent"]
        )
    )

    sufficient_clean_tail = (
        len(clean_tail)
        >= REQUIRED_CONSECUTIVE_CLEAN_FRAMES
    )

    historical_critical_frames = [
        int(row["frame_index_0based"])
        for row in rows
        if row["severity"]
        == "CRITICAL_STRUCTURAL_VIOLATION"
    ]

    historical_warning_frames = [
        int(row["frame_index_0based"])
        for row in rows
        if row["severity"]
        == "TRANSIENT_NOMINAL_BOND_WINDOW_WARNING"
    ]

    continuation_checks = {
        "execution_running": execution_running,
        "latest_frame_clean": latest_frame_clean,
        "latest_no_lost_edges": latest_no_lost_edges,
        "latest_no_gained_edges": latest_no_gained_edges,
        "latest_no_degree_failures": (
            latest_no_degree_failures
        ),
        "latest_single_component": (
            latest_single_component
        ),
        "latest_no_hard_contacts": (
            latest_no_hard_contacts
        ),
        "latest_required_V7A_edges_present": (
            latest_required_v7a_edges_present
        ),
        "latest_V6B_forbidden_pair_absent": (
            latest_v6b_forbidden_pair_absent
        ),
        "sufficient_consecutive_clean_frames": (
            sufficient_clean_tail
        ),
    }

    continue_execution = all(
        continuation_checks.values()
    )

    report = {
        "system": "QM_F06_UPPER_V7A_R1",
        "purpose": (
            "Live execution-continuation assessment. "
            "This report does not authorize final "
            "post-QM acceptance, RESP, or MD."
        ),
        "execution_directory": str(
            execution_directory
        ),
        "execution_running": execution_running,
        "complete_frame_count": len(rows),
        "latest_frame_index_0based": int(
            latest["frame_index_0based"]
        ),
        "latest_energy_Eh": float(
            latest["energy_Eh"]
        ),
        "latest_severity": latest["severity"],
        "required_consecutive_clean_frames": (
            REQUIRED_CONSECUTIVE_CLEAN_FRAMES
        ),
        "consecutive_clean_tail_frames": clean_tail,
        "historical_warning_frames": (
            historical_warning_frames
        ),
        "historical_critical_frames": (
            historical_critical_frames
        ),
        "continuation_checks": (
            continuation_checks
        ),
        "continue_current_ORCA_execution": (
            continue_execution
        ),
        "post_QM_acceptance_authorized": False,
        "RESP_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_JSON.write_text(
        json.dumps(
            report,
            indent=2,
        )
        + "\n"
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "EXECUTION-CONTINUATION GATE"
    )
    print("=" * 112)

    for key, value in continuation_checks.items():
        status = "PASS" if value else "FAIL"
        print(f"{key:52s}: {status}")

    print()
    print(
        "Latest complete frame: "
        f"{report['latest_frame_index_0based']}"
    )
    print(
        "Latest energy Eh: "
        f"{report['latest_energy_Eh']:.12f}"
    )
    print(
        "Consecutive clean tail: "
        f"{clean_tail}"
    )
    print(
        "Historical warning frames: "
        f"{historical_warning_frames}"
    )
    print(
        "Historical critical frames: "
        f"{historical_critical_frames}"
    )
    print()

    if continue_execution:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "CURRENT_EXECUTION_CONTINUE"
        )
    else:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "CURRENT_EXECUTION_REVIEW_REQUIRED"
        )

    print(f"Decision: {decision}")
    print(f"Report: {OUTPUT_JSON}")
    print()
    print(
        "Continue current ORCA execution: "
        f"{continue_execution}"
    )
    print(
        "Post-QM acceptance authorized: False"
    )
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
