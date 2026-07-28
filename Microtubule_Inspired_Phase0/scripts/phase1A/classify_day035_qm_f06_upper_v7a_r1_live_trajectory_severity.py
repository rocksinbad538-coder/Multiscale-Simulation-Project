#!/usr/bin/env python3

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import csv
import json


ROOT = Path(__file__).resolve().parents[2]

LIVE_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_live_trajectory"
)

LIVE_REPORT = (
    LIVE_DIR
    / "QM_F06_UPPER_V7A_R1_LIVE_TRAJECTORY_AUDIT.json"
)

LIVE_CSV = (
    LIVE_DIR
    / "QM_F06_UPPER_V7A_R1_live_trajectory_audit.csv"
)

OUTPUT_REPORT = (
    LIVE_DIR
    / "QM_F06_UPPER_V7A_R1_LIVE_TRAJECTORY_SEVERITY.json"
)

OUTPUT_CSV = (
    LIVE_DIR
    / "QM_F06_UPPER_V7A_R1_live_trajectory_severity.csv"
)


def parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def parse_int(value: str) -> int:
    return int(value.strip())


def parse_float_or_none(value: str) -> float | None:
    value = value.strip()

    if not value:
        return None

    return float(value)


def main() -> None:
    for path in (LIVE_REPORT, LIVE_CSV):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing live-trajectory artifact: {path}"
            )

    report = json.loads(
        LIVE_REPORT.read_text(
            encoding="utf-8"
        )
    )

    with LIVE_CSV.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    detailed_by_frame = {
        int(record["frame_index_0based"]): record
        for record in report.get(
            "violation_frames",
            [],
        )
    }

    severity_rows = []
    critical_frames = []
    warning_frames = []
    clean_frames = []

    for row in rows:
        frame_index = int(
            row["frame_index_0based"]
        )

        detailed = detailed_by_frame.get(
            frame_index,
            {},
        )

        nominal_bond_failures = (
            detailed.get(
                "nominal_bond_failures",
                [],
            )
        )

        hard_topology_failure = any((
            parse_int(
                row["lost_edge_count"]
            ) > 0,
            parse_int(
                row["gained_edge_count"]
            ) > 0,
            parse_int(
                row["degree_failure_count"]
            ) > 0,
            parse_int(
                row["connected_component_count"]
            ) != 1,
            parse_int(
                row["hard_contact_count"]
            ) > 0,
            not parse_bool(
                row["required_V7A_edges_present"]
            ),
            not parse_bool(
                row["V6B_forbidden_pair_absent"]
            ),
            not parse_bool(
                row["identity_pass"]
            ),
            not parse_bool(
                row["composition_pass"]
            ),
        ))

        nominal_window_warning = (
            len(nominal_bond_failures) > 0
            and not hard_topology_failure
        )

        if hard_topology_failure:
            severity = "CRITICAL_STRUCTURAL_VIOLATION"
            action = "REVIEW_CURRENT_EXECUTION"
            critical_frames.append(frame_index)

        elif nominal_window_warning:
            severity = "TRANSIENT_NOMINAL_BOND_WINDOW_WARNING"
            action = "CONTINUE_AND_MONITOR"
            warning_frames.append(frame_index)

        else:
            severity = "STRUCTURALLY_CLEAN"
            action = "CONTINUE"
            clean_frames.append(frame_index)

        maximum_nominal_excess = None
        limiting_bond = None

        for failure in nominal_bond_failures:
            distance_A = float(
                failure["distance_A"]
            )

            minimum_A = float(
                failure["minimum_A"]
            )

            maximum_A = float(
                failure["maximum_A"]
            )

            if distance_A > maximum_A:
                excess = distance_A - maximum_A
            elif distance_A < minimum_A:
                excess = minimum_A - distance_A
            else:
                excess = 0.0

            if (
                maximum_nominal_excess is None
                or excess > maximum_nominal_excess
            ):
                maximum_nominal_excess = excess
                limiting_bond = (
                    f"{failure['first_atom']}--"
                    f"{failure['second_atom']}"
                )

        severity_rows.append({
            "frame_index_0based": frame_index,
            "energy_Eh": row["energy_Eh"],
            "original_frame_pass": (
                row["frame_pass"]
            ),
            "nominal_bond_failure_count": (
                len(nominal_bond_failures)
            ),
            "lost_edge_count": (
                row["lost_edge_count"]
            ),
            "gained_edge_count": (
                row["gained_edge_count"]
            ),
            "degree_failure_count": (
                row["degree_failure_count"]
            ),
            "connected_component_count": (
                row["connected_component_count"]
            ),
            "hard_contact_count": (
                row["hard_contact_count"]
            ),
            "required_V7A_edges_present": (
                row["required_V7A_edges_present"]
            ),
            "V6B_forbidden_pair_absent": (
                row["V6B_forbidden_pair_absent"]
            ),
            "severity": severity,
            "recommended_action": action,
            "limiting_nominal_bond": (
                limiting_bond or ""
            ),
            "maximum_nominal_window_excess_A": (
                maximum_nominal_excess
                if maximum_nominal_excess is not None
                else ""
            ),
        })

    if critical_frames:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "LIVE_TRAJECTORY_CRITICAL_"
            "STRUCTURAL_REVIEW_REQUIRED"
        )
        continue_execution = False

    elif warning_frames:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "LIVE_TRAJECTORY_TRANSIENT_"
            "BOND_WINDOW_WARNING_"
            "CONTINUE_AND_MONITOR"
        )
        continue_execution = True

    else:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "LIVE_TRAJECTORY_STRUCTURALLY_CLEAN_"
            "CONTINUE"
        )
        continue_execution = True

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                severity_rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(severity_rows)

    output = {
        "model": "QM_F06_UPPER_V7A_R1",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "source_live_audit": str(
            LIVE_REPORT
        ),
        "decision": decision,
        "classification_policy": {
            "critical_structural_violation": (
                "Lost or gained geometric edge, incorrect "
                "coordination, fragmentation, hard contact, "
                "identity/composition failure, missing required "
                "V7-A edge, or reappearance of the V6-B "
                "forbidden pair."
            ),
            "transient_nominal_bond_window_warning": (
                "One or more nominal bonds are outside the "
                "preferred acceptance window while geometric "
                "connectivity, coordination, component integrity, "
                "hard-contact gates, and required-edge gates "
                "remain intact."
            ),
            "final_acceptance_policy": (
                "The final post-QM gate remains unchanged and "
                "requires every nominal bond to be inside its "
                "approved final acceptance window."
            ),
        },
        "summary": {
            "frame_count": len(rows),
            "clean_frames": clean_frames,
            "warning_frames": warning_frames,
            "critical_frames": critical_frames,
        },
        "authorizations": {
            "continue_current_ORCA_execution": (
                continue_execution
            ),
            "post_QM_acceptance_authorized": False,
            "RESP_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            output,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "LIVE TRAJECTORY SEVERITY CLASSIFICATION"
    )
    print("=" * 112)

    for row in severity_rows:
        print(
            f"frame={row['frame_index_0based']:>3} | "
            f"E={row['energy_Eh']} | "
            f"severity={row['severity']} | "
            f"bond={row['limiting_nominal_bond'] or 'none'} | "
            f"excess_A="
            f"{row['maximum_nominal_window_excess_A'] or 0}"
        )

    print()
    print("Clean frames:", clean_frames)
    print("Warning frames:", warning_frames)
    print("Critical frames:", critical_frames)

    print()
    print("Decision:", decision)
    print("CSV:", OUTPUT_CSV)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "Continue current ORCA execution:",
        continue_execution,
    )
    print(
        "Final post-QM acceptance authorized: False"
    )
    print("RESP authorized: False")
    print("MD authorized: False")

    if critical_frames:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
