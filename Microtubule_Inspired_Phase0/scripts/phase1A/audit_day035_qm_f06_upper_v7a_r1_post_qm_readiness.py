#!/usr/bin/env python3

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json
import re


ROOT = Path(__file__).resolve().parents[2]

EXEC_PARENT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_FILE = (
    EXEC_PARENT
    / "LATEST_V7A_R1_EXECUTION.txt"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_post_qm"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_POST_QM_READINESS.json"
)

EXPECTED_ATOM_COUNT = 52


def read_optional_text(path: Path) -> str:
    if not path.is_file():
        return ""

    return path.read_text(
        encoding="utf-8",
        errors="replace",
    )


def parse_integer_file(
    path: Path,
) -> int | None:
    if not path.is_file():
        return None

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).strip()

    try:
        return int(text)
    except ValueError:
        return None


def parse_single_xyz(
    path: Path,
) -> dict:
    result = {
        "present": False,
        "complete": False,
        "declared_atom_count": None,
        "parsed_atom_count": 0,
        "comment": "",
    }

    if not path.is_file() or path.stat().st_size == 0:
        return result

    result["present"] = True

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if not lines:
        return result

    try:
        declared = int(lines[0].strip())
    except ValueError:
        return result

    result["declared_atom_count"] = declared

    if len(lines) >= 2:
        result["comment"] = lines[1]

    parsed = 0

    for line in lines[2:2 + declared]:
        fields = line.split()

        if len(fields) < 4:
            break

        try:
            float(fields[1])
            float(fields[2])
            float(fields[3])
        except ValueError:
            break

        parsed += 1

    result["parsed_atom_count"] = parsed
    result["complete"] = (
        declared == EXPECTED_ATOM_COUNT
        and parsed == declared
    )

    return result


def parse_xyz_trajectory(
    path: Path,
) -> dict:
    result = {
        "present": False,
        "complete_frame_count": 0,
        "atom_counts": [],
        "comments": [],
        "incomplete_tail": False,
    }

    if not path.is_file() or path.stat().st_size == 0:
        return result

    result["present"] = True

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    cursor = 0

    while cursor < len(lines):
        try:
            atom_count = int(
                lines[cursor].strip()
            )
        except (ValueError, IndexError):
            result["incomplete_tail"] = True
            break

        frame_end = cursor + atom_count + 2

        if frame_end > len(lines):
            result["incomplete_tail"] = True
            break

        atom_lines = lines[
            cursor + 2:
            frame_end
        ]

        valid = True

        for line in atom_lines:
            fields = line.split()

            if len(fields) < 4:
                valid = False
                break

            try:
                float(fields[1])
                float(fields[2])
                float(fields[3])
            except ValueError:
                valid = False
                break

        if not valid:
            result["incomplete_tail"] = True
            break

        result["complete_frame_count"] += 1
        result["atom_counts"].append(
            atom_count
        )
        result["comments"].append(
            lines[cursor + 1]
            if cursor + 1 < len(lines)
            else ""
        )

        cursor = frame_end

    return result


def last_float_match(
    pattern: str,
    text: str,
) -> float | None:
    matches = re.findall(
        pattern,
        text,
        flags=(
            re.IGNORECASE
            | re.MULTILINE
        ),
    )

    if not matches:
        return None

    return float(matches[-1])


def main() -> None:
    if not LATEST_FILE.is_file():
        raise RuntimeError(
            f"Missing execution pointer: {LATEST_FILE}"
        )

    execution_relative = (
        LATEST_FILE
        .read_text(encoding="utf-8")
        .strip()
    )

    execution_dir = (
        ROOT
        / execution_relative
    )

    output_path = (
        execution_dir
        / "v7a_r1.out"
    )

    final_xyz_path = (
        execution_dir
        / "v7a_r1.xyz"
    )

    trajectory_path = (
        execution_dir
        / "v7a_r1_trj.xyz"
    )

    opt_path = (
        execution_dir
        / "v7a_r1.opt"
    )

    stderr_path = (
        execution_dir
        / "v7a_r1.stderr"
    )

    exit_status_path = (
        execution_dir
        / "v7a_r1.exit_status"
    )

    shell_status_path = (
        execution_dir
        / "v7a_r1.orca_shell_status"
    )

    classification_path = (
        execution_dir
        / "v7a_r1.termination_classification"
    )

    output_text = read_optional_text(
        output_path
    )

    stderr_text = read_optional_text(
        stderr_path
    )

    exit_status = parse_integer_file(
        exit_status_path
    )

    shell_status = parse_integer_file(
        shell_status_path
    )

    classification = (
        read_optional_text(
            classification_path
        ).strip()
    )

    normal_termination_count = (
        output_text.count(
            "ORCA TERMINATED NORMALLY"
        )
    )

    error_termination_count = (
        output_text.count(
            "ORCA finished by error termination"
        )
    )

    optimization_converged_count = (
        output_text.count(
            "THE OPTIMIZATION HAS CONVERGED"
        )
    )

    geometry_cycle_count = len(
        re.findall(
            r"GEOMETRY OPTIMIZATION CYCLE",
            output_text,
        )
    )

    scf_convergence_count = len(
        re.findall(
            r"SCF CONVERGED AFTER",
            output_text,
        )
    )

    final_energy = last_float_match(
        r"FINAL SINGLE POINT ENERGY\s+"
        r"(-?\d+\.\d+)",
        output_text,
    )

    final_geometry = parse_single_xyz(
        final_xyz_path
    )

    trajectory = parse_xyz_trajectory(
        trajectory_path
    )

    running = (
        exit_status is None
        and classification == ""
    )

    terminated = not running

    valid_classified_completion = (
        exit_status == 0
        and shell_status == 0
        and classification
        == "ORCA_NORMAL_TERMINATION"
    )

    formal_post_qm_ready = all((
        terminated,
        valid_classified_completion,
        normal_termination_count >= 1,
        error_termination_count == 0,
        optimization_converged_count >= 1,
        final_energy is not None,
        final_geometry["complete"],
        trajectory["present"],
        (
            trajectory[
                "complete_frame_count"
            ]
            >= 1
        ),
        all(
            count == EXPECTED_ATOM_COUNT
            for count in trajectory[
                "atom_counts"
            ]
        ),
        opt_path.is_file(),
        opt_path.stat().st_size > 0,
        stderr_text.strip() == "",
    ))

    if running:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "POST_QM_NOT_READY_EXECUTION_RUNNING"
        )
    elif formal_post_qm_ready:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "POST_QM_READINESS_GATE_PASS_"
            "STRUCTURAL_AUDIT_AUTHORIZED"
        )
    else:
        decision = (
            "QM_F06_UPPER_V7A_R1_"
            "POST_QM_READINESS_GATE_FAIL_"
            "EXECUTION_REVIEW_REQUIRED"
        )

    gates = {
        "execution_finished": terminated,
        "classified_status_zero": (
            exit_status == 0
        ),
        "ORCA_shell_status_zero": (
            shell_status == 0
        ),
        "termination_classification_normal": (
            classification
            == "ORCA_NORMAL_TERMINATION"
        ),
        "normal_ORCA_termination_marker": (
            normal_termination_count >= 1
        ),
        "no_ORCA_error_termination_marker": (
            error_termination_count == 0
        ),
        "optimization_convergence_marker": (
            optimization_converged_count
            >= 1
        ),
        "final_energy_available": (
            final_energy is not None
        ),
        "final_XYZ_complete_52_atoms": (
            final_geometry["complete"]
        ),
        "trajectory_available": (
            trajectory["present"]
        ),
        "trajectory_has_complete_frames": (
            trajectory[
                "complete_frame_count"
            ]
            >= 1
        ),
        "all_complete_frames_have_52_atoms": (
            bool(trajectory["atom_counts"])
            and all(
                count == EXPECTED_ATOM_COUNT
                for count in trajectory[
                    "atom_counts"
                ]
            )
        ),
        "optimization_file_available": (
            opt_path.is_file()
            and opt_path.stat().st_size > 0
        ),
        "stderr_empty": (
            stderr_text.strip() == ""
        ),
    }

    report = {
        "model": "QM_F06_UPPER_V7A_R1",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir
        ),
        "decision": decision,
        "execution_state": (
            "RUNNING"
            if running
            else "FINISHED"
        ),
        "gates": gates,
        "summary": {
            "exit_status": exit_status,
            "ORCA_shell_status": shell_status,
            "termination_classification": (
                classification or None
            ),
            "normal_termination_markers": (
                normal_termination_count
            ),
            "error_termination_markers": (
                error_termination_count
            ),
            "optimization_convergence_markers": (
                optimization_converged_count
            ),
            "geometry_cycle_count": (
                geometry_cycle_count
            ),
            "SCF_convergence_count": (
                scf_convergence_count
            ),
            "latest_final_energy_Eh": (
                final_energy
            ),
            "final_geometry": final_geometry,
            "trajectory": trajectory,
            "stderr_character_count": len(
                stderr_text
            ),
        },
        "authorizations": {
            "post_QM_structural_audit_authorized": (
                formal_post_qm_ready
            ),
            "RESP_input_preparation_authorized": (
                False
            ),
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": (
                False
            ),
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "POST-QM READINESS GATE"
    )
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:58s}: "
            f"{'PASS' if value else 'WAIT/FAIL'}"
        )

    print()
    print(
        "Execution state:",
        report["execution_state"],
    )
    print("Exit status:", exit_status)
    print(
        "ORCA shell status:",
        shell_status,
    )
    print(
        "Termination classification:",
        classification or "PENDING",
    )
    print(
        "Geometry cycles:",
        geometry_cycle_count,
    )
    print(
        "SCF convergences:",
        scf_convergence_count,
    )
    print(
        "Latest energy Eh:",
        final_energy,
    )
    print(
        "Final XYZ present:",
        final_geometry["present"],
    )
    print(
        "Final XYZ complete:",
        final_geometry["complete"],
    )
    print(
        "Trajectory complete frames:",
        trajectory[
            "complete_frame_count"
        ],
    )
    print(
        "Trajectory incomplete tail:",
        trajectory["incomplete_tail"],
    )

    print()
    print("Decision:", decision)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "Post-QM structural audit authorized:",
        formal_post_qm_ready,
    )
    print(
        "RESP input preparation authorized: False"
    )
    print("RESP execution authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
