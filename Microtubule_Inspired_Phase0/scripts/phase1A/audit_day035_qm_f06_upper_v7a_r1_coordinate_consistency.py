#!/usr/bin/env python3

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import math


ROOT = Path(__file__).resolve().parents[2]

POST_QM_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_post_qm"
)

POST_QM_REPORT = (
    POST_QM_DIR
    / "QM_F06_UPPER_V7A_R1_POST_QM_AUDIT.json"
)

POST_QM_FINAL_XYZ = (
    POST_QM_DIR
    / "QM_F06_UPPER_V7A_R1_FINAL.xyz"
)

EXEC_PARENT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_EXECUTION_FILE = (
    EXEC_PARENT
    / "LATEST_V7A_R1_EXECUTION.txt"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_coordinate_consistency"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_COORDINATE_CONSISTENCY.json"
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_coordinate_differences.csv"
)

EXPECTED_POST_QM_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "POST_QM_GATE_PASS_"
    "RESP_INPUT_PREPARATION_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52

COORDINATE_IDENTITY_TOLERANCE_A = 5.0e-12

PASS_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "COORDINATE_CONSISTENCY_PASS_"
    "COORDINATE_ADOPTION_AUTHORIZED"
)

FAIL_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "COORDINATE_CONSISTENCY_FAIL_"
    "REVIEW_REQUIRED"
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


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="strict",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"XYZ file is incomplete: {path}"
        )

    try:
        atom_count = int(
            lines[0].strip()
        )
    except ValueError as error:
        raise RuntimeError(
            f"Invalid XYZ atom count: {path}"
        ) from error

    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ atom record "
                f"at index {index}: {path}"
            )

        coordinates = tuple(
            float(value)
            for value in fields[1:4]
        )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": coordinates,
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"XYZ declared {atom_count} atoms "
            f"but only {len(atoms)} were parsed."
        )

    return atoms


def read_trajectory(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    frames = []

    cursor = 0

    while cursor < len(lines):

        try:
            atom_count = int(
                lines[cursor].strip()
            )
        except (ValueError, IndexError):
            break

        end = cursor + atom_count + 2

        if end > len(lines):
            break

        atoms = []

        for index, line in enumerate(
            lines[cursor + 2:end]
        ):
            fields = line.split()

            if len(fields) < 4:
                atoms = []
                break

            atoms.append({
                "index_0based": index,
                "element": fields[0],
                "xyz_A": tuple(
                    float(value)
                    for value in fields[1:4]
                ),
            })

        if len(atoms) != atom_count:
            break

        frames.append({
            "frame_index_0based": len(frames),
            "comment": lines[cursor + 1],
            "atoms": atoms,
        })

        cursor = end

    return frames


def coordinate_difference(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:

    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def write_csv(
    path: Path,
    rows: list[dict],
) -> None:

    if not rows:
        path.write_text(
            "",
            encoding="utf-8",
        )
        return

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
        )

        writer.writeheader()
        writer.writerows(rows)

def read_last_orca_coordinates(
    path: Path,
) -> list[dict]:

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    header = (
        "CARTESIAN COORDINATES (ANGSTROEM)"
    )

    start_indices = [
        index
        for index, line in enumerate(lines)
        if header in line
    ]

    if not start_indices:
        raise RuntimeError(
            "No ORCA coordinate blocks were found."
        )

    last_atoms = None

    for start in start_indices:

        cursor = start + 1

        #
        # Skip blank lines and dashed separators
        #
        while (
            cursor < len(lines)
            and (
                not lines[cursor].strip()
                or set(lines[cursor].strip()) == {"-"}
            )
        ):
            cursor += 1

        atoms = []

        while cursor < len(lines):

            line = lines[cursor].strip()

            if not line:
                break

            fields = line.split()

            #
            # Typical coordinate line:
            #
            # B   -0.123456   1.234567   2.345678
            #
            if (
                len(fields) != 4
            ):
                break

            try:

                xyz = tuple(
                    float(value)
                    for value in fields[1:4]
                )

            except ValueError:
                break

            atoms.append({
                "index_0based": len(atoms),
                "element": fields[0],
                "xyz_A": xyz,
            })

            cursor += 1

        if len(atoms) == EXPECTED_ATOM_COUNT:
            last_atoms = atoms

    if last_atoms is None:
        raise RuntimeError(
            "No complete ORCA coordinate block "
            f"containing {EXPECTED_ATOM_COUNT} atoms "
            "was found."
        )

    return last_atoms

def compare_coordinate_sets(
    reference_name: str,
    candidate_name: str,
    reference_atoms: list[dict],
    candidate_atoms: list[dict],
) -> dict:

    result = {
        "reference_name": reference_name,
        "candidate_name": candidate_name,
        "reference_atom_count": len(reference_atoms),
        "candidate_atom_count": len(candidate_atoms),
        "atom_count_match": False,
        "element_order_match": False,
        "rmsd_A": None,
        "maximum_displacement_A": None,
        "maximum_displacement_atom_index_0based": None,
        "maximum_displacement_element": None,
        "maximum_displacement_dx_A": None,
        "maximum_displacement_dy_A": None,
        "maximum_displacement_dz_A": None,
        "coordinate_identity_tolerance_A": (
            COORDINATE_IDENTITY_TOLERANCE_A
        ),
        "coordinate_identity_pass": False,
        "difference_rows": [],
        "failure_reasons": [],
    }

    if len(reference_atoms) != len(candidate_atoms):
        result["failure_reasons"].append(
            "atom_count_mismatch"
        )
        return result

    result["atom_count_match"] = True

    reference_elements = [
        atom["element"]
        for atom in reference_atoms
    ]

    candidate_elements = [
        atom["element"]
        for atom in candidate_atoms
    ]

    if reference_elements != candidate_elements:
        result["failure_reasons"].append(
            "element_order_mismatch"
        )
        return result

    result["element_order_match"] = True

    squared_displacements = []
    maximum_row = None

    for reference_atom, candidate_atom in zip(
        reference_atoms,
        candidate_atoms,
    ):
        index = reference_atom["index_0based"]
        element = reference_atom["element"]

        reference_xyz = reference_atom["xyz_A"]
        candidate_xyz = candidate_atom["xyz_A"]

        dx = candidate_xyz[0] - reference_xyz[0]
        dy = candidate_xyz[1] - reference_xyz[1]
        dz = candidate_xyz[2] - reference_xyz[2]

        displacement = coordinate_difference(
            reference_xyz,
            candidate_xyz,
        )

        squared_displacements.append(
            displacement ** 2
        )

        row = {
            "comparison": (
                f"{reference_name}_vs_{candidate_name}"
            ),
            "reference_name": reference_name,
            "candidate_name": candidate_name,
            "atom_index_0based": index,
            "element": element,
            "reference_x_A": reference_xyz[0],
            "reference_y_A": reference_xyz[1],
            "reference_z_A": reference_xyz[2],
            "candidate_x_A": candidate_xyz[0],
            "candidate_y_A": candidate_xyz[1],
            "candidate_z_A": candidate_xyz[2],
            "dx_A": dx,
            "dy_A": dy,
            "dz_A": dz,
            "displacement_A": displacement,
            "within_identity_tolerance": (
                displacement
                <= COORDINATE_IDENTITY_TOLERANCE_A
            ),
        }

        result["difference_rows"].append(
            row
        )

        if (
            maximum_row is None
            or displacement
            > maximum_row["displacement_A"]
        ):
            maximum_row = row

    rmsd = math.sqrt(
        sum(squared_displacements)
        / len(squared_displacements)
    )

    maximum_displacement = (
        maximum_row["displacement_A"]
        if maximum_row is not None
        else 0.0
    )

    result["rmsd_A"] = rmsd
    result["maximum_displacement_A"] = (
        maximum_displacement
    )

    if maximum_row is not None:
        result[
            "maximum_displacement_atom_index_0based"
        ] = maximum_row["atom_index_0based"]

        result[
            "maximum_displacement_element"
        ] = maximum_row["element"]

        result[
            "maximum_displacement_dx_A"
        ] = maximum_row["dx_A"]

        result[
            "maximum_displacement_dy_A"
        ] = maximum_row["dy_A"]

        result[
            "maximum_displacement_dz_A"
        ] = maximum_row["dz_A"]

    result["coordinate_identity_pass"] = (
        maximum_displacement
        <= COORDINATE_IDENTITY_TOLERANCE_A
    )

    if not result["coordinate_identity_pass"]:
        result["failure_reasons"].append(
            "coordinate_identity_tolerance_exceeded"
        )

    return result

def main() -> None:
    generated_at_utc = datetime.now(
        timezone.utc
    ).isoformat()

    if not POST_QM_REPORT.is_file():
        raise RuntimeError(
            "Post-QM structural audit report is not available: "
            f"{POST_QM_REPORT}"
        )

    if (
        not POST_QM_FINAL_XYZ.is_file()
        or POST_QM_FINAL_XYZ.stat().st_size == 0
    ):
        raise RuntimeError(
            "Post-QM audited final XYZ is not available: "
            f"{POST_QM_FINAL_XYZ}"
        )

    if not LATEST_EXECUTION_FILE.is_file():
        raise RuntimeError(
            "Latest R1 execution pointer is missing: "
            f"{LATEST_EXECUTION_FILE}"
        )

    post_qm = json.loads(
        POST_QM_REPORT.read_text(
            encoding="utf-8"
        )
    )

    post_qm_decision = post_qm.get(
        "decision"
    )

    post_qm_authorized = (
        post_qm.get(
            "authorizations",
            {},
        ).get(
            "RESP_input_preparation_authorized"
        )
        is True
    )

    if (
        post_qm_decision
        != EXPECTED_POST_QM_DECISION
        or not post_qm_authorized
    ):
        raise RuntimeError(
            "Coordinate-consistency audit is not "
            "authorized by the post-QM gate: "
            f"decision={post_qm_decision}; "
            f"RESP_input_preparation_authorized="
            f"{post_qm_authorized}"
        )

    execution_relative = (
        LATEST_EXECUTION_FILE
        .read_text(
            encoding="utf-8"
        )
        .strip()
    )

    if not execution_relative:
        raise RuntimeError(
            "Latest R1 execution pointer is empty: "
            f"{LATEST_EXECUTION_FILE}"
        )

    execution_dir = (
        ROOT
        / execution_relative
    )

    source_orca_xyz = (
        execution_dir
        / "v7a_r1.xyz"
    )

    source_trajectory = (
        execution_dir
        / "v7a_r1_trj.xyz"
    )

    source_orca_output = (
        execution_dir
        / "v7a_r1.out"
    )

    source_paths = {
        "post_qm_final_xyz": POST_QM_FINAL_XYZ,
        "orca_xyz": source_orca_xyz,
        "orca_trajectory": source_trajectory,
        "orca_output": source_orca_output,
    }

    for source_name, path in source_paths.items():
        if not path.is_file():
            raise RuntimeError(
                f"Missing coordinate source "
                f"{source_name}: {path}"
            )

        if path.stat().st_size == 0:
            raise RuntimeError(
                f"Empty coordinate source "
                f"{source_name}: {path}"
            )

    post_qm_atoms = read_xyz(
        POST_QM_FINAL_XYZ
    )

    orca_xyz_atoms = read_xyz(
        source_orca_xyz
    )

    trajectory_frames = read_trajectory(
        source_trajectory
    )

    if not trajectory_frames:
        raise RuntimeError(
            "No complete trajectory frames were parsed: "
            f"{source_trajectory}"
        )

    last_trajectory_frame = (
        trajectory_frames[-1]
    )

    trajectory_atoms = (
        last_trajectory_frame["atoms"]
    )

    orca_output_atoms = (
        read_last_orca_coordinates(
            source_orca_output
        )
    )

    atom_counts = {
        "post_qm_final_xyz": len(
            post_qm_atoms
        ),
        "orca_xyz": len(
            orca_xyz_atoms
        ),
        "last_trajectory_frame": len(
            trajectory_atoms
        ),
        "last_orca_output_coordinate_block": len(
            orca_output_atoms
        ),
    }

    element_orders = {
        "post_qm_final_xyz": [
            atom["element"]
            for atom in post_qm_atoms
        ],
        "orca_xyz": [
            atom["element"]
            for atom in orca_xyz_atoms
        ],
        "last_trajectory_frame": [
            atom["element"]
            for atom in trajectory_atoms
        ],
        "last_orca_output_coordinate_block": [
            atom["element"]
            for atom in orca_output_atoms
        ],
    }

    gates = {}

    gates[
        "post_QM_gate_authorizes_consistency_audit"
    ] = (
        post_qm_decision
        == EXPECTED_POST_QM_DECISION
        and post_qm_authorized
    )

    gates[
        "all_coordinate_sources_available"
    ] = all(
        path.is_file()
        and path.stat().st_size > 0
        for path in source_paths.values()
    )

    gates[
        "all_coordinate_sources_have_52_atoms"
    ] = all(
        count == EXPECTED_ATOM_COUNT
        for count in atom_counts.values()
    )

    gates[
        "all_coordinate_sources_have_identical_element_order"
    ] = all(
        order
        == element_orders["post_qm_final_xyz"]
        for order in element_orders.values()
    )

    comparison_orca_xyz = (
        compare_coordinate_sets(
            "post_qm_final_xyz",
            "orca_xyz",
            post_qm_atoms,
            orca_xyz_atoms,
        )
    )

    comparison_trajectory = (
        compare_coordinate_sets(
            "post_qm_final_xyz",
            "last_trajectory_frame",
            post_qm_atoms,
            trajectory_atoms,
        )
    )

    comparison_orca_output = (
        compare_coordinate_sets(
            "post_qm_final_xyz",
            "last_orca_output_coordinate_block",
            post_qm_atoms,
            orca_output_atoms,
        )
    )

    comparisons = {
        "post_qm_final_vs_orca_xyz": (
            comparison_orca_xyz
        ),
        "post_qm_final_vs_last_trajectory_frame": (
            comparison_trajectory
        ),
        "post_qm_final_vs_last_orca_output_block": (
            comparison_orca_output
        ),
    }

    gates[
        "post_qm_final_matches_orca_xyz"
    ] = comparison_orca_xyz[
        "coordinate_identity_pass"
    ]

    gates[
        "post_qm_final_matches_last_trajectory_frame"
    ] = comparison_trajectory[
        "coordinate_identity_pass"
    ]

    gates[
        "post_qm_final_matches_last_orca_output_block"
    ] = comparison_orca_output[
        "coordinate_identity_pass"
    ]

    all_gates_pass = all(
        gates.values()
    )

    decision = (
        PASS_DECISION
        if all_gates_pass
        else FAIL_DECISION
    )

    failure_reasons = [
        gate_name
        for gate_name, passed in gates.items()
        if not passed
    ]

    difference_rows = []

    for comparison in comparisons.values():
        difference_rows.extend(
            comparison["difference_rows"]
        )

    write_csv(
        OUTPUT_CSV,
        difference_rows,
    )

    comparison_summaries = {}

    for name, comparison in comparisons.items():
        comparison_summaries[name] = {
            key: value
            for key, value in comparison.items()
            if key != "difference_rows"
        }

    report = {
        "schema_version": 1,
        "generated_at_utc": generated_at_utc,
        "audit": (
            "QM F06 UPPER V7A R1 "
            "final-coordinate consistency"
        ),
        "execution": {
            "latest_execution_pointer": str(
                LATEST_EXECUTION_FILE
            ),
            "execution_relative_path": (
                execution_relative
            ),
            "execution_directory": str(
                execution_dir
            ),
        },
        "upstream_post_QM_gate": {
            "report": str(
                POST_QM_REPORT
            ),
            "decision": (
                post_qm_decision
            ),
            "RESP_input_preparation_authorized": (
                post_qm_authorized
            ),
            "expected_decision": (
                EXPECTED_POST_QM_DECISION
            ),
        },
        "configuration": {
            "expected_atom_count": (
                EXPECTED_ATOM_COUNT
            ),
            "coordinate_identity_tolerance_A": (
                COORDINATE_IDENTITY_TOLERANCE_A
            ),
            "comparison_reference": (
                "post_qm_final_xyz"
            ),
            "coordinate_alignment_applied": False,
            "coordinate_alignment_method": None,
        },
        "sources": {
            "post_qm_final_xyz": {
                "path": str(
                    POST_QM_FINAL_XYZ
                ),
                "sha256": sha256(
                    POST_QM_FINAL_XYZ
                ),
            },
            "orca_xyz": {
                "path": str(
                    source_orca_xyz
                ),
                "sha256": sha256(
                    source_orca_xyz
                ),
            },
            "orca_trajectory": {
                "path": str(
                    source_trajectory
                ),
                "sha256": sha256(
                    source_trajectory
                ),
                "complete_frame_count": len(
                    trajectory_frames
                ),
                "selected_frame_index_0based": (
                    last_trajectory_frame[
                        "frame_index_0based"
                    ]
                ),
                "selected_frame_comment": (
                    last_trajectory_frame[
                        "comment"
                    ]
                ),
            },
            "orca_output": {
                "path": str(
                    source_orca_output
                ),
                "sha256": sha256(
                    source_orca_output
                ),
            },
        },
        "atom_counts": atom_counts,
        "element_orders": element_orders,
        "gates": gates,
        "all_gates_pass": all_gates_pass,
        "failure_reasons": failure_reasons,
        "comparisons": comparison_summaries,
        "decision": decision,
        "authorizations": {
            "coordinate_adoption_authorized": (
                all_gates_pass
            ),
            "RESP_input_preparation_authorized": False,
            "RESP_execution_authorized": False,
            "MD_authorized": False,
        },
        "outputs": {
            "coordinate_consistency_report": str(
                OUTPUT_REPORT
            ),
            "coordinate_differences_csv": str(
                OUTPUT_CSV
            ),
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        "Coordinate consistency decision:",
        decision,
    )

    print(
        "All gates pass:",
        all_gates_pass,
    )

    print(
        "Coordinate adoption authorized:",
        report[
            "authorizations"
        ][
            "coordinate_adoption_authorized"
        ],
    )

    print(
        "Post-QM final vs ORCA XYZ max difference (A):",
        comparison_orca_xyz[
            "maximum_displacement_A"
        ],
    )

    print(
        "Post-QM final vs last trajectory frame "
        "max difference (A):",
        comparison_trajectory[
            "maximum_displacement_A"
        ],
    )

    print(
        "Post-QM final vs last ORCA output block "
        "max difference (A):",
        comparison_orca_output[
            "maximum_displacement_A"
        ],
    )

    print(
        "Report:",
        OUTPUT_REPORT,
    )

    print(
        "Differences CSV:",
        OUTPUT_CSV,
    )


if __name__ == "__main__":
    main()
