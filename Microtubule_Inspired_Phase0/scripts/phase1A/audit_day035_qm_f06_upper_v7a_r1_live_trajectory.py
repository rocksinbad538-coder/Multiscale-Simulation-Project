#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import math


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
    / "day035_qm_f06_upper_v7a_r1_live_trajectory"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_live_trajectory_audit.csv"
)

OUTPUT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_LIVE_TRAJECTORY_AUDIT.json"
)

EXPECTED_ATOM_COUNT = 52

EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_DEGREE = {
    "B": 3,
    "N": 3,
    "H": 1,
}

BOND_WINDOWS = {
    ("B", "N"): (1.35, 1.75),
    ("B", "H"): (0.90, 1.35),
    ("H", "N"): (0.80, 1.25),
}

GEOMETRIC_BOND_MAXIMUM = {
    ("B", "N"): 1.90,
    ("B", "H"): 1.35,
    ("H", "N"): 1.25,
}

HARD_CONTACT_MINIMUM = {
    ("B", "B"): 1.20,
    ("B", "N"): 1.20,
    ("N", "N"): 1.20,
    ("B", "H"): 0.75,
    ("H", "N"): 0.70,
    ("H", "H"): 0.65,
}

REQUIRED_V7A_EDGES = {
    frozenset((
        "A:UPPER:13:1",
        "A:UPPER:14:2",
    )),
    frozenset((
        "A:UPPER:13:1",
        "A:UPPER:14:0",
    )),
    frozenset((
        "A:UPPER:13:1",
        "HCAPV7:UPPER:A13_1:A11_1",
    )),
    frozenset((
        "A:UPPER:14:0",
        "HCAPV7:UPPER:A14_0:A13_M1",
    )),
    frozenset((
        "A:UPPER:14:0",
        "H4:UPPER:0045:0",
    )),
}

FORBIDDEN_V6B_PAIR = frozenset((
    "A:UPPER:14:2",
    "P:1641",
))


def canonical_elements(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def read_csv(path: Path) -> list[dict]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_trajectory(
    path: Path,
) -> list[dict]:
    if not path.is_file() or path.stat().st_size == 0:
        return []

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

        frame_end = cursor + atom_count + 2

        if frame_end > len(lines):
            break

        comment = lines[cursor + 1]
        atoms = []
        valid = True

        for index, line in enumerate(
            lines[cursor + 2:frame_end]
        ):
            fields = line.split()

            if len(fields) < 4:
                valid = False
                break

            try:
                xyz = tuple(
                    float(value)
                    for value in fields[1:4]
                )
            except ValueError:
                valid = False
                break

            atoms.append({
                "index_0based": index,
                "element": fields[0],
                "xyz_A": xyz,
            })

        if not valid:
            break

        frames.append({
            "frame_index_0based": len(frames),
            "comment": comment,
            "atom_count": atom_count,
            "atoms": atoms,
        })

        cursor = frame_end

    return frames


def extract_energy(comment: str) -> float | None:
    fields = comment.split()

    for index, token in enumerate(fields):
        if token == "E" and index + 1 < len(fields):
            try:
                return float(fields[index + 1])
            except ValueError:
                return None

    return None


def main() -> None:
    execution_relative = (
        LATEST_FILE
        .read_text(encoding="utf-8")
        .strip()
    )

    execution_dir = ROOT / execution_relative

    trajectory_path = (
        execution_dir
        / "v7a_r1_trj.xyz"
    )

    map_path = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_constraint_map.csv"
    )

    edges_path = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_nominal_edges.csv"
    )

    for path in (
        trajectory_path,
        map_path,
        edges_path,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing required input: {path}"
            )

    map_rows = read_csv(map_path)

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    atom_ids = [
        row["atom_id"]
        for row in map_rows
    ]

    elements = {
        row["atom_id"]: row["element"]
        for row in map_rows
    }

    edge_rows = read_csv(edges_path)

    nominal_edges = {
        frozenset((
            row["first_atom"],
            row["second_atom"],
        ))
        for row in edge_rows
    }

    frames = read_trajectory(
        trajectory_path
    )

    if not frames:
        raise RuntimeError(
            "No complete trajectory frames"
        )

    frame_records = []
    violation_frames = []

    for frame in frames:
        atoms = frame["atoms"]

        identity_pass = (
            len(atoms) == len(map_rows)
            and all(
                atoms[index]["element"]
                == map_rows[index]["element"]
                for index in range(len(map_rows))
            )
        )

        coordinates = {
            atom_ids[index]:
            atoms[index]["xyz_A"]
            for index in range(
                min(len(atoms), len(atom_ids))
            )
        }

        composition = Counter(
            atom["element"]
            for atom in atoms
        )

        geometric_edges = set()
        hard_contacts = []
        minimum_nonnominal_margin = None

        for first_index in range(len(atom_ids)):
            first = atom_ids[first_index]
            first_element = elements[first]

            for second_index in range(
                first_index + 1,
                len(atom_ids),
            ):
                second = atom_ids[second_index]
                second_element = elements[second]

                pair = frozenset((first, second))
                pair_class = canonical_elements(
                    first_element,
                    second_element,
                )

                value = distance(
                    coordinates[first],
                    coordinates[second],
                )

                hard_minimum = (
                    HARD_CONTACT_MINIMUM.get(
                        pair_class
                    )
                )

                if (
                    hard_minimum is not None
                    and value < hard_minimum
                ):
                    hard_contacts.append({
                        "first_atom": first,
                        "second_atom": second,
                        "distance_A": value,
                    })

                geometric_maximum = (
                    GEOMETRIC_BOND_MAXIMUM.get(
                        pair_class
                    )
                )

                if (
                    geometric_maximum is not None
                    and value <= geometric_maximum
                ):
                    geometric_edges.add(pair)

                if (
                    pair not in nominal_edges
                    and geometric_maximum is not None
                ):
                    margin = (
                        value - geometric_maximum
                    )

                    if (
                        minimum_nonnominal_margin
                        is None
                        or margin
                        < minimum_nonnominal_margin
                    ):
                        minimum_nonnominal_margin = margin

        lost_edges = (
            nominal_edges - geometric_edges
        )

        gained_edges = (
            geometric_edges - nominal_edges
        )

        adjacency = defaultdict(set)

        for pair in geometric_edges:
            first, second = tuple(pair)

            adjacency[first].add(second)
            adjacency[second].add(first)

        degree_failures = []

        for atom_id in atom_ids:
            expected = EXPECTED_DEGREE[
                elements[atom_id]
            ]

            observed = len(
                adjacency[atom_id]
            )

            if observed != expected:
                degree_failures.append({
                    "atom_id": atom_id,
                    "element": elements[atom_id],
                    "degree": observed,
                    "expected_degree": expected,
                })

        components = []
        visited = set()

        for atom_id in atom_ids:
            if atom_id in visited:
                continue

            stack = [atom_id]
            component = set()

            while stack:
                current = stack.pop()

                if current in component:
                    continue

                component.add(current)
                visited.add(current)

                stack.extend(
                    adjacency[current]
                    - component
                )

            components.append(component)

        nominal_bond_failures = []

        for row in edge_rows:
            first = row["first_atom"]
            second = row["second_atom"]

            pair_class = canonical_elements(
                row["first_element"],
                row["second_element"],
            )

            minimum, maximum = (
                BOND_WINDOWS[pair_class]
            )

            value = distance(
                coordinates[first],
                coordinates[second],
            )

            if not (
                minimum <= value <= maximum
            ):
                nominal_bond_failures.append({
                    "first_atom": first,
                    "second_atom": second,
                    "distance_A": value,
                    "minimum_A": minimum,
                    "maximum_A": maximum,
                })

        required_edges_present = (
            REQUIRED_V7A_EDGES
            <= geometric_edges
        )

        forbidden_pair_absent = (
            FORBIDDEN_V6B_PAIR
            not in geometric_edges
        )

        frame_pass = all((
            identity_pass,
            len(atoms) == EXPECTED_ATOM_COUNT,
            dict(composition)
            == EXPECTED_COMPOSITION,
            len(nominal_bond_failures) == 0,
            len(lost_edges) == 0,
            len(gained_edges) == 0,
            len(degree_failures) == 0,
            len(components) == 1,
            len(hard_contacts) == 0,
            required_edges_present,
            forbidden_pair_absent,
        ))

        record = {
            "frame_index_0based": (
                frame["frame_index_0based"]
            ),
            "energy_Eh": (
                extract_energy(frame["comment"])
            ),
            "atom_count": len(atoms),
            "identity_pass": identity_pass,
            "composition_pass": (
                dict(composition)
                == EXPECTED_COMPOSITION
            ),
            "geometric_edge_count": (
                len(geometric_edges)
            ),
            "nominal_bond_failure_count": (
                len(nominal_bond_failures)
            ),
            "lost_edge_count": len(lost_edges),
            "gained_edge_count": len(gained_edges),
            "degree_failure_count": (
                len(degree_failures)
            ),
            "connected_component_count": (
                len(components)
            ),
            "hard_contact_count": (
                len(hard_contacts)
            ),
            "minimum_nonnominal_margin_A": (
                minimum_nonnominal_margin
            ),
            "required_V7A_edges_present": (
                required_edges_present
            ),
            "V6B_forbidden_pair_absent": (
                forbidden_pair_absent
            ),
            "frame_pass": frame_pass,
        }

        frame_records.append(record)

        if not frame_pass:
            violation_frames.append({
                "frame_index_0based": (
                    frame["frame_index_0based"]
                ),
                "nominal_bond_failures": (
                    nominal_bond_failures
                ),
                "lost_edges": [
                    sorted(pair)
                    for pair in lost_edges
                ],
                "gained_edges": [
                    sorted(pair)
                    for pair in gained_edges
                ],
                "degree_failures": (
                    degree_failures
                ),
                "hard_contacts": hard_contacts,
            })

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                frame_records[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(frame_records)

    decision = (
        "QM_F06_UPPER_V7A_R1_"
        "LIVE_TRAJECTORY_DIAGNOSTIC_PASS"
        if not violation_frames
        else
        "QM_F06_UPPER_V7A_R1_"
        "LIVE_TRAJECTORY_DIAGNOSTIC_"
        "VIOLATIONS_DETECTED"
    )

    report = {
        "model": "QM_F06_UPPER_V7A_R1",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir
        ),
        "trajectory": str(
            trajectory_path
        ),
        "decision": decision,
        "summary": {
            "complete_frame_count": (
                len(frames)
            ),
            "passing_frame_count": (
                len(frames)
                - len(violation_frames)
            ),
            "violation_frame_count": (
                len(violation_frames)
            ),
            "first_violation_frame": (
                violation_frames[0][
                    "frame_index_0based"
                ]
                if violation_frames
                else None
            ),
        },
        "violation_frames": violation_frames,
        "authorizations": {
            "continue_current_ORCA_execution": (
                len(violation_frames) == 0
            ),
            "post_QM_acceptance_authorized": False,
            "RESP_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_JSON.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "LIVE TRAJECTORY STRUCTURAL DIAGNOSTIC"
    )
    print("=" * 112)

    print(
        "Complete frames:",
        len(frames),
    )
    print(
        "Passing frames:",
        len(frames)
        - len(violation_frames),
    )
    print(
        "Violation frames:",
        len(violation_frames),
    )

    for row in frame_records:
        print(
            f"frame={row['frame_index_0based']:3d} | "
            f"E={row['energy_Eh']} | "
            f"edges={row['geometric_edge_count']:2d} | "
            f"lost={row['lost_edge_count']:2d} | "
            f"gained={row['gained_edge_count']:2d} | "
            f"degree_fail={row['degree_failure_count']:2d} | "
            f"components={row['connected_component_count']:2d} | "
            f"hard={row['hard_contact_count']:2d} | "
            f"{'PASS' if row['frame_pass'] else 'FAIL'}"
        )

    print()
    print("Decision:", decision)
    print("CSV:", OUTPUT_CSV)
    print("Report:", OUTPUT_JSON)
    print()
    print(
        "Continue current ORCA execution:",
        len(violation_frames) == 0,
    )
    print(
        "Post-QM acceptance authorized: False"
    )
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
