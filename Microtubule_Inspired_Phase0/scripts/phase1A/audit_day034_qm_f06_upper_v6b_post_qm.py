#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

EXEC_PARENT = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_executions"
)

LATEST_POINTER = (
    EXEC_PARENT
    / "LATEST_V6B_EXECUTION.txt"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day034_qm_f06_upper_v6b_post_qm"
)

FINAL_XYZ_OUT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_FINAL.xyz"
)

ATOM_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_atom_displacements.csv"
)

BOND_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_bond_audit.csv"
)

CONTACT_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_contacts.csv"
)

TRAJECTORY_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_trajectory_connectivity.csv"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_POST_QM_AUDIT.json"
)


EXPECTED_ATOM_COUNT = 48

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 19,
})

FIXED_ATOM_TOLERANCE_A = 5.0e-4

BN_MIN_A = 1.25
BN_MAX_A = 1.90

BH_MIN_A = 0.90
BH_MAX_A = 1.35

NH_MIN_A = 0.80
NH_MAX_A = 1.25

CLOSURE_MIN_A = 1.40
CLOSURE_MAX_A = 1.70

MINIMUM_CLOSURE_MARGIN_A = 0.02

GLOBAL_NEAR_CONTACT_MARGIN_A = 0.06

HH_HARD_CONTACT_A = 0.70
HX_HARD_CONTACT_A = 0.75
HEAVY_HEAVY_HARD_CONTACT_A = 1.10

CLOSURE_PAIRS = {
    tuple(sorted((
        "BR4:UPPER:00:3",
        "S:1739",
    ))),
    tuple(sorted((
        "BR4:UPPER:00:4",
        "BR4:UPPER:14:1",
    ))),
}


def require_file(path: Path) -> None:
    if not path.is_file():
        raise RuntimeError(
            f"Missing file: {path}"
        )

    if path.stat().st_size == 0:
        raise RuntimeError(
            f"Empty file: {path}"
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


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path) -> list[dict]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0].strip())
    coordinate_lines = lines[2:2 + count]

    if len(coordinate_lines) != count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    atoms = []

    for index, line in enumerate(
        coordinate_lines
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line: {line}"
            )

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": tuple(
                map(float, fields[1:4])
            ),
        })

    return atoms


def read_xyz_trajectory(
    path: Path,
) -> list[list[dict]]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    frames = []
    cursor = 0

    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue

        try:
            count = int(
                lines[cursor].strip()
            )
        except ValueError:
            break

        start = cursor + 2
        end = start + count

        if end > len(lines):
            break

        frame = []

        for index, line in enumerate(
            lines[start:end]
        ):
            fields = line.split()

            if len(fields) < 4:
                frame = []
                break

            frame.append({
                "index": index,
                "element": fields[0],
                "xyz_A": tuple(
                    map(float, fields[1:4])
                ),
            })

        if len(frame) != count:
            break

        frames.append(frame)
        cursor = end

    return frames


def distance(a, b) -> float:
    return math.sqrt(
        sum(
            (x - y) ** 2
            for x, y in zip(a, b)
        )
    )


def canonical_pair(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def nominal_range(
    first_element: str,
    second_element: str,
    pair: tuple[str, str],
) -> tuple[float, float, str]:
    elements = frozenset((
        first_element,
        second_element,
    ))

    if pair in CLOSURE_PAIRS:
        return (
            CLOSURE_MIN_A,
            CLOSURE_MAX_A,
            "CLOSURE_B-N",
        )

    if elements == frozenset(("B", "N")):
        return BN_MIN_A, BN_MAX_A, "B-N"

    if elements == frozenset(("B", "H")):
        return BH_MIN_A, BH_MAX_A, "B-H"

    if elements == frozenset(("N", "H")):
        return NH_MIN_A, NH_MAX_A, "N-H"

    raise RuntimeError(
        "Unsupported nominal bond class: "
        f"{first_element}-{second_element}"
    )


def geometric_bond(
    first_element: str,
    second_element: str,
    value: float,
) -> bool:
    elements = frozenset((
        first_element,
        second_element,
    ))

    if elements == frozenset(("B", "N")):
        return BN_MIN_A <= value <= BN_MAX_A

    if elements == frozenset(("B", "H")):
        return BH_MIN_A <= value <= BH_MAX_A

    if elements == frozenset(("N", "H")):
        return NH_MIN_A <= value <= NH_MAX_A

    return False


def geometric_upper_limit(
    first_element: str,
    second_element: str,
) -> float | None:
    elements = frozenset((
        first_element,
        second_element,
    ))

    if elements == frozenset(("B", "N")):
        return BN_MAX_A

    if elements == frozenset(("B", "H")):
        return BH_MAX_A

    if elements == frozenset(("N", "H")):
        return NH_MAX_A

    return None


def hard_contact_threshold(
    first_element: str,
    second_element: str,
) -> float:
    if (
        first_element == "H"
        and second_element == "H"
    ):
        return HH_HARD_CONTACT_A

    if (
        first_element == "H"
        or second_element == "H"
    ):
        return HX_HARD_CONTACT_A

    return HEAVY_HEAVY_HARD_CONTACT_A


def connected_component_count(
    atom_ids: list[str],
    edges: set[tuple[str, str]],
) -> int:
    adjacency = {
        atom_id: set()
        for atom_id in atom_ids
    }

    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    unvisited = set(atom_ids)
    components = 0

    while unvisited:
        components += 1
        seed = next(iter(unvisited))
        queue = deque([seed])
        unvisited.remove(seed)

        while queue:
            current = queue.popleft()

            for neighbor in adjacency[current]:
                if neighbor in unvisited:
                    unvisited.remove(neighbor)
                    queue.append(neighbor)

    return components


def frame_connectivity(
    frame: list[dict],
    atom_ids: list[str],
    nominal_edges: set[tuple[str, str]],
) -> dict:
    geometric_edges = set()
    distances = {}

    for first_index in range(len(frame)):
        for second_index in range(
            first_index + 1,
            len(frame),
        ):
            first_id = atom_ids[first_index]
            second_id = atom_ids[second_index]

            value = distance(
                frame[first_index]["xyz_A"],
                frame[second_index]["xyz_A"],
            )

            pair = canonical_pair(
                first_id,
                second_id,
            )

            distances[pair] = value

            if geometric_bond(
                frame[first_index]["element"],
                frame[second_index]["element"],
                value,
            ):
                geometric_edges.add(pair)

    gained = geometric_edges - nominal_edges
    lost = nominal_edges - geometric_edges

    adjacency = {
        atom_id: set()
        for atom_id in atom_ids
    }

    for first, second in geometric_edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    overcoordinated = []

    for atom_id in atom_ids:
        index = atom_ids.index(atom_id)
        element = frame[index]["element"]
        maximum_degree = 1 if element == "H" else 3

        if len(adjacency[atom_id]) > maximum_degree:
            overcoordinated.append(atom_id)

    return {
        "geometric_edges": geometric_edges,
        "distances": distances,
        "gained": gained,
        "lost": lost,
        "overcoordinated": overcoordinated,
        "components": connected_component_count(
            atom_ids,
            geometric_edges,
        ),
    }


def main() -> None:
    require_file(LATEST_POINTER)

    execution_relative = (
        LATEST_POINTER.read_text(
            encoding="utf-8"
        ).strip()
    )

    execution_dir = (
        ROOT / execution_relative
    ).resolve()

    start_xyz_path = (
        execution_dir / "v6b_start.xyz"
    )

    final_xyz_path = (
        execution_dir / "v6b.xyz"
    )

    trajectory_path = (
        execution_dir / "v6b_trj.xyz"
    )

    orca_output_path = (
        execution_dir / "v6b.out"
    )

    exit_status_path = (
        execution_dir / "v6b.exit_status"
    )

    constraint_map_path = (
        execution_dir
        / "QM_F06_UPPER_V6B_constraint_map.csv"
    )

    nominal_edges_path = (
        execution_dir
        / "QM_F06_UPPER_V6A_nominal_edges.csv"
    )

    pre_qm_report_path = (
        execution_dir
        / "QM_F06_UPPER_V6B_PRE_QM_AUDIT.json"
    )

    execution_manifest_path = (
        execution_dir
        / "QM_F06_UPPER_V6B_EXECUTION_MANIFEST.json"
    )

    required_files = [
        start_xyz_path,
        final_xyz_path,
        trajectory_path,
        orca_output_path,
        exit_status_path,
        constraint_map_path,
        nominal_edges_path,
        pre_qm_report_path,
        execution_manifest_path,
    ]

    for path in required_files:
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    map_rows = read_csv(
        constraint_map_path
    )

    retained_rows = [
        row
        for row in map_rows
        if row[
            "v6a_retained"
        ].strip().lower() == "true"
    ]

    retained_rows.sort(
        key=lambda row: int(
            row["v6b_index_0based"]
        )
    )

    atom_ids = [
        row["atom_id"]
        for row in retained_rows
    ]

    elements_from_map = [
        row["element"]
        for row in retained_rows
    ]

    fixed_by_id = {
        row["atom_id"]: (
            row["v6b_fixed"]
            .strip()
            .lower()
            == "true"
        )
        for row in retained_rows
    }

    start_atoms = read_xyz(
        start_xyz_path
    )

    final_atoms = read_xyz(
        final_xyz_path
    )

    trajectory_frames = (
        read_xyz_trajectory(
            trajectory_path
        )
    )

    identity_gate = (
        len(start_atoms)
        == len(final_atoms)
        == len(atom_ids)
        and [
            atom["element"]
            for atom in start_atoms
        ] == elements_from_map
        and [
            atom["element"]
            for atom in final_atoms
        ] == elements_from_map
    )

    atom_count_gate = (
        len(final_atoms)
        == EXPECTED_ATOM_COUNT
    )

    composition = Counter(
        atom["element"]
        for atom in final_atoms
    )

    composition_gate = (
        composition
        == EXPECTED_COMPOSITION
    )

    atom_records = []
    fixed_failures = []

    for index, atom_id in enumerate(atom_ids):
        value = distance(
            start_atoms[index]["xyz_A"],
            final_atoms[index]["xyz_A"],
        )

        fixed = fixed_by_id[atom_id]

        fixed_pass = (
            not fixed
            or value <= FIXED_ATOM_TOLERANCE_A
        )

        record = {
            "index_0based": index,
            "atom_id": atom_id,
            "element": final_atoms[index][
                "element"
            ],
            "fixed": fixed,
            "mobile": not fixed,
            "displacement_A": value,
            "fixed_tolerance_A": (
                FIXED_ATOM_TOLERANCE_A
            ),
            "fixed_displacement_pass": (
                fixed_pass
            ),
        }

        atom_records.append(record)

        if not fixed_pass:
            fixed_failures.append(record)

    with ATOM_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                atom_records[0]
            ),
        )

        writer.writeheader()
        writer.writerows(atom_records)

    nominal_edge_rows = read_csv(
        nominal_edges_path
    )

    element_by_id = {
        atom_id: elements_from_map[index]
        for index, atom_id in enumerate(
            atom_ids
        )
    }

    index_by_id = {
        atom_id: index
        for index, atom_id in enumerate(
            atom_ids
        )
    }

    nominal_edges = set()

    for row in nominal_edge_rows:
        first = row["first_atom"]
        second = row["second_atom"]

        if (
            first not in index_by_id
            or second not in index_by_id
        ):
            raise RuntimeError(
                "Nominal edge references "
                "an absent atom: "
                f"{first}--{second}"
            )

        nominal_edges.add(
            canonical_pair(first, second)
        )

    bond_records = []
    bond_failures = []
    closure_records = []
    closure_failures = []

    for first, second in sorted(
        nominal_edges
    ):
        first_index = index_by_id[first]
        second_index = index_by_id[second]

        value = distance(
            final_atoms[first_index]["xyz_A"],
            final_atoms[second_index]["xyz_A"],
        )

        minimum, maximum, bond_class = (
            nominal_range(
                element_by_id[first],
                element_by_id[second],
                (first, second),
            )
        )

        margin = min(
            value - minimum,
            maximum - value,
        )

        passed = (
            minimum <= value <= maximum
        )

        if (
            (first, second)
            in CLOSURE_PAIRS
        ):
            passed = (
                passed
                and margin
                >= MINIMUM_CLOSURE_MARGIN_A
            )

        record = {
            "first_atom": first,
            "first_element": (
                element_by_id[first]
            ),
            "second_atom": second,
            "second_element": (
                element_by_id[second]
            ),
            "bond_class": bond_class,
            "distance_A": value,
            "minimum_A": minimum,
            "maximum_A": maximum,
            "margin_A": margin,
            "closure_edge": (
                (first, second)
                in CLOSURE_PAIRS
            ),
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

        if (
            (first, second)
            in CLOSURE_PAIRS
        ):
            closure_records.append(record)

            if not passed:
                closure_failures.append(
                    record
                )

    with BOND_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                bond_records[0]
            ),
        )

        writer.writeheader()
        writer.writerows(bond_records)

    final_connectivity = frame_connectivity(
        final_atoms,
        atom_ids,
        nominal_edges,
    )

    geometric_edges = (
        final_connectivity[
            "geometric_edges"
        ]
    )

    gained_edges = sorted(
        final_connectivity["gained"]
    )

    lost_edges = sorted(
        final_connectivity["lost"]
    )

    overcoordinated = (
        final_connectivity[
            "overcoordinated"
        ]
    )

    component_count = (
        final_connectivity[
            "components"
        ]
    )

    contact_records = []
    hard_contacts = []
    global_near_contact_conflicts = []

    for first_index in range(
        len(final_atoms)
    ):
        for second_index in range(
            first_index + 1,
            len(final_atoms),
        ):
            first = atom_ids[first_index]
            second = atom_ids[second_index]
            pair = canonical_pair(
                first,
                second,
            )

            if pair in nominal_edges:
                continue

            first_element = (
                final_atoms[first_index][
                    "element"
                ]
            )

            second_element = (
                final_atoms[second_index][
                    "element"
                ]
            )

            value = distance(
                final_atoms[first_index][
                    "xyz_A"
                ],
                final_atoms[second_index][
                    "xyz_A"
                ],
            )

            hard_threshold = (
                hard_contact_threshold(
                    first_element,
                    second_element,
                )
            )

            hard_contact = (
                value < hard_threshold
            )

            upper_limit = (
                geometric_upper_limit(
                    first_element,
                    second_element,
                )
            )

            if upper_limit is None:
                clearance = None
                near_conflict = False
            else:
                clearance = (
                    value - upper_limit
                )

                near_conflict = (
                    0.0
                    < clearance
                    < GLOBAL_NEAR_CONTACT_MARGIN_A
                )

            record = {
                "first_atom": first,
                "first_element": first_element,
                "second_atom": second,
                "second_element": second_element,
                "distance_A": value,
                "hard_contact_threshold_A": (
                    hard_threshold
                ),
                "hard_contact": hard_contact,
                "bond_upper_limit_A": (
                    ""
                    if upper_limit is None
                    else upper_limit
                ),
                "nonnominal_clearance_A": (
                    ""
                    if clearance is None
                    else clearance
                ),
                "global_near_contact_conflict": (
                    near_conflict
                ),
            }

            contact_records.append(record)

            if hard_contact:
                hard_contacts.append(record)

            if near_conflict:
                global_near_contact_conflicts.append(
                    record
                )

    with CONTACT_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                contact_records[0]
            ),
        )

        writer.writeheader()
        writer.writerows(contact_records)

    trajectory_records = []
    trajectory_violation_frames = []

    for frame_index, frame in enumerate(
        trajectory_frames
    ):
        if len(frame) != len(atom_ids):
            trajectory_records.append({
                "frame_index": frame_index,
                "atom_count": len(frame),
                "gained_edges": "",
                "lost_edges": "",
                "overcoordinated_atoms": "",
                "connected_components": "",
                "topology_matches_nominal": False,
            })

            trajectory_violation_frames.append(
                frame_index
            )
            continue

        connectivity = frame_connectivity(
            frame,
            atom_ids,
            nominal_edges,
        )

        topology_matches = (
            len(connectivity["gained"]) == 0
            and len(connectivity["lost"]) == 0
            and len(
                connectivity[
                    "overcoordinated"
                ]
            ) == 0
            and connectivity[
                "components"
            ] == 1
        )

        record = {
            "frame_index": frame_index,
            "atom_count": len(frame),
            "gained_edges": "|".join(
                f"{first}--{second}"
                for first, second in sorted(
                    connectivity["gained"]
                )
            ),
            "lost_edges": "|".join(
                f"{first}--{second}"
                for first, second in sorted(
                    connectivity["lost"]
                )
            ),
            "overcoordinated_atoms": (
                "|".join(
                    sorted(
                        connectivity[
                            "overcoordinated"
                        ]
                    )
                )
            ),
            "connected_components": (
                connectivity[
                    "components"
                ]
            ),
            "topology_matches_nominal": (
                topology_matches
            ),
        }

        trajectory_records.append(record)

        if not topology_matches:
            trajectory_violation_frames.append(
                frame_index
            )

    with TRAJECTORY_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                trajectory_records[0]
            ),
        )

        writer.writeheader()
        writer.writerows(
            trajectory_records
        )

    orca_text = (
        orca_output_path.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    exit_status = int(
        exit_status_path.read_text(
            encoding="utf-8"
        ).strip()
    )

    normal_termination_gate = (
        exit_status == 0
        and "THE OPTIMIZATION HAS CONVERGED"
        in orca_text
        and "ORCA TERMINATED NORMALLY"
        in orca_text
        and "ORCA finished by error termination"
        not in orca_text
        and "SCF NOT CONVERGED"
        not in orca_text
    )

    final_energy_lines = [
        line.strip()
        for line in orca_text.splitlines()
        if "FINAL SINGLE POINT ENERGY"
        in line
    ]

    final_energy_Eh = None

    if final_energy_lines:
        final_energy_Eh = float(
            final_energy_lines[-1].split()[-1]
        )

    fixed_geometry_gate = (
        len(fixed_failures) == 0
    )

    bond_gate = (
        len(bond_failures) == 0
    )

    closure_gate = (
        len(closure_records)
        == len(CLOSURE_PAIRS)
        and len(closure_failures) == 0
    )

    reconnectivity_gate = (
        len(gained_edges) == 0
        and len(lost_edges) == 0
    )

    overcoordination_gate = (
        len(overcoordinated) == 0
    )

    connected_gate = (
        component_count == 1
    )

    hard_contact_gate = (
        len(hard_contacts) == 0
    )

    global_near_contact_gate = (
        len(
            global_near_contact_conflicts
        ) == 0
    )

    trajectory_available_gate = (
        len(trajectory_frames) > 0
    )

    gates = {
        "normal_ORCA_termination_and_optimization_convergence": (
            normal_termination_gate
        ),
        "atom_identity_and_order": (
            identity_gate
        ),
        "atom_count": atom_count_gate,
        "composition": composition_gate,
        "fixed_atoms_preserved": (
            fixed_geometry_gate
        ),
        "all_nominal_bonds_in_range": (
            bond_gate
        ),
        "closure_BN_geometry": (
            closure_gate
        ),
        "no_final_geometric_reconnectivity": (
            reconnectivity_gate
        ),
        "no_final_overcoordinated_atoms": (
            overcoordination_gate
        ),
        "single_final_connected_component": (
            connected_gate
        ),
        "no_final_hard_contacts": (
            hard_contact_gate
        ),
        "global_nonnominal_near_contact_gate": (
            global_near_contact_gate
        ),
        "trajectory_available_and_parseable": (
            trajectory_available_gate
        ),
    }

    passed = all(gates.values())

    if passed:
        decision = (
            "QM_F06_UPPER_V6B_POST_QM_GATE_PASS_"
            "RESP_INPUT_PREPARATION_AUTHORIZED"
        )
    else:
        decision = (
            "QM_F06_UPPER_V6B_POST_QM_GATE_FAIL_"
            "STRUCTURAL_REVIEW_REQUIRED"
        )

    with FINAL_XYZ_OUT.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"{len(final_atoms)}\n"
        )

        handle.write(
            "QM_F06 UPPER V6-B final "
            "post-QM audited geometry; "
            f"gate_pass={passed}; "
            f"energy_Eh={final_energy_Eh}\n"
        )

        for atom in final_atoms:
            x_value, y_value, z_value = (
                atom["xyz_A"]
            )

            handle.write(
                f"{atom['element']:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    files_for_hash = {
        "start_xyz": start_xyz_path,
        "final_orca_xyz": final_xyz_path,
        "audited_final_xyz": FINAL_XYZ_OUT,
        "trajectory": trajectory_path,
        "orca_output": orca_output_path,
        "exit_status": exit_status_path,
        "constraint_map": constraint_map_path,
        "nominal_edges": nominal_edges_path,
        "pre_qm_report": pre_qm_report_path,
        "execution_manifest": execution_manifest_path,
        "atom_audit": ATOM_AUDIT_CSV,
        "bond_audit": BOND_AUDIT_CSV,
        "contact_audit": CONTACT_AUDIT_CSV,
        "trajectory_audit": TRAJECTORY_AUDIT_CSV,
    }

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "ORCA": {
            "exit_status": exit_status,
            "optimization_converged": (
                "THE OPTIMIZATION HAS CONVERGED"
                in orca_text
            ),
            "terminated_normally": (
                "ORCA TERMINATED NORMALLY"
                in orca_text
            ),
            "final_energy_Eh": (
                final_energy_Eh
            ),
        },
        "geometry": {
            "atom_count": len(final_atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "fixed_atom_count": sum(
                fixed_by_id.values()
            ),
            "mobile_atom_count": sum(
                not value
                for value in fixed_by_id.values()
            ),
            "maximum_fixed_displacement_A": (
                max(
                    record["displacement_A"]
                    for record in atom_records
                    if record["fixed"]
                )
            ),
            "fixed_displacement_tolerance_A": (
                FIXED_ATOM_TOLERANCE_A
            ),
        },
        "topology": {
            "nominal_edge_count": len(
                nominal_edges
            ),
            "geometric_edge_count": len(
                geometric_edges
            ),
            "gained_edge_count": len(
                gained_edges
            ),
            "gained_edges": [
                f"{first}--{second}"
                for first, second in gained_edges
            ],
            "lost_edge_count": len(
                lost_edges
            ),
            "lost_edges": [
                f"{first}--{second}"
                for first, second in lost_edges
            ],
            "overcoordinated_atom_count": (
                len(overcoordinated)
            ),
            "overcoordinated_atoms": (
                sorted(overcoordinated)
            ),
            "connected_component_count": (
                component_count
            ),
        },
        "closure_bonds": closure_records,
        "contact_screen": {
            "hard_contact_count": len(
                hard_contacts
            ),
            "global_nonnominal_near_contact_margin_A": (
                GLOBAL_NEAR_CONTACT_MARGIN_A
            ),
            "global_nonnominal_near_contact_conflict_count": (
                len(
                    global_near_contact_conflicts
                )
            ),
            "global_nonnominal_near_contact_conflicts": (
                global_near_contact_conflicts
            ),
        },
        "trajectory_diagnostic": {
            "complete_frame_count": len(
                trajectory_frames
            ),
            "topology_violation_frame_count": (
                len(
                    trajectory_violation_frames
                )
            ),
            "topology_violation_frames": (
                trajectory_violation_frames
            ),
            "note": (
                "Trajectory topology is diagnostic. "
                "RESP authorization is governed by "
                "the converged final geometry."
            ),
        },
        "failure_counts": {
            "fixed_atom_failures": len(
                fixed_failures
            ),
            "bond_failures": len(
                bond_failures
            ),
            "closure_failures": len(
                closure_failures
            ),
            "gained_edges": len(
                gained_edges
            ),
            "lost_edges": len(
                lost_edges
            ),
            "overcoordinated_atoms": len(
                overcoordinated
            ),
            "hard_contacts": len(
                hard_contacts
            ),
            "global_near_contact_conflicts": (
                len(
                    global_near_contact_conflicts
                )
            ),
        },
        "gates": gates,
        "structural_pass": passed,
        "files": {
            key: str(
                path.relative_to(ROOT)
            )
            for key, path
            in files_for_hash.items()
        },
        "sha256": {
            key: sha256(path)
            for key, path
            in files_for_hash.items()
        },
        "authorization": {
            "RESP_input_preparation_authorized": (
                passed
            ),
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": (
                False
            ),
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

    print("=" * 108)
    print("QM_F06 UPPER V6-B FORMAL POST-QM AUDIT")
    print("=" * 108)

    for name, value in gates.items():
        print(
            f"{name:58s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Final energy Eh:", final_energy_Eh)
    print("Atoms:", len(final_atoms))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print(
        "Maximum fixed displacement A:",
        report["geometry"][
            "maximum_fixed_displacement_A"
        ],
    )

    print()
    print("Nominal edges:", len(nominal_edges))
    print("Geometric edges:", len(geometric_edges))
    print("Bond failures:", len(bond_failures))
    print("Closure failures:", len(closure_failures))
    print("Gained edges:", len(gained_edges))
    print("Lost edges:", len(lost_edges))
    print(
        "Overcoordinated atoms:",
        len(overcoordinated),
    )
    print(
        "Connected components:",
        component_count,
    )
    print("Hard contacts:", len(hard_contacts))
    print(
        "Global near-contact conflicts:",
        len(global_near_contact_conflicts),
    )

    print()
    print("Closure B-N bonds:")

    for record in closure_records:
        print(
            f"  {record['first_atom']:28s} -- "
            f"{record['second_atom']:28s} "
            f"{float(record['distance_A']):.6f} Å | "
            f"margin={float(record['margin_A']):.6f} Å | "
            f"{'PASS' if record['pass'] else 'FAIL'}"
        )

    print()
    print(
        "Trajectory complete frames:",
        len(trajectory_frames),
    )
    print(
        "Trajectory diagnostic violation frames:",
        len(trajectory_violation_frames),
    )

    if trajectory_violation_frames:
        print(
            "First trajectory violation frames:",
            trajectory_violation_frames[:20],
        )

    if fixed_failures:
        print()
        print("Fixed-atom failures:")

        for record in fixed_failures:
            print(
                f"  {record['atom_id']:28s} "
                f"{float(record['displacement_A']):.8f} Å"
            )

    if bond_failures:
        print()
        print("Bond failures:")

        for record in bond_failures:
            print(
                f"  {record['first_atom']:28s} -- "
                f"{record['second_atom']:28s} "
                f"{float(record['distance_A']):.6f} Å"
            )

    if gained_edges:
        print()
        print("Nonnominal geometric edges:")

        for first, second in gained_edges:
            print(f"  {first} -- {second}")

    if lost_edges:
        print()
        print("Lost nominal edges:")

        for first, second in lost_edges:
            print(f"  {first} -- {second}")

    if global_near_contact_conflicts:
        print()
        print("Global near-contact conflicts:")

        for record in (
            global_near_contact_conflicts
        ):
            print(
                f"  {record['first_atom']:28s} -- "
                f"{record['second_atom']:28s} "
                f"{float(record['distance_A']):.6f} Å | "
                f"clearance="
                f"{float(record['nonnominal_clearance_A']):.6f} Å"
            )

    print()
    print("Decision:", decision)
    print("Final XYZ:", FINAL_XYZ_OUT)
    print("Atom audit:", ATOM_AUDIT_CSV)
    print("Bond audit:", BOND_AUDIT_CSV)
    print("Contact audit:", CONTACT_AUDIT_CSV)
    print(
        "Trajectory audit:",
        TRAJECTORY_AUDIT_CSV,
    )
    print("Report:", REPORT_PATH)

    print()
    print(
        "RESP input preparation authorized:",
        passed,
    )
    print("RESP execution authorized: False")
    print(
        "Force-field adoption authorized: False"
    )
    print("MD authorized: False")


if __name__ == "__main__":
    main()
