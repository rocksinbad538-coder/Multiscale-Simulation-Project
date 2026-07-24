#!/usr/bin/env python3
"""
Post-QM geometry auditor for QM_F06 UPPER V5-B.

The auditor:
- locates the execution selected by LATEST_V5B_EXECUTION.txt;
- requires normal ORCA termination and geometry convergence for a final audit;
- supports --allow-incomplete for non-authorizing live diagnostics;
- extracts the latest complete geometry;
- checks composition and atom identity;
- verifies displacement of fixed atoms;
- imports the validated pre-QM bond graph;
- checks connectivity, nominal valence, bond preservation and cap ownership;
- screens non-topological hard contacts;
- verifies the complete V5-B new-cap inventory, ownership and bond ranges;
- writes deterministic JSON, CSV, XYZ and SHA256 artifacts.

RESP, force-field adoption and MD remain blocked unless every final gate passes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_PARENT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_V5B_EXECUTION.txt"
)

PRE_QM_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit"
)

PRE_QM_REPORT = (
    PRE_QM_DIR
    / "QM_F06_UPPER_V5B_PRE_QM_AUDIT.json"
)

PRE_QM_VALENCE = (
    PRE_QM_DIR
    / "QM_F06_UPPER_V5B_valence.csv"
)

PRE_QM_CAPS = (
    PRE_QM_DIR
    / "QM_F06_UPPER_V5B_new_cap_audit.csv"
)

CONSTRAINT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_orca_input"
)

CONSTRAINT_MAP = (
    CONSTRAINT_DIR
    / "QM_F06_UPPER_V5B_constraint_map.csv"
)

CONSTRUCTION_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction"
)

PROVENANCE_MAP = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

CONSTRUCTION_CAPS = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_new_artificial_caps.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5b_post_qm"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_POST_QM_AUDIT.json"
)

FINAL_XYZ_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_FINAL.xyz"
)

ATOM_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_atom_displacements.csv"
)

BOND_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_bond_audit.csv"
)

CAP_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_new_cap_audit.csv"
)

CONTACT_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_contacts.csv"
)

EXPECTED_ATOM_COUNT = 52

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 23,
})

FIXED_DISPLACEMENT_TOLERANCE_A = 5.0e-4

BN_MIN_A = 1.25
BN_MAX_A = 1.85

BH_MIN_A = 0.95
BH_MAX_A = 1.35

NH_MIN_A = 0.85
NH_MAX_A = 1.25

HH_HARD_CONTACT_A = 1.20
HX_HARD_CONTACT_A = 0.85
HEAVY_HEAVY_HARD_CONTACT_A = 1.20

CAP_OWNER_MARGIN_A = 0.10

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


def distance(first, second) -> float:
    return math.sqrt(sum(
        (a - b) ** 2
        for a, b in zip(first, second)
    ))


def angle(first, center, second) -> float:
    vector_1 = tuple(
        a - b
        for a, b in zip(first, center)
    )

    vector_2 = tuple(
        a - b
        for a, b in zip(second, center)
    )

    norm_1 = math.sqrt(sum(
        value * value
        for value in vector_1
    ))

    norm_2 = math.sqrt(sum(
        value * value
        for value in vector_2
    ))

    if norm_1 == 0.0 or norm_2 == 0.0:
        raise RuntimeError(
            "Cannot calculate angle from a zero-length vector."
        )

    cosine = sum(
        a * b
        for a, b in zip(vector_1, vector_2)
    ) / (norm_1 * norm_2)

    cosine = max(-1.0, min(1.0, cosine))

    return math.degrees(math.acos(cosine))


def canonical_pair(first: str, second: str):
    return tuple(sorted((first, second)))


def read_xyz(path: Path) -> list[dict]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    declared = int(lines[0].strip())

    atoms = []

    for index, line in enumerate(
        lines[2:2 + declared]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index}: {line}"
            )

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    if len(atoms) != declared:
        raise RuntimeError(
            f"Incomplete XYZ: expected {declared}, "
            f"found {len(atoms)}"
        )

    return atoms


def read_xyz_trajectory(path: Path) -> list[list[dict]]:
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
            declared = int(lines[cursor].strip())
        except ValueError:
            break

        frame_end = cursor + declared + 2

        if frame_end > len(lines):
            break

        atoms = []

        for index, line in enumerate(
            lines[cursor + 2:frame_end]
        ):
            fields = line.split()

            if len(fields) < 4:
                atoms = []
                break

            atoms.append({
                "index": index,
                "element": fields[0],
                "xyz_A": tuple(
                    float(value)
                    for value in fields[1:4]
                ),
            })

        if len(atoms) != declared:
            break

        frames.append(atoms)
        cursor = frame_end

    return frames


def select_latest_geometry(
    execution_dir: Path,
) -> tuple[Path, list[dict], int | None]:
    final_xyz = execution_dir / "v5b.xyz"
    trajectory = execution_dir / "v5b_trj.xyz"
    start_xyz = execution_dir / "v5b_start.xyz"

    if final_xyz.is_file() and final_xyz.stat().st_size > 0:
        atoms = read_xyz(final_xyz)

        if len(atoms) == EXPECTED_ATOM_COUNT:
            return final_xyz, atoms, None

    if trajectory.is_file() and trajectory.stat().st_size > 0:
        frames = read_xyz_trajectory(trajectory)

        if frames:
            return trajectory, frames[-1], len(frames) - 1

    return start_xyz, read_xyz(start_xyz), 0


def parse_orca_status(text: str) -> dict:
    energies = [
        float(value)
        for value in re.findall(
            r"FINAL SINGLE POINT ENERGY\s+"
            r"([-+]?\d+\.\d+)",
            text,
        )
    ]

    geometry_cycles = len(re.findall(
        r"GEOMETRY OPTIMIZATION CYCLE",
        text,
        re.IGNORECASE,
    ))

    scf_convergences = len(re.findall(
        r"SCF CONVERGED",
        text,
        re.IGNORECASE,
    ))

    normal_termination = (
        "ORCA TERMINATED NORMALLY" in text
    )

    geometry_converged = (
        "THE OPTIMIZATION HAS CONVERGED" in text
    )

    error_patterns = {
        "error_termination": (
            "ORCA finished by error termination" in text
        ),
        "aborting": bool(re.search(
            r"\baborting\b",
            text,
            re.IGNORECASE,
        )),
        "segmentation_fault": bool(re.search(
            r"segmentation fault",
            text,
            re.IGNORECASE,
        )),
        "scf_not_converged": bool(re.search(
            r"SCF NOT CONVERGED|SCF FAILED",
            text,
            re.IGNORECASE,
        )),
    }

    return {
        "normal_termination": normal_termination,
        "geometry_converged": geometry_converged,
        "geometry_cycles_started": geometry_cycles,
        "scf_convergences": scf_convergences,
        "completed_energy_points": len(energies),
        "last_completed_energy_hartree": (
            energies[-1] if energies else None
        ),
        "error_patterns": error_patterns,
        "error_detected": any(
            error_patterns.values()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help=(
            "Audit the latest complete frame without "
            "authorizing final acceptance."
        ),
    )

    args = parser.parse_args()

    for path in (
        LATEST_POINTER,
        PRE_QM_REPORT,
        PRE_QM_VALENCE,
        PRE_QM_CAPS,
        CONSTRAINT_MAP,
        PROVENANCE_MAP,
    ):
        require_file(path)

    execution_relative = Path(
        LATEST_POINTER.read_text(
            encoding="utf-8"
        ).strip()
    )

    execution_dir = ROOT / execution_relative

    output_path = execution_dir / "v5b.out"
    stderr_path = execution_dir / "v5b.stderr"
    start_xyz_path = execution_dir / "v5b_start.xyz"
    manifest_path = (
        execution_dir
        / "QM_F06_UPPER_V5B_EXECUTION_MANIFEST.json"
    )

    for path in (
        output_path,
        start_xyz_path,
        manifest_path,
    ):
        require_file(path)

    output_text = output_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    status = parse_orca_status(output_text)

    execution_complete = (
        status["normal_termination"]
        and status["geometry_converged"]
        and not status["error_detected"]
    )

    if not execution_complete and not args.allow_incomplete:
        raise RuntimeError(
            "Final audit blocked: ORCA has not completed "
            "with normal termination and converged geometry. "
            "Use --allow-incomplete for a diagnostic-only audit."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    start_atoms = read_xyz(start_xyz_path)

    (
        geometry_source,
        final_atoms,
        trajectory_frame_index,
    ) = select_latest_geometry(execution_dir)

    constraint_rows = read_csv(CONSTRAINT_MAP)
    provenance_rows = read_csv(PROVENANCE_MAP)
    valence_rows = read_csv(PRE_QM_VALENCE)
    pre_qm_cap_rows = read_csv(PRE_QM_CAPS)
    construction_cap_rows = read_csv(CONSTRUCTION_CAPS)

    if not (
        len(start_atoms)
        == len(final_atoms)
        == len(constraint_rows)
        == len(provenance_rows)
        == EXPECTED_ATOM_COUNT
    ):
        raise RuntimeError(
            "Atom-count mismatch among V4 artifacts."
        )

    mapped = {}

    for index, (
        start,
        final,
        constraint,
        provenance,
    ) in enumerate(zip(
        start_atoms,
        final_atoms,
        constraint_rows,
        provenance_rows,
        strict=True,
    )):
        if int(constraint["index_0based"]) != index:
            raise RuntimeError(
                f"Constraint-map index mismatch at {index}"
            )

        if int(provenance["index_0based"]) != index:
            raise RuntimeError(
                f"Provenance-map index mismatch at {index}"
            )

        atom_id = constraint["atom_id"]

        if provenance["atom_id"] != atom_id:
            raise RuntimeError(
                f"Atom-ID mismatch at index {index}"
            )

        if not (
            start["element"]
            == final["element"]
            == constraint["element"]
            == provenance["element"]
        ):
            raise RuntimeError(
                f"Element mismatch for {atom_id}"
            )

        mapped[atom_id] = {
            "index": index,
            "element": final["element"],
            "start_xyz_A": start["xyz_A"],
            "final_xyz_A": final["xyz_A"],
            "fixed": parse_bool(
                constraint["v5b_fixed"]
            ),
            "mobile": parse_bool(
                constraint["v5b_mobile"]
            ),
            "artificial_cap": parse_bool(
                provenance["artificial_cap"]
            ),
            "atom_role": provenance["atom_role"],
            "node_type": provenance["node_type"],
        }

    composition = Counter(
        record["element"]
        for record in mapped.values()
    )

    atom_count_gate = (
        len(mapped) == EXPECTED_ATOM_COUNT
    )

    composition_gate = (
        composition == EXPECTED_COMPOSITION
    )

    adjacency = defaultdict(set)
    bonded_pairs = set()

    for row in valence_rows:
        atom_id = row["atom_id"]

        if atom_id not in mapped:
            raise RuntimeError(
                f"Pre-QM valence atom absent from V4: {atom_id}"
            )

        neighbors = [
            value
            for value in row["neighbors"].split("|")
            if value
        ]

        for neighbor in neighbors:
            if neighbor not in mapped:
                raise RuntimeError(
                    "Pre-QM bonded neighbor absent from final "
                    f"model: {atom_id} -- {neighbor}"
                )

            pair = canonical_pair(atom_id, neighbor)
            bonded_pairs.add(pair)

    for first, second in bonded_pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)

    atom_ids = set(mapped)

    visited = set()
    queue = deque([next(iter(atom_ids))])

    while queue:
        current = queue.popleft()

        if current in visited:
            continue

        visited.add(current)

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                queue.append(neighbor)

    connected_gate = (
        visited == atom_ids
    )

    displacement_records = []
    fixed_failures = []

    sum_squared_mobile = 0.0
    mobile_count = 0

    for atom_id, record in sorted(
        mapped.items(),
        key=lambda item: item[1]["index"],
    ):
        displacement_A = distance(
            record["start_xyz_A"],
            record["final_xyz_A"],
        )

        fixed_pass = (
            not record["fixed"]
            or displacement_A
            <= FIXED_DISPLACEMENT_TOLERANCE_A
        )

        if record["fixed"] and not fixed_pass:
            fixed_failures.append(atom_id)

        if record["mobile"]:
            sum_squared_mobile += displacement_A ** 2
            mobile_count += 1

        displacement_records.append({
            "index_0based": record["index"],
            "atom_id": atom_id,
            "element": record["element"],
            "fixed": record["fixed"],
            "mobile": record["mobile"],
            "artificial_cap": record["artificial_cap"],
            "displacement_A": displacement_A,
            "fixed_displacement_tolerance_A": (
                FIXED_DISPLACEMENT_TOLERANCE_A
            ),
            "fixed_displacement_pass": fixed_pass,
        })

    fixed_displacement_gate = (
        len(fixed_failures) == 0
    )

    mobile_rms_displacement_A = (
        math.sqrt(
            sum_squared_mobile / mobile_count
        )
        if mobile_count
        else 0.0
    )

    with ATOM_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                displacement_records[0]
            ),
        )
        writer.writeheader()
        writer.writerows(displacement_records)

    bond_records = []
    bond_failures = []

    for first, second in sorted(bonded_pairs):
        first_record = mapped[first]
        second_record = mapped[second]

        value = distance(
            first_record["final_xyz_A"],
            second_record["final_xyz_A"],
        )

        elements = {
            first_record["element"],
            second_record["element"],
        }

        if elements == {"B", "N"}:
            bond_class = "B-N"
            minimum_A = BN_MIN_A
            maximum_A = BN_MAX_A
        elif elements == {"B", "H"}:
            bond_class = "B-H"
            minimum_A = BH_MIN_A
            maximum_A = BH_MAX_A
        elif elements == {"N", "H"}:
            bond_class = "N-H"
            minimum_A = NH_MIN_A
            maximum_A = NH_MAX_A
        else:
            bond_class = "OTHER"
            minimum_A = 0.0
            maximum_A = float("inf")

        passed = minimum_A <= value <= maximum_A

        record = {
            "first_atom": first,
            "first_element": first_record["element"],
            "second_atom": second,
            "second_element": second_record["element"],
            "bond_class": bond_class,
            "distance_A": value,
            "minimum_A": minimum_A,
            "maximum_A": maximum_A,
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

    with BOND_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(bond_records[0]),
        )
        writer.writeheader()
        writer.writerows(bond_records)

    bond_gate = (
        len(bond_failures) == 0
    )

    cap_records = []
    cap_failures = []

    for row in pre_qm_cap_rows:
        cap_id = row["cap_id"]
        center_id = row["center_atom"]

        if cap_id not in mapped or center_id not in mapped:
            raise RuntimeError(
                f"Missing cap ownership pair: "
                f"{cap_id} -- {center_id}"
            )

        cap = mapped[cap_id]
        center = mapped[center_id]

        owner_distance_A = distance(
            cap["final_xyz_A"],
            center["final_xyz_A"],
        )

        heavy_candidates = sorted(
            (
                distance(
                    cap["final_xyz_A"],
                    record["final_xyz_A"],
                ),
                atom_id,
                record["element"],
            )
            for atom_id, record in mapped.items()
            if record["element"] in {"B", "N"}
        )

        (
            nearest_distance_A,
            nearest_heavy_id,
            nearest_heavy_element,
        ) = heavy_candidates[0]

        second_distance_A = (
            heavy_candidates[1][0]
            if len(heavy_candidates) > 1
            else None
        )

        ownership_pass = (
            nearest_heavy_id == center_id
            and (
                second_distance_A is None
                or second_distance_A
                - nearest_distance_A
                >= CAP_OWNER_MARGIN_A
            )
        )

        if center["element"] == "B":
            bond_range_pass = (
                BH_MIN_A
                <= owner_distance_A
                <= BH_MAX_A
            )
        elif center["element"] == "N":
            bond_range_pass = (
                NH_MIN_A
                <= owner_distance_A
                <= NH_MAX_A
            )
        else:
            bond_range_pass = False

        passed = (
            ownership_pass
            and bond_range_pass
        )

        record = {
            "cap_id": cap_id,
            "center_atom": center_id,
            "center_element": center["element"],
            "owner_distance_A": owner_distance_A,
            "nearest_heavy_atom": nearest_heavy_id,
            "nearest_heavy_element": nearest_heavy_element,
            "nearest_heavy_distance_A": nearest_distance_A,
            "second_heavy_distance_A": second_distance_A,
            "ownership_pass": ownership_pass,
            "bond_range_pass": bond_range_pass,
            "pass": passed,
        }

        cap_records.append(record)

        if not passed:
            cap_failures.append(record)

    with CAP_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(cap_records[0]),
        )
        writer.writeheader()
        writer.writerows(cap_records)

    construction_cap_ids = sorted(
        row["cap_id"]
        for row in construction_cap_rows
    )

    pre_qm_cap_ids = sorted(
        row["cap_id"]
        for row in pre_qm_cap_rows
    )

    post_qm_cap_ids = sorted(
        row["cap_id"]
        for row in cap_records
    )

    new_cap_inventory_gate = (
        len(construction_cap_ids) > 0
        and construction_cap_ids == pre_qm_cap_ids
        and pre_qm_cap_ids == post_qm_cap_ids
    )

    new_cap_ownership_gate = (
        len(cap_records) == len(construction_cap_ids)
        and all(
            record["ownership_pass"]
            for record in cap_records
        )
    )

    new_cap_bond_range_gate = (
        len(cap_records) == len(construction_cap_ids)
        and all(
            record["bond_range_pass"]
            for record in cap_records
        )
    )

    cap_gate = (
        new_cap_inventory_gate
        and new_cap_ownership_gate
        and new_cap_bond_range_gate
        and len(cap_failures) == 0
    )

    graph_distances = {}

    for origin in atom_ids:
        distances = {origin: 0}
        local_queue = deque([origin])

        while local_queue:
            current = local_queue.popleft()

            for neighbor in adjacency[current]:
                if neighbor not in distances:
                    distances[neighbor] = (
                        distances[current] + 1
                    )
                    local_queue.append(neighbor)

        for target, separation in distances.items():
            graph_distances[
                canonical_pair(origin, target)
            ] = separation

    contact_records = []
    hard_contacts = []

    ordered = sorted(
        mapped.items(),
        key=lambda item: item[1]["index"],
    )

    for position, (first_id, first) in enumerate(ordered):
        for second_id, second in ordered[position + 1:]:
            value = distance(
                first["final_xyz_A"],
                second["final_xyz_A"],
            )

            separation = graph_distances.get(
                canonical_pair(first_id, second_id)
            )

            if separation in {1, 2}:
                classification = "TOPOLOGICAL_EXCLUSION"
                is_hard = False
            else:
                if (
                    first["element"] == "H"
                    and second["element"] == "H"
                ):
                    threshold_A = HH_HARD_CONTACT_A
                elif (
                    first["element"] == "H"
                    or second["element"] == "H"
                ):
                    threshold_A = HX_HARD_CONTACT_A
                else:
                    threshold_A = (
                        HEAVY_HEAVY_HARD_CONTACT_A
                    )

                is_hard = value < threshold_A

                classification = (
                    "HARD_CONTACT"
                    if is_hard
                    else "NONBONDED_OK"
                )

            record = {
                "first_atom": first_id,
                "first_element": first["element"],
                "second_atom": second_id,
                "second_element": second["element"],
                "distance_A": value,
                "graph_separation": separation,
                "classification": classification,
                "hard_contact": is_hard,
            }

            contact_records.append(record)

            if is_hard:
                hard_contacts.append(record)

    with CONTACT_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(contact_records[0]),
        )
        writer.writeheader()
        writer.writerows(contact_records)

    contact_gate = (
        len(hard_contacts) == 0
    )

    with FINAL_XYZ_PATH.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(f"{len(final_atoms)}\n")
        handle.write(
            "QM_F06 UPPER V5-B latest audited frame; "
            f"final={execution_complete}\n"
        )

        for atom in final_atoms:
            x_value, y_value, z_value = atom["xyz_A"]

            handle.write(
                f"{atom['element']:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    structural_gates = {
        "atom_count": atom_count_gate,
        "composition": composition_gate,
        "single_connected_component": connected_gate,
        "fixed_atoms_unchanged": fixed_displacement_gate,
        "all_validated_bonds_preserved": bond_gate,
        "new_cap_inventory_consistent": (
            new_cap_inventory_gate
        ),
        "all_new_caps_owned": (
            new_cap_ownership_gate
        ),
        "all_new_cap_bonds_in_range": (
            new_cap_bond_range_gate
        ),
        "no_unresolved_hard_contacts": contact_gate,
    }

    structural_pass = all(
        structural_gates.values()
    )

    final_acceptance = (
        execution_complete
        and structural_pass
    )

    if final_acceptance:
        decision = (
            "QM_F06_UPPER_V5B_POST_QM_GATE_PASS_"
            "GEOMETRIC_REFERENCE_ACCEPTED_"
            "RESP_PREPARATION_AUTHORIZED"
        )
    elif args.allow_incomplete:
        decision = (
            "QM_F06_UPPER_V5B_POST_QM_INTERIM_"
            "DIAGNOSTIC_ONLY_NO_AUTHORIZATION"
        )
    else:
        decision = (
            "QM_F06_UPPER_V5B_POST_QM_GATE_FAIL_"
            "GEOMETRY_REVIEW_REQUIRED"
        )

    files_for_hash = {
        "orca_output": output_path,
        "start_xyz": start_xyz_path,
        "geometry_source": geometry_source,
        "constraint_map": CONSTRAINT_MAP,
        "provenance_map": PROVENANCE_MAP,
        "construction_caps": CONSTRUCTION_CAPS,
        "pre_qm_report": PRE_QM_REPORT,
        "pre_qm_valence": PRE_QM_VALENCE,
        "pre_qm_caps": PRE_QM_CAPS,
        "execution_manifest": manifest_path,
        "final_xyz": FINAL_XYZ_PATH,
        "atom_audit_csv": ATOM_AUDIT_CSV,
        "bond_audit_csv": BOND_AUDIT_CSV,
        "cap_audit_csv": CAP_AUDIT_CSV,
        "contact_audit_csv": CONTACT_AUDIT_CSV,
    }

    if stderr_path.is_file():
        files_for_hash["orca_stderr"] = stderr_path

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "diagnostic_mode": bool(
            args.allow_incomplete
            and not execution_complete
        ),
        "orca_status": status,
        "geometry": {
            "source": str(
                geometry_source.relative_to(ROOT)
            ),
            "trajectory_frame_index": (
                trajectory_frame_index
            ),
            "atom_count": len(final_atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "mobile_rms_displacement_A": (
                mobile_rms_displacement_A
            ),
            "fixed_displacement_failure_count": (
                len(fixed_failures)
            ),
            "fixed_displacement_failures": (
                fixed_failures
            ),
        },
        "structural_gates": structural_gates,
        "structural_pass": structural_pass,
        "bond_failure_count": len(
            bond_failures
        ),
        "cap_failure_count": len(
            cap_failures
        ),
        "hard_contact_count": len(
            hard_contacts
        ),
        "new_cap_inventory": {
            "construction_cap_count": len(
                construction_cap_ids
            ),
            "pre_qm_cap_count": len(
                pre_qm_cap_ids
            ),
            "post_qm_cap_count": len(
                post_qm_cap_ids
            ),
            "construction_cap_ids": (
                construction_cap_ids
            ),
            "pre_qm_cap_ids": pre_qm_cap_ids,
            "post_qm_cap_ids": post_qm_cap_ids,
            "inventory_consistent": (
                new_cap_inventory_gate
            ),
            "all_owned": (
                new_cap_ownership_gate
            ),
            "all_bonds_in_range": (
                new_cap_bond_range_gate
            ),
            "pass": cap_gate,
        },
        "final_acceptance": final_acceptance,
        "files": {
            key: str(path.relative_to(ROOT))
            for key, path in files_for_hash.items()
        },
        "files_sha256": {
            key: sha256(path)
            for key, path in files_for_hash.items()
        },
        "authorization": {
            "geometric_reference_accepted": (
                final_acceptance
            ),
            "RESP_input_preparation_authorized": (
                final_acceptance
            ),
            "RESP_execution_authorized": False,
            "electronic_reference_accepted": False,
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

    print("=" * 78)
    print("QM_F06 UPPER V5-B POST-QM AUDIT")
    print("=" * 78)

    print("Execution:", execution_dir)
    print("Geometry source:", geometry_source)
    print(
        "Trajectory frame index:",
        trajectory_frame_index,
    )
    print()

    print(
        "ORCA normal termination:",
        status["normal_termination"],
    )
    print(
        "Geometry converged:",
        status["geometry_converged"],
    )
    print(
        "Error detected:",
        status["error_detected"],
    )
    print(
        "Geometry cycles:",
        status["geometry_cycles_started"],
    )
    print(
        "SCF convergences:",
        status["scf_convergences"],
    )
    print(
        "Energy points:",
        status["completed_energy_points"],
    )
    print(
        "Last energy Eh:",
        status["last_completed_energy_hartree"],
    )

    print()

    for gate, passed in structural_gates.items():
        print(
            f"{gate:40s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Fixed displacement failures:", len(
        fixed_failures
    ))
    print("Bond failures:", len(bond_failures))
    print("Cap failures:", len(cap_failures))
    print("Hard contacts:", len(hard_contacts))
    print(
        "New-cap inventory counts "
        "(construction/pre-QM/post-QM):",
        len(construction_cap_ids),
        len(pre_qm_cap_ids),
        len(post_qm_cap_ids),
    )
    print(
        "New-cap inventories consistent:",
        new_cap_inventory_gate,
    )
    print(
        "All new caps owned:",
        new_cap_ownership_gate,
    )
    print(
        "All new cap bonds in range:",
        new_cap_bond_range_gate,
    )
    print(
        "Mobile RMS displacement A:",
        mobile_rms_displacement_A,
    )

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Final/latest XYZ:", FINAL_XYZ_PATH)
    print("Atom audit:", ATOM_AUDIT_CSV)
    print("Bond audit:", BOND_AUDIT_CSV)
    print("Cap audit:", CAP_AUDIT_CSV)
    print("Contact audit:", CONTACT_AUDIT_CSV)
    print()
    print(
        "RESP input preparation authorized:",
        final_acceptance,
    )
    print("RESP execution authorized: False")
    print("Force-field adoption authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
