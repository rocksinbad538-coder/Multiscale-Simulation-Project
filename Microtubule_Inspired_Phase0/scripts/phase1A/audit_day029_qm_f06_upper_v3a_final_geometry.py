#!/usr/bin/env python3
"""
Final structural auditor for QM_F06 UPPER Boundary V3-A.

This script:
- locates the latest restart4 execution;
- requires normal ORCA termination and geometry convergence;
- reads the initial and optimized 30-atom geometries;
- validates composition, atom ordering, connectivity, local valence,
  bond lengths, nonbonded hard contacts, fixed-atom preservation,
  and mobile-atom displacements;
- writes reproducible JSON/CSV/XYZ audit artifacts;
- does not authorize RESP, force-field adoption, or MD.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_ROOT = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_upper_boundary_v3a_restart4_executions"
)

WORKFLOW_ROOT = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3_workflow"
)

RESTART_ROOT = WORKFLOW_ROOT / "restart4"

INITIAL_XYZ = RESTART_ROOT / "v3a_restart4_start.xyz"
INPUT_PATH = RESTART_ROOT / "v3a_restart4.inp"
STATE_PATH = WORKFLOW_ROOT / "v3_workflow_state.json"

AUDIT_ROOT = ROOT / (
    "runs/phase1A/day029_qm_f06_upper_v3a_final_audit"
)

EXPECTED_COMPOSITION = {
    "B": 8,
    "N": 7,
    "H": 15,
}

# Covalent graph thresholds used only to identify chemically plausible bonds.
# They are deliberately broader than ideal equilibrium distances.
BOND_WINDOWS = {
    frozenset(("B", "N")): (1.20, 1.75),
    frozenset(("B", "H")): (0.95, 1.35),
    frozenset(("N", "H")): (0.85, 1.25),
}

# Hard-contact thresholds for pairs that are not assigned as covalent.
NONBONDED_HARD_CONTACT = {
    frozenset(("H", "H")): 1.20,
    frozenset(("B", "H")): 1.30,
    frozenset(("N", "H")): 1.25,
    frozenset(("B", "B")): 1.75,
    frozenset(("N", "N")): 1.75,
    frozenset(("B", "N")): 1.80,
}

EXPECTED_VALENCE = {
    "B": {3},
    "N": {3},
    "H": {1},
}

FIXED_DISPLACEMENT_TOLERANCE_A = 1.0e-5


@dataclass(frozen=True)
class Atom:
    index: int
    element: str
    x: float
    y: float
    z: float

    @property
    def xyz(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def distance(a: Atom, b: Atom) -> float:
    return math.dist(a.xyz, b.xyz)


def read_xyz(path: Path) -> tuple[str, list[Atom]]:
    require_file(path)
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(f"Malformed XYZ: {path}")

    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise RuntimeError(
            f"Invalid atom count in {path}"
        ) from exc

    atom_lines = lines[2:2 + atom_count]

    if len(atom_lines) != atom_count:
        raise RuntimeError(
            f"Incomplete XYZ: expected {atom_count}, "
            f"found {len(atom_lines)}"
        )

    atoms: list[Atom] = []

    for index, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed atom line {index}: {line}"
            )

        atoms.append(
            Atom(
                index=index,
                element=fields[0],
                x=float(fields[1]),
                y=float(fields[2]),
                z=float(fields[3]),
            )
        )

    return lines[1].strip(), atoms


def write_xyz(
    path: Path,
    comment: str,
    atoms: Iterable[Atom],
) -> None:
    atoms = list(atoms)

    lines = [
        str(len(atoms)),
        comment,
    ]

    for atom in atoms:
        lines.append(
            f"{atom.element:2s} "
            f"{atom.x: .10f} "
            f"{atom.y: .10f} "
            f"{atom.z: .10f}"
        )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_fixed_indices(input_path: Path) -> list[int]:
    require_file(input_path)
    text = input_path.read_text(encoding="utf-8")

    indices = [
        int(value)
        for value in re.findall(
            r"\{\s*C\s+(\d+)\s+C\s*\}",
            text,
        )
    ]

    if not indices:
        raise RuntimeError(
            "No fixed Cartesian constraints found."
        )

    if len(indices) != len(set(indices)):
        raise RuntimeError(
            "Duplicate fixed-atom constraint indices."
        )

    return sorted(indices)


def latest_execution() -> Path:
    candidates = sorted(
        path
        for path in EXECUTION_ROOT.glob("restart4_*")
        if path.is_dir()
    )

    if not candidates:
        raise RuntimeError(
            "No restart4 execution directory found."
        )

    return candidates[-1]


def select_final_xyz(execution_dir: Path) -> Path:
    candidates = [
        execution_dir / "restart4.xyz",
        execution_dir / "restart4_trj.xyz",
    ]

    final_xyz = candidates[0]

    if final_xyz.is_file() and final_xyz.stat().st_size > 0:
        return final_xyz

    trajectory = candidates[1]

    if not trajectory.is_file():
        raise RuntimeError(
            "No optimized XYZ or trajectory found."
        )

    return extract_last_complete_frame(
        trajectory,
        execution_dir / "restart4_last_frame.xyz",
    )


def extract_last_complete_frame(
    trajectory: Path,
    output_path: Path,
) -> Path:
    lines = trajectory.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    frames: list[tuple[str, list[str]]] = []
    cursor = 0

    while cursor < len(lines):
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError:
            cursor += 1
            continue

        end = cursor + 2 + atom_count

        if end > len(lines):
            break

        comment = lines[cursor + 1]
        atom_lines = lines[cursor + 2:end]

        if all(len(line.split()) >= 4 for line in atom_lines):
            frames.append((comment, atom_lines))

        cursor = end

    if not frames:
        raise RuntimeError(
            "No complete trajectory frame available."
        )

    comment, atom_lines = frames[-1]

    output_path.write_text(
        "\n".join(
            [
                str(len(atom_lines)),
                comment,
                *atom_lines,
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return output_path


def build_graph(
    atoms: list[Atom],
) -> tuple[
    list[dict],
    dict[int, set[int]],
]:
    bonds: list[dict] = []
    adjacency: dict[int, set[int]] = defaultdict(set)

    for i, atom_i in enumerate(atoms):
        for atom_j in atoms[i + 1:]:
            pair = frozenset(
                (atom_i.element, atom_j.element)
            )

            if pair not in BOND_WINDOWS:
                continue

            lower, upper = BOND_WINDOWS[pair]
            value = distance(atom_i, atom_j)

            if lower <= value <= upper:
                bonds.append(
                    {
                        "i": atom_i.index,
                        "j": atom_j.index,
                        "element_i": atom_i.element,
                        "element_j": atom_j.element,
                        "distance_A": value,
                    }
                )

                adjacency[atom_i.index].add(atom_j.index)
                adjacency[atom_j.index].add(atom_i.index)

    for atom in atoms:
        adjacency.setdefault(atom.index, set())

    return bonds, adjacency


def is_connected(
    adjacency: dict[int, set[int]],
) -> bool:
    if not adjacency:
        return False

    start = next(iter(adjacency))
    visited = {start}
    queue = deque([start])

    while queue:
        current = queue.popleft()

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    return len(visited) == len(adjacency)


def audit_valence(
    atoms: list[Atom],
    adjacency: dict[int, set[int]],
) -> list[dict]:
    failures: list[dict] = []

    for atom in atoms:
        observed = len(adjacency[atom.index])
        allowed = EXPECTED_VALENCE.get(
            atom.element,
            set(),
        )

        if observed not in allowed:
            failures.append(
                {
                    "index": atom.index,
                    "element": atom.element,
                    "observed_valence": observed,
                    "allowed_valence": sorted(allowed),
                    "neighbors": sorted(
                        adjacency[atom.index]
                    ),
                }
            )

    return failures


def audit_nonbonded_contacts(
    atoms: list[Atom],
    adjacency: dict[int, set[int]],
) -> list[dict]:
    contacts: list[dict] = []

    for i, atom_i in enumerate(atoms):
        for atom_j in atoms[i + 1:]:
            if atom_j.index in adjacency[atom_i.index]:
                continue

            pair = frozenset(
                (atom_i.element, atom_j.element)
            )

            threshold = NONBONDED_HARD_CONTACT.get(pair)

            if threshold is None:
                continue

            value = distance(atom_i, atom_j)

            if value < threshold:
                contacts.append(
                    {
                        "i": atom_i.index,
                        "j": atom_j.index,
                        "element_i": atom_i.element,
                        "element_j": atom_j.element,
                        "distance_A": value,
                        "hard_contact_threshold_A": threshold,
                    }
                )

    return contacts


def displacement_rows(
    initial: list[Atom],
    final: list[Atom],
    fixed_indices: set[int],
) -> list[dict]:
    rows: list[dict] = []

    if len(initial) != len(final):
        raise RuntimeError(
            "Initial and final atom counts differ."
        )

    for atom_initial, atom_final in zip(initial, final):
        if atom_initial.element != atom_final.element:
            raise RuntimeError(
                "Element ordering changed at atom "
                f"{atom_initial.index}: "
                f"{atom_initial.element} -> "
                f"{atom_final.element}"
            )

        displacement = distance(
            atom_initial,
            atom_final,
        )

        rows.append(
            {
                "index": atom_initial.index,
                "element": atom_initial.element,
                "fixed": (
                    atom_initial.index in fixed_indices
                ),
                "displacement_A": displacement,
                "initial_x_A": atom_initial.x,
                "initial_y_A": atom_initial.y,
                "initial_z_A": atom_initial.z,
                "final_x_A": atom_final.x,
                "final_y_A": atom_final.y,
                "final_z_A": atom_final.z,
            }
        )

    return rows


def rms(values: Iterable[float]) -> float:
    values = list(values)

    if not values:
        return 0.0

    return math.sqrt(
        sum(value * value for value in values)
        / len(values)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help=(
            "Allow a diagnostic dry audit of the latest frame "
            "without accepting it as a final geometry."
        ),
    )
    args = parser.parse_args()

    execution_dir = latest_execution()
    output_path = execution_dir / "restart4.out"
    stderr_path = execution_dir / "restart4.stderr"

    require_file(output_path)
    require_file(INITIAL_XYZ)
    require_file(INPUT_PATH)

    output_text = output_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    normal_termination = (
        "ORCA TERMINATED NORMALLY" in output_text
    )

    geometry_converged = (
        "THE OPTIMIZATION HAS CONVERGED"
        in output_text
    )

    execution_active_or_incomplete = not (
        normal_termination
        and geometry_converged
    )

    if (
        execution_active_or_incomplete
        and not args.allow_incomplete
    ):
        raise RuntimeError(
            "Final audit blocked: ORCA has not yet "
            "terminated normally with a converged geometry."
        )

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    audit_dir = AUDIT_ROOT / f"audit_{timestamp}"
    audit_dir.mkdir(parents=True, exist_ok=False)

    final_xyz = select_final_xyz(execution_dir)

    initial_comment, initial_atoms = read_xyz(
        INITIAL_XYZ
    )
    final_comment, final_atoms = read_xyz(
        final_xyz
    )

    fixed_indices = set(
        parse_fixed_indices(INPUT_PATH)
    )

    composition = Counter(
        atom.element
        for atom in final_atoms
    )

    bonds, adjacency = build_graph(final_atoms)
    valence_failures = audit_valence(
        final_atoms,
        adjacency,
    )
    hard_contacts = audit_nonbonded_contacts(
        final_atoms,
        adjacency,
    )

    displacements = displacement_rows(
        initial_atoms,
        final_atoms,
        fixed_indices,
    )

    fixed_displacements = [
        row["displacement_A"]
        for row in displacements
        if row["fixed"]
    ]

    mobile_displacements = [
        row["displacement_A"]
        for row in displacements
        if not row["fixed"]
    ]

    fixed_failures = [
        row
        for row in displacements
        if (
            row["fixed"]
            and row["displacement_A"]
            > FIXED_DISPLACEMENT_TOLERANCE_A
        )
    ]

    energy_values = [
        float(value)
        for value in re.findall(
            r"FINAL SINGLE POINT ENERGY\s+"
            r"(-?\d+\.\d+)",
            output_text,
        )
    ]

    checks = {
        "normal_termination": normal_termination,
        "geometry_converged": geometry_converged,
        "atom_count_30": len(final_atoms) == 30,
        "composition_B8_N7_H15": (
            dict(composition)
            == EXPECTED_COMPOSITION
        ),
        "element_order_preserved": all(
            a.element == b.element
            for a, b in zip(
                initial_atoms,
                final_atoms,
            )
        ),
        "fixed_atom_count_20": (
            len(fixed_indices) == 20
        ),
        "graph_connected": is_connected(adjacency),
        "valence_failures_zero": (
            len(valence_failures) == 0
        ),
        "hard_contacts_zero": (
            len(hard_contacts) == 0
        ),
        "fixed_displacement_failures_zero": (
            len(fixed_failures) == 0
        ),
        "final_energy_available": bool(
            energy_values
        ),
    }

    final_gate_pass = all(checks.values())

    decision = (
        "QM_F06_UPPER_V3A_FINAL_STRUCTURAL_AUDIT_PASS"
        if final_gate_pass
        else (
            "QM_F06_UPPER_V3A_INTERIM_DIAGNOSTIC_ONLY"
            if args.allow_incomplete
            else
            "QM_F06_UPPER_V3A_FINAL_STRUCTURAL_AUDIT_FAIL"
        )
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "orca_output": str(
            output_path.relative_to(ROOT)
        ),
        "orca_stderr": (
            str(stderr_path.relative_to(ROOT))
            if stderr_path.is_file()
            else None
        ),
        "initial_xyz": str(
            INITIAL_XYZ.relative_to(ROOT)
        ),
        "final_xyz_source": str(
            final_xyz.relative_to(ROOT)
        ),
        "initial_comment": initial_comment,
        "final_comment": final_comment,
        "checks": checks,
        "final_gate_pass": final_gate_pass,
        "diagnostic_only": bool(
            args.allow_incomplete
            and not final_gate_pass
        ),
        "atom_count": len(final_atoms),
        "composition": dict(composition),
        "bond_count": len(bonds),
        "graph_connected": is_connected(adjacency),
        "fixed_atom_count": len(fixed_indices),
        "mobile_atom_count": (
            len(final_atoms) - len(fixed_indices)
        ),
        "fixed_displacement_max_A": (
            max(fixed_displacements)
            if fixed_displacements
            else None
        ),
        "fixed_displacement_rms_A": rms(
            fixed_displacements
        ),
        "mobile_displacement_max_A": (
            max(mobile_displacements)
            if mobile_displacements
            else None
        ),
        "mobile_displacement_rms_A": rms(
            mobile_displacements
        ),
        "valence_failures": valence_failures,
        "hard_contacts": hard_contacts,
        "fixed_displacement_failures": (
            fixed_failures
        ),
        "final_energy_hartree": (
            energy_values[-1]
            if energy_values
            else None
        ),
        "completed_energy_points": len(
            energy_values
        ),
        "authorization": {
            "geometry_reference_accepted": (
                final_gate_pass
            ),
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    report_path = (
        audit_dir
        / "QM_F06_UPPER_V3A_STRUCTURAL_AUDIT.json"
    )

    report_path.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    bond_csv = audit_dir / "bond_inventory.csv"

    with bond_csv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "i",
                "j",
                "element_i",
                "element_j",
                "distance_A",
            ],
        )
        writer.writeheader()
        writer.writerows(bonds)

    displacement_csv = (
        audit_dir / "atom_displacements.csv"
    )

    with displacement_csv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                displacements[0].keys()
            ),
        )
        writer.writeheader()
        writer.writerows(displacements)

    audited_xyz = (
        audit_dir
        / "QM_F06_UPPER_V3A_AUDITED_GEOMETRY.xyz"
    )

    write_xyz(
        audited_xyz,
        (
            f"{decision}; "
            f"source={final_xyz.name}"
        ),
        final_atoms,
    )

    manifest_entries = {}

    for artifact in (
        report_path,
        bond_csv,
        displacement_csv,
        audited_xyz,
        output_path,
        INPUT_PATH,
        INITIAL_XYZ,
    ):
        manifest_entries[
            str(artifact.relative_to(ROOT))
        ] = sha256(artifact)

    manifest = {
        "decision": (
            "QM_F06_UPPER_V3A_AUDIT_MANIFEST"
        ),
        "final_gate_pass": final_gate_pass,
        "files_sha256": manifest_entries,
    }

    manifest_path = (
        audit_dir
        / "QM_F06_UPPER_V3A_AUDIT_MANIFEST.json"
    )

    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V3-A STRUCTURAL AUDIT")
    print("=" * 78)
    print("Decision:", decision)
    print("Execution:", execution_dir)
    print("Final XYZ source:", final_xyz)
    print("Normal termination:", normal_termination)
    print("Geometry converged:", geometry_converged)
    print("Atoms:", len(final_atoms))
    print("Composition:", dict(composition))
    print("Bonds:", len(bonds))
    print("Connected:", is_connected(adjacency))
    print("Valence failures:", len(valence_failures))
    print("Hard contacts:", len(hard_contacts))
    print(
        "Fixed displacement failures:",
        len(fixed_failures),
    )
    print(
        "Final energy (Eh):",
        (
            energy_values[-1]
            if energy_values
            else None
        ),
    )
    print("Final structural gate:", final_gate_pass)
    print("Report:", report_path)
    print("Manifest:", manifest_path)


if __name__ == "__main__":
    main()
