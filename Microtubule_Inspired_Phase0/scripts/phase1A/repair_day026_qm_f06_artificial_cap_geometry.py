#!/usr/bin/env python3
"""
Targeted steric repair of artificial H caps in QM_F06 fragments.

Only artificial caps participating in topology-aware severe clashes or
strong compressions are moved.

Constraints
-----------
- Parent B/N atom remains fixed.
- X-H bond length remains fixed.
- All original R2 atoms remain fixed.
- Candidate cap directions are sampled within a cone around the original
  cut-bond direction.
- The selected direction maximizes the minimum normalized nonbonded
  distance involving the cap.
- No QM calculation or geometry optimization is executed.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

MAX_CONE_ANGLE_DEG = 75.0
CONE_ANGLE_STEP_DEG = 2.5
AZIMUTH_STEP_DEG = 5.0

TARGET_MIN_RATIO = 0.80


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def vector(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        b - a
        for a, b in zip(first, second, strict=True)
    )


def norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(sum(component * component for component in v))


def normalize(
    v: tuple[float, float, float],
) -> tuple[float, float, float]:
    length = norm(v)

    if length <= 1.0e-14:
        raise RuntimeError("Cannot normalize zero-length vector.")

    return tuple(component / length for component in v)


def dot(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return sum(
        a * b
        for a, b in zip(first, second, strict=True)
    )


def cross(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def add(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        a + b
        for a, b in zip(first, second, strict=True)
    )


def scale(
    v: tuple[float, float, float],
    factor: float,
) -> tuple[float, float, float]:
    return tuple(component * factor for component in v)


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return norm(vector(first, second))


def coordinates(
    atom: dict[str, str],
) -> tuple[float, float, float]:
    return (
        float(atom["x_angstrom"]),
        float(atom["y_angstrom"]),
        float(atom["z_angstrom"]),
    )


def orthonormal_basis(
    axis: tuple[float, float, float],
) -> tuple[
    tuple[float, float, float],
    tuple[float, float, float],
]:
    axis = normalize(axis)

    reference = (
        (1.0, 0.0, 0.0)
        if abs(axis[0]) < 0.85
        else (0.0, 1.0, 0.0)
    )

    first = normalize(cross(axis, reference))
    second = normalize(cross(axis, first))

    return first, second


def candidate_direction(
    axis: tuple[float, float, float],
    basis_1: tuple[float, float, float],
    basis_2: tuple[float, float, float],
    cone_angle_rad: float,
    azimuth_rad: float,
) -> tuple[float, float, float]:
    radial = add(
        scale(basis_1, math.cos(azimuth_rad)),
        scale(basis_2, math.sin(azimuth_rad)),
    )

    direction = add(
        scale(axis, math.cos(cone_angle_rad)),
        scale(radial, math.sin(cone_angle_rad)),
    )

    return normalize(direction)


def write_xyz(
    path: Path,
    atoms: list[dict[str, Any]],
    comment: str,
) -> None:
    center = (
        sum(float(row["x_angstrom"]) for row in atoms) / len(atoms),
        sum(float(row["y_angstrom"]) for row in atoms) / len(atoms),
        sum(float(row["z_angstrom"]) for row in atoms) / len(atoms),
    )

    lines = [str(len(atoms)), comment]

    for atom in atoms:
        lines.append(
            f"{atom['element']:<2s} "
            f"{float(atom['x_angstrom']) - center[0]: .10f} "
            f"{float(atom['y_angstrom']) - center[1]: .10f} "
            f"{float(atom['z_angstrom']) - center[2]: .10f}"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    combined_audit = read_csv(
        F06_DIR
        / "QM_F06_topology_aware_steric_audit_combined.csv"
    )

    repair_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []
    all_cap_gates_pass = True

    for end in ENDS:
        label = f"QM_F06_{end}_CAPPED"

        atoms = read_csv(
            F06_DIR / f"{label}_atoms.csv"
        )
        caps = read_csv(
            F06_DIR / f"{label}_caps.csv"
        )
        valence_rows = read_csv(
            F06_DIR / f"{label}_valence_audit.csv"
        )

        atom_lookup = {
            row["atom_id"]: row
            for row in atoms
        }

        adjacency: dict[str, set[str]] = defaultdict(set)

        for row in valence_rows:
            atom_id = row["atom_id"]

            neighbors = [
                item
                for item in row["neighbors"].split("|")
                if item
            ]

            adjacency[atom_id].update(neighbors)

        conflicting_caps = sorted(
            {
                atom_id
                for row in combined_audit
                if row["fragment"] == label
                and row["requires_geometry_repair"].lower() == "true"
                and row["involves_artificial_cap"].lower() == "true"
                for atom_id in (row["atom_1"], row["atom_2"])
                if atom_id.startswith("HCAP:")
            }
        )

        cap_lookup = {
            row["cap_id"]: row
            for row in caps
        }

        repaired_caps = 0

        for cap_id in conflicting_caps:
            cap_meta = cap_lookup[cap_id]
            cap_atom = atom_lookup[cap_id]

            parent_id = cap_meta["parent_inside_node"]
            parent_atom = atom_lookup[parent_id]

            parent_xyz = coordinates(parent_atom)
            original_cap_xyz = coordinates(cap_atom)

            original_axis = normalize(
                vector(parent_xyz, original_cap_xyz)
            )

            basis_1, basis_2 = orthonormal_basis(
                original_axis
            )

            bond_length = float(
                cap_meta["target_XH_distance_angstrom"]
            )

            best_score: tuple[float, float, float] | None = None
            best_xyz: tuple[float, float, float] | None = None
            best_angle = 0.0
            best_azimuth = 0.0
            best_min_partner = ""

            angle_steps = int(
                round(
                    MAX_CONE_ANGLE_DEG
                    / CONE_ANGLE_STEP_DEG
                )
            )

            azimuth_steps = int(
                round(360.0 / AZIMUTH_STEP_DEG)
            )

            for angle_index in range(angle_steps + 1):
                cone_deg = (
                    angle_index * CONE_ANGLE_STEP_DEG
                )
                cone_rad = math.radians(cone_deg)

                for azimuth_index in range(azimuth_steps):
                    azimuth_deg = (
                        azimuth_index * AZIMUTH_STEP_DEG
                    )
                    azimuth_rad = math.radians(
                        azimuth_deg
                    )

                    direction = candidate_direction(
                        original_axis,
                        basis_1,
                        basis_2,
                        cone_rad,
                        azimuth_rad,
                    )

                    candidate_xyz = add(
                        parent_xyz,
                        scale(direction, bond_length),
                    )

                    minimum_ratio = float("inf")
                    minimum_distance = float("inf")
                    minimum_partner = ""

                    for other_id, other_atom in atom_lookup.items():
                        if other_id in {cap_id, parent_id}:
                            continue

                        # Exclude 1–3 contacts through the parent atom.
                        if other_id in adjacency[parent_id]:
                            continue

                        measured = distance(
                            candidate_xyz,
                            coordinates(other_atom),
                        )

                        vdw_sum = (
                            VDW_RADII["H"]
                            + VDW_RADII[other_atom["element"]]
                        )

                        ratio = measured / vdw_sum

                        if ratio < minimum_ratio:
                            minimum_ratio = ratio
                            minimum_distance = measured
                            minimum_partner = other_id

                    # Primary objective: maximize minimum ratio.
                    # Secondary objective: minimize angular displacement.
                    score = (
                        minimum_ratio,
                        -cone_deg,
                        minimum_distance,
                    )

                    if best_score is None or score > best_score:
                        best_score = score
                        best_xyz = candidate_xyz
                        best_angle = cone_deg
                        best_azimuth = azimuth_deg
                        best_min_partner = minimum_partner

            if best_score is None or best_xyz is None:
                raise RuntimeError(
                    f"No candidate geometry generated for {cap_id}"
                )

            old_xyz = original_cap_xyz
            old_direction = normalize(
                vector(parent_xyz, old_xyz)
            )
            new_direction = normalize(
                vector(parent_xyz, best_xyz)
            )

            angular_change = math.degrees(
                math.acos(
                    max(
                        -1.0,
                        min(1.0, dot(old_direction, new_direction)),
                    )
                )
            )

            cap_atom["x_angstrom"] = f"{best_xyz[0]:.10f}"
            cap_atom["y_angstrom"] = f"{best_xyz[1]:.10f}"
            cap_atom["z_angstrom"] = f"{best_xyz[2]:.10f}"

            passed = best_score[0] >= TARGET_MIN_RATIO
            all_cap_gates_pass = (
                all_cap_gates_pass and passed
            )

            repair_rows.append(
                {
                    "fragment": label,
                    "cap_id": cap_id,
                    "parent_atom": parent_id,
                    "old_x_angstrom": f"{old_xyz[0]:.10f}",
                    "old_y_angstrom": f"{old_xyz[1]:.10f}",
                    "old_z_angstrom": f"{old_xyz[2]:.10f}",
                    "new_x_angstrom": f"{best_xyz[0]:.10f}",
                    "new_y_angstrom": f"{best_xyz[1]:.10f}",
                    "new_z_angstrom": f"{best_xyz[2]:.10f}",
                    "bond_length_angstrom": f"{bond_length:.10f}",
                    "angular_change_deg": f"{angular_change:.6f}",
                    "search_cone_angle_deg": f"{best_angle:.6f}",
                    "search_azimuth_deg": f"{best_azimuth:.6f}",
                    "minimum_post_repair_vdw_ratio": (
                        f"{best_score[0]:.10f}"
                    ),
                    "minimum_post_repair_partner": (
                        best_min_partner
                    ),
                    "target_minimum_ratio": (
                        TARGET_MIN_RATIO
                    ),
                    "cap_repair_gate_pass": passed,
                }
            )

            repaired_caps += 1

        repaired_atoms_path = (
            F06_DIR / f"{label}_REPAIRED_atoms.csv"
        )
        repaired_xyz_path = (
            F06_DIR / f"{label}_REPAIRED.xyz"
        )

        write_csv(repaired_atoms_path, atoms)

        write_xyz(
            repaired_xyz_path,
            atoms,
            (
                f"{label}_REPAIRED; only conflicting artificial caps "
                "were angularly repositioned; original R2 atoms fixed; "
                "unoptimized; no QM calculation executed"
            ),
        )

        report_sections.extend(
            [
                f"## {label}",
                "",
                (
                    f"- Conflicting artificial caps identified: "
                    f"**{len(conflicting_caps)}**"
                ),
                f"- Caps repositioned: **{repaired_caps}**",
                "- Original R2 atoms moved: **0**",
                "- X–H bond lengths changed: **NO**",
                "",
            ]
        )

    write_csv(
        F06_DIR / "QM_F06_artificial_cap_repair_manifest.csv",
        repair_rows,
    )

    decision = (
        "QM_F06_ARTIFICIAL_CAP_GEOMETRY_REPAIR_COMPLETED"
        if all_cap_gates_pass
        else "QM_F06_CAP_SEARCH_INSUFFICIENT_FRAGMENT_REDESIGN_REQUIRED"
    )

    (
        F06_DIR / "QM_F06_ARTIFICIAL_CAP_REPAIR_REPORT.md"
    ).write_text(
        "\n".join(
            [
                "# QM_F06 Artificial-Cap Geometry Repair — Day026",
                "",
                "## Scope",
                "",
                (
                    "Only artificial caps involved in topology-aware "
                    "hard steric failures were repositioned."
                ),
                "",
                (
                    "All original R2 atoms, bridge atoms, attachment "
                    "centers and existing hydrogen atoms remained fixed."
                ),
                "",
                f"## Decision: **{decision}**",
                "",
                *report_sections,
                "## Search constraints",
                "",
                (
                    f"- Maximum angular displacement sampled: "
                    f"**{MAX_CONE_ANGLE_DEG:.1f}°**"
                ),
                (
                    f"- Required minimum distance/vdW-sum ratio: "
                    f"**{TARGET_MIN_RATIO:.2f}**"
                ),
                "- X–H distances: **FIXED**",
                "",
                "## Authorization state",
                "",
                "- Targeted cap repair: **COMPLETED**",
                "- Original fragment geometry optimization: **NOT RUN**",
                "- QM calculation executed: **NO**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (
        F06_DIR / "QM_F06_artificial_cap_repair_summary.json"
    ).write_text(
        json.dumps(
            {
                "decision": decision,
                "all_cap_repair_gates_pass": (
                    all_cap_gates_pass
                ),
                "original_r2_atoms_moved": 0,
                "qm_calculation_executed": False,
                "required_next_step": (
                    "REPEAT_TOPOLOGY_AWARE_AUDIT_ON_REPAIRED_GEOMETRIES"
                    if all_cap_gates_pass
                    else
                    "EXPAND_FRAGMENT_OR_REDESIGN_BOUNDARY"
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("Targeted QM_F06 artificial-cap repair completed.")
    print(f"Decision: {decision}")
    print(f"Caps repaired: {len(repair_rows)}")


if __name__ == "__main__":
    main()
