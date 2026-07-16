#!/usr/bin/env python3
"""
Targeted angular repair of HCAPV2:LOWER:02.

Only the artificial cap HCAPV2:LOWER:02 is repositioned.

Constraints:
- parent atom A:LOWER:11:-3 remains fixed;
- N-H distance remains exactly 1.01 Å;
- all real R2 atoms remain fixed;
- all other caps remain fixed;
- candidate orientations are screened against topology-aware
  nonbonded contacts;
- no QM calculation is executed.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict, deque
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

V2_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2"
)

ATOMS_PATH = V2_DIR / "QM_F06_LOWER_BOUNDARY_V2_atoms.csv"

FULL_EDGES_PATH = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = V2_DIR / "cap02_repair"

TARGET_CAP = "HCAPV2:LOWER:02"
PARENT_ATOM = "A:LOWER:11:-3"
TARGET_BOND_LENGTH = 1.01

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

TARGET_MINIMUM_RATIO = 0.80


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows in {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []

    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def xyz(row: dict[str, str]) -> tuple[float, float, float]:
    return (
        float(row["x_angstrom"]),
        float(row["y_angstrom"]),
        float(row["z_angstrom"]),
    )


def add(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        a + b
        for a, b in zip(first, second, strict=True)
    )


def subtract(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        a - b
        for a, b in zip(first, second, strict=True)
    )


def scale(
    vector: tuple[float, float, float],
    factor: float,
) -> tuple[float, float, float]:
    return tuple(value * factor for value in vector)


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


def norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(dot(vector, vector))


def normalize(
    vector: tuple[float, float, float],
) -> tuple[float, float, float]:
    length = norm(vector)

    if length <= 1.0e-14:
        raise RuntimeError("Cannot normalize zero vector.")

    return scale(vector, 1.0 / length)


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return norm(subtract(first, second))


def shortest_path_length(
    adjacency: dict[str, set[str]],
    source: str,
    target: str,
) -> int | None:
    if source == target:
        return 0

    visited = {source}
    queue = deque([(source, 0)])

    while queue:
        node, depth = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor == target:
                return depth + 1

            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return None


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = read_csv(ATOMS_PATH)
    atoms = {row["atom_id"]: row for row in rows}

    if TARGET_CAP not in atoms:
        raise RuntimeError(f"Target cap not found: {TARGET_CAP}")

    if PARENT_ATOM not in atoms:
        raise RuntimeError(f"Parent not found: {PARENT_ATOM}")

    atom_ids = list(atoms)
    atom_set = set(atom_ids)

    coordinates = {
        atom_id: xyz(row)
        for atom_id, row in atoms.items()
    }

    elements = {
        atom_id: row["element"]
        for atom_id, row in atoms.items()
    }

    adjacency: dict[str, set[str]] = defaultdict(set)

    for row in read_csv(FULL_EDGES_PATH):
        first = row["source_node"]
        second = row["target_node"]

        if first in atom_set and second in atom_set:
            adjacency[first].add(second)
            adjacency[second].add(first)

    for atom_id, row in atoms.items():
        if row["artificial_cap"].lower() == "true":
            parent = row["parent_inside_node"]

            if not parent:
                raise RuntimeError(
                    f"Artificial cap has no parent: {atom_id}"
                )

            adjacency[atom_id].add(parent)
            adjacency[parent].add(atom_id)

    parent_xyz = coordinates[PARENT_ATOM]
    old_cap_xyz = coordinates[TARGET_CAP]

    original_direction = normalize(
        subtract(old_cap_xyz, parent_xyz)
    )

    # Construct a stable orthonormal basis around the original N-H vector.
    reference = (0.0, 0.0, 1.0)

    if abs(dot(original_direction, reference)) > 0.90:
        reference = (0.0, 1.0, 0.0)

    basis_u = normalize(cross(original_direction, reference))
    basis_v = normalize(cross(original_direction, basis_u))

    best: dict[str, Any] | None = None

    for cone_angle_deg in range(5, 176, 5):
        theta = math.radians(cone_angle_deg)

        for azimuth_deg in range(0, 360, 5):
            phi = math.radians(azimuth_deg)

            candidate_direction = add(
                scale(original_direction, math.cos(theta)),
                add(
                    scale(
                        basis_u,
                        math.sin(theta) * math.cos(phi),
                    ),
                    scale(
                        basis_v,
                        math.sin(theta) * math.sin(phi),
                    ),
                ),
            )

            candidate_direction = normalize(candidate_direction)

            candidate_xyz = add(
                parent_xyz,
                scale(candidate_direction, TARGET_BOND_LENGTH),
            )

            minimum_ratio = float("inf")
            minimum_partner = ""

            hard_failure = False

            for other_id in atom_ids:
                if other_id in {TARGET_CAP, PARENT_ATOM}:
                    continue

                graph_separation = shortest_path_length(
                    adjacency,
                    TARGET_CAP,
                    other_id,
                )

                # Exclude 1-2, 1-3 and 1-4 intramolecular pairs.
                if graph_separation in {1, 2, 3}:
                    continue

                measured = distance(
                    candidate_xyz,
                    coordinates[other_id],
                )

                ratio = measured / (
                    VDW_RADII["H"]
                    + VDW_RADII[elements[other_id]]
                )

                if ratio < minimum_ratio:
                    minimum_ratio = ratio
                    minimum_partner = other_id

                if ratio < TARGET_MINIMUM_RATIO:
                    hard_failure = True

            candidate = {
                "cone_angle_deg": cone_angle_deg,
                "azimuth_deg": azimuth_deg,
                "candidate_xyz": candidate_xyz,
                "minimum_ratio": minimum_ratio,
                "minimum_partner": minimum_partner,
                "gate_pass": not hard_failure,
            }

            if best is None:
                best = candidate
            elif candidate["gate_pass"] and not best["gate_pass"]:
                best = candidate
            elif (
                candidate["gate_pass"] == best["gate_pass"]
                and candidate["minimum_ratio"]
                > best["minimum_ratio"]
            ):
                best = candidate

    if best is None:
        raise RuntimeError("No candidate orientations were generated.")

    if not best["gate_pass"]:
        raise RuntimeError(
            "No cap orientation reached the required minimum "
            f"vdW ratio of {TARGET_MINIMUM_RATIO:.2f}. "
            f"Best ratio: {best['minimum_ratio']:.6f}"
        )

    repaired_rows: list[dict[str, Any]] = []

    for row in rows:
        updated = dict(row)

        if row["atom_id"] == TARGET_CAP:
            new_xyz = best["candidate_xyz"]

            updated["x_angstrom"] = f"{new_xyz[0]:.10f}"
            updated["y_angstrom"] = f"{new_xyz[1]:.10f}"
            updated["z_angstrom"] = f"{new_xyz[2]:.10f}"
            updated["coordinate_source"] = (
                "TARGETED_CAP02_ANGULAR_REPAIR"
            )

        repaired_rows.append(updated)

    repaired_atoms_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_REPAIRED_atoms.csv"
    )
    write_csv(repaired_atoms_path, repaired_rows)

    repaired_xyz_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_REPAIRED.xyz"
    )

    xyz_lines = [
        str(len(repaired_rows)),
        (
            "QM_F06_LOWER_BOUNDARY_V2_REPAIRED; "
            "only HCAPV2:LOWER:02 angularly repositioned; "
            "all real atoms fixed; unoptimized"
        ),
    ]

    for row in repaired_rows:
        xyz_lines.append(
            f"{row['element']:<2s} "
            f"{float(row['x_angstrom']): .10f} "
            f"{float(row['y_angstrom']): .10f} "
            f"{float(row['z_angstrom']): .10f}"
        )

    repaired_xyz_path.write_text(
        "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )

    new_cap_xyz = best["candidate_xyz"]
    measured_bond = distance(parent_xyz, new_cap_xyz)

    manifest = [
        {
            "cap_id": TARGET_CAP,
            "parent_atom": PARENT_ATOM,
            "old_x_angstrom": f"{old_cap_xyz[0]:.10f}",
            "old_y_angstrom": f"{old_cap_xyz[1]:.10f}",
            "old_z_angstrom": f"{old_cap_xyz[2]:.10f}",
            "new_x_angstrom": f"{new_cap_xyz[0]:.10f}",
            "new_y_angstrom": f"{new_cap_xyz[1]:.10f}",
            "new_z_angstrom": f"{new_cap_xyz[2]:.10f}",
            "bond_length_angstrom": f"{measured_bond:.10f}",
            "angular_change_deg": best["cone_angle_deg"],
            "azimuth_deg": best["azimuth_deg"],
            "minimum_post_repair_vdw_ratio": (
                f"{best['minimum_ratio']:.10f}"
            ),
            "minimum_post_repair_partner": (
                best["minimum_partner"]
            ),
            "target_minimum_ratio": TARGET_MINIMUM_RATIO,
            "repair_gate_pass": best["gate_pass"],
        }
    ]

    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_cap02_repair_manifest.csv",
        manifest,
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_CAP02_REPAIR_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2 Cap02 Repair — Day027",
                "",
                "## Scope",
                "",
                "- Repaired atom: `HCAPV2:LOWER:02`",
                "- Parent: `A:LOWER:11:-3`",
                "- Real R2 atoms moved: **0**",
                "- Other artificial caps moved: **0**",
                "- N–H bond length: **FIXED at 1.01 Å**",
                "",
                "## Search result",
                "",
                (
                    "- Angular displacement: "
                    f"**{best['cone_angle_deg']}°**"
                ),
                (
                    "- Azimuth: "
                    f"**{best['azimuth_deg']}°**"
                ),
                (
                    "- Minimum post-repair vdW ratio: "
                    f"**{best['minimum_ratio']:.6f}**"
                ),
                (
                    "- Closest screened partner: "
                    f"`{best['minimum_partner']}`"
                ),
                "",
                "## Decision",
                "",
                "**QM_F06_LOWER_BOUNDARY_V2_CAP02_REPAIR_COMPLETED**",
                "",
                "- QM calculation executed: **NO**",
                "- Pre-QM gate repeated: **NO — required next**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": (
            "QM_F06_LOWER_BOUNDARY_V2_CAP02_REPAIR_COMPLETED"
        ),
        "cap_repaired": TARGET_CAP,
        "real_atoms_moved": 0,
        "other_caps_moved": 0,
        "bond_length_angstrom": measured_bond,
        "minimum_post_repair_vdw_ratio": (
            best["minimum_ratio"]
        ),
        "repair_gate_pass": best["gate_pass"],
        "qm_executed": False,
        "required_next_step": (
            "REPEAT_BOUNDARY_V2_PRE_QM_AUDIT"
        ),
    }

    (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_cap02_repair_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Boundary V2 Cap02 repair completed.")
    print("Cap:", TARGET_CAP)
    print("Bond length:", f"{measured_bond:.10f} Å")
    print(
        "Minimum vdW ratio:",
        f"{best['minimum_ratio']:.10f}",
    )
    print("Closest partner:", best["minimum_partner"])
    print("Real atoms moved: 0")
    print("QM executed: False")
    print("Output:", repaired_xyz_path)


if __name__ == "__main__":
    main()
