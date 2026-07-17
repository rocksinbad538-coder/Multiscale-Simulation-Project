#!/usr/bin/env python3
"""
Construct QM_F06 UPPER Boundary V2.

Changes relative to the repaired 22-atom UPPER fragment:

- remove HCAP:UPPER:07;
- restore three real heavy atoms:
    A:UPPER:11:3
    A:UPPER:13:3
    A:UPPER:14:2
- restore real R2 hydrogen:
    H4:UPPER:0046:0
- add three new artificial caps for cut edges:
    E:2897
    E:2898
    E:2913

No geometry optimization or QM calculation is executed.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

SOURCE_ATOMS = F06_DIR / (
    "QM_F06_UPPER_CAPPED_REPAIRED_atoms.csv"
)

SOURCE_XYZ = F06_DIR / (
    "QM_F06_UPPER_CAPPED_REPAIRED.xyz"
)

FULL_COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

FULL_NODES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

FULL_EDGES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V2"
)

OUTPUT_ATOMS = OUTPUT_DIR / (
    "QM_F06_UPPER_BOUNDARY_V2_atoms.csv"
)

OUTPUT_CAPS = OUTPUT_DIR / (
    "QM_F06_UPPER_BOUNDARY_V2_caps.csv"
)

OUTPUT_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_BOUNDARY_V2.xyz"
)

REMOVE_CAP = "HCAP:UPPER:07"

REAL_ADDITIONS = (
    "A:UPPER:11:3",
    "A:UPPER:13:3",
    "A:UPPER:14:2",
    "H4:UPPER:0046:0",
)

NEW_CAP_EDGE_IDS = (
    "E:2897",
    "E:2898",
    "E:2913",
)

XH_TARGETS = {
    "B": 1.19,
    "N": 1.01,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    return rows


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(f"Cannot write empty table: {path}")

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

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


def vector(first, second):
    return tuple(
        second_value - first_value
        for first_value, second_value
        in zip(first, second, strict=True)
    )


def norm(value) -> float:
    return math.sqrt(sum(component**2 for component in value))


def add(first, second):
    return tuple(
        first_value + second_value
        for first_value, second_value
        in zip(first, second, strict=True)
    )


def scale(value, factor):
    return tuple(component * factor for component in value)


def parse_xyz(path: Path):
    require_file(path)

    lines = path.read_text(encoding="utf-8").splitlines()
    expected = int(lines[0].strip())

    rows = [
        (
            fields[0],
            float(fields[1]),
            float(fields[2]),
            float(fields[3]),
        )
        for line in lines[2:]
        if line.strip()
        for fields in [line.split()]
    ]

    if len(rows) != expected:
        raise RuntimeError(
            f"XYZ atom-count mismatch: {expected} vs {len(rows)}"
        )

    return rows


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    source_atoms = read_csv(SOURCE_ATOMS)
    source_xyz = parse_xyz(SOURCE_XYZ)
    coordinate_rows = read_csv(FULL_COORDINATES)
    node_rows = read_csv(FULL_NODES)
    edge_rows = read_csv(FULL_EDGES)

    if len(source_atoms) != 22 or len(source_xyz) != 22:
        raise RuntimeError("Expected repaired 22-atom UPPER source.")

    node_map = {
        row["node_id"]: row
        for row in node_rows
    }

    coordinate_map = {
        row["node_id"]: row
        for row in coordinate_rows
    }

    edge_map = {
        row["edge_id"]: row
        for row in edge_rows
    }

    source_ids = [row["atom_id"] for row in source_atoms]

    if REMOVE_CAP not in source_ids:
        raise RuntimeError(f"Cap to remove is absent: {REMOVE_CAP}")

    source_coordinate_map = {}

    for atom, xyz in zip(
        source_atoms,
        source_xyz,
        strict=True,
    ):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Source element mismatch: {atom['atom_id']}"
            )

        source_coordinate_map[atom["atom_id"]] = xyz[1:]

    retained_rows = []
    final_coordinates = {}

    for atom in source_atoms:
        atom_id = atom["atom_id"]

        if atom_id == REMOVE_CAP:
            continue

        retained = dict(atom)
        retained["boundary_v2_source"] = (
            "RETAINED_UPPER_REPAIRED_FRAGMENT_ATOM"
        )

        retained_rows.append(retained)
        final_coordinates[atom_id] = source_coordinate_map[atom_id]

    # The original repaired fragment coordinates may be translated
    # relative to the Day024 global frame. Determine the translation
    # using real atoms common to both data sets.
    translation_vectors = []

    for atom in source_atoms:
        atom_id = atom["atom_id"]

        if atom.get("artificial_cap", "").lower() == "true":
            continue

        if atom_id not in coordinate_map:
            continue

        global_coord = (
            10.0 * float(coordinate_map[atom_id]["x_nm"]),
            10.0 * float(coordinate_map[atom_id]["y_nm"]),
            10.0 * float(coordinate_map[atom_id]["z_nm"]),
        )

        local_coord = source_coordinate_map[atom_id]

        translation_vectors.append(
            tuple(
                local - global_value
                for local, global_value
                in zip(local_coord, global_coord, strict=True)
            )
        )

    if not translation_vectors:
        raise RuntimeError("Unable to determine coordinate translation.")

    translation = tuple(
        sum(vector_value[index] for vector_value in translation_vectors)
        / len(translation_vectors)
        for index in range(3)
    )

    max_translation_residual = 0.0

    for vector_value in translation_vectors:
        residual = norm(
            tuple(
                vector_value[index] - translation[index]
                for index in range(3)
            )
        )
        max_translation_residual = max(
            max_translation_residual,
            residual,
        )

    if max_translation_residual > 1.0e-5:
        raise RuntimeError(
            "Source/global coordinates are not related by a uniform "
            f"translation; max residual={max_translation_residual:.6e} Å"
        )

    for atom_id in REAL_ADDITIONS:
        if atom_id in final_coordinates:
            raise RuntimeError(
                f"Real addition already present: {atom_id}"
            )

        if atom_id not in node_map or atom_id not in coordinate_map:
            raise RuntimeError(
                f"Missing topology or coordinates for {atom_id}"
            )

        node = node_map[atom_id]
        coordinate = coordinate_map[atom_id]

        global_coord = (
            10.0 * float(coordinate["x_nm"]),
            10.0 * float(coordinate["y_nm"]),
            10.0 * float(coordinate["z_nm"]),
        )

        local_coord = add(global_coord, translation)

        role = (
            "REAL_R2_BOUNDARY_EXPANSION_HYDROGEN"
            if node["element"] == "H"
            else "REAL_R2_BOUNDARY_EXPANSION_ATOM"
        )

        retained_rows.append(
            {
                "atom_id": atom_id,
                "element": node["element"],
                "atom_role": role,
                "node_type": node["node_type"],
                "source_edge_id": "",
                "parent_inside_node": "",
                "x_angstrom": f"{local_coord[0]:.10f}",
                "y_angstrom": f"{local_coord[1]:.10f}",
                "z_angstrom": f"{local_coord[2]:.10f}",
                "artificial_cap": False,
                "boundary_v2_source": (
                    "RESTORED_REAL_R2_COORDINATION"
                ),
            }
        )

        final_coordinates[atom_id] = local_coord

    cap_rows = []

    for cap_number, edge_id in enumerate(
        NEW_CAP_EDGE_IDS,
        start=1,
    ):
        if edge_id not in edge_map:
            raise RuntimeError(f"Missing cut edge: {edge_id}")

        edge = edge_map[edge_id]
        source = edge["source_node"]
        target = edge["target_node"]

        source_inside = source in final_coordinates
        target_inside = target in final_coordinates

        if source_inside == target_inside:
            raise RuntimeError(
                f"Expected exactly one inside atom for {edge_id}: "
                f"{source}, {target}"
            )

        inside = source if source_inside else target
        outside = target if source_inside else source

        inside_element = node_map[inside]["element"]
        outside_element = node_map[outside]["element"]

        if {
            inside_element,
            outside_element,
        } != {"B", "N"}:
            raise RuntimeError(
                f"Expected B-N cut for {edge_id}; observed "
                f"{inside_element}-{outside_element}"
            )

        outside_global = (
            10.0 * float(coordinate_map[outside]["x_nm"]),
            10.0 * float(coordinate_map[outside]["y_nm"]),
            10.0 * float(coordinate_map[outside]["z_nm"]),
        )

        outside_local = add(outside_global, translation)
        inside_coord = final_coordinates[inside]

        direction = vector(inside_coord, outside_local)
        direction_norm = norm(direction)

        if direction_norm <= 1.0e-12:
            raise RuntimeError(
                f"Zero cut-vector length for {edge_id}"
            )

        target_distance = XH_TARGETS[inside_element]
        unit_direction = scale(
            direction,
            1.0 / direction_norm,
        )

        cap_coord = add(
            inside_coord,
            scale(unit_direction, target_distance),
        )

        cap_id = f"HCAPV2:UPPER:{cap_number:02d}"

        retained_rows.append(
            {
                "atom_id": cap_id,
                "element": "H",
                "atom_role": "ARTIFICIAL_BOUNDARY_CAP_V2",
                "node_type": "QM_BOUNDARY_CAP_H",
                "source_edge_id": edge_id,
                "parent_inside_node": inside,
                "x_angstrom": f"{cap_coord[0]:.10f}",
                "y_angstrom": f"{cap_coord[1]:.10f}",
                "z_angstrom": f"{cap_coord[2]:.10f}",
                "artificial_cap": True,
                "boundary_v2_source": (
                    "NEW_CAP_ALONG_ORIGINAL_BN_CUT_VECTOR"
                ),
            }
        )

        final_coordinates[cap_id] = cap_coord

        cap_rows.append(
            {
                "cap_id": cap_id,
                "parent_inside_node": inside,
                "parent_element": inside_element,
                "removed_outside_node": outside,
                "removed_outside_element": outside_element,
                "source_edge_id": edge_id,
                "target_XH_distance_angstrom": target_distance,
                "placement_method": (
                    "ALONG_ORIGINAL_INSIDE_TO_OUTSIDE_BN_VECTOR"
                ),
            }
        )

    if len(retained_rows) != 28:
        raise RuntimeError(
            f"Expected 28 final atoms; found {len(retained_rows)}"
        )

    atom_ids = [row["atom_id"] for row in retained_rows]

    if len(atom_ids) != len(set(atom_ids)):
        raise RuntimeError("Duplicate atom IDs in final fragment.")

    for row in retained_rows:
        atom_id = row["atom_id"]
        coord = final_coordinates[atom_id]

        row["x_angstrom"] = f"{coord[0]:.10f}"
        row["y_angstrom"] = f"{coord[1]:.10f}"
        row["z_angstrom"] = f"{coord[2]:.10f}"

    write_csv(OUTPUT_ATOMS, retained_rows)
    write_csv(OUTPUT_CAPS, cap_rows)

    xyz_lines = [
        str(len(retained_rows)),
        (
            "QM_F06_UPPER_BOUNDARY_V2; HCAP:UPPER:07 replaced "
            "by real R2 coordination; unoptimized; no QM executed"
        ),
    ]

    for row in retained_rows:
        atom_id = row["atom_id"]
        x, y, z = final_coordinates[atom_id]

        xyz_lines.append(
            f"{row['element']:<2s} "
            f"{x: .10f} {y: .10f} {z: .10f}"
        )

    OUTPUT_XYZ.write_text(
        "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )

    element_counts = Counter(
        row["element"]
        for row in retained_rows
    )

    role_counts = Counter(
        row["atom_role"]
        for row in retained_rows
    )

    summary = {
        "decision": (
            "QM_F06_UPPER_BOUNDARY_V2_CONSTRUCTED_"
            "PRE_QM_AUDIT_REQUIRED"
        ),
        "source_atoms": 22,
        "removed_cap": REMOVE_CAP,
        "real_heavy_atoms_added": 3,
        "real_hydrogens_added": 1,
        "new_artificial_caps_added": 3,
        "final_atom_count": len(retained_rows),
        "element_counts": dict(element_counts),
        "role_counts": dict(role_counts),
        "translation_angstrom": translation,
        "translation_fit_max_residual_angstrom": (
            max_translation_residual
        ),
        "geometry_optimized": False,
        "qm_calculation_executed": False,
        "pre_qm_audit_required": True,
        "orca_input_preparation_authorized": False,
        "orca_execution_authorized": False,
    }

    (
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_build_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_BUILD_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER Boundary V2 Build — Day028",
                "",
                "## Boundary correction",
                "",
                f"- Removed cap: `{REMOVE_CAP}`",
                "",
                "## Real R2 atoms restored",
                "",
                *[f"- `{atom_id}`" for atom_id in REAL_ADDITIONS],
                "",
                "## New peripheral caps",
                "",
                *[
                    (
                        f"- `{row['cap_id']}` on `{row['source_edge_id']}`: "
                        f"`{row['parent_inside_node']} — "
                        f"{row['removed_outside_node']}`"
                    )
                    for row in cap_rows
                ],
                "",
                "## Result",
                "",
                f"- Final atoms: **{len(retained_rows)}**",
                f"- Element counts: `{dict(element_counts)}`",
                f"- Role counts: `{dict(role_counts)}`",
                "- Geometry optimized: **NO**",
                "- QM calculation executed: **NO**",
                "",
                "## Decision",
                "",
                (
                    "**QM_F06_UPPER_BOUNDARY_V2_CONSTRUCTED_"
                    "PRE_QM_AUDIT_REQUIRED**"
                ),
                "",
                "## Authorization state",
                "",
                "- Boundary construction: **COMPLETED**",
                "- Pre-QM audit: **PENDING**",
                "- ORCA input preparation: **NOT AUTHORIZED**",
                "- ORCA execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER Boundary V2 construction completed.")
    print("Final atoms:", len(retained_rows))
    print("Element counts:", dict(element_counts))
    print("Real heavy atoms added: 3")
    print("Real hydrogens added: 1")
    print("New artificial caps: 3")
    print("QM executed: False")
    print("Output:", OUTPUT_XYZ)


if __name__ == "__main__":
    main()
