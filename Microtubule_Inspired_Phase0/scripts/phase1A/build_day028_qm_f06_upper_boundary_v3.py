#!/usr/bin/env python3
"""
Construct QM_F06 UPPER Boundary V3 from Boundary V2.

Correction:
- remove two redundant artificial caps directed toward the same omitted
  real nitrogen center;
- restore A:UPPER:10:4 and A:UPPER:8:4;
- terminate the two new peripheral B-N cuts E:2878 and E:2879 with
  artificial hydrogen caps placed along the original cut-bond vectors.

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

V2_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V2"
)

SOURCE_ATOMS = V2_DIR / "QM_F06_UPPER_BOUNDARY_V2_atoms.csv"
SOURCE_XYZ = V2_DIR / "QM_F06_UPPER_BOUNDARY_V2.xyz"

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
    "QM_F06_UPPER_BOUNDARY_V3"
)

OUTPUT_ATOMS = OUTPUT_DIR / "QM_F06_UPPER_BOUNDARY_V3_atoms.csv"
OUTPUT_CAPS = OUTPUT_DIR / "QM_F06_UPPER_BOUNDARY_V3_caps.csv"
OUTPUT_XYZ = OUTPUT_DIR / "QM_F06_UPPER_BOUNDARY_V3.xyz"

CAPS_TO_REMOVE = {
    "HCAP:UPPER:05",
    "HCAPV2:UPPER:02",
}

REAL_ADDITIONS = (
    "A:UPPER:10:4",
    "A:UPPER:8:4",
)

NEW_CAP_EDGE_IDS = (
    "E:2878",
    "E:2879",
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
        raise RuntimeError(f"No records in {path}")

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

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def parse_xyz(
    path: Path,
) -> list[tuple[str, float, float, float]]:
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


def add(first, second):
    return tuple(
        a + b
        for a, b in zip(first, second, strict=True)
    )


def subtract(first, second):
    return tuple(
        a - b
        for a, b in zip(first, second, strict=True)
    )


def scale(value, factor):
    return tuple(component * factor for component in value)


def norm(value) -> float:
    return math.sqrt(
        sum(component * component for component in value)
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    source_atoms = read_csv(SOURCE_ATOMS)
    source_xyz = parse_xyz(SOURCE_XYZ)
    node_rows = read_csv(FULL_NODES)
    edge_rows = read_csv(FULL_EDGES)
    coordinate_rows = read_csv(FULL_COORDINATES)

    if len(source_atoms) != 28 or len(source_xyz) != 28:
        raise RuntimeError(
            "Expected the 28-atom UPPER Boundary V2 source."
        )

    node_map = {
        row["node_id"]: row
        for row in node_rows
    }

    edge_map = {
        row["edge_id"]: row
        for row in edge_rows
    }

    coordinate_map = {
        row["node_id"]: row
        for row in coordinate_rows
    }

    source_coordinates: dict[
        str,
        tuple[float, float, float],
    ] = {}

    for atom, xyz in zip(
        source_atoms,
        source_xyz,
        strict=True,
    ):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Element mismatch for {atom['atom_id']}"
            )

        source_coordinates[atom["atom_id"]] = xyz[1:]

    source_ids = set(source_coordinates)

    missing_caps = CAPS_TO_REMOVE - source_ids

    if missing_caps:
        raise RuntimeError(
            f"Caps selected for removal are absent: {missing_caps}"
        )

    for atom_id in REAL_ADDITIONS:
        if atom_id in source_ids:
            raise RuntimeError(
                f"Real addition already present: {atom_id}"
            )

        if atom_id not in node_map:
            raise RuntimeError(
                f"Missing topology for {atom_id}"
            )

        if atom_id not in coordinate_map:
            raise RuntimeError(
                f"Missing validated coordinates for {atom_id}"
            )

    retained_rows: list[dict[str, Any]] = []
    final_coordinates: dict[
        str,
        tuple[float, float, float],
    ] = {}

    for atom in source_atoms:
        atom_id = atom["atom_id"]

        if atom_id in CAPS_TO_REMOVE:
            continue

        row = dict(atom)
        row["boundary_v3_source"] = (
            "RETAINED_UPPER_BOUNDARY_V2_ATOM"
        )

        retained_rows.append(row)
        final_coordinates[atom_id] = source_coordinates[atom_id]

    # Determine translation from validated Day024 global coordinates
    # into the local Boundary V2 coordinate frame.
    translation_samples = []

    for atom in source_atoms:
        atom_id = atom["atom_id"]

        if atom["artificial_cap"].lower() == "true":
            continue

        if atom_id not in coordinate_map:
            continue

        global_coord = (
            10.0 * float(coordinate_map[atom_id]["x_nm"]),
            10.0 * float(coordinate_map[atom_id]["y_nm"]),
            10.0 * float(coordinate_map[atom_id]["z_nm"]),
        )

        translation_samples.append(
            subtract(
                source_coordinates[atom_id],
                global_coord,
            )
        )

    if not translation_samples:
        raise RuntimeError(
            "Could not determine the V2 coordinate translation."
        )

    translation = tuple(
        sum(sample[index] for sample in translation_samples)
        / len(translation_samples)
        for index in range(3)
    )

    translation_residuals = [
        norm(
            subtract(sample, translation)
        )
        for sample in translation_samples
    ]

    maximum_translation_residual = max(
        translation_residuals
    )

    if maximum_translation_residual > 1.0e-5:
        raise RuntimeError(
            "Boundary V2 and Day024 coordinates are not related "
            "by a uniform translation. "
            f"Maximum residual: "
            f"{maximum_translation_residual:.6e} Å"
        )

    for atom_id in REAL_ADDITIONS:
        node = node_map[atom_id]
        coordinate = coordinate_map[atom_id]

        global_coord = (
            10.0 * float(coordinate["x_nm"]),
            10.0 * float(coordinate["y_nm"]),
            10.0 * float(coordinate["z_nm"]),
        )

        local_coord = add(
            global_coord,
            translation,
        )

        retained_rows.append(
            {
                "atom_id": atom_id,
                "element": node["element"],
                "atom_role": (
                    "REAL_R2_BOUNDARY_V3_EXPANSION_ATOM"
                ),
                "node_type": node["node_type"],
                "source_edge_id": "",
                "parent_inside_node": "",
                "x_angstrom": f"{local_coord[0]:.10f}",
                "y_angstrom": f"{local_coord[1]:.10f}",
                "z_angstrom": f"{local_coord[2]:.10f}",
                "artificial_cap": False,
                "boundary_v2_source": "",
                "boundary_v3_source": (
                    "RESTORED_REAL_R2_COORDINATION"
                ),
            }
        )

        final_coordinates[atom_id] = local_coord

    cap_rows: list[dict[str, Any]] = []

    for number, edge_id in enumerate(
        NEW_CAP_EDGE_IDS,
        start=1,
    ):
        if edge_id not in edge_map:
            raise RuntimeError(
                f"Missing graph edge {edge_id}"
            )

        edge = edge_map[edge_id]
        source = edge["source_node"]
        target = edge["target_node"]

        source_inside = source in final_coordinates
        target_inside = target in final_coordinates

        if source_inside == target_inside:
            raise RuntimeError(
                f"Expected one inside and one outside atom for "
                f"{edge_id}: {source}, {target}"
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

        if outside not in coordinate_map:
            raise RuntimeError(
                f"Missing coordinate for external atom {outside}"
            )

        outside_global = (
            10.0 * float(coordinate_map[outside]["x_nm"]),
            10.0 * float(coordinate_map[outside]["y_nm"]),
            10.0 * float(coordinate_map[outside]["z_nm"]),
        )

        outside_local = add(
            outside_global,
            translation,
        )

        inside_coord = final_coordinates[inside]
        cut_vector = subtract(
            outside_local,
            inside_coord,
        )
        cut_length = norm(cut_vector)

        if cut_length <= 1.0e-12:
            raise RuntimeError(
                f"Zero cut vector for {edge_id}"
            )

        target_distance = XH_TARGETS[inside_element]

        cap_coord = add(
            inside_coord,
            scale(
                cut_vector,
                target_distance / cut_length,
            ),
        )

        cap_id = f"HCAPV3:UPPER:{number:02d}"

        retained_rows.append(
            {
                "atom_id": cap_id,
                "element": "H",
                "atom_role": "ARTIFICIAL_BOUNDARY_CAP_V3",
                "node_type": "QM_BOUNDARY_CAP_H",
                "source_edge_id": edge_id,
                "parent_inside_node": inside,
                "x_angstrom": f"{cap_coord[0]:.10f}",
                "y_angstrom": f"{cap_coord[1]:.10f}",
                "z_angstrom": f"{cap_coord[2]:.10f}",
                "artificial_cap": True,
                "boundary_v2_source": "",
                "boundary_v3_source": (
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
                "target_XH_distance_angstrom": (
                    target_distance
                ),
                "placement_method": (
                    "ALONG_ORIGINAL_INSIDE_TO_OUTSIDE_BN_VECTOR"
                ),
            }
        )

    if len(retained_rows) != 30:
        raise RuntimeError(
            f"Expected 30 final atoms; found {len(retained_rows)}"
        )

    final_ids = [
        row["atom_id"]
        for row in retained_rows
    ]

    if len(final_ids) != len(set(final_ids)):
        raise RuntimeError(
            "Duplicate atom IDs in Boundary V3."
        )

    for row in retained_rows:
        atom_id = row["atom_id"]
        x, y, z = final_coordinates[atom_id]

        row["x_angstrom"] = f"{x:.10f}"
        row["y_angstrom"] = f"{y:.10f}"
        row["z_angstrom"] = f"{z:.10f}"

    write_csv(
        OUTPUT_ATOMS,
        retained_rows,
    )

    write_csv(
        OUTPUT_CAPS,
        cap_rows,
    )

    xyz_lines = [
        str(len(retained_rows)),
        (
            "QM_F06_UPPER_BOUNDARY_V3; shared omitted "
            "A:UPPER:10:4/A:UPPER:8:4 coordination restored; "
            "unoptimized; no QM executed"
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

    neutral_valence_electrons = (
        3 * element_counts["B"]
        + 5 * element_counts["N"]
        + element_counts["H"]
    )

    summary = {
        "decision": (
            "QM_F06_UPPER_BOUNDARY_V3_CONSTRUCTED_"
            "PRE_QM_AUDIT_REQUIRED"
        ),
        "source_fragment": "QM_F06_UPPER_BOUNDARY_V2",
        "source_atom_count": 28,
        "removed_caps": sorted(CAPS_TO_REMOVE),
        "real_atoms_added": list(REAL_ADDITIONS),
        "new_cap_edge_ids": list(NEW_CAP_EDGE_IDS),
        "new_artificial_caps_added": 2,
        "final_atom_count": len(retained_rows),
        "element_counts": dict(element_counts),
        "role_counts": dict(role_counts),
        "neutral_valence_electrons": (
            neutral_valence_electrons
        ),
        "provisional_charge": 0,
        "provisional_multiplicity": (
            1
            if neutral_valence_electrons % 2 == 0
            else 2
        ),
        "translation_angstrom": translation,
        "translation_fit_max_residual_angstrom": (
            maximum_translation_residual
        ),
        "geometry_optimized": False,
        "qm_calculation_executed": False,
        "pre_qm_audit_required": True,
        "orca_input_preparation_authorized": False,
        "orca_execution_authorized": False,
    }

    (
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V3_build_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V3_BUILD_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER Boundary V3 Build — Day028",
                "",
                "## Removed redundant caps",
                "",
                *[
                    f"- `{cap_id}`"
                    for cap_id in sorted(CAPS_TO_REMOVE)
                ],
                "",
                "## Restored real R2 atoms",
                "",
                *[
                    f"- `{atom_id}`"
                    for atom_id in REAL_ADDITIONS
                ],
                "",
                "## New peripheral caps",
                "",
                *[
                    (
                        f"- `{row['cap_id']}` on "
                        f"`{row['source_edge_id']}`: "
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
                (
                    "- Neutral valence-electron count: "
                    f"**{neutral_valence_electrons}**"
                ),
                "- Provisional charge/multiplicity: **0 / 1**",
                "- Geometry optimized: **NO**",
                "- QM calculation executed: **NO**",
                "",
                "## Decision",
                "",
                (
                    "**QM_F06_UPPER_BOUNDARY_V3_CONSTRUCTED_"
                    "PRE_QM_AUDIT_REQUIRED**"
                ),
                "",
                "## Authorization state",
                "",
                "- Boundary V3 construction: **COMPLETED**",
                "- Pre-QM audit: **PENDING**",
                "- ORCA input preparation: **NOT AUTHORIZED**",
                "- ORCA execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER Boundary V3 construction completed.")
    print("Final atoms:", len(retained_rows))
    print("Element counts:", dict(element_counts))
    print(
        "Neutral valence electrons:",
        neutral_valence_electrons,
    )
    print("Real atoms added:", len(REAL_ADDITIONS))
    print("Caps removed:", len(CAPS_TO_REMOVE))
    print("New artificial caps:", len(cap_rows))
    print("QM executed: False")
    print("Output:", OUTPUT_XYZ)


if __name__ == "__main__":
    main()
