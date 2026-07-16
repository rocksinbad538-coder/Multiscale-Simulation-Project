#!/usr/bin/env python3
"""
Build QM_F06 LOWER boundary-expanded V2 fragment.

Scientific correction
---------------------
Remove artificial cap HCAP:LOWER:07, which migrated toward the bridge
during Stage 2, and restore its real R2 environment:

- A:LOWER:13:-3
- A:LOWER:11:-3
- A:LOWER:14:-2
- H4:LOWER:0016:0

The remaining peripheral B-N cuts are saturated with new H caps along
the original R2 cut-bond vectors.

Existing V1 atoms use the converged Stage-2 geometry.
New real atoms use validated Day024 R2 coordinates.
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

WORKFLOW_STATE = F06_DIR / (
    "orca_workflow/QM_F06_LOWER_CAPPED_REPAIRED/"
    "workflow_state.json"
)

V1_ATOMS = F06_DIR / (
    "QM_F06_LOWER_CAPPED_REPAIRED_atoms.csv"
)

FULL_COORDS = ROOT / (
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

V1_INTERNAL_EDGES = F06_DIR / (
    "QM_F06_LOWER_internal_edges.csv"
)

V1_BOUNDARY_AUDIT = F06_DIR / (
    "QM_F06_LOWER_boundary_edge_audit.csv"
)

V1_CAPS = F06_DIR / (
    "QM_F06_LOWER_CAPPED_caps.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2"
)

REMOVE_ATOMS = {
    "HCAP:LOWER:07",
}

ADD_REAL_ATOMS = {
    "A:LOWER:13:-3",
    "A:LOWER:11:-3",
    "A:LOWER:14:-2",
    "H4:LOWER:0016:0",
}

NEW_CAP_CUTS = (
    ("A:LOWER:11:-3", "A:LOWER:10:-2"),
    ("A:LOWER:11:-3", "A:LOWER:10:-4"),
    ("A:LOWER:14:-2", "A:LOWER:13:-1"),
)

TARGET_XH = {
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
        raise RuntimeError(f"No rows in {path}")

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
                {field: row.get(field, "") for field in fields}
            )


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    require_file(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    count = int(lines[0].strip())
    rows = [
        (
            parts[0],
            float(parts[1]),
            float(parts[2]),
            float(parts[3]),
        )
        for line in lines[2:]
        if line.strip()
        for parts in [line.split()]
    ]

    if len(rows) != count:
        raise RuntimeError(
            f"XYZ count mismatch: header={count}, rows={len(rows)}"
        )

    return rows


def vector(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        b - a
        for a, b in zip(first, second, strict=True)
    )


def norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(sum(value * value for value in v))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    state = json.loads(
        WORKFLOW_STATE.read_text(encoding="utf-8")
    )

    if state.get("stage2_validation_pass") is not True:
        raise RuntimeError("Stage-2 validation did not pass.")

    stage2_xyz_path = ROOT / state["stage2_optimized_xyz"]
    stage2_xyz = read_xyz(stage2_xyz_path)

    v1_rows = read_csv(V1_ATOMS)

    if len(v1_rows) != len(stage2_xyz):
        raise RuntimeError("V1 atom order and Stage-2 XYZ mismatch.")

    stage2_coords: dict[str, tuple[float, float, float]] = {}

    for atom, xyz in zip(v1_rows, stage2_xyz, strict=True):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Element mismatch for {atom['atom_id']}"
            )

        stage2_coords[atom["atom_id"]] = xyz[1:]

    full_coord_rows = read_csv(FULL_COORDS)
    full_coords = {
        row["node_id"]: (
            float(row["x_nm"]) * 10.0,
            float(row["y_nm"]) * 10.0,
            float(row["z_nm"]) * 10.0,
        )
        for row in full_coord_rows
    }

    full_node_rows = read_csv(FULL_NODES)
    full_nodes = {
        row["node_id"]: row
        for row in full_node_rows
    }

    full_edge_rows = read_csv(FULL_EDGES)
    full_edges = {
        tuple(sorted((row["source_node"], row["target_node"]))): row
        for row in full_edge_rows
    }

    output_atoms: list[dict[str, Any]] = []

    # Retain all Stage-2 V1 atoms except invalid cap07.
    for row in v1_rows:
        atom_id = row["atom_id"]

        if atom_id in REMOVE_ATOMS:
            continue

        x, y, z = stage2_coords[atom_id]

        output_atoms.append(
            {
                "atom_id": atom_id,
                "element": row["element"],
                "atom_role": row["atom_role"],
                "node_type": row["node_type"],
                "coordinate_source": "CONVERGED_STAGE2_V1",
                "x_angstrom": f"{x:.10f}",
                "y_angstrom": f"{y:.10f}",
                "z_angstrom": f"{z:.10f}",
                "artificial_cap": row["artificial_cap"],
                "parent_inside_node": row.get(
                    "parent_inside_node",
                    "",
                ),
                "source_edge_id": row.get("source_edge_id", ""),
            }
        )

    # Add real R2 atoms from validated Day024 coordinates.
    for atom_id in sorted(ADD_REAL_ATOMS):
        node = full_nodes[atom_id]
        x, y, z = full_coords[atom_id]

        output_atoms.append(
            {
                "atom_id": atom_id,
                "element": node["element"],
                "atom_role": "REAL_R2_BOUNDARY_EXPANSION_ATOM",
                "node_type": node["node_type"],
                "coordinate_source": "VALIDATED_DAY024_R2_COORDINATE",
                "x_angstrom": f"{x:.10f}",
                "y_angstrom": f"{y:.10f}",
                "z_angstrom": f"{z:.10f}",
                "artificial_cap": False,
                "parent_inside_node": "",
                "source_edge_id": "",
            }
        )

    atom_lookup = {
        row["atom_id"]: row
        for row in output_atoms
    }

    # Add new peripheral caps.
    cap_rows: list[dict[str, Any]] = []

    for index, (inside, outside) in enumerate(
        NEW_CAP_CUTS,
        start=1,
    ):
        if inside not in atom_lookup:
            raise RuntimeError(f"Inside atom absent: {inside}")

        inside_element = full_nodes[inside]["element"]
        target = TARGET_XH[inside_element]

        inside_xyz = tuple(
            float(atom_lookup[inside][key])
            for key in (
                "x_angstrom",
                "y_angstrom",
                "z_angstrom",
            )
        )
        outside_xyz = full_coords[outside]

        direction = vector(inside_xyz, outside_xyz)
        length = norm(direction)

        if length <= 1.0e-12:
            raise RuntimeError(
                f"Invalid cut vector for {inside} — {outside}"
            )

        unit = tuple(component / length for component in direction)
        cap_xyz = tuple(
            coordinate + target * component
            for coordinate, component in zip(
                inside_xyz,
                unit,
                strict=True,
            )
        )

        cap_id = f"HCAPV2:LOWER:{index:02d}"
        edge = full_edges[tuple(sorted((inside, outside)))]

        cap_record = {
            "atom_id": cap_id,
            "element": "H",
            "atom_role": "ARTIFICIAL_BOUNDARY_CAP_V2",
            "node_type": "QM_BOUNDARY_CAP_H",
            "coordinate_source": "ORIGINAL_CUT_BOND_VECTOR",
            "x_angstrom": f"{cap_xyz[0]:.10f}",
            "y_angstrom": f"{cap_xyz[1]:.10f}",
            "z_angstrom": f"{cap_xyz[2]:.10f}",
            "artificial_cap": True,
            "parent_inside_node": inside,
            "source_edge_id": edge["edge_id"],
        }

        output_atoms.append(cap_record)
        atom_lookup[cap_id] = cap_record

        cap_rows.append(
            {
                "cap_id": cap_id,
                "parent_inside_node": inside,
                "parent_element": inside_element,
                "removed_outside_node": outside,
                "removed_outside_element": (
                    full_nodes[outside]["element"]
                ),
                "source_edge_id": edge["edge_id"],
                "target_XH_distance_angstrom": target,
                "placement_method": (
                    "ALONG_ORIGINAL_INSIDE_TO_OUTSIDE_BN_VECTOR"
                ),
            }
        )

    expected_atoms = 28

    if len(output_atoms) != expected_atoms:
        raise RuntimeError(
            f"Expected {expected_atoms} atoms; "
            f"found {len(output_atoms)}"
        )

    output_atoms_path = (
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2_atoms.csv"
    )
    write_csv(output_atoms_path, output_atoms)

    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2_caps.csv",
        cap_rows,
    )

    xyz_path = (
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2.xyz"
    )

    xyz_lines = [
        str(len(output_atoms)),
        (
            "QM_F06_LOWER_BOUNDARY_V2; Stage-2 bridge geometry; "
            "CAP07 replaced by real R2 coordination; unoptimized"
        ),
    ]

    for row in output_atoms:
        xyz_lines.append(
            f"{row['element']:<2s} "
            f"{float(row['x_angstrom']): .10f} "
            f"{float(row['y_angstrom']): .10f} "
            f"{float(row['z_angstrom']): .10f}"
        )

    xyz_path.write_text(
        "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )

    element_counts = Counter(
        row["element"]
        for row in output_atoms
    )

    role_counts = Counter(
        row["atom_role"]
        for row in output_atoms
    )

    report_path = (
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2_BUILD_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2 Build — Day027",
                "",
                "## Correction",
                "",
                "- Removed artificial cap: `HCAP:LOWER:07`",
                "- Restored real atom: `A:LOWER:13:-3`",
                (
                    "- Restored first-shell atoms: "
                    "`A:LOWER:11:-3`, `A:LOWER:14:-2`"
                ),
                "- Restored real passivant: `H4:LOWER:0016:0`",
                "- New peripheral caps: **3**",
                "",
                "## Result",
                "",
                f"- Total atoms: **{len(output_atoms)}**",
                f"- Element counts: `{dict(element_counts)}`",
                f"- Role counts: `{dict(role_counts)}`",
                "- Geometry optimized: **NO**",
                "- QM calculation executed: **NO**",
                "",
                "## Coordinate policy",
                "",
                (
                    "Atoms retained from V1 use the converged Stage-2 "
                    "geometry. Newly restored R2 atoms use validated "
                    "Day024 coordinates. New caps are placed along the "
                    "original R2 cut-bond vectors."
                ),
                "",
                "## Required next step",
                "",
                (
                    "Reconstruct connectivity and audit valence, bond "
                    "distances, boundary contacts and bridge proximity "
                    "before preparing a V2 QM optimization."
                ),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": (
            "QM_F06_LOWER_BOUNDARY_V2_BUILT_PRE_QM_AUDIT_REQUIRED"
        ),
        "atom_count": len(output_atoms),
        "element_counts": dict(element_counts),
        "removed_atoms": sorted(REMOVE_ATOMS),
        "added_real_atoms": sorted(ADD_REAL_ATOMS),
        "new_cap_count": len(cap_rows),
        "geometry_optimized": False,
        "qm_calculation_executed": False,
        "required_next_step": (
            "AUDIT_QM_F06_LOWER_BOUNDARY_V2"
        ),
    }

    (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_build_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Boundary V2 construction completed.")
    print("Atoms:", len(output_atoms))
    print("Elements:", dict(element_counts))
    print("New caps:", len(cap_rows))
    print("QM executed: False")
    print("Output:", xyz_path)


if __name__ == "__main__":
    main()
