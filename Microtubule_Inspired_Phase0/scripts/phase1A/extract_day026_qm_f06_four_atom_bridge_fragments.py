#!/usr/bin/env python3
"""
Extract representative QM_F06 four-atom B-N-B-N bridge fragments.

For each R2 end (LOWER and UPPER), the script:

1. Selects one bridge path deterministically.
2. Includes the complete four-atom bridge.
3. Includes both bridge attachment centers.
4. Includes direct neighbors of the bridge atoms.
5. Includes the first coordination shell of both attachment centers.
6. Preserves all existing hydrogen atoms already present in the graph.
7. Writes XYZ coordinates in angstrom.
8. Reports internal and cut/boundary edges.

No artificial capping atoms are added.
No geometry optimization, force-field assignment or QM calculation is run.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

NODES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

EDGES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

BRIDGE_TYPE = "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
ENDS = ("LOWER", "UPPER")


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No records found in {path}")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for output: {path}")

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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def natural_bridge_key(row: dict[str, str]) -> tuple[str, int, int]:
    path_id = row.get("bridge_path_id", "")
    bridge_index = int(float(row.get("bridge_index") or 0))
    bridge_position = int(float(row.get("bridge_position") or 0))
    return path_id, bridge_index, bridge_position


def write_xyz(
    path: Path,
    atom_rows: list[dict[str, Any]],
    comment: str,
) -> None:
    # Convert nm to angstrom and center the extracted fragment.
    x_values = [float(row["x_nm"]) * 10.0 for row in atom_rows]
    y_values = [float(row["y_nm"]) * 10.0 for row in atom_rows]
    z_values = [float(row["z_nm"]) * 10.0 for row in atom_rows]

    cx = sum(x_values) / len(x_values)
    cy = sum(y_values) / len(y_values)
    cz = sum(z_values) / len(z_values)

    lines = [str(len(atom_rows)), comment]

    for row, x, y, z in zip(
        atom_rows,
        x_values,
        y_values,
        z_values,
        strict=True,
    ):
        lines.append(
            f"{row['element']:<2s} "
            f"{x - cx: .10f} "
            f"{y - cy: .10f} "
            f"{z - cz: .10f}"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    coordinate_rows = read_csv(COORDINATES)
    node_rows = read_csv(NODES)
    edge_rows = read_csv(EDGES)

    coordinates = {
        row["node_id"]: row
        for row in coordinate_rows
    }
    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    if set(coordinates) != set(nodes):
        raise RuntimeError(
            "Coordinate and graph node IDs are not identical."
        )

    adjacency: dict[str, set[str]] = defaultdict(set)

    for edge in edge_rows:
        source = edge["source_node"]
        target = edge["target_node"]

        if source not in nodes or target not in nodes:
            raise RuntimeError(
                f"Edge references unknown node: {source} -- {target}"
            )

        adjacency[source].add(target)
        adjacency[target].add(source)

    manifest_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []

    for end in ENDS:
        bridge_candidates = [
            row
            for row in node_rows
            if row["node_type"] == BRIDGE_TYPE
            and row["end"] == end
        ]

        if not bridge_candidates:
            raise RuntimeError(
                f"No bridge atoms found for end {end}"
            )

        path_ids = sorted(
            {
                row["bridge_path_id"]
                for row in bridge_candidates
                if row["bridge_path_id"]
            }
        )

        if not path_ids:
            raise RuntimeError(
                f"No bridge_path_id values found for {end}"
            )

        selected_path_id = path_ids[0]

        bridge_rows = sorted(
            [
                row
                for row in bridge_candidates
                if row["bridge_path_id"] == selected_path_id
            ],
            key=natural_bridge_key,
        )

        bridge_ids = {
            row["node_id"]
            for row in bridge_rows
        }

        if len(bridge_ids) != 4:
            raise RuntimeError(
                f"{end}/{selected_path_id}: expected 4 bridge atoms; "
                f"found {len(bridge_ids)}"
            )

        attachment_ids: set[str] = set()

        for row in bridge_rows:
            for field in (
                "attached_seed_node",
                "attached_annulus_node",
            ):
                node_id = row.get(field, "").strip()
                if node_id:
                    attachment_ids.add(node_id)

        # Fallback based on graph connectivity if metadata is incomplete.
        for node_id in bridge_ids:
            for neighbor in adjacency[node_id]:
                if neighbor not in bridge_ids:
                    neighbor_type = nodes[neighbor]["node_type"]
                    if neighbor_type not in {
                        "BRIDGE_PASSIVANT_H",
                    }:
                        attachment_ids.add(neighbor)

        if len(attachment_ids) != 2:
            raise RuntimeError(
                f"{end}/{selected_path_id}: expected exactly two "
                f"attachment centers; found {len(attachment_ids)}: "
                f"{sorted(attachment_ids)}"
            )

        selected_ids = set(bridge_ids)
        selected_ids.update(attachment_ids)

        # Direct neighbors of bridge atoms: existing passivants and attachments.
        for node_id in bridge_ids:
            selected_ids.update(adjacency[node_id])

        # First coordination shell around both attachment centers.
        for node_id in attachment_ids:
            selected_ids.update(adjacency[node_id])

        internal_edges: list[dict[str, Any]] = []
        boundary_edges: list[dict[str, Any]] = []

        for edge in edge_rows:
            source = edge["source_node"]
            target = edge["target_node"]

            source_inside = source in selected_ids
            target_inside = target in selected_ids

            if source_inside and target_inside:
                internal_edges.append(edge)
            elif source_inside != target_inside:
                inside_node = source if source_inside else target
                outside_node = target if source_inside else source

                boundary_edges.append(
                    {
                        **edge,
                        "inside_node": inside_node,
                        "outside_node": outside_node,
                        "inside_element": nodes[inside_node]["element"],
                        "outside_element": nodes[outside_node]["element"],
                        "inside_node_type": nodes[inside_node]["node_type"],
                        "outside_node_type": nodes[outside_node]["node_type"],
                        "capping_required_review": True,
                    }
                )

        atom_rows: list[dict[str, Any]] = []

        for node_id in sorted(selected_ids):
            coord = coordinates[node_id]
            meta = nodes[node_id]

            atom_rows.append(
                {
                    "node_id": node_id,
                    "element": meta["element"],
                    "node_type": meta["node_type"],
                    "end": meta["end"],
                    "bridge_path_id": meta.get("bridge_path_id", ""),
                    "bridge_index": meta.get("bridge_index", ""),
                    "bridge_position": meta.get("bridge_position", ""),
                    "is_bridge_atom": node_id in bridge_ids,
                    "is_attachment_center": node_id in attachment_ids,
                    "x_nm": coord["x_nm"],
                    "y_nm": coord["y_nm"],
                    "z_nm": coord["z_nm"],
                    "coordinate_source": coord["coordinate_source"],
                    "energy_minimized": coord["energy_minimized"],
                    "MD_relaxed": coord["MD_relaxed"],
                }
            )

        label = f"QM_F06_{end}"
        xyz_path = OUTPUT_DIR / f"{label}.xyz"
        atoms_path = OUTPUT_DIR / f"{label}_atoms.csv"
        edges_path = OUTPUT_DIR / f"{label}_internal_edges.csv"
        boundary_path = OUTPUT_DIR / f"{label}_boundary_edges.csv"

        write_xyz(
            xyz_path,
            atom_rows,
            (
                f"{label}; bridge_path_id={selected_path_id}; "
                "coordinates from Day024 Gate3P2; unoptimized; "
                "no artificial capping atoms"
            ),
        )
        write_csv(atoms_path, atom_rows)
        write_csv(edges_path, internal_edges)
        write_csv(boundary_path, boundary_edges)

        element_counts: dict[str, int] = defaultdict(int)
        type_counts: dict[str, int] = defaultdict(int)

        for row in atom_rows:
            element_counts[row["element"]] += 1
            type_counts[row["node_type"]] += 1

        manifest_rows.extend(
            [
                {
                    "fragment": label,
                    "role": "xyz_geometry",
                    "file": str(xyz_path.relative_to(ROOT)),
                    "sha256": sha256(xyz_path),
                },
                {
                    "fragment": label,
                    "role": "atom_manifest",
                    "file": str(atoms_path.relative_to(ROOT)),
                    "sha256": sha256(atoms_path),
                },
                {
                    "fragment": label,
                    "role": "internal_edges",
                    "file": str(edges_path.relative_to(ROOT)),
                    "sha256": sha256(edges_path),
                },
                {
                    "fragment": label,
                    "role": "boundary_edges",
                    "file": str(boundary_path.relative_to(ROOT)),
                    "sha256": sha256(boundary_path),
                },
            ]
        )

        report_sections.append(
            "\n".join(
                [
                    f"## {label}",
                    "",
                    f"- Selected bridge path: `{selected_path_id}`",
                    f"- Bridge atoms: **{len(bridge_ids)}**",
                    f"- Attachment centers: **{len(attachment_ids)}**",
                    f"- Total extracted atoms: **{len(atom_rows)}**",
                    f"- Internal edges: **{len(internal_edges)}**",
                    f"- Boundary/cut edges: **{len(boundary_edges)}**",
                    f"- Element counts: `{dict(sorted(element_counts.items()))}`",
                    f"- Node-type counts: `{dict(sorted(type_counts.items()))}`",
                    "- Artificial capping atoms added: **NO**",
                    "- Geometry optimized: **NO**",
                    "- QM calculation authorized or executed: **NO**",
                    "",
                ]
            )
        )

    manifest_path = OUTPUT_DIR / "QM_F06_extraction_manifest.csv"
    write_csv(manifest_path, manifest_rows)

    report_path = OUTPUT_DIR / "QM_F06_EXTRACTION_REPORT.md"
    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 Four-Atom Bridge Fragment Extraction — Day026",
                "",
                "## Scope",
                "",
                (
                    "Representative LOWER and UPPER B–N–B–N bridge "
                    "fragments were extracted from the validated Day024 "
                    "Gate3P2 coordinate set and the Gate3M graph topology."
                ),
                "",
                (
                    "Each fragment includes the complete four-atom bridge, "
                    "both attachment centers, direct bridge neighbors and "
                    "the first coordination shell of both attachment centers."
                ),
                "",
                (
                    "Boundary edges are reported explicitly. No artificial "
                    "hydrogen capping, geometry optimization, force-field "
                    "assignment or QM calculation was performed."
                ),
                "",
                *report_sections,
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": "QM_F06_REPRESENTATIVE_FRAGMENTS_EXTRACTED_CAP_AUDIT_REQUIRED",
        "ends_extracted": list(ENDS),
        "artificial_capping_added": False,
        "geometry_optimized": False,
        "qm_calculation_executed": False,
        "required_next_step": (
            "AUDIT_BOUNDARY_EDGES_AND_DEFINE_CHEMICALLY_VALID_CAPPING"
        ),
    }

    json_path = OUTPUT_DIR / "QM_F06_extraction_summary.json"
    json_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 extraction completed.")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Manifest records: {len(manifest_rows)}")
    print(
        "Next step: audit boundary edges before adding any capping atoms."
    )


if __name__ == "__main__":
    main()
