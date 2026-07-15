#!/usr/bin/env python3
"""
Build chemically capped QM_F06 bridge fragments.

Workflow
--------
1. Load the previously extracted LOWER and UPPER QM_F06 fragments.
2. Incorporate existing hydrogen atoms identified by the boundary audit.
3. Replace each peripheral cut B-N bond with one artificial H cap attached
   to the atom retained inside the fragment.
4. Place each cap along the original inside-to-outside bond vector.
5. Use target bond lengths already present in the validated R2 geometry:
      B-H = 1.19 angstrom
      N-H = 1.01 angstrom
6. Generate capped XYZ, atom manifests, cap manifests and a validation report.

No optimization, charge assignment, force-field fitting or QM calculation
is performed or authorized by this script.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
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

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")

TARGET_XH_ANGSTROM = {
    "B": 1.19,
    "N": 1.01,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def coordinates_angstrom(row: dict[str, str]) -> tuple[float, float, float]:
    return (
        float(row["x_nm"]) * 10.0,
        float(row["y_nm"]) * 10.0,
        float(row["z_nm"]) * 10.0,
    )


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second, strict=True)
        )
    )


def place_cap(
    inside_xyz: tuple[float, float, float],
    outside_xyz: tuple[float, float, float],
    target_length: float,
) -> tuple[float, float, float]:
    vector = tuple(
        outside - inside
        for inside, outside in zip(
            inside_xyz,
            outside_xyz,
            strict=True,
        )
    )

    norm = math.sqrt(sum(component * component for component in vector))

    if norm <= 1.0e-12:
        raise RuntimeError("Zero-length cut-bond vector encountered.")

    unit = tuple(component / norm for component in vector)

    return tuple(
        inside + target_length * direction
        for inside, direction in zip(
            inside_xyz,
            unit,
            strict=True,
        )
    )


def write_xyz(
    path: Path,
    atoms: list[dict[str, Any]],
    comment: str,
) -> None:
    x_values = [float(row["x_angstrom"]) for row in atoms]
    y_values = [float(row["y_angstrom"]) for row in atoms]
    z_values = [float(row["z_angstrom"]) for row in atoms]

    center = (
        sum(x_values) / len(atoms),
        sum(y_values) / len(atoms),
        sum(z_values) / len(atoms),
    )

    lines = [str(len(atoms)), comment]

    for row in atoms:
        lines.append(
            f"{row['element']:<2s} "
            f"{float(row['x_angstrom']) - center[0]: .10f} "
            f"{float(row['y_angstrom']) - center[1]: .10f} "
            f"{float(row['z_angstrom']) - center[2]: .10f}"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    coordinate_rows = read_csv(COORDINATES)
    node_rows = read_csv(NODES)

    coordinates = {
        row["node_id"]: row
        for row in coordinate_rows
    }

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    manifest_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []

    for end in ENDS:
        original_atoms = read_csv(
            F06_DIR / f"QM_F06_{end}_atoms.csv"
        )

        boundary_audit = read_csv(
            F06_DIR / f"QM_F06_{end}_boundary_edge_audit.csv"
        )

        final_atoms: list[dict[str, Any]] = []
        selected_ids: set[str] = set()

        # Preserve every atom from the original extraction.
        for atom in original_atoms:
            node_id = atom["node_id"]
            xyz = coordinates_angstrom(coordinates[node_id])

            selected_ids.add(node_id)

            final_atoms.append(
                {
                    "atom_id": node_id,
                    "element": atom["element"],
                    "atom_role": "ORIGINAL_FRAGMENT_ATOM",
                    "node_type": atom["node_type"],
                    "source_edge_id": "",
                    "parent_inside_node": "",
                    "x_angstrom": f"{xyz[0]:.10f}",
                    "y_angstrom": f"{xyz[1]:.10f}",
                    "z_angstrom": f"{xyz[2]:.10f}",
                    "artificial_cap": False,
                }
            )

        cap_rows: list[dict[str, Any]] = []
        included_existing_h: list[str] = []

        cap_index = 0

        for boundary in boundary_audit:
            action = boundary["preliminary_action"]
            inside_node = boundary["inside_node"]
            outside_node = boundary["outside_node"]

            if action == "INCLUDE_EXISTING_HYDROGEN":
                if outside_node in selected_ids:
                    continue

                outside_meta = nodes[outside_node]

                if outside_meta["element"] != "H":
                    raise RuntimeError(
                        f"Expected outside H atom, found "
                        f"{outside_meta['element']} for {outside_node}"
                    )

                xyz = coordinates_angstrom(coordinates[outside_node])

                final_atoms.append(
                    {
                        "atom_id": outside_node,
                        "element": "H",
                        "atom_role": "EXISTING_R2_HYDROGEN_ADDED",
                        "node_type": outside_meta["node_type"],
                        "source_edge_id": boundary["edge_id"],
                        "parent_inside_node": inside_node,
                        "x_angstrom": f"{xyz[0]:.10f}",
                        "y_angstrom": f"{xyz[1]:.10f}",
                        "z_angstrom": f"{xyz[2]:.10f}",
                        "artificial_cap": False,
                    }
                )

                selected_ids.add(outside_node)
                included_existing_h.append(outside_node)
                continue

            if action != "CANDIDATE_BN_CUT_FOR_HYDROGEN_CAPPING":
                raise RuntimeError(
                    f"Unsupported boundary action for {end}: {action}"
                )

            inside_element = nodes[inside_node]["element"]
            outside_element = nodes[outside_node]["element"]

            if inside_element not in TARGET_XH_ANGSTROM:
                raise RuntimeError(
                    f"Cannot cap inside element {inside_element}: "
                    f"{inside_node}"
                )

            if {inside_element, outside_element} != {"B", "N"}:
                raise RuntimeError(
                    f"Expected B-N cut edge, found "
                    f"{inside_element}-{outside_element}"
                )

            inside_xyz = coordinates_angstrom(coordinates[inside_node])
            outside_xyz = coordinates_angstrom(coordinates[outside_node])

            target_length = TARGET_XH_ANGSTROM[inside_element]

            cap_xyz = place_cap(
                inside_xyz,
                outside_xyz,
                target_length,
            )

            cap_index += 1
            cap_id = f"HCAP:{end}:{cap_index:02d}"

            final_atoms.append(
                {
                    "atom_id": cap_id,
                    "element": "H",
                    "atom_role": "ARTIFICIAL_BOUNDARY_CAP",
                    "node_type": "QM_BOUNDARY_CAP_H",
                    "source_edge_id": boundary["edge_id"],
                    "parent_inside_node": inside_node,
                    "x_angstrom": f"{cap_xyz[0]:.10f}",
                    "y_angstrom": f"{cap_xyz[1]:.10f}",
                    "z_angstrom": f"{cap_xyz[2]:.10f}",
                    "artificial_cap": True,
                }
            )

            cap_rows.append(
                {
                    "fragment": f"QM_F06_{end}",
                    "cap_id": cap_id,
                    "source_edge_id": boundary["edge_id"],
                    "parent_inside_node": inside_node,
                    "parent_element": inside_element,
                    "removed_outside_node": outside_node,
                    "removed_outside_element": outside_element,
                    "original_BN_distance_angstrom": (
                        boundary["original_distance_angstrom"]
                    ),
                    "target_XH_distance_angstrom": (
                        f"{target_length:.8f}"
                    ),
                    "placement_rule": (
                        "ALONG_ORIGINAL_INSIDE_TO_OUTSIDE_BOND_VECTOR"
                    ),
                    "geometry_optimized": False,
                    "capping_method_status": (
                        "GEOMETRIC_INITIAL_GUESS_REQUIRES_QM_VALIDATION"
                    ),
                }
            )

        label = f"QM_F06_{end}_CAPPED"

        atoms_path = F06_DIR / f"{label}_atoms.csv"
        caps_path = F06_DIR / f"{label}_caps.csv"
        xyz_path = F06_DIR / f"{label}.xyz"

        write_csv(atoms_path, final_atoms)
        write_csv(caps_path, cap_rows)

        write_xyz(
            xyz_path,
            final_atoms,
            (
                f"{label}; includes existing R2 H atoms and artificial "
                "boundary H caps; unoptimized; no QM calculation executed"
            ),
        )

        # Validation: all artificial X-H distances must equal target values.
        validation_rows: list[dict[str, Any]] = []

        atom_lookup = {
            row["atom_id"]: row
            for row in final_atoms
        }

        for cap in cap_rows:
            parent = atom_lookup[cap["parent_inside_node"]]
            capped_h = atom_lookup[cap["cap_id"]]

            parent_xyz = (
                float(parent["x_angstrom"]),
                float(parent["y_angstrom"]),
                float(parent["z_angstrom"]),
            )
            cap_xyz = (
                float(capped_h["x_angstrom"]),
                float(capped_h["y_angstrom"]),
                float(capped_h["z_angstrom"]),
            )

            measured = distance(parent_xyz, cap_xyz)
            target = float(cap["target_XH_distance_angstrom"])
            error = abs(measured - target)

            validation_rows.append(
                {
                    "fragment": f"QM_F06_{end}",
                    "cap_id": cap["cap_id"],
                    "parent_inside_node": cap["parent_inside_node"],
                    "parent_element": cap["parent_element"],
                    "target_distance_angstrom": f"{target:.10f}",
                    "measured_distance_angstrom": f"{measured:.10f}",
                    "absolute_error_angstrom": f"{error:.12e}",
                    "distance_gate_pass": error <= 1.0e-8,
                }
            )

        validation_path = F06_DIR / (
            f"{label}_cap_distance_validation.csv"
        )
        write_csv(validation_path, validation_rows)

        failures = [
            row
            for row in validation_rows
            if str(row["distance_gate_pass"]).lower() != "true"
        ]

        if failures:
            raise RuntimeError(
                f"{end}: artificial cap distance validation failed."
            )

        element_counts = Counter(
            row["element"]
            for row in final_atoms
        )

        report_sections.extend(
            [
                f"## {label}",
                "",
                f"- Original extracted atoms: **{len(original_atoms)}**",
                (
                    f"- Existing R2 hydrogen atoms incorporated: "
                    f"**{len(included_existing_h)}**"
                ),
                f"- Artificial caps added: **{len(cap_rows)}**",
                f"- Final atoms: **{len(final_atoms)}**",
                f"- Element counts: `{dict(sorted(element_counts.items()))}`",
                "- Cap-distance validation failures: **0**",
                "- Geometry optimized: **NO**",
                "- QM calculation executed: **NO**",
                "",
            ]
        )

        for role, path in (
            ("capped_xyz", xyz_path),
            ("capped_atom_manifest", atoms_path),
            ("cap_manifest", caps_path),
            ("cap_distance_validation", validation_path),
        ):
            manifest_rows.append(
                {
                    "fragment": label,
                    "role": role,
                    "file": str(path.relative_to(ROOT)),
                    "sha256": sha256(path),
                }
            )

    manifest_path = F06_DIR / "QM_F06_capped_fragment_manifest.csv"
    write_csv(manifest_path, manifest_rows)

    report_path = F06_DIR / "QM_F06_CAPPED_FRAGMENT_REPORT.md"
    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 Chemically Capped Fragments — Day026",
                "",
                "## Decision",
                "",
                (
                    "The original 14-atom bridge fragments required only "
                    "minimal boundary completion. Existing R2 hydrogen atoms "
                    "were restored, and peripheral B-N cuts were saturated "
                    "with artificial H caps along the original bond vectors."
                ),
                "",
                (
                    "No additional full coordination shell was required "
                    "because none of the cut edges touched the bridge core "
                    "or either attachment center."
                ),
                "",
                "## Capping geometry",
                "",
                "- B-H target distance: **1.19 Å**",
                "- N-H target distance: **1.01 Å**",
                (
                    "- Placement: along each original inside-to-outside "
                    "B-N bond vector."
                ),
                "",
                *report_sections,
                "## Authorization state",
                "",
                "- Fragment construction: **COMPLETED**",
                "- Artificial caps: **GEOMETRIC INITIAL GUESSES ONLY**",
                "- Geometry optimization authorized: **NO**",
                "- QM calculation executed: **NO**",
                "",
                "## Required next step",
                "",
                (
                    "Audit capped-fragment valence, atom count, minimum "
                    "interatomic distances and net charge/multiplicity "
                    "requirements before preparing electronic-structure inputs."
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    summary = {
        "decision": (
            "QM_F06_CAPPED_GEOMETRIES_BUILT_PRE_QM_VALIDATION_REQUIRED"
        ),
        "ends": list(ENDS),
        "target_bh_angstrom": TARGET_XH_ANGSTROM["B"],
        "target_nh_angstrom": TARGET_XH_ANGSTROM["N"],
        "geometry_optimized": False,
        "qm_calculation_executed": False,
        "required_next_step": (
            "AUDIT_CAPPED_FRAGMENT_VALENCE_CLASHES_CHARGE_AND_MULTIPLICITY"
        ),
    }

    (
        F06_DIR / "QM_F06_capped_fragment_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 capped-fragment construction completed.")
    print(f"Output directory: {F06_DIR}")
    print(
        "Next step: pre-QM valence, clash, charge and multiplicity audit."
    )


if __name__ == "__main__":
    main()
