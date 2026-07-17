#!/usr/bin/env python3
"""
Chemical audit of the only remaining topology-aware hard contact in the
QM_F06 LOWER Boundary V2-B optimized geometry:

    BR4:LOWER:00:3 — H4:LOWER:0017:0

The script determines:

1. covalent parent and graph neighborhood of the hydrogen;
2. graph separation between the contact atoms;
3. contact evolution from V2-A to V2-B;
4. local bonded distances around both atoms;
5. local angles involving the residual-contact atoms;
6. local torsions spanning the residual-contact region;
7. nearest nonbonded partners of both atoms;
8. whether the contact satisfies a conservative unintended-covalent
   threshold;
9. whether the contact is improving, worsening or unchanged.

No coordinates are changed and no QM calculation is executed.
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

STATE_PATH = V2_DIR / (
    "orca_v2_workflow/v2_workflow_state.json"
)

ATOM_MANIFEST = V2_DIR / (
    "cap02_repair/"
    "QM_F06_LOWER_BOUNDARY_V2_REPAIRED_atoms.csv"
)

FULL_EDGES_PATH = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_postprocessing/"
    "residual_real_contact_audit"
)

BRIDGE_ATOM = "BR4:LOWER:00:3"
HYDROGEN_ATOM = "H4:LOWER:0017:0"

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

# Conservative geometric thresholds used only to flag possible
# unintended covalent connectivity. They do not assign a bond.
COVALENT_THRESHOLDS = {
    tuple(sorted(("B", "H"))): 1.45,
    tuple(sorted(("N", "H"))): 1.30,
    tuple(sorted(("B", "N"))): 1.85,
    tuple(sorted(("B", "B"))): 1.90,
    tuple(sorted(("N", "N"))): 1.70,
    tuple(sorted(("H", "H"))): 0.90,
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
        return

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


def read_xyz(
    path: Path,
) -> list[tuple[str, float, float, float]]:
    require_file(path)

    lines = path.read_text(encoding="utf-8").splitlines()
    expected = int(lines[0].strip())

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

    if len(rows) != expected:
        raise RuntimeError(
            f"XYZ count mismatch in {path}: "
            f"{expected} vs {len(rows)}"
        )

    return rows


def distance(first, second) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second, strict=True)
        )
    )


def vector(origin, target):
    return tuple(
        target_value - origin_value
        for origin_value, target_value
        in zip(origin, target, strict=True)
    )


def dot(first, second) -> float:
    return sum(
        a * b
        for a, b in zip(first, second, strict=True)
    )


def cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def norm(value) -> float:
    return math.sqrt(dot(value, value))


def angle_degrees(first, center, third) -> float:
    vector_1 = vector(center, first)
    vector_2 = vector(center, third)

    denominator = norm(vector_1) * norm(vector_2)

    if denominator <= 1.0e-14:
        raise RuntimeError("Undefined angle from zero-length vector.")

    cosine = max(
        -1.0,
        min(1.0, dot(vector_1, vector_2) / denominator),
    )

    return math.degrees(math.acos(cosine))


def dihedral_degrees(first, second, third, fourth) -> float:
    bond_1 = vector(first, second)
    bond_2 = vector(second, third)
    bond_3 = vector(third, fourth)

    normal_1 = cross(bond_1, bond_2)
    normal_2 = cross(bond_2, bond_3)

    normal_1_norm = norm(normal_1)
    normal_2_norm = norm(normal_2)
    bond_2_norm = norm(bond_2)

    if min(normal_1_norm, normal_2_norm, bond_2_norm) <= 1.0e-14:
        raise RuntimeError("Undefined dihedral from collinear geometry.")

    unit_normal_1 = tuple(
        value / normal_1_norm for value in normal_1
    )
    unit_normal_2 = tuple(
        value / normal_2_norm for value in normal_2
    )
    unit_bond_2 = tuple(
        value / bond_2_norm for value in bond_2
    )

    helper = cross(unit_normal_1, unit_bond_2)

    x_value = dot(unit_normal_1, unit_normal_2)
    y_value = dot(helper, unit_normal_2)

    return math.degrees(math.atan2(y_value, x_value))


def canonical_edge(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def shortest_path(
    adjacency: dict[str, set[str]],
    source: str,
    target: str,
) -> list[str] | None:
    if source == target:
        return [source]

    visited = {source}
    queue = deque([(source, [source])])

    while queue:
        node, path = queue.popleft()

        for neighbor in sorted(adjacency[node]):
            if neighbor == target:
                return path + [neighbor]

            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    require_file(STATE_PATH)

    state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    required_state = {
        "v2a_executed": state.get("v2a_executed") is True,
        "v2b_executed": state.get("v2b_executed") is True,
        "v2b_validation_pass": (
            state.get("v2b_validation_pass") is True
        ),
        "v2a_xyz_present": bool(
            state.get("v2a_optimized_xyz")
        ),
        "v2b_xyz_present": bool(
            state.get("v2b_optimized_xyz")
        ),
    }

    if not all(required_state.values()):
        raise RuntimeError(
            f"Incomplete V2 workflow state: {required_state}"
        )

    atom_rows = read_csv(ATOM_MANIFEST)

    if len(atom_rows) != 28:
        raise RuntimeError(
            f"Expected 28 atoms; found {len(atom_rows)}"
        )

    atom_ids = [row["atom_id"] for row in atom_rows]
    atom_set = set(atom_ids)

    if BRIDGE_ATOM not in atom_set:
        raise RuntimeError(f"Missing bridge atom: {BRIDGE_ATOM}")

    if HYDROGEN_ATOM not in atom_set:
        raise RuntimeError(
            f"Missing hydrogen atom: {HYDROGEN_ATOM}"
        )

    atoms = {
        row["atom_id"]: row
        for row in atom_rows
    }

    elements = {
        row["atom_id"]: row["element"]
        for row in atom_rows
    }

    v2a_xyz_path = ROOT / state["v2a_optimized_xyz"]
    v2b_xyz_path = ROOT / state["v2b_optimized_xyz"]

    v2a_xyz_rows = read_xyz(v2a_xyz_path)
    v2b_xyz_rows = read_xyz(v2b_xyz_path)

    if (
        len(v2a_xyz_rows) != len(atom_ids)
        or len(v2b_xyz_rows) != len(atom_ids)
    ):
        raise RuntimeError("XYZ/manifest atom-count mismatch.")

    v2a_coords = {}
    v2b_coords = {}

    for atom_id, atom, v2a, v2b in zip(
        atom_ids,
        atom_rows,
        v2a_xyz_rows,
        v2b_xyz_rows,
        strict=True,
    ):
        if atom["element"] != v2a[0]:
            raise RuntimeError(
                f"V2-A element mismatch for {atom_id}"
            )

        if atom["element"] != v2b[0]:
            raise RuntimeError(
                f"V2-B element mismatch for {atom_id}"
            )

        v2a_coords[atom_id] = v2a[1:]
        v2b_coords[atom_id] = v2b[1:]

    edge_records: dict[
        tuple[str, str],
        dict[str, str],
    ] = {}

    for row in read_csv(FULL_EDGES_PATH):
        first = row["source_node"]
        second = row["target_node"]

        if first in atom_set and second in atom_set:
            edge_records[
                canonical_edge(first, second)
            ] = {
                "edge_origin": "REAL_R2_GRAPH_EDGE",
                "source_edge_id": row["edge_id"],
                "edge_type": row["edge_type"],
            }

    for row in atom_rows:
        if row["artificial_cap"].lower() != "true":
            continue

        cap = row["atom_id"]
        parent = row["parent_inside_node"]

        if parent not in atom_set:
            raise RuntimeError(
                f"Cap parent absent: {cap} -> {parent}"
            )

        edge_records[
            canonical_edge(cap, parent)
        ] = {
            "edge_origin": "ARTIFICIAL_CAP_EDGE",
            "source_edge_id": row["source_edge_id"],
            "edge_type": "QM_BOUNDARY_CAP",
        }

    adjacency: dict[str, set[str]] = defaultdict(set)

    for first, second in edge_records:
        adjacency[first].add(second)
        adjacency[second].add(first)

    hydrogen_neighbors = sorted(
        adjacency[HYDROGEN_ATOM]
    )

    if len(hydrogen_neighbors) != 1:
        raise RuntimeError(
            f"Expected one covalent parent for {HYDROGEN_ATOM}; "
            f"found {hydrogen_neighbors}"
        )

    hydrogen_parent = hydrogen_neighbors[0]

    path = shortest_path(
        adjacency,
        BRIDGE_ATOM,
        HYDROGEN_ATOM,
    )

    if path is None:
        raise RuntimeError(
            "Residual-contact atoms are disconnected."
        )

    graph_separation = len(path) - 1

    residual_pair = tuple(
        sorted(
            (
                elements[BRIDGE_ATOM],
                elements[HYDROGEN_ATOM],
            )
        )
    )

    residual_threshold = COVALENT_THRESHOLDS[
        residual_pair
    ]

    v2a_contact_distance = distance(
        v2a_coords[BRIDGE_ATOM],
        v2a_coords[HYDROGEN_ATOM],
    )

    v2b_contact_distance = distance(
        v2b_coords[BRIDGE_ATOM],
        v2b_coords[HYDROGEN_ATOM],
    )

    vdw_sum = (
        VDW_RADII[elements[BRIDGE_ATOM]]
        + VDW_RADII[elements[HYDROGEN_ATOM]]
    )

    v2a_ratio = v2a_contact_distance / vdw_sum
    v2b_ratio = v2b_contact_distance / vdw_sum

    distance_change = (
        v2b_contact_distance - v2a_contact_distance
    )

    compression_improved = distance_change > 0.0

    possible_unintended_covalent_contact = (
        v2b_contact_distance <= residual_threshold
    )

    # Local bonded environment around both contact atoms.
    local_atoms = {
        BRIDGE_ATOM,
        HYDROGEN_ATOM,
        hydrogen_parent,
        *adjacency[BRIDGE_ATOM],
        *adjacency[hydrogen_parent],
    }

    local_bond_rows = []

    for edge in sorted(edge_records):
        first, second = edge

        if (
            first not in local_atoms
            and second not in local_atoms
        ):
            continue

        if (
            first not in {
                BRIDGE_ATOM,
                HYDROGEN_ATOM,
                hydrogen_parent,
            }
            and second not in {
                BRIDGE_ATOM,
                HYDROGEN_ATOM,
                hydrogen_parent,
            }
        ):
            continue

        v2a_length = distance(
            v2a_coords[first],
            v2a_coords[second],
        )

        v2b_length = distance(
            v2b_coords[first],
            v2b_coords[second],
        )

        local_bond_rows.append(
            {
                "atom_1": first,
                "element_1": elements[first],
                "atom_2": second,
                "element_2": elements[second],
                "edge_origin": (
                    edge_records[edge]["edge_origin"]
                ),
                "v2a_distance_angstrom": (
                    f"{v2a_length:.10f}"
                ),
                "v2b_distance_angstrom": (
                    f"{v2b_length:.10f}"
                ),
                "distance_change_angstrom": (
                    f"{v2b_length - v2a_length:.10f}"
                ),
                "touches_bridge_atom": (
                    BRIDGE_ATOM in edge
                ),
                "touches_hydrogen_parent": (
                    hydrogen_parent in edge
                ),
                "is_hydrogen_parent_bond": (
                    HYDROGEN_ATOM in edge
                ),
            }
        )

    # H-parent-X angles define whether H rotates toward or away from
    # the bridge-contact atom.
    local_angle_rows = []

    for neighbor in sorted(adjacency[hydrogen_parent]):
        if neighbor == HYDROGEN_ATOM:
            continue

        v2a_angle = angle_degrees(
            v2a_coords[HYDROGEN_ATOM],
            v2a_coords[hydrogen_parent],
            v2a_coords[neighbor],
        )

        v2b_angle = angle_degrees(
            v2b_coords[HYDROGEN_ATOM],
            v2b_coords[hydrogen_parent],
            v2b_coords[neighbor],
        )

        local_angle_rows.append(
            {
                "atom_1": HYDROGEN_ATOM,
                "center_atom": hydrogen_parent,
                "atom_3": neighbor,
                "angle_pattern": (
                    f"{elements[HYDROGEN_ATOM]}-"
                    f"{elements[hydrogen_parent]}-"
                    f"{elements[neighbor]}"
                ),
                "v2a_angle_deg": f"{v2a_angle:.10f}",
                "v2b_angle_deg": f"{v2b_angle:.10f}",
                "angle_change_deg": (
                    f"{v2b_angle - v2a_angle:.10f}"
                ),
            }
        )

    # Direct orientation angle:
    # H-parent -> H vector against H-parent -> bridge-atom vector.
    v2a_parent_h_bridge_angle = angle_degrees(
        v2a_coords[HYDROGEN_ATOM],
        v2a_coords[hydrogen_parent],
        v2a_coords[BRIDGE_ATOM],
    )

    v2b_parent_h_bridge_angle = angle_degrees(
        v2b_coords[HYDROGEN_ATOM],
        v2b_coords[hydrogen_parent],
        v2b_coords[BRIDGE_ATOM],
    )

    orientation_rows = [
        {
            "hydrogen": HYDROGEN_ATOM,
            "parent": hydrogen_parent,
            "bridge_atom": BRIDGE_ATOM,
            "v2a_parent_H_vs_parent_bridge_angle_deg": (
                f"{v2a_parent_h_bridge_angle:.10f}"
            ),
            "v2b_parent_H_vs_parent_bridge_angle_deg": (
                f"{v2b_parent_h_bridge_angle:.10f}"
            ),
            "angle_change_deg": (
                f"{v2b_parent_h_bridge_angle - v2a_parent_h_bridge_angle:.10f}"
            ),
            "v2a_contact_distance_angstrom": (
                f"{v2a_contact_distance:.10f}"
            ),
            "v2b_contact_distance_angstrom": (
                f"{v2b_contact_distance:.10f}"
            ),
        }
    ]

    # Torsions across first-shell paths containing either the bridge
    # atom or hydrogen-parent bond.
    torsion_paths: set[tuple[str, str, str, str]] = set()

    for first in atom_ids:
        for second in adjacency[first]:
            for third in adjacency[second]:
                if third == first:
                    continue

                for fourth in adjacency[third]:
                    if fourth in {first, second}:
                        continue

                    path_4 = (
                        first,
                        second,
                        third,
                        fourth,
                    )

                    reverse = tuple(reversed(path_4))
                    canonical = min(path_4, reverse)

                    if (
                        BRIDGE_ATOM in canonical
                        and (
                            HYDROGEN_ATOM in canonical
                            or hydrogen_parent in canonical
                        )
                    ):
                        torsion_paths.add(canonical)

    local_torsion_rows = []

    for first, second, third, fourth in sorted(
        torsion_paths
    ):
        try:
            v2a_torsion = dihedral_degrees(
                v2a_coords[first],
                v2a_coords[second],
                v2a_coords[third],
                v2a_coords[fourth],
            )

            v2b_torsion = dihedral_degrees(
                v2b_coords[first],
                v2b_coords[second],
                v2b_coords[third],
                v2b_coords[fourth],
            )
        except RuntimeError:
            continue

        raw_change = v2b_torsion - v2a_torsion

        wrapped_change = (
            (raw_change + 180.0) % 360.0
        ) - 180.0

        local_torsion_rows.append(
            {
                "atom_1": first,
                "atom_2": second,
                "atom_3": third,
                "atom_4": fourth,
                "torsion_pattern": "-".join(
                    elements[atom]
                    for atom in (
                        first,
                        second,
                        third,
                        fourth,
                    )
                ),
                "v2a_torsion_deg": (
                    f"{v2a_torsion:.10f}"
                ),
                "v2b_torsion_deg": (
                    f"{v2b_torsion:.10f}"
                ),
                "wrapped_change_deg": (
                    f"{wrapped_change:.10f}"
                ),
            }
        )

    # Rank nonbonded neighbors for both target atoms.
    nearest_rows = []

    for target in (BRIDGE_ATOM, HYDROGEN_ATOM):
        candidates = []

        for other in atom_ids:
            if other == target:
                continue

            if other in adjacency[target]:
                continue

            target_path = shortest_path(
                adjacency,
                target,
                other,
            )

            separation = (
                len(target_path) - 1
                if target_path is not None
                else None
            )

            if separation in {2, 3}:
                continue

            measured = distance(
                v2b_coords[target],
                v2b_coords[other],
            )

            ratio = measured / (
                VDW_RADII[elements[target]]
                + VDW_RADII[elements[other]]
            )

            candidates.append(
                {
                    "target_atom": target,
                    "target_element": elements[target],
                    "partner_atom": other,
                    "partner_element": elements[other],
                    "graph_separation": (
                        separation
                        if separation is not None
                        else "DISCONNECTED"
                    ),
                    "v2b_distance_angstrom": (
                        f"{measured:.10f}"
                    ),
                    "distance_over_vdw_sum": (
                        f"{ratio:.10f}"
                    ),
                    "partner_role": (
                        atoms[other]["atom_role"]
                    ),
                    "partner_node_type": (
                        atoms[other]["node_type"]
                    ),
                }
            )

        candidates.sort(
            key=lambda row: float(
                row["v2b_distance_angstrom"]
            )
        )

        nearest_rows.extend(candidates[:10])

    contact_rows = [
        {
            "bridge_atom": BRIDGE_ATOM,
            "bridge_element": elements[BRIDGE_ATOM],
            "hydrogen_atom": HYDROGEN_ATOM,
            "hydrogen_parent": hydrogen_parent,
            "hydrogen_parent_element": (
                elements[hydrogen_parent]
            ),
            "graph_separation": graph_separation,
            "shortest_graph_path": "|".join(path),
            "v2a_distance_angstrom": (
                f"{v2a_contact_distance:.10f}"
            ),
            "v2b_distance_angstrom": (
                f"{v2b_contact_distance:.10f}"
            ),
            "distance_change_angstrom": (
                f"{distance_change:.10f}"
            ),
            "v2a_distance_over_vdw_sum": (
                f"{v2a_ratio:.10f}"
            ),
            "v2b_distance_over_vdw_sum": (
                f"{v2b_ratio:.10f}"
            ),
            "compression_improved": compression_improved,
            "covalent_threshold_angstrom": (
                f"{residual_threshold:.10f}"
            ),
            "possible_unintended_covalent_contact": (
                possible_unintended_covalent_contact
            ),
            "hard_contact_below_0p70": (
                v2b_ratio < 0.70
            ),
            "involves_artificial_cap": False,
        }
    ]

    write_csv(
        OUTPUT_DIR / "residual_contact_summary.csv",
        contact_rows,
    )

    write_csv(
        OUTPUT_DIR / "residual_contact_local_bonds.csv",
        local_bond_rows,
    )

    write_csv(
        OUTPUT_DIR / "residual_contact_local_angles.csv",
        local_angle_rows,
    )

    write_csv(
        OUTPUT_DIR / "residual_contact_orientation.csv",
        orientation_rows,
    )

    write_csv(
        OUTPUT_DIR / "residual_contact_local_torsions.csv",
        local_torsion_rows,
    )

    write_csv(
        OUTPUT_DIR / "residual_contact_nearest_nonbonded.csv",
        nearest_rows,
    )

    # This audit does not yet declare the contact acceptable merely
    # because it is above the covalent threshold. It assigns an
    # evidence-based intermediate classification.
    if possible_unintended_covalent_contact:
        classification = (
            "POSSIBLE_UNINTENDED_COVALENT_CONTACT_BLOCKING"
        )
        geometry_acceptance_authorized = False
    elif graph_separation <= 3:
        classification = (
            "LOCAL_INTRAMOLECULAR_PAIR_EXCLUDED_FROM_HARD_GATE"
        )
        geometry_acceptance_authorized = True
    elif (
        graph_separation == 4
        and not possible_unintended_covalent_contact
    ):
        classification = (
            "REAL_1_5_INTRAMOLECULAR_COMPRESSION_"
            "REQUIRES_ELECTRONIC_AND_COMPARATIVE_VALIDATION"
        )
        geometry_acceptance_authorized = True
    else:
        classification = (
            "REAL_LONG_RANGE_COMPRESSION_"
            "REQUIRES_FURTHER_STRUCTURAL_REVIEW"
        )
        geometry_acceptance_authorized = False

    report_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_V2B_RESIDUAL_REAL_CONTACT_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER V2-B Residual Real-Contact Audit — Day028",
                "",
                "## Contact",
                "",
                f"- Bridge atom: `{BRIDGE_ATOM}`",
                f"- Hydrogen: `{HYDROGEN_ATOM}`",
                f"- Hydrogen covalent parent: `{hydrogen_parent}`",
                (
                    f"- Graph separation: "
                    f"**{graph_separation} bonds**"
                ),
                (
                    f"- Shortest path: "
                    f"`{' → '.join(path)}`"
                ),
                "",
                "## Geometry evolution",
                "",
                (
                    f"- V2-A distance: "
                    f"**{v2a_contact_distance:.6f} Å**"
                ),
                (
                    f"- V2-B distance: "
                    f"**{v2b_contact_distance:.6f} Å**"
                ),
                (
                    f"- Distance change: "
                    f"**{distance_change:+.6f} Å**"
                ),
                (
                    f"- V2-A vdW ratio: "
                    f"**{v2a_ratio:.6f}**"
                ),
                (
                    f"- V2-B vdW ratio: "
                    f"**{v2b_ratio:.6f}**"
                ),
                (
                    "- Compression improved during V2-B: "
                    f"**{'YES' if compression_improved else 'NO'}**"
                ),
                "",
                "## Covalent-connectivity test",
                "",
                (
                    f"- Conservative B–H threshold: "
                    f"**{residual_threshold:.3f} Å**"
                ),
                (
                    "- Possible unintended covalent contact: "
                    f"**{'YES' if possible_unintended_covalent_contact else 'NO'}**"
                ),
                "",
                "## Hydrogen orientation",
                "",
                (
                    "- V2-A parent–H versus parent–bridge angle: "
                    f"**{v2a_parent_h_bridge_angle:.6f}°**"
                ),
                (
                    "- V2-B parent–H versus parent–bridge angle: "
                    f"**{v2b_parent_h_bridge_angle:.6f}°**"
                ),
                "",
                "## Classification",
                "",
                f"**{classification}**",
                "",
                "## Interpretation",
                "",
                (
                    "The pair is composed entirely of real R2 atoms "
                    "and is not an artificial-cap contact. Its graph "
                    "separation and conservative covalent-distance "
                    "test must therefore be considered before treating "
                    "the vdW-ratio criterion as a structural failure."
                ),
                "",
                "## Authorization state",
                "",
                (
                    "- LOWER V2-B geometry acceptance retained: "
                    f"**{'YES' if geometry_acceptance_authorized else 'NO'}**"
                ),
                (
                    "- Electronic-property protocol definition: "
                    f"**{'AUTHORIZED' if geometry_acceptance_authorized else 'NOT AUTHORIZED'}**"
                ),
                "- ESP/RESP execution: **NOT AUTHORIZED**",
                "- Force-field parameter adoption: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": classification,
        "bridge_atom": BRIDGE_ATOM,
        "hydrogen_atom": HYDROGEN_ATOM,
        "hydrogen_parent": hydrogen_parent,
        "graph_separation": graph_separation,
        "shortest_path": path,
        "v2a_distance_angstrom": v2a_contact_distance,
        "v2b_distance_angstrom": v2b_contact_distance,
        "distance_change_angstrom": distance_change,
        "v2a_vdw_ratio": v2a_ratio,
        "v2b_vdw_ratio": v2b_ratio,
        "compression_improved": compression_improved,
        "possible_unintended_covalent_contact": (
            possible_unintended_covalent_contact
        ),
        "v2a_parent_H_bridge_angle_deg": (
            v2a_parent_h_bridge_angle
        ),
        "v2b_parent_H_bridge_angle_deg": (
            v2b_parent_h_bridge_angle
        ),
        "geometry_acceptance_retained": (
            geometry_acceptance_authorized
        ),
        "electronic_property_protocol_definition_authorized": (
            geometry_acceptance_authorized
        ),
        "esp_resp_execution_authorized": False,
        "qm_executed_by_this_script": False,
    }

    (
        OUTPUT_DIR
        / "residual_real_contact_audit_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Residual real-contact audit completed.")
    print("Bridge atom:", BRIDGE_ATOM)
    print("Hydrogen:", HYDROGEN_ATOM)
    print("Hydrogen parent:", hydrogen_parent)
    print("Graph separation:", graph_separation)
    print(
        "V2-A distance:",
        f"{v2a_contact_distance:.10f} Å",
    )
    print(
        "V2-B distance:",
        f"{v2b_contact_distance:.10f} Å",
    )
    print(
        "V2-B vdW ratio:",
        f"{v2b_ratio:.10f}",
    )
    print(
        "Possible unintended covalent contact:",
        possible_unintended_covalent_contact,
    )
    print("Decision:", classification)
    print(
        "Electronic-property protocol definition authorized:",
        geometry_acceptance_authorized,
    )
    print("QM executed: False")
    print("Report:", report_path)


if __name__ == "__main__":
    main()
