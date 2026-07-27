#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import json

ROOT = Path(__file__).resolve().parents[2]

CANONICAL_EDGES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

CANONICAL_NODES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

V6B_MAP = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_orca_input/"
    "QM_F06_UPPER_V6B_constraint_map.csv"
)

SELECTION_REPORT = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7_expansion_selection/"
    "QM_F06_UPPER_V7_EXPANSION_SELECTION.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7_boundary_provenance"
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_boundary_cut_provenance.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_BOUNDARY_PROVENANCE.json"
)

REMOVED_CAP = "HCAPV2:UPPER:03"

REQUIRED_ADDED = {
    "A:UPPER:13:1",
    "A:UPPER:14:0",
}


def read_csv(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def boundary_inventory(
    inside_atoms,
    heavy_graph,
):
    cuts = []
    external_degree = defaultdict(int)

    for inside_atom in sorted(inside_atoms):
        for outside_atom in sorted(
            heavy_graph.get(inside_atom, set())
        ):
            if outside_atom in inside_atoms:
                continue

            pair = tuple(sorted((
                inside_atom,
                outside_atom,
            )))

            cuts.append((
                pair,
                inside_atom,
                outside_atom,
            ))

            external_degree[inside_atom] += 1

    unique = {}

    for pair, inside_atom, outside_atom in cuts:
        unique[pair] = (
            inside_atom,
            outside_atom,
        )

    return unique, dict(external_degree)


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    edge_rows = read_csv(CANONICAL_EDGES)
    node_rows = read_csv(CANONICAL_NODES)
    map_rows = read_csv(V6B_MAP)

    selection = json.loads(
        SELECTION_REPORT.read_text(
            encoding="utf-8"
        )
    )

    selected_atoms = set(
        selection[
            "selected_candidate"
        ]["added_atoms"].split("|")
    )

    if selected_atoms != REQUIRED_ADDED:
        raise RuntimeError(
            "Unexpected V7 selected subset: "
            f"{sorted(selected_atoms)}"
        )

    element = {
        row["node_id"]: row["element"]
        for row in node_rows
    }

    node_type = {
        row["node_id"]: row["node_type"]
        for row in node_rows
    }

    heavy_graph = defaultdict(set)
    heavy_edges = set()

    for row in edge_rows:
        if (
            row["source_element"] == "H"
            or row["target_element"] == "H"
        ):
            continue

        first = row["source_node"]
        second = row["target_node"]

        heavy_graph[first].add(second)
        heavy_graph[second].add(first)

        heavy_edges.add(
            tuple(sorted((first, second)))
        )

    v6b_real_atoms = {
        row["atom_id"]
        for row in map_rows
        if (
            row["element"] != "H"
            and row["atom_id"] != REMOVED_CAP
        )
    }

    v7_real_atoms = (
        v6b_real_atoms
        | selected_atoms
    )

    v6b_cuts, v6b_external = (
        boundary_inventory(
            v6b_real_atoms,
            heavy_graph,
        )
    )

    v7_cuts, v7_external = (
        boundary_inventory(
            v7_real_atoms,
            heavy_graph,
        )
    )

    all_pairs = sorted(
        set(v6b_cuts)
        | set(v7_cuts)
    )

    records = []

    for pair in all_pairs:
        in_v6b = pair in v6b_cuts
        in_v7 = pair in v7_cuts

        if in_v7:
            inside_atom, outside_atom = (
                v7_cuts[pair]
            )
        else:
            inside_atom, outside_atom = (
                v6b_cuts[pair]
            )

        if in_v6b and in_v7:
            provenance = "INHERITED_V6B_CUT"
        elif in_v6b and not in_v7:
            provenance = "CLOSED_BY_V7"
        elif not in_v6b and in_v7:
            provenance = "NEW_V7_CUT"
        else:
            raise RuntimeError(
                f"Invalid cut classification: {pair}"
            )

        records.append({
            "first_atom": pair[0],
            "second_atom": pair[1],
            "inside_atom": inside_atom,
            "inside_element": element.get(
                inside_atom,
                "",
            ),
            "inside_node_type": node_type.get(
                inside_atom,
                "",
            ),
            "outside_atom": outside_atom,
            "outside_element": element.get(
                outside_atom,
                "",
            ),
            "outside_node_type": node_type.get(
                outside_atom,
                "",
            ),
            "present_in_v6b": in_v6b,
            "present_in_v7": in_v7,
            "provenance": provenance,
        })

    fieldnames = [
        "first_atom",
        "second_atom",
        "inside_atom",
        "inside_element",
        "inside_node_type",
        "outside_atom",
        "outside_element",
        "outside_node_type",
        "present_in_v6b",
        "present_in_v7",
        "provenance",
    ]

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(records)

    v6b_multicut = {
        atom_id: count
        for atom_id, count
        in v6b_external.items()
        if count > 1
    }

    v7_multicut = {
        atom_id: count
        for atom_id, count
        in v7_external.items()
        if count > 1
    }

    inherited_multicut = {
        atom_id: count
        for atom_id, count
        in v7_multicut.items()
        if atom_id in v6b_multicut
    }

    new_multicut = {
        atom_id: count
        for atom_id, count
        in v7_multicut.items()
        if atom_id not in v6b_multicut
    }

    required_internal_edges = {
        tuple(sorted((
            "A:UPPER:13:1",
            "A:UPPER:14:2",
        ))),
        tuple(sorted((
            "A:UPPER:13:1",
            "A:UPPER:14:0",
        ))),
    }

    required_edges_present = (
        required_internal_edges
        <= heavy_edges
        and all(
            first in v7_real_atoms
            and second in v7_real_atoms
            for first, second
            in required_internal_edges
        )
    )

    removed_cap_absent = all(
        row["atom_id"] != REMOVED_CAP
        for row in map_rows
    )

    a13_external_degree = (
        v7_external.get(
            "A:UPPER:13:1",
            0,
        )
    )

    inherited_A8_4_multicut = (
        "A:UPPER:8:4" in inherited_multicut
        and "A:UPPER:8:4" not in new_multicut
    )

    gates = {
        "selection_matches_expected_subset": (
            selected_atoms == REQUIRED_ADDED
        ),
        "required_canonical_edges_internal": (
            required_edges_present
        ),
        "A13_1_external_degree_at_most_one": (
            a13_external_degree <= 1
        ),
        "no_new_multicut_center": (
            len(new_multicut) == 0
        ),
        "A8_4_multicut_is_inherited": (
            inherited_A8_4_multicut
        ),
        "obsolete_cap_not_in_v6b_real_set": (
            REMOVED_CAP not in v6b_real_atoms
        ),
    }

    passed = all(gates.values())

    report = {
        "model": "QM_F06_UPPER_V7",
        "decision": (
            "QM_F06_UPPER_V7_BOUNDARY_"
            "PROVENANCE_GATE_PASS_"
            "CONSTRUCTION_AUTHORIZED"
            if passed
            else
            "QM_F06_UPPER_V7_BOUNDARY_"
            "PROVENANCE_GATE_FAIL_REVIEW_REQUIRED"
        ),
        "selected_added_atoms": sorted(
            selected_atoms
        ),
        "v6b_real_atom_count": len(
            v6b_real_atoms
        ),
        "v7_real_atom_count": len(
            v7_real_atoms
        ),
        "v6b_boundary_cut_count": len(
            v6b_cuts
        ),
        "v7_boundary_cut_count": len(
            v7_cuts
        ),
        "v6b_multicut_centers": (
            v6b_multicut
        ),
        "v7_multicut_centers": (
            v7_multicut
        ),
        "inherited_multicut_centers": (
            inherited_multicut
        ),
        "new_multicut_centers": (
            new_multicut
        ),
        "A_UPPER_13_1_external_degree": (
            a13_external_degree
        ),
        "gates": gates,
        "construction_authorized": passed,
        "orca_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 108)
    print("QM_F06 UPPER V7 BOUNDARY-CUT PROVENANCE AUDIT")
    print("=" * 108)

    for name, value in gates.items():
        print(
            f"{name:48s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(
        "V6-B real atoms:",
        len(v6b_real_atoms),
    )
    print(
        "V7 real atoms:",
        len(v7_real_atoms),
    )
    print(
        "V6-B boundary cuts:",
        len(v6b_cuts),
    )
    print(
        "V7 boundary cuts:",
        len(v7_cuts),
    )

    print()
    print(
        "V6-B multi-cut centers:",
        v6b_multicut or "none",
    )
    print(
        "V7 multi-cut centers:",
        v7_multicut or "none",
    )
    print(
        "Inherited multi-cut centers:",
        inherited_multicut or "none",
    )
    print(
        "New multi-cut centers:",
        new_multicut or "none",
    )

    print()
    print(
        "A:UPPER:13:1 external degree:",
        a13_external_degree,
    )

    print()
    print("V7 new cuts:")

    new_cut_records = [
        row
        for row in records
        if row["provenance"] == "NEW_V7_CUT"
    ]

    for row in new_cut_records:
        print(
            f"  {row['inside_atom']:28s} -- "
            f"{row['outside_atom']:28s}"
        )

    print()
    print("Cuts closed by V7:")

    closed_records = [
        row
        for row in records
        if row["provenance"] == "CLOSED_BY_V7"
    ]

    for row in closed_records:
        print(
            f"  {row['inside_atom']:28s} -- "
            f"{row['outside_atom']:28s}"
        )

    print()
    print("Decision:", report["decision"])
    print("CSV:", OUTPUT_CSV)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "V7 construction authorized:",
        passed,
    )
    print("ORCA authorized: False")
    print("RESP authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
