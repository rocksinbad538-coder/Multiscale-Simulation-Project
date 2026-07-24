#!/usr/bin/env python3
"""
Define the minimum canonical expansion scope for QM_F06 UPPER V5.

The audit starts from the accepted V4 atom inventory and the canonical
R2 heavy-atom graph. It evaluates the neighborhoods implicated by the
V4 post-QM reorganization:

1. Canonical B-N edge E:2915:
   A:UPPER:13:3 -- A:UPPER:14:4

2. Hydrogen-transfer region:
   P:1641 -- S:1739 and their canonical neighbors

3. Ambiguous cap region:
   S:1738 and BR4:UPPER:00:4

The script:
- resolves graph-column names defensively;
- reconstructs the canonical heavy-atom adjacency;
- identifies V4 real atoms and artificial caps;
- evaluates first- and second-shell expansion candidates;
- counts boundary cuts for each candidate model;
- reports which critical local bonds become internal;
- does not construct V5 or authorize ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CANONICAL_EDGE_TABLE = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

CANONICAL_NODE_TABLE = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

V4_MAP = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction/"
    "QM_F06_UPPER_V4_atom_role_provenance_map.csv"
)

V4_TERMINATION_POINTER = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_executions/"
    "LATEST_V4_EXECUTION.txt"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_expansion_scope"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_EXPANSION_SCOPE.json"
)

NODE_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_candidate_nodes.csv"
)

CUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_candidate_boundary_cuts.csv"
)

MODEL_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_candidate_models.csv"
)

CRITICAL_ATOMS = {
    "A:UPPER:13:3",
    "A:UPPER:14:4",
    "P:1641",
    "S:1739",
    "S:1738",
    "BR4:UPPER:00:4",
}

CRITICAL_EDGES = {
    tuple(sorted((
        "A:UPPER:13:3",
        "A:UPPER:14:4",
    ))): "CANONICAL_EDGE_E2915",

    tuple(sorted((
        "P:1641",
        "S:1739",
    ))): "HYDROGEN_TRANSFER_REGION",

    tuple(sorted((
        "S:1738",
        "BR4:UPPER:00:4",
    ))): "AMBIGUOUS_CAP_REGION",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        headers = reader.fieldnames or []

    return headers, rows


def find_column(
    headers: list[str],
    candidates: tuple[str, ...],
) -> str:
    lowered = {
        header.lower(): header
        for header in headers
    }

    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]

    raise RuntimeError(
        "Could not resolve required column. "
        f"Candidates={candidates}; headers={headers}"
    )


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def canonical_pair(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def graph_shell(
    seeds: set[str],
    adjacency: dict[str, set[str]],
    depth: int,
) -> dict[str, int]:
    distances = {
        seed: 0
        for seed in seeds
    }

    queue = deque(seeds)

    while queue:
        current = queue.popleft()
        current_depth = distances[current]

        if current_depth >= depth:
            continue

        for neighbor in adjacency.get(current, set()):
            if neighbor not in distances:
                distances[neighbor] = current_depth + 1
                queue.append(neighbor)

    return distances


def boundary_cuts(
    model_atoms: set[str],
    adjacency: dict[str, set[str]],
) -> set[tuple[str, str]]:
    cuts = set()

    for atom_id in model_atoms:
        for neighbor in adjacency.get(atom_id, set()):
            if neighbor not in model_atoms:
                cuts.add(
                    canonical_pair(atom_id, neighbor)
                )

    return cuts


def main() -> None:
    for path in (
        CANONICAL_EDGE_TABLE,
        CANONICAL_NODE_TABLE,
        V4_MAP,
        V4_TERMINATION_POINTER,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    edge_headers, edge_rows = read_csv(
        CANONICAL_EDGE_TABLE
    )

    node_headers, node_rows = read_csv(
        CANONICAL_NODE_TABLE
    )

    _, v4_rows = read_csv(V4_MAP)

    edge_first_key = find_column(
        edge_headers,
        (
            "node_a",
            "atom_id_1",
            "source",
            "source_node",
            "atom1",
            "first_atom",
            "first_node",
        ),
    )

    edge_second_key = find_column(
        edge_headers,
        (
            "node_b",
            "atom_id_2",
            "target",
            "target_node",
            "atom2",
            "second_atom",
            "second_node",
        ),
    )

    edge_id_key = None

    for candidate in (
        "edge_id",
        "id",
        "bond_id",
    ):
        if candidate in edge_headers:
            edge_id_key = candidate
            break

    edge_element_first_key = None
    edge_element_second_key = None

    for candidate in (
        "element_a",
        "element_1",
        "first_element",
        "source_element",
    ):
        if candidate in edge_headers:
            edge_element_first_key = candidate
            break

    for candidate in (
        "element_b",
        "element_2",
        "second_element",
        "target_element",
    ):
        if candidate in edge_headers:
            edge_element_second_key = candidate
            break

    node_id_key = find_column(
        node_headers,
        (
            "node_id",
            "atom_id",
            "id",
        ),
    )

    node_element_key = find_column(
        node_headers,
        (
            "element",
            "atom_element",
        ),
    )

    node_type_key = None

    for candidate in (
        "node_type",
        "atom_role",
        "type",
    ):
        if candidate in node_headers:
            node_type_key = candidate
            break

    node_metadata = {}

    for row in node_rows:
        atom_id = row[node_id_key]
        node_metadata[atom_id] = {
            "element": row[node_element_key],
            "node_type": (
                row[node_type_key]
                if node_type_key
                else ""
            ),
        }

    adjacency = defaultdict(set)
    edge_metadata = {}

    for row in edge_rows:
        first = row[edge_first_key]
        second = row[edge_second_key]

        if not first or not second:
            continue

        first_element = (
            row[edge_element_first_key]
            if edge_element_first_key
            else node_metadata.get(
                first,
                {},
            ).get("element", "")
        )

        second_element = (
            row[edge_element_second_key]
            if edge_element_second_key
            else node_metadata.get(
                second,
                {},
            ).get("element", "")
        )

        # V5 scope is defined from the canonical heavy-atom graph.
        if (
            first_element == "H"
            or second_element == "H"
        ):
            continue

        pair = canonical_pair(first, second)

        adjacency[first].add(second)
        adjacency[second].add(first)

        edge_metadata[pair] = {
            "edge_id": (
                row[edge_id_key]
                if edge_id_key
                else ""
            ),
            "first": first,
            "second": second,
            "first_element": first_element,
            "second_element": second_element,
        }

    missing_critical = (
        CRITICAL_ATOMS
        - set(adjacency)
    )

    if missing_critical:
        raise RuntimeError(
            "Critical atoms missing from canonical graph: "
            f"{sorted(missing_critical)}"
        )

    v4_real_atoms = {
        row["atom_id"]
        for row in v4_rows
        if not parse_bool(row["artificial_cap"])
    }

    v4_caps = {
        row["atom_id"]
        for row in v4_rows
        if parse_bool(row["artificial_cap"])
    }

    critical_existing = (
        CRITICAL_ATOMS
        & v4_real_atoms
    )

    if critical_existing != CRITICAL_ATOMS:
        raise RuntimeError(
            "Not all critical atoms are present as real V4 atoms: "
            f"{sorted(CRITICAL_ATOMS - critical_existing)}"
        )

    shell_1 = graph_shell(
        CRITICAL_ATOMS,
        adjacency,
        depth=1,
    )

    shell_2 = graph_shell(
        CRITICAL_ATOMS,
        adjacency,
        depth=2,
    )

    first_shell_new = {
        atom_id
        for atom_id, depth in shell_1.items()
        if depth == 1
        and atom_id not in v4_real_atoms
    }

    second_shell_new = {
        atom_id
        for atom_id, depth in shell_2.items()
        if depth == 2
        and atom_id not in v4_real_atoms
    }

    candidate_models = {
        "V4_REAL_BASELINE": set(v4_real_atoms),

        "V5_CRITICAL_FIRST_SHELL": (
            set(v4_real_atoms)
            | first_shell_new
        ),

        "V5_CRITICAL_SECOND_SHELL": (
            set(v4_real_atoms)
            | first_shell_new
            | second_shell_new
        ),
    }

    model_records = []
    cut_records = []

    baseline_cuts = boundary_cuts(
        candidate_models["V4_REAL_BASELINE"],
        adjacency,
    )

    for model_name, model_atoms in candidate_models.items():
        cuts = boundary_cuts(
            model_atoms,
            adjacency,
        )

        critical_internal = {}

        for pair, label in CRITICAL_EDGES.items():
            critical_internal[label] = (
                pair[0] in model_atoms
                and pair[1] in model_atoms
                and pair in edge_metadata
            )

        newly_added = (
            model_atoms
            - v4_real_atoms
        )

        new_cuts_relative_to_v4 = (
            cuts - baseline_cuts
        )

        cuts_closed_relative_to_v4 = (
            baseline_cuts - cuts
        )

        model_records.append({
            "model": model_name,
            "real_atom_count": len(model_atoms),
            "new_real_atoms": len(newly_added),
            "boundary_cut_count": len(cuts),
            "new_cuts_relative_to_v4": len(
                new_cuts_relative_to_v4
            ),
            "cuts_closed_relative_to_v4": len(
                cuts_closed_relative_to_v4
            ),
            "E2915_internal": critical_internal[
                "CANONICAL_EDGE_E2915"
            ],
            "P1641_S1739_internal": critical_internal[
                "HYDROGEN_TRANSFER_REGION"
            ],
            "S1738_BR4_00_4_internal": critical_internal[
                "AMBIGUOUS_CAP_REGION"
            ],
            "new_atom_ids": "|".join(
                sorted(newly_added)
            ),
        })

        for first, second in sorted(cuts):
            metadata = edge_metadata.get(
                canonical_pair(first, second),
                {},
            )

            inside = (
                first
                if first in model_atoms
                else second
            )

            outside = (
                second
                if first in model_atoms
                else first
            )

            cut_records.append({
                "model": model_name,
                "edge_id": metadata.get(
                    "edge_id",
                    "",
                ),
                "inside_atom": inside,
                "inside_element": node_metadata.get(
                    inside,
                    {},
                ).get("element", ""),
                "outside_atom": outside,
                "outside_element": node_metadata.get(
                    outside,
                    {},
                ).get("element", ""),
                "preexisting_v4_cut": (
                    canonical_pair(first, second)
                    in baseline_cuts
                ),
                "critical_region": (
                    first in CRITICAL_ATOMS
                    or second in CRITICAL_ATOMS
                ),
            })

    node_records = []

    all_candidate_atoms = (
        set(v4_real_atoms)
        | first_shell_new
        | second_shell_new
    )

    for atom_id in sorted(all_candidate_atoms):
        metadata = node_metadata.get(
            atom_id,
            {},
        )

        node_records.append({
            "atom_id": atom_id,
            "element": metadata.get(
                "element",
                "",
            ),
            "node_type": metadata.get(
                "node_type",
                "",
            ),
            "in_v4_real_model": (
                atom_id in v4_real_atoms
            ),
            "critical_atom": (
                atom_id in CRITICAL_ATOMS
            ),
            "first_shell_new": (
                atom_id in first_shell_new
            ),
            "second_shell_new": (
                atom_id in second_shell_new
            ),
            "canonical_heavy_degree": len(
                adjacency.get(atom_id, set())
            ),
            "canonical_heavy_neighbors": "|".join(
                sorted(adjacency.get(atom_id, set()))
            ),
        })

    with NODE_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(node_records[0]),
        )
        writer.writeheader()
        writer.writerows(node_records)

    with CUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(cut_records[0]),
        )
        writer.writeheader()
        writer.writerows(cut_records)

    with MODEL_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(model_records[0]),
        )
        writer.writeheader()
        writer.writerows(model_records)

    latest_execution = Path(
        V4_TERMINATION_POINTER.read_text(
            encoding="utf-8"
        ).strip()
    )

    termination_record = (
        ROOT
        / latest_execution
        / "QM_F06_UPPER_V4_TERMINATION_RECORD.json"
    )

    require_file(termination_record)

    report = {
        "decision": (
            "QM_F06_UPPER_V5_EXPANSION_SCOPE_AUDITED_"
            "MODEL_SELECTION_PENDING"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "v4_failure_basis": {
            "termination_record": str(
                termination_record.relative_to(ROOT)
            ),
            "canonical_bond_rupture": (
                "A:UPPER:13:3--A:UPPER:14:4"
            ),
            "hydrogen_transfer_region": (
                "S:1739--H4:UPPER:0203:0--P:1641"
            ),
            "ambiguous_cap_region": (
                "S:1738--HCAPV4:UPPER:S1738:BR4_14_1--"
                "BR4:UPPER:00:4"
            ),
        },
        "critical_atoms": sorted(CRITICAL_ATOMS),
        "v4": {
            "real_atom_count": len(v4_real_atoms),
            "artificial_cap_count": len(v4_caps),
            "baseline_heavy_boundary_cuts": len(
                baseline_cuts
            ),
        },
        "first_shell_new_atoms": sorted(
            first_shell_new
        ),
        "second_shell_new_atoms": sorted(
            second_shell_new
        ),
        "candidate_models": model_records,
        "files": {
            "candidate_nodes": str(
                NODE_CSV.relative_to(ROOT)
            ),
            "candidate_boundary_cuts": str(
                CUT_CSV.relative_to(ROOT)
            ),
            "candidate_models": str(
                MODEL_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "canonical_edge_table": sha256(
                CANONICAL_EDGE_TABLE
            ),
            "canonical_node_table": sha256(
                CANONICAL_NODE_TABLE
            ),
            "v4_map": sha256(V4_MAP),
            "v4_termination_record": sha256(
                termination_record
            ),
            "candidate_nodes": sha256(NODE_CSV),
            "candidate_boundary_cuts": sha256(
                CUT_CSV
            ),
            "candidate_models": sha256(
                MODEL_CSV
            ),
        },
        "authorization": {
            "v5_model_selected": False,
            "v5_geometry_construction_authorized": False,
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 90)
    print("QM_F06 UPPER V5 EXPANSION-SCOPE AUDIT")
    print("=" * 90)

    print(
        "V4 real atoms:",
        len(v4_real_atoms),
    )
    print(
        "V4 artificial caps:",
        len(v4_caps),
    )
    print(
        "V4 heavy boundary cuts:",
        len(baseline_cuts),
    )
    print()

    print(
        "First-shell new atoms:",
        sorted(first_shell_new),
    )
    print(
        "Second-shell new atoms:",
        sorted(second_shell_new),
    )
    print()

    print("CANDIDATE MODELS")

    for row in model_records:
        print(
            f"{row['model']:30s} "
            f"real={row['real_atom_count']:3d} "
            f"new={row['new_real_atoms']:3d} "
            f"cuts={row['boundary_cut_count']:3d} "
            f"new_cuts={row['new_cuts_relative_to_v4']:3d} "
            f"closed={row['cuts_closed_relative_to_v4']:3d}"
        )

    print()
    print(
        "Decision:",
        report["decision"],
    )
    print("Report:", REPORT_PATH)
    print("Nodes:", NODE_CSV)
    print("Cuts:", CUT_CSV)
    print("Models:", MODEL_CSV)
    print()
    print("V5 model selected: False")
    print("V5 construction authorized: False")
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
