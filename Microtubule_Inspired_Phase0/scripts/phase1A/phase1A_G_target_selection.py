#!/usr/bin/env python3
"""
DAY040 / D040-A2

Authoritative force-field target selection for Phase 1A-G.

The block audits the accepted and validated GROMACS topology family,
its include graph, molecule declarations, coordinate candidates, atom
sections and possible charge-bearing targets.

No topology is modified.
No charges are written into a force field.
No MD calculation is executed.
"""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

PREFLIGHT_JSON = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_force_field_preflight"
    / "QM_F06_UPPER_V7A_R1_PHASE1A_G_FORCE_FIELD_PREFLIGHT.json"
)

CLOSURE_POINTER = (
    ROOT
    / "runs/phase1A"
    / "LATEST_PHASE1A_F_CHARGE_MODEL_CLOSURE.txt"
)

EXPECTED_PREFLIGHT_DECISION = (
    "D040_A1_FORCE_FIELD_INTEGRATION_PREFLIGHT_PASS_"
    "TARGET_SELECTION_REVIEW_AUTHORIZED"
)

EXPECTED_CLOSURE_DECISION = (
    "D039_A16_PHASE1A_F_CLOSED_"
    "WORKING_CHARGE_MODEL_ADOPTED_FOR_"
    "FORCE_FIELD_INTEGRATION"
)

TARGET_DIR = (
    ROOT
    / "parameters/phase1A/accepted"
    / "hybrid_hbnBonded_kang2000_improperGeo100_validated"
)

TARGET_TOP = (
    TARGET_DIR
    / "hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_"
    "hbnBonded_kang2000_improperGeo100.top"
)

TARGET_HBN_ITP = (
    TARGET_DIR
    / "hbn_bonded_candidate_kang2000_improperGeo100.itp"
)

TARGET_PYRENE_ITP = TARGET_DIR / "pyrene.itp"
TARGET_WATER_ITP = TARGET_DIR / "tip4p2005.itp"

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_target_selection"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_G_TARGET_SELECTION.json"
)

INCLUDE_GRAPH_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GROMACS_INCLUDE_GRAPH.csv"
)

MOLECULE_TYPES_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GROMACS_MOLECULE_TYPES.csv"
)

ATOM_SECTION_SUMMARY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GROMACS_ATOM_SECTION_SUMMARY.csv"
)


def read_lines(path: Path) -> list[str]:
    require_file(path)
    return path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()


def strip_gromacs_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse_includes(path: Path) -> list[str]:
    includes = []

    pattern = re.compile(
        r'^\s*#include\s+["<]([^">]+)[">]'
    )

    for line in read_lines(path):
        match = pattern.match(line)

        if match:
            includes.append(match.group(1))

    return includes


def resolve_include(
    parent: Path,
    include_value: str,
) -> tuple[Path | None, str]:
    local_candidate = (
        parent.parent
        / include_value
    ).resolve()

    if local_candidate.is_file():
        return local_candidate, "LOCAL_RELATIVE"

    root_candidate = (
        ROOT
        / include_value
    ).resolve()

    if root_candidate.is_file():
        return root_candidate, "PROJECT_ROOT_RELATIVE"

    basename_matches = [
        path.resolve()
        for path in ROOT.rglob(
            Path(include_value).name
        )
        if path.is_file()
    ]

    unique_matches = sorted(
        set(basename_matches),
        key=str,
    )

    if len(unique_matches) == 1:
        return unique_matches[0], "UNIQUE_BASENAME_MATCH"

    if len(unique_matches) > 1:
        return None, (
            "AMBIGUOUS_BASENAME_MATCH:"
            + str(len(unique_matches))
        )

    return None, "UNRESOLVED"


def parse_sections(path: Path) -> dict[str, list[list[str]]]:
    sections: dict[str, list[list[str]]] = {}
    active_section: str | None = None

    for raw_line in read_lines(path):
        stripped = strip_gromacs_comment(
            raw_line
        )

        if not stripped:
            continue

        section_match = re.match(
            r"^\[\s*([^\]]+?)\s*\]$",
            stripped,
        )

        if section_match:
            active_section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            sections.setdefault(
                active_section,
                [],
            )

            continue

        if stripped.startswith("#"):
            continue

        if active_section is not None:
            sections[
                active_section
            ].append(
                stripped.split()
            )

    return sections


def parse_molecule_types(
    path: Path,
) -> list[dict]:
    sections = parse_sections(path)
    rows = []

    for tokens in sections.get(
        "moleculetype",
        [],
    ):
        if len(tokens) >= 2:
            rows.append(
                {
                    "source_path": str(
                        path.relative_to(ROOT)
                    ),
                    "molecule_type": tokens[0],
                    "nrexcl": tokens[1],
                }
            )

    return rows


def summarize_atoms(
    path: Path,
) -> dict:
    sections = parse_sections(path)
    atom_rows = sections.get(
        "atoms",
        [],
    )

    parsed = []

    for tokens in atom_rows:
        if len(tokens) < 7:
            continue

        try:
            atom_index = int(tokens[0])
            charge = float(tokens[6])
        except ValueError:
            continue

        parsed.append(
            {
                "atom_index": atom_index,
                "atom_type": tokens[1],
                "residue_number": tokens[2],
                "residue_name": tokens[3],
                "atom_name": tokens[4],
                "charge_group": tokens[5],
                "charge_e": charge,
                "mass": (
                    float(tokens[7])
                    if len(tokens) >= 8
                    else None
                ),
            }
        )

    return {
        "source_path": str(
            path.relative_to(ROOT)
        ),
        "atom_count": len(parsed),
        "charge_sum_e": sum(
            row["charge_e"]
            for row in parsed
        ),
        "atom_types": sorted(
            {
                row["atom_type"]
                for row in parsed
            }
        ),
        "residue_names": sorted(
            {
                row["residue_name"]
                for row in parsed
            }
        ),
        "atom_names": [
            row["atom_name"]
            for row in parsed
        ],
        "atoms": parsed,
        "has_atoms_section": bool(
            atom_rows
        ),
    }


def candidate_coordinate_files() -> list[Path]:
    candidates = []

    priority_roots = (
        ROOT / "parameters/phase1A/accepted",
        ROOT / "structures/phase1A/accepted",
        ROOT / "runs/phase1A/day021_accepted_hydrated_topology_audit",
    )

    for search_root in priority_roots:
        if not search_root.exists():
            continue

        for suffix in (
            "*.gro",
            "*.pdb",
        ):
            candidates.extend(
                path
                for path in search_root.rglob(
                    suffix
                )
                if path.is_file()
            )

    return sorted(
        set(candidates),
        key=lambda path: str(
            path.relative_to(ROOT)
        ),
    )


print("=" * 100)
print("DAY040 / D040-A2 — AUTHORITATIVE FORCE-FIELD TARGET SELECTION")
print("=" * 100)


print("\n[1] UPSTREAM AUTHORIZATION")

require_file(PREFLIGHT_JSON)
require_file(CLOSURE_POINTER)

preflight = load_json(
    PREFLIGHT_JSON
)

if (
    preflight.get("decision")
    != EXPECTED_PREFLIGHT_DECISION
):
    raise RuntimeError(
        "Unexpected D040-A1 decision.\n"
        f"Observed: {preflight.get('decision')}"
    )

if (
    preflight.get(
        "authorizations",
        {},
    ).get(
        "force_field_target_selection_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Force-field target selection is not authorized"
    )

closure_dir = (
    ROOT
    / CLOSURE_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

closure_json = (
    closure_dir
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_CLOSURE.json"
)

charges_csv = (
    closure_dir
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.csv"
)

require_file(closure_json)
require_file(charges_csv)

closure = load_json(
    closure_json
)

if (
    closure.get("decision")
    != EXPECTED_CLOSURE_DECISION
):
    raise RuntimeError(
        "Unexpected Phase 1A-F closure decision"
    )

print("D040_A1_decision_gate = PASS")
print("Phase1A_F_closure_gate = PASS")
print("target_selection_review_gate = PASS")
print("topology_modification_blocked_gate = PASS")


print("\n[2] AUTHORITATIVE TARGET FAMILY")

required_target_files = (
    TARGET_TOP,
    TARGET_HBN_ITP,
    TARGET_PYRENE_ITP,
    TARGET_WATER_ITP,
)

for path in required_target_files:
    require_file(path)
    print(
        f"FOUND bytes={path.stat().st_size:8d} "
        f"sha256={sha256(path)} "
        f"{path}"
    )

print(
    "selected_engine = GROMACS"
)
print(
    "selected_force_field_family = "
    "HYBRID_HBN_BONDED_KANG2000_IMPROPER_GEO100"
)
print(
    "selected_status = ACCEPTED_AND_VALIDATED"
)
print(
    f"selected_topology = "
    f"{TARGET_TOP.relative_to(ROOT)}"
)
print(
    f"selected_hbn_component = "
    f"{TARGET_HBN_ITP.relative_to(ROOT)}"
)


print("\n[3] INCLUDE GRAPH")

queue = [TARGET_TOP]
visited: set[Path] = set()
include_records = []
all_resolved_files: set[Path] = set()

while queue:
    parent = queue.pop(0).resolve()

    if parent in visited:
        continue

    visited.add(parent)
    all_resolved_files.add(parent)

    for include_value in parse_includes(
        parent
    ):
        resolved, resolution = (
            resolve_include(
                parent,
                include_value,
            )
        )

        include_record = {
            "parent": str(
                parent.relative_to(ROOT)
            ),
            "include_value": include_value,
            "resolution": resolution,
            "resolved_path": (
                str(
                    resolved.relative_to(ROOT)
                )
                if resolved is not None
                and ROOT in resolved.parents
                else (
                    str(resolved)
                    if resolved is not None
                    else ""
                )
            ),
            "resolved_exists": (
                resolved is not None
                and resolved.is_file()
            ),
        }

        include_records.append(
            include_record
        )

        if (
            resolved is not None
            and resolved.is_file()
            and ROOT in resolved.parents
            and resolved.suffix.lower()
            in {
                ".top",
                ".itp",
            }
        ):
            queue.append(resolved)

unresolved_records = [
    row
    for row in include_records
    if not row["resolved_exists"]
]

print(
    f"include_record_count = "
    f"{len(include_records)}"
)
print(
    f"resolved_project_file_count = "
    f"{len(all_resolved_files)}"
)
print(
    f"unresolved_include_count = "
    f"{len(unresolved_records)}"
)

for row in include_records:
    print(
        f"parent={row['parent']} "
        f"include={row['include_value']} "
        f"resolution={row['resolution']} "
        f"resolved={row['resolved_path']}"
    )


print("\n[4] MOLECULE TYPES")

molecule_type_rows = []

for path in sorted(
    all_resolved_files,
    key=str,
):
    if path.suffix.lower() not in {
        ".top",
        ".itp",
    }:
        continue

    molecule_type_rows.extend(
        parse_molecule_types(path)
    )

print(
    f"molecule_type_record_count = "
    f"{len(molecule_type_rows)}"
)

for row in molecule_type_rows:
    print(
        f"molecule_type={row['molecule_type']} "
        f"nrexcl={row['nrexcl']} "
        f"source={row['source_path']}"
    )


print("\n[5] ATOM-SECTION AUDIT")

atom_summaries = []

for path in sorted(
    all_resolved_files
    | {
        TARGET_HBN_ITP.resolve(),
        TARGET_PYRENE_ITP.resolve(),
        TARGET_WATER_ITP.resolve(),
    },
    key=str,
):
    if path.suffix.lower() not in {
        ".top",
        ".itp",
    }:
        continue

    summary = summarize_atoms(path)

    if summary[
        "has_atoms_section"
    ]:
        atom_summaries.append(
            summary
        )

        print(
            f"path={summary['source_path']} "
            f"atoms={summary['atom_count']} "
            f"charge_sum_e={summary['charge_sum_e']:.12g} "
            f"residues={summary['residue_names']} "
            f"atom_types={summary['atom_types']}"
        )


print("\n[6] ADOPTED CHARGE CONTRACT")

with charges_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    charge_rows = list(
        csv.DictReader(handle)
    )

charge_rows.sort(
    key=lambda row: int(
        row[
            "real_atom_sequence_index"
        ]
    )
)

if len(charge_rows) != 37:
    raise RuntimeError(
        "Expected 37 adopted working charges"
    )

charge_atom_ids = [
    row["atom_id"]
    for row in charge_rows
]

charge_elements = [
    row["element"]
    for row in charge_rows
]

charge_original_indices = [
    int(
        row[
            "original_atom_index_0based"
        ]
    )
    for row in charge_rows
]

charge_sum_e = sum(
    float(
        row[
            "adopted_working_charge_e"
        ]
    )
    for row in charge_rows
)

print(
    f"adopted_charge_count = "
    f"{len(charge_rows)}"
)
print(
    f"adopted_charge_sum_e = "
    f"{charge_sum_e:.16g}"
)
print(
    f"adopted_composition = "
    f"B:{charge_elements.count('B')} "
    f"N:{charge_elements.count('N')} "
    f"H:{charge_elements.count('H')}"
)
print(
    f"original_index_min = "
    f"{min(charge_original_indices)}"
)
print(
    f"original_index_max = "
    f"{max(charge_original_indices)}"
)
print(
    f"unique_atom_id_count = "
    f"{len(set(charge_atom_ids))}"
)


print("\n[7] CANDIDATE COORDINATE FILES")

coordinate_candidates = (
    candidate_coordinate_files()
)

print(
    f"coordinate_candidate_count = "
    f"{len(coordinate_candidates)}"
)

for path in coordinate_candidates:
    name_lower = path.name.lower()

    score = 0
    reasons = []

    if "accepted" in str(
        path.relative_to(ROOT)
    ).lower():
        score += 2
        reasons.append("ACCEPTED_PATH")

    if "hydr" in name_lower:
        score += 2
        reasons.append("HYDRATED_NAME")

    if "solvat" in name_lower:
        score += 3
        reasons.append("SOLVATED_NAME")

    if "gap45" in name_lower:
        score += 2
        reasons.append("GAP45_NAME")

    if "pyr5shift" in name_lower:
        score += 2
        reasons.append("PYR5SHIFT_NAME")

    if score > 0:
        print(
            f"score={score:2d} "
            f"bytes={path.stat().st_size:8d} "
            f"path={path.relative_to(ROOT)} "
            f"reasons={reasons}"
        )


print("\n[8] TARGET-SELECTION INTERPRETATION")

hbn_summary_matches = [
    summary
    for summary in atom_summaries
    if (
        Path(
            summary["source_path"]
        ).name
        == TARGET_HBN_ITP.name
    )
]

if len(hbn_summary_matches) != 1:
    raise RuntimeError(
        "Could not identify the HBN atom section uniquely"
    )

hbn_summary = hbn_summary_matches[0]

selected_component_has_B = any(
    atom_type.lower().startswith("b")
    or atom_type.lower() == "b"
    for atom_type in hbn_summary[
        "atom_types"
    ]
)

selected_component_has_N = any(
    atom_type.lower().startswith("n")
    or atom_type.lower() == "n"
    for atom_type in hbn_summary[
        "atom_types"
    ]
)

selected_component_atom_count = (
    hbn_summary["atom_count"]
)

print(
    "authoritative_topology_family_gate = PASS"
)
print(
    f"hbn_component_atom_count = "
    f"{selected_component_atom_count}"
)
print(
    f"hbn_component_charge_sum_e = "
    f"{hbn_summary['charge_sum_e']:.16g}"
)
print(
    f"hbn_component_has_B_type = "
    f"{selected_component_has_B}"
)
print(
    f"hbn_component_has_N_type = "
    f"{selected_component_has_N}"
)
print(
    "integration_target_interpretation = "
    "HBN_BONDED_COMPONENT_REQUIRES_LOCAL_37_ATOM_MAPPING"
)
print(
    "direct_global_charge_replacement_authorized = False"
)


print("\n[9] WRITE OUTPUTS")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

with INCLUDE_GRAPH_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    fieldnames = [
        "parent",
        "include_value",
        "resolution",
        "resolved_path",
        "resolved_exists",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()
    writer.writerows(
        include_records
    )

with MOLECULE_TYPES_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    fieldnames = [
        "source_path",
        "molecule_type",
        "nrexcl",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()
    writer.writerows(
        molecule_type_rows
    )

with ATOM_SECTION_SUMMARY_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    fieldnames = [
        "source_path",
        "atom_count",
        "charge_sum_e",
        "residue_names",
        "atom_types",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for summary in atom_summaries:
        writer.writerow(
            {
                "source_path": (
                    summary[
                        "source_path"
                    ]
                ),
                "atom_count": (
                    summary[
                        "atom_count"
                    ]
                ),
                "charge_sum_e": (
                    summary[
                        "charge_sum_e"
                    ]
                ),
                "residue_names": json.dumps(
                    summary[
                        "residue_names"
                    ]
                ),
                "atom_types": json.dumps(
                    summary[
                        "atom_types"
                    ]
                ),
            }
        )


print("\n[10] GATES")

gates = {
    "D040_A1_decision_gate": True,
    "Phase1A_F_closure_gate": True,
    "selected_engine_GROMACS_gate": True,
    "selected_target_directory_gate": (
        TARGET_DIR.is_dir()
    ),
    "selected_topology_gate": (
        TARGET_TOP.is_file()
    ),
    "selected_HBN_component_gate": (
        TARGET_HBN_ITP.is_file()
    ),
    "selected_pyrene_component_gate": (
        TARGET_PYRENE_ITP.is_file()
    ),
    "selected_water_component_gate": (
        TARGET_WATER_ITP.is_file()
    ),
    "include_graph_created_gate": (
        INCLUDE_GRAPH_CSV.is_file()
        and INCLUDE_GRAPH_CSV.stat().st_size > 0
    ),
    "molecule_type_inventory_created_gate": (
        MOLECULE_TYPES_CSV.is_file()
        and MOLECULE_TYPES_CSV.stat().st_size > 0
    ),
    "atom_section_inventory_created_gate": (
        ATOM_SECTION_SUMMARY_CSV.is_file()
        and ATOM_SECTION_SUMMARY_CSV.stat().st_size > 0
    ),
    "adopted_charge_count_37_gate": (
        len(charge_rows) == 37
    ),
    "adopted_charge_neutrality_gate": (
        abs(charge_sum_e) <= 1.0e-10
    ),
    "selected_component_has_atoms_gate": (
        selected_component_atom_count > 0
    ),
    "selected_component_has_B_gate": (
        selected_component_has_B
    ),
    "selected_component_has_N_gate": (
        selected_component_has_N
    ),
    "no_topology_modified_gate": True,
    "no_charge_mapping_executed_gate": True,
    "force_field_adoption_blocked_gate": True,
    "MD_execution_blocked_gate": True,
}

for gate_name, value in gates.items():
    print(
        f"{gate_name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[11] WRITE REPORT")

decision = (
    "D040_A2_FORCE_FIELD_TARGET_SELECTION_PASS_"
    "LOCAL_37_ATOM_MAPPING_DESIGN_AUTHORIZED"
    if all_gates_pass
    else
    "D040_A2_FORCE_FIELD_TARGET_SELECTION_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "selected_target": {
        "engine": "GROMACS",
        "force_field_family": (
            "HYBRID_HBN_BONDED_KANG2000_IMPROPER_GEO100"
        ),
        "status": (
            "ACCEPTED_AND_VALIDATED"
        ),
        "directory": str(
            TARGET_DIR.relative_to(ROOT)
        ),
        "topology": {
            "path": str(
                TARGET_TOP.relative_to(ROOT)
            ),
            "sha256": sha256(TARGET_TOP),
        },
        "hbn_component": {
            "path": str(
                TARGET_HBN_ITP.relative_to(ROOT)
            ),
            "sha256": sha256(
                TARGET_HBN_ITP
            ),
            "atom_count": (
                selected_component_atom_count
            ),
            "charge_sum_e": (
                hbn_summary[
                    "charge_sum_e"
                ]
            ),
            "atom_types": (
                hbn_summary[
                    "atom_types"
                ]
            ),
            "residue_names": (
                hbn_summary[
                    "residue_names"
                ]
            ),
        },
        "pyrene_component": {
            "path": str(
                TARGET_PYRENE_ITP.relative_to(ROOT)
            ),
            "sha256": sha256(
                TARGET_PYRENE_ITP
            ),
        },
        "water_component": {
            "path": str(
                TARGET_WATER_ITP.relative_to(ROOT)
            ),
            "sha256": sha256(
                TARGET_WATER_ITP
            ),
        },
    },
    "adopted_charge_contract": {
        "path": str(
            charges_csv.relative_to(ROOT)
        ),
        "sha256": sha256(
            charges_csv
        ),
        "atom_count": len(
            charge_rows
        ),
        "charge_sum_e": (
            charge_sum_e
        ),
        "atom_ids": charge_atom_ids,
        "original_indices_0based": (
            charge_original_indices
        ),
    },
    "include_graph": (
        include_records
    ),
    "molecule_types": (
        molecule_type_rows
    ),
    "atom_section_summaries": [
        {
            key: value
            for key, value in summary.items()
            if key != "atoms"
        }
        for summary in atom_summaries
    ],
    "coordinate_candidates": [
        {
            "path": str(
                path.relative_to(ROOT)
            ),
            "bytes": int(
                path.stat().st_size
            ),
            "sha256": sha256(path),
        }
        for path in coordinate_candidates
    ],
    "gates": gates,
    "authorizations": {
        "authoritative_force_field_target_selected": (
            all_gates_pass
        ),
        "local_37_atom_mapping_design_authorized": (
            all_gates_pass
        ),
        "charge_to_topology_mapping_execution_authorized": False,
        "topology_modification_authorized": False,
        "force_field_adoption_authorized": False,
        "energy_execution_authorized": False,
        "minimization_execution_authorized": False,
        "validation_MD_execution_authorized": False,
        "production_MD_authorized": False,
    },
    "next_required_block": {
        "name": (
            "D040_A3_LOCAL_37_ATOM_MAPPING_DESIGN"
        ),
        "required_actions": [
            (
                "Identify the exact 37 HBN topology atoms "
                "corresponding to the adopted QM real atoms."
            ),
            (
                "Use atom IDs, original indices, coordinates "
                "and local connectivity; do not map by element alone."
            ),
            (
                "Produce a one-to-one mapping table with "
                "distance and topology-consistency diagnostics."
            ),
            (
                "Keep topology modification blocked until "
                "the mapping is complete and reviewed."
            ),
        ],
    },
}

REPORT_JSON.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"report_path = {REPORT_JSON}")
print(
    f"report_sha256 = "
    f"{sha256(REPORT_JSON)}"
)


print("\n[12] DECISION")

print(f"decision={decision}")
print(
    "authoritative_force_field_target_selected="
    f"{all_gates_pass}"
)
print(
    "local_37_atom_mapping_design_authorized="
    f"{all_gates_pass}"
)
print(
    "charge_to_topology_mapping_execution_authorized=False"
)
print(
    "topology_modification_authorized=False"
)
print(
    "force_field_adoption_authorized=False"
)
print(
    "validation_MD_execution_authorized=False"
)
print(
    "production_MD_authorized=False"
)
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
