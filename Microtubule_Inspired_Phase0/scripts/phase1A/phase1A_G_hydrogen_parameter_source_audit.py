#!/usr/bin/env python3
"""
DAY040 / D040-A5

Hydrogen parameter-source audit for Phase 1A-G.

The block searches the project repository for existing GROMACS
definitions or explicit interactions relevant to edge-passivated hBN:

- hydrogen atom types;
- B-H and N-H bonds;
- B-B-H, N-B-H, B-N-H and N-N-H angles;
- H-containing dihedrals and impropers;
- parameter provenance in accepted, candidate and historical files.

No parameter is adopted.
No accepted topology is modified.
No coordinates are modified.
No GROMACS calculation is executed.
"""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

A4_REPORT = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_hydrogen_augmentation_design"
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_AUGMENTATION_DESIGN.json"
)

A4_PARENT_MAP = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_hydrogen_augmentation_design"
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARENT_MAPPING.csv"
)

A4_REQUIREMENTS = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_hydrogen_augmentation_design"
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_REQUIREMENTS.csv"
)

EXPECTED_A4_DECISION = (
    "D040_A4_HYDROGEN_AUGMENTATION_DESIGN_PASS_"
    "PARAMETER_SOURCE_AUDIT_AUTHORIZED"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_hydrogen_parameter_source_audit"
)

ATOMTYPE_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_ATOMTYPE_CANDIDATES.csv"
)

BONDED_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_BONDED_PARAMETER_CANDIDATES.csv"
)

TEXT_HITS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_TEXT_HITS.csv"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_SOURCE_AUDIT.json"
)

SEARCH_EXTENSIONS = {
    ".itp",
    ".top",
    ".rtp",
    ".prm",
    ".str",
    ".ff",
    ".txt",
    ".md",
    ".json",
    ".csv",
}

EXCLUDED_DIR_NAMES = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "node_modules",
}

GROMACS_PARAMETER_SECTIONS = {
    "atomtypes",
    "bondtypes",
    "angletypes",
    "dihedraltypes",
    "nonbond_params",
}

EXPLICIT_INTERACTION_SECTIONS = {
    "atoms",
    "bonds",
    "angles",
    "dihedrals",
}

HYDROGEN_TYPE_PATTERNS = (
    re.compile(r"^H", re.IGNORECASE),
    re.compile(r".*H$", re.IGNORECASE),
)

TEXT_PATTERNS = {
    "B_H": re.compile(
        r"(?i)(\bB[\w+-]*\s*[-–]\s*H[\w+-]*\b|\bB-H\b|\bB H\b)"
    ),
    "N_H": re.compile(
        r"(?i)(\bN[\w+-]*\s*[-–]\s*H[\w+-]*\b|\bN-H\b|\bN H\b)"
    ),
    "HBN_PASSIVATION": re.compile(
        r"(?i)(hbn.*passivat|passivat.*hbn|edge.*hydrogen|hydrogen.*edge)"
    ),
    "HYDROGEN_ATOMTYPE": re.compile(
        r"(?i)(hydrogen atom type|atomtype.*hydrogen|hydrogen.*atomtype)"
    ),
}


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def is_excluded(path: Path) -> bool:
    return any(
        part in EXCLUDED_DIR_NAMES
        for part in path.relative_to(ROOT).parts
    )


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse_sections(path: Path) -> dict[str, list[dict]]:
    sections: dict[str, list[dict]] = {}
    active_section: str | None = None

    for line_number, raw_line in enumerate(
        path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines(),
        start=1,
    ):
        content = strip_comment(raw_line)

        if not content:
            continue

        section_match = re.match(
            r"^\[\s*([^\]]+?)\s*\]$",
            content,
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

        if content.startswith("#"):
            continue

        if active_section is not None:
            sections[
                active_section
            ].append(
                {
                    "line_number": line_number,
                    "raw_line": raw_line,
                    "tokens": content.split(),
                }
            )

    return sections


def path_status(path: Path) -> str:
    lower = relative(path).lower()

    if "/accepted/" in f"/{lower}/":
        return "ACCEPTED"

    if "/candidates/" in f"/{lower}/":
        return "CANDIDATE"

    if "/runs/" in f"/{lower}/":
        return "HISTORICAL_RUN"

    if "/parameters/" in f"/{lower}/":
        return "PARAMETER_LIBRARY"

    return "OTHER"


def infer_element_from_type(atom_type: str) -> str:
    normalized = atom_type.strip().upper()

    if normalized.startswith("B"):
        return "B"

    if normalized.startswith("N"):
        return "N"

    if normalized.startswith("H"):
        return "H"

    return "UNKNOWN"


def candidate_files() -> list[Path]:
    files = []

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue

        if is_excluded(path):
            continue

        if path.suffix.lower() not in SEARCH_EXTENSIONS:
            continue

        files.append(path)

    return sorted(
        files,
        key=lambda path: relative(path),
    )


print("=" * 100)
print("DAY040 / D040-A5 — HYDROGEN PARAMETER-SOURCE AUDIT")
print("=" * 100)


print("\n[1] UPSTREAM AUTHORIZATION")

for path in (
    A4_REPORT,
    A4_PARENT_MAP,
    A4_REQUIREMENTS,
):
    require_file(path)
    print(
        f"FOUND bytes={path.stat().st_size:8d} "
        f"sha256={sha256(path)} "
        f"{path}"
    )

a4_report = load_json(
    A4_REPORT
)

if (
    a4_report.get("decision")
    != EXPECTED_A4_DECISION
):
    raise RuntimeError(
        "Unexpected D040-A4 decision.\n"
        f"Observed: {a4_report.get('decision')}"
    )

if (
    a4_report.get(
        "authorizations",
        {},
    ).get(
        "hydrogen_parameter_source_audit_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Hydrogen parameter-source audit is not authorized"
    )

print("D040_A4_decision_gate = PASS")
print("parameter_source_audit_authorization_gate = PASS")
print("parameter_adoption_blocked_gate = PASS")
print("topology_modification_blocked_gate = PASS")


print("\n[2] REQUIRED INTERACTION CONTRACT")

with A4_PARENT_MAP.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    parent_rows = list(
        csv.DictReader(handle)
    )

with A4_REQUIREMENTS.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    requirement_rows = list(
        csv.DictReader(handle)
    )

B_H_count = sum(
    row["bond_class"] == "B-H"
    for row in parent_rows
)

N_H_count = sum(
    row["bond_class"] == "N-H"
    for row in parent_rows
)

print(f"required_H_atom_count = {len(parent_rows)}")
print(f"required_B_H_instances = {B_H_count}")
print(f"required_N_H_instances = {N_H_count}")
print(
    "required_interactions = "
    "H_ATOMTYPE,B_H_BOND,N_H_BOND,H_CONTAINING_ANGLES,"
    "H_CONTAINING_IMPROPERS"
)


print("\n[3] REPOSITORY SEARCH INVENTORY")

files = candidate_files()

print(f"search_file_count = {len(files)}")

extension_counts = {}

for path in files:
    extension = path.suffix.lower()
    extension_counts[
        extension
    ] = extension_counts.get(
        extension,
        0,
    ) + 1

for extension in sorted(
    extension_counts
):
    print(
        f"{extension or '<none>'} = "
        f"{extension_counts[extension]}"
    )


print("\n[4] GROMACS ATOMTYPE CANDIDATES")

atomtype_records = []
bonded_records = []
text_hit_records = []

for path in files:
    try:
        text = path.read_text(
            encoding="utf-8",
            errors="replace",
        )
    except Exception:
        continue

    status = path_status(path)

    for pattern_name, pattern in TEXT_PATTERNS.items():
        for line_number, raw_line in enumerate(
            text.splitlines(),
            start=1,
        ):
            if pattern.search(raw_line):
                text_hit_records.append(
                    {
                        "path": relative(path),
                        "status": status,
                        "pattern": pattern_name,
                        "line_number": line_number,
                        "line": raw_line.strip(),
                    }
                )

    if path.suffix.lower() not in {
        ".itp",
        ".top",
        ".rtp",
        ".prm",
        ".str",
    }:
        continue

    sections = parse_sections(path)

    for row in sections.get(
        "atomtypes",
        [],
    ):
        tokens = row["tokens"]

        if len(tokens) < 6:
            continue

        atom_type = tokens[0]

        is_H_candidate = any(
            pattern.match(atom_type)
            for pattern in HYDROGEN_TYPE_PATTERNS
        )

        if not is_H_candidate:
            continue

        atomtype_records.append(
            {
                "path": relative(path),
                "status": status,
                "section": "atomtypes",
                "line_number": row["line_number"],
                "atom_type": atom_type,
                "inferred_element": (
                    infer_element_from_type(
                        atom_type
                    )
                ),
                "token_count": len(tokens),
                "tokens_json": json.dumps(
                    tokens
                ),
                "raw_line": row["raw_line"].strip(),
            }
        )

    for row in sections.get(
        "atoms",
        [],
    ):
        tokens = row["tokens"]

        if len(tokens) < 7:
            continue

        atom_type = tokens[1]
        atom_name = tokens[4]

        is_H_candidate = (
            any(
                pattern.match(atom_type)
                for pattern in HYDROGEN_TYPE_PATTERNS
            )
            or atom_name.upper().startswith("H")
        )

        if not is_H_candidate:
            continue

        atomtype_records.append(
            {
                "path": relative(path),
                "status": status,
                "section": "atoms",
                "line_number": row["line_number"],
                "atom_type": atom_type,
                "inferred_element": (
                    infer_element_from_type(
                        atom_type
                    )
                ),
                "token_count": len(tokens),
                "tokens_json": json.dumps(
                    tokens
                ),
                "raw_line": row["raw_line"].strip(),
            }
        )

print(
    f"hydrogen_atomtype_candidate_count = "
    f"{len(atomtype_records)}"
)

for record in atomtype_records[:100]:
    print(
        f"status={record['status']:<16s} "
        f"section={record['section']:<10s} "
        f"type={record['atom_type']:<10s} "
        f"path={record['path']} "
        f"line={record['line_number']}"
    )


print("\n[5] BONDED PARAMETER CANDIDATES")

for path in files:
    if path.suffix.lower() not in {
        ".itp",
        ".top",
        ".rtp",
        ".prm",
        ".str",
    }:
        continue

    try:
        sections = parse_sections(path)
    except Exception:
        continue

    status = path_status(path)

    for section_name in (
        "bondtypes",
        "angletypes",
        "dihedraltypes",
        "bonds",
        "angles",
        "dihedrals",
    ):
        for row in sections.get(
            section_name,
            [],
        ):
            tokens = row["tokens"]

            if section_name in {
                "bondtypes",
                "bonds",
            }:
                type_token_count = 2
            elif section_name in {
                "angletypes",
                "angles",
            }:
                type_token_count = 3
            else:
                type_token_count = 4

            if len(tokens) < type_token_count:
                continue

            leading_tokens = tokens[
                :type_token_count
            ]

            normalized = [
                token.upper()
                for token in leading_tokens
            ]

            has_H = any(
                token.startswith("H")
                for token in normalized
            )

            has_B = any(
                token.startswith("B")
                for token in normalized
            )

            has_N = any(
                token.startswith("N")
                for token in normalized
            )

            relevant = (
                has_H
                and (
                    has_B
                    or has_N
                )
            )

            if not relevant:
                continue

            if has_B and has_N:
                chemistry_class = (
                    "B_N_H_MIXED"
                )
            elif has_B:
                chemistry_class = (
                    "B_H"
                )
            elif has_N:
                chemistry_class = (
                    "N_H"
                )
            else:
                chemistry_class = (
                    "H_ONLY_OTHER"
                )

            bonded_records.append(
                {
                    "path": relative(path),
                    "status": status,
                    "section": section_name,
                    "line_number": row["line_number"],
                    "chemistry_class": chemistry_class,
                    "leading_tokens_json": (
                        json.dumps(
                            leading_tokens
                        )
                    ),
                    "all_tokens_json": (
                        json.dumps(tokens)
                    ),
                    "raw_line": (
                        row["raw_line"].strip()
                    ),
                }
            )

print(
    f"bonded_parameter_candidate_count = "
    f"{len(bonded_records)}"
)

section_counts = {}

for record in bonded_records:
    key = (
        record["section"],
        record["chemistry_class"],
    )

    section_counts[key] = (
        section_counts.get(
            key,
            0,
        )
        + 1
    )

for key in sorted(
    section_counts
):
    print(
        f"section={key[0]:<14s} "
        f"class={key[1]:<12s} "
        f"count={section_counts[key]}"
    )

for record in bonded_records[:150]:
    print(
        f"status={record['status']:<16s} "
        f"section={record['section']:<14s} "
        f"class={record['chemistry_class']:<12s} "
        f"path={record['path']} "
        f"line={record['line_number']} "
        f"record={record['raw_line']}"
    )


print("\n[6] TEXTUAL PROVENANCE HITS")

print(
    f"text_hit_count = "
    f"{len(text_hit_records)}"
)

pattern_counts = {}

for record in text_hit_records:
    pattern_counts[
        record["pattern"]
    ] = (
        pattern_counts.get(
            record["pattern"],
            0,
        )
        + 1
    )

for pattern_name in sorted(
    pattern_counts
):
    print(
        f"{pattern_name} = "
        f"{pattern_counts[pattern_name]}"
    )

for record in text_hit_records[:150]:
    print(
        f"pattern={record['pattern']:<20s} "
        f"status={record['status']:<16s} "
        f"path={record['path']} "
        f"line={record['line_number']} "
        f"text={record['line']}"
    )


print("\n[7] ACCEPTED-SOURCE FILTER")

accepted_atomtypes = [
    record
    for record in atomtype_records
    if record["status"] == "ACCEPTED"
]

accepted_bonded = [
    record
    for record in bonded_records
    if record["status"] == "ACCEPTED"
]

accepted_text_hits = [
    record
    for record in text_hit_records
    if record["status"] == "ACCEPTED"
]

accepted_B_H_bonds = [
    record
    for record in accepted_bonded
    if (
        record["section"]
        in {
            "bondtypes",
            "bonds",
        }
        and record[
            "chemistry_class"
        ]
        in {
            "B_H",
            "B_N_H_MIXED",
        }
    )
]

accepted_N_H_bonds = [
    record
    for record in accepted_bonded
    if (
        record["section"]
        in {
            "bondtypes",
            "bonds",
        }
        and record[
            "chemistry_class"
        ]
        in {
            "N_H",
            "B_N_H_MIXED",
        }
    )
]

accepted_H_angles = [
    record
    for record in accepted_bonded
    if record["section"]
    in {
        "angletypes",
        "angles",
    }
]

accepted_H_dihedrals = [
    record
    for record in accepted_bonded
    if record["section"]
    in {
        "dihedraltypes",
        "dihedrals",
    }
]

print(
    f"accepted_H_atomtype_candidate_count = "
    f"{len(accepted_atomtypes)}"
)
print(
    f"accepted_B_H_bond_candidate_count = "
    f"{len(accepted_B_H_bonds)}"
)
print(
    f"accepted_N_H_bond_candidate_count = "
    f"{len(accepted_N_H_bonds)}"
)
print(
    f"accepted_H_angle_candidate_count = "
    f"{len(accepted_H_angles)}"
)
print(
    f"accepted_H_dihedral_candidate_count = "
    f"{len(accepted_H_dihedrals)}"
)
print(
    f"accepted_text_provenance_hit_count = "
    f"{len(accepted_text_hits)}"
)


print("\n[8] SCIENTIFIC INTERPRETATION")

complete_accepted_parameter_set_found = (
    len(accepted_atomtypes) > 0
    and len(accepted_B_H_bonds) > 0
    and len(accepted_N_H_bonds) > 0
    and len(accepted_H_angles) > 0
)

partial_project_parameter_evidence_found = (
    len(atomtype_records) > 0
    or len(bonded_records) > 0
    or len(text_hit_records) > 0
)

if complete_accepted_parameter_set_found:
    interpretation = (
        "COMPLETE_ACCEPTED_HYDROGEN_PARAMETER_SET_FOUND_"
        "SCIENTIFIC_COMPARISON_REQUIRED"
    )
    next_block = (
        "D040_A6_ACCEPTED_HYDROGEN_PARAMETER_COMPARISON"
    )
elif partial_project_parameter_evidence_found:
    interpretation = (
        "PARTIAL_PROJECT_PARAMETER_EVIDENCE_FOUND_"
        "EXTERNAL_OR_PRIMARY_SOURCE_REVIEW_REQUIRED"
    )
    next_block = (
        "D040_A6_HYDROGEN_PARAMETER_PROVENANCE_REVIEW"
    )
else:
    interpretation = (
        "NO_PROJECT_PARAMETER_SOURCE_FOUND_"
        "PRIMARY_LITERATURE_PARAMETERIZATION_REQUIRED"
    )
    next_block = (
        "D040_A6_PRIMARY_SOURCE_PARAMETERIZATION_PLAN"
    )

print(
    "complete_accepted_parameter_set_found = "
    f"{complete_accepted_parameter_set_found}"
)
print(
    "partial_project_parameter_evidence_found = "
    f"{partial_project_parameter_evidence_found}"
)
print(
    f"interpretation = {interpretation}"
)
print(
    f"next_block = {next_block}"
)


print("\n[9] WRITE OUTPUTS")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

atomtype_fieldnames = [
    "path",
    "status",
    "section",
    "line_number",
    "atom_type",
    "inferred_element",
    "token_count",
    "tokens_json",
    "raw_line",
]

with ATOMTYPE_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=atomtype_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        atomtype_records
    )

bonded_fieldnames = [
    "path",
    "status",
    "section",
    "line_number",
    "chemistry_class",
    "leading_tokens_json",
    "all_tokens_json",
    "raw_line",
]

with BONDED_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=bonded_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        bonded_records
    )

text_hit_fieldnames = [
    "path",
    "status",
    "pattern",
    "line_number",
    "line",
]

with TEXT_HITS_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=text_hit_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        text_hit_records
    )


print("\n[10] GATES")

gates = {
    "D040_A4_decision_gate": True,
    "required_parent_map_loaded_gate": (
        len(parent_rows) == 6
    ),
    "required_B_H_count_gate": (
        B_H_count == 5
    ),
    "required_N_H_count_gate": (
        N_H_count == 1
    ),
    "repository_search_inventory_gate": (
        len(files) > 0
    ),
    "atomtype_output_created_gate": (
        ATOMTYPE_CSV.is_file()
        and ATOMTYPE_CSV.stat().st_size > 0
    ),
    "bonded_output_created_gate": (
        BONDED_CSV.is_file()
        and BONDED_CSV.stat().st_size > 0
    ),
    "text_hits_output_created_gate": (
        TEXT_HITS_CSV.is_file()
        and TEXT_HITS_CSV.stat().st_size > 0
    ),
    "no_parameter_adopted_gate": True,
    "no_topology_modified_gate": True,
    "no_coordinates_modified_gate": True,
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
    "D040_A5_HYDROGEN_PARAMETER_SOURCE_AUDIT_PASS_"
    + next_block
    + "_AUTHORIZED"
    if all_gates_pass
    else
    "D040_A5_HYDROGEN_PARAMETER_SOURCE_AUDIT_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_identity": {
        "A4_report": {
            "path": relative(
                A4_REPORT
            ),
            "sha256": sha256(
                A4_REPORT
            ),
        },
        "A4_parent_map": {
            "path": relative(
                A4_PARENT_MAP
            ),
            "sha256": sha256(
                A4_PARENT_MAP
            ),
        },
        "A4_requirements": {
            "path": relative(
                A4_REQUIREMENTS
            ),
            "sha256": sha256(
                A4_REQUIREMENTS
            ),
        },
    },
    "required_contract": {
        "hydrogen_count": len(
            parent_rows
        ),
        "B_H_instances": B_H_count,
        "N_H_instances": N_H_count,
    },
    "search_inventory": {
        "file_count": len(files),
        "extension_counts": (
            extension_counts
        ),
    },
    "candidate_counts": {
        "all_H_atomtype_candidates": (
            len(atomtype_records)
        ),
        "all_bonded_candidates": (
            len(bonded_records)
        ),
        "all_text_hits": (
            len(text_hit_records)
        ),
        "accepted_H_atomtype_candidates": (
            len(accepted_atomtypes)
        ),
        "accepted_B_H_bond_candidates": (
            len(accepted_B_H_bonds)
        ),
        "accepted_N_H_bond_candidates": (
            len(accepted_N_H_bonds)
        ),
        "accepted_H_angle_candidates": (
            len(accepted_H_angles)
        ),
        "accepted_H_dihedral_candidates": (
            len(accepted_H_dihedrals)
        ),
    },
    "scientific_interpretation": {
        "complete_accepted_parameter_set_found": (
            complete_accepted_parameter_set_found
        ),
        "partial_project_parameter_evidence_found": (
            partial_project_parameter_evidence_found
        ),
        "interpretation": interpretation,
        "parameter_adoption_status": (
            "NOT_AUTHORIZED"
        ),
    },
    "accepted_candidates": {
        "atomtypes": accepted_atomtypes,
        "B_H_bonds": accepted_B_H_bonds,
        "N_H_bonds": accepted_N_H_bonds,
        "angles": accepted_H_angles,
        "dihedrals": (
            accepted_H_dihedrals
        ),
        "text_hits": accepted_text_hits,
    },
    "gates": gates,
    "authorizations": {
        "hydrogen_parameter_provenance_review_authorized": (
            all_gates_pass
        ),
        "primary_literature_parameter_search_authorized": (
            all_gates_pass
            and not complete_accepted_parameter_set_found
        ),
        "parameter_comparison_authorized": (
            all_gates_pass
            and complete_accepted_parameter_set_found
        ),
        "new_atom_type_definition_authorized": False,
        "bonded_parameter_modification_authorized": False,
        "hydrogen_coordinate_insertion_authorized": False,
        "charge_to_topology_mapping_execution_authorized": False,
        "topology_modification_authorized": False,
        "force_field_adoption_authorized": False,
        "energy_execution_authorized": False,
        "minimization_execution_authorized": False,
        "validation_MD_execution_authorized": False,
        "production_MD_authorized": False,
    },
    "next_required_block": {
        "name": next_block,
        "required_actions": (
            [
                (
                    "Compare all accepted H atom types and bonded "
                    "terms for units, functional form and provenance."
                ),
                (
                    "Determine whether they are chemically applicable "
                    "to B-H and N-H edge passivation in this model."
                ),
            ]
            if complete_accepted_parameter_set_found
            else [
                (
                    "Identify primary literature or authoritative "
                    "force-field sources for edge-passivated hBN."
                ),
                (
                    "Establish H nonbonded, B-H, N-H, angle and "
                    "improper terms with explicit provenance."
                ),
                (
                    "Do not infer missing constants from geometry alone."
                ),
            ]
        ),
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

print(f"atomtype_csv = {ATOMTYPE_CSV}")
print(
    f"atomtype_csv_sha256 = "
    f"{sha256(ATOMTYPE_CSV)}"
)
print(f"bonded_csv = {BONDED_CSV}")
print(
    f"bonded_csv_sha256 = "
    f"{sha256(BONDED_CSV)}"
)
print(f"text_hits_csv = {TEXT_HITS_CSV}")
print(
    f"text_hits_csv_sha256 = "
    f"{sha256(TEXT_HITS_CSV)}"
)
print(f"report_json = {REPORT_JSON}")
print(
    f"report_json_sha256 = "
    f"{sha256(REPORT_JSON)}"
)


print("\n[12] DECISION")

print(f"decision={decision}")
print(
    "hydrogen_parameter_provenance_review_authorized="
    f"{all_gates_pass}"
)
print(
    "primary_literature_parameter_search_authorized="
    f"{all_gates_pass and not complete_accepted_parameter_set_found}"
)
print(
    "parameter_comparison_authorized="
    f"{all_gates_pass and complete_accepted_parameter_set_found}"
)
print(
    "new_atom_type_definition_authorized=False"
)
print(
    "bonded_parameter_modification_authorized=False"
)
print(
    "hydrogen_coordinate_insertion_authorized=False"
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
