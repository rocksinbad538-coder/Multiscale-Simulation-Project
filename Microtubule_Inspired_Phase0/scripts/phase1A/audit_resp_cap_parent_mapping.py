#!/usr/bin/env python3
"""
DAY038 / D038-F2

Authoritative artificial-cap to real-atom mapping inventory for
QM_F06_UPPER_V7A_R1.

The script searches project artifacts for explicit evidence connecting
each artificial QM boundary cap to one or more real atoms.

No charge redistribution is performed.
RESP Stage 2 remains blocked.
"""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_PARENT = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

EXPECTED_EQUIVALENCE_DECISION = (
    "D038_F1_RESP_STAGE2_EQUIVALENCE_DESIGN_PASS_"
    "NO_STAGE2_EQUIVALENCE_AUTOMATICALLY_AUTHORIZED"
)

SEARCH_ROOTS = (
    ROOT / "runs/phase1A/day035_qm_f06_upper_v7a_r1_resp_input_design",
    ROOT / "runs/phase1A/day035_qm_f06_upper_v7a_r1_coordinate_adoption",
    ROOT / "runs/phase1A/day036_qm_f06_upper_v7a_r1_resp_preparation",
    ROOT / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions",
    ROOT / "scripts/phase1A",
)

TEXT_SUFFIXES = {
    ".json",
    ".csv",
    ".tsv",
    ".txt",
    ".md",
    ".map",
    ".dat",
    ".in",
    ".inp",
    ".xyz",
    ".pdb",
    ".mol2",
    ".ac",
    ".py",
}


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()

    if normalized in {"true", "1", "yes", "y"}:
        return True

    if normalized in {"false", "0", "no", "n", ""}:
        return False

    raise RuntimeError(
        f"Unrecognized Boolean value: {value!r}"
    )


def extract_encoded_reference_tokens(
    atom_id: str,
) -> list[str]:
    """
    Extract identifier-like tokens from cap names.

    These are candidate references only. They are not treated as
    verified cap-parent assignments.
    """

    tokens = re.findall(
        r"\b(?:P|S|A|BR)\d+(?:_\d+)*\b",
        atom_id,
    )

    unique_tokens = []

    for token in tokens:
        if token not in unique_tokens:
            unique_tokens.append(token)

    return unique_tokens


print("=" * 100)
print("DAY038 / D038-F2 — ARTIFICIAL-CAP PARENT-MAPPING INVENTORY")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_POINTER)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

equivalence_design_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE2_EQUIVALENCE_DESIGN.json"
)

require_file(transferability_csv)
require_file(equivalence_design_json)

print(f"execution_dir          = {execution_dir}")
print(f"transferability_csv    = {transferability_csv}")
print(f"equivalence_design     = {equivalence_design_json}")


print("\n[2] UPSTREAM AUTHORIZATION")

equivalence_design = load_json(
    equivalence_design_json
)

if (
    equivalence_design.get("decision")
    != EXPECTED_EQUIVALENCE_DECISION
):
    raise RuntimeError(
        "Unexpected Stage 2 equivalence-design decision.\n"
        f"Observed: {equivalence_design.get('decision')}"
    )

authorizations = equivalence_design.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "RESP_stage2_protocol_design_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Stage 2 protocol design is not authorized"
    )

if (
    authorizations.get(
        "RESP_stage2_execution_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected Stage 2 execution authorization"
    )

print("upstream_decision_gate          = PASS")
print("protocol_design_authorized_gate = PASS")
print("stage2_execution_blocked_gate   = PASS")


print("\n[3] LOAD ARTIFICIAL CAPS")

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    atom_rows = list(
        csv.DictReader(handle)
    )

if len(atom_rows) != 52:
    raise RuntimeError(
        f"Expected 52 atoms, observed {len(atom_rows)}"
    )

caps = []

real_atoms = []

for row in atom_rows:
    row["atom_index_0based_int"] = int(
        row["atom_index_0based"]
    )
    row["RESP_stage1_charge_e_float"] = float(
        row["RESP_stage1_charge_e"]
    )
    row["artificial_cap_bool"] = parse_bool(
        row["artificial_cap"]
    )

    if row["artificial_cap_bool"]:
        caps.append(row)
    else:
        real_atoms.append(row)

if len(caps) != 15:
    raise RuntimeError(
        f"Expected 15 artificial caps, observed {len(caps)}"
    )

cap_charge_sum = sum(
    row["RESP_stage1_charge_e_float"]
    for row in caps
)

real_charge_sum = sum(
    row["RESP_stage1_charge_e_float"]
    for row in real_atoms
)

print(f"artificial_cap_count = {len(caps)}")
print(f"real_atom_count       = {len(real_atoms)}")
print(f"cap_charge_sum_e      = {cap_charge_sum:.16g}")
print(f"real_charge_sum_e     = {real_charge_sum:.16g}")
print("cap_inventory_gate    = PASS")


print("\n[4] BUILD SEARCH CORPUS")

candidate_files: list[Path] = []

for search_root in SEARCH_ROOTS:
    if not search_root.exists():
        print(f"MISSING_SEARCH_ROOT  {search_root}")
        continue

    for path in search_root.rglob("*"):
        if (
            path.is_file()
            and path.suffix.lower() in TEXT_SUFFIXES
        ):
            candidate_files.append(path)

candidate_files = sorted(
    set(candidate_files)
)

print(f"search_file_count = {len(candidate_files)}")

text_cache: dict[Path, list[str]] = {}

for path in candidate_files:
    try:
        text_cache[path] = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()
    except Exception:
        continue

print(f"readable_text_file_count = {len(text_cache)}")


print("\n[5] CAP-BY-CAP EVIDENCE SEARCH")

mapping_records = []

for cap in caps:
    atom_id = cap["atom_id"].strip()

    encoded_tokens = extract_encoded_reference_tokens(
        atom_id
    )

    exact_occurrences = []
    token_occurrences = []

    for path, lines in text_cache.items():
        relative_path = str(
            path.relative_to(ROOT)
        )

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            if atom_id in line:
                exact_occurrences.append(
                    {
                        "path": relative_path,
                        "line_number": line_number,
                        "line": line.strip(),
                    }
                )

            matched_tokens = [
                token
                for token in encoded_tokens
                if token in line
            ]

            if matched_tokens:
                token_occurrences.append(
                    {
                        "path": relative_path,
                        "line_number": line_number,
                        "matched_tokens": matched_tokens,
                        "line": line.strip(),
                    }
                )

    mapping_status = (
        "EXPLICIT_CAP_ID_EVIDENCE_FOUND"
        if exact_occurrences
        else
        "NO_EXPLICIT_CAP_ID_EVIDENCE"
    )

    record = {
        "cap_atom_index_0based": (
            cap["atom_index_0based_int"]
        ),
        "cap_atom_index_1based": int(
            cap["atom_index_1based"]
        ),
        "cap_atom_id": atom_id,
        "cap_role": cap["atom_role"],
        "cap_RESP_stage1_charge_e": (
            cap["RESP_stage1_charge_e_float"]
        ),
        "cap_coordinates_A": {
            "x": float(cap["x_A"]),
            "y": float(cap["y_A"]),
            "z": float(cap["z_A"]),
        },
        "encoded_reference_tokens_unverified": (
            encoded_tokens
        ),
        "exact_cap_id_occurrence_count": len(
            exact_occurrences
        ),
        "encoded_token_occurrence_count": len(
            token_occurrences
        ),
        "exact_cap_id_occurrences": (
            exact_occurrences
        ),
        "encoded_token_occurrences": (
            token_occurrences
        ),
        "mapping_status": mapping_status,
        "verified_parent_atom_indices_0based": [],
        "verified_parent_atom_ids": [],
        "charge_redistribution_authorized": False,
    }

    mapping_records.append(record)

    print(
        f"\ncap_atom={record['cap_atom_index_0based']:>2} "
        f"id={atom_id} "
        f"charge={record['cap_RESP_stage1_charge_e']: .6f}"
    )
    print(
        f"  encoded_tokens={encoded_tokens}"
    )
    print(
        f"  exact_cap_id_occurrences="
        f"{len(exact_occurrences)}"
    )
    print(
        f"  token_occurrences="
        f"{len(token_occurrences)}"
    )
    print(
        f"  mapping_status={mapping_status}"
    )

    for occurrence in exact_occurrences[:10]:
        print(
            f"    EXACT {occurrence['path']}:"
            f"{occurrence['line_number']} "
            f"{occurrence['line']}"
        )

    if len(exact_occurrences) > 10:
        print(
            f"    ... {len(exact_occurrences) - 10} "
            "additional exact occurrences"
        )


print("\n[6] EVIDENCE SUMMARY")

explicit_evidence_count = sum(
    record["mapping_status"]
    == "EXPLICIT_CAP_ID_EVIDENCE_FOUND"
    for record in mapping_records
)

verified_parent_count = sum(
    bool(
        record[
            "verified_parent_atom_indices_0based"
        ]
    )
    for record in mapping_records
)

print(
    f"caps_with_explicit_id_evidence = "
    f"{explicit_evidence_count}/{len(mapping_records)}"
)
print(
    f"caps_with_verified_parent_mapping = "
    f"{verified_parent_count}/{len(mapping_records)}"
)
print(
    f"unresolved_parent_mapping_count = "
    f"{len(mapping_records) - verified_parent_count}"
)


print("\n[7] WRITE OUTPUTS")

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_MAPPING_INVENTORY.json"
)

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_MAPPING_INVENTORY.csv"
)

csv_fieldnames = [
    "cap_atom_index_0based",
    "cap_atom_index_1based",
    "cap_atom_id",
    "cap_role",
    "cap_RESP_stage1_charge_e",
    "encoded_reference_tokens_unverified",
    "exact_cap_id_occurrence_count",
    "encoded_token_occurrence_count",
    "mapping_status",
    "verified_parent_atom_indices_0based",
    "verified_parent_atom_ids",
    "charge_redistribution_authorized",
]

with output_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=csv_fieldnames,
    )

    writer.writeheader()

    for record in mapping_records:
        writer.writerow(
            {
                "cap_atom_index_0based": (
                    record["cap_atom_index_0based"]
                ),
                "cap_atom_index_1based": (
                    record["cap_atom_index_1based"]
                ),
                "cap_atom_id": (
                    record["cap_atom_id"]
                ),
                "cap_role": (
                    record["cap_role"]
                ),
                "cap_RESP_stage1_charge_e": (
                    record[
                        "cap_RESP_stage1_charge_e"
                    ]
                ),
                "encoded_reference_tokens_unverified": (
                    json.dumps(
                        record[
                            "encoded_reference_tokens_unverified"
                        ]
                    )
                ),
                "exact_cap_id_occurrence_count": (
                    record[
                        "exact_cap_id_occurrence_count"
                    ]
                ),
                "encoded_token_occurrence_count": (
                    record[
                        "encoded_token_occurrence_count"
                    ]
                ),
                "mapping_status": (
                    record["mapping_status"]
                ),
                "verified_parent_atom_indices_0based": (
                    json.dumps(
                        record[
                            "verified_parent_atom_indices_0based"
                        ]
                    )
                ),
                "verified_parent_atom_ids": (
                    json.dumps(
                        record[
                            "verified_parent_atom_ids"
                        ]
                    )
                ),
                "charge_redistribution_authorized": False,
            }
        )


gates = {
    "upstream_decision_gate": True,
    "protocol_design_authorized_gate": True,
    "atom_count_gate": len(atom_rows) == 52,
    "artificial_cap_count_gate": len(caps) == 15,
    "real_atom_count_gate": len(real_atoms) == 37,
    "finite_charge_gate": all(
        abs(
            row["RESP_stage1_charge_e_float"]
        ) < float("inf")
        for row in atom_rows
    ),
    "search_corpus_gate": len(text_cache) > 0,
    "no_charge_redistribution_performed_gate": all(
        record["charge_redistribution_authorized"]
        is False
        for record in mapping_records
    ),
}

all_gates_pass = all(
    gates.values()
)

decision = (
    "D038_F2_CAP_PARENT_MAPPING_INVENTORY_PASS_"
    "PARENT_ASSIGNMENTS_REQUIRE_EVIDENCE_REVIEW"
    if all_gates_pass
    else
    "D038_F2_CAP_PARENT_MAPPING_INVENTORY_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_execution_directory": str(
        execution_dir.resolve()
    ),
    "source_identity": {
        "transferability_csv": str(
            transferability_csv.resolve()
        ),
        "transferability_csv_sha256": sha256(
            transferability_csv
        ),
        "equivalence_design_json": str(
            equivalence_design_json.resolve()
        ),
        "equivalence_design_json_sha256": sha256(
            equivalence_design_json
        ),
    },
    "charge_summary": {
        "artificial_cap_count": len(caps),
        "real_atom_count": len(real_atoms),
        "artificial_cap_RESP_charge_sum_e": (
            cap_charge_sum
        ),
        "real_atom_RESP_charge_sum_e": (
            real_charge_sum
        ),
    },
    "evidence_summary": {
        "search_file_count": len(
            candidate_files
        ),
        "readable_text_file_count": len(
            text_cache
        ),
        "caps_with_explicit_id_evidence": (
            explicit_evidence_count
        ),
        "caps_with_verified_parent_mapping": (
            verified_parent_count
        ),
        "unresolved_parent_mapping_count": (
            len(mapping_records)
            - verified_parent_count
        ),
    },
    "cap_mapping_records": mapping_records,
    "gates": gates,
    "authorizations": {
        "cap_parent_mapping_review_authorized": (
            all_gates_pass
        ),
        "cap_charge_redistribution_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "mapping_csv": str(
            output_csv.resolve()
        ),
    },
}

output_json.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"output_csv = {output_csv}")
print(
    f"output_csv_sha256 = "
    f"{sha256(output_csv)}"
)
print(f"output_json = {output_json}")
print(
    f"output_json_sha256 = "
    f"{sha256(output_json)}"
)


print("\n[8] GATES")

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )


print("\n[9] DECISION")

print(f"decision = {decision}")
print(
    "cap_parent_mapping_review_authorized = "
    f"{all_gates_pass}"
)
print("cap_charge_redistribution_authorized = False")
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
