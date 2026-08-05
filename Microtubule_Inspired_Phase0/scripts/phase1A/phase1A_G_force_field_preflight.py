#!/usr/bin/env python3
"""
DAY040 / D040-A1

Phase 1A-G force-field integration preflight.

This block inventories candidate topology and force-field files without
modifying them. It binds the preflight to the formally closed Phase 1A-F
working charge model.

No topology is modified.
No force field is adopted.
No MD is executed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

CLOSURE_POINTER = (
    ROOT
    / "runs/phase1A"
    / "LATEST_PHASE1A_F_CHARGE_MODEL_CLOSURE.txt"
)

EXPECTED_CLOSURE_DECISION = (
    "D039_A16_PHASE1A_F_CLOSED_"
    "WORKING_CHARGE_MODEL_ADOPTED_FOR_"
    "FORCE_FIELD_INTEGRATION"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_force_field_preflight"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_G_FORCE_FIELD_PREFLIGHT.json"
)

CANDIDATE_EXTENSIONS = {
    ".top",
    ".itp",
    ".gro",
    ".pdb",
    ".mol2",
    ".frcmod",
    ".prmtop",
    ".parm7",
    ".inpcrd",
    ".rst7",
    ".psf",
    ".rtf",
    ".prm",
    ".str",
    ".xml",
    ".off",
    ".lib",
    ".prep",
    ".in",
    ".mdp",
    ".lmp",
    ".data",
}

EXCLUDED_DIRECTORY_NAMES = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "node_modules",
}


def file_record(path: Path) -> dict:
    return {
        "path": str(path.relative_to(ROOT)),
        "absolute_path": str(path.resolve()),
        "suffix": path.suffix.lower(),
        "bytes": int(path.stat().st_size),
        "sha256": sha256(path),
    }


def detect_format(path: Path) -> str:
    suffix = path.suffix.lower()

    if suffix in {
        ".top",
        ".itp",
        ".gro",
        ".mdp",
    }:
        return "GROMACS"

    if suffix in {
        ".prmtop",
        ".parm7",
        ".inpcrd",
        ".rst7",
        ".mol2",
        ".frcmod",
        ".off",
        ".lib",
        ".prep",
    }:
        return "AMBER"

    if suffix in {
        ".psf",
        ".rtf",
        ".prm",
        ".str",
    }:
        return "CHARMM"

    if suffix == ".xml":
        return "OPENMM_OR_GENERIC_XML"

    if suffix in {
        ".lmp",
        ".data",
    }:
        return "LAMMPS_OR_GENERIC_DATA"

    if suffix == ".pdb":
        return "STRUCTURE_PDB"

    if suffix == ".in":
        return "GENERIC_INPUT"

    return "OTHER"


def safe_text_probe(path: Path) -> dict:
    try:
        text = path.read_text(
            encoding="utf-8",
            errors="replace",
        )
    except Exception as error:
        return {
            "readable_as_text": False,
            "error": str(error),
        }

    lowered = text.lower()

    probes = {
        "contains_atomtypes": (
            "[ atomtypes ]" in lowered
            or "mass" in lowered
        ),
        "contains_atoms_section": (
            "[ atoms ]" in lowered
            or "%flag charge" in lowered
        ),
        "contains_bonds": (
            "[ bonds ]" in lowered
            or "%flag bonds" in lowered
        ),
        "contains_angles": (
            "[ angles ]" in lowered
            or "%flag angles" in lowered
        ),
        "contains_dihedrals": (
            "[ dihedrals ]" in lowered
            or "%flag dihedrals" in lowered
        ),
        "mentions_boron": (
            " boron" in lowered
            or "\nb " in lowered
            or " b " in lowered
        ),
        "mentions_nitrogen": (
            " nitrogen" in lowered
            or "\nn " in lowered
            or " n " in lowered
        ),
        "mentions_QM_F06": (
            "qm_f06" in lowered
        ),
        "mentions_upper_v7a": (
            "upper_v7a" in lowered
        ),
    }

    return {
        "readable_as_text": True,
        "line_count": len(
            text.splitlines()
        ),
        "character_count": len(text),
        "probes": probes,
    }


print("=" * 100)
print("DAY040 / D040-A1 — FORCE-FIELD INTEGRATION PREFLIGHT")
print("=" * 100)


print("\n[1] PHASE 1A-F CLOSURE BINDING")

require_file(CLOSURE_POINTER)

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

charges_json = (
    closure_dir
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.json"
)

manifest_json = (
    closure_dir
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_MANIFEST.json"
)

for path in (
    closure_json,
    charges_csv,
    charges_json,
    manifest_json,
):
    require_file(path)
    print(f"FOUND  {path}")

closure_report = load_json(
    closure_json
)

if (
    closure_report.get("decision")
    != EXPECTED_CLOSURE_DECISION
):
    raise RuntimeError(
        "Unexpected Phase 1A-F closure decision.\n"
        f"Observed: {closure_report.get('decision')}"
    )

authorizations = closure_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "Phase1A_G_force_field_integration_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Phase 1A-G force-field integration is not authorized"
    )

if (
    authorizations.get(
        "force_field_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected force-field adoption authorization"
    )

print("Phase1A_F_closure_gate = PASS")
print("Phase1A_G_integration_authorization_gate = PASS")
print("force_field_adoption_blocked_gate = PASS")


print("\n[2] INVENTORY CANDIDATE FILES")

candidate_files = []

for path in ROOT.rglob("*"):
    if not path.is_file():
        continue

    relative_parts = path.relative_to(
        ROOT
    ).parts

    if any(
        part in EXCLUDED_DIRECTORY_NAMES
        for part in relative_parts
    ):
        continue

    if path.suffix.lower() not in CANDIDATE_EXTENSIONS:
        continue

    candidate_files.append(path)

candidate_files.sort(
    key=lambda path: str(
        path.relative_to(ROOT)
    )
)

print(
    f"candidate_file_count = "
    f"{len(candidate_files)}"
)

format_counts = {}

records = []

for path in candidate_files:
    detected_format = detect_format(
        path
    )

    format_counts[
        detected_format
    ] = (
        format_counts.get(
            detected_format,
            0,
        )
        + 1
    )

    record = file_record(path)
    record["detected_format"] = (
        detected_format
    )
    record["text_probe"] = (
        safe_text_probe(path)
    )

    records.append(record)

for format_name in sorted(
    format_counts
):
    print(
        f"{format_name} = "
        f"{format_counts[format_name]}"
    )


print("\n[3] HIGH-RELEVANCE CANDIDATES")

high_relevance = []

for record in records:
    path_lower = record[
        "path"
    ].lower()

    probe = record[
        "text_probe"
    ]

    probe_flags = (
        probe.get(
            "probes",
            {},
        )
        if probe.get(
            "readable_as_text"
        )
        else {}
    )

    score = 0
    reasons = []

    if "phase1a" in path_lower:
        score += 2
        reasons.append(
            "PATH_CONTAINS_PHASE1A"
        )

    if "force" in path_lower:
        score += 2
        reasons.append(
            "PATH_CONTAINS_FORCE"
        )

    if "topolog" in path_lower:
        score += 3
        reasons.append(
            "PATH_CONTAINS_TOPOLOGY"
        )

    if "qm_f06" in path_lower:
        score += 5
        reasons.append(
            "PATH_CONTAINS_QM_F06"
        )

    if "upper_v7a" in path_lower:
        score += 5
        reasons.append(
            "PATH_CONTAINS_UPPER_V7A"
        )

    if probe_flags.get(
        "mentions_QM_F06"
    ):
        score += 5
        reasons.append(
            "CONTENT_MENTIONS_QM_F06"
        )

    if probe_flags.get(
        "mentions_upper_v7a"
    ):
        score += 5
        reasons.append(
            "CONTENT_MENTIONS_UPPER_V7A"
        )

    if probe_flags.get(
        "contains_atoms_section"
    ):
        score += 2
        reasons.append(
            "CONTAINS_ATOMS_SECTION"
        )

    if probe_flags.get(
        "contains_atomtypes"
    ):
        score += 2
        reasons.append(
            "CONTAINS_ATOMTYPES"
        )

    if score > 0:
        high_relevance.append(
            {
                "score": score,
                "reasons": reasons,
                **record,
            }
        )

high_relevance.sort(
    key=lambda item: (
        -item["score"],
        item["path"],
    )
)

print(
    f"high_relevance_candidate_count = "
    f"{len(high_relevance)}"
)

for record in high_relevance[:50]:
    print(
        f"score={record['score']:2d} "
        f"format={record['detected_format']:<24s} "
        f"bytes={record['bytes']:8d} "
        f"path={record['path']}"
    )


print("\n[4] ENGINE-SPECIFIC INVENTORY")

engine_groups = {}

for record in records:
    engine_groups.setdefault(
        record["detected_format"],
        [],
    ).append(
        record["path"]
    )

for engine_name in sorted(
    engine_groups
):
    print(
        f"\n{engine_name}"
    )

    for path in engine_groups[
        engine_name
    ][:40]:
        print(f"  {path}")


print("\n[5] PREFLIGHT GATES")

closure_binding_gate = True

adopted_charge_artifact_gate = all(
    path.is_file()
    and path.stat().st_size > 0
    for path in (
        charges_csv,
        charges_json,
        manifest_json,
    )
)

candidate_inventory_gate = (
    len(candidate_files) > 0
)

recognized_force_field_format_gate = any(
    format_name in format_counts
    for format_name in (
        "GROMACS",
        "AMBER",
        "CHARMM",
        "OPENMM_OR_GENERIC_XML",
        "LAMMPS_OR_GENERIC_DATA",
    )
)

no_files_modified_gate = True
force_field_adoption_blocked_gate = True
MD_execution_blocked_gate = True

gates = {
    "closure_binding_gate": (
        closure_binding_gate
    ),
    "adopted_charge_artifact_gate": (
        adopted_charge_artifact_gate
    ),
    "candidate_inventory_gate": (
        candidate_inventory_gate
    ),
    "recognized_force_field_format_gate": (
        recognized_force_field_format_gate
    ),
    "no_files_modified_gate": (
        no_files_modified_gate
    ),
    "force_field_adoption_blocked_gate": (
        force_field_adoption_blocked_gate
    ),
    "MD_execution_blocked_gate": (
        MD_execution_blocked_gate
    ),
}

for gate_name, value in gates.items():
    print(
        f"{gate_name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[6] WRITE REPORT")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

decision = (
    "D040_A1_FORCE_FIELD_INTEGRATION_PREFLIGHT_PASS_"
    "TARGET_SELECTION_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D040_A1_FORCE_FIELD_INTEGRATION_PREFLIGHT_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "Phase1A_F_closure": {
        "closure_json": str(
            closure_json.resolve()
        ),
        "closure_json_sha256": sha256(
            closure_json
        ),
        "working_charges_csv": str(
            charges_csv.resolve()
        ),
        "working_charges_csv_sha256": sha256(
            charges_csv
        ),
        "working_charges_json": str(
            charges_json.resolve()
        ),
        "working_charges_json_sha256": sha256(
            charges_json
        ),
    },
    "candidate_file_count": len(
        candidate_files
    ),
    "format_counts": format_counts,
    "candidate_files": records,
    "high_relevance_candidates": (
        high_relevance
    ),
    "gates": gates,
    "authorizations": {
        "force_field_target_selection_review_authorized": (
            all_gates_pass
        ),
        "charge_to_topology_mapping_authorized": False,
        "topology_modification_authorized": False,
        "force_field_adoption_authorized": False,
        "energy_execution_authorized": False,
        "minimization_execution_authorized": False,
        "validation_MD_execution_authorized": False,
        "production_MD_authorized": False,
    },
    "next_required_block": {
        "name": (
            "D040_A2_FORCE_FIELD_TARGET_SELECTION"
        ),
        "required_review": [
            (
                "Select the authoritative topology "
                "or structure target."
            ),
            (
                "Identify the simulation engine and "
                "force-field family."
            ),
            (
                "Verify whether B and N atom types and "
                "bonded parameters already exist."
            ),
            (
                "Establish the exact 37-atom mapping "
                "before modifying any topology."
            ),
        ],
    },
}

REPORT_PATH.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"report_path = {REPORT_PATH}")
print(
    f"report_sha256 = "
    f"{sha256(REPORT_PATH)}"
)


print("\n[7] DECISION")

print(f"decision={decision}")
print(
    "force_field_target_selection_review_authorized="
    f"{all_gates_pass}"
)
print(
    "charge_to_topology_mapping_authorized=False"
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
