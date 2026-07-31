#!/usr/bin/env python3
"""
Independent final audit of the QM_F06 UPPER V7-A R1 ORCA
ESP/CHELPG single-point calculation.

The audit validates:
- upstream ESP-input preflight;
- execution/input identity;
- ORCA shell status;
- normal termination and SCF convergence;
- absence of fatal markers and stderr output;
- final single-point energy;
- CHELPG charge block completeness and total charge;
- stable electronic-structure outputs.

This audit does not authorize RESP execution. ORCA CHELPG charges are
retained as diagnostic fitted charges; conversion or extraction of the
raw QM ESP dataset for AmberTools RESP requires a separate validated
stage.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import math
import re


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_PARENT = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_r1_esp_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_ESP_UPPER_V7A_R1_EXECUTION.txt"
)

DAY036_PREPARATION_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_r1_resp_preparation"
)

PREFLIGHT_REPORT = (
    DAY036_PREPARATION_DIR
    / "QM_F06_UPPER_V7A_R1_ESP_INPUT_PREFLIGHT.json"
)

OUTPUT_AUDIT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_r1_esp_audit"
)

OUTPUT_REPORT = (
    OUTPUT_AUDIT_DIR
    / "QM_F06_UPPER_V7A_R1_ESP_RESULT_AUDIT.json"
)

OUTPUT_CHARGES = (
    OUTPUT_AUDIT_DIR
    / "QM_F06_UPPER_V7A_R1_CHELPG_CHARGES.csv"
)

OUTPUT_SUMMARY = (
    OUTPUT_AUDIT_DIR
    / "QM_F06_UPPER_V7A_R1_ESP_RESULT_SUMMARY.md"
)

EXPECTED_PREFLIGHT_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "ESP_INPUT_PREFLIGHT_PASS_"
    "ESP_EXECUTION_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

TOTAL_CHARGE_TOLERANCE_E = 1.0e-5

ENERGY_RE = re.compile(
    r"FINAL SINGLE POINT ENERGY\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
)

SCF_RE = re.compile(
    r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES",
    re.IGNORECASE,
)

CHELPG_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]+)\s*:\s*"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s*$"
)

TOTAL_CHARGE_RE = re.compile(
    r"Total charge\s*:\s*"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)",
    re.IGNORECASE,
)

FATAL_PATTERNS = {
    "orca_error": re.compile(
        r"\bORCA\b.*\bERROR\b|\bERROR\b.*\bORCA\b",
        re.IGNORECASE,
    ),
    "abort": re.compile(
        r"\bABORT(?:ED)?\b",
        re.IGNORECASE,
    ),
    "fatal": re.compile(
        r"\bFATAL\b",
        re.IGNORECASE,
    ),
    "segmentation_fault": re.compile(
        r"segmentation fault|segfault",
        re.IGNORECASE,
    ),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(
            f"Required JSON is missing: {path}"
        )

    return json.loads(
        path.read_text(encoding="utf-8")
    )


def resolve_execution_directory() -> Path:
    if not LATEST_POINTER.is_file():
        raise RuntimeError(
            f"Missing execution pointer: {LATEST_POINTER}"
        )

    raw = LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()

    if not raw:
        raise RuntimeError(
            f"Empty execution pointer: {LATEST_POINTER}"
        )

    path = Path(raw)

    if not path.is_absolute():
        path = ROOT / path

    path = path.resolve()

    if not path.is_dir():
        raise RuntimeError(
            f"Execution directory does not exist: {path}"
        )

    return path


def parse_chelpg_charges(
    text: str,
) -> tuple[list[dict], float | None]:
    lines = text.splitlines()

    headers = [
        index
        for index, line in enumerate(lines)
        if line.strip() == "CHELPG Charges"
    ]

    if not headers:
        return [], None

    start = headers[-1]
    rows = []
    reported_total = None
    parsing_started = False

    for line in lines[start + 1:start + 150]:
        total_match = TOTAL_CHARGE_RE.search(line)

        if total_match:
            reported_total = float(
                total_match.group(1)
            )

        match = CHELPG_ROW_RE.match(line)

        if match:
            parsing_started = True
            index, element, charge = match.groups()

            rows.append({
                "atom_index_0based": int(index),
                "atom_index_1based": int(index) + 1,
                "element": element,
                "CHELPG_charge_e": float(charge),
            })
            continue

        if (
            parsing_started
            and rows
            and not line.strip()
        ):
            # Do not stop immediately; ORCA may print a blank
            # before the total-charge line.
            continue

        if (
            parsing_started
            and rows
            and (
                "Total charge" in line
                or "CHELPG" in line
            )
        ):
            continue

    return rows, reported_total


def write_csv(
    path: Path,
    rows: list[dict],
) -> None:
    if not rows:
        raise RuntimeError(
            "No CHELPG rows available for CSV output."
        )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT_AUDIT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    preflight = read_json(
        PREFLIGHT_REPORT
    )

    execution_dir = (
        resolve_execution_directory()
    )

    input_path = (
        execution_dir
        / "esp_upper_v7a_r1.inp"
    )

    output_path = (
        execution_dir
        / "esp_upper_v7a_r1.out"
    )

    stderr_path = (
        execution_dir
        / "esp_upper_v7a_r1.stderr"
    )

    status_path = (
        execution_dir
        / "esp_upper_v7a_r1.exit_status"
    )

    gbw_path = (
        execution_dir
        / "esp_upper_v7a_r1.gbw"
    )

    densities_path = (
        execution_dir
        / "esp_upper_v7a_r1.densities"
    )

    property_path = (
        execution_dir
        / "esp_upper_v7a_r1.property.txt"
    )

    for path in (
        input_path,
        output_path,
        stderr_path,
        status_path,
    ):
        if not path.is_file():
            raise RuntimeError(
                f"Required execution artifact is missing: {path}"
            )

    text = output_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    stderr_text = stderr_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    try:
        shell_status = int(
            status_path.read_text(
                encoding="utf-8"
            ).strip()
        )
    except ValueError as error:
        raise RuntimeError(
            f"Invalid ORCA shell status: {status_path}"
        ) from error

    energy_matches = [
        float(match.group(1))
        for match in ENERGY_RE.finditer(text)
    ]

    scf_matches = [
        int(match.group(1))
        for match in SCF_RE.finditer(text)
    ]

    chelpg_rows, reported_total = (
        parse_chelpg_charges(text)
    )

    computed_total = (
        math.fsum(
            row["CHELPG_charge_e"]
            for row in chelpg_rows
        )
        if chelpg_rows
        else None
    )

    composition = Counter(
        row["element"]
        for row in chelpg_rows
    )

    fatal_markers = {
        name: bool(pattern.search(text))
        for name, pattern
        in FATAL_PATTERNS.items()
    }

    preflight_input_hash = (
        preflight.get(
            "summary",
            {},
        ).get("ESP_input_sha256")
    )

    # Older preflight schema may not expose the hash directly.
    # The copied input must at least match the authorized Day036 input.
    authorized_input = (
        DAY036_PREPARATION_DIR
        / "QM_F06_UPPER_V7A_R1_ESP_ORCA.inp"
    )

    gates = {
        "preflight_decision_matches": (
            preflight.get("decision")
            == EXPECTED_PREFLIGHT_DECISION
        ),
        "ESP_execution_authorized_upstream": (
            preflight.get(
                "authorizations",
                {},
            ).get(
                "ESP_execution_authorized"
            )
            is True
        ),
        "execution_input_exists": (
            input_path.is_file()
            and input_path.stat().st_size > 0
        ),
        "execution_input_matches_authorized_input": (
            authorized_input.is_file()
            and sha256(input_path)
            == sha256(authorized_input)
        ),
        "preflight_hash_consistent_if_reported": (
            preflight_input_hash is None
            or preflight_input_hash
            == sha256(input_path)
        ),
        "ORCA_shell_status_zero": (
            shell_status == 0
        ),
        "stderr_empty": (
            stderr_text == ""
        ),
        "normal_ORCA_termination_marker": (
            "ORCA TERMINATED NORMALLY"
            in text
        ),
        "SCF_converged": (
            len(scf_matches) >= 1
        ),
        "final_single_point_energy_available": (
            len(energy_matches) >= 1
        ),
        "no_fatal_markers": (
            not any(
                fatal_markers.values()
            )
        ),
        "CHELPG_charge_block_present": (
            "CHELPG Charges" in text
        ),
        "CHELPG_charge_count_52": (
            len(chelpg_rows)
            == EXPECTED_ATOM_COUNT
        ),
        "CHELPG_atom_indices_complete": (
            len(chelpg_rows)
            == EXPECTED_ATOM_COUNT
            and [
                row["atom_index_0based"]
                for row in chelpg_rows
            ]
            == list(
                range(EXPECTED_ATOM_COUNT)
            )
        ),
        "CHELPG_composition_B17_N14_H21": (
            dict(composition)
            == EXPECTED_COMPOSITION
        ),
        "CHELPG_computed_total_charge_zero": (
            computed_total is not None
            and abs(computed_total)
            <= TOTAL_CHARGE_TOLERANCE_E
        ),
        "CHELPG_reported_total_charge_zero": (
            reported_total is None
            or abs(reported_total)
            <= TOTAL_CHARGE_TOLERANCE_E
        ),
        "GBW_available": (
            gbw_path.is_file()
            and gbw_path.stat().st_size > 0
        ),
        "densities_available": (
            densities_path.is_file()
            and densities_path.stat().st_size > 0
        ),
        "property_file_available": (
            property_path.is_file()
            and property_path.stat().st_size > 0
        ),
    }

    passed = all(gates.values())

    decision = (
        "QM_F06_UPPER_V7A_R1_"
        "ESP_RESULT_AUDIT_PASS_"
        "RESP_DATA_EXTRACTION_DESIGN_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V7A_R1_"
        "ESP_RESULT_AUDIT_FAIL_"
        "RESP_WORKFLOW_BLOCKED"
    )

    write_csv(
        OUTPUT_CHARGES,
        chelpg_rows,
    )

    report = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "model": "QM_F06_UPPER_V7A_R1",
        "stage": "DAY036_ESP_RESULT_AUDIT",
        "execution_directory": str(
            execution_dir
        ),
        "decision": decision,
        "gates": gates,
        "summary": {
            "ORCA_shell_status": shell_status,
            "SCF_convergence_markers": len(
                scf_matches
            ),
            "last_SCF_cycle_count": (
                scf_matches[-1]
                if scf_matches
                else None
            ),
            "final_single_point_energy_Eh": (
                energy_matches[-1]
                if energy_matches
                else None
            ),
            "CHELPG_charge_count": len(
                chelpg_rows
            ),
            "CHELPG_computed_total_charge_e": (
                computed_total
            ),
            "CHELPG_reported_total_charge_e": (
                reported_total
            ),
            "composition": dict(
                sorted(composition.items())
            ),
            "fatal_markers": fatal_markers,
        },
        "files": {
            "input": {
                "path": str(input_path),
                "sha256": sha256(input_path),
            },
            "output": {
                "path": str(output_path),
                "sha256": sha256(output_path),
            },
            "GBW": {
                "path": str(gbw_path),
                "sha256": sha256(gbw_path),
            },
            "densities": {
                "path": str(densities_path),
                "sha256": sha256(densities_path),
            },
            "property_file": {
                "path": str(property_path),
                "sha256": sha256(property_path),
            },
            "CHELPG_charges": {
                "path": str(
                    OUTPUT_CHARGES
                ),
            },
        },
        "scientific_scope": {
            "CHELPG_charges_are_diagnostic": True,
            "raw_ESP_dataset_for_Amber_RESP_verified": False,
            "note": (
                "CHELPG is itself an ESP-fitted charge model. "
                "A separate validated extraction or conversion "
                "stage is required before AmberTools RESP input "
                "generation can be authorized."
            ),
        },
        "authorizations": {
            "RESP_data_extraction_design_authorized": (
                passed
            ),
            "RESP_input_generation_authorized": False,
            "RESP_execution_authorized": False,
            "RESP_validation_authorized": False,
            "charge_adoption_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    summary_text = f"""# QM_F06 UPPER V7-A R1 ESP Result Audit

## Decision

`{decision}`

## ORCA result

- Shell status: {shell_status}
- Normal termination: {"yes" if "ORCA TERMINATED NORMALLY" in text else "no"}
- SCF convergence cycles: {scf_matches[-1] if scf_matches else "unavailable"}
- Final energy: {energy_matches[-1] if energy_matches else "unavailable"} Eh
- CHELPG charge rows: {len(chelpg_rows)}
- Computed CHELPG total charge: {computed_total}
- Reported CHELPG total charge: {reported_total}

## Scientific interpretation

The ORCA CHELPG calculation is accepted as a completed and internally
consistent electronic-property calculation.

CHELPG charges are retained as diagnostic ESP-fitted charges. This
result does not by itself prove that a raw ESP dataset compatible with
AmberTools RESP has been generated.

The next authorized activity is the design and validation of the
ORCA-to-RESP ESP-data extraction or conversion stage.

RESP input generation, RESP execution, charge adoption, force-field
adoption and molecular dynamics remain blocked.
"""

    OUTPUT_SUMMARY.write_text(
        summary_text,
        encoding="utf-8",
    )

    print("=" * 108)
    print("QM_F06 UPPER V7-A R1 ESP RESULT AUDIT")
    print("=" * 108)

    for name, value in gates.items():
        print(
            f"{name:62s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(
        "Final single-point energy Eh:",
        energy_matches[-1]
        if energy_matches
        else None,
    )
    print(
        "SCF converged after cycles:",
        scf_matches[-1]
        if scf_matches
        else None,
    )
    print(
        "CHELPG charge count:",
        len(chelpg_rows),
    )
    print(
        "CHELPG computed total charge e:",
        computed_total,
    )
    print("Decision:", decision)
    print("Report:", OUTPUT_REPORT)
    print("Charges:", OUTPUT_CHARGES)
    print(
        "RESP data-extraction design authorized:",
        passed,
    )
    print("RESP execution authorized: False")
    print("Force-field adoption authorized: False")
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
