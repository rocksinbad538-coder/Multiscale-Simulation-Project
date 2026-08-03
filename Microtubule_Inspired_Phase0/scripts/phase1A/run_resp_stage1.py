#!/usr/bin/env python3
"""
DAY038 / D038-D1

Reproducible AmberTools RESP Stage 1 execution for
QM_F06_UPPER_V7A_R1.

Scientific policy
-----------------
- Uses the authorized ORCA-derived Amber ESP dataset.
- Preserves the adopted 52-atom ordering.
- Uses total molecular charge = 0.
- Enforces no non-singleton atom-charge equivalences.
- Uses the standard RESP Stage 1 hyperbolic restraint weight:
  qwt = 0.0005.
- Does not adopt fitted charges.
- Does not execute RESP Stage 2.
- Does not infer connectivity or GAFF atom types.

Outputs are isolated in a timestamped execution directory.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from resp_common import (
    load_json,
    require_file,
    sha256,
)


ROOT = Path(__file__).resolve().parents[2]

DAY036_PREP = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_resp_preparation"
)

DAY038_GENERATION = (
    ROOT
    / "runs/phase1A/day038_resp_generation"
)

EXECUTION_PARENT = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
)

PREFLIGHT = (
    DAY038_GENERATION
    / "DAY038_RESP_PREFLIGHT.json"
)

AMBER_ESP = (
    DAY038_GENERATION
    / "candidate_from_orca_vpot.esp"
)

ATOM_CLASSES = (
    DAY036_PREP
    / "QM_F06_UPPER_V7A_R1_RESP_ATOM_CLASSES.csv"
)

EQUIVALENCE_GROUPS = (
    DAY036_PREP
    / "QM_F06_UPPER_V7A_R1_RESP_EQUIVALENCE_GROUPS.csv"
)

RESP_PREPARATION = (
    DAY036_PREP
    / "QM_F06_UPPER_V7A_R1_RESP_PREPARATION.json"
)

RESP_PROTOCOL = (
    DAY036_PREP
    / "QM_F06_UPPER_V7A_R1_RESP_PROTOCOL.md"
)

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

AUTHORIZED_AMBER_ESP_SHA256 = (
    "3a69a39b7848efce3d2b2467d8b7f28f"
    "dc1925eacca479a5c9a39019b741aac6"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_GRID_POINT_COUNT = 24835
EXPECTED_NET_CHARGE = 0
RESP_STAGE1_QWT = 0.0005

ATOMIC_NUMBERS = {
    "H": 1,
    "B": 5,
    "N": 7,
}


def utc_compact() -> str:
    return datetime.now(timezone.utc).strftime(
        "%Y%m%dT%H%M%SZ"
    )


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def require_true(value: object, label: str) -> None:
    if value is not True:
        raise RuntimeError(
            f"Required gate is not true: {label}"
        )


def parse_qout(path: Path) -> list[float]:
    require_file(path)

    values: list[float] = []

    for token in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).split():
        try:
            value = float(
                token.replace("D", "E").replace("d", "e")
            )
        except ValueError:
            continue

        if not math.isfinite(value):
            raise RuntimeError(
                f"Non-finite value in RESP qout: {token!r}"
            )

        values.append(value)

    return values


def build_resp1_input(
    elements: list[str],
) -> str:
    lines = [
        "QM_F06_UPPER_V7A_R1 RESP Stage 1",
        " &cntrl",
        "  nmol = 1,",
        "  ihfree = 1,",
        "  iqopt = 1,",
        f"  qwt = {RESP_STAGE1_QWT:.7f},",
        " /",
        "    1.0",
        "QM_F06_UPPER_V7A_R1 neutral 52-atom model",
        (
            f"{EXPECTED_NET_CHARGE:5d}"
            f"{len(elements):5d}"
        ),
    ]

    for element in elements:
        atomic_number = ATOMIC_NUMBERS[element]

        # ivary = 0:
        # each atom is an independently fitted variable.
        lines.append(
            f"{atomic_number:5d}{0:5d}"
        )

    # Number of additional group-charge constraints.
    # Day036 authorizes none, so RESP requires a terminal zero.
    lines.append(f"{0:5d}")

    return "\n".join(lines) + "\n"


print("=" * 100)
print("DAY038 / D038-D1 — REPRODUCIBLE RESP STAGE 1")
print("=" * 100)


print("\n[1] REQUIRED INPUTS")

for path in (
    PREFLIGHT,
    AMBER_ESP,
    ATOM_CLASSES,
    EQUIVALENCE_GROUPS,
    RESP_PREPARATION,
    RESP_PROTOCOL,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] AMBERTOOLS EXECUTABLE")

resp_executable = shutil.which("resp")

if resp_executable is None:
    raise RuntimeError(
        "AmberTools RESP executable was not found in PATH"
    )

resp_executable_path = Path(
    resp_executable
).resolve()

print(f"resp = {resp_executable_path}")


print("\n[3] DAY038 PREFLIGHT")

preflight = load_json(PREFLIGHT)

if preflight.get("decision") != "D038_RESP_PREFLIGHT_PASS":
    raise RuntimeError(
        "Unexpected Day038 RESP preflight decision:\n"
        f"{preflight.get('decision')!r}"
    )

gates = preflight.get("gates", {})

for gate_name in (
    "authorized_vpot",
    "atom_count_52",
    "grid_point_count_24835",
    "amber_matches_vpot",
):
    require_true(
        gates.get(gate_name),
        f"preflight.gates.{gate_name}",
    )

if (
    preflight["vpot"]["sha256"]
    != AUTHORIZED_VPOT_SHA256
):
    raise RuntimeError(
        "Authorized VPOT hash mismatch in preflight"
    )

observed_esp_sha256 = sha256(AMBER_ESP)

if (
    observed_esp_sha256
    != AUTHORIZED_AMBER_ESP_SHA256
):
    raise RuntimeError(
        "Amber ESP SHA256 mismatch.\n"
        f"Expected: {AUTHORIZED_AMBER_ESP_SHA256}\n"
        f"Observed: {observed_esp_sha256}"
    )

if (
    preflight["amber_esp"]["sha256"]
    != observed_esp_sha256
):
    raise RuntimeError(
        "Amber ESP hash does not match the Day038 preflight"
    )

print("preflight_decision_gate = PASS")
print("authorized_vpot_gate    = PASS")
print("authorized_esp_gate     = PASS")


print("\n[4] DAY036 PREPARATION POLICY")

preparation = load_json(RESP_PREPARATION)

classification = preparation.get(
    "classification",
    {},
)

if (
    classification.get(
        "enforced_nonsingleton_equivalence_groups"
    )
    != 0
):
    raise RuntimeError(
        "Day036 policy does not authorize a fully "
        "unconstrained-equivalence Stage 1 fit"
    )

electronic_structure = preparation.get(
    "electronic_structure",
    {},
)

if (
    electronic_structure.get("net_charge")
    != EXPECTED_NET_CHARGE
):
    raise RuntimeError(
        "Unexpected net charge in Day036 preparation"
    )

preparation_gates = preparation.get("gates", {})

require_true(
    preparation_gates.get(
        "map_element_order_matches_XYZ"
    ),
    "Day036 map element order",
)

require_true(
    preparation_gates.get(
        "charge_and_multiplicity_defined"
    ),
    "Day036 charge and multiplicity",
)

print("net_charge_zero_gate              = PASS")
print("atom_order_gate                   = PASS")
print("no_enforced_equivalence_gate      = PASS")


print("\n[5] ATOM ORDER AND ELEMENTS")

atom_rows = read_csv_rows(ATOM_CLASSES)

if len(atom_rows) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        "RESP atom-class row count mismatch.\n"
        f"Expected: {EXPECTED_ATOM_COUNT}\n"
        f"Observed: {len(atom_rows)}"
    )

observed_indices = [
    int(row["atom_index_0based"])
    for row in atom_rows
]

if observed_indices != list(
    range(EXPECTED_ATOM_COUNT)
):
    raise RuntimeError(
        "RESP atom-class indices are not exactly 0..51"
    )

elements = [
    row["element"].strip()
    for row in atom_rows
]

unsupported_elements = sorted(
    set(elements) - set(ATOMIC_NUMBERS)
)

if unsupported_elements:
    raise RuntimeError(
        "Unsupported elements in RESP input: "
        f"{unsupported_elements}"
    )

composition = {
    element: elements.count(element)
    for element in sorted(set(elements))
}

expected_composition = {
    "B": 17,
    "H": 21,
    "N": 14,
}

if composition != expected_composition:
    raise RuntimeError(
        "RESP atom composition mismatch.\n"
        f"Expected: {expected_composition}\n"
        f"Observed: {composition}"
    )

print(f"atom_count  = {len(elements)}")
print(f"composition = {composition}")
print("element_order_gate = PASS")


print("\n[6] EQUIVALENCE POLICY")

equivalence_rows = read_csv_rows(
    EQUIVALENCE_GROUPS
)

if len(equivalence_rows) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        "RESP equivalence row count mismatch"
    )

enforced_rows = [
    row
    for row in equivalence_rows
    if row["equivalence_enforced"].strip().lower()
    in {"true", "1", "yes"}
]

if enforced_rows:
    raise RuntimeError(
        "One or more atom equivalences are enforced, "
        "contrary to the Day036 policy"
    )

group_ids = [
    row["enforced_group_id"].strip()
    for row in equivalence_rows
]

if len(set(group_ids)) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        "Expected one singleton enforced group per atom"
    )

print(f"equivalence_rows       = {len(equivalence_rows)}")
print(f"enforced_rows          = {len(enforced_rows)}")
print(f"unique_enforced_groups = {len(set(group_ids))}")
print("equivalence_policy_gate = PASS")


print("\n[7] CREATE ISOLATED EXECUTION")

timestamp = utc_compact()

execution_dir = (
    EXECUTION_PARENT
    / f"resp_stage1_upper_v7a_r1_{timestamp}"
)

if execution_dir.exists():
    raise RuntimeError(
        f"Execution directory already exists: {execution_dir}"
    )

execution_dir.mkdir(
    parents=True,
    exist_ok=False,
)

esp_copy = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP.esp"
)

resp_input = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_resp1.in"
)

resp_output = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_resp1.out"
)

resp_punch = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_resp1.pch"
)

resp_qout = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_resp1.qout"
)

stdout_path = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_resp1.stdout"
)

stderr_path = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_resp1.stderr"
)

report_path = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1.json"
)

shutil.copy2(
    AMBER_ESP,
    esp_copy,
)

resp_input.write_text(
    build_resp1_input(elements),
    encoding="utf-8",
)

if sha256(esp_copy) != observed_esp_sha256:
    raise RuntimeError(
        "Copied ESP SHA256 mismatch"
    )

print(f"execution_dir = {execution_dir}")
print(f"resp_input    = {resp_input}")
print(f"esp_copy      = {esp_copy}")


print("\n[8] RESP COMMAND")

command = [
    str(resp_executable_path),
    "-O",
    "-i",
    resp_input.name,
    "-o",
    resp_output.name,
    "-p",
    resp_punch.name,
    "-t",
    resp_qout.name,
    "-e",
    esp_copy.name,
]

print(" ".join(command))


print("\n[9] EXECUTE RESP STAGE 1")

result = subprocess.run(
    command,
    cwd=execution_dir,
    text=True,
    capture_output=True,
)

stdout_path.write_text(
    result.stdout,
    encoding="utf-8",
)

stderr_path.write_text(
    result.stderr,
    encoding="utf-8",
)

print(f"returncode = {result.returncode}")
print(f"stdout_bytes = {stdout_path.stat().st_size}")
print(f"stderr_bytes = {stderr_path.stat().st_size}")


print("\n[10] OUTPUT INVENTORY")

output_paths = (
    resp_output,
    resp_punch,
    resp_qout,
    stdout_path,
    stderr_path,
)

output_inventory = {}

for path in output_paths:
    exists = path.is_file()
    bytes_count = (
        path.stat().st_size
        if exists
        else 0
    )

    output_inventory[path.name] = {
        "exists": exists,
        "bytes": bytes_count,
        "sha256": (
            sha256(path)
            if exists and bytes_count > 0
            else None
        ),
    }

    print(
        f"{'FOUND' if exists else 'MISSING'}  "
        f"bytes={bytes_count:>10}  {path.name}"
    )


print("\n[11] CANDIDATE CHARGE PARSE")

candidate_charges: list[float] = []
charge_sum = None
charge_minimum = None
charge_maximum = None

if (
    resp_qout.is_file()
    and resp_qout.stat().st_size > 0
):
    candidate_charges = parse_qout(
        resp_qout
    )

    if len(candidate_charges) == EXPECTED_ATOM_COUNT:
        charge_sum = sum(candidate_charges)
        charge_minimum = min(candidate_charges)
        charge_maximum = max(candidate_charges)

print(
    f"candidate_charge_count = "
    f"{len(candidate_charges)}"
)
print(f"candidate_charge_sum   = {charge_sum}")
print(f"candidate_charge_min   = {charge_minimum}")
print(f"candidate_charge_max   = {charge_maximum}")


execution_process_gate = (
    result.returncode == 0
)

required_outputs_gate = all(
    path.is_file() and path.stat().st_size > 0
    for path in (
        resp_output,
        resp_punch,
        resp_qout,
    )
)

charge_count_gate = (
    len(candidate_charges)
    == EXPECTED_ATOM_COUNT
)

charge_sum_gate = (
    charge_sum is not None
    and abs(
        charge_sum - EXPECTED_NET_CHARGE
    )
    <= 5.0e-5
)


print("\n[12] EXECUTION GATES")
print(
    "execution_process_gate = "
    + (
        "PASS"
        if execution_process_gate
        else "FAIL"
    )
)
print(
    "required_outputs_gate  = "
    + (
        "PASS"
        if required_outputs_gate
        else "FAIL"
    )
)
print(
    "charge_count_gate      = "
    + (
        "PASS"
        if charge_count_gate
        else "FAIL"
    )
)
print(
    "charge_sum_gate        = "
    + (
        "PASS"
        if charge_sum_gate
        else "FAIL"
    )
)


all_execution_gates_pass = (
    execution_process_gate
    and required_outputs_gate
    and charge_count_gate
    and charge_sum_gate
)

decision = (
    "D038_D1_RESP_STAGE1_EXECUTION_PASS_"
    "CANDIDATE_CHARGES_GENERATED_NOT_ADOPTED"
    if all_execution_gates_pass
    else
    "D038_D1_RESP_STAGE1_EXECUTION_REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "execution_directory": str(
        execution_dir.resolve()
    ),
    "command": command,
    "returncode": result.returncode,
    "protocol": {
        "stage": "RESP_STAGE1",
        "nmol": 1,
        "ihfree": 1,
        "iqopt": 1,
        "qwt": RESP_STAGE1_QWT,
        "net_charge": EXPECTED_NET_CHARGE,
        "atom_count": EXPECTED_ATOM_COUNT,
        "grid_point_count": (
            EXPECTED_GRID_POINT_COUNT
        ),
        "equivalence_policy": (
            "NO_NONSINGLETON_EQUIVALENCE_ENFORCED"
        ),
    },
    "source_identity": {
        "day038_preflight": str(
            PREFLIGHT.resolve()
        ),
        "day038_preflight_sha256": sha256(
            PREFLIGHT
        ),
        "source_amber_esp": str(
            AMBER_ESP.resolve()
        ),
        "source_amber_esp_sha256": (
            observed_esp_sha256
        ),
        "execution_amber_esp_sha256": (
            sha256(esp_copy)
        ),
        "atom_classes_sha256": sha256(
            ATOM_CLASSES
        ),
        "equivalence_groups_sha256": sha256(
            EQUIVALENCE_GROUPS
        ),
        "resp_preparation_sha256": sha256(
            RESP_PREPARATION
        ),
        "resp_protocol_sha256": sha256(
            RESP_PROTOCOL
        ),
        "resp_executable": str(
            resp_executable_path
        ),
        "resp_executable_sha256": sha256(
            resp_executable_path
        ),
    },
    "composition": composition,
    "candidate_charges": {
        "count": len(candidate_charges),
        "sum_e": charge_sum,
        "minimum_e": charge_minimum,
        "maximum_e": charge_maximum,
        "values_e": candidate_charges,
    },
    "outputs": output_inventory,
    "gates": {
        "preflight_decision_gate": True,
        "authorized_vpot_gate": True,
        "authorized_esp_gate": True,
        "net_charge_zero_gate": True,
        "atom_order_gate": True,
        "element_order_gate": True,
        "no_enforced_equivalence_gate": True,
        "execution_process_gate": (
            execution_process_gate
        ),
        "required_outputs_gate": (
            required_outputs_gate
        ),
        "charge_count_gate": (
            charge_count_gate
        ),
        "charge_sum_gate": (
            charge_sum_gate
        ),
    },
    "authorizations": {
        "candidate_charge_analysis_authorized": (
            all_execution_gates_pass
        ),
        "RESP_stage2_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
}

report_path.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

latest_pointer = (
    EXECUTION_PARENT
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

latest_pointer.write_text(
    str(execution_dir.relative_to(ROOT))
    + "\n",
    encoding="utf-8",
)


print("\n[13] REPORT")
print(f"report_path = {report_path}")
print(f"report_sha256 = {sha256(report_path)}")
print(f"latest_pointer = {latest_pointer}")


print("\n[14] D038-D1 DECISION")
print(f"decision = {decision}")
print(
    "RESP Stage 2 remains blocked."
)
print(
    "Charge adoption remains blocked."
)
print("=" * 100)


# D038_DRIVER_RETURNS_NONZERO_ON_FAILED_GATES
if not all_execution_gates_pass:
    raise SystemExit(2)
