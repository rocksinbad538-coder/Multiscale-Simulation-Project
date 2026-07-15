#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "07_r1_frozen_solute_nvt_20ps_preparation"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_grompp.log"
)

MDOUT = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_mdout.mdp"
)

TPR = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.tpr"
)

PREPARATION_SUMMARY = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_preparation_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_nvt_grompp_warning_authorization_gates.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_nvt_grompp_warning_authorization_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_NVT_GROMPP_WARNING_AUTHORIZATION_DAY023.md"
)

EXPECTED_WARNING = (
    "Some temperature coupling groups do not use "
    "temperature coupling. We will assume their "
    "temperature is not more than 300.000 K. If their "
    "temperature is higher, the energy error and the "
    "Verlet buffer might be underestimated."
)


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def normalize(text: str) -> str:
    return " ".join(
        text.split()
    )


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; "
            f"found {len(rows)}"
        )

    return rows[0]


def parse_mdp(
    path: Path,
) -> dict[str, str]:
    require_file(path)

    result: dict[str, str] = {}

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        active = raw_line.split(
            ";",
            1,
        )[0].strip()

        if (
            not active
            or active.startswith("#")
            or "=" not in active
        ):
            continue

        key, value = active.split(
            "=",
            1,
        )

        normalized_key = (
            key.strip()
            .lower()
            .replace("-", "_")
        )

        result[
            normalized_key
        ] = " ".join(
            value.strip().split()
        )

    return result


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    key: row.get(
                        key,
                        "",
                    )
                    for key in fields
                }
            )


def parse_dof(
    text: str,
    group: str,
) -> float | None:
    match = re.search(
        (
            r"Number of degrees of freedom in "
            r"T-Coupling group "
            + re.escape(group)
            + r" is\s+([-+0-9.eE]+)"
        ),
        text,
    )

    if match is None:
        return None

    return float(
        match.group(1)
    )


def main() -> None:
    for required in (
        GROMPP_LOG,
        MDOUT,
        TPR,
        PREPARATION_SUMMARY,
    ):
        require_file(required)

    log_text = GROMPP_LOG.read_text(
        encoding="utf-8",
        errors="replace",
    )

    normalized_log = normalize(
        log_text
    )

    mdp = parse_mdp(
        MDOUT
    )

    preparation = read_single_csv_row(
        PREPARATION_SUMMARY
    )

    warning_count = len(
        re.findall(
            r"(?m)^WARNING\s+\d+\s+\[",
            log_text,
        )
    )

    expected_warning_present = (
        normalize(
            EXPECTED_WARNING
        )
        in normalized_log
    )

    fatal_error_present = bool(
        re.search(
            r"(?mi)^Fatal error:",
            log_text,
        )
    )

    maxwarn_one_present = bool(
        re.search(
            r"(?:^|\s)-maxwarn\s+1(?:\s|$)",
            log_text,
        )
    )

    capl_unbound_note = (
        "In moleculetype 'CAPL' 163 atoms are not bound"
        in log_text
    )

    capu_unbound_note = (
        "In moleculetype 'CAPU' 163 atoms are not bound"
        in log_text
    )

    vcm_note_present = (
        "2110 atoms are not part of any of the VCM groups"
        in log_text
    )

    dof_sol = parse_dof(
        log_text,
        "SOL",
    )

    dof_hbn_pyr = parse_dof(
        log_text,
        "HBN_PYR",
    )

    dof_caps = parse_dof(
        log_text,
        "CAPS",
    )

    tc_groups = mdp.get(
        "tc_grps",
        "",
    ).split()

    tau_t = [
        float(value)
        for value in mdp.get(
            "tau_t",
            "",
        ).split()
    ]

    ref_t = [
        float(value)
        for value in mdp.get(
            "ref_t",
            "",
        ).split()
    ]

    freeze_groups = mdp.get(
        "freezegrps",
        "",
    ).split()

    freeze_dimensions = mdp.get(
        "freezedim",
        "",
    ).split()

    gates = {
        "TPR_exists_and_is_nonempty": (
            TPR.exists()
            and TPR.stat().st_size > 0
        ),
        "exactly_one_grompp_warning": (
            warning_count == 1
        ),
        "warning_text_is_expected": (
            expected_warning_present
        ),
        "no_fatal_error_after_authorization": (
            not fatal_error_present
        ),
        "grompp_used_exactly_maxwarn_1": (
            maxwarn_one_present
        ),
        "temperature_groups_are_correct": (
            tc_groups
            == [
                "SOL",
                "HBN_PYR",
                "CAPS",
            ]
        ),
        "tau_t_values_are_correct": (
            tau_t
            == [
                0.1,
                -1.0,
                -1.0,
            ]
        ),
        "reference_temperatures_are_300K": (
            ref_t
            == [
                300.0,
                300.0,
                300.0,
            ]
        ),
        "freeze_groups_are_correct": (
            freeze_groups
            == [
                "HBN_PYR",
                "CAPS",
            ]
        ),
        "all_six_freeze_dimensions_are_active": (
            freeze_dimensions
            == [
                "Y",
                "Y",
                "Y",
                "Y",
                "Y",
                "Y",
            ]
        ),
        "SOL_has_positive_degrees_of_freedom": (
            dof_sol is not None
            and dof_sol > 0.0
        ),
        "HBN_PYR_has_zero_degrees_of_freedom": (
            dof_hbn_pyr == 0.0
        ),
        "CAPS_has_zero_degrees_of_freedom": (
            dof_caps == 0.0
        ),
        "expected_CAPL_unbound_note_is_present": (
            capl_unbound_note
        ),
        "expected_CAPU_unbound_note_is_present": (
            capu_unbound_note
        ),
        "expected_VCM_note_is_present": (
            vcm_note_present
        ),
        "preparation_decision_is_pass": (
            preparation.get(
                "decision"
            )
            ==
            "R1_FROZEN_SOLUTE_NVT_20PS_PREPARED"
        ),
        "preparation_authorized_execution": (
            parse_bool(
                preparation.get(
                    "NVT_execution_authorized",
                    "false",
                )
            )
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    authorized = (
        len(failed_gates) == 0
    )

    decision = (
        "R1_NVT_GROMPP_WARNING_AUTHORIZED"
        if authorized
        else
        "R1_NVT_GROMPP_WARNING_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R1_FROZEN_SOLUTE_NVT_20PS"
        if authorized
        else
        "RESOLVE_R1_GROMPP_WARNING_AUDIT_FAILURES"
    )

    gate_rows = [
        {
            "gate": name,
            "pass": passed,
        }
        for name, passed in gates.items()
    ]

    write_csv(
        GATE_CSV,
        gate_rows,
    )

    summary = {
        "decision": decision,
        "grompp_warning_count": (
            warning_count
        ),
        "expected_warning_present": (
            expected_warning_present
        ),
        "fatal_error_present": (
            fatal_error_present
        ),
        "maxwarn_one_present": (
            maxwarn_one_present
        ),
        "SOL_degrees_of_freedom": (
            dof_sol
        ),
        "HBN_PYR_degrees_of_freedom": (
            dof_hbn_pyr
        ),
        "CAPS_degrees_of_freedom": (
            dof_caps
        ),
        "temperature_groups": (
            " ".join(tc_groups)
        ),
        "tau_t_ps": (
            " ".join(
                str(value)
                for value in tau_t
            )
        ),
        "reference_temperatures_K": (
            " ".join(
                str(value)
                for value in ref_t
            )
        ),
        "freeze_groups": (
            " ".join(
                freeze_groups
            )
        ),
        "freeze_dimensions": (
            " ".join(
                freeze_dimensions
            )
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "NVT_execution_authorized": (
            authorized
        ),
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R1 NVT Grompp Warning Authorization

## Warning reviewed

The only GROMACS warning is:

> {EXPECTED_WARNING}

This warning is authorized only for this R1 frozen-solute screening
protocol because:

- HBN_PYR has {dof_hbn_pyr} degrees of freedom;
- CAPS has {dof_caps} degrees of freedom;
- both groups are frozen in all Cartesian dimensions;
- both groups use tau-t = -1 and are not thermostatted;
- SOL is the only mobile and thermostatted group;
- no additional warning or fatal error is present.

The use of `-maxwarn 1` is therefore restricted to this exact,
programmatically audited warning.

## Coupling state

- T-coupling groups:
  **{' '.join(tc_groups)}**
- tau-t:
  **{' '.join(str(value) for value in tau_t)} ps**
- ref-t:
  **{' '.join(str(value) for value in ref_t)} K**
- Freeze groups:
  **{' '.join(freeze_groups)}**
- Freeze dimensions:
  **{' '.join(freeze_dimensions)}**

## Degrees of freedom

- SOL:
  **{dof_sol}**
- HBN_PYR:
  **{dof_hbn_pyr}**
- CAPS:
  **{dof_caps}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- NVT execution authorized:
  **{'YES' if authorized else 'NO'}**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 NVT grompp-warning audit completed."
    )

    print(
        "Warning count / expected warning: "
        f"{warning_count} / "
        f"{'YES' if expected_warning_present else 'NO'}"
    )

    print(
        "Fatal error present: "
        f"{'YES' if fatal_error_present else 'NO'}"
    )

    print(
        "T-coupling groups: "
        + " ".join(tc_groups)
    )

    print(
        "tau-t values: "
        + " ".join(
            str(value)
            for value in tau_t
        )
        + " ps"
    )

    print(
        "Degrees of freedom SOL/HBN_PYR/CAPS: "
        f"{dof_sol}/"
        f"{dof_hbn_pyr}/"
        f"{dof_caps}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(
                failed_gates
            )
        )
    )

    print(
        "NVT execution authorized: "
        f"{'YES' if authorized else 'NO'}"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(GATE_CSV)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not authorized:
        raise RuntimeError(
            "The R1 grompp warning could not be authorized."
        )


if __name__ == "__main__":
    main()
