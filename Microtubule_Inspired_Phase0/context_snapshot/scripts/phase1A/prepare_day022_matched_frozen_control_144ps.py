#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

EXECUTION = PROTOCOL / "execution"

CONTROL_ROOT = (
    PROTOCOL
    / "matched_frozen_control_144ps"
)

INPUT_ROOT = (
    CONTROL_ROOT
    / "prepared_input"
)

STATIC_ROOT = (
    CONTROL_ROOT
    / "static_validation"
)

SOURCE_STAGE = "02_nvt_k10000_1ps"

SOURCE_TPR = (
    EXECUTION
    / SOURCE_STAGE
    / f"{SOURCE_STAGE}.tpr"
)

SOURCE_MDP = (
    PROTOCOL
    / "protocol_inputs/mdp"
    / f"{SOURCE_STAGE}.mdp"
)

START_GRO = (
    EXECUTION
    / "01_em_k10000"
    / "01_em_k10000.gro"
)

TOP = (
    PROTOCOL
    / "protocol_inputs/topology"
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    PROTOCOL
    / "protocol_inputs"
    / "mobile_release_index.ndx"
)

AUDIT_MODULE_PATH = (
    ROOT
    / "scripts/phase1A/"
    "audit_day022_matched_frozen_control_inputs.py"
)

TARGET_STAGE = "matched_frozen_control_144ps"

TARGET_MDP = (
    INPUT_ROOT
    / f"{TARGET_STAGE}.mdp"
)

TARGET_TPR = (
    STATIC_ROOT
    / f"{TARGET_STAGE}.tpr"
)

TARGET_MDOUT = (
    STATIC_ROOT
    / f"{TARGET_STAGE}_mdout.mdp"
)

TARGET_PROCESSED_TOP = (
    STATIC_ROOT
    / f"{TARGET_STAGE}_processed.top"
)

GROMPP_LOG = (
    STATIC_ROOT
    / f"{TARGET_STAGE}_grompp.log"
)

SOURCE_DUMP_STDOUT = (
    STATIC_ROOT
    / "source_stage02_tpr_dump_stdout.txt"
)

SOURCE_DUMP_STDERR = (
    STATIC_ROOT
    / "source_stage02_tpr_dump_stderr.txt"
)

CONTROL_DUMP_STDOUT = (
    STATIC_ROOT
    / "matched_frozen_control_tpr_dump_stdout.txt"
)

CONTROL_DUMP_STDERR = (
    STATIC_ROOT
    / "matched_frozen_control_tpr_dump_stderr.txt"
)

SUMMARY_CSV = (
    STATIC_ROOT
    / "matched_frozen_control_144ps_static_summary.csv"
)

REPORT_MD = (
    STATIC_ROOT
    / "MATCHED_FROZEN_CONTROL_144PS_STATIC_AUTHORIZATION_DAY022.md"
)

EXPECTED_SOURCE_WATER_HASH = (
    "d11ac84fa2d4a2a9a594a91ac0a6e0714"
    "dc3caa9da33f0ed12617371feff722d"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26

SOLUTE_COUNT = (
    HBN_COUNT
    + PYR_COUNT * PYR_ATOMS
)

TOTAL_ATOMS = 68320

DT_PS = 0.0005
NSTEPS = 288000
DURATION_PS = 144.0
XTC_INTERVAL_STEPS = 1000
EXPECTED_XTC_FRAMES = 289


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


def locate_gmx() -> str:
    executable = shutil.which("gmx")

    if executable:
        return executable

    fallback = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if fallback.exists():
        return str(fallback)

    raise RuntimeError(
        "Could not locate GROMACS"
    )


def load_audit_module():
    specification = (
        importlib.util.spec_from_file_location(
            "matched_control_input_audit_module",
            AUDIT_MODULE_PATH,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            f"Could not load audit module: "
            f"{AUDIT_MODULE_PATH}"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(module)

    return module


def canonical_key(value: str) -> str:
    return (
        value.strip()
        .lower()
        .replace("_", "-")
    )


def parse_mdp(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}

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
            or "=" not in active
        ):
            continue

        key, value = active.split(
            "=",
            1,
        )

        values[
            canonical_key(key)
        ] = value.strip()

    return values


def rewrite_mdp() -> None:
    source_text = SOURCE_MDP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    overrides = {
        "integrator": "md",
        "dt": "0.0005",
        "nsteps": "288000",
        "continuation": "no",
        "gen-vel": "yes",
        "gen-seed": "20260706",
        "freezegrps": "HBN_PYR",
        "freezedim": "Y Y Y",
        "nstxout": "0",
        "nstvout": "0",
        "nstfout": "0",
        "nstxout-compressed": "1000",
        "nstcalcenergy": "100",
        "nstenergy": "200",
        "nstlog": "1000",
    }

    display_keys = {
        "integrator": "integrator",
        "dt": "dt",
        "nsteps": "nsteps",
        "continuation": "continuation",
        "gen-vel": "gen-vel",
        "gen-seed": "gen-seed",
        "freezegrps": "freezegrps",
        "freezedim": "freezedim",
        "nstxout": "nstxout",
        "nstvout": "nstvout",
        "nstfout": "nstfout",
        "nstxout-compressed": (
            "nstxout-compressed"
        ),
        "nstcalcenergy": "nstcalcenergy",
        "nstenergy": "nstenergy",
        "nstlog": "nstlog",
    }

    used: set[str] = set()
    output_lines: list[str] = []

    for raw_line in source_text.splitlines():
        active = raw_line.split(
            ";",
            1,
        )[0].strip()

        if (
            not active
            or "=" not in active
        ):
            output_lines.append(
                raw_line
            )
            continue

        key, _ = active.split(
            "=",
            1,
        )

        canonical = canonical_key(
            key
        )

        if canonical == "define":
            continue

        if canonical in overrides:
            if canonical not in used:
                output_lines.append(
                    f"{display_keys[canonical]} = "
                    f"{overrides[canonical]}"
                )

                used.add(
                    canonical
                )

            continue

        output_lines.append(
            raw_line
        )

    for canonical, value in overrides.items():
        if canonical in used:
            continue

        output_lines.append(
            f"{display_keys[canonical]} = {value}"
        )

    output_lines.extend(
        [
            "",
            (
                "; Day022 matched frozen control"
            ),
            (
                "; Same initial coordinates and "
                "velocity-generation seed as Stage02"
            ),
            (
                "; HBN and PYR frozen in all "
                "Cartesian dimensions"
            ),
            (
                "; 288000 steps x 0.0005 ps = "
                "144 ps"
            ),
            (
                "; compressed coordinates every "
                "0.5 ps"
            ),
            (
                "; expected frames including "
                "t=0 and t=144 ps: 289"
            ),
        ]
    )

    TARGET_MDP.write_text(
        "\n".join(
            output_lines
        )
        + "\n",
        encoding="utf-8",
    )


def count_position_restraints(
    path: Path,
) -> int:
    section = ""
    count = 0

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        active = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not active:
            continue

        if (
            active.startswith("[")
            and active.endswith("]")
        ):
            section = (
                active[1:-1]
                .strip()
                .lower()
            )
            continue

        if section != "position_restraints":
            continue

        fields = active.split()

        try:
            int(fields[0])
        except (
            IndexError,
            ValueError,
        ):
            continue

        count += 1

    return count


def run_grompp(
    gmx: str,
) -> int:
    command = [
        gmx,
        "grompp",
        "-f",
        str(TARGET_MDP),
        "-c",
        str(START_GRO),
        "-r",
        str(START_GRO),
        "-p",
        str(TOP),
        "-n",
        str(NDX),
        "-o",
        str(TARGET_TPR),
        "-po",
        str(TARGET_MDOUT),
        "-pp",
        str(TARGET_PROCESSED_TOP),
        "-maxwarn",
        "0",
    ]

    completed = subprocess.run(
        command,
        cwd=TOP.parent,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    GROMPP_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    return completed.returncode


def dump_tpr(
    gmx: str,
    tpr: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> str:
    completed = subprocess.run(
        [
            gmx,
            "dump",
            "-s",
            str(tpr),
        ],
        cwd=tpr.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    stdout_path.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    stderr_path.write_text(
        completed.stderr,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            f"gmx dump failed for {tpr}. "
            f"See {stderr_path}"
        )

    if not completed.stdout.strip():
        raise RuntimeError(
            f"gmx dump returned empty stdout "
            f"for {tpr}"
        )

    return completed.stdout


def nonzero_fraction(
    velocities: np.ndarray,
) -> float:
    nonzero = np.any(
        np.abs(
            velocities
        )
        > 0.0,
        axis=1,
    )

    return float(
        np.mean(
            nonzero
        )
    )


def write_csv(
    row: dict[str, object],
) -> None:
    with SUMMARY_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                row.keys()
            ),
        )

        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    for required in (
        SOURCE_TPR,
        SOURCE_MDP,
        START_GRO,
        TOP,
        NDX,
        AUDIT_MODULE_PATH,
    ):
        require_file(
            required
        )

    INPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    STATIC_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    gmx = locate_gmx()
    audit = load_audit_module()

    rewrite_mdp()

    target_input = parse_mdp(
        TARGET_MDP
    )

    preparation_failures: list[str] = []

    required_input_values = {
        "integrator": "md",
        "dt": "0.0005",
        "nsteps": "288000",
        "continuation": "no",
        "gen-vel": "yes",
        "gen-seed": "20260706",
        "freezegrps": "HBN_PYR",
        "freezedim": "Y Y Y",
        "nstxout-compressed": "1000",
    }

    for field, expected in (
        required_input_values.items()
    ):
        actual = target_input.get(
            field,
            "",
        )

        if (
            " ".join(
                actual.split()
            ).lower()
            !=
            " ".join(
                expected.split()
            ).lower()
        ):
            preparation_failures.append(
                f"{field}: expected "
                f"{expected!r}, found {actual!r}"
            )

    if "define" in target_input:
        preparation_failures.append(
            "A topology preprocessor define "
            "remains active"
        )

    if preparation_failures:
        raise RuntimeError(
            "Matched-control MDP contract failed:\n"
            + "\n".join(
                preparation_failures
            )
        )

    grompp_return_code = run_grompp(
        gmx
    )

    output_files_present = all(
        path.exists()
        and path.stat().st_size > 0
        for path in (
            TARGET_TPR,
            TARGET_MDOUT,
            TARGET_PROCESSED_TOP,
        )
    )

    position_restraint_entries = (
        count_position_restraints(
            TARGET_PROCESSED_TOP
        )
        if TARGET_PROCESSED_TOP.exists()
        else -1
    )

    processed_mdp = (
        parse_mdp(
            TARGET_MDOUT
        )
        if TARGET_MDOUT.exists()
        else {}
    )

    source_dump = dump_tpr(
        gmx,
        SOURCE_TPR,
        SOURCE_DUMP_STDOUT,
        SOURCE_DUMP_STDERR,
    )

    control_dump = dump_tpr(
        gmx,
        TARGET_TPR,
        CONTROL_DUMP_STDOUT,
        CONTROL_DUMP_STDERR,
    )

    (
        source_positions,
        source_velocities,
        _,
    ) = audit.parse_tpr_vectors(
        source_dump
    )

    (
        control_positions,
        control_velocities,
        _,
    ) = audit.parse_tpr_vectors(
        control_dump
    )

    source_water_velocities = (
        source_velocities[
            SOLUTE_COUNT:
        ]
    )

    control_water_velocities = (
        control_velocities[
            SOLUTE_COUNT:
        ]
    )

    source_coordinate_hash = (
        audit.array_hash(
            source_positions
        )
    )

    control_coordinate_hash = (
        audit.array_hash(
            control_positions
        )
    )

    source_water_hash = (
        audit.array_hash(
            source_water_velocities
        )
    )

    control_water_hash = (
        audit.array_hash(
            control_water_velocities
        )
    )

    coordinate_exact_match = bool(
        np.array_equal(
            source_positions,
            control_positions,
        )
    )

    water_velocity_exact_match = bool(
        np.array_equal(
            source_water_velocities,
            control_water_velocities,
        )
    )

    coordinate_max_abs_difference = float(
        np.max(
            np.abs(
                source_positions
                - control_positions
            )
        )
    )

    water_velocity_max_abs_difference = float(
        np.max(
            np.abs(
                source_water_velocities
                - control_water_velocities
            )
        )
    )

    source_water_nonzero_fraction = (
        nonzero_fraction(
            source_water_velocities
        )
    )

    control_water_nonzero_fraction = (
        nonzero_fraction(
            control_water_velocities
        )
    )

    failures: list[str] = []

    if grompp_return_code != 0:
        failures.append(
            "grompp return code is not zero"
        )

    if not output_files_present:
        failures.append(
            "one or more static outputs are missing"
        )

    if position_restraint_entries != 0:
        failures.append(
            "processed topology contains "
            f"{position_restraint_entries} "
            "position-restraint entries"
        )

    processed_checks = {
        "processed freezegrps is HBN_PYR": (
            processed_mdp.get(
                "freezegrps",
                "",
            ).strip().lower()
            == "hbn_pyr"
        ),
        "processed freezedim is Y Y Y": (
            " ".join(
                processed_mdp.get(
                    "freezedim",
                    "",
                ).split()
            ).lower()
            == "y y y"
        ),
        "processed nsteps is 288000": (
            processed_mdp.get(
                "nsteps",
                "",
            ).strip()
            == "288000"
        ),
        "processed dt is 0.0005 ps": (
            abs(
                float(
                    processed_mdp.get(
                        "dt",
                        "nan",
                    )
                )
                - DT_PS
            )
            < 1.0e-12
        ),
        "processed gen-vel is yes": (
            processed_mdp.get(
                "gen-vel",
                "",
            ).strip().lower()
            == "yes"
        ),
        "processed gen-seed is 20260706": (
            processed_mdp.get(
                "gen-seed",
                "",
            ).strip()
            == "20260706"
        ),
    }

    failures.extend(
        label
        for label, passed
        in processed_checks.items()
        if not passed
    )

    if not coordinate_exact_match:
        failures.append(
            "control and Stage02 initial "
            "coordinates are not identical"
        )

    if (
        source_coordinate_hash
        != control_coordinate_hash
    ):
        failures.append(
            "initial coordinate hashes differ"
        )

    if (
        source_water_hash
        != EXPECTED_SOURCE_WATER_HASH
    ):
        failures.append(
            "recomputed Stage02 water-velocity "
            "hash differs from the audited hash"
        )

    if not water_velocity_exact_match:
        failures.append(
            "control and Stage02 water "
            "velocities are not identical"
        )

    if source_water_hash != control_water_hash:
        failures.append(
            "water-velocity hashes differ"
        )

    if (
        abs(
            source_water_nonzero_fraction
            - 0.75
        )
        > 1.0e-12
    ):
        failures.append(
            "Stage02 water nonzero-velocity "
            "fraction is not 0.75"
        )

    if (
        abs(
            control_water_nonzero_fraction
            - 0.75
        )
        > 1.0e-12
    ):
        failures.append(
            "control water nonzero-velocity "
            "fraction is not 0.75"
        )

    execution_authorized = (
        len(failures) == 0
    )

    decision = (
        "PASS"
        if execution_authorized
        else "BLOCKED"
    )

    row = {
        "control_stage": TARGET_STAGE,
        "duration_ps": DURATION_PS,
        "dt_ps": DT_PS,
        "nsteps": NSTEPS,
        "coordinate_interval_ps": (
            DT_PS
            * XTC_INTERVAL_STEPS
        ),
        "expected_xtc_frames": (
            EXPECTED_XTC_FRAMES
        ),
        "initial_coordinate_source": (
            relative(
                START_GRO
            )
        ),
        "source_stage02_tpr": (
            relative(
                SOURCE_TPR
            )
        ),
        "control_tpr": relative(
            TARGET_TPR
        ),
        "grompp_return_code": (
            grompp_return_code
        ),
        "position_restraint_entries": (
            position_restraint_entries
        ),
        "freezegrps": (
            processed_mdp.get(
                "freezegrps",
                "",
            )
        ),
        "freezedim": (
            processed_mdp.get(
                "freezedim",
                "",
            )
        ),
        "source_coordinate_sha256": (
            source_coordinate_hash
        ),
        "control_coordinate_sha256": (
            control_coordinate_hash
        ),
        "coordinate_exact_match": (
            coordinate_exact_match
        ),
        "coordinate_max_abs_difference_nm": (
            coordinate_max_abs_difference
        ),
        "expected_stage02_water_velocity_sha256": (
            EXPECTED_SOURCE_WATER_HASH
        ),
        "source_water_velocity_sha256": (
            source_water_hash
        ),
        "control_water_velocity_sha256": (
            control_water_hash
        ),
        "water_velocity_exact_match": (
            water_velocity_exact_match
        ),
        "water_velocity_max_abs_difference_nm_per_ps": (
            water_velocity_max_abs_difference
        ),
        "source_water_nonzero_fraction": (
            source_water_nonzero_fraction
        ),
        "control_water_nonzero_fraction": (
            control_water_nonzero_fraction
        ),
        "static_decision": decision,
        "matched_control_execution_authorized": (
            execution_authorized
        ),
        "failure_reasons": (
            " | ".join(
                failures
            )
        ),
    }

    write_csv(
        row
    )

    REPORT_MD.write_text(
        f"""# Matched Frozen-Control 144 ps Static Authorization

## Control definition

- Initial coordinates:
  `{relative(START_GRO)}`
- Reference mobile initial state:
  `{relative(SOURCE_TPR)}`
- HBN/PYR treatment: **frozen in X, Y, and Z**
- Position restraints: **{position_restraint_entries}**
- Duration: **{DURATION_PS:.1f} ps**
- Time step: **{DT_PS:.4f} ps**
- Expected trajectory frames: **{EXPECTED_XTC_FRAMES}**

## Initial-state identity

- Coordinate exact match: **{coordinate_exact_match}**
- Coordinate maximum absolute difference:
  **{coordinate_max_abs_difference:.12f} nm**
- Stage02 water-velocity SHA256:
  `{source_water_hash}`
- Control water-velocity SHA256:
  `{control_water_hash}`
- Water-velocity exact match:
  **{water_velocity_exact_match}**
- Water-velocity maximum absolute difference:
  **{water_velocity_max_abs_difference:.12f} nm ps^-1**

## Decision

- Static decision: **{decision}**
- Matched frozen-control execution authorized:
  **{'YES' if execution_authorized else 'NO'}**
- Failure reasons:
  **{'NONE' if not failures else ' | '.join(failures)}**

No MD execution was performed by this preparation workflow.
""",
        encoding="utf-8",
    )

    print(
        "Day022 matched frozen-control "
        "144ps preparation completed."
    )

    print(
        "Control duration / steps / dt: "
        f"{DURATION_PS:.1f} ps / "
        f"{NSTEPS} / {DT_PS:.4f} ps"
    )

    print(
        "Expected compressed trajectory frames: "
        f"{EXPECTED_XTC_FRAMES}"
    )

    print(
        "Static grompp return code: "
        f"{grompp_return_code}"
    )

    print(
        "Active position restraints: "
        f"{position_restraint_entries}/0"
    )

    print(
        "Processed freezegrps / freezedim: "
        f"{processed_mdp.get('freezegrps', 'MISSING')} / "
        f"{processed_mdp.get('freezedim', 'MISSING')}"
    )

    print(
        "Initial coordinate exact match: "
        f"{coordinate_exact_match}"
    )

    print(
        "Initial coordinate maximum difference: "
        f"{coordinate_max_abs_difference:.12f} nm"
    )

    print(
        "Stage02 water velocity SHA256: "
        f"{source_water_hash}"
    )

    print(
        "Control water velocity SHA256: "
        f"{control_water_hash}"
    )

    print(
        "Water velocity exact match: "
        f"{water_velocity_exact_match}"
    )

    print(
        "Water velocity maximum difference: "
        f"{water_velocity_max_abs_difference:.12f} nm/ps"
    )

    print(
        "Stage02/control water nonzero fractions: "
        f"{source_water_nonzero_fraction:.6f}/"
        f"{control_water_nonzero_fraction:.6f}"
    )

    print(
        f"Static decision: {decision}"
    )

    print(
        "Matched frozen-control execution authorized: "
        f"{'YES' if execution_authorized else 'NO'}"
    )

    if failures:
        print(
            "Failure reasons: "
            + " | ".join(
                failures
            )
        )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if failures:
        raise RuntimeError(
            "Matched frozen-control execution "
            "remains blocked."
        )


if __name__ == "__main__":
    main()
