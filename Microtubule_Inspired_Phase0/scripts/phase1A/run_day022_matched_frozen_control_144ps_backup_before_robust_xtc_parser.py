#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
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

CONTROL_ROOT = (
    PROTOCOL
    / "matched_frozen_control_144ps"
)

STATIC_ROOT = (
    CONTROL_ROOT
    / "static_validation"
)

EXECUTION_ROOT = (
    CONTROL_ROOT
    / "execution"
)

STAGE = "matched_frozen_control_144ps"

STATIC_SUMMARY = (
    STATIC_ROOT
    / "matched_frozen_control_144ps_static_summary.csv"
)

STATIC_TPR = (
    STATIC_ROOT
    / f"{STAGE}.tpr"
)

START_GRO = (
    PROTOCOL
    / "execution/01_em_k10000/"
    "01_em_k10000.gro"
)

RUN_TPR = (
    EXECUTION_ROOT
    / f"{STAGE}.tpr"
)

RUN_XTC = (
    EXECUTION_ROOT
    / f"{STAGE}.xtc"
)

RUN_GRO = (
    EXECUTION_ROOT
    / f"{STAGE}.gro"
)

RUN_EDR = (
    EXECUTION_ROOT
    / f"{STAGE}.edr"
)

RUN_LOG = (
    EXECUTION_ROOT
    / f"{STAGE}.log"
)

RUN_CPT = (
    EXECUTION_ROOT
    / f"{STAGE}.cpt"
)

MDRUN_CONSOLE = (
    EXECUTION_ROOT
    / "run_matched_frozen_control_144ps_mdrun_console.log"
)

XTC_CHECK_LOG = (
    EXECUTION_ROOT
    / "run_matched_frozen_control_144ps_xtc_check.log"
)

TEMPERATURE_XVG = (
    EXECUTION_ROOT
    / "matched_frozen_control_temperature.xvg"
)

POTENTIAL_XVG = (
    EXECUTION_ROOT
    / "matched_frozen_control_potential.xvg"
)

PRESSURE_XVG = (
    EXECUTION_ROOT
    / "matched_frozen_control_pressure.xvg"
)

SUMMARY_CSV = (
    EXECUTION_ROOT
    / "matched_frozen_control_144ps_execution_summary.csv"
)

REPORT_MD = (
    EXECUTION_ROOT
    / "MATCHED_FROZEN_CONTROL_144PS_EXECUTION_DAY022.md"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26

SOLUTE_COUNT = (
    HBN_COUNT
    + PYR_COUNT * PYR_ATOMS
)

TOTAL_ATOMS = 68320

EXPECTED_FRAMES = 289
EXPECTED_DURATION_PS = 144.0
EXPECTED_INTERVAL_PS = 0.5

EXPECTED_WATER_VELOCITY_HASH = (
    "d11ac84fa2d4a2a9a594a91ac0a6e0714"
    "dc3caa9da33f0ed12617371feff722d"
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


def truthy(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def read_single_row(
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


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


def validate_static_authorization() -> dict[str, str]:
    row = read_single_row(
        STATIC_SUMMARY
    )

    failures = []

    if (
        row.get(
            "static_decision",
            "",
        ).strip().upper()
        != "PASS"
    ):
        failures.append(
            "static decision is not PASS"
        )

    if not truthy(
        row.get(
            "matched_control_execution_authorized",
            "",
        )
    ):
        failures.append(
            "execution authorization is not TRUE"
        )

    if not truthy(
        row.get(
            "coordinate_exact_match",
            "",
        )
    ):
        failures.append(
            "initial coordinates are not identical"
        )

    if not truthy(
        row.get(
            "water_velocity_exact_match",
            "",
        )
    ):
        failures.append(
            "initial water velocities are not identical"
        )

    source_hash = row.get(
        "source_water_velocity_sha256",
        "",
    ).strip()

    control_hash = row.get(
        "control_water_velocity_sha256",
        "",
    ).strip()

    expected_hash = row.get(
        "expected_stage02_water_velocity_sha256",
        "",
    ).strip()

    if source_hash != EXPECTED_WATER_VELOCITY_HASH:
        failures.append(
            "source water-velocity hash differs "
            "from the audited hash"
        )

    if control_hash != EXPECTED_WATER_VELOCITY_HASH:
        failures.append(
            "control water-velocity hash differs "
            "from the audited hash"
        )

    if expected_hash != EXPECTED_WATER_VELOCITY_HASH:
        failures.append(
            "stored expected water-velocity hash differs"
        )

    if row.get(
        "freezegrps",
        "",
    ).strip().lower() != "hbn_pyr":
        failures.append(
            "freezegrps is not HBN_PYR"
        )

    if (
        " ".join(
            row.get(
                "freezedim",
                "",
            ).split()
        ).lower()
        != "y y y"
    ):
        failures.append(
            "freezedim is not Y Y Y"
        )

    if int(
        float(
            row.get(
                "position_restraint_entries",
                "-1",
            )
        )
    ) != 0:
        failures.append(
            "position-restraint entries are not zero"
        )

    if failures:
        raise RuntimeError(
            "Static authorization validation failed:\n"
            + "\n".join(failures)
        )

    return row


def prepare_execution_directory() -> tuple[str, str]:
    EXECUTION_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    protected_outputs = (
        RUN_XTC,
        RUN_GRO,
        RUN_EDR,
        RUN_LOG,
        RUN_CPT,
        SUMMARY_CSV,
    )

    existing = [
        path
        for path in protected_outputs
        if path.exists()
    ]

    if existing:
        raise RuntimeError(
            "Execution outputs already exist; "
            "refusing to overwrite:\n"
            + "\n".join(
                str(path)
                for path in existing
            )
        )

    shutil.copy2(
        STATIC_TPR,
        RUN_TPR,
    )

    static_hash = sha256_file(
        STATIC_TPR
    )

    run_hash = sha256_file(
        RUN_TPR
    )

    if static_hash != run_hash:
        raise RuntimeError(
            "Copied TPR hash differs from "
            "the authorized static TPR"
        )

    return static_hash, run_hash


def run_mdrun(gmx: str) -> int:
    command = [
        gmx,
        "mdrun",
        "-deffnm",
        STAGE,
    ]

    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    print()
    print(
        "===== MATCHED FROZEN CONTROL MDRUN ====="
    )
    print(
        "Command: "
        + " ".join(command)
    )
    print()

    with MDRUN_CONSOLE.open(
        "w",
        encoding="utf-8",
    ) as log_handle:
        process = subprocess.Popen(
            command,
            cwd=EXECUTION_ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )

        if process.stdout is None:
            raise RuntimeError(
                "Could not capture mdrun output"
            )

        for line in process.stdout:
            print(
                line,
                end="",
                flush=True,
            )

            log_handle.write(line)
            log_handle.flush()

        return_code = process.wait()

    return return_code


def run_xtc_check(
    gmx: str,
) -> tuple[int, float, float]:
    completed = subprocess.run(
        [
            gmx,
            "check",
            "-f",
            str(RUN_XTC),
        ],
        cwd=EXECUTION_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    XTC_CHECK_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "gmx check failed. "
            f"See {XTC_CHECK_LOG}"
        )

    frame_match = re.search(
        r"Coords\s+(\d+)\s+"
        r"([-+0-9.eE]+)",
        completed.stdout,
    )

    last_match = re.search(
        r"Last frame\s+\d+\s+time\s+"
        r"([-+0-9.eE]+)",
        completed.stdout,
    )

    if (
        frame_match is None
        or last_match is None
    ):
        raise RuntimeError(
            "Could not parse XTC frame count "
            "or final time"
        )

    return (
        int(
            frame_match.group(1)
        ),
        float(
            frame_match.group(2)
        ),
        float(
            last_match.group(1)
        ),
    )


def extract_energy(
    gmx: str,
    term: str,
    output: Path,
) -> np.ndarray:
    completed = subprocess.run(
        [
            gmx,
            "energy",
            "-f",
            str(RUN_EDR),
            "-o",
            str(output),
        ],
        cwd=EXECUTION_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        input=f"{term}\n0\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    if (
        completed.returncode != 0
        or not output.exists()
        or output.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Could not extract energy term: {term}"
        )

    values = []

    for raw_line in output.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.strip()

        if (
            not line
            or line.startswith("#")
            or line.startswith("@")
        ):
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        values.append(
            (
                float(fields[0]),
                float(fields[1]),
            )
        )

    if not values:
        raise RuntimeError(
            f"No numeric values in {output}"
        )

    array = np.array(
        values,
        dtype=float,
    )

    if not np.all(
        np.isfinite(array)
    ):
        raise RuntimeError(
            f"Non-finite values in {output}"
        )

    return array


def read_gro_positions(
    path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    natoms = int(
        lines[1].strip()
    )

    if natoms != TOTAL_ATOMS:
        raise RuntimeError(
            f"Expected {TOTAL_ATOMS} atoms in "
            f"{path}; found {natoms}"
        )

    positions = np.empty(
        (natoms, 3),
        dtype=float,
    )

    for atom_index, line in enumerate(
        lines[2 : 2 + natoms]
    ):
        positions[
            atom_index
        ] = (
            float(line[20:28]),
            float(line[28:36]),
            float(line[36:44]),
        )

    box_values = [
        float(value)
        for value in lines[
            2 + natoms
        ].split()
    ]

    if len(box_values) < 3:
        raise RuntimeError(
            f"Invalid box in {path}"
        )

    box = np.array(
        box_values[:3],
        dtype=float,
    )

    return positions, box


def frozen_solute_displacement() -> tuple[
    float,
    float,
    float,
]:
    initial_positions, initial_box = (
        read_gro_positions(
            START_GRO
        )
    )

    final_positions, final_box = (
        read_gro_positions(
            RUN_GRO
        )
    )

    if not np.allclose(
        initial_box,
        final_box,
        atol=1.0e-6,
        rtol=0.0,
    ):
        raise RuntimeError(
            "NVT control box changed unexpectedly"
        )

    displacement = (
        final_positions
        - initial_positions
    )

    displacement -= (
        final_box
        * np.round(
            displacement / final_box
        )
    )

    solute_distance = np.linalg.norm(
        displacement[
            :SOLUTE_COUNT
        ],
        axis=1,
    )

    water_distance = np.linalg.norm(
        displacement[
            SOLUTE_COUNT:
        ],
        axis=1,
    )

    solute_rms = float(
        np.sqrt(
            np.mean(
                solute_distance ** 2
            )
        )
    )

    solute_max = float(
        solute_distance.max()
    )

    water_rms = float(
        np.sqrt(
            np.mean(
                water_distance ** 2
            )
        )
    )

    return (
        solute_rms,
        solute_max,
        water_rms,
    )


def instability_count() -> tuple[int, list[str]]:
    patterns = (
        r"\bFatal error\b",
        r"\bSegmentation fault\b",
        r"\bnon-finite\b",
        r"\bNaN\b",
        r"\bLINCS WARNING\b",
        r"\bSETTLE.*error\b",
        r"\bconstraint failure\b",
    )

    texts = []

    for path in (
        MDRUN_CONSOLE,
        RUN_LOG,
    ):
        if path.exists():
            texts.append(
                path.read_text(
                    encoding="utf-8",
                    errors="replace",
                )
            )

    combined = "\n".join(texts)

    matched = [
        pattern
        for pattern in patterns
        if re.search(
            pattern,
            combined,
            flags=re.IGNORECASE,
        )
    ]

    return (
        len(matched),
        matched,
    )


def write_summary(
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
        STATIC_SUMMARY,
        STATIC_TPR,
        START_GRO,
    ):
        require_file(required)

    authorization = (
        validate_static_authorization()
    )

    static_tpr_hash, run_tpr_hash = (
        prepare_execution_directory()
    )

    gmx = locate_gmx()

    return_code = run_mdrun(
        gmx
    )

    required_outputs = (
        RUN_XTC,
        RUN_GRO,
        RUN_EDR,
        RUN_LOG,
        RUN_CPT,
    )

    missing_outputs = [
        path
        for path in required_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    failures = []

    if return_code != 0:
        failures.append(
            f"mdrun return code is {return_code}"
        )

    if missing_outputs:
        failures.append(
            "missing outputs: "
            + ", ".join(
                path.name
                for path in missing_outputs
            )
        )

    if failures:
        raise RuntimeError(
            "Matched frozen-control execution failed:\n"
            + "\n".join(failures)
        )

    frames, timestep_ps, final_time_ps = (
        run_xtc_check(
            gmx
        )
    )

    temperature = extract_energy(
        gmx,
        "Temperature",
        TEMPERATURE_XVG,
    )

    potential = extract_energy(
        gmx,
        "Potential",
        POTENTIAL_XVG,
    )

    pressure = extract_energy(
        gmx,
        "Pressure",
        PRESSURE_XVG,
    )

    (
        solute_rms_nm,
        solute_max_nm,
        water_rms_nm,
    ) = frozen_solute_displacement()

    signature_count, signatures = (
        instability_count()
    )

    temperature_values = temperature[:, 1]
    potential_values = potential[:, 1]
    pressure_values = pressure[:, 1]

    validation_failures = []

    if frames != EXPECTED_FRAMES:
        validation_failures.append(
            f"expected {EXPECTED_FRAMES} frames; "
            f"found {frames}"
        )

    if not math.isclose(
        timestep_ps,
        EXPECTED_INTERVAL_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        validation_failures.append(
            f"expected {EXPECTED_INTERVAL_PS} ps "
            f"frame interval; found {timestep_ps}"
        )

    if not math.isclose(
        final_time_ps,
        EXPECTED_DURATION_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        validation_failures.append(
            f"expected {EXPECTED_DURATION_PS} ps; "
            f"found {final_time_ps}"
        )

    if signature_count != 0:
        validation_failures.append(
            "instability signatures detected: "
            + " | ".join(signatures)
        )

    if (
        float(
            temperature_values.mean()
        )
        < 290.0
        or float(
            temperature_values.mean()
        )
        > 310.0
    ):
        validation_failures.append(
            "mean temperature outside 290-310 K"
        )

    if solute_max_nm > 0.001:
        validation_failures.append(
            "frozen HBN/PYR maximum displacement "
            "exceeds 0.001 nm"
        )

    if water_rms_nm <= 0.01:
        validation_failures.append(
            "water atoms did not exhibit expected motion"
        )

    decision = (
        "PASS"
        if not validation_failures
        else "REVIEW"
    )

    row = {
        "stage": STAGE,
        "duration_ps": final_time_ps,
        "frames": frames,
        "coordinate_interval_ps": (
            timestep_ps
        ),
        "mdrun_return_code": (
            return_code
        ),
        "temperature_mean_K": float(
            temperature_values.mean()
        ),
        "temperature_std_K": float(
            temperature_values.std()
        ),
        "temperature_min_K": float(
            temperature_values.min()
        ),
        "temperature_max_K": float(
            temperature_values.max()
        ),
        "temperature_last_K": float(
            temperature_values[-1]
        ),
        "potential_first_kJ_mol": float(
            potential_values[0]
        ),
        "potential_last_kJ_mol": float(
            potential_values[-1]
        ),
        "pressure_mean_bar": float(
            pressure_values.mean()
        ),
        "pressure_min_bar": float(
            pressure_values.min()
        ),
        "pressure_max_bar": float(
            pressure_values.max()
        ),
        "frozen_solute_rms_displacement_nm": (
            solute_rms_nm
        ),
        "frozen_solute_max_displacement_nm": (
            solute_max_nm
        ),
        "water_atom_rms_displacement_nm": (
            water_rms_nm
        ),
        "instability_signature_count": (
            signature_count
        ),
        "static_tpr_sha256": (
            static_tpr_hash
        ),
        "execution_tpr_sha256": (
            run_tpr_hash
        ),
        "initial_water_velocity_sha256": (
            authorization[
                "control_water_velocity_sha256"
            ]
        ),
        "execution_decision": decision,
        "matched_comparison_authorized": (
            decision == "PASS"
        ),
        "electronic_recalculation_authorized": (
            False
        ),
        "validation_reasons": (
            " | ".join(
                validation_failures
            )
        ),
    }

    write_summary(row)

    REPORT_MD.write_text(
        f"""# Matched Frozen-Control 144 ps Execution

## Execution

- Duration: **{final_time_ps:.1f} ps**
- Frames: **{frames}**
- Coordinate interval: **{timestep_ps:.3f} ps**
- `mdrun` return code: **{return_code}**
- Execution decision: **{decision}**

## Thermodynamics

- Temperature mean/std/min/max:
  {row['temperature_mean_K']:.4f}/
  {row['temperature_std_K']:.4f}/
  {row['temperature_min_K']:.4f}/
  {row['temperature_max_K']:.4f} K
- Potential energy first/last:
  {row['potential_first_kJ_mol']:.4f}/
  {row['potential_last_kJ_mol']:.4f} kJ mol^-1
- Pressure mean/min/max:
  {row['pressure_mean_bar']:.4f}/
  {row['pressure_min_bar']:.4f}/
  {row['pressure_max_bar']:.4f} bar

## Frozen-solute verification

- HBN/PYR RMS displacement:
  {solute_rms_nm:.12f} nm
- HBN/PYR maximum displacement:
  {solute_max_nm:.12f} nm
- Water-atom RMS displacement:
  {water_rms_nm:.6f} nm
- Instability signatures:
  {signature_count}

## Scope

If execution passes, the authorized comparison is:

- matched frozen control: **44-144 ps**
- mobile Stage08: **0-100 ps**

Electronic recalculation remains unauthorized until the matched
water comparison and representative-snapshot review are complete.
""",
        encoding="utf-8",
    )

    print()
    print(
        "===== MATCHED FROZEN CONTROL VALIDATION ====="
    )

    print(
        "Trajectory frames / duration / interval: "
        f"{frames} / {final_time_ps:.3f} ps / "
        f"{timestep_ps:.3f} ps"
    )

    print(
        "Temperature mean/std/min/max/last: "
        f"{row['temperature_mean_K']:.4f}/"
        f"{row['temperature_std_K']:.4f}/"
        f"{row['temperature_min_K']:.4f}/"
        f"{row['temperature_max_K']:.4f}/"
        f"{row['temperature_last_K']:.4f} K"
    )

    print(
        "Potential energy first/last: "
        f"{row['potential_first_kJ_mol']:.4f}/"
        f"{row['potential_last_kJ_mol']:.4f} kJ/mol"
    )

    print(
        "Pressure mean/min/max: "
        f"{row['pressure_mean_bar']:.4f}/"
        f"{row['pressure_min_bar']:.4f}/"
        f"{row['pressure_max_bar']:.4f} bar"
    )

    print(
        "Frozen HBN/PYR RMS/max displacement: "
        f"{solute_rms_nm:.12f}/"
        f"{solute_max_nm:.12f} nm"
    )

    print(
        "Water atom RMS displacement: "
        f"{water_rms_nm:.6f} nm"
    )

    print(
        "Instability signatures: "
        f"{signature_count}"
    )

    print(
        f"Execution decision: {decision}"
    )

    print(
        "Matched 44-144 ps comparison authorized: "
        f"{'YES' if decision == 'PASS' else 'NO'}"
    )

    print(
        "Electronic recalculation authorized: NO"
    )

    if validation_failures:
        print(
            "Validation reasons: "
            + " | ".join(
                validation_failures
            )
        )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if validation_failures:
        raise RuntimeError(
            "Matched frozen-control execution "
            "requires review."
        )


if __name__ == "__main__":
    main()
