#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

SOURCE_SCRIPT = (
    ROOT
    / "scripts/phase1A/"
    "run_day022_matched_frozen_control_144ps.py"
)


def load_runner_module():
    if (
        not SOURCE_SCRIPT.exists()
        or SOURCE_SCRIPT.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing runner script: {SOURCE_SCRIPT}"
        )

    specification = (
        importlib.util.spec_from_file_location(
            "matched_frozen_control_runner",
            SOURCE_SCRIPT,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            f"Could not load runner module: {SOURCE_SCRIPT}"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(
        module
    )

    return module


def main() -> None:
    module = load_runner_module()

    required_files = (
        module.STATIC_SUMMARY,
        module.STATIC_TPR,
        module.RUN_TPR,
        module.RUN_XTC,
        module.RUN_GRO,
        module.RUN_EDR,
        module.RUN_LOG,
        module.RUN_CPT,
        module.MDRUN_CONSOLE,
        module.START_GRO,
    )

    for path in required_files:
        module.require_file(
            path
        )

    authorization = (
        module.validate_static_authorization()
    )

    static_tpr_hash = module.sha256_file(
        module.STATIC_TPR
    )

    execution_tpr_hash = module.sha256_file(
        module.RUN_TPR
    )

    if static_tpr_hash != execution_tpr_hash:
        raise RuntimeError(
            "The execution TPR differs from the "
            "statically authorized TPR."
        )

    console_text = (
        module.MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    completion_checks = {
        "expected step count": (
            "288000 steps"
            in console_text
        ),
        "final coordinates written": (
            "Writing final coordinates"
            in console_text
        ),
        "performance reported": (
            "Performance:"
            in console_text
        ),
    }

    missing_completion_checks = [
        label
        for label, passed
        in completion_checks.items()
        if not passed
    ]

    if missing_completion_checks:
        raise RuntimeError(
            "The completed MD execution could not be "
            "verified from its console log:\n"
            + "\n".join(
                missing_completion_checks
            )
        )

    gmx = module.locate_gmx()

    (
        frames,
        coordinate_interval_ps,
        final_time_ps,
    ) = module.run_xtc_check(
        gmx
    )

    temperature = module.extract_energy(
        gmx,
        "Temperature",
        module.TEMPERATURE_XVG,
    )

    potential = module.extract_energy(
        gmx,
        "Potential",
        module.POTENTIAL_XVG,
    )

    pressure = module.extract_energy(
        gmx,
        "Pressure",
        module.PRESSURE_XVG,
    )

    (
        frozen_solute_rms_nm,
        frozen_solute_max_nm,
        water_rms_nm,
    ) = module.frozen_solute_displacement()

    (
        instability_signature_count,
        instability_signatures,
    ) = module.instability_count()

    temperature_values = np.asarray(
        temperature[:, 1],
        dtype=float,
    )

    potential_values = np.asarray(
        potential[:, 1],
        dtype=float,
    )

    pressure_values = np.asarray(
        pressure[:, 1],
        dtype=float,
    )

    numeric_arrays = (
        temperature_values,
        potential_values,
        pressure_values,
    )

    if not all(
        np.all(
            np.isfinite(array)
        )
        for array in numeric_arrays
    ):
        raise RuntimeError(
            "Non-finite thermodynamic values were found."
        )

    failures: list[str] = []

    if frames != module.EXPECTED_FRAMES:
        failures.append(
            f"expected {module.EXPECTED_FRAMES} frames; "
            f"found {frames}"
        )

    if not math.isclose(
        coordinate_interval_ps,
        module.EXPECTED_INTERVAL_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        failures.append(
            "expected coordinate interval "
            f"{module.EXPECTED_INTERVAL_PS:.6f} ps; "
            f"found {coordinate_interval_ps:.6f} ps"
        )

    if not math.isclose(
        final_time_ps,
        module.EXPECTED_DURATION_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        failures.append(
            "expected duration "
            f"{module.EXPECTED_DURATION_PS:.6f} ps; "
            f"found {final_time_ps:.6f} ps"
        )

    temperature_mean = float(
        temperature_values.mean()
    )

    if not (
        290.0
        <= temperature_mean
        <= 310.0
    ):
        failures.append(
            "mean temperature is outside 290-310 K"
        )

    if instability_signature_count != 0:
        failures.append(
            "instability signatures detected: "
            + " | ".join(
                instability_signatures
            )
        )

    if frozen_solute_max_nm > 0.001:
        failures.append(
            "frozen HBN/PYR maximum displacement "
            "exceeds 0.001 nm"
        )

    if water_rms_nm <= 0.01:
        failures.append(
            "water atoms did not exhibit the expected motion"
        )

    decision = (
        "PASS"
        if not failures
        else "REVIEW"
    )

    row = {
        "stage": module.STAGE,
        "validation_mode": (
            "POSTPROCESS_RECOVERY_NO_MDRUN"
        ),
        "existing_md_reused": True,
        "new_mdrun_executed": False,
        "mdrun_completion_verified": True,
        "duration_ps": final_time_ps,
        "frames": frames,
        "coordinate_interval_ps": (
            coordinate_interval_ps
        ),
        "temperature_mean_K": (
            temperature_mean
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
            frozen_solute_rms_nm
        ),
        "frozen_solute_max_displacement_nm": (
            frozen_solute_max_nm
        ),
        "water_atom_rms_displacement_nm": (
            water_rms_nm
        ),
        "instability_signature_count": (
            instability_signature_count
        ),
        "static_tpr_sha256": (
            static_tpr_hash
        ),
        "execution_tpr_sha256": (
            execution_tpr_hash
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
                failures
            )
        ),
    }

    module.write_summary(
        row
    )

    module.REPORT_MD.write_text(
        f"""# Matched Frozen-Control 144 ps Recovered Validation

## Recovery scope

- Existing MD outputs reused: **YES**
- New `mdrun` executed: **NO**
- Original execution completion verified: **YES**
- Authorized and execution TPR hashes identical: **YES**

## Trajectory

- Frames: **{frames}**
- Duration: **{final_time_ps:.3f} ps**
- Coordinate interval: **{coordinate_interval_ps:.3f} ps**
- Execution decision: **{decision}**

## Thermodynamics

- Temperature mean/std/min/max/last:
  {row['temperature_mean_K']:.4f}/
  {row['temperature_std_K']:.4f}/
  {row['temperature_min_K']:.4f}/
  {row['temperature_max_K']:.4f}/
  {row['temperature_last_K']:.4f} K
- Potential energy first/last:
  {row['potential_first_kJ_mol']:.4f}/
  {row['potential_last_kJ_mol']:.4f} kJ mol^-1
- Pressure mean/min/max:
  {row['pressure_mean_bar']:.4f}/
  {row['pressure_min_bar']:.4f}/
  {row['pressure_max_bar']:.4f} bar

## Frozen-solute verification

- HBN/PYR RMS displacement:
  {frozen_solute_rms_nm:.12f} nm
- HBN/PYR maximum displacement:
  {frozen_solute_max_nm:.12f} nm
- Water-atom RMS displacement:
  {water_rms_nm:.6f} nm
- Instability signatures:
  {instability_signature_count}

## Authorization

- Matched comparison, frozen 44-144 ps versus
  mobile Stage08 0-100 ps:
  **{'YES' if decision == 'PASS' else 'NO'}**
- Electronic recalculation:
  **NO**
- Validation reasons:
  **{'NONE' if not failures else ' | '.join(failures)}**
""",
        encoding="utf-8",
    )

    print()
    print(
        "===== MATCHED FROZEN CONTROL "
        "RECOVERED VALIDATION ====="
    )

    print(
        "Existing MD reused / new mdrun executed: "
        "YES / NO"
    )

    print(
        "Trajectory frames / duration / interval: "
        f"{frames} / "
        f"{final_time_ps:.3f} ps / "
        f"{coordinate_interval_ps:.3f} ps"
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
        f"{row['potential_last_kJ_mol']:.4f} "
        "kJ/mol"
    )

    print(
        "Pressure mean/min/max: "
        f"{row['pressure_mean_bar']:.4f}/"
        f"{row['pressure_min_bar']:.4f}/"
        f"{row['pressure_max_bar']:.4f} bar"
    )

    print(
        "Frozen HBN/PYR RMS/max displacement: "
        f"{frozen_solute_rms_nm:.12f}/"
        f"{frozen_solute_max_nm:.12f} nm"
    )

    print(
        "Water atom RMS displacement: "
        f"{water_rms_nm:.6f} nm"
    )

    print(
        "Instability signatures: "
        f"{instability_signature_count}"
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

    if failures:
        print(
            "Validation reasons: "
            + " | ".join(
                failures
            )
        )

    print(
        "Wrote: "
        + module.relative(
            module.REPORT_MD
        )
    )

    if failures:
        raise RuntimeError(
            "Recovered validation requires review."
        )


if __name__ == "__main__":
    main()
