#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PROTOCOL_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

INPUT_ROOT = (
    PROTOCOL_ROOT
    / "protocol_inputs"
)

TOPOLOGY_ROOT = (
    INPUT_ROOT
    / "topology"
)

MDP = (
    INPUT_ROOT
    / "mdp/00_em_k100000.mdp"
)

TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    INPUT_ROOT
    / "mobile_release_index.ndx"
)

START_GRO = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute/"
    "nvt_100ps_frozenSolute.gro"
)

STATIC_CSV = (
    PROTOCOL_ROOT
    / "static_validation/"
    "mobile_release_static_validation.csv"
)

INPUT_MANIFEST = (
    INPUT_ROOT
    / "mobile_release_protocol_manifest.csv"
)

RUN_ROOT = (
    PROTOCOL_ROOT
    / "execution/00_em_k100000"
)

PREFIX_NAME = "00_em_k100000"

TPR = RUN_ROOT / f"{PREFIX_NAME}.tpr"
PROCESSED_TOP = RUN_ROOT / f"{PREFIX_NAME}_processed.top"
MDOUT = RUN_ROOT / f"{PREFIX_NAME}_mdout.mdp"
GROMPP_LOG = RUN_ROOT / f"{PREFIX_NAME}_grompp.log"
MDRUN_CONSOLE_LOG = RUN_ROOT / f"{PREFIX_NAME}_mdrun_console.log"
MDRUN_LOG = RUN_ROOT / f"{PREFIX_NAME}.log"
FINAL_GRO = RUN_ROOT / f"{PREFIX_NAME}.gro"
EDR = RUN_ROOT / f"{PREFIX_NAME}.edr"

COMMANDS_TXT = RUN_ROOT / "commands.txt"
SUMMARY_CSV = RUN_ROOT / "stage00_em_summary.csv"
REPORT_MD = RUN_ROOT / "STAGE00_EM_K100000_DAY021.md"
OUTPUT_MANIFEST = RUN_ROOT / "stage00_output_manifest.csv"

EXPECTED_ATOMS = 68320
EXPECTED_EMTOL = 1000.0
MAX_ALLOWED_RESTRAINED_DISPLACEMENT_NM = 0.05

HBN_RANGE = (1, 1680)
PYR_RANGE = (1681, 1784)

FLOAT_PATTERN = (
    r"[+-]?"
    r"(?:\d+(?:\.\d*)?|\.\d+)"
    r"(?:[Ee][+-]?\d+)?"
)


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(PROJECT_ROOT)
        )
    except ValueError:
        return str(path.resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def find_gmx() -> Path:
    preferred = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if preferred.exists():
        return preferred

    discovered = shutil.which("gmx")

    if discovered is None:
        raise RuntimeError(
            "Could not locate the GROMACS gmx executable"
        )

    return Path(discovered)


def require_inputs() -> None:
    required = (
        MDP,
        TOP,
        NDX,
        START_GRO,
        STATIC_CSV,
        INPUT_MANIFEST,
    )

    missing = [
        path
        for path in required
        if not path.exists()
        or path.stat().st_size == 0
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty required inputs:\n"
            + "\n".join(str(path) for path in missing)
        )


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def verify_static_validation() -> None:
    with STATIC_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            "Static validation CSV is empty"
        )

    failed = []

    for row in rows:
        if not parse_bool(
            row.get("grompp_pass", "")
        ):
            failed.append(
                f"{row.get('stage', '?')}: grompp"
            )

        if not parse_bool(
            row.get(
                "restraint_entry_validation",
                "",
            )
        ):
            failed.append(
                f"{row.get('stage', '?')}: restraints"
            )

    if failed:
        raise RuntimeError(
            "Protocol static validation is not complete:\n"
            + "\n".join(failed)
        )


def verify_input_manifest() -> int:
    with INPUT_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            "Protocol input manifest is empty"
        )

    mismatches = []

    for row in rows:
        path = PROJECT_ROOT / row["path"]
        expected_hash = row["sha256"]

        if not path.exists():
            mismatches.append(
                f"missing: {row['path']}"
            )
            continue

        current_hash = sha256(path)

        if current_hash != expected_hash:
            mismatches.append(
                f"hash mismatch: {row['path']}"
            )

    if mismatches:
        raise RuntimeError(
            "Protocol input-manifest verification failed:\n"
            + "\n".join(mismatches)
        )

    return len(rows)


def ensure_clean_run_target() -> None:
    RUN_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    existing_outputs = [
        path
        for path in (
            TPR,
            MDRUN_LOG,
            FINAL_GRO,
            EDR,
        )
        if path.exists()
    ]

    if existing_outputs:
        raise RuntimeError(
            "Stage00 execution outputs already exist. "
            "No files were overwritten:\n"
            + "\n".join(
                str(path)
                for path in existing_outputs
            )
        )


def run_streaming(
    command: list[str],
    cwd: Path,
    log_path: Path,
) -> int:
    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    with log_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )

        assert process.stdout is not None

        for line in process.stdout:
            handle.write(line)
            handle.flush()
            print(
                line,
                end="",
                flush=True,
            )

        return process.wait()


def read_gro(
    path: Path,
) -> tuple[
    list[tuple[str, str, int, float, float, float]],
    list[float],
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    natoms = int(lines[1].strip())

    atom_lines = lines[2 : 2 + natoms]

    if len(atom_lines) != natoms:
        raise RuntimeError(
            f"GRO atom-count mismatch in {path}"
        )

    atoms = []

    for line in atom_lines:
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom record in {path}"
            )

        residue_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        atom_number = int(line[15:20])

        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])

        atoms.append(
            (
                residue_name,
                atom_name,
                atom_number,
                x,
                y,
                z,
            )
        )

    box = [
        float(value)
        for value in lines[2 + natoms].split()
    ]

    return atoms, box


def group_displacement(
    start_atoms,
    final_atoms,
    first_atom: int,
    last_atom: int,
) -> dict[str, float | int]:
    squared_sum = 0.0
    max_displacement = -1.0
    max_atom = -1

    for atom_index in range(
        first_atom - 1,
        last_atom,
    ):
        start = start_atoms[atom_index]
        final = final_atoms[atom_index]

        if start[:3] != final[:3]:
            raise RuntimeError(
                "Atom identity/order mismatch at "
                f"atom {atom_index + 1}"
            )

        dx = final[3] - start[3]
        dy = final[4] - start[4]
        dz = final[5] - start[5]

        displacement = math.sqrt(
            dx * dx
            + dy * dy
            + dz * dz
        )

        squared_sum += displacement * displacement

        if displacement > max_displacement:
            max_displacement = displacement
            max_atom = atom_index + 1

    atom_count = (
        last_atom
        - first_atom
        + 1
    )

    return {
        "atom_count": atom_count,
        "rms_displacement_nm": math.sqrt(
            squared_sum / atom_count
        ),
        "max_displacement_nm": (
            max_displacement
        ),
        "max_displacement_atom": (
            max_atom
        ),
    }


def last_float_match(
    text: str,
    pattern: str,
) -> float | None:
    matches = re.findall(
        pattern,
        text,
        flags=re.IGNORECASE,
    )

    if not matches:
        return None

    value = matches[-1]

    if isinstance(value, tuple):
        value = value[0]

    return float(value)


def parse_mdrun_log(
    path: Path,
) -> dict[str, object]:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    potential_energy = last_float_match(
        text,
        rf"Potential Energy\s*=\s*({FLOAT_PATTERN})",
    )

    maximum_force = last_float_match(
        text,
        rf"Maximum force\s*=\s*({FLOAT_PATTERN})",
    )

    norm_force = last_float_match(
        text,
        rf"Norm of force\s*=\s*({FLOAT_PATTERN})",
    )

    atom_matches = re.findall(
        r"Maximum force\s*=\s*"
        + FLOAT_PATTERN
        + r"\s+on atom\s+(\d+)",
        text,
        flags=re.IGNORECASE,
    )

    maximum_force_atom = (
        int(atom_matches[-1])
        if atom_matches
        else None
    )

    convergence_lines = [
        line.strip()
        for line in text.splitlines()
        if (
            "Steepest Descents" in line
            and (
                "converged" in line.lower()
                or "machine precision" in line.lower()
            )
        )
    ]

    convergence_message = (
        convergence_lines[-1]
        if convergence_lines
        else ""
    )

    converged = (
        bool(convergence_message)
        and "converged to fmax"
        in convergence_message.lower()
        and "did not converge"
        not in convergence_message.lower()
    )

    step_match = re.search(
        r"in\s+(\d+)\s+steps",
        convergence_message,
        flags=re.IGNORECASE,
    )

    steps = (
        int(step_match.group(1))
        if step_match
        else None
    )

    return {
        "potential_energy_kj_mol": (
            potential_energy
        ),
        "maximum_force_kj_mol_nm": (
            maximum_force
        ),
        "maximum_force_atom": (
            maximum_force_atom
        ),
        "norm_force_kj_mol_nm": (
            norm_force
        ),
        "convergence_message": (
            convergence_message
        ),
        "converged": converged,
        "steps": steps,
    }


def finite_coordinates(atoms) -> bool:
    return all(
        math.isfinite(value)
        for atom in atoms
        for value in atom[3:6]
    )


def write_csv(
    path: Path,
    row: dict[str, object],
) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(row.keys()),
        )

        writer.writeheader()
        writer.writerow(row)


def write_output_manifest() -> None:
    paths = sorted(
        path
        for path in RUN_ROOT.iterdir()
        if path.is_file()
        and path != OUTPUT_MANIFEST
    )

    rows = [
        {
            "path": relative(path),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in paths
    ]

    with OUTPUT_MANIFEST.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "path",
                "size_bytes",
                "sha256",
            ],
        )

        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    require_inputs()
    verify_static_validation()

    verified_manifest_files = (
        verify_input_manifest()
    )

    ensure_clean_run_target()

    gmx = find_gmx()

    grompp_command = [
        str(gmx),
        "grompp",
        "-f",
        str(MDP),
        "-c",
        str(START_GRO),
        "-r",
        str(START_GRO),
        "-p",
        str(TOP),
        "-n",
        str(NDX),
        "-o",
        str(TPR),
        "-po",
        str(MDOUT),
        "-pp",
        str(PROCESSED_TOP),
        "-maxwarn",
        "0",
    ]

    mdrun_command = [
        str(gmx),
        "mdrun",
        "-deffnm",
        PREFIX_NAME,
    ]

    COMMANDS_TXT.write_text(
        " ".join(grompp_command)
        + "\n\n"
        + " ".join(mdrun_command)
        + "\n",
        encoding="utf-8",
    )

    log(
        "Day021 Stage00 strongly restrained "
        "minimization"
    )

    log(
        f"Verified protocol-manifest files: "
        f"{verified_manifest_files}"
    )

    log(
        "Static protocol validation: PASS"
    )

    log(
        "Accepted source files modified: NO"
    )

    log()
    log(
        "===== GROMPP STAGE00 ====="
    )

    grompp_return_code = run_streaming(
        grompp_command,
        TOPOLOGY_ROOT,
        GROMPP_LOG,
    )

    grompp_pass = (
        grompp_return_code == 0
        and TPR.exists()
        and TPR.stat().st_size > 0
    )

    if not grompp_pass:
        raise RuntimeError(
            "Stage00 grompp failed. "
            f"See {GROMPP_LOG}"
        )

    log()
    log(
        "===== MDRUN STAGE00 ====="
    )

    mdrun_return_code = run_streaming(
        mdrun_command,
        RUN_ROOT,
        MDRUN_CONSOLE_LOG,
    )

    required_outputs = (
        MDRUN_LOG,
        FINAL_GRO,
        EDR,
    )

    outputs_present = all(
        path.exists()
        and path.stat().st_size > 0
        for path in required_outputs
    )

    blocked_reasons = []

    if mdrun_return_code != 0:
        blocked_reasons.append(
            f"mdrun return code {mdrun_return_code}"
        )

    if not outputs_present:
        blocked_reasons.append(
            "missing Stage00 output files"
        )

    if blocked_reasons:
        decision = "BLOCKED"

        summary = {
            "stage": PREFIX_NAME,
            "decision": decision,
            "grompp_return_code": (
                grompp_return_code
            ),
            "mdrun_return_code": (
                mdrun_return_code
            ),
            "blocked_reasons": (
                " | ".join(blocked_reasons)
            ),
        }

        write_csv(
            SUMMARY_CSV,
            summary,
        )

        write_output_manifest()

        log()
        log(
            f"Stage00 decision: {decision}"
        )

        log(
            "Reasons: "
            + " | ".join(blocked_reasons)
        )

        return 2

    metrics = parse_mdrun_log(
        MDRUN_LOG
    )

    start_atoms, start_box = read_gro(
        START_GRO
    )

    final_atoms, final_box = read_gro(
        FINAL_GRO
    )

    atom_count_pass = (
        len(start_atoms)
        == len(final_atoms)
        == EXPECTED_ATOMS
    )

    coordinates_finite = (
        finite_coordinates(final_atoms)
    )

    box_finite = all(
        math.isfinite(value)
        for value in final_box
    )

    box_positive = (
        len(final_box) >= 3
        and all(
            value > 0.0
            for value in final_box[:3]
        )
    )

    hbn_displacement = group_displacement(
        start_atoms,
        final_atoms,
        *HBN_RANGE,
    )

    pyr_displacement = group_displacement(
        start_atoms,
        final_atoms,
        *PYR_RANGE,
    )

    potential_energy = metrics[
        "potential_energy_kj_mol"
    ]

    maximum_force = metrics[
        "maximum_force_kj_mol_nm"
    ]

    metrics_finite = all(
        value is not None
        and math.isfinite(float(value))
        for value in (
            potential_energy,
            maximum_force,
            metrics[
                "norm_force_kj_mol_nm"
            ],
        )
    )

    force_pass = (
        maximum_force is not None
        and float(maximum_force)
        <= EXPECTED_EMTOL * 1.000001
    )

    displacement_pass = (
        float(
            hbn_displacement[
                "max_displacement_nm"
            ]
        )
        <= MAX_ALLOWED_RESTRAINED_DISPLACEMENT_NM
        and float(
            pyr_displacement[
                "max_displacement_nm"
            ]
        )
        <= MAX_ALLOWED_RESTRAINED_DISPLACEMENT_NM
    )

    revise_reasons = []

    if not bool(metrics["converged"]):
        revise_reasons.append(
            "steepest-descent convergence criterion not reached"
        )

    if not force_pass:
        revise_reasons.append(
            "final maximum force exceeds emtol"
        )

    if not displacement_pass:
        revise_reasons.append(
            "restrained-solute displacement exceeds 0.05 nm"
        )

    if not atom_count_pass:
        blocked_reasons.append(
            "atom-count mismatch"
        )

    if not coordinates_finite:
        blocked_reasons.append(
            "non-finite final coordinates"
        )

    if not metrics_finite:
        blocked_reasons.append(
            "missing or non-finite minimization metrics"
        )

    if not box_finite or not box_positive:
        blocked_reasons.append(
            "invalid final simulation box"
        )

    if blocked_reasons:
        decision = "BLOCKED"
    elif revise_reasons:
        decision = "REVISE"
    else:
        decision = "PASS"

    summary = {
        "stage": PREFIX_NAME,
        "decision": decision,
        "grompp_return_code": (
            grompp_return_code
        ),
        "mdrun_return_code": (
            mdrun_return_code
        ),
        "atom_count": len(final_atoms),
        "atom_count_pass": (
            atom_count_pass
        ),
        "coordinates_finite": (
            coordinates_finite
        ),
        "box_finite_positive": (
            box_finite
            and box_positive
        ),
        "converged": (
            metrics["converged"]
        ),
        "convergence_message": (
            metrics[
                "convergence_message"
            ]
        ),
        "minimization_steps": (
            metrics["steps"]
        ),
        "potential_energy_kj_mol": (
            potential_energy
        ),
        "maximum_force_kj_mol_nm": (
            maximum_force
        ),
        "maximum_force_atom": (
            metrics[
                "maximum_force_atom"
            ]
        ),
        "norm_force_kj_mol_nm": (
            metrics[
                "norm_force_kj_mol_nm"
            ]
        ),
        "emtol_kj_mol_nm": (
            EXPECTED_EMTOL
        ),
        "force_pass": force_pass,
        "HBN_rms_displacement_nm": (
            hbn_displacement[
                "rms_displacement_nm"
            ]
        ),
        "HBN_max_displacement_nm": (
            hbn_displacement[
                "max_displacement_nm"
            ]
        ),
        "HBN_max_displacement_atom": (
            hbn_displacement[
                "max_displacement_atom"
            ]
        ),
        "PYR_rms_displacement_nm": (
            pyr_displacement[
                "rms_displacement_nm"
            ]
        ),
        "PYR_max_displacement_nm": (
            pyr_displacement[
                "max_displacement_nm"
            ]
        ),
        "PYR_max_displacement_atom": (
            pyr_displacement[
                "max_displacement_atom"
            ]
        ),
        "restrained_displacement_pass": (
            displacement_pass
        ),
        "revise_reasons": (
            " | ".join(revise_reasons)
        ),
        "blocked_reasons": (
            " | ".join(blocked_reasons)
        ),
    }

    write_csv(
        SUMMARY_CSV,
        summary,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Stage00 Strongly Restrained Minimization\n\n"
        )

        handle.write(
            f"- Decision: **{decision}**\n"
        )

        handle.write(
            f"- Convergence: "
            f"`{metrics['convergence_message']}`\n"
        )

        handle.write(
            f"- Potential energy: "
            f"{potential_energy} kJ/mol\n"
        )

        handle.write(
            f"- Maximum force: "
            f"{maximum_force} kJ mol^-1 nm^-1\n"
        )

        handle.write(
            f"- HBN RMS/max displacement: "
            f"{hbn_displacement['rms_displacement_nm']:.8f}/"
            f"{hbn_displacement['max_displacement_nm']:.8f} nm\n"
        )

        handle.write(
            f"- PYR RMS/max displacement: "
            f"{pyr_displacement['rms_displacement_nm']:.8f}/"
            f"{pyr_displacement['max_displacement_nm']:.8f} nm\n"
        )

        handle.write(
            f"- Atom count and finite-coordinate checks: "
            f"{'PASS' if atom_count_pass and coordinates_finite else 'FAIL'}\n"
        )

        if revise_reasons:
            handle.write(
                "\n## Revision reasons\n\n"
            )

            for reason in revise_reasons:
                handle.write(
                    f"- {reason}\n"
                )

        if blocked_reasons:
            handle.write(
                "\n## Blocking reasons\n\n"
            )

            for reason in blocked_reasons:
                handle.write(
                    f"- {reason}\n"
                )

        handle.write(
            "\nNo subsequent release stage was executed.\n"
        )

    write_output_manifest()

    log()
    log(
        "===== STAGE00 ANALYSIS ====="
    )

    log(
        f"Minimization convergence: "
        f"{'PASS' if metrics['converged'] else 'FAIL'}"
    )

    log(
        f"Convergence message: "
        f"{metrics['convergence_message']}"
    )

    log(
        f"Potential energy: "
        f"{potential_energy} kJ/mol"
    )

    log(
        f"Maximum force: "
        f"{maximum_force} kJ mol^-1 nm^-1"
    )

    log(
        "HBN RMS/max displacement: "
        f"{hbn_displacement['rms_displacement_nm']:.8f}/"
        f"{hbn_displacement['max_displacement_nm']:.8f} nm"
    )

    log(
        "PYR RMS/max displacement: "
        f"{pyr_displacement['rms_displacement_nm']:.8f}/"
        f"{pyr_displacement['max_displacement_nm']:.8f} nm"
    )

    log(
        f"Atom count: "
        f"{len(final_atoms)}/{EXPECTED_ATOMS}"
    )

    log(
        f"Finite-coordinate and box checks: "
        f"{'PASS' if coordinates_finite and box_finite and box_positive else 'FAIL'}"
    )

    log(
        f"Stage00 decision: {decision}"
    )

    log(
        "Subsequent stages executed: NO"
    )

    log(
        f"Wrote: {relative(RUN_ROOT)}"
    )

    return (
        0
        if decision == "PASS"
        else 2
    )


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            f"FATAL: {error}",
            file=sys.stderr,
            flush=True,
        )
        raise SystemExit(2)
