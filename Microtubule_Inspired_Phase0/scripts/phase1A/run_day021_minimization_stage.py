#!/usr/bin/env python3

from __future__ import annotations

import argparse
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

TOPOLOGY_ROOT = (
    PROTOCOL_ROOT
    / "protocol_inputs/topology"
)

TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    PROTOCOL_ROOT
    / "protocol_inputs/mobile_release_index.ndx"
)

STATIC_CSV = (
    PROTOCOL_ROOT
    / "static_validation/mobile_release_static_validation.csv"
)

INPUT_MANIFEST = (
    PROTOCOL_ROOT
    / "protocol_inputs/mobile_release_protocol_manifest.csv"
)

EXPECTED_ATOMS = 68320

HBN_RANGE = (1, 1680)
PYR_RANGE = (1681, 1784)

FLOAT_PATTERN = (
    r"[+-]?"
    r"(?:\d+(?:\.\d*)?|\.\d+)"
    r"(?:[Ee][+-]?\d+)?"
)


def resolve_path(value: str) -> Path:
    path = Path(value)

    if path.is_absolute():
        return path

    return PROJECT_ROOT / path


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(
                PROJECT_ROOT
            )
        )
    except ValueError:
        return str(path.resolve())


def sha256(path: Path) -> str:
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


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


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


def verify_required_inputs(
    paths: tuple[Path, ...],
) -> None:
    missing = [
        path
        for path in paths
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty required inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def verify_previous_stage(
    summary_path: Path,
) -> None:
    with summary_path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 1:
        raise RuntimeError(
            "Previous-stage summary must contain "
            "exactly one data row"
        )

    decision = rows[0].get(
        "decision",
        "",
    ).strip().upper()

    if decision != "PASS":
        raise RuntimeError(
            "Previous stage is not authorized "
            f"for continuation: {decision or 'UNKNOWN'}"
        )


def verify_static_stage(
    stage_name: str,
) -> None:
    with STATIC_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    matching = [
        row
        for row in rows
        if row.get("stage") == stage_name
    ]

    if len(matching) != 1:
        raise RuntimeError(
            f"Static-validation row not uniquely "
            f"identified for {stage_name}"
        )

    row = matching[0]

    if not parse_bool(
        row.get("grompp_pass", "")
    ):
        raise RuntimeError(
            f"Static grompp validation failed "
            f"for {stage_name}"
        )

    if not parse_bool(
        row.get(
            "restraint_entry_validation",
            "",
        )
    ):
        raise RuntimeError(
            f"Static restraint validation failed "
            f"for {stage_name}"
        )


def verify_input_manifest() -> int:
    with INPUT_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if not rows:
        raise RuntimeError(
            "Protocol input manifest is empty"
        )

    mismatches: list[str] = []

    for row in rows:
        path = PROJECT_ROOT / row["path"]

        if not path.exists():
            mismatches.append(
                f"missing: {row['path']}"
            )
            continue

        current_hash = sha256(path)

        if current_hash != row["sha256"]:
            mismatches.append(
                f"hash mismatch: {row['path']}"
            )

    if mismatches:
        raise RuntimeError(
            "Protocol manifest verification failed:\n"
            + "\n".join(mismatches)
        )

    return len(rows)


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
    list[
        tuple[
            str,
            str,
            int,
            float,
            float,
            float,
        ]
    ],
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

    natoms = int(
        lines[1].strip()
    )

    atom_lines = lines[
        2 : 2 + natoms
    ]

    if len(atom_lines) != natoms:
        raise RuntimeError(
            f"GRO atom-count mismatch: {path}"
        )

    atoms = []

    for line in atom_lines:
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom line: {path}"
            )

        atoms.append(
            (
                line[5:10].strip(),
                line[10:15].strip(),
                int(line[15:20]),
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            )
        )

    box_line = lines[
        2 + natoms
    ].split()

    box = [
        float(value)
        for value in box_line
    ]

    return atoms, box


def finite_coordinates(
    atoms,
) -> bool:
    return all(
        math.isfinite(value)
        for atom in atoms
        for value in atom[3:6]
    )


def displacement_metrics(
    source_atoms,
    target_atoms,
    atom_range: tuple[int, int],
) -> dict[str, object]:
    first_atom, last_atom = (
        atom_range
    )

    squared_sum = 0.0
    maximum = -1.0
    maximum_atom = -1

    for atom_number in range(
        first_atom,
        last_atom + 1,
    ):
        index = atom_number - 1

        source = source_atoms[index]
        target = target_atoms[index]

        if source[:3] != target[:3]:
            raise RuntimeError(
                "Atom identity/order mismatch at "
                f"atom {atom_number}"
            )

        dx = target[3] - source[3]
        dy = target[4] - source[4]
        dz = target[5] - source[5]

        displacement = math.sqrt(
            dx * dx
            + dy * dy
            + dz * dz
        )

        squared_sum += (
            displacement * displacement
        )

        if displacement > maximum:
            maximum = displacement
            maximum_atom = atom_number

    count = (
        last_atom
        - first_atom
        + 1
    )

    return {
        "atom_count": count,
        "rms_nm": math.sqrt(
            squared_sum / count
        ),
        "max_nm": maximum,
        "max_atom": maximum_atom,
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


def parse_minimization_log(
    path: Path,
) -> dict[str, object]:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    potential = last_float_match(
        text,
        rf"Potential Energy\s*=\s*"
        rf"({FLOAT_PATTERN})",
    )

    maximum_force = last_float_match(
        text,
        rf"Maximum force\s*=\s*"
        rf"({FLOAT_PATTERN})",
    )

    norm_force = last_float_match(
        text,
        rf"Norm of force\s*=\s*"
        rf"({FLOAT_PATTERN})",
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
                or "machine precision"
                in line.lower()
            )
        )
    ]

    convergence_message = (
        convergence_lines[-1]
        if convergence_lines
        else ""
    )

    converged = (
        "converged to fmax"
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
        "potential_energy": potential,
        "maximum_force": maximum_force,
        "maximum_force_atom": (
            maximum_force_atom
        ),
        "norm_force": norm_force,
        "convergence_message": (
            convergence_message
        ),
        "converged": converged,
        "steps": steps,
    }


def atom_group(
    atom_number: int | None,
) -> str:
    if atom_number is None:
        return "UNKNOWN"

    if 1 <= atom_number <= 1680:
        return "HBN"

    if 1681 <= atom_number <= 1784:
        return "PYR"

    if 1785 <= atom_number <= 68320:
        return "SOL"

    return "OUT_OF_RANGE"


def write_single_row_csv(
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


def write_manifest(
    run_root: Path,
    manifest_path: Path,
) -> None:
    paths = sorted(
        path
        for path in run_root.iterdir()
        if (
            path.is_file()
            and path != manifest_path
            and not path.name.endswith(
                "_wrapper.log"
            )
        )
    )

    rows = [
        {
            "path": relative(path),
            "size_bytes": (
                path.stat().st_size
            ),
            "sha256": sha256(path),
        }
        for path in paths
    ]

    with manifest_path.open(
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


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage-name",
        required=True,
    )

    parser.add_argument(
        "--mdp",
        required=True,
    )

    parser.add_argument(
        "--start-gro",
        required=True,
    )

    parser.add_argument(
        "--reference-gro",
        required=True,
    )

    parser.add_argument(
        "--previous-summary",
        required=True,
    )

    parser.add_argument(
        "--run-root",
        required=True,
    )

    parser.add_argument(
        "--expected-emtol",
        type=float,
        default=1000.0,
    )

    parser.add_argument(
        "--max-incremental-solute-displacement-nm",
        type=float,
        default=0.05,
    )

    parser.add_argument(
        "--max-cumulative-solute-displacement-nm",
        type=float,
        default=0.05,
    )

    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()

    stage_name = arguments.stage_name
    mdp = resolve_path(arguments.mdp)
    start_gro = resolve_path(
        arguments.start_gro
    )
    reference_gro = resolve_path(
        arguments.reference_gro
    )
    previous_summary = resolve_path(
        arguments.previous_summary
    )
    run_root = resolve_path(
        arguments.run_root
    )

    run_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    prefix = run_root / stage_name

    tpr = prefix.with_suffix(".tpr")
    final_gro = prefix.with_suffix(".gro")
    edr = prefix.with_suffix(".edr")
    mdrun_log = prefix.with_suffix(".log")

    processed_top = (
        run_root
        / f"{stage_name}_processed.top"
    )

    mdout = (
        run_root
        / f"{stage_name}_mdout.mdp"
    )

    grompp_log = (
        run_root
        / f"{stage_name}_grompp.log"
    )

    mdrun_console = (
        run_root
        / f"{stage_name}_mdrun_console.log"
    )

    commands_path = (
        run_root
        / "commands.txt"
    )

    summary_path = (
        run_root
        / f"{stage_name}_summary.csv"
    )

    report_path = (
        run_root
        / f"{stage_name.upper()}_DAY021.md"
    )

    output_manifest = (
        run_root
        / f"{stage_name}_output_manifest.csv"
    )

    verify_required_inputs(
        (
            mdp,
            start_gro,
            reference_gro,
            previous_summary,
            TOP,
            NDX,
            STATIC_CSV,
            INPUT_MANIFEST,
        )
    )

    verify_previous_stage(
        previous_summary
    )

    verify_static_stage(
        stage_name
    )

    manifest_file_count = (
        verify_input_manifest()
    )

    conflicting_outputs = [
        path
        for path in (
            tpr,
            final_gro,
            edr,
            mdrun_log,
        )
        if path.exists()
    ]

    if conflicting_outputs:
        raise RuntimeError(
            "Stage outputs already exist; "
            "nothing was overwritten:\n"
            + "\n".join(
                str(path)
                for path in conflicting_outputs
            )
        )

    gmx = find_gmx()

    grompp_command = [
        str(gmx),
        "grompp",
        "-f",
        str(mdp),
        "-c",
        str(start_gro),
        "-r",
        str(reference_gro),
        "-p",
        str(TOP),
        "-n",
        str(NDX),
        "-o",
        str(tpr),
        "-po",
        str(mdout),
        "-pp",
        str(processed_top),
        "-maxwarn",
        "0",
    ]

    mdrun_command = [
        str(gmx),
        "mdrun",
        "-deffnm",
        stage_name,
    ]

    commands_path.write_text(
        " ".join(grompp_command)
        + "\n\n"
        + " ".join(mdrun_command)
        + "\n",
        encoding="utf-8",
    )

    print(
        f"Day021 minimization stage: "
        f"{stage_name}",
        flush=True,
    )

    print(
        f"Previous stage authorization: PASS",
        flush=True,
    )

    print(
        f"Verified protocol-manifest files: "
        f"{manifest_file_count}",
        flush=True,
    )

    print(
        f"Starting coordinates: "
        f"{relative(start_gro)}",
        flush=True,
    )

    print(
        f"Restraint reference coordinates: "
        f"{relative(reference_gro)}",
        flush=True,
    )

    print(
        "\n===== GROMPP =====",
        flush=True,
    )

    grompp_return_code = run_streaming(
        grompp_command,
        TOPOLOGY_ROOT,
        grompp_log,
    )

    if (
        grompp_return_code != 0
        or not tpr.exists()
    ):
        raise RuntimeError(
            "grompp failed; no minimization started"
        )

    print(
        "\n===== MDRUN =====",
        flush=True,
    )

    mdrun_return_code = run_streaming(
        mdrun_command,
        run_root,
        mdrun_console,
    )

    required_outputs = (
        final_gro,
        edr,
        mdrun_log,
    )

    missing_outputs = [
        path
        for path in required_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if (
        mdrun_return_code != 0
        or missing_outputs
    ):
        reasons = []

        if mdrun_return_code != 0:
            reasons.append(
                f"mdrun return code "
                f"{mdrun_return_code}"
            )

        if missing_outputs:
            reasons.append(
                "missing outputs: "
                + ", ".join(
                    path.name
                    for path in missing_outputs
                )
            )

        write_single_row_csv(
            summary_path,
            {
                "stage": stage_name,
                "decision": "BLOCKED",
                "reasons": " | ".join(
                    reasons
                ),
            },
        )

        print(
            "\nStage decision: BLOCKED",
            flush=True,
        )

        return 2

    metrics = parse_minimization_log(
        mdrun_log
    )

    start_atoms, start_box = read_gro(
        start_gro
    )

    reference_atoms, reference_box = (
        read_gro(reference_gro)
    )

    final_atoms, final_box = read_gro(
        final_gro
    )

    atom_count_pass = (
        len(start_atoms)
        == len(reference_atoms)
        == len(final_atoms)
        == EXPECTED_ATOMS
    )

    finite_pass = all(
        (
            finite_coordinates(
                final_atoms
            ),
            all(
                math.isfinite(value)
                for value in final_box
            ),
        )
    )

    incremental_hbn = (
        displacement_metrics(
            start_atoms,
            final_atoms,
            HBN_RANGE,
        )
    )

    incremental_pyr = (
        displacement_metrics(
            start_atoms,
            final_atoms,
            PYR_RANGE,
        )
    )

    cumulative_hbn = (
        displacement_metrics(
            reference_atoms,
            final_atoms,
            HBN_RANGE,
        )
    )

    cumulative_pyr = (
        displacement_metrics(
            reference_atoms,
            final_atoms,
            PYR_RANGE,
        )
    )

    box_delta = max(
        abs(final - initial)
        for final, initial in zip(
            final_box,
            start_box,
        )
    )

    maximum_force = metrics[
        "maximum_force"
    ]

    potential_energy = metrics[
        "potential_energy"
    ]

    force_pass = (
        maximum_force is not None
        and math.isfinite(
            float(maximum_force)
        )
        and float(maximum_force)
        <= (
            arguments.expected_emtol
            * 1.000001
        )
    )

    incremental_displacement_pass = (
        float(
            incremental_hbn["max_nm"]
        )
        <= arguments.max_incremental_solute_displacement_nm
        and float(
            incremental_pyr["max_nm"]
        )
        <= arguments.max_incremental_solute_displacement_nm
    )

    cumulative_displacement_pass = (
        float(
            cumulative_hbn["max_nm"]
        )
        <= arguments.max_cumulative_solute_displacement_nm
        and float(
            cumulative_pyr["max_nm"]
        )
        <= arguments.max_cumulative_solute_displacement_nm
    )

    blocked_reasons: list[str] = []
    revise_reasons: list[str] = []

    if not atom_count_pass:
        blocked_reasons.append(
            "atom-count mismatch"
        )

    if not finite_pass:
        blocked_reasons.append(
            "non-finite coordinates or box"
        )

    if (
        potential_energy is None
        or not math.isfinite(
            float(potential_energy)
        )
    ):
        blocked_reasons.append(
            "invalid potential energy"
        )

    if not bool(metrics["converged"]):
        revise_reasons.append(
            "minimization did not reach Fmax criterion"
        )

    if not force_pass:
        revise_reasons.append(
            "maximum force exceeds emtol"
        )

    if not incremental_displacement_pass:
        revise_reasons.append(
            "incremental restrained-solute "
            "displacement exceeds threshold"
        )

    if not cumulative_displacement_pass:
        revise_reasons.append(
            "cumulative restrained-solute "
            "displacement exceeds threshold"
        )

    if box_delta > 1.0e-6:
        revise_reasons.append(
            "simulation box changed during EM"
        )

    if blocked_reasons:
        decision = "BLOCKED"
    elif revise_reasons:
        decision = "REVISE"
    else:
        decision = "PASS"

    maximum_force_atom = metrics[
        "maximum_force_atom"
    ]

    maximum_force_group = atom_group(
        maximum_force_atom
    )

    summary = {
        "stage": stage_name,
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
        "finite_coordinate_box_pass": (
            finite_pass
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
            maximum_force_atom
        ),
        "maximum_force_group": (
            maximum_force_group
        ),
        "norm_force_kj_mol_nm": (
            metrics["norm_force"]
        ),
        "emtol_kj_mol_nm": (
            arguments.expected_emtol
        ),
        "force_pass": force_pass,
        "HBN_incremental_rms_nm": (
            incremental_hbn["rms_nm"]
        ),
        "HBN_incremental_max_nm": (
            incremental_hbn["max_nm"]
        ),
        "PYR_incremental_rms_nm": (
            incremental_pyr["rms_nm"]
        ),
        "PYR_incremental_max_nm": (
            incremental_pyr["max_nm"]
        ),
        "HBN_cumulative_rms_nm": (
            cumulative_hbn["rms_nm"]
        ),
        "HBN_cumulative_max_nm": (
            cumulative_hbn["max_nm"]
        ),
        "PYR_cumulative_rms_nm": (
            cumulative_pyr["rms_nm"]
        ),
        "PYR_cumulative_max_nm": (
            cumulative_pyr["max_nm"]
        ),
        "incremental_displacement_pass": (
            incremental_displacement_pass
        ),
        "cumulative_displacement_pass": (
            cumulative_displacement_pass
        ),
        "maximum_box_change_nm": (
            box_delta
        ),
        "revise_reasons": (
            " | ".join(revise_reasons)
        ),
        "blocked_reasons": (
            " | ".join(blocked_reasons)
        ),
    }

    write_single_row_csv(
        summary_path,
        summary,
    )

    with report_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"# Day021 {stage_name}\n\n"
        )

        handle.write(
            f"- Decision: **{decision}**\n"
        )

        handle.write(
            f"- Convergence: "
            f"{metrics['convergence_message']}\n"
        )

        handle.write(
            f"- Potential energy: "
            f"{potential_energy} kJ/mol\n"
        )

        handle.write(
            f"- Maximum force: "
            f"{maximum_force} kJ mol^-1 nm^-1 "
            f"on atom {maximum_force_atom} "
            f"({maximum_force_group})\n"
        )

        handle.write(
            f"- HBN incremental RMS/max: "
            f"{incremental_hbn['rms_nm']:.8f}/"
            f"{incremental_hbn['max_nm']:.8f} nm\n"
        )

        handle.write(
            f"- PYR incremental RMS/max: "
            f"{incremental_pyr['rms_nm']:.8f}/"
            f"{incremental_pyr['max_nm']:.8f} nm\n"
        )

        handle.write(
            f"- HBN cumulative RMS/max: "
            f"{cumulative_hbn['rms_nm']:.8f}/"
            f"{cumulative_hbn['max_nm']:.8f} nm\n"
        )

        handle.write(
            f"- PYR cumulative RMS/max: "
            f"{cumulative_pyr['rms_nm']:.8f}/"
            f"{cumulative_pyr['max_nm']:.8f} nm\n"
        )

        handle.write(
            "\nNo subsequent protocol stage was executed.\n"
        )

    write_manifest(
        run_root,
        output_manifest,
    )

    print(
        "\n===== MINIMIZATION ANALYSIS =====",
        flush=True,
    )

    print(
        f"Convergence: "
        f"{'PASS' if metrics['converged'] else 'FAIL'}",
        flush=True,
    )

    print(
        f"Convergence message: "
        f"{metrics['convergence_message']}",
        flush=True,
    )

    print(
        f"Potential energy: "
        f"{potential_energy} kJ/mol",
        flush=True,
    )

    print(
        f"Maximum force: "
        f"{maximum_force} kJ mol^-1 nm^-1 "
        f"on atom {maximum_force_atom} "
        f"({maximum_force_group})",
        flush=True,
    )

    print(
        "HBN incremental RMS/max: "
        f"{incremental_hbn['rms_nm']:.8f}/"
        f"{incremental_hbn['max_nm']:.8f} nm",
        flush=True,
    )

    print(
        "PYR incremental RMS/max: "
        f"{incremental_pyr['rms_nm']:.8f}/"
        f"{incremental_pyr['max_nm']:.8f} nm",
        flush=True,
    )

    print(
        "HBN cumulative RMS/max: "
        f"{cumulative_hbn['rms_nm']:.8f}/"
        f"{cumulative_hbn['max_nm']:.8f} nm",
        flush=True,
    )

    print(
        "PYR cumulative RMS/max: "
        f"{cumulative_pyr['rms_nm']:.8f}/"
        f"{cumulative_pyr['max_nm']:.8f} nm",
        flush=True,
    )

    print(
        f"Atom count: "
        f"{len(final_atoms)}/{EXPECTED_ATOMS}",
        flush=True,
    )

    print(
        f"Finite-coordinate and box checks: "
        f"{'PASS' if finite_pass else 'FAIL'}",
        flush=True,
    )

    print(
        f"Stage decision: {decision}",
        flush=True,
    )

    print(
        "Subsequent stages executed: NO",
        flush=True,
    )

    print(
        f"Wrote: {relative(run_root)}",
        flush=True,
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
