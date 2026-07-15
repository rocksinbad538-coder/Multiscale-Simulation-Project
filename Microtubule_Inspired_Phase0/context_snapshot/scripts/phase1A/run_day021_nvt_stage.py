#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import re
import shutil
import statistics
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

ENERGY_TERMS = {
    "Temperature": "temperature",
    "Potential": "potential",
    "Total-Energy": "total_energy",
    "Pressure": "pressure",
}


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


def require_files(
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
    allow_revise: bool = False,
) -> str:
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
            "exactly one row"
        )

    decision = rows[0].get(
        "decision",
        "",
    ).strip().upper()

    allowed = {"PASS"}

    if allow_revise:
        allowed.add("REVISE")

    if decision not in allowed:
        raise RuntimeError(
            "Previous stage is not authorized "
            f"for continuation: {decision or 'UNKNOWN'}"
        )

    return decision


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

    matches = [
        row
        for row in rows
        if row.get("stage") == stage_name
    ]

    if len(matches) != 1:
        raise RuntimeError(
            "Static-validation row was not uniquely "
            f"identified for {stage_name}"
        )

    row = matches[0]

    if not parse_bool(
        row.get("grompp_pass", "")
    ):
        raise RuntimeError(
            f"Static grompp validation failed for {stage_name}"
        )

    if not parse_bool(
        row.get(
            "restraint_entry_validation",
            "",
        )
    ):
        raise RuntimeError(
            f"Static restraint validation failed for {stage_name}"
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

        if sha256(path) != row["sha256"]:
            mismatches.append(
                f"hash mismatch: {row['path']}"
            )

    if mismatches:
        raise RuntimeError(
            "Protocol input-manifest verification failed:\n"
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
    list[dict[str, object]],
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

    atoms: list[
        dict[str, object]
    ] = []

    for line in atom_lines:
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom record: {path}"
            )

        atom = {
            "residue": line[5:10].strip(),
            "name": line[10:15].strip(),
            "number": int(line[15:20]),
            "x": float(line[20:28]),
            "y": float(line[28:36]),
            "z": float(line[36:44]),
            "vx": None,
            "vy": None,
            "vz": None,
        }

        if len(line) >= 68:
            atom["vx"] = float(
                line[44:52]
            )
            atom["vy"] = float(
                line[52:60]
            )
            atom["vz"] = float(
                line[60:68]
            )

        atoms.append(atom)

    box = [
        float(value)
        for value in lines[
            2 + natoms
        ].split()
    ]

    return atoms, box


def finite_coordinates(
    atoms: list[dict[str, object]],
) -> bool:
    return all(
        math.isfinite(
            float(atom[key])
        )
        for atom in atoms
        for key in (
            "x",
            "y",
            "z",
        )
    )


def displacement_metrics(
    source_atoms: list[dict[str, object]],
    target_atoms: list[dict[str, object]],
    atom_range: tuple[int, int],
) -> dict[str, object]:
    first_atom, last_atom = atom_range

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

        identity_source = (
            source["residue"],
            source["name"],
            source["number"],
        )

        identity_target = (
            target["residue"],
            target["name"],
            target["number"],
        )

        if identity_source != identity_target:
            raise RuntimeError(
                "Atom identity/order mismatch at "
                f"atom {atom_number}"
            )

        dx = (
            float(target["x"])
            - float(source["x"])
        )
        dy = (
            float(target["y"])
            - float(source["y"])
        )
        dz = (
            float(target["z"])
            - float(source["z"])
        )

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
        "count": count,
        "rms_nm": math.sqrt(
            squared_sum / count
        ),
        "max_nm": maximum,
        "max_atom": maximum_atom,
    }


def velocity_summary(
    atoms: list[dict[str, object]],
) -> dict[str, dict[str, float | int]]:
    accumulator: dict[
        str,
        dict[str, float | int],
    ] = {}

    for atom in atoms:
        residue = str(
            atom["residue"]
        )

        vx = atom["vx"]
        vy = atom["vy"]
        vz = atom["vz"]

        if (
            vx is None
            or vy is None
            or vz is None
        ):
            raise RuntimeError(
                "Final GRO does not contain velocities"
            )

        speed = math.sqrt(
            float(vx) ** 2
            + float(vy) ** 2
            + float(vz) ** 2
        )

        group = accumulator.setdefault(
            residue,
            {
                "count": 0,
                "nonzero": 0,
                "speed2_sum": 0.0,
                "max_speed": 0.0,
            },
        )

        group["count"] = (
            int(group["count"])
            + 1
        )

        if speed > 1.0e-8:
            group["nonzero"] = (
                int(group["nonzero"])
                + 1
            )

        group["speed2_sum"] = (
            float(group["speed2_sum"])
            + speed * speed
        )

        group["max_speed"] = max(
            float(group["max_speed"]),
            speed,
        )

    result: dict[
        str,
        dict[str, float | int],
    ] = {}

    for residue, group in accumulator.items():
        count = int(
            group["count"]
        )
        nonzero = int(
            group["nonzero"]
        )

        result[residue] = {
            "count": count,
            "nonzero": nonzero,
            "nonzero_fraction": (
                nonzero / count
            ),
            "rms_speed_nm_ps": math.sqrt(
                float(group["speed2_sum"])
                / count
            ),
            "max_speed_nm_ps": float(
                group["max_speed"]
            ),
        }

    return result


def extract_energy_series(
    gmx: Path,
    edr: Path,
    term: str,
    output_xvg: Path,
    output_log: Path,
    cwd: Path,
) -> list[
    tuple[float, float]
]:
    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    result = subprocess.run(
        [
            str(gmx),
            "energy",
            "-f",
            str(edr),
            "-o",
            str(output_xvg),
        ],
        cwd=cwd,
        env=environment,
        input=f"{term}\n0\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    output_log.write_text(
        result.stdout,
        encoding="utf-8",
    )

    if (
        result.returncode != 0
        or not output_xvg.exists()
        or output_xvg.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Could not extract energy term: {term}"
        )

    series: list[
        tuple[float, float]
    ] = []

    for raw_line in output_xvg.read_text(
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

        time_ps = float(fields[0])
        value = float(fields[1])

        if not (
            math.isfinite(time_ps)
            and math.isfinite(value)
        ):
            raise RuntimeError(
                f"Non-finite value in {term} series"
            )

        series.append(
            (
                time_ps,
                value,
            )
        )

    if not series:
        raise RuntimeError(
            f"No data points extracted for {term}"
        )

    return series


def series_statistics(
    series: list[
        tuple[float, float]
    ],
) -> dict[str, float | int]:
    values = [
        value
        for _, value in series
    ]

    return {
        "points": len(values),
        "first": values[0],
        "last": values[-1],
        "mean": statistics.fmean(
            values
        ),
        "std": (
            statistics.pstdev(values)
            if len(values) > 1
            else 0.0
        ),
        "min": min(values),
        "max": max(values),
    }


def instability_signatures(
    paths: tuple[Path, ...],
) -> list[str]:
    serious_patterns = {
        "LINCS_WARNING": re.compile(
            r"LINCS WARNING",
            re.IGNORECASE,
        ),
        "SETTLE_FAILURE": re.compile(
            r"(cannot be settled|SETTLE.*(?:error|failed))",
            re.IGNORECASE,
        ),
        "SHAKE_FAILURE": re.compile(
            r"SHAKE.*(?:failed|did not converge)",
            re.IGNORECASE,
        ),
        "FATAL_ERROR": re.compile(
            r"Fatal error",
            re.IGNORECASE,
        ),
    }

    nonfinite_pattern = re.compile(
        r"(?<![A-Za-z])"
        r"(?:nan|[-+]?inf(?:inity)?)"
        r"(?![A-Za-z])",
        re.IGNORECASE,
    )

    harmless_epsilon_rf = re.compile(
        r"^\s*epsilon-rf\s*=\s*"
        r"(?:inf|infinity)\s*$",
        re.IGNORECASE,
    )

    detected: list[str] = []

    for path in paths:
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            for label, pattern in (
                serious_patterns.items()
            ):
                if pattern.search(line):
                    detected.append(
                        f"{label}:{path.name}:"
                        f"{line_number}"
                    )

            if harmless_epsilon_rf.match(line):
                continue

            if nonfinite_pattern.search(line):
                detected.append(
                    f"NONFINITE_VALUE:"
                    f"{path.name}:{line_number}"
                )

    return sorted(
        set(detected)
    )


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
    files = sorted(
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
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in files
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
        "--checkpoint",
        default="",
    )
    parser.add_argument(
        "--temperature-lower",
        type=float,
        default=270.0,
    )
    parser.add_argument(
        "--temperature-upper",
        type=float,
        default=330.0,
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
    parser.add_argument(
        "--require-mobile-velocities",
        action="store_true",
    )

    parser.add_argument(
        "--allow-previous-revise",
        action="store_true",
        help=(
            "Permit an adaptive hold stage to continue "
            "from a scientifically REVISE previous stage."
        ),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_arguments()

    stage_name = args.stage_name

    mdp = resolve_path(args.mdp)
    start_gro = resolve_path(
        args.start_gro
    )
    reference_gro = resolve_path(
        args.reference_gro
    )
    previous_summary = resolve_path(
        args.previous_summary
    )
    run_root = resolve_path(
        args.run_root
    )

    checkpoint = (
        resolve_path(args.checkpoint)
        if args.checkpoint
        else None
    )

    run_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    tpr = run_root / f"{stage_name}.tpr"
    final_gro = run_root / f"{stage_name}.gro"
    edr = run_root / f"{stage_name}.edr"
    cpt = run_root / f"{stage_name}.cpt"
    mdrun_log = run_root / f"{stage_name}.log"
    xtc = run_root / f"{stage_name}.xtc"

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

    required = (
        mdp,
        start_gro,
        reference_gro,
        previous_summary,
        TOP,
        NDX,
        STATIC_CSV,
        INPUT_MANIFEST,
    )

    if checkpoint is not None:
        required = (
            *required,
            checkpoint,
        )

    require_files(required)

    previous_decision = verify_previous_stage(
        previous_summary,
        allow_revise=args.allow_previous_revise,
    )
    verify_static_stage(
        stage_name
    )

    manifest_count = (
        verify_input_manifest()
    )

    conflicts = [
        path
        for path in (
            tpr,
            final_gro,
            edr,
            cpt,
            mdrun_log,
        )
        if path.exists()
    ]

    if conflicts:
        raise RuntimeError(
            "Stage outputs already exist; "
            "nothing was overwritten:\n"
            + "\n".join(
                str(path)
                for path in conflicts
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

    if checkpoint is not None:
        grompp_command.extend(
            [
                "-t",
                str(checkpoint),
            ]
        )

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
        f"Day021 NVT stage: {stage_name}",
        flush=True,
    )
    print(
        "Previous stage authorization: "
        f"{previous_decision}"
        + (
            " -> adaptive continuation"
            if previous_decision == "REVISE"
            else ""
        ),
        flush=True,
    )
    print(
        f"Verified protocol-manifest files: "
        f"{manifest_count}",
        flush=True,
    )
    print(
        f"Starting coordinates: "
        f"{relative(start_gro)}",
        flush=True,
    )
    print(
        f"Checkpoint input: "
        f"{relative(checkpoint) if checkpoint else 'NONE'}",
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
        or tpr.stat().st_size == 0
    ):
        raise RuntimeError(
            "grompp failed; NVT was not started"
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
        cpt,
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
                "blocked_reasons": (
                    " | ".join(reasons)
                ),
            },
        )

        print(
            "\nStage decision: BLOCKED",
            flush=True,
        )

        return 2

    energy_statistics: dict[
        str,
        dict[str, float | int],
    ] = {}

    for term, stem in ENERGY_TERMS.items():
        xvg_path = (
            run_root
            / f"{stage_name}_{stem}.xvg"
        )
        energy_log = (
            run_root
            / f"{stage_name}_{stem}_energy.log"
        )

        series = extract_energy_series(
            gmx,
            edr,
            term,
            xvg_path,
            energy_log,
            run_root,
        )

        energy_statistics[stem] = (
            series_statistics(series)
        )

    start_atoms, start_box = read_gro(
        start_gro
    )
    reference_atoms, _ = read_gro(
        reference_gro
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

    box_change = max(
        abs(final - initial)
        for final, initial in zip(
            final_box,
            start_box,
        )
    )

    hbn_incremental = (
        displacement_metrics(
            start_atoms,
            final_atoms,
            HBN_RANGE,
        )
    )
    pyr_incremental = (
        displacement_metrics(
            start_atoms,
            final_atoms,
            PYR_RANGE,
        )
    )
    hbn_cumulative = (
        displacement_metrics(
            reference_atoms,
            final_atoms,
            HBN_RANGE,
        )
    )
    pyr_cumulative = (
        displacement_metrics(
            reference_atoms,
            final_atoms,
            PYR_RANGE,
        )
    )

    velocities = velocity_summary(
        final_atoms
    )

    temperature = energy_statistics[
        "temperature"
    ]

    temperature_pass = all(
        (
            float(temperature["mean"])
            >= args.temperature_lower,
            float(temperature["mean"])
            <= args.temperature_upper,
            float(temperature["last"])
            >= args.temperature_lower,
            float(temperature["last"])
            <= args.temperature_upper,
            float(temperature["min"])
            >= 250.0,
            float(temperature["max"])
            <= 350.0,
        )
    )

    incremental_pass = all(
        (
            float(
                hbn_incremental["max_nm"]
            )
            <= args.max_incremental_solute_displacement_nm,
            float(
                pyr_incremental["max_nm"]
            )
            <= args.max_incremental_solute_displacement_nm,
        )
    )

    cumulative_pass = all(
        (
            float(
                hbn_cumulative["max_nm"]
            )
            <= args.max_cumulative_solute_displacement_nm,
            float(
                pyr_cumulative["max_nm"]
            )
            <= args.max_cumulative_solute_displacement_nm,
        )
    )

    velocity_pass = True

    if args.require_mobile_velocities:
        required_residues = (
            "HBN",
            "PYR",
            "SOL",
        )

        velocity_pass = all(
            residue in velocities
            and float(
                velocities[residue][
                    "nonzero_fraction"
                ]
            )
            >= 0.95
            for residue in required_residues
        )

    signatures = instability_signatures(
        (
            mdrun_log,
            mdrun_console,
        )
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

    if signatures:
        blocked_reasons.append(
            "instability signatures: "
            + ",".join(signatures)
        )

    if not temperature_pass:
        revise_reasons.append(
            "temperature outside acceptance window"
        )

    if not incremental_pass:
        revise_reasons.append(
            "incremental solute displacement "
            "exceeds threshold"
        )

    if not cumulative_pass:
        revise_reasons.append(
            "cumulative solute displacement "
            "exceeds threshold"
        )

    if not velocity_pass:
        revise_reasons.append(
            "mobile velocities were not established "
            "for HBN/PYR/SOL"
        )

    if box_change > 1.0e-6:
        revise_reasons.append(
            "box changed during NVT"
        )

    if blocked_reasons:
        decision = "BLOCKED"
    elif revise_reasons:
        decision = "REVISE"
    else:
        decision = "PASS"

    summary: dict[str, object] = {
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
        "temperature_points": (
            temperature["points"]
        ),
        "temperature_first_K": (
            temperature["first"]
        ),
        "temperature_last_K": (
            temperature["last"]
        ),
        "temperature_mean_K": (
            temperature["mean"]
        ),
        "temperature_std_K": (
            temperature["std"]
        ),
        "temperature_min_K": (
            temperature["min"]
        ),
        "temperature_max_K": (
            temperature["max"]
        ),
        "temperature_pass": (
            temperature_pass
        ),
        "potential_first_kj_mol": (
            energy_statistics[
                "potential"
            ]["first"]
        ),
        "potential_last_kj_mol": (
            energy_statistics[
                "potential"
            ]["last"]
        ),
        "potential_mean_kj_mol": (
            energy_statistics[
                "potential"
            ]["mean"]
        ),
        "total_energy_first_kj_mol": (
            energy_statistics[
                "total_energy"
            ]["first"]
        ),
        "total_energy_last_kj_mol": (
            energy_statistics[
                "total_energy"
            ]["last"]
        ),
        "pressure_mean_bar": (
            energy_statistics[
                "pressure"
            ]["mean"]
        ),
        "pressure_min_bar": (
            energy_statistics[
                "pressure"
            ]["min"]
        ),
        "pressure_max_bar": (
            energy_statistics[
                "pressure"
            ]["max"]
        ),
        "HBN_incremental_rms_nm": (
            hbn_incremental["rms_nm"]
        ),
        "HBN_incremental_max_nm": (
            hbn_incremental["max_nm"]
        ),
        "PYR_incremental_rms_nm": (
            pyr_incremental["rms_nm"]
        ),
        "PYR_incremental_max_nm": (
            pyr_incremental["max_nm"]
        ),
        "HBN_cumulative_rms_nm": (
            hbn_cumulative["rms_nm"]
        ),
        "HBN_cumulative_max_nm": (
            hbn_cumulative["max_nm"]
        ),
        "PYR_cumulative_rms_nm": (
            pyr_cumulative["rms_nm"]
        ),
        "PYR_cumulative_max_nm": (
            pyr_cumulative["max_nm"]
        ),
        "incremental_displacement_pass": (
            incremental_pass
        ),
        "cumulative_displacement_pass": (
            cumulative_pass
        ),
        "maximum_box_change_nm": (
            box_change
        ),
        "velocity_validation_required": (
            args.require_mobile_velocities
        ),
        "velocity_pass": (
            velocity_pass
        ),
        "HBN_nonzero_velocity_fraction": (
            velocities.get(
                "HBN",
                {},
            ).get(
                "nonzero_fraction",
                "",
            )
        ),
        "PYR_nonzero_velocity_fraction": (
            velocities.get(
                "PYR",
                {},
            ).get(
                "nonzero_fraction",
                "",
            )
        ),
        "SOL_nonzero_velocity_fraction": (
            velocities.get(
                "SOL",
                {},
            ).get(
                "nonzero_fraction",
                "",
            )
        ),
        "instability_signatures": (
            ",".join(signatures)
        ),
        "revise_reasons": (
            " | ".join(revise_reasons)
        ),
        "blocked_reasons": (
            " | ".join(blocked_reasons)
        ),
        "checkpoint_path": (
            relative(cpt)
        ),
        "trajectory_path": (
            relative(xtc)
            if xtc.exists()
            else ""
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
            f"- Temperature mean/last: "
            f"{temperature['mean']:.4f}/"
            f"{temperature['last']:.4f} K\n"
        )
        handle.write(
            f"- Temperature range: "
            f"{temperature['min']:.4f}–"
            f"{temperature['max']:.4f} K\n"
        )
        handle.write(
            f"- Potential energy first/last: "
            f"{energy_statistics['potential']['first']}/"
            f"{energy_statistics['potential']['last']} kJ/mol\n"
        )
        handle.write(
            f"- HBN incremental RMS/max: "
            f"{hbn_incremental['rms_nm']:.8f}/"
            f"{hbn_incremental['max_nm']:.8f} nm\n"
        )
        handle.write(
            f"- PYR incremental RMS/max: "
            f"{pyr_incremental['rms_nm']:.8f}/"
            f"{pyr_incremental['max_nm']:.8f} nm\n"
        )
        handle.write(
            f"- HBN cumulative RMS/max: "
            f"{hbn_cumulative['rms_nm']:.8f}/"
            f"{hbn_cumulative['max_nm']:.8f} nm\n"
        )
        handle.write(
            f"- PYR cumulative RMS/max: "
            f"{pyr_cumulative['rms_nm']:.8f}/"
            f"{pyr_cumulative['max_nm']:.8f} nm\n"
        )
        handle.write(
            f"- Instability signatures: "
            f"{signatures or 'none'}\n"
        )
        handle.write(
            "\nNo subsequent protocol stage was executed.\n"
        )

    write_manifest(
        run_root,
        output_manifest,
    )

    print(
        "\n===== NVT ANALYSIS =====",
        flush=True,
    )
    print(
        "Temperature mean/std/min/max/last: "
        f"{temperature['mean']:.4f}/"
        f"{temperature['std']:.4f}/"
        f"{temperature['min']:.4f}/"
        f"{temperature['max']:.4f}/"
        f"{temperature['last']:.4f} K",
        flush=True,
    )
    print(
        "Potential energy first/last: "
        f"{energy_statistics['potential']['first']}/"
        f"{energy_statistics['potential']['last']} kJ/mol",
        flush=True,
    )
    print(
        "Pressure mean/min/max: "
        f"{energy_statistics['pressure']['mean']:.4f}/"
        f"{energy_statistics['pressure']['min']:.4f}/"
        f"{energy_statistics['pressure']['max']:.4f} bar",
        flush=True,
    )
    print(
        "HBN incremental RMS/max: "
        f"{hbn_incremental['rms_nm']:.8f}/"
        f"{hbn_incremental['max_nm']:.8f} nm",
        flush=True,
    )
    print(
        "PYR incremental RMS/max: "
        f"{pyr_incremental['rms_nm']:.8f}/"
        f"{pyr_incremental['max_nm']:.8f} nm",
        flush=True,
    )
    print(
        "HBN cumulative RMS/max: "
        f"{hbn_cumulative['rms_nm']:.8f}/"
        f"{hbn_cumulative['max_nm']:.8f} nm",
        flush=True,
    )
    print(
        "PYR cumulative RMS/max: "
        f"{pyr_cumulative['rms_nm']:.8f}/"
        f"{pyr_cumulative['max_nm']:.8f} nm",
        flush=True,
    )

    for residue in (
        "HBN",
        "PYR",
        "SOL",
    ):
        group = velocities.get(
            residue
        )

        if group:
            print(
                f"{residue} nonzero velocities: "
                f"{group['nonzero']}/"
                f"{group['count']} "
                f"({group['nonzero_fraction']:.6f})",
                flush=True,
            )

    print(
        f"Instability signatures: "
        f"{signatures or 'NONE'}",
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
