#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_md_confined_water_input_audit"
)

INVENTORY_CSV = OUTPUT_ROOT / "file_inventory.csv"
STRUCTURE_CSV = OUTPUT_ROOT / "structure_composition.csv"
INDEX_CSV = OUTPUT_ROOT / "index_groups.csv"
PARAMETER_CSV = OUTPUT_ROOT / "md_parameter_summary.csv"
REPORT_MD = OUTPUT_ROOT / "MD_CONFINED_WATER_INPUT_AUDIT_DAY020.md"

TEXT_SUFFIXES = {
    ".mdp",
    ".log",
    ".top",
    ".itp",
    ".ndx",
    ".txt",
    ".out",
    ".err",
}

PARAMETER_KEYS = (
    "integrator",
    "dt",
    "nsteps",
    "continuation",
    "constraints",
    "constraint_algorithm",
    "cutoff-scheme",
    "nstlist",
    "rlist",
    "coulombtype",
    "rcoulomb",
    "vdwtype",
    "rvdw",
    "tcoupl",
    "tc-grps",
    "tau-t",
    "ref-t",
    "pcoupl",
    "ref-p",
    "compressibility",
    "pbc",
    "comm-mode",
    "comm-grps",
    "freezegrps",
    "freezedim",
    "define",
    "gen-vel",
    "gen-temp",
    "gen-seed",
    "nstxout",
    "nstvout",
    "nstfout",
    "nstxout-compressed",
    "compressed-x-precision",
)

WATER_RESIDUE_NAMES = {
    "SOL",
    "HOH",
    "WAT",
    "TIP3",
    "TIP3P",
    "SPC",
    "SPCE",
}

ION_RESIDUE_NAMES = {
    "NA",
    "CL",
    "K",
    "CA",
    "MG",
    "ZN",
    "SOD",
    "CLA",
    "POT",
}


def log(message: str = "") -> None:
    print(message, flush=True)


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        rows = [{"status": "no_rows"}]

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


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def human_size(size_bytes: int) -> str:
    units = (
        "B",
        "KiB",
        "MiB",
        "GiB",
        "TiB",
    )

    value = float(size_bytes)

    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.3f} {unit}"

        value /= 1024.0

    return f"{size_bytes} B"


def discover_files() -> list[Path]:
    if not RUN_ROOT.exists():
        raise RuntimeError(
            f"Accepted run directory does not exist:\n"
            f"{RUN_ROOT}"
        )

    return sorted(
        path
        for path in RUN_ROOT.rglob("*")
        if path.is_file()
    )


def build_inventory(
    files: list[Path],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for path in files:
        stat = path.stat()

        rows.append(
            {
                "path": relative(path),
                "suffix": path.suffix.lower(),
                "size_bytes": stat.st_size,
                "size_human": human_size(
                    stat.st_size
                ),
            }
        )

    return rows


def select_primary_file(
    files: list[Path],
    suffixes: tuple[str, ...],
) -> Path | None:
    candidates = [
        path
        for path in files
        if path.suffix.lower() in suffixes
    ]

    if not candidates:
        return None

    return max(
        candidates,
        key=lambda path: path.stat().st_size,
    )


def find_gromacs() -> str | None:
    configured = os.environ.get(
        "GMX_BIN"
    )

    if configured:
        resolved = shutil.which(
            configured
        )

        if resolved is not None:
            return resolved

    for candidate in (
        "gmx",
        "gmx_mpi",
        "gmx_d",
        "gmx_mpi_d",
    ):
        resolved = shutil.which(
            candidate
        )

        if resolved is not None:
            return resolved

    return None


def run_command(
    command: list[str],
    output_path: Path,
    timeout_seconds: int = 180,
) -> tuple[int, str]:
    try:
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )

        output = result.stdout

        output_path.write_text(
            output,
            encoding="utf-8",
        )

        return result.returncode, output

    except subprocess.TimeoutExpired as exc:
        output = (
            f"Command timed out after "
            f"{timeout_seconds} seconds.\n"
            f"Command: {' '.join(command)}\n"
        )

        if exc.stdout:
            output += str(exc.stdout)

        output_path.write_text(
            output,
            encoding="utf-8",
        )

        return 124, output


def parse_mdp_lines(
    text: str,
    source: Path,
) -> list[dict[str, object]]:
    parameters: dict[str, str] = {}

    for raw_line in text.splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line or "=" not in line:
            continue

        key, raw_value = line.split(
            "=",
            1,
        )

        key = key.strip().lower()
        parameter_value = (
            raw_value.strip()
        )

        if key in PARAMETER_KEYS:
            parameters[key] = (
                parameter_value
            )

    return [
        {
            "source": relative(source),
            "parameter": key,
            "value": parameters.get(
                key,
                "",
            ),
        }
        for key in PARAMETER_KEYS
    ]


def parse_tpr_dump_parameters(
    text: str,
    source_label: str,
) -> list[dict[str, object]]:
    normalized_keys = {
        key.lower(): key
        for key in PARAMETER_KEYS
    }

    found: dict[str, str] = {}

    for raw_line in text.splitlines():
        line = raw_line.strip()

        match = re.match(
            r"([A-Za-z0-9_-]+)\s*=\s*(.*)$",
            line,
        )

        if match is None:
            continue

        key = match.group(1).lower()
        parameter_value = (
            match.group(2).strip()
        )

        if key in normalized_keys:
            found[key] = parameter_value

    return [
        {
            "source": source_label,
            "parameter": key,
            "value": found.get(
                key,
                "",
            ),
        }
        for key in PARAMETER_KEYS
    ]


def parse_gro(
    path: Path,
) -> tuple[
    list[dict[str, object]],
    int,
    int,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Invalid GRO file: {path}"
        )

    try:
        declared_atoms = int(
            lines[1].strip()
        )
    except ValueError as exc:
        raise RuntimeError(
            f"Could not parse atom count "
            f"from {path}"
        ) from exc

    atom_lines = lines[
        2 : 2 + declared_atoms
    ]

    if len(atom_lines) != declared_atoms:
        raise RuntimeError(
            f"GRO atom count mismatch in {path}: "
            f"declared {declared_atoms}, "
            f"read {len(atom_lines)}"
        )

    residue_atoms: Counter[
        tuple[int, str]
    ] = Counter()

    residue_names: dict[
        tuple[int, str],
        str,
    ] = {}

    atom_name_counts: Counter[
        tuple[str, str]
    ] = Counter()

    for line in atom_lines:
        if len(line) < 20:
            continue

        try:
            residue_number = int(
                line[0:5]
            )
        except ValueError:
            continue

        residue_name = (
            line[5:10].strip()
        )

        atom_name = (
            line[10:15].strip()
        )

        residue_key = (
            residue_number,
            residue_name,
        )

        residue_atoms[
            residue_key
        ] += 1

        residue_names[
            residue_key
        ] = residue_name

        atom_name_counts[
            (
                residue_name,
                atom_name,
            )
        ] += 1

    residue_count_by_name: Counter[
        str
    ] = Counter(
        residue_names.values()
    )

    atom_count_by_residue: Counter[
        str
    ] = Counter()

    for (
        residue_name,
        _atom_name,
    ), count in atom_name_counts.items():
        atom_count_by_residue[
            residue_name
        ] += count

    rows: list[dict[str, object]] = []

    for residue_name in sorted(
        residue_count_by_name
    ):
        if residue_name in (
            WATER_RESIDUE_NAMES
        ):
            category = "water"
        elif residue_name in (
            ION_RESIDUE_NAMES
        ):
            category = "ion"
        else:
            category = "solute_or_other"

        rows.append(
            {
                "source": relative(path),
                "residue_name": residue_name,
                "category": category,
                "residue_count": (
                    residue_count_by_name[
                        residue_name
                    ]
                ),
                "atom_count": (
                    atom_count_by_residue[
                        residue_name
                    ]
                ),
            }
        )

    return (
        rows,
        declared_atoms,
        len(residue_names),
    )


def parse_ndx(
    path: Path,
) -> list[dict[str, object]]:
    groups: list[
        dict[str, object]
    ] = []

    current_name: str | None = None
    current_indices: list[int] = []

    def flush() -> None:
        nonlocal current_name
        nonlocal current_indices

        if current_name is None:
            return

        groups.append(
            {
                "source": relative(path),
                "group_name": current_name,
                "atom_count": len(
                    current_indices
                ),
                "minimum_atom_index": (
                    min(current_indices)
                    if current_indices
                    else ""
                ),
                "maximum_atom_index": (
                    max(current_indices)
                    if current_indices
                    else ""
                ),
            }
        )

        current_name = None
        current_indices = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if (
            line.startswith("[")
            and line.endswith("]")
        ):
            flush()

            current_name = (
                line[1:-1].strip()
            )

            continue

        if current_name is None:
            continue

        for token in line.split():
            try:
                current_indices.append(
                    int(token)
                )
            except ValueError:
                pass

    flush()

    return groups


def search_freeze_settings(
    files: list[Path],
) -> list[
    tuple[str, str]
]:
    matches: list[
        tuple[str, str]
    ] = []

    patterns = (
        "freezegrps",
        "freezedim",
        "freeze groups",
        "frozen atoms",
    )

    for path in files:
        if path.suffix.lower() not in (
            TEXT_SUFFIXES
        ):
            continue

        try:
            text = path.read_text(
                encoding="utf-8",
                errors="replace",
            )
        except OSError:
            continue

        for line in text.splitlines():
            lower = line.lower()

            if any(
                pattern in lower
                for pattern in patterns
            ):
                matches.append(
                    (
                        relative(path),
                        line.strip(),
                    )
                )

    return matches


def parse_trajectory_check(
    text: str,
) -> dict[str, object]:
    last_frame: int | str = ""
    last_time_ps: float | str = ""
    first_time_ps: float | str = ""
    timestep_ps: float | str = ""

    frame_matches = re.findall(
        r"(?:Reading|Last)\s+frame\s+"
        r"(\d+)\s+time\s+"
        r"([-+0-9.eE]+)",
        text,
    )

    if frame_matches:
        first_time_ps = float(
            frame_matches[0][1]
        )

        last_frame = int(
            frame_matches[-1][0]
        )

        last_time_ps = float(
            frame_matches[-1][1]
        )

    timestep_match = re.search(
        r"Step\s+(\d+)\s+([-.+0-9eE]+)",
        text,
    )

    if timestep_match:
        timestep_ps = float(
            timestep_match.group(2)
        )

    return {
        "last_frame": last_frame,
        "first_time_ps": first_time_ps,
        "last_time_ps": last_time_ps,
        "reported_timestep_ps": (
            timestep_ps
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    log(
        "Day020 MD/confined-water "
        "input audit"
    )

    files = discover_files()

    inventory_rows = build_inventory(
        files
    )

    write_csv(
        INVENTORY_CSV,
        inventory_rows,
    )

    trajectories = [
        path
        for path in files
        if path.suffix.lower()
        in (
            ".xtc",
            ".trr",
            ".tng",
        )
    ]

    tpr_files = [
        path
        for path in files
        if path.suffix.lower()
        == ".tpr"
    ]

    gro_files = [
        path
        for path in files
        if path.suffix.lower()
        == ".gro"
    ]

    mdp_files = [
        path
        for path in files
        if path.suffix.lower()
        == ".mdp"
    ]

    ndx_files = [
        path
        for path in files
        if path.suffix.lower()
        == ".ndx"
    ]

    top_files = [
        path
        for path in files
        if path.suffix.lower()
        in (
            ".top",
            ".itp",
        )
    ]

    primary_trajectory = (
        select_primary_file(
            files,
            (
                ".xtc",
                ".trr",
                ".tng",
            ),
        )
    )

    primary_tpr = (
        select_primary_file(
            files,
            (".tpr",),
        )
    )

    primary_gro = (
        select_primary_file(
            files,
            (".gro",),
        )
    )

    gromacs = find_gromacs()

    gromacs_version = ""
    trajectory_check = {
        "last_frame": "",
        "first_time_ps": "",
        "last_time_ps": "",
        "reported_timestep_ps": "",
    }

    tpr_dump_text = ""

    if gromacs is not None:
        return_code, version_text = (
            run_command(
                [
                    gromacs,
                    "--version",
                ],
                OUTPUT_ROOT
                / "gromacs_version.txt",
            )
        )

        if return_code == 0:
            gromacs_version = (
                version_text.strip()
            )

        if primary_trajectory is not None:
            (
                trajectory_return_code,
                trajectory_output,
            ) = run_command(
                [
                    gromacs,
                    "check",
                    "-f",
                    str(
                        primary_trajectory
                    ),
                ],
                OUTPUT_ROOT
                / "gromacs_check_trajectory.txt",
            )

            trajectory_check = (
                parse_trajectory_check(
                    trajectory_output
                )
            )

            trajectory_check[
                "return_code"
            ] = trajectory_return_code

        if primary_tpr is not None:
            run_command(
                [
                    gromacs,
                    "check",
                    "-s1",
                    str(primary_tpr),
                ],
                OUTPUT_ROOT
                / "gromacs_check_tpr.txt",
            )

            (
                _dump_return_code,
                tpr_dump_text,
            ) = run_command(
                [
                    gromacs,
                    "dump",
                    "-s",
                    str(primary_tpr),
                ],
                OUTPUT_ROOT
                / "gromacs_dump_tpr.txt",
            )

    parameter_rows: list[
        dict[str, object]
    ] = []

    for mdp_file in mdp_files:
        parameter_rows.extend(
            parse_mdp_lines(
                mdp_file.read_text(
                    encoding="utf-8",
                    errors="replace",
                ),
                mdp_file,
            )
        )

    if tpr_dump_text:
        parameter_rows.extend(
            parse_tpr_dump_parameters(
                tpr_dump_text,
                (
                    "gmx dump: "
                    + relative(
                        primary_tpr
                    )
                ),
            )
        )

    write_csv(
        PARAMETER_CSV,
        parameter_rows,
    )

    structure_rows: list[
        dict[str, object]
    ] = []

    declared_atoms: int | str = ""
    total_residues: int | str = ""

    if primary_gro is not None:
        (
            structure_rows,
            declared_atoms,
            total_residues,
        ) = parse_gro(
            primary_gro
        )

    write_csv(
        STRUCTURE_CSV,
        structure_rows,
    )

    index_rows: list[
        dict[str, object]
    ] = []

    for ndx_file in ndx_files:
        index_rows.extend(
            parse_ndx(
                ndx_file
            )
        )

    write_csv(
        INDEX_CSV,
        index_rows,
    )

    freeze_matches = (
        search_freeze_settings(
            files
        )
    )

    freeze_detected = bool(
        freeze_matches
    )

    water_residue_count = sum(
        int(row["residue_count"])
        for row in structure_rows
        if row["category"] == "water"
    )

    water_atom_count = sum(
        int(row["atom_count"])
        for row in structure_rows
        if row["category"] == "water"
    )

    ion_count = sum(
        int(row["residue_count"])
        for row in structure_rows
        if row["category"] == "ion"
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 MD and Confined-Water "
            "Input Audit\n\n"
        )

        handle.write(
            "## Accepted run\n\n"
        )

        handle.write(
            f"- Run directory: "
            f"`{relative(RUN_ROOT)}`.\n"
        )

        handle.write(
            f"- Files discovered: "
            f"{len(files)}.\n"
        )

        handle.write(
            f"- Trajectory files: "
            f"{len(trajectories)}.\n"
        )

        handle.write(
            f"- TPR files: "
            f"{len(tpr_files)}.\n"
        )

        handle.write(
            f"- GRO files: "
            f"{len(gro_files)}.\n"
        )

        handle.write(
            f"- MDP files: "
            f"{len(mdp_files)}.\n"
        )

        handle.write(
            f"- NDX files: "
            f"{len(ndx_files)}.\n"
        )

        handle.write(
            f"- TOP/ITP files: "
            f"{len(top_files)}.\n\n"
        )

        handle.write(
            "## Primary files selected\n\n"
        )

        handle.write(
            f"- Trajectory: "
            f"`{relative(primary_trajectory) if primary_trajectory else 'not found'}`.\n"
        )

        handle.write(
            f"- Run input: "
            f"`{relative(primary_tpr) if primary_tpr else 'not found'}`.\n"
        )

        handle.write(
            f"- Structure: "
            f"`{relative(primary_gro) if primary_gro else 'not found'}`.\n\n"
        )

        handle.write(
            "## GROMACS inspection\n\n"
        )

        handle.write(
            f"- Executable: "
            f"`{gromacs or 'not found'}`.\n"
        )

        handle.write(
            f"- GROMACS available: "
            f"{gromacs is not None}.\n"
        )

        handle.write(
            f"- Last reported frame: "
            f"{trajectory_check.get('last_frame', '')}.\n"
        )

        handle.write(
            f"- First trajectory time: "
            f"{trajectory_check.get('first_time_ps', '')} ps.\n"
        )

        handle.write(
            f"- Last trajectory time: "
            f"{trajectory_check.get('last_time_ps', '')} ps.\n\n"
        )

        handle.write(
            "## Structure composition\n\n"
        )

        handle.write(
            f"- Declared atoms in selected GRO: "
            f"{declared_atoms}.\n"
        )

        handle.write(
            f"- Residues in selected GRO: "
            f"{total_residues}.\n"
        )

        handle.write(
            f"- Water residues: "
            f"{water_residue_count}.\n"
        )

        handle.write(
            f"- Water atoms: "
            f"{water_atom_count}.\n"
        )

        handle.write(
            f"- Ions: {ion_count}.\n\n"
        )

        handle.write(
            "## Frozen-solute evidence\n\n"
        )

        handle.write(
            f"- Freeze-related settings found: "
            f"{freeze_detected}.\n"
        )

        for source, line in freeze_matches:
            handle.write(
                f"- `{source}`: `{line}`\n"
            )

        if not freeze_matches:
            handle.write(
                "- No textual freeze setting was "
                "found. Inspect the TPR dump and run "
                "log before assuming that the solute "
                "was mobile.\n"
            )

        handle.write(
            "\n## Scientifically valid analyses "
            "for the accepted trajectory\n\n"
        )

        handle.write(
            "- Confined-water density and spatial "
            "heterogeneity.\n"
        )

        handle.write(
            "- Water orientation relative to the "
            "tube axis and local chromophore geometry.\n"
        )

        handle.write(
            "- Water–solute contacts and hydrogen-bond "
            "occupancy.\n"
        )

        handle.write(
            "- Snapshot-resolved electrostatic "
            "environment and disorder.\n"
        )

        handle.write(
            "- Short-time water correlation functions, "
            "provided their sampling limitations are "
            "reported explicitly.\n\n"
        )

        handle.write(
            "## Analyses not supported by the "
            "frozen-solute trajectory\n\n"
        )

        handle.write(
            "- Solute RMSD or RMSF as thermal-stability "
            "metrics.\n"
        )

        handle.write(
            "- Scaffold or chromophore conformational "
            "stability.\n"
        )

        handle.write(
            "- Coupled water–solute structural dynamics.\n"
        )

        handle.write(
            "- Converged long-time diffusion or residence "
            "times from the current 100 ps window.\n"
        )

        handle.write(
            "- A microscopic spectral density derived "
            "from coupled structural fluctuations.\n\n"
        )

        handle.write(
            "## Required next decision\n\n"
        )

        handle.write(
            "Use the audited groups and topology to "
            "define the confined-water structural "
            "analysis, then construct a controlled "
            "restraint-release sequence for a mobile-solute "
            "trajectory. The mobile trajectory must be "
            "validated before RMSD, RMSF, structural "
            "stability, or coupled bath dynamics are "
            "reported.\n"
        )

    log("")
    log(
        "Day020 MD/confined-water "
        "input audit completed."
    )

    log(
        f"Files discovered: {len(files)}"
    )

    log(
        f"Trajectories: {len(trajectories)}"
    )

    log(
        f"TPR files: {len(tpr_files)}"
    )

    log(
        f"GRO files: {len(gro_files)}"
    )

    log(
        f"MDP files: {len(mdp_files)}"
    )

    log(
        f"NDX files: {len(ndx_files)}"
    )

    log(
        f"TOP/ITP files: {len(top_files)}"
    )

    log(
        f"GROMACS executable: "
        f"{gromacs or 'not found'}"
    )

    log(
        f"Primary trajectory: "
        f"{relative(primary_trajectory) if primary_trajectory else 'not found'}"
    )

    log(
        f"Primary TPR: "
        f"{relative(primary_tpr) if primary_tpr else 'not found'}"
    )

    log(
        f"Primary GRO: "
        f"{relative(primary_gro) if primary_gro else 'not found'}"
    )

    log(
        f"Last trajectory frame: "
        f"{trajectory_check.get('last_frame', '')}"
    )

    log(
        f"Last trajectory time: "
        f"{trajectory_check.get('last_time_ps', '')} ps"
    )

    log(
        f"Declared atoms: "
        f"{declared_atoms}"
    )

    log(
        f"Water residues: "
        f"{water_residue_count}"
    )

    log(
        f"Ions: {ion_count}"
    )

    log(
        f"Freeze evidence detected: "
        f"{freeze_detected}"
    )

    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
