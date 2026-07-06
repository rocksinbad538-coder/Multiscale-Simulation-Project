#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
import re
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

MOBILE_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_"
    "improperGeo100_validated"
)

BASELINE_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_"
    "improperGeo100_hydrated_baseline"
)

FROZEN_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032"
)

ACCEPTED_RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

FINAL_GRO = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.gro"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

FILE_INVENTORY_CSV = (
    OUTPUT_ROOT
    / "mobile_protocol_input_inventory.csv"
)

RESTRAINTS_CSV = (
    OUTPUT_ROOT
    / "existing_position_restraints.csv"
)

MDP_SETTINGS_CSV = (
    OUTPUT_ROOT
    / "existing_mdp_relevant_settings.csv"
)

INDEX_GROUPS_CSV = (
    OUTPUT_ROOT
    / "existing_index_groups.csv"
)

VELOCITY_CSV = (
    OUTPUT_ROOT
    / "final_gro_velocity_audit.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "MOBILE_RESTRAINT_PROTOCOL_INPUT_AUDIT_DAY021.md"
)

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)

INCLUDE_PATTERN = re.compile(
    r'^\s*#include\s+[<"]([^>"]+)[>"]'
)

IFDEF_PATTERN = re.compile(
    r"^\s*#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)$"
)

RELEVANT_SUFFIXES = {
    ".top",
    ".itp",
    ".mdp",
    ".ndx",
    ".gro",
    ".tpr",
    ".cpt",
}

RELEVANT_MDP_KEYS = (
    "integrator",
    "dt",
    "nsteps",
    "continuation",
    "gen-vel",
    "gen-temp",
    "gen-seed",
    "ld-seed",
    "define",
    "freezegrps",
    "freezedim",
    "constraints",
    "constraint-algorithm",
    "tcoupl",
    "tc-grps",
    "tau-t",
    "ref-t",
    "pcoupl",
    "ref-p",
    "tau-p",
    "nstxout-compressed",
    "nstenergy",
    "nstlog",
    "emtol",
    "emstep",
)


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


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    empty_row: dict[str, object],
) -> None:
    if not rows:
        rows = [empty_row]

    fields: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)

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
                    field: row.get(field, "")
                    for field in fields
                }
            )


def candidate_roots() -> list[Path]:
    roots = [
        MOBILE_ROOT,
        BASELINE_ROOT,
        FROZEN_ROOT,
        ACCEPTED_RUN_ROOT,
    ]

    return [
        path
        for path in roots
        if path.exists()
    ]


def inventory_files() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()

    for root in candidate_roots():
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue

            if path.suffix.lower() not in RELEVANT_SUFFIXES:
                continue

            resolved = path.resolve()

            if resolved in seen:
                continue

            seen.add(resolved)

            rows.append(
                {
                    "root": relative(root),
                    "path": relative(path),
                    "suffix": path.suffix.lower(),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )

    return rows


def parse_topology_directives(
    path: Path,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
]:
    restraints: list[dict[str, object]] = []
    directives: list[dict[str, object]] = []

    current_section = ""
    restraint_entries = 0
    restraint_start_line: int | None = None

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    def close_restraint_section() -> None:
        nonlocal restraint_entries
        nonlocal restraint_start_line

        if (
            current_section == "position_restraints"
            and restraint_start_line is not None
        ):
            restraints.append(
                {
                    "path": relative(path),
                    "section_start_line": restraint_start_line,
                    "entry_count": restraint_entries,
                }
            )

        restraint_entries = 0
        restraint_start_line = None

    for line_number, raw_line in enumerate(
        lines,
        start=1,
    ):
        include_match = INCLUDE_PATTERN.match(raw_line)

        if include_match:
            directives.append(
                {
                    "path": relative(path),
                    "line_number": line_number,
                    "directive": "include",
                    "value": include_match.group(1),
                }
            )

        conditional_match = IFDEF_PATTERN.match(raw_line)

        if conditional_match:
            directives.append(
                {
                    "path": relative(path),
                    "line_number": line_number,
                    "directive": conditional_match.group(1),
                    "value": conditional_match.group(2).strip(),
                }
            )

        line = strip_comment(raw_line)

        if not line or line.startswith("#"):
            continue

        section_match = SECTION_PATTERN.match(line)

        if section_match:
            close_restraint_section()

            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            if current_section == "position_restraints":
                restraint_start_line = line_number

            continue

        if current_section == "position_restraints":
            fields = line.split()

            try:
                int(fields[0])
            except (IndexError, ValueError):
                continue

            restraint_entries += 1

    close_restraint_section()

    return restraints, directives


def parse_mdp(path: Path) -> dict[str, str]:
    settings: dict[str, str] = {}

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = strip_comment(raw_line)

        if not line or "=" not in line:
            continue

        key, value = line.split("=", 1)

        normalized_key = (
            key.strip()
            .lower()
            .replace("_", "-")
        )

        settings[normalized_key] = value.strip()

    return settings


def mdp_rows(
    paths: list[Path],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for path in sorted(paths):
        settings = parse_mdp(path)

        row: dict[str, object] = {
            "path": relative(path),
        }

        for key in RELEVANT_MDP_KEYS:
            row[key] = settings.get(key, "")

        rows.append(row)

    return rows


def parse_ndx(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    current_group: str | None = None
    atom_count = 0

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.strip()

        section_match = SECTION_PATTERN.match(line)

        if section_match:
            if current_group is not None:
                rows.append(
                    {
                        "path": relative(path),
                        "group": current_group,
                        "atom_count": atom_count,
                    }
                )

            current_group = section_match.group(1).strip()
            atom_count = 0
            continue

        if current_group is None or not line:
            continue

        atom_count += sum(
            1
            for token in line.split()
            if token.isdigit()
        )

    if current_group is not None:
        rows.append(
            {
                "path": relative(path),
                "group": current_group,
                "atom_count": atom_count,
            }
        )

    return rows


def parse_final_gro_velocities() -> list[dict[str, object]]:
    lines = FINAL_GRO.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    natoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]

    if len(atom_lines) != natoms:
        raise RuntimeError(
            "GRO atom count does not match file length"
        )

    counts: Counter[str] = Counter()
    nonzero: Counter[str] = Counter()
    speed_sum: Counter[str] = Counter()
    speed2_sum: Counter[str] = Counter()
    max_speed: dict[str, float] = {}

    for atom_line in atom_lines:
        residue_name = atom_line[5:10].strip()

        if len(atom_line) < 68:
            raise RuntimeError(
                "Final GRO does not contain complete velocities"
            )

        vx = float(atom_line[44:52])
        vy = float(atom_line[52:60])
        vz = float(atom_line[60:68])

        speed = math.sqrt(
            vx * vx + vy * vy + vz * vz
        )

        counts[residue_name] += 1
        speed_sum[residue_name] += speed
        speed2_sum[residue_name] += speed * speed

        if speed > 1.0e-10:
            nonzero[residue_name] += 1

        max_speed[residue_name] = max(
            max_speed.get(residue_name, 0.0),
            speed,
        )

    rows: list[dict[str, object]] = []

    for residue_name in sorted(counts):
        count = counts[residue_name]

        rows.append(
            {
                "residue_name": residue_name,
                "atom_count": count,
                "nonzero_velocity_atoms": nonzero[residue_name],
                "zero_velocity_atoms": (
                    count - nonzero[residue_name]
                ),
                "nonzero_fraction": (
                    nonzero[residue_name] / count
                ),
                "mean_speed_nm_per_ps": (
                    speed_sum[residue_name] / count
                ),
                "rms_speed_nm_per_ps": math.sqrt(
                    speed2_sum[residue_name] / count
                ),
                "max_speed_nm_per_ps": (
                    max_speed[residue_name]
                ),
            }
        )

    return rows


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not FINAL_GRO.exists():
        raise RuntimeError(
            f"Missing final GRO: {FINAL_GRO}"
        )

    file_rows = inventory_files()

    topology_paths = [
        PROJECT_ROOT / row["path"]
        for row in file_rows
        if row["suffix"] in {".top", ".itp"}
    ]

    mdp_paths = [
        PROJECT_ROOT / row["path"]
        for row in file_rows
        if row["suffix"] == ".mdp"
    ]

    ndx_paths = [
        PROJECT_ROOT / row["path"]
        for row in file_rows
        if row["suffix"] == ".ndx"
    ]

    restraint_rows: list[dict[str, object]] = []
    directive_rows: list[dict[str, object]] = []

    for path in topology_paths:
        restraints, directives = (
            parse_topology_directives(path)
        )

        restraint_rows.extend(restraints)
        directive_rows.extend(directives)

    mdp_setting_rows = mdp_rows(mdp_paths)

    index_rows: list[dict[str, object]] = []

    for path in ndx_paths:
        index_rows.extend(parse_ndx(path))

    velocity_rows = parse_final_gro_velocities()

    write_csv(
        FILE_INVENTORY_CSV,
        file_rows,
        {"status": "NO_FILES_FOUND"},
    )

    write_csv(
        RESTRAINTS_CSV,
        restraint_rows,
        {"status": "NO_POSITION_RESTRAINTS_FOUND"},
    )

    write_csv(
        MDP_SETTINGS_CSV,
        mdp_setting_rows,
        {"status": "NO_MDP_FILES_FOUND"},
    )

    write_csv(
        INDEX_GROUPS_CSV,
        index_rows,
        {"status": "NO_INDEX_FILES_FOUND"},
    )

    write_csv(
        VELOCITY_CSV,
        velocity_rows,
        {"status": "NO_VELOCITY_DATA"},
    )

    hbn_velocity = next(
        (
            row
            for row in velocity_rows
            if row["residue_name"] == "HBN"
        ),
        None,
    )

    pyr_velocity = next(
        (
            row
            for row in velocity_rows
            if row["residue_name"] == "PYR"
        ),
        None,
    )

    sol_velocity = next(
        (
            row
            for row in velocity_rows
            if row["residue_name"] == "SOL"
        ),
        None,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Mobile Restraint-Protocol Input Audit\n\n"
        )

        handle.write(
            f"- Relevant files: {len(file_rows)}\n"
        )

        handle.write(
            f"- TOP/ITP files: {len(topology_paths)}\n"
        )

        handle.write(
            f"- Existing position-restraint sections: "
            f"{len(restraint_rows)}\n"
        )

        handle.write(
            f"- Existing MDP files: {len(mdp_paths)}\n"
        )

        handle.write(
            f"- Existing NDX files: {len(ndx_paths)}\n"
        )

        handle.write(
            f"- Conditional/include directives: "
            f"{len(directive_rows)}\n\n"
        )

        handle.write(
            "## Final GRO velocities\n\n"
        )

        for row in velocity_rows:
            handle.write(
                f"- {row['residue_name']}: "
                f"{row['nonzero_velocity_atoms']}/"
                f"{row['atom_count']} atoms with nonzero velocity; "
                f"RMS speed "
                f"{row['rms_speed_nm_per_ps']:.8f} nm/ps.\n"
            )

        handle.write(
            "\n## Protocol implication\n\n"
        )

        if (
            hbn_velocity is not None
            and pyr_velocity is not None
            and sol_velocity is not None
        ):
            hbn_zero = (
                int(hbn_velocity["nonzero_velocity_atoms"]) == 0
            )

            pyr_zero = (
                int(pyr_velocity["nonzero_velocity_atoms"]) == 0
            )

            sol_mobile = (
                int(sol_velocity["nonzero_velocity_atoms"]) > 0
            )

            if hbn_zero and pyr_zero and sol_mobile:
                handle.write(
                    "The final frozen-solute GRO contains zero "
                    "velocities for HBN and PYR but propagated "
                    "water velocities. Direct continuation into a "
                    "mobile-solute trajectory would therefore create "
                    "a non-equilibrated kinetic partition. A staged "
                    "protocol must explicitly regenerate or "
                    "re-equilibrate velocities before full release.\n"
                )
            else:
                handle.write(
                    "Velocity handling requires manual review because "
                    "the final GRO does not show the expected frozen-"
                    "solute/nonzero-water partition.\n"
                )

    print(
        "Day021 mobile restraint-protocol input audit completed."
    )

    print(
        f"Relevant files: {len(file_rows)}"
    )

    print(
        f"Existing position-restraint sections: "
        f"{len(restraint_rows)}"
    )

    print(
        f"Existing MDP files: {len(mdp_paths)}"
    )

    print(
        f"Existing NDX files: {len(ndx_paths)}"
    )

    print(
        "Final GRO velocity summary:"
    )

    for row in velocity_rows:
        print(
            f"  {row['residue_name']}: "
            f"{row['nonzero_velocity_atoms']}/"
            f"{row['atom_count']} nonzero; "
            f"RMS={row['rms_speed_nm_per_ps']:.8f} nm/ps"
        )

    print(
        f"Wrote: {relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
