#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import re
from collections import OrderedDict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

MOBILE_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_"
    "improperGeo100_validated"
)

MOBILE_TOP = (
    MOBILE_ROOT
    / "hbn_pyrene_4_hydratable_gap45_"
    "pyr5shift_clean032_hbnBonded_"
    "kang2000_improperGeo100.top"
)

AUDIT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_mobile_restraint_protocol"
)

RESTRAINTS_CSV = (
    AUDIT_ROOT
    / "existing_position_restraints.csv"
)

MDP_SETTINGS_CSV = (
    AUDIT_ROOT
    / "existing_mdp_relevant_settings.csv"
)

OUTPUT_PATH = (
    AUDIT_ROOT
    / "active_restraint_controls_compact.txt"
)

INCLUDE_PATTERN = re.compile(
    r'^\s*#include\s+[<"]([^>"]+)[>"]'
)

PREPROCESSOR_PATTERN = re.compile(
    r"^\s*#\s*"
    r"(ifdef|ifndef|if|elif|else|endif|define|undef)"
    r"\b(.*)$"
)

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)


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


def strip_comment(line: str) -> str:
    return line.split(
        ";",
        1,
    )[0].strip()


def resolve_include(
    parent: Path,
    include_token: str,
) -> Path | None:
    candidates = (
        parent.parent / include_token,
        MOBILE_ROOT / include_token,
        PROJECT_ROOT / include_token,
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    return None


def collect_include_closure(
    root_topology: Path,
) -> tuple[list[Path], list[str]]:
    visited: set[Path] = set()
    ordered: list[Path] = []
    unresolved: list[str] = []

    def recurse(path: Path) -> None:
        resolved = path.resolve()

        if resolved in visited:
            return

        visited.add(resolved)
        ordered.append(resolved)

        for raw_line in resolved.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines():
            match = INCLUDE_PATTERN.match(
                raw_line
            )

            if match is None:
                continue

            token = match.group(1)

            included = resolve_include(
                resolved,
                token,
            )

            if included is None:
                unresolved.append(
                    f"{relative(resolved)}::{token}"
                )

                continue

            recurse(included)

    recurse(root_topology)

    return ordered, unresolved


def parse_file(
    path: Path,
) -> dict[str, object]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    includes: list[
        tuple[int, str]
    ] = []

    preprocessors: list[
        tuple[int, str, str]
    ] = []

    restraint_sections: list[
        dict[str, object]
    ] = []

    current_section = ""
    current_restraint: (
        dict[str, object] | None
    ) = None

    for line_number, raw_line in enumerate(
        lines,
        start=1,
    ):
        include_match = INCLUDE_PATTERN.match(
            raw_line
        )

        if include_match:
            includes.append(
                (
                    line_number,
                    include_match.group(1),
                )
            )

        preprocessor_match = (
            PREPROCESSOR_PATTERN.match(
                raw_line
            )
        )

        if preprocessor_match:
            preprocessors.append(
                (
                    line_number,
                    preprocessor_match.group(1),
                    preprocessor_match.group(2).strip(),
                )
            )

        line = strip_comment(raw_line)

        if not line:
            continue

        if line.startswith("#"):
            continue

        section_match = SECTION_PATTERN.match(
            line
        )

        if section_match:
            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            if (
                current_section
                == "position_restraints"
            ):
                current_restraint = {
                    "start_line": line_number,
                    "entries": [],
                }

                restraint_sections.append(
                    current_restraint
                )
            else:
                current_restraint = None

            continue

        if (
            current_section
            != "position_restraints"
            or current_restraint is None
        ):
            continue

        fields = line.split()

        if len(fields) < 5:
            continue

        try:
            atom_index = int(fields[0])
            function_type = int(fields[1])
            force_x = float(fields[2])
            force_y = float(fields[3])
            force_z = float(fields[4])
        except ValueError:
            continue

        current_restraint[
            "entries"
        ].append(
            (
                atom_index,
                function_type,
                force_x,
                force_y,
                force_z,
            )
        )

    return {
        "includes": includes,
        "preprocessors": preprocessors,
        "restraint_sections": restraint_sections,
    }


def format_force(value: float) -> str:
    if value.is_integer():
        return str(int(value))

    return f"{value:.8g}"


def load_csv(
    path: Path,
) -> list[dict[str, str]]:
    if not path.exists():
        return []

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


def main() -> None:
    AUDIT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not MOBILE_TOP.exists():
        raise RuntimeError(
            f"Missing mobile topology: "
            f"{MOBILE_TOP}"
        )

    closure, unresolved = (
        collect_include_closure(
            MOBILE_TOP
        )
    )

    parsed = OrderedDict(
        (
            path,
            parse_file(path),
        )
        for path in closure
    )

    all_restraint_rows = load_csv(
        RESTRAINTS_CSV
    )

    mdp_rows = load_csv(
        MDP_SETTINGS_CSV
    )

    active_restraint_section_count = 0
    active_restraint_entry_count = 0

    with OUTPUT_PATH.open(
        "w",
        encoding="utf-8",
    ) as output:
        output.write(
            "DAY021 ACTIVE RESTRAINT CONTROLS\n"
        )

        output.write(
            "=" * 88 + "\n"
        )

        output.write(
            f"mobile_topology: "
            f"{relative(MOBILE_TOP)}\n"
        )

        output.write(
            f"include_closure_files: "
            f"{len(closure)}\n"
        )

        output.write(
            f"unresolved_includes: "
            f"{len(unresolved)}\n"
        )

        for item in unresolved:
            output.write(
                f"  unresolved: {item}\n"
            )

        output.write(
            "\nACTIVE INCLUDE CLOSURE\n"
        )

        output.write(
            "-" * 88 + "\n"
        )

        for path, information in parsed.items():
            output.write(
                f"\nFILE: {relative(path)}\n"
            )

            output.write(
                f"sha256: {sha256(path)}\n"
            )

            includes = information[
                "includes"
            ]

            preprocessors = information[
                "preprocessors"
            ]

            restraints = information[
                "restraint_sections"
            ]

            output.write(
                f"include_count: "
                f"{len(includes)}\n"
            )

            for line_number, token in includes:
                output.write(
                    f"  line {line_number}: "
                    f"#include \"{token}\"\n"
                )

            output.write(
                f"preprocessor_directive_count: "
                f"{len(preprocessors)}\n"
            )

            for (
                line_number,
                directive,
                value,
            ) in preprocessors:
                suffix = (
                    f" {value}"
                    if value
                    else ""
                )

                output.write(
                    f"  line {line_number}: "
                    f"#{directive}{suffix}\n"
                )

            output.write(
                f"position_restraint_section_count: "
                f"{len(restraints)}\n"
            )

            for section_index, section in enumerate(
                restraints,
                start=1,
            ):
                entries = section[
                    "entries"
                ]

                active_restraint_section_count += 1
                active_restraint_entry_count += len(
                    entries
                )

                force_triplets = sorted(
                    {
                        (
                            entry[2],
                            entry[3],
                            entry[4],
                        )
                        for entry in entries
                    }
                )

                atom_indices = [
                    entry[0]
                    for entry in entries
                ]

                output.write(
                    f"  section {section_index}:\n"
                )

                output.write(
                    f"    start_line: "
                    f"{section['start_line']}\n"
                )

                output.write(
                    f"    entry_count: "
                    f"{len(entries)}\n"
                )

                if atom_indices:
                    output.write(
                        f"    atom_index_range: "
                        f"{min(atom_indices)}-"
                        f"{max(atom_indices)}\n"
                    )

                output.write(
                    f"    unique_force_triplets: "
                    f"{len(force_triplets)}\n"
                )

                for force_triplet in (
                    force_triplets[:20]
                ):
                    output.write(
                        "      "
                        + " ".join(
                            format_force(value)
                            for value
                            in force_triplet
                        )
                        + "\n"
                    )

                if len(force_triplets) > 20:
                    output.write(
                        "      ... truncated ...\n"
                    )

                if entries:
                    first = entries[0]
                    last = entries[-1]

                    output.write(
                        "    first_entry: "
                        + " ".join(
                            str(value)
                            for value in first
                        )
                        + "\n"
                    )

                    output.write(
                        "    last_entry: "
                        + " ".join(
                            str(value)
                            for value in last
                        )
                        + "\n"
                    )

        output.write(
            "\nGLOBAL RESTRAINT INVENTORY\n"
        )

        output.write(
            "-" * 88 + "\n"
        )

        output.write(
            f"all_discovered_restraint_sections: "
            f"{len(all_restraint_rows)}\n"
        )

        for index, row in enumerate(
            all_restraint_rows,
            start=1,
        ):
            output.write(
                f"{index:02d}. "
                f"{row.get('path', '')} "
                f"line={row.get('section_start_line', '')} "
                f"entries={row.get('entry_count', '')}\n"
            )

        output.write(
            "\nEXISTING MDP SETTINGS\n"
        )

        output.write(
            "-" * 88 + "\n"
        )

        output.write(
            f"mdp_file_count: "
            f"{len(mdp_rows)}\n"
        )

        for row in mdp_rows:
            output.write(
                f"\nMDP: "
                f"{row.get('path', '')}\n"
            )

            for key, value in row.items():
                if key == "path":
                    continue

                if not value:
                    continue

                output.write(
                    f"  {key}: {value}\n"
                )

        output.write(
            "\nSUMMARY\n"
        )

        output.write(
            "-" * 88 + "\n"
        )

        output.write(
            f"active_restraint_sections: "
            f"{active_restraint_section_count}\n"
        )

        output.write(
            f"active_restraint_entries: "
            f"{active_restraint_entry_count}\n"
        )

        output.write(
            f"active_include_closure_complete: "
            f"{len(unresolved) == 0}\n"
        )

    print(
        "Day021 active restraint-control "
        "inspection completed."
    )

    print(
        f"Active include files: "
        f"{len(closure)}"
    )

    print(
        f"Unresolved includes: "
        f"{len(unresolved)}"
    )

    print(
        f"Active position-restraint sections: "
        f"{active_restraint_section_count}"
    )

    print(
        f"Active position-restraint entries: "
        f"{active_restraint_entry_count}"
    )

    print(
        f"All discovered restraint sections: "
        f"{len(all_restraint_rows)}"
    )

    print(
        f"Existing MDP files: "
        f"{len(mdp_rows)}"
    )

    print(
        f"Wrote: {relative(OUTPUT_PATH)}"
    )


if __name__ == "__main__":
    main()
