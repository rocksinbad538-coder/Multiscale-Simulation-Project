#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import itertools
import re
from collections import OrderedDict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

AUDIT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit"
)

ACCEPTED_DUMP = (
    AUDIT_ROOT
    / "accepted_tpr_dump.txt"
)

REBUILT_DUMP = (
    AUDIT_ROOT
    / "rebuilt_tpr_dump.txt"
)

GMX_CHECK_LOG = (
    AUDIT_ROOT
    / "accepted_vs_rebuilt_tpr_check.log"
)

OUTPUT_PATH = (
    AUDIT_ROOT
    / "day021_tpr_dump_compact_comparison.txt"
)

MAX_DIFFERENCES_PER_SECTION = 80
MAX_LINE_LENGTH = 500

SECTION_HEADER_PATTERN = re.compile(
    r"^[^\s].*:\s*$"
)

ASSIGNMENT_PATTERN = re.compile(
    r"^\s*([^=]+?)\s*=\s*(.*?)\s*$"
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_text(
    lines: list[str],
) -> str:
    return "\n".join(
        line.rstrip()
        for line in lines
    ) + "\n"


def section_hash(
    lines: list[str],
) -> str:
    return sha256_bytes(
        canonical_text(lines).encode(
            "utf-8"
        )
    )


def truncate_line(
    line: str,
) -> str:
    if len(line) <= MAX_LINE_LENGTH:
        return line

    return (
        line[:MAX_LINE_LENGTH]
        + " ... [truncated]"
    )


def parse_sections(
    path: Path,
) -> OrderedDict[str, list[str]]:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    lines = text.splitlines()

    sections: OrderedDict[
        str,
        list[str]
    ] = OrderedDict()

    current_name = "__preamble__"
    sections[current_name] = []

    duplicate_counts: dict[str, int] = {}

    for line in lines:
        if SECTION_HEADER_PATTERN.match(
            line
        ):
            base_name = line.rstrip(
                ":"
            ).strip()

            duplicate_counts[
                base_name
            ] = (
                duplicate_counts.get(
                    base_name,
                    0,
                )
                + 1
            )

            occurrence = (
                duplicate_counts[
                    base_name
                ]
            )

            current_name = (
                base_name
                if occurrence == 1
                else (
                    f"{base_name}"
                    f" [occurrence {occurrence}]"
                )
            )

            sections[current_name] = [
                line
            ]

        else:
            sections[
                current_name
            ].append(line)

    return sections


def normalized_section_name(
    name: str,
) -> str:
    return re.sub(
        r"\s+\[occurrence\s+\d+\]$",
        "",
        name,
    ).strip().lower()


def is_dynamic_section(
    name: str,
) -> bool:
    normalized = (
        normalized_section_name(name)
    )

    dynamic_prefixes = (
        "x (",
        "v (",
        "box (",
        "box_rel",
        "boxv",
        "vir_prev",
        "pres_prev",
        "fvir",
        "svir",
    )

    return normalized.startswith(
        dynamic_prefixes
    )


def find_section(
    sections: OrderedDict[
        str,
        list[str],
    ],
    prefix: str,
) -> tuple[str, list[str]] | None:
    normalized_prefix = prefix.lower()

    for name, lines in sections.items():
        if normalized_section_name(
            name
        ).startswith(
            normalized_prefix
        ):
            return name, lines

    return None


def parse_assignments(
    lines: list[str],
) -> OrderedDict[str, str]:
    assignments: OrderedDict[
        str,
        str
    ] = OrderedDict()

    duplicate_counts: dict[str, int] = {}

    for line in lines:
        match = ASSIGNMENT_PATTERN.match(
            line
        )

        if match is None:
            continue

        key = " ".join(
            match.group(1).split()
        )

        value = " ".join(
            match.group(2).split()
        )

        duplicate_counts[key] = (
            duplicate_counts.get(
                key,
                0,
            )
            + 1
        )

        occurrence = duplicate_counts[
            key
        ]

        stored_key = (
            key
            if occurrence == 1
            else (
                f"{key}"
                f" [occurrence {occurrence}]"
            )
        )

        assignments[
            stored_key
        ] = value

    return assignments


def first_line_differences(
    accepted_lines: list[str],
    rebuilt_lines: list[str],
    limit: int,
) -> list[
    tuple[int, str, str]
]:
    differences: list[
        tuple[int, str, str]
    ] = []

    for line_number, (
        accepted_line,
        rebuilt_line,
    ) in enumerate(
        itertools.zip_longest(
            accepted_lines,
            rebuilt_lines,
            fillvalue="<LINE MISSING>",
        ),
        start=1,
    ):
        accepted_normalized = (
            accepted_line.rstrip()
        )

        rebuilt_normalized = (
            rebuilt_line.rstrip()
        )

        if (
            accepted_normalized
            == rebuilt_normalized
        ):
            continue

        differences.append(
            (
                line_number,
                truncate_line(
                    accepted_normalized
                ),
                truncate_line(
                    rebuilt_normalized
                ),
            )
        )

        if len(differences) >= limit:
            break

    return differences


def write_line(
    handle,
    text: str = "",
) -> None:
    handle.write(
        text + "\n"
    )


def main() -> None:
    required = (
        ACCEPTED_DUMP,
        REBUILT_DUMP,
        GMX_CHECK_LOG,
    )

    missing = [
        path
        for path in required
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    accepted_sections = parse_sections(
        ACCEPTED_DUMP
    )

    rebuilt_sections = parse_sections(
        REBUILT_DUMP
    )

    accepted_names = list(
        accepted_sections.keys()
    )

    rebuilt_names = list(
        rebuilt_sections.keys()
    )

    ordered_names: list[str] = []
    seen_names: set[str] = set()

    for name in (
        accepted_names
        + rebuilt_names
    ):
        if name in seen_names:
            continue

        seen_names.add(name)
        ordered_names.append(name)

    gmx_check_text = (
        GMX_CHECK_LOG.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    with OUTPUT_PATH.open(
        "w",
        encoding="utf-8",
    ) as output:
        write_line(
            output,
            "DAY021 COMPACT TPR DUMP COMPARISON",
        )

        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            f"accepted_dump: "
            f"{ACCEPTED_DUMP.relative_to(PROJECT_ROOT)}",
        )

        write_line(
            output,
            f"rebuilt_dump: "
            f"{REBUILT_DUMP.relative_to(PROJECT_ROOT)}",
        )

        write_line(
            output,
            f"accepted_dump_size_bytes: "
            f"{ACCEPTED_DUMP.stat().st_size}",
        )

        write_line(
            output,
            f"rebuilt_dump_size_bytes: "
            f"{REBUILT_DUMP.stat().st_size}",
        )

        write_line(
            output,
            f"accepted_section_count: "
            f"{len(accepted_sections)}",
        )

        write_line(
            output,
            f"rebuilt_section_count: "
            f"{len(rebuilt_sections)}",
        )

        write_line(output)
        write_line(
            output,
            "SECTION HASH AND SIZE COMPARISON",
        )

        write_line(
            output,
            "-" * 96,
        )

        unequal_static_sections: list[
            str
        ] = []

        unequal_dynamic_sections: list[
            str
        ] = []

        missing_sections: list[str] = []

        for name in ordered_names:
            accepted_lines = (
                accepted_sections.get(name)
            )

            rebuilt_lines = (
                rebuilt_sections.get(name)
            )

            write_line(output)
            write_line(
                output,
                f"[SECTION] {name}",
            )

            if accepted_lines is None:
                write_line(
                    output,
                    "accepted_present: False",
                )

                write_line(
                    output,
                    "rebuilt_present: True",
                )

                missing_sections.append(
                    name
                )

                continue

            if rebuilt_lines is None:
                write_line(
                    output,
                    "accepted_present: True",
                )

                write_line(
                    output,
                    "rebuilt_present: False",
                )

                missing_sections.append(
                    name
                )

                continue

            accepted_text = canonical_text(
                accepted_lines
            )

            rebuilt_text = canonical_text(
                rebuilt_lines
            )

            accepted_hash = sha256_bytes(
                accepted_text.encode(
                    "utf-8"
                )
            )

            rebuilt_hash = sha256_bytes(
                rebuilt_text.encode(
                    "utf-8"
                )
            )

            equal = (
                accepted_hash
                == rebuilt_hash
            )

            dynamic = is_dynamic_section(
                name
            )

            write_line(
                output,
                "accepted_present: True",
            )

            write_line(
                output,
                "rebuilt_present: True",
            )

            write_line(
                output,
                f"dynamic_state_section: {dynamic}",
            )

            write_line(
                output,
                f"accepted_lines: "
                f"{len(accepted_lines)}",
            )

            write_line(
                output,
                f"rebuilt_lines: "
                f"{len(rebuilt_lines)}",
            )

            write_line(
                output,
                f"accepted_bytes: "
                f"{len(accepted_text.encode('utf-8'))}",
            )

            write_line(
                output,
                f"rebuilt_bytes: "
                f"{len(rebuilt_text.encode('utf-8'))}",
            )

            write_line(
                output,
                f"accepted_sha256: "
                f"{accepted_hash}",
            )

            write_line(
                output,
                f"rebuilt_sha256: "
                f"{rebuilt_hash}",
            )

            write_line(
                output,
                f"exactly_equal: {equal}",
            )

            if not equal:
                if dynamic:
                    unequal_dynamic_sections.append(
                        name
                    )
                else:
                    unequal_static_sections.append(
                        name
                    )

        write_line(output)
        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            "INPUTREC PARAMETER COMPARISON",
        )

        write_line(
            output,
            "-" * 96,
        )

        accepted_inputrec = find_section(
            accepted_sections,
            "inputrec",
        )

        rebuilt_inputrec = find_section(
            rebuilt_sections,
            "inputrec",
        )

        if (
            accepted_inputrec is None
            or rebuilt_inputrec is None
        ):
            write_line(
                output,
                "inputrec section not found "
                "in one or both dumps.",
            )

        else:
            accepted_parameters = (
                parse_assignments(
                    accepted_inputrec[1]
                )
            )

            rebuilt_parameters = (
                parse_assignments(
                    rebuilt_inputrec[1]
                )
            )

            parameter_keys = list(
                dict.fromkeys(
                    list(
                        accepted_parameters.keys()
                    )
                    + list(
                        rebuilt_parameters.keys()
                    )
                )
            )

            parameter_differences = []

            for key in parameter_keys:
                accepted_value = (
                    accepted_parameters.get(
                        key,
                        "<MISSING>",
                    )
                )

                rebuilt_value = (
                    rebuilt_parameters.get(
                        key,
                        "<MISSING>",
                    )
                )

                if (
                    accepted_value
                    != rebuilt_value
                ):
                    parameter_differences.append(
                        (
                            key,
                            accepted_value,
                            rebuilt_value,
                        )
                    )

            write_line(
                output,
                f"accepted_parameter_count: "
                f"{len(accepted_parameters)}",
            )

            write_line(
                output,
                f"rebuilt_parameter_count: "
                f"{len(rebuilt_parameters)}",
            )

            write_line(
                output,
                f"differing_parameter_count: "
                f"{len(parameter_differences)}",
            )

            for (
                key,
                accepted_value,
                rebuilt_value,
            ) in parameter_differences:
                write_line(output)
                write_line(
                    output,
                    f"parameter: {key}",
                )

                write_line(
                    output,
                    f"accepted: "
                    f"{truncate_line(accepted_value)}",
                )

                write_line(
                    output,
                    f"rebuilt: "
                    f"{truncate_line(rebuilt_value)}",
                )

        write_line(output)
        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            "FIRST DIFFERENCES IN STATIC SECTIONS",
        )

        write_line(
            output,
            "-" * 96,
        )

        if not unequal_static_sections:
            write_line(
                output,
                "No unequal static sections detected.",
            )

        for name in (
            unequal_static_sections
        ):
            accepted_lines = (
                accepted_sections[
                    name
                ]
            )

            rebuilt_lines = (
                rebuilt_sections[
                    name
                ]
            )

            differences = (
                first_line_differences(
                    accepted_lines,
                    rebuilt_lines,
                    MAX_DIFFERENCES_PER_SECTION,
                )
            )

            write_line(output)
            write_line(
                output,
                f"[STATIC SECTION] {name}",
            )

            write_line(
                output,
                f"reported_differences: "
                f"{len(differences)}",
            )

            for (
                line_number,
                accepted_line,
                rebuilt_line,
            ) in differences:
                write_line(output)
                write_line(
                    output,
                    f"relative_line: "
                    f"{line_number}",
                )

                write_line(
                    output,
                    f"accepted: "
                    f"{accepted_line}",
                )

                write_line(
                    output,
                    f"rebuilt: "
                    f"{rebuilt_line}",
                )

        write_line(output)
        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            "DYNAMIC STATE SECTION SUMMARY",
        )

        write_line(
            output,
            "-" * 96,
        )

        if not unequal_dynamic_sections:
            write_line(
                output,
                "No unequal dynamic-state "
                "sections detected.",
            )

        for name in (
            unequal_dynamic_sections
        ):
            accepted_lines = (
                accepted_sections[
                    name
                ]
            )

            rebuilt_lines = (
                rebuilt_sections[
                    name
                ]
            )

            differences = (
                first_line_differences(
                    accepted_lines,
                    rebuilt_lines,
                    5,
                )
            )

            write_line(output)
            write_line(
                output,
                f"[DYNAMIC SECTION] {name}",
            )

            write_line(
                output,
                f"accepted_sha256: "
                f"{section_hash(accepted_lines)}",
            )

            write_line(
                output,
                f"rebuilt_sha256: "
                f"{section_hash(rebuilt_lines)}",
            )

            write_line(
                output,
                "Only the first five differences "
                "are shown because these sections "
                "contain per-atom state arrays.",
            )

            for (
                line_number,
                accepted_line,
                rebuilt_line,
            ) in differences:
                write_line(output)
                write_line(
                    output,
                    f"relative_line: "
                    f"{line_number}",
                )

                write_line(
                    output,
                    f"accepted: "
                    f"{accepted_line}",
                )

                write_line(
                    output,
                    f"rebuilt: "
                    f"{rebuilt_line}",
                )

        write_line(output)
        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            "MISSING SECTION SUMMARY",
        )

        write_line(
            output,
            "-" * 96,
        )

        write_line(
            output,
            f"missing_section_count: "
            f"{len(missing_sections)}",
        )

        for name in missing_sections:
            write_line(
                output,
                f"- {name}",
            )

        write_line(output)
        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            "COMPLETE GMX CHECK OUTPUT",
        )

        write_line(
            output,
            "-" * 96,
        )

        output.write(
            gmx_check_text.rstrip()
            + "\n"
        )

        write_line(output)
        write_line(
            output,
            "=" * 96,
        )

        write_line(
            output,
            "COMPACT COMPARISON SUMMARY",
        )

        write_line(
            output,
            "-" * 96,
        )

        write_line(
            output,
            f"unequal_static_section_count: "
            f"{len(unequal_static_sections)}",
        )

        write_line(
            output,
            f"unequal_dynamic_section_count: "
            f"{len(unequal_dynamic_sections)}",
        )

        write_line(
            output,
            f"missing_section_count: "
            f"{len(missing_sections)}",
        )

        write_line(
            output,
            "This file intentionally excludes "
            "complete coordinate and velocity arrays.",
        )

    if (
        not OUTPUT_PATH.exists()
        or OUTPUT_PATH.stat().st_size == 0
    ):
        raise RuntimeError(
            "Compact comparison was not generated"
        )

    print(
        "Day021 compact TPR comparison completed."
    )

    print(
        "Output: "
        f"{OUTPUT_PATH.relative_to(PROJECT_ROOT)}"
    )

    print(
        "Size: "
        f"{OUTPUT_PATH.stat().st_size} bytes"
    )

    print(
        "SHA256: "
        f"{sha256_bytes(OUTPUT_PATH.read_bytes())}"
    )


if __name__ == "__main__":
    main()
