#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

AUDIT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_accepted_hydrated_topology_audit"
)

FILES = {
    "accepted_dump": (
        AUDIT_ROOT
        / "accepted_tpr_dump.txt"
    ),
    "rebuilt_dump": (
        AUDIT_ROOT
        / "rebuilt_tpr_dump.txt"
    ),
    "gmx_check": (
        AUDIT_ROOT
        / "accepted_vs_rebuilt_tpr_check.log"
    ),
    "processed_top": (
        AUDIT_ROOT
        / "candidate_processed.top"
    ),
    "mdout": (
        AUDIT_ROOT
        / "candidate_mdout.mdp"
    ),
}

OUTPUT_PATH = (
    AUDIT_ROOT
    / "day021_tpr_dump_diagnostic_full.txt"
)


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


def write_line(
    handle,
    text: str = "",
) -> None:
    handle.write(
        text + "\n"
    )


def main() -> None:
    missing = [
        path
        for path in FILES.values()
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty diagnostic inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    with OUTPUT_PATH.open(
        "w",
        encoding="utf-8",
    ) as output:
        write_line(
            output,
            "DAY021 TPR DUMP DIAGNOSTIC",
        )

        write_line(
            output,
            "=" * 88,
        )

        write_line(
            output,
            f"project_root: {PROJECT_ROOT}",
        )

        write_line(
            output,
            f"audit_root: {AUDIT_ROOT}",
        )

        for label, path in FILES.items():
            write_line(output)
            write_line(
                output,
                f"===== {label} =====",
            )

            write_line(
                output,
                f"path: {path}",
            )

            write_line(
                output,
                f"size_bytes: {path.stat().st_size}",
            )

            write_line(
                output,
                f"sha256: {sha256(path)}",
            )

            text = path.read_text(
                encoding="utf-8",
                errors="replace",
            )

            lines = text.splitlines()

            write_line(
                output,
                f"line_count: {len(lines)}",
            )

            if label.endswith("_dump"):
                write_line(output)
                write_line(
                    output,
                    "Top-level candidate headers:",
                )

                header_rows = []

                for line_number, line in enumerate(
                    lines,
                    start=1,
                ):
                    if not line:
                        continue

                    if line[0].isspace():
                        continue

                    if (
                        line.endswith(":")
                        or re.match(
                            (
                                r"^[A-Za-z_]"
                                r"[A-Za-z0-9_ ]*"
                                r"\s*\([^)]*\):$"
                            ),
                            line,
                        )
                    ):
                        header_rows.append(
                            (
                                line_number,
                                line,
                            )
                        )

                for (
                    line_number,
                    line,
                ) in header_rows:
                    write_line(
                        output,
                        f"{line_number:8d}: {line}",
                    )

                write_line(output)
                write_line(
                    output,
                    "Key section matches:",
                )

                patterns = (
                    "inputrec",
                    "topology",
                    "moltype",
                    "atoms",
                    "idef",
                    "x (",
                    "v (",
                    "box (",
                )

                for pattern in patterns:
                    matches = [
                        (
                            line_number,
                            line.strip(),
                        )
                        for line_number, line
                        in enumerate(
                            lines,
                            start=1,
                        )
                        if pattern.lower()
                        in line.lower()
                    ]

                    write_line(output)

                    write_line(
                        output,
                        (
                            f"[{pattern}] "
                            f"matches={len(matches)}"
                        ),
                    )

                    for (
                        line_number,
                        line,
                    ) in matches:
                        write_line(
                            output,
                            (
                                f"{line_number:8d}: "
                                f"{line}"
                            ),
                        )

                write_line(output)
                write_line(
                    output,
                    "Complete TPR dump:",
                )

                write_line(
                    output,
                    "-" * 88,
                )

                for line_number, line in enumerate(
                    lines,
                    start=1,
                ):
                    write_line(
                        output,
                        f"{line_number:8d}: {line}",
                    )

            elif label == "gmx_check":
                write_line(output)
                write_line(
                    output,
                    "Complete gmx check output:",
                )

                write_line(
                    output,
                    "-" * 88,
                )

                for line_number, line in enumerate(
                    lines,
                    start=1,
                ):
                    write_line(
                        output,
                        f"{line_number:8d}: {line}",
                    )

            else:
                write_line(output)
                write_line(
                    output,
                    "Complete file contents:",
                )

                write_line(
                    output,
                    "-" * 88,
                )

                for line_number, line in enumerate(
                    lines,
                    start=1,
                ):
                    write_line(
                        output,
                        f"{line_number:8d}: {line}",
                    )

    if (
        not OUTPUT_PATH.exists()
        or OUTPUT_PATH.stat().st_size == 0
    ):
        raise RuntimeError(
            "Diagnostic output was not generated"
        )

    print(
        "Day021 TPR dump diagnostic export completed."
    )

    print(
        f"Output: "
        f"{OUTPUT_PATH.relative_to(PROJECT_ROOT)}"
    )

    print(
        f"Size: "
        f"{OUTPUT_PATH.stat().st_size} bytes"
    )

    print(
        f"SHA256: "
        f"{sha256(OUTPUT_PATH)}"
    )


if __name__ == "__main__":
    main()
