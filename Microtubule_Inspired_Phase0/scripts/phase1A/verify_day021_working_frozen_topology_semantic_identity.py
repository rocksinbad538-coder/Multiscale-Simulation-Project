#!/usr/bin/env python3

from __future__ import annotations

import csv
import difflib
import hashlib
import os
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

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

ACCEPTED_TOPOLOGY_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032"
)

ACCEPTED_TOP = (
    ACCEPTED_TOPOLOGY_ROOT
    / "hbn_pyrene_4_hydratable_gap45_"
    "pyr5shift_clean032.top"
)

WORKING_TOPOLOGY_ROOT = (
    PROJECT_ROOT
    / "parameters/phase1A/"
    "hybrid_hydrated_gromacs"
)

WORKING_TOP = (
    WORKING_TOPOLOGY_ROOT
    / "hbn_pyrene_4_hydratable_gap45_"
    "pyr5shift_clean032.top"
)

PRIOR_IDENTITY_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit/"
    "original_frozen_topology_identity"
)

AUDIT_MDP = (
    PRIOR_IDENTITY_ROOT
    / "nvt_100ps_frozenSolute_"
    "accepted_seed.mdp"
)

ACCEPTED_CLEAN_DUMP = (
    PRIOR_IDENTITY_ROOT
    / "accepted_tpr_clean_dump.txt"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit/"
    "working_frozen_topology_semantic_identity"
)

WORKING_REBUILT_TPR = (
    OUTPUT_ROOT
    / "working_frozen_topology_rebuilt.tpr"
)

WORKING_PROCESSED_TOP = (
    OUTPUT_ROOT
    / "working_frozen_topology_processed.top"
)

WORKING_MDOUT = (
    OUTPUT_ROOT
    / "working_frozen_topology_mdout.mdp"
)

WORKING_DUMP = (
    OUTPUT_ROOT
    / "working_frozen_topology_dump.txt"
)

WORKING_DUMP_STDERR = (
    OUTPUT_ROOT
    / "working_frozen_topology_dump_stderr.log"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "grompp_working_frozen_topology.log"
)

RAW_TOP_DIFF = (
    OUTPUT_ROOT
    / "accepted_vs_working_top_raw.diff"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "working_frozen_topology_semantic_identity.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "WORKING_FROZEN_TOPOLOGY_"
    "SEMANTIC_IDENTITY_DAY021.md"
)


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(
                PROJECT_ROOT
            )
        )
    except ValueError:
        return str(path.resolve())


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


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


def text_hash(lines: list[str]) -> str:
    normalized = (
        "\n".join(
            line.rstrip()
            for line in lines
        )
        + "\n"
    )

    return sha256_bytes(
        normalized.encode("utf-8")
    )


def canonical_topology_source(
    path: Path,
) -> list[str]:
    result: list[str] = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        without_comment = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not without_comment:
            continue

        result.append(
            " ".join(
                without_comment.split()
            )
        )

    return result


def extract_between(
    text: str,
    start_name: str,
    end_name: str,
) -> list[str]:
    lines = text.splitlines()

    start_index: int | None = None
    end_index: int | None = None

    for index, line in enumerate(lines):
        if line.strip() == start_name:
            start_index = index
            break

    if start_index is None:
        raise RuntimeError(
            f"Section start not found: {start_name}"
        )

    for index in range(
        start_index + 1,
        len(lines),
    ):
        if lines[index].strip() == end_name:
            end_index = index
            break

    if end_index is None:
        raise RuntimeError(
            f"Section end not found: {end_name}"
        )

    return lines[start_index:end_index]


def find_gmx() -> Path:
    preferred = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if preferred.exists():
        return preferred

    discovered = shutil.which("gmx")

    if discovered is None:
        raise RuntimeError(
            "Could not locate GROMACS gmx"
        )

    return Path(discovered)


def run_combined(
    command: list[str],
    cwd: Path,
    log_path: Path,
) -> subprocess.CompletedProcess[str]:
    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    result = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    log_path.write_text(
        result.stdout,
        encoding="utf-8",
    )

    return result


def run_dump(
    gmx: Path,
    tpr_path: Path,
) -> None:
    environment = {
        **os.environ,
        "GMX_MAXBACKUP": "-1",
    }

    result = subprocess.run(
        [
            str(gmx),
            "dump",
            "-s",
            str(tpr_path),
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    WORKING_DUMP.write_text(
        result.stdout,
        encoding="utf-8",
    )

    WORKING_DUMP_STDERR.write_text(
        result.stderr,
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RuntimeError(
            "gmx dump failed for the "
            "working-copy reconstruction"
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


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = (
        FINAL_GRO,
        ACCEPTED_TOP,
        WORKING_TOP,
        AUDIT_MDP,
        ACCEPTED_CLEAN_DUMP,
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
            "Missing or empty required inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    gmx = find_gmx()

    accepted_raw_hash = sha256_file(
        ACCEPTED_TOP
    )

    working_raw_hash = sha256_file(
        WORKING_TOP
    )

    raw_top_equal = (
        accepted_raw_hash
        == working_raw_hash
    )

    accepted_canonical = (
        canonical_topology_source(
            ACCEPTED_TOP
        )
    )

    working_canonical = (
        canonical_topology_source(
            WORKING_TOP
        )
    )

    canonical_top_equal = (
        accepted_canonical
        == working_canonical
    )

    raw_diff_lines = list(
        difflib.unified_diff(
            ACCEPTED_TOP.read_text(
                encoding="utf-8",
                errors="replace",
            ).splitlines(),
            WORKING_TOP.read_text(
                encoding="utf-8",
                errors="replace",
            ).splitlines(),
            fromfile=relative(ACCEPTED_TOP),
            tofile=relative(WORKING_TOP),
            lineterm="",
        )
    )

    RAW_TOP_DIFF.write_text(
        "\n".join(raw_diff_lines) + "\n",
        encoding="utf-8",
    )

    grompp_result = run_combined(
        [
            str(gmx),
            "grompp",
            "-f",
            str(AUDIT_MDP),
            "-c",
            str(FINAL_GRO),
            "-p",
            str(WORKING_TOP),
            "-o",
            str(WORKING_REBUILT_TPR),
            "-po",
            str(WORKING_MDOUT),
            "-pp",
            str(WORKING_PROCESSED_TOP),
            "-maxwarn",
            "0",
        ],
        WORKING_TOPOLOGY_ROOT,
        GROMPP_LOG,
    )

    grompp_pass = (
        grompp_result.returncode == 0
        and WORKING_REBUILT_TPR.exists()
        and WORKING_REBUILT_TPR.stat().st_size > 0
    )

    if not grompp_pass:
        raise RuntimeError(
            "Working-copy topology reconstruction "
            f"failed. See {GROMPP_LOG}"
        )

    run_dump(
        gmx,
        WORKING_REBUILT_TPR,
    )

    accepted_text = (
        ACCEPTED_CLEAN_DUMP.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    working_text = (
        WORKING_DUMP.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    accepted_inputrec = extract_between(
        accepted_text,
        "inputrec:",
        "qm-opts:",
    )

    working_inputrec = extract_between(
        working_text,
        "inputrec:",
        "qm-opts:",
    )

    accepted_topology = extract_between(
        accepted_text,
        "topology:",
        "box (3x3):",
    )

    working_topology = extract_between(
        working_text,
        "topology:",
        "box (3x3):",
    )

    accepted_box = extract_between(
        accepted_text,
        "box (3x3):",
        "box_rel (3x3):",
    )

    working_box = extract_between(
        working_text,
        "box (3x3):",
        "box_rel (3x3):",
    )

    inputrec_equal = (
        text_hash(accepted_inputrec)
        == text_hash(working_inputrec)
    )

    processed_topology_equal = (
        text_hash(accepted_topology)
        == text_hash(working_topology)
    )

    box_equal = (
        text_hash(accepted_box)
        == text_hash(working_box)
    )

    semantic_identity = all(
        (
            grompp_pass,
            inputrec_equal,
            processed_topology_equal,
            box_equal,
        )
    )

    summary = {
        "accepted_TOP": (
            relative(ACCEPTED_TOP)
        ),
        "working_TOP": (
            relative(WORKING_TOP)
        ),
        "accepted_TOP_sha256": (
            accepted_raw_hash
        ),
        "working_TOP_sha256": (
            working_raw_hash
        ),
        "raw_TOP_exactly_equal": (
            raw_top_equal
        ),
        "comment_whitespace_normalized_TOP_equal": (
            canonical_top_equal
        ),
        "raw_diff_line_count": (
            len(raw_diff_lines)
        ),
        "grompp_return_code": (
            grompp_result.returncode
        ),
        "grompp_pass": (
            grompp_pass
        ),
        "inputrec_exactly_equal_to_accepted_TPR": (
            inputrec_equal
        ),
        "processed_topology_exactly_equal_to_accepted_TPR": (
            processed_topology_equal
        ),
        "box_exactly_equal_to_accepted_TPR": (
            box_equal
        ),
        "working_copy_semantic_identity": (
            "PASS"
            if semantic_identity
            else "FAIL"
        ),
        "canonical_source_decision": (
            relative(ACCEPTED_TOP)
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
            "# Day021 Working Frozen-Topology "
            "Semantic Identity\n\n"
        )

        handle.write(
            "## Source-level comparison\n\n"
        )

        handle.write(
            f"- Raw TOP equality: "
            f"{raw_top_equal}.\n"
        )

        handle.write(
            f"- Comment/whitespace-normalized "
            f"equality: {canonical_top_equal}.\n"
        )

        handle.write(
            f"- Unified-diff lines: "
            f"{len(raw_diff_lines)}.\n\n"
        )

        handle.write(
            "## Processed TPR comparison\n\n"
        )

        handle.write(
            f"- `grompp`: "
            f"{'PASS' if grompp_pass else 'FAIL'}.\n"
        )

        handle.write(
            f"- `inputrec` identity: "
            f"{inputrec_equal}.\n"
        )

        handle.write(
            f"- Processed topology identity: "
            f"{processed_topology_equal}.\n"
        )

        handle.write(
            f"- Box identity: "
            f"{box_equal}.\n\n"
        )

        handle.write(
            "## Decision\n\n"
        )

        handle.write(
            f"- Working-copy semantic identity: "
            f"{'PASS' if semantic_identity else 'FAIL'}.\n"
        )

        handle.write(
            "- The accepted-directory topology "
            "remains the canonical provenance source "
            "regardless of working-copy equivalence.\n"
        )

    required_outputs = (
        WORKING_REBUILT_TPR,
        WORKING_PROCESSED_TOP,
        WORKING_MDOUT,
        WORKING_DUMP,
        WORKING_DUMP_STDERR,
        GROMPP_LOG,
        RAW_TOP_DIFF,
        SUMMARY_CSV,
        REPORT_MD,
    )

    missing_outputs = [
        path
        for path in required_outputs
        if not path.exists()
    ]

    if missing_outputs:
        raise RuntimeError(
            "Missing expected outputs:\n"
            + "\n".join(
                str(path)
                for path in missing_outputs
            )
        )

    log(
        "Day021 working frozen-topology "
        "semantic-identity verification completed."
    )

    log(
        f"Raw TOP equality: "
        f"{'PASS' if raw_top_equal else 'FAIL'}"
    )

    log(
        f"Comment/whitespace-normalized "
        f"TOP equality: "
        f"{'PASS' if canonical_top_equal else 'FAIL'}"
    )

    log(
        f"Raw diff lines: "
        f"{len(raw_diff_lines)}"
    )

    log(
        f"grompp reconstruction: "
        f"{'PASS' if grompp_pass else 'FAIL'}"
    )

    log(
        f"Inputrec identity: "
        f"{'PASS' if inputrec_equal else 'FAIL'}"
    )

    log(
        f"Processed topology identity: "
        f"{'PASS' if processed_topology_equal else 'FAIL'}"
    )

    log(
        f"Box identity: "
        f"{'PASS' if box_equal else 'FAIL'}"
    )

    log(
        f"Working-copy semantic identity: "
        f"{'PASS' if semantic_identity else 'FAIL'}"
    )

    log(
        "Canonical source: "
        f"{relative(ACCEPTED_TOP)}"
    )

    log(
        f"Wrote: {relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
