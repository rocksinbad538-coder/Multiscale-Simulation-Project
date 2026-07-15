#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import os
import re
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

ACCEPTED_TPR = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.tpr"
)

FINAL_GRO = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.gro"
)

ORIGINAL_MDP = (
    ACCEPTED_RUN_ROOT
    / "nvt_100ps_frozenSolute.mdp"
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

ACCEPTED_HBN_ITP = (
    ACCEPTED_TOPOLOGY_ROOT
    / "hbn_fixed_dummy.itp"
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

WORKING_HBN_ITP = (
    WORKING_TOPOLOGY_ROOT
    / "hbn_fixed_dummy.itp"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day021_accepted_hydrated_topology_audit/"
    "original_frozen_topology_identity"
)

AUDIT_MDP = (
    OUTPUT_ROOT
    / "nvt_100ps_frozenSolute_"
    "accepted_seed.mdp"
)

REBUILT_TPR = (
    OUTPUT_ROOT
    / "original_frozen_topology_"
    "rebuilt_from_final_gro.tpr"
)

PROCESSED_TOP = (
    OUTPUT_ROOT
    / "original_frozen_topology_processed.top"
)

MDOUT_MDP = (
    OUTPUT_ROOT
    / "original_frozen_topology_mdout.mdp"
)

ACCEPTED_DUMP = (
    OUTPUT_ROOT
    / "accepted_tpr_clean_dump.txt"
)

REBUILT_DUMP = (
    OUTPUT_ROOT
    / "rebuilt_original_frozen_"
    "topology_clean_dump.txt"
)

ACCEPTED_DUMP_STDERR = (
    OUTPUT_ROOT
    / "accepted_tpr_dump_stderr.log"
)

REBUILT_DUMP_STDERR = (
    OUTPUT_ROOT
    / "rebuilt_tpr_dump_stderr.log"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "grompp_original_frozen_topology.log"
)

GMX_CHECK_LOG = (
    OUTPUT_ROOT
    / "accepted_vs_original_frozen_"
    "rebuilt_gmx_check.log"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "original_frozen_topology_"
    "identity_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "ORIGINAL_FROZEN_TOPOLOGY_"
    "IDENTITY_DAY021.md"
)

ACCEPTED_LD_SEED = 2099210158
EXPECTED_TOTAL_ATOMS = 68320
EXPECTED_SOLUTE_ATOMS = 1784

LD_SEED_PATTERN = re.compile(
    r"^\s*ld[-_]seed\s*=.*$",
    flags=re.IGNORECASE,
)

NATOMS_PATTERN = re.compile(
    r"\bnatoms\s*=\s*(\d+)"
)

COORDINATE_PATTERN = re.compile(
    r"^\s*x\[\s*(\d+)\s*\]="
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


def text_hash(lines: list[str]) -> str:
    normalized = (
        "\n".join(
            line.rstrip()
            for line in lines
        )
        + "\n"
    )

    return hashlib.sha256(
        normalized.encode("utf-8")
    ).hexdigest()


def require_inputs() -> None:
    required = (
        ACCEPTED_TPR,
        FINAL_GRO,
        ORIGINAL_MDP,
        ACCEPTED_TOP,
        ACCEPTED_HBN_ITP,
        WORKING_TOP,
        WORKING_HBN_ITP,
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


def find_gmx() -> Path:
    preferred = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if preferred.exists():
        return preferred

    discovered = shutil.which("gmx")

    if discovered is None:
        raise RuntimeError(
            "Could not locate the GROMACS "
            "gmx executable"
        )

    return Path(discovered)


def write_audit_mdp() -> None:
    lines = ORIGINAL_MDP.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    output_lines: list[str] = []
    replacement_count = 0

    for line in lines:
        if LD_SEED_PATTERN.match(line):
            output_lines.append(
                f"ld-seed = {ACCEPTED_LD_SEED}"
            )

            replacement_count += 1
        else:
            output_lines.append(line)

    if replacement_count == 0:
        output_lines.append("")
        output_lines.append(
            f"ld-seed = {ACCEPTED_LD_SEED}"
        )

    if replacement_count > 1:
        raise RuntimeError(
            "Multiple ld-seed entries were "
            "found in the accepted MDP"
        )

    AUDIT_MDP.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )


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
    output_path: Path,
    stderr_path: Path,
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

    output_path.write_text(
        result.stdout,
        encoding="utf-8",
    )

    stderr_path.write_text(
        result.stderr,
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"gmx dump failed for {tpr_path}"
        )


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

    return lines[
        start_index:end_index
    ]


def extract_coordinate_section(
    text: str,
) -> list[str]:
    lines = text.splitlines()

    start_index: int | None = None
    end_index: int | None = None

    for index, line in enumerate(lines):
        stripped = line.strip()

        if (
            stripped.startswith("x (")
            and stripped.endswith("):")
        ):
            start_index = index
            break

    if start_index is None:
        raise RuntimeError(
            "Coordinate section was not found"
        )

    for index in range(
        start_index + 1,
        len(lines),
    ):
        stripped = lines[index].strip()

        if (
            stripped.startswith("v (")
            and stripped.endswith("):")
        ):
            end_index = index
            break

    if end_index is None:
        raise RuntimeError(
            "Velocity section boundary "
            "was not found"
        )

    return lines[
        start_index:end_index
    ]


def parse_natoms(text: str) -> int:
    match = NATOMS_PATTERN.search(text)

    if match is None:
        raise RuntimeError(
            "Could not parse natoms from dump"
        )

    return int(match.group(1))


def parse_assignments(
    lines: list[str],
) -> dict[str, str]:
    result: dict[str, str] = {}

    for line in lines:
        if "=" not in line:
            continue

        key, value = line.split(
            "=",
            1,
        )

        normalized_key = " ".join(
            key.split()
        )

        normalized_value = " ".join(
            value.split()
        )

        result[
            normalized_key
        ] = normalized_value

    return result


def first_coordinate_difference(
    accepted_lines: list[str],
    rebuilt_lines: list[str],
) -> int | None:
    accepted_coordinates = [
        line.rstrip()
        for line in accepted_lines
        if COORDINATE_PATTERN.match(line)
    ]

    rebuilt_coordinates = [
        line.rstrip()
        for line in rebuilt_lines
        if COORDINATE_PATTERN.match(line)
    ]

    if (
        len(accepted_coordinates)
        != len(rebuilt_coordinates)
    ):
        raise RuntimeError(
            "Coordinate-array length mismatch: "
            f"{len(accepted_coordinates)} versus "
            f"{len(rebuilt_coordinates)}"
        )

    for (
        accepted_line,
        rebuilt_line,
    ) in zip(
        accepted_coordinates,
        rebuilt_coordinates,
    ):
        if accepted_line == rebuilt_line:
            continue

        accepted_match = (
            COORDINATE_PATTERN.match(
                accepted_line
            )
        )

        if accepted_match is None:
            raise RuntimeError(
                "Could not parse coordinate index"
            )

        return int(
            accepted_match.group(1)
        )

    return None


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

    require_inputs()

    gmx = find_gmx()

    accepted_top_matches_working = (
        sha256(ACCEPTED_TOP)
        == sha256(WORKING_TOP)
    )

    accepted_hbn_matches_working = (
        sha256(ACCEPTED_HBN_ITP)
        == sha256(WORKING_HBN_ITP)
    )

    write_audit_mdp()

    grompp_result = run_combined(
        [
            str(gmx),
            "grompp",
            "-f",
            str(AUDIT_MDP),
            "-c",
            str(FINAL_GRO),
            "-p",
            str(ACCEPTED_TOP),
            "-o",
            str(REBUILT_TPR),
            "-po",
            str(MDOUT_MDP),
            "-pp",
            str(PROCESSED_TOP),
            "-maxwarn",
            "0",
        ],
        ACCEPTED_TOPOLOGY_ROOT,
        GROMPP_LOG,
    )

    grompp_pass = (
        grompp_result.returncode == 0
        and REBUILT_TPR.exists()
        and REBUILT_TPR.stat().st_size > 0
    )

    if not grompp_pass:
        raise RuntimeError(
            "Original frozen-topology "
            "grompp reconstruction failed. "
            f"See {GROMPP_LOG}"
        )

    run_dump(
        gmx,
        ACCEPTED_TPR,
        ACCEPTED_DUMP,
        ACCEPTED_DUMP_STDERR,
    )

    run_dump(
        gmx,
        REBUILT_TPR,
        REBUILT_DUMP,
        REBUILT_DUMP_STDERR,
    )

    accepted_text = (
        ACCEPTED_DUMP.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    rebuilt_text = (
        REBUILT_DUMP.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    accepted_natoms = parse_natoms(
        accepted_text
    )

    rebuilt_natoms = parse_natoms(
        rebuilt_text
    )

    accepted_inputrec = extract_between(
        accepted_text,
        "inputrec:",
        "qm-opts:",
    )

    rebuilt_inputrec = extract_between(
        rebuilt_text,
        "inputrec:",
        "qm-opts:",
    )

    accepted_topology = extract_between(
        accepted_text,
        "topology:",
        "box (3x3):",
    )

    rebuilt_topology = extract_between(
        rebuilt_text,
        "topology:",
        "box (3x3):",
    )

    accepted_box = extract_between(
        accepted_text,
        "box (3x3):",
        "box_rel (3x3):",
    )

    rebuilt_box = extract_between(
        rebuilt_text,
        "box (3x3):",
        "box_rel (3x3):",
    )

    accepted_coordinates = (
        extract_coordinate_section(
            accepted_text
        )
    )

    rebuilt_coordinates = (
        extract_coordinate_section(
            rebuilt_text
        )
    )

    inputrec_equal = (
        text_hash(accepted_inputrec)
        == text_hash(rebuilt_inputrec)
    )

    topology_equal = (
        text_hash(accepted_topology)
        == text_hash(rebuilt_topology)
    )

    box_equal = (
        text_hash(accepted_box)
        == text_hash(rebuilt_box)
    )

    first_coordinate_difference_index = (
        first_coordinate_difference(
            accepted_coordinates,
            rebuilt_coordinates,
        )
    )

    frozen_solute_prefix_pass = (
        first_coordinate_difference_index
        == EXPECTED_SOLUTE_ATOMS
    )

    accepted_parameters = (
        parse_assignments(
            accepted_inputrec
        )
    )

    rebuilt_parameters = (
        parse_assignments(
            rebuilt_inputrec
        )
    )

    parameter_keys = sorted(
        set(accepted_parameters)
        | set(rebuilt_parameters)
    )

    differing_parameters = [
        key
        for key in parameter_keys
        if (
            accepted_parameters.get(
                key,
                "<MISSING>",
            )
            != rebuilt_parameters.get(
                key,
                "<MISSING>",
            )
        )
    ]

    gmx_check_result = run_combined(
        [
            str(gmx),
            "check",
            "-s1",
            str(ACCEPTED_TPR),
            "-s2",
            str(REBUILT_TPR),
        ],
        PROJECT_ROOT,
        GMX_CHECK_LOG,
    )

    atom_count_pass = (
        accepted_natoms
        == rebuilt_natoms
        == EXPECTED_TOTAL_ATOMS
    )

    exact_original_topology_identity = all(
        (
            grompp_pass,
            atom_count_pass,
            inputrec_equal,
            topology_equal,
            box_equal,
        )
    )

    provenance_pass = all(
        (
            exact_original_topology_identity,
            frozen_solute_prefix_pass,
        )
    )

    summary = {
        "accepted_TPR": (
            relative(ACCEPTED_TPR)
        ),
        "candidate_original_TOP": (
            relative(ACCEPTED_TOP)
        ),
        "candidate_original_HBN_ITP": (
            relative(ACCEPTED_HBN_ITP)
        ),
        "working_copy_TOP": (
            relative(WORKING_TOP)
        ),
        "working_copy_HBN_ITP": (
            relative(WORKING_HBN_ITP)
        ),
        "accepted_TOP_matches_working_copy": (
            accepted_top_matches_working
        ),
        "accepted_HBN_ITP_matches_working_copy": (
            accepted_hbn_matches_working
        ),
        "accepted_TPR_atom_count": (
            accepted_natoms
        ),
        "rebuilt_TPR_atom_count": (
            rebuilt_natoms
        ),
        "grompp_return_code": (
            grompp_result.returncode
        ),
        "grompp_pass": grompp_pass,
        "inputrec_exactly_equal": (
            inputrec_equal
        ),
        "differing_inputrec_parameter_count": (
            len(differing_parameters)
        ),
        "differing_inputrec_parameters": (
            ",".join(differing_parameters)
        ),
        "topology_exactly_equal": (
            topology_equal
        ),
        "box_exactly_equal": (
            box_equal
        ),
        "first_coordinate_difference_atom_index": (
            first_coordinate_difference_index
        ),
        "expected_first_water_atom_index": (
            EXPECTED_SOLUTE_ATOMS
        ),
        "frozen_HBN_PYR_coordinate_prefix_pass": (
            frozen_solute_prefix_pass
        ),
        "gmx_check_return_code": (
            gmx_check_result.returncode
        ),
        "exact_original_topology_identity": (
            exact_original_topology_identity
        ),
        "provenance_validation": (
            "PASS"
            if provenance_pass
            else "FAIL"
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
            "# Day021 Original Frozen-Topology "
            "Identity Verification\n\n"
        )

        handle.write(
            "## Candidate\n\n"
        )

        handle.write(
            f"- TOP: `{relative(ACCEPTED_TOP)}`\n"
        )

        handle.write(
            f"- HBN ITP: "
            f"`{relative(ACCEPTED_HBN_ITP)}`\n"
        )

        handle.write(
            f"- Explicit reconstructed "
            f"`ld-seed`: {ACCEPTED_LD_SEED}\n\n"
        )

        handle.write(
            "## Exact TPR-section comparison\n\n"
        )

        handle.write(
            f"- Atom count: "
            f"{accepted_natoms} versus "
            f"{rebuilt_natoms}.\n"
        )

        handle.write(
            f"- `inputrec` exact equality: "
            f"{inputrec_equal}.\n"
        )

        handle.write(
            f"- Differing `inputrec` parameters: "
            f"{len(differing_parameters)}.\n"
        )

        handle.write(
            f"- `topology` exact equality: "
            f"{topology_equal}.\n"
        )

        handle.write(
            f"- `box` exact equality: "
            f"{box_equal}.\n"
        )

        handle.write(
            f"- First coordinate difference: "
            f"atom index "
            f"{first_coordinate_difference_index}.\n"
        )

        handle.write(
            f"- Frozen HBN+PYR coordinate prefix: "
            f"{frozen_solute_prefix_pass}.\n\n"
        )

        handle.write(
            "## Interpretation\n\n"
        )

        if provenance_pass:
            handle.write(
                "The unbonded HBN topology stored "
                "under the accepted parameter set is "
                "the exact topology embedded in the "
                "accepted frozen-solute TPR. The first "
                "coordinate difference occurs at atom "
                "1784, the first water atom, because "
                "the reconstruction uses the final "
                "100 ps GRO while the accepted TPR "
                "stores the initial production state.\n"
            )
        else:
            handle.write(
                "Exact provenance has not been "
                "established. No restraint-release "
                "workflow should be initiated until "
                "the failed comparison is resolved.\n"
            )

        handle.write(
            "\n## Decision\n\n"
        )

        handle.write(
            f"- Exact original-topology identity: "
            f"{'PASS' if exact_original_topology_identity else 'FAIL'}.\n"
        )

        handle.write(
            f"- Provenance validation: "
            f"{'PASS' if provenance_pass else 'FAIL'}.\n"
        )

        handle.write(
            "- Mobile bonded-HBN topology remains a "
            "separate model and must be validated "
            "independently before mobile production.\n"
        )

    required_outputs = (
        AUDIT_MDP,
        REBUILT_TPR,
        PROCESSED_TOP,
        MDOUT_MDP,
        ACCEPTED_DUMP,
        REBUILT_DUMP,
        ACCEPTED_DUMP_STDERR,
        REBUILT_DUMP_STDERR,
        GROMPP_LOG,
        GMX_CHECK_LOG,
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
            "Missing required outputs:\n"
            + "\n".join(
                str(path)
                for path in missing_outputs
            )
        )

    log(
        "Day021 original frozen-topology "
        "identity verification completed."
    )

    log(
        f"Accepted-directory TOP matches "
        f"working-copy TOP: "
        f"{accepted_top_matches_working}"
    )

    log(
        f"Accepted-directory HBN ITP matches "
        f"working-copy HBN ITP: "
        f"{accepted_hbn_matches_working}"
    )

    log(
        f"Accepted/rebuilt atoms: "
        f"{accepted_natoms}/{rebuilt_natoms}"
    )

    log(
        f"Inputrec exact equality: "
        f"{'PASS' if inputrec_equal else 'FAIL'}"
    )

    log(
        f"Differing inputrec parameters: "
        f"{len(differing_parameters)}"
    )

    if differing_parameters:
        log(
            "Differing parameter names: "
            + ", ".join(
                differing_parameters
            )
        )

    log(
        f"Topology exact equality: "
        f"{'PASS' if topology_equal else 'FAIL'}"
    )

    log(
        f"Box exact equality: "
        f"{'PASS' if box_equal else 'FAIL'}"
    )

    log(
        f"First coordinate difference "
        f"atom index: "
        f"{first_coordinate_difference_index}"
    )

    log(
        f"Frozen HBN+PYR coordinate prefix: "
        f"{'PASS' if frozen_solute_prefix_pass else 'FAIL'}"
    )

    log(
        f"Exact original-topology identity: "
        f"{'PASS' if exact_original_topology_identity else 'FAIL'}"
    )

    log(
        f"Provenance validation: "
        f"{'PASS' if provenance_pass else 'FAIL'}"
    )

    log(
        f"Wrote: {relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
