#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

AUDIT_ROOT = (
    DAY023_ROOT
    / "00_input_audit"
)

R0_ACCEPTED = (
    ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

PROTOCOL_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

STAGE02_START_GRO = (
    PROTOCOL_ROOT
    / "execution/01_em_k10000/"
    "01_em_k10000.gro"
)

STAGE02_TPR = (
    PROTOCOL_ROOT
    / "execution/02_nvt_k10000_1ps/"
    "02_nvt_k10000_1ps.tpr"
)

STAGE08_TPR = (
    PROTOCOL_ROOT
    / "execution/08_nvt_mobile_100ps/"
    "08_nvt_mobile_100ps.tpr"
)

STAGE08_XTC = (
    PROTOCOL_ROOT
    / "execution/08_nvt_mobile_100ps/"
    "08_nvt_mobile_100ps.xtc"
)

MATCHED_CONTROL_TPR = (
    PROTOCOL_ROOT
    / "matched_frozen_control_144ps/execution/"
    "matched_frozen_control_144ps.tpr"
)

MATCHED_CONTROL_XTC = (
    PROTOCOL_ROOT
    / "matched_frozen_control_144ps/execution/"
    "matched_frozen_control_144ps.xtc"
)

VALIDATED_TOPOLOGY = (
    ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_improperGeo100_validated/"
    "hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_"
    "hbnBonded_kang2000_improperGeo100.top"
)

INPUT_INDEX = (
    PROTOCOL_ROOT
    / "protocol_inputs/mobile_release_index.ndx"
)

INPUT_MDP_STAGE02 = (
    PROTOCOL_ROOT
    / "protocol_inputs/mdp/"
    "02_nvt_k10000_1ps.mdp"
)

INPUT_MDP_STAGE08 = (
    PROTOCOL_ROOT
    / "protocol_inputs/mdp/"
    "08_nvt_mobile_100ps.mdp"
)

INVENTORY_CSV = (
    AUDIT_ROOT
    / "day023_confinement_input_inventory.csv"
)

REPORT_MD = (
    AUDIT_ROOT
    / "R1_DESIGN_INPUT_AUDIT_DAY023.md"
)

CONTRACT_MD = (
    AUDIT_ROOT
    / "R1_CONFINEMENT_DESIGN_CONTRACT_DRAFT_DAY023.md"
)

EXPECTED_ATOMS = 68320

ALLOWED_EXTENSIONS = {
    ".gro",
    ".tpr",
    ".xtc",
    ".trr",
    ".cpt",
    ".top",
    ".itp",
    ".mdp",
    ".ndx",
    ".csv",
    ".md",
    ".log",
}


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def run_command(
    command: list[str],
    cwd: Path | None = None,
) -> tuple[int, str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    return (
        completed.returncode,
        completed.stdout,
    )


def git_value(
    arguments: list[str],
) -> str:
    return_code, output = run_command(
        ["git", *arguments],
        cwd=ROOT,
    )

    if return_code != 0:
        return "UNAVAILABLE"

    return output.strip()


def gro_atom_count(
    path: Path,
) -> int | None:
    if (
        not path.exists()
        or path.suffix.lower() != ".gro"
    ):
        return None

    try:
        with path.open(
            "r",
            encoding="utf-8",
            errors="replace",
        ) as handle:
            handle.readline()
            return int(
                handle.readline().strip()
            )
    except (
        OSError,
        ValueError,
    ):
        return None


def classify_role(
    path: Path,
) -> str:
    name = path.name.lower()

    if path.suffix.lower() == ".xtc":
        return "trajectory"

    if path.suffix.lower() == ".tpr":
        return "run_input"

    if path.suffix.lower() == ".cpt":
        return "checkpoint"

    if path.suffix.lower() == ".gro":
        if any(
            token in name
            for token in (
                "initial",
                "start",
                "input",
                "em",
                "0ps",
                "conf",
            )
        ):
            return "coordinate_candidate"

        return "coordinate_file"

    if path.suffix.lower() in {
        ".top",
        ".itp",
    }:
        return "topology"

    if path.suffix.lower() == ".mdp":
        return "md_parameters"

    if path.suffix.lower() == ".ndx":
        return "index"

    if path.suffix.lower() in {
        ".csv",
        ".md",
    }:
        return "analysis_or_report"

    return "supporting_file"


def file_record(
    group: str,
    path: Path,
    role: str | None = None,
) -> dict[str, Any]:
    exists = (
        path.exists()
        and path.is_file()
    )

    return {
        "group": group,
        "role": (
            role
            if role is not None
            else classify_role(path)
        ),
        "path": relative(path),
        "exists": exists,
        "size_bytes": (
            path.stat().st_size
            if exists
            else 0
        ),
        "size_MiB": (
            path.stat().st_size
            / (1024.0 ** 2)
            if exists
            else 0.0
        ),
        "gro_atom_count": (
            gro_atom_count(path)
            if exists
            else ""
        ),
    }


def discover_files(
    root: Path,
    maximum: int = 5000,
) -> list[Path]:
    if not root.exists():
        return []

    files = []

    for path in root.rglob("*"):
        if not path.is_file():
            continue

        if path.suffix.lower() not in ALLOWED_EXTENSIONS:
            continue

        files.append(path)

        if len(files) >= maximum:
            break

    return sorted(
        files,
        key=lambda item: str(item),
    )


def find_gmx() -> str | None:
    executable = shutil.which("gmx")

    if executable:
        return executable

    fallback = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if fallback.exists():
        return str(fallback)

    return None


def parse_xtc_check(
    gmx: str,
    xtc: Path,
) -> dict[str, Any]:
    return_code, output = run_command(
        [
            gmx,
            "check",
            "-f",
            str(xtc),
        ],
        cwd=xtc.parent,
    )

    result: dict[str, Any] = {
        "return_code": return_code,
        "atoms": "",
        "frames": "",
        "interval_ps": "",
        "duration_ps": "",
    }

    if return_code != 0:
        return result

    atom_match = re.search(
        r"(?mi)^\s*#\s*Atoms\s+(\d+)\s*$",
        output,
    )

    coords_match = re.search(
        r"(?mi)^\s*Coords\s+(\d+)\s+"
        r"([-+0-9.eE]+)\s*$",
        output,
    )

    if atom_match is not None:
        result["atoms"] = int(
            atom_match.group(1)
        )

    if coords_match is not None:
        frames = int(
            coords_match.group(1)
        )

        interval = float(
            coords_match.group(2)
        )

        result["frames"] = frames
        result["interval_ps"] = interval
        result["duration_ps"] = (
            (frames - 1)
            * interval
        )

    return result


def write_csv(
    rows: list[dict[str, Any]],
) -> None:
    fields = [
        "group",
        "role",
        "path",
        "exists",
        "size_bytes",
        "size_MiB",
        "gro_atom_count",
    ]

    with INVENTORY_CSV.open(
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
            writer.writerow(row)


def status_word(
    condition: bool,
) -> str:
    return (
        "PASS"
        if condition
        else "MISSING"
    )


def main() -> None:
    AUDIT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    records: list[dict[str, Any]] = []

    authoritative_inputs = {
        "validated_topology": VALIDATED_TOPOLOGY,
        "stage02_start_gro": STAGE02_START_GRO,
        "stage02_tpr": STAGE02_TPR,
        "stage08_tpr": STAGE08_TPR,
        "stage08_xtc": STAGE08_XTC,
        "matched_control_tpr": MATCHED_CONTROL_TPR,
        "matched_control_xtc": MATCHED_CONTROL_XTC,
        "input_index": INPUT_INDEX,
        "stage02_mdp": INPUT_MDP_STAGE02,
        "stage08_mdp": INPUT_MDP_STAGE08,
    }

    for role, path in authoritative_inputs.items():
        records.append(
            file_record(
                "known_authoritative_inputs",
                path,
                role,
            )
        )

    accepted_files = discover_files(
        R0_ACCEPTED
    )

    for path in accepted_files:
        records.append(
            file_record(
                "accepted_frozen_r0",
                path,
            )
        )

    script_candidates = sorted(
        {
            *(
                ROOT
                / "scripts/phase1A"
            ).glob("*hbn*.py"),
            *(
                ROOT
                / "scripts/phase1A"
            ).glob("*structure*.py"),
            *(
                ROOT
                / "scripts/phase1A"
            ).glob("*topolog*.py"),
            *(
                ROOT
                / "scripts/phase1A"
            ).glob("*cap*.py"),
            *(
                ROOT
                / "scripts/phase1A"
            ).glob("*water*.py"),
        },
        key=lambda item: str(item),
    )

    for path in script_candidates:
        if path.is_file():
            records.append(
                file_record(
                    "reusable_script_candidates",
                    path,
                )
            )

    write_csv(
        records
    )

    accepted_gro = [
        path
        for path in accepted_files
        if path.suffix.lower() == ".gro"
    ]

    accepted_tpr = [
        path
        for path in accepted_files
        if path.suffix.lower() == ".tpr"
    ]

    accepted_xtc = [
        path
        for path in accepted_files
        if path.suffix.lower() == ".xtc"
    ]

    accepted_cpt = [
        path
        for path in accepted_files
        if path.suffix.lower() == ".cpt"
    ]

    correct_atom_gro = [
        path
        for path in accepted_gro
        if gro_atom_count(path)
        == EXPECTED_ATOMS
    ]

    cap_scripts = [
        path
        for path in script_candidates
        if "cap" in path.name.lower()
    ]

    gmx = find_gmx()

    xtc_results = []

    if gmx is not None:
        for xtc in accepted_xtc:
            metadata = parse_xtc_check(
                gmx,
                xtc,
            )

            xtc_results.append(
                {
                    "path": relative(xtc),
                    **metadata,
                }
            )

    hydrated_start_resolved = bool(
        len(accepted_xtc) == 1
        and len(accepted_tpr) >= 1
        and len(xtc_results) == 1
        and xtc_results[0][
            "frames"
        ] not in (
            "",
            None,
        )
    )

    if hydrated_start_resolved:
        next_step = (
            "EXTRACT_AND_AUDIT_ACCEPTED_R0_T0_"
            "HYDRATED_STATE"
        )
    else:
        next_step = (
            "RESOLVE_ACCEPTED_R0_HYDRATED_"
            "START_STATE"
        )

    git_root = git_value(
        [
            "rev-parse",
            "--show-toplevel",
        ]
    )

    git_branch = git_value(
        [
            "branch",
            "--show-current",
        ]
    )

    git_head = git_value(
        [
            "rev-parse",
            "HEAD",
        ]
    )

    accepted_file_lines = "\n".join(
        (
            f"- `{relative(path)}`"
            f" — {classify_role(path)}"
            + (
                f"; atoms={gro_atom_count(path)}"
                if path.suffix.lower() == ".gro"
                else ""
            )
        )
        for path in accepted_files
    )

    if not accepted_file_lines:
        accepted_file_lines = (
            "- No accepted R0 files discovered."
        )

    xtc_lines = "\n".join(
        (
            f"- `{item['path']}`: "
            f"return code={item['return_code']}; "
            f"atoms={item['atoms']}; "
            f"frames={item['frames']}; "
            f"interval={item['interval_ps']} ps; "
            f"duration={item['duration_ps']} ps"
        )
        for item in xtc_results
    )

    if not xtc_lines:
        xtc_lines = (
            "- No accepted R0 XTC metadata available."
        )

    script_lines = "\n".join(
        f"- `{relative(path)}`"
        for path in script_candidates
    )

    if not script_lines:
        script_lines = (
            "- No reusable structure/topology scripts discovered."
        )

    REPORT_MD.write_text(
        f"""# Day023 R1 Design Input Audit

## Repository state

- Project root: `{relative(ROOT)}`
- Git root: `{git_root}`
- Branch: `{git_branch}`
- HEAD: `{git_head}`

## Known authoritative inputs

- Validated topology:
  **{status_word(VALIDATED_TOPOLOGY.exists())}**
  `{relative(VALIDATED_TOPOLOGY)}`
- Stage02 starting coordinates:
  **{status_word(STAGE02_START_GRO.exists())}**
  `{relative(STAGE02_START_GRO)}`
- Stage02 TPR:
  **{status_word(STAGE02_TPR.exists())}**
- Stage08 mobile TPR/XTC:
  **{status_word(STAGE08_TPR.exists() and STAGE08_XTC.exists())}**
- Matched frozen-control TPR/XTC:
  **{status_word(MATCHED_CONTROL_TPR.exists() and MATCHED_CONTROL_XTC.exists())}**
- Index file:
  **{status_word(INPUT_INDEX.exists())}**

## Accepted frozen R0 directory

- Directory:
  `{relative(R0_ACCEPTED)}`
- Directory status:
  **{status_word(R0_ACCEPTED.exists())}**
- GRO files: **{len(accepted_gro)}**
- GRO files with {EXPECTED_ATOMS} atoms:
  **{len(correct_atom_gro)}**
- TPR files: **{len(accepted_tpr)}**
- XTC files: **{len(accepted_xtc)}**
- CPT files: **{len(accepted_cpt)}**

### Accepted R0 files

{accepted_file_lines}

## Accepted trajectory metadata

{xtc_lines}

## Existing reusable workflow candidates

{script_lines}

## Cap-generation status

- Existing cap-related scripts:
  **{len(cap_scripts)}**
- Authoritative hydrated R1 starting state resolved:
  **{'YES' if hydrated_start_resolved else 'NO'}**

## Decision

- R0 remains the validated open-tube reference.
- No existing R0 file is modified by this audit.
- No MD or QM calculation was executed.
- Next required step:
  **{next_step}**

R1 must be generated from the earliest hydrated state of the accepted
frozen-solute R0 trajectory, not from the depleted mobile branch state.
""",
        encoding="utf-8",
    )

    CONTRACT_MD.write_text(
        """# R1 Confinement Design Contract — Draft

## Purpose

R1 is a fully capped positive control for persistent lumen-water
confinement. It is not assumed to be the final device architecture.

## Invariants inherited from R0

R1 must preserve, unless an explicit audit demonstrates otherwise:

- the validated h-BN scaffold;
- all four pyrene chromophores;
- atom ordering for the inherited R0 atoms;
- the accepted force-field parameters;
- TIP4P/2005 water;
- the simulation box;
- the 300 K reference condition;
- the established lumen-axis convention;
- the existing structural-analysis definitions.

## Starting-state requirement

R1 must start from the earliest hydrated frame of the accepted
frozen-solute R0 trajectory.

The depleted Stage02/mobile-branch state must not be used as the
primary R1 positive-control starting state.

## Cap requirements

The initial R1 caps must:

- close both axial exits;
- avoid atomic overlaps with HBN, PYR, and water;
- preserve the accessible lumen interior;
- introduce no uncontrolled net charge;
- have explicitly documented composition and bonding;
- be generated reproducibly by a script under `scripts/phase1A/`.

## Gate sequence

1. Identify and extract the authoritative hydrated R0 starting state.
2. Quantify lumen axis, end planes, radius, and accessible volume.
3. Generate an R1 cap prototype.
4. Audit composition, bonding, overlaps, charge, and geometry.
5. Minimize and perform static preprocessing.
6. Run only a short frozen-solute confinement screening.
7. Authorize longer or mobile simulations only after screening passes.

## Current prohibitions

- No long mobile MD.
- No multitemperature production.
- No new QM calculations.
- No replacement of the accepted R0 baseline.
- No interpretation of R1 as a final experimental architecture before
  R2 and R3 are evaluated.

## Provisional scientific objective

R1 should demonstrate that eliminating axial escape prevents the
progressive depletion observed in R0. Its role is to validate the
confinement methodology and establish a positive-control reference.
""",
        encoding="utf-8",
    )

    print(
        "Day023 confinement-design input audit completed."
    )

    print(
        "Git branch / HEAD: "
        f"{git_branch} / {git_head}"
    )

    print(
        "Accepted R0 directory: "
        f"{status_word(R0_ACCEPTED.exists())}"
    )

    print(
        "Accepted R0 GRO/TPR/XTC/CPT counts: "
        f"{len(accepted_gro)}/"
        f"{len(accepted_tpr)}/"
        f"{len(accepted_xtc)}/"
        f"{len(accepted_cpt)}"
    )

    print(
        "Accepted GRO files with 68320 atoms: "
        f"{len(correct_atom_gro)}"
    )

    print(
        "Validated topology / Stage02 start / "
        "Stage08 trajectory: "
        f"{status_word(VALIDATED_TOPOLOGY.exists())}/"
        f"{status_word(STAGE02_START_GRO.exists())}/"
        f"{status_word(STAGE08_XTC.exists())}"
    )

    print(
        "Existing reusable script candidates: "
        f"{len(script_candidates)}"
    )

    print(
        "Existing cap-related scripts: "
        f"{len(cap_scripts)}"
    )

    print(
        "Authoritative hydrated R1 start state resolved: "
        f"{'YES' if hydrated_start_resolved else 'NO'}"
    )

    for item in xtc_results:
        print(
            "Accepted XTC metadata: "
            f"{item['frames']} frames / "
            f"{item['interval_ps']} ps / "
            f"{item['duration_ps']} ps / "
            f"{item['atoms']} atoms"
        )

    print(
        f"Required next step: {next_step}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    print(
        f"Wrote: {relative(CONTRACT_MD)}"
    )

    print(
        f"Wrote: {relative(INVENTORY_CSV)}"
    )


if __name__ == "__main__":
    main()
