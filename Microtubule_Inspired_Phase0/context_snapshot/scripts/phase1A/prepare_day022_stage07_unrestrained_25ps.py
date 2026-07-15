#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

STAGE06 = "06_nvt_unrestrained_10ps"
STAGE07 = "07_nvt_unrestrained_25ps"

INPUT_ROOT = PROTOCOL / "protocol_inputs"
MDP_ROOT = INPUT_ROOT / "mdp"
TOPOLOGY_ROOT = INPUT_ROOT / "topology"

RUN06 = PROTOCOL / "execution" / STAGE06

SOURCE_MDP = MDP_ROOT / f"{STAGE06}.mdp"
TARGET_MDP = MDP_ROOT / f"{STAGE07}.mdp"

START_GRO = RUN06 / f"{STAGE06}.gro"
START_CPT = RUN06 / f"{STAGE06}.cpt"

TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    INPUT_ROOT
    / "mobile_release_index.ndx"
)

OPERATIONAL_SUMMARY = (
    RUN06
    / f"{STAGE06}_summary.csv"
)

STRUCTURAL_SUMMARY = (
    RUN06
    / f"{STAGE06}_structural_summary.csv"
)

PERSISTENCE_SUMMARY = (
    RUN06
    / "stage06_hbn_improper_persistence_summary.csv"
)

STATIC_ROOT = (
    PROTOCOL
    / "static_validation"
    / STAGE07
)

STATIC_TPR = STATIC_ROOT / f"{STAGE07}.tpr"
STATIC_TOP = STATIC_ROOT / f"{STAGE07}_processed.top"
STATIC_MDOUT = STATIC_ROOT / f"{STAGE07}_mdout.mdp"
STATIC_LOG = STATIC_ROOT / f"{STAGE07}_grompp.log"

PREPARATION_SUMMARY = (
    STATIC_ROOT
    / f"{STAGE07}_preparation_summary.csv"
)

PREPARATION_REPORT = (
    STATIC_ROOT
    / "STAGE07_PREPARATION_DAY022.md"
)


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly one row in {path}"
        )

    return rows[0]


def require_inputs() -> None:
    required = (
        SOURCE_MDP,
        START_GRO,
        START_CPT,
        TOP,
        NDX,
        OPERATIONAL_SUMMARY,
        STRUCTURAL_SUMMARY,
        PERSISTENCE_SUMMARY,
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
            "Missing or empty required files:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def validate_stage06_authorization() -> None:
    operational = read_single_row(
        OPERATIONAL_SUMMARY
    )

    structural = read_single_row(
        STRUCTURAL_SUMMARY
    )

    persistence = read_single_row(
        PERSISTENCE_SUMMARY
    )

    failures = []

    if (
        operational.get("decision", "")
        .strip()
        .upper()
        != "PASS"
    ):
        failures.append(
            "Stage06 operational decision is not PASS"
        )

    if (
        structural.get(
            "structural_screen",
            "",
        )
        .strip()
        .upper()
        != "STABLE_CANDIDATE"
    ):
        failures.append(
            "Stage06 structural screen is not STABLE_CANDIDATE"
        )

    if (
        persistence.get(
            "targeted_decision",
            "",
        )
        .strip()
        .upper()
        != "PASS_WITH_MONITORING"
    ):
        failures.append(
            "Stage06 persistence decision is not PASS_WITH_MONITORING"
        )

    if (
        persistence.get(
            "authorized_next_step",
            "",
        )
        .strip()
        != "25PS_EXTENDED_UNRESTRAINED_VALIDATION"
    ):
        failures.append(
            "The 25 ps extended validation is not authorized"
        )

    if (
        persistence.get(
            "blocked_reasons",
            "",
        )
        .strip()
    ):
        failures.append(
            "Stage06 contains blocking reasons"
        )

    if failures:
        raise RuntimeError(
            "Stage06 authorization failed:\n"
            + "\n".join(failures)
        )


def set_mdp_value(
    text: str,
    key: str,
    value: str,
) -> str:
    pattern = re.compile(
        rf"(?mi)^(\s*{re.escape(key)}\s*=\s*).*$"
    )

    if pattern.search(text):
        return pattern.sub(
            rf"\g<1>{value}",
            text,
            count=1,
        )

    return (
        text.rstrip()
        + f"\n{key} = {value}\n"
    )


def create_stage07_mdp() -> None:
    text = SOURCE_MDP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    text = set_mdp_value(
        text,
        "nsteps",
        "50000",
    )

    text = set_mdp_value(
        text,
        "dt",
        "0.0005",
    )

    text = set_mdp_value(
        text,
        "continuation",
        "yes",
    )

    text = set_mdp_value(
        text,
        "gen-vel",
        "no",
    )

    text = set_mdp_value(
        text,
        "nstxout-compressed",
        "1000",
    )

    text = set_mdp_value(
        text,
        "nstenergy",
        "200",
    )

    text = set_mdp_value(
        text,
        "nstlog",
        "1000",
    )

    text = set_mdp_value(
        text,
        "nstcalcenergy",
        "100",
    )

    parsed = {}

    for raw_line in text.splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line or "=" not in line:
            continue

        key, value = line.split(
            "=",
            1,
        )

        parsed[
            key.strip().lower()
        ] = value.strip().lower()

    expected = {
        "integrator": "md",
        "nsteps": "50000",
        "dt": "0.0005",
        "continuation": "yes",
        "gen-vel": "no",
        "nstxout-compressed": "1000",
    }

    failures = []

    for key, expected_value in expected.items():
        actual = parsed.get(
            key,
            "",
        )

        if actual != expected_value:
            failures.append(
                f"{key}: expected {expected_value}, "
                f"found {actual or 'MISSING'}"
            )

    define_value = parsed.get(
        "define",
        "",
    )

    if "posres" in define_value:
        failures.append(
            f"Position restraints remain active: {define_value}"
        )

    if failures:
        raise RuntimeError(
            "Stage07 MDP validation failed:\n"
            + "\n".join(failures)
        )

    text = (
        text.rstrip()
        + "\n\n"
        + "; Day022 extended unrestrained mobile validation\n"
        + "; 50000 steps x 0.0005 ps = 25 ps\n"
        + "; compressed coordinates every 0.5 ps\n"
    )

    if TARGET_MDP.exists():
        existing = TARGET_MDP.read_text(
            encoding="utf-8",
            errors="replace",
        )

        if existing != text:
            raise RuntimeError(
                "Stage07 MDP already exists with different content"
            )
    else:
        TARGET_MDP.write_text(
            text,
            encoding="utf-8",
        )


def count_position_restraints(
    path: Path,
) -> int:
    section = ""
    count = 0

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line:
            continue

        if (
            line.startswith("[")
            and line.endswith("]")
        ):
            section = (
                line[1:-1]
                .strip()
                .lower()
            )
            continue

        if section != "position_restraints":
            continue

        fields = line.split()

        try:
            int(fields[0])
        except (
            IndexError,
            ValueError,
        ):
            continue

        count += 1

    return count


def run_grompp() -> tuple[int, int, str]:
    STATIC_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    gmx = shutil.which("gmx")

    if gmx is None:
        default = Path(
            "/usr/local/gromacs/bin/gmx"
        )

        if default.exists():
            gmx = str(default)

    if gmx is None:
        raise RuntimeError(
            "Could not locate the GROMACS executable"
        )

    command = [
        gmx,
        "grompp",
        "-f",
        str(TARGET_MDP),
        "-c",
        str(START_GRO),
        "-r",
        str(START_GRO),
        "-t",
        str(START_CPT),
        "-p",
        str(TOP),
        "-n",
        str(NDX),
        "-o",
        str(STATIC_TPR),
        "-po",
        str(STATIC_MDOUT),
        "-pp",
        str(STATIC_TOP),
        "-maxwarn",
        "0",
    ]

    completed = subprocess.run(
        command,
        cwd=TOPOLOGY_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    STATIC_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    restraint_count = (
        count_position_restraints(
            STATIC_TOP
        )
        if STATIC_TOP.exists()
        else -1
    )

    return (
        completed.returncode,
        restraint_count,
        completed.stdout,
    )


def write_outputs(
    return_code: int,
    restraint_count: int,
) -> None:
    outputs_present = all(
        path.exists()
        and path.stat().st_size > 0
        for path in (
            STATIC_TPR,
            STATIC_TOP,
            STATIC_MDOUT,
        )
    )

    grompp_pass = (
        return_code == 0
        and outputs_present
    )

    restraint_pass = (
        restraint_count == 0
    )

    overall_pass = (
        grompp_pass
        and restraint_pass
    )

    row = {
        "stage": STAGE07,
        "duration_ps": 25.0,
        "dt_ps": 0.0005,
        "nsteps": 50000,
        "coordinate_interval_ps": 0.5,
        "grompp_return_code": return_code,
        "grompp_pass": grompp_pass,
        "position_restraint_entries": (
            restraint_count
        ),
        "position_restraint_validation": (
            restraint_pass
        ),
        "preparation_decision": (
            "PASS"
            if overall_pass
            else "BLOCKED"
        ),
        "long_mobile_production_authorized": False,
        "mdp": relative(TARGET_MDP),
        "static_tpr": relative(STATIC_TPR),
        "processed_topology": (
            relative(STATIC_TOP)
        ),
        "grompp_log": relative(STATIC_LOG),
    }

    with PREPARATION_SUMMARY.open(
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

    PREPARATION_REPORT.write_text(
        f"""# Day022 Stage07 Preparation

- Previous checkpoint: `{STAGE06}`
- Stage: `{STAGE07}`
- Duration: 25 ps
- Time step: 0.0005 ps
- Steps: 50000
- Coordinate interval: 0.5 ps
- Static grompp: **{'PASS' if grompp_pass else 'FAIL'}**
- Active position restraints: **{restraint_count}**
- Preparation decision: **{row['preparation_decision']}**
- Long mobile production authorized: **NO**
""",
        encoding="utf-8",
    )

    if not overall_pass:
        raise RuntimeError(
            "Stage07 preparation is BLOCKED. "
            f"Review {STATIC_LOG}"
        )


def main() -> None:
    require_inputs()
    validate_stage06_authorization()
    create_stage07_mdp()

    return_code, restraint_count, _ = (
        run_grompp()
    )

    write_outputs(
        return_code,
        restraint_count,
    )

    print(
        "Day022 Stage07 preparation completed."
    )

    print(
        "Stage06 authorization: PASS_WITH_MONITORING"
    )

    print(
        "Stage07: 07_nvt_unrestrained_25ps"
    )

    print(
        "Stage07 duration / steps / dt: "
        "25 ps / 50000 / 0.0005 ps"
    )

    print(
        "Coordinate sampling interval: 0.5 ps"
    )

    print(
        "Stage07 static grompp: PASS"
    )

    print(
        f"Stage07 active position restraints: "
        f"{restraint_count}/0"
    )

    print(
        "Stage07 preparation decision: PASS"
    )

    print(
        "Long mobile production authorized: NO"
    )

    print(
        f"Wrote: {relative(PREPARATION_REPORT)}"
    )


if __name__ == "__main__":
    main()
