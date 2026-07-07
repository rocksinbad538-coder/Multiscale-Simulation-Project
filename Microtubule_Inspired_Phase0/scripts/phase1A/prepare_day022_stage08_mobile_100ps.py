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

STAGE07 = "07_nvt_unrestrained_25ps"
STAGE08 = "08_nvt_mobile_100ps"

INPUT_ROOT = PROTOCOL / "protocol_inputs"
MDP_ROOT = INPUT_ROOT / "mdp"
TOPOLOGY_ROOT = INPUT_ROOT / "topology"

RUN07 = PROTOCOL / "execution" / STAGE07

SOURCE_MDP = MDP_ROOT / f"{STAGE07}.mdp"
TARGET_MDP = MDP_ROOT / f"{STAGE08}.mdp"

START_GRO = RUN07 / f"{STAGE07}.gro"
START_CPT = RUN07 / f"{STAGE07}.cpt"

TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    INPUT_ROOT
    / "mobile_release_index.ndx"
)

OPERATIONAL_SUMMARY = (
    RUN07
    / f"{STAGE07}_summary.csv"
)

STRUCTURAL_SUMMARY = (
    RUN07
    / f"{STAGE07}_structural_summary.csv"
)

TIME_RESOLVED_SUMMARY = (
    RUN07
    / "time_resolved_stability/"
    "stage07_time_resolved_stability_summary.csv"
)

STATIC_ROOT = (
    PROTOCOL
    / "static_validation"
    / STAGE08
)

STATIC_TPR = STATIC_ROOT / f"{STAGE08}.tpr"
STATIC_TOP = STATIC_ROOT / f"{STAGE08}_processed.top"
STATIC_MDOUT = STATIC_ROOT / f"{STAGE08}_mdout.mdp"
STATIC_LOG = STATIC_ROOT / f"{STAGE08}_grompp.log"

MASTER_STATIC_CSV = (
    PROTOCOL
    / "static_validation/"
    "mobile_release_static_validation.csv"
)

PREPARATION_SUMMARY = (
    STATIC_ROOT
    / f"{STAGE08}_preparation_summary.csv"
)

AUTHORIZATION_REPORT = (
    STATIC_ROOT
    / "STAGE08_100PS_MOBILE_PRODUCTION_AUTHORIZATION_DAY022.md"
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
            f"Expected one row in {path}, found {len(rows)}"
        )

    return rows[0]


def is_true(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def require_inputs() -> None:
    required = (
        SOURCE_MDP,
        START_GRO,
        START_CPT,
        TOP,
        NDX,
        OPERATIONAL_SUMMARY,
        STRUCTURAL_SUMMARY,
        TIME_RESOLVED_SUMMARY,
        MASTER_STATIC_CSV,
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


def validate_stage07() -> None:
    operational = read_single_row(
        OPERATIONAL_SUMMARY
    )

    structural = read_single_row(
        STRUCTURAL_SUMMARY
    )

    temporal = read_single_row(
        TIME_RESOLVED_SUMMARY
    )

    failures = []

    if (
        operational.get("decision", "")
        .strip()
        .upper()
        != "PASS"
    ):
        failures.append(
            "Stage07 operational decision is not PASS"
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
            "Stage07 structural screen is not STABLE_CANDIDATE"
        )

    if (
        temporal.get(
            "time_resolved_decision",
            "",
        )
        .strip()
        .upper()
        != "PASS_WITH_MONITORING"
    ):
        failures.append(
            "Stage07 temporal decision is not PASS_WITH_MONITORING"
        )

    if (
        temporal.get(
            "authorized_next_step",
            "",
        )
        .strip()
        != "100PS_MOBILE_PRODUCTION_CANDIDATE"
    ):
        failures.append(
            "The 100 ps mobile production is not authorized"
        )

    if temporal.get(
        "blocked_reasons",
        "",
    ).strip():
        failures.append(
            "Stage07 has blocking reasons"
        )

    quantitative_checks = {
        "no impropers persistent above 20 degrees in >=80% frames": (
            int(
                float(
                    temporal.get(
                        "HBN_impropers_above20_persistent80_count",
                        "-1",
                    )
                )
            )
            == 0
        ),
        "minimum contact >= 0.14 nm": (
            float(
                temporal.get(
                    "minimum_intergroup_contact_nm",
                    "0",
                )
            )
            >= 0.14
        ),
        "PYR maximum aligned RMS <= 0.08 nm": (
            float(
                temporal.get(
                    "PYR_max_aligned_rms_nm",
                    "999",
                )
            )
            <= 0.08
        ),
    }

    failures.extend(
        label
        for label, passed
        in quantitative_checks.items()
        if not passed
    )

    if failures:
        raise RuntimeError(
            "Stage07 production authorization failed:\n"
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


def create_mdp() -> None:
    text = SOURCE_MDP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    settings = {
        "nsteps": "200000",
        "dt": "0.0005",
        "continuation": "yes",
        "gen-vel": "no",
        "nstxout": "0",
        "nstvout": "0",
        "nstfout": "0",
        "nstxout-compressed": "1000",
        "nstcalcenergy": "100",
        "nstenergy": "200",
        "nstlog": "1000",
    }

    for key, value in settings.items():
        text = set_mdp_value(
            text,
            key,
            value,
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

    required = {
        "integrator": "md",
        "nsteps": "200000",
        "dt": "0.0005",
        "continuation": "yes",
        "gen-vel": "no",
        "nstxout-compressed": "1000",
    }

    failures = []

    for key, expected in required.items():
        actual = parsed.get(
            key,
            "",
        )

        if actual != expected:
            failures.append(
                f"{key}: expected {expected}, "
                f"found {actual or 'MISSING'}"
            )

    if "posres" in parsed.get(
        "define",
        "",
    ):
        failures.append(
            "Position restraints remain active"
        )

    if failures:
        raise RuntimeError(
            "Stage08 MDP contract failed:\n"
            + "\n".join(failures)
        )

    text = (
        text.rstrip()
        + "\n\n"
        + "; Day022 mobile NVT production at 300 K\n"
        + "; 200000 steps x 0.0005 ps = 100 ps\n"
        + "; compressed coordinates every 0.5 ps\n"
        + "; expected trajectory frames including endpoints: 201\n"
    )

    if TARGET_MDP.exists():
        existing = TARGET_MDP.read_text(
            encoding="utf-8",
            errors="replace",
        )

        if existing != text:
            raise RuntimeError(
                "Stage08 MDP already exists with different content"
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


def run_grompp() -> tuple[int, int]:
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
            "Could not locate GROMACS"
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
    )


def register_static_validation(
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

    if (
        return_code != 0
        or not outputs_present
        or restraint_count != 0
    ):
        raise RuntimeError(
            "Stage08 static validation failed. "
            f"Review {STATIC_LOG}"
        )

    with MASTER_STATIC_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    new_row = {
        "stage": STAGE08,
        "kind": "nvt",
        "restraint_k": "0",
        "grompp_return_code": "0",
        "grompp_pass": "True",
        "position_restraint_entries": "0",
        "expected_position_restraint_entries": "0",
        "restraint_entry_validation": "True",
        "tpr_path": relative(STATIC_TPR),
        "processed_topology_path": relative(
            STATIC_TOP
        ),
        "grompp_log": relative(STATIC_LOG),
    }

    for key in new_row:
        if key not in fields:
            fields.append(key)

    rows = [
        row
        for row in rows
        if row.get(
            "stage",
            "",
        ).strip() != STAGE08
    ]

    rows.append(new_row)

    temporary = MASTER_STATIC_CSV.with_suffix(
        ".csv.tmp"
    )

    with temporary.open(
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
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )

    temporary.replace(
        MASTER_STATIC_CSV
    )

    with MASTER_STATIC_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if row.get(
                "stage",
                "",
            ).strip() == STAGE08
        ]

    if len(matches) != 1:
        raise RuntimeError(
            "Stage08 was not uniquely registered"
        )


def write_outputs(
    return_code: int,
    restraint_count: int,
) -> None:
    row = {
        "stage": STAGE08,
        "duration_ps": 100.0,
        "dt_ps": 0.0005,
        "nsteps": 200000,
        "coordinate_interval_ps": 0.5,
        "expected_xtc_frames": 201,
        "continued_from": STAGE07,
        "grompp_return_code": return_code,
        "grompp_pass": (
            return_code == 0
        ),
        "active_position_restraints": (
            restraint_count
        ),
        "preparation_decision": "PASS",
        "mobile_100ps_authorized": True,
        "longer_mobile_production_authorized": False,
        "multitemperature_production_authorized": False,
        "mdp": relative(TARGET_MDP),
        "tpr": relative(STATIC_TPR),
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

    AUTHORIZATION_REPORT.write_text(
        f"""# Day022 Stage08 Mobile Production Authorization

- Previous validated checkpoint: `{STAGE07}`
- Stage07 time-resolved decision: **PASS_WITH_MONITORING**
- Authorized stage: `{STAGE08}`
- Duration: **100 ps**
- Temperature ensemble: **NVT at 300 K**
- Time step: **0.0005 ps**
- Coordinate interval: **0.5 ps**
- Expected trajectory frames: **201**
- Active position restraints: **0**
- Static grompp: **PASS**
- Stage08 execution authorized: **YES**
- Production longer than 100 ps authorized: **NO**
- Multitemperature production authorized: **NO**
""",
        encoding="utf-8",
    )


def main() -> None:
    require_inputs()
    validate_stage07()
    create_mdp()

    return_code, restraint_count = run_grompp()

    register_static_validation(
        return_code,
        restraint_count,
    )

    write_outputs(
        return_code,
        restraint_count,
    )

    print(
        "Day022 Stage08 mobile-production preparation completed."
    )
    print(
        "Stage07 time-resolved decision: PASS_WITH_MONITORING"
    )
    print(
        "Stage08: 08_nvt_mobile_100ps"
    )
    print(
        "Duration / steps / dt: "
        "100 ps / 200000 / 0.0005 ps"
    )
    print(
        "Coordinate sampling interval: 0.5 ps"
    )
    print(
        "Expected compressed trajectory frames: 201"
    )
    print(
        "Stage08 static grompp: PASS"
    )
    print(
        f"Stage08 active position restraints: "
        f"{restraint_count}/0"
    )
    print(
        "Stage08 master static-registration rows: 1"
    )
    print(
        "Stage08 execution authorization: PASS"
    )
    print(
        "Production longer than 100 ps authorized: NO"
    )
    print(
        "Multitemperature production authorized: NO"
    )
    print(
        f"Wrote: {relative(AUTHORIZATION_REPORT)}"
    )


if __name__ == "__main__":
    main()
