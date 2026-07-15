#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
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

INPUT_ROOT = PROTOCOL / "protocol_inputs"
MDP_ROOT = INPUT_ROOT / "mdp"
TOPOLOGY_ROOT = INPUT_ROOT / "topology"

STAGE05 = "05_nvt_unrestrained_2ps"
STAGE06 = "06_nvt_unrestrained_10ps"

RUN05 = (
    PROTOCOL
    / "execution"
    / STAGE05
)

SOURCE_MDP = (
    MDP_ROOT
    / f"{STAGE05}.mdp"
)

TARGET_MDP = (
    MDP_ROOT
    / f"{STAGE06}.mdp"
)

START_GRO = (
    RUN05
    / f"{STAGE05}.gro"
)

START_CPT = (
    RUN05
    / f"{STAGE05}.cpt"
)

TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    INPUT_ROOT
    / "mobile_release_index.ndx"
)

OPERATIONAL_SUMMARY = (
    RUN05
    / f"{STAGE05}_summary.csv"
)

HBN_STRUCTURAL_SUMMARY = (
    RUN05
    / f"{STAGE05}_structural_summary.csv"
)

INTEGRATED_SUMMARY = (
    RUN05
    / "stage05_integrated_mobile_pilot_summary.csv"
)

IMPROPER_SUMMARY = (
    RUN05
    / "stage05_hbn_improper_phase_diagnostic.csv"
)

AUTHORIZATION_REPORT = (
    RUN05
    / "STAGE05_SHORT_MOBILE_VALIDATION_AUTHORIZATION_DAY021.md"
)

STATIC_CSV = (
    PROTOCOL
    / "static_validation/mobile_release_static_validation.csv"
)

INPUT_MANIFEST = (
    INPUT_ROOT
    / "mobile_release_protocol_manifest.csv"
)

STATIC_ROOT = (
    PROTOCOL
    / "static_validation"
    / STAGE06
)

STATIC_TPR = (
    STATIC_ROOT
    / f"{STAGE06}.tpr"
)

STATIC_TOP = (
    STATIC_ROOT
    / f"{STAGE06}_processed.top"
)

STATIC_MDOUT = (
    STATIC_ROOT
    / f"{STAGE06}_mdout.mdp"
)

STATIC_LOG = (
    STATIC_ROOT
    / f"{STAGE06}_grompp.log"
)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(ROOT)
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


def require_files() -> None:
    required = (
        SOURCE_MDP,
        START_GRO,
        START_CPT,
        TOP,
        NDX,
        OPERATIONAL_SUMMARY,
        HBN_STRUCTURAL_SUMMARY,
        INTEGRATED_SUMMARY,
        IMPROPER_SUMMARY,
        STATIC_CSV,
        INPUT_MANIFEST,
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


def read_single_row(
    path: Path,
) -> tuple[dict[str, str], list[str]]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}"
        )

    return rows[0], fields


def write_single_row(
    path: Path,
    row: dict[str, object],
    fields: list[str],
) -> None:
    for key in row:
        if key not in fields:
            fields.append(key)

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

        writer.writerow(
            {
                field: row.get(field, "")
                for field in fields
            }
        )


def number(
    row: dict[str, str],
    field: str,
) -> float:
    try:
        return float(row[field])
    except (KeyError, ValueError) as error:
        raise RuntimeError(
            f"Invalid field {field}"
        ) from error


def validate_stage05() -> None:
    operational, _ = read_single_row(
        OPERATIONAL_SUMMARY
    )

    structural, _ = read_single_row(
        HBN_STRUCTURAL_SUMMARY
    )

    integrated, _ = read_single_row(
        INTEGRATED_SUMMARY
    )

    improper, _ = read_single_row(
        IMPROPER_SUMMARY
    )

    failures = []

    if (
        operational.get("decision", "")
        .strip()
        .upper()
        != "PASS"
    ):
        failures.append(
            "Stage05 operational decision is not PASS"
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
            "Stage05 HBN structural screen is not stable"
        )

    if integrated.get(
        "blocked_reasons",
        "",
    ).strip():
        failures.append(
            "Stage05 integrated validation has blocking reasons"
        )

    if (
        improper.get(
            "targeted_decision",
            "",
        )
        .strip()
        .upper()
        != "PASS"
    ):
        failures.append(
            "HBN improper-phase diagnostic is not PASS"
        )

    quantitative_checks = {
        "planarity q99 <= 20 deg": (
            number(
                improper,
                "current_planarity_q99_deg",
            )
            <= 20.0
        ),
        "planarity maximum <= 40 deg": (
            number(
                improper,
                "current_planarity_max_deg",
            )
            <= 40.0
        ),
        "calibrated equilibrium q99 <= 20 deg": (
            number(
                improper,
                "calibrated_equilibrium_q99_deg",
            )
            <= 20.0
        ),
        "calibrated equilibrium maximum <= 40 deg": (
            number(
                improper,
                "calibrated_equilibrium_max_deg",
            )
            <= 40.0
        ),
        "stage-change maximum <= 60 deg": (
            number(
                improper,
                "stage_change_max_deg",
            )
            <= 60.0
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
            "Stage05 short-validation authorization failed:\n"
            + "\n".join(failures)
        )


def parse_mdp(text: str) -> dict[str, str]:
    values = {}

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

        values[
            key.strip().lower()
        ] = value.strip()

    return values


def create_stage06_mdp() -> None:
    source_text = SOURCE_MDP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    values = parse_mdp(
        source_text
    )

    expected = {
        "integrator": "md",
        "nsteps": "4000",
        "dt": "0.0005",
        "continuation": "yes",
        "gen-vel": "no",
        "comm-mode": "linear",
        "comm-grps": "system",
    }

    failures = []

    for key, expected_value in expected.items():
        actual = values.get(
            key,
            "",
        ).lower()

        if actual != expected_value:
            failures.append(
                f"{key}: expected {expected_value}, "
                f"found {actual or 'MISSING'}"
            )

    if "posres" in source_text.lower():
        failures.append(
            "Stage05 MDP unexpectedly contains POSRES"
        )

    if failures:
        raise RuntimeError(
            "Stage05 MDP contract failed:\n"
            + "\n".join(failures)
        )

    target_text, replacements = re.subn(
        r"(?mi)^(\s*nsteps\s*=\s*)\S+.*$",
        r"\g<1>20000",
        source_text,
        count=1,
    )

    if replacements != 1:
        raise RuntimeError(
            "Could not replace nsteps uniquely"
        )

    target_text = (
        target_text.rstrip()
        + "\n\n"
        + "; Day021 short unrestrained mobile validation\n"
        + "; 20000 steps x 0.0005 ps = 10 ps\n"
    )

    if TARGET_MDP.exists():
        existing = TARGET_MDP.read_text(
            encoding="utf-8",
            errors="replace",
        )

        if existing != target_text:
            raise RuntimeError(
                "Stage06 MDP already exists with different content"
            )
    else:
        TARGET_MDP.write_text(
            target_text,
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


def run_static_grompp() -> dict[str, object]:
    STATIC_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    gmx_path = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if not gmx_path.exists():
        discovered = shutil.which("gmx")

        if discovered is None:
            raise RuntimeError(
                "Could not locate gmx"
            )

        gmx_path = Path(discovered)

    command = [
        str(gmx_path),
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

    result = subprocess.run(
        command,
        cwd=TOPOLOGY_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    STATIC_LOG.write_text(
        result.stdout,
        encoding="utf-8",
    )

    outputs_present = all(
        path.exists()
        and path.stat().st_size > 0
        for path in (
            STATIC_TPR,
            STATIC_MDOUT,
            STATIC_TOP,
        )
    )

    restraint_count = (
        count_position_restraints(
            STATIC_TOP
        )
        if STATIC_TOP.exists()
        else -1
    )

    return {
        "stage": STAGE06,
        "kind": "nvt",
        "restraint_k": 0,
        "grompp_return_code": result.returncode,
        "grompp_pass": (
            result.returncode == 0
            and outputs_present
        ),
        "position_restraint_entries": (
            restraint_count
        ),
        "expected_position_restraint_entries": 0,
        "restraint_entry_validation": (
            restraint_count == 0
        ),
        "tpr_path": relative(STATIC_TPR),
        "processed_topology_path": (
            relative(STATIC_TOP)
        ),
        "grompp_log": relative(STATIC_LOG),
    }


def update_static_csv(
    new_row: dict[str, object],
) -> None:
    with STATIC_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    rows = [
        row
        for row in rows
        if row.get("stage") != STAGE06
    ]

    rows.append(
        {
            key: str(value)
            for key, value in new_row.items()
        }
    )

    for key in new_row:
        if key not in fields:
            fields.append(key)

    with STATIC_CSV.open(
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


def update_manifest() -> int:
    with INPUT_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    for row in rows:
        path = ROOT / row["path"]

        if not path.exists():
            raise RuntimeError(
                f"Manifest file missing: {path}"
            )

        if sha256(path) != row["sha256"]:
            raise RuntimeError(
                f"Manifest hash mismatch: {path}"
            )

    rows = [
        row
        for row in rows
        if row["path"] != relative(TARGET_MDP)
    ]

    rows.append(
        {
            "path": relative(TARGET_MDP),
            "size_bytes": TARGET_MDP.stat().st_size,
            "sha256": sha256(TARGET_MDP),
        }
    )

    rows.sort(
        key=lambda row: row["path"]
    )

    with INPUT_MANIFEST.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "path",
                "size_bytes",
                "sha256",
            ],
        )

        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


def reconcile_stage05() -> None:
    integrated, integrated_fields = (
        read_single_row(
            INTEGRATED_SUMMARY
        )
    )

    integrated.update(
        {
            "pilot_readiness": (
                "READY_FOR_SHORT_MOBILE_VALIDATION"
            ),
            "authorized_next_step": STAGE06,
            "long_mobile_production_authorized": "False",
            "review_reasons": "",
            "blocked_reasons": "",
            "improper_phase_reconciliation_status": "PASS",
            "improper_phase_reconciliation_reason": (
                "The approximately 180-degree mismatch was "
                "a dihedral phase-convention artifact; "
                "calibrated planarity and equilibrium "
                "deviations passed targeted thresholds"
            ),
        }
    )

    write_single_row(
        INTEGRATED_SUMMARY,
        integrated,
        integrated_fields,
    )

    operational, operational_fields = (
        read_single_row(
            OPERATIONAL_SUMMARY
        )
    )

    operational.update(
        {
            "integrated_validation_status": "PASS",
            "short_mobile_validation_authorized": "True",
            "long_mobile_production_authorized": "False",
            "next_action": STAGE06,
        }
    )

    write_single_row(
        OPERATIONAL_SUMMARY,
        operational,
        operational_fields,
    )

    AUTHORIZATION_REPORT.write_text(
        f"""# Day021 Stage05 Short Mobile Validation Authorization

- Stage05 operational decision: **PASS**
- HBN structural screen: **STABLE_CANDIDATE**
- HBN improper-phase targeted decision: **PASS**
- Integrated pilot decision: **READY_FOR_SHORT_MOBILE_VALIDATION**
- Authorized next stage: `{STAGE06}`
- Long mobile production authorized: **NO**

The previous approximately 180-degree equilibrium mismatch was
caused by a dihedral phase-convention artifact. Calibrated
planarity and equilibrium-deviation metrics passed the targeted
acceptance criteria.
""",
        encoding="utf-8",
    )


def main() -> None:
    require_files()
    validate_stage05()
    create_stage06_mdp()

    static_result = run_static_grompp()

    if not static_result["grompp_pass"]:
        raise RuntimeError(
            f"Stage06 static grompp failed. See {STATIC_LOG}"
        )

    if not static_result[
        "restraint_entry_validation"
    ]:
        raise RuntimeError(
            "Stage06 active-restraint count failed"
        )

    update_static_csv(
        static_result
    )

    manifest_count = update_manifest()

    reconcile_stage05()

    print(
        "Day021 Stage06 short-validation preparation completed."
    )
    print(
        "Stage05 scientific decision: PASS"
    )
    print(
        "Stage05 readiness: READY_FOR_SHORT_MOBILE_VALIDATION"
    )
    print(
        f"Stage06 MDP: {relative(TARGET_MDP)}"
    )
    print(
        "Stage06 duration: 10 ps"
    )
    print(
        "Stage06 static grompp: PASS"
    )
    print(
        "Stage06 active position restraints: 0/0"
    )
    print(
        f"Protocol-manifest files: {manifest_count}"
    )
    print(
        "Long mobile production authorized: NO"
    )
    print(
        f"Wrote: {relative(AUTHORIZATION_REPORT)}"
    )


if __name__ == "__main__":
    main()
