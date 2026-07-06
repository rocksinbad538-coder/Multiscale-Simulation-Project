#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PROTOCOL_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

INPUT_ROOT = PROTOCOL_ROOT / "protocol_inputs"
MDP_ROOT = INPUT_ROOT / "mdp"
TOPOLOGY_ROOT = INPUT_ROOT / "topology"

STAGE03 = "03_nvt_k1000_2ps"
STAGE03B = "03b_nvt_k1000_hold_2ps"

RUN03 = (
    PROTOCOL_ROOT
    / "execution"
    / STAGE03
)

DIAGNOSTIC_CSV = (
    RUN03
    / "stage03_release_response_diagnostic.csv"
)

SUMMARY03 = (
    RUN03
    / f"{STAGE03}_summary.csv"
)

REPORT03 = (
    RUN03
    / f"{STAGE03.upper()}_DAY021.md"
)

MANIFEST03 = (
    RUN03
    / f"{STAGE03}_output_manifest.csv"
)

SOURCE_MDP = (
    MDP_ROOT
    / f"{STAGE03}.mdp"
)

TARGET_MDP = (
    MDP_ROOT
    / f"{STAGE03B}.mdp"
)

STAGE03_GRO = (
    RUN03
    / f"{STAGE03}.gro"
)

STAGE03_CPT = (
    RUN03
    / f"{STAGE03}.cpt"
)

REFERENCE_GRO = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute/"
    "nvt_100ps_frozenSolute.gro"
)

TOP = (
    TOPOLOGY_ROOT
    / "hbn_pyrene_mobile_release.top"
)

NDX = (
    INPUT_ROOT
    / "mobile_release_index.ndx"
)

STATIC_ROOT = (
    PROTOCOL_ROOT
    / "static_validation"
    / STAGE03B
)

STATIC_TPR = (
    STATIC_ROOT
    / f"{STAGE03B}.tpr"
)

STATIC_PROCESSED_TOP = (
    STATIC_ROOT
    / f"{STAGE03B}_processed.top"
)

STATIC_MDOUT = (
    STATIC_ROOT
    / f"{STAGE03B}_mdout.mdp"
)

STATIC_LOG = (
    STATIC_ROOT
    / f"{STAGE03B}_grompp.log"
)

STATIC_CSV = (
    PROTOCOL_ROOT
    / "static_validation/"
    "mobile_release_static_validation.csv"
)

INPUT_MANIFEST = (
    INPUT_ROOT
    / "mobile_release_protocol_manifest.csv"
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


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def require_files() -> None:
    required = (
        DIAGNOSTIC_CSV,
        SUMMARY03,
        REPORT03,
        SOURCE_MDP,
        STAGE03_GRO,
        STAGE03_CPT,
        REFERENCE_GRO,
        TOP,
        NDX,
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
            "Missing or empty inputs:\n"
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
        fields = list(
            reader.fieldnames or []
        )

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


def verify_diagnostic() -> dict[str, str]:
    row, _ = read_single_row(
        DIAGNOSTIC_CSV
    )

    if (
        row.get("decision", "")
        .strip()
        .upper()
        != "REVISE"
    ):
        raise RuntimeError(
            "Stage03 diagnostic is not REVISE"
        )

    zero_fields = (
        "harmful_nonfinite_matches",
        "serious_instability_matches",
        "numeric_xvg_nonfinite_values",
    )

    failed = [
        field
        for field in zero_fields
        if int(
            float(
                row.get(field, "-1")
            )
        )
        != 0
    ]

    if failed:
        raise RuntimeError(
            "Stage03 contains blocking diagnostic "
            "signals:\n"
            + "\n".join(failed)
        )

    return row


def reconcile_stage03_summary() -> None:
    row, fields = read_single_row(
        SUMMARY03
    )

    row["decision"] = "REVISE"
    row["instability_signatures"] = ""
    row["blocked_reasons"] = ""

    row["reconciliation_status"] = "PASS"

    row["reconciliation_reason"] = (
        "epsilon-rf = inf is harmless PME "
        "metadata; local HBN displacement "
        "requires adaptive hold at k=1000"
    )

    row["rerun_required"] = "False"

    row["next_action"] = (
        "03b_nvt_k1000_hold_2ps"
    )

    write_single_row(
        SUMMARY03,
        row,
        fields,
    )

    report_text = REPORT03.read_text(
        encoding="utf-8",
        errors="replace",
    )

    report_text = report_text.replace(
        "- Decision: **BLOCKED**",
        "- Decision: **REVISE**",
        1,
    )

    report_text = report_text.replace(
        "- Instability signatures: "
        "['NONFINITE_VALUE:"
        "03_nvt_k1000_2ps.log:194']",
        "- Instability signatures: none",
        1,
    )

    section = """
## Stage03 reconciliation

The apparent non-finite value was the normal GROMACS
parameter record `epsilon-rf = inf`, not a runtime
observable. No LINCS, SETTLE, fatal, or numeric XVG
failures were present.

The stage remains **REVISE**, rather than PASS, because
five HBN atoms exceeded 0.08 nm and two exceeded 0.10 nm.
The authorized next action is an additional 2 ps hold at
the same restraint strength, k = 1000 kJ mol^-1 nm^-2.

- Stage03 rerun required: **NO**
- Stage04 authorized: **NO**
"""

    if "## Stage03 reconciliation" not in report_text:
        report_text = (
            report_text.rstrip()
            + "\n\n"
            + section.strip()
            + "\n"
        )

    REPORT03.write_text(
        report_text,
        encoding="utf-8",
    )


def rewrite_stage03_manifest() -> None:
    files = sorted(
        path
        for path in RUN03.iterdir()
        if (
            path.is_file()
            and path != MANIFEST03
            and not path.name.endswith(
                "_wrapper.log"
            )
        )
    )

    with MANIFEST03.open(
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

        for path in files:
            writer.writerow(
                {
                    "path": relative(path),
                    "size_bytes": (
                        path.stat().st_size
                    ),
                    "sha256": sha256(path),
                }
            )


def create_hold_mdp() -> None:
    text = SOURCE_MDP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    text = text.replace(
        "; Day021 mobile-release protocol",
        "; Day021 adaptive k=1000 hold stage",
        1,
    )

    text = (
        text.rstrip()
        + "\n\n"
        + "; Stage03b: additional 2 ps at unchanged "
        + "restraint strength after Stage03 REVISE\n"
    )

    required_records = (
        "integrator                   = md",
        "dt                           = 0.00025",
        "nsteps                       = 8000",
        "define                       = "
        "-DPOSRES_HBN_K1000 "
        "-DPOSRES_PYR_K1000",
        "gen-vel                      = no",
        "continuation                 = yes",
    )

    missing = [
        record
        for record in required_records
        if record not in text
    ]

    if missing:
        raise RuntimeError(
            "Source MDP lacks expected records:\n"
            + "\n".join(missing)
        )

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


def run_static_grompp() -> dict[str, object]:
    STATIC_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for path in (
        STATIC_TPR,
        STATIC_PROCESSED_TOP,
        STATIC_MDOUT,
    ):
        if path.exists():
            path.unlink()

    gmx = (
        Path("/usr/local/gromacs/bin/gmx")
        if Path(
            "/usr/local/gromacs/bin/gmx"
        ).exists()
        else Path(
            shutil.which("gmx") or ""
        )
    )

    if not gmx.exists():
        raise RuntimeError(
            "Could not locate gmx"
        )

    command = [
        str(gmx),
        "grompp",
        "-f",
        str(TARGET_MDP),
        "-c",
        str(STAGE03_GRO),
        "-r",
        str(REFERENCE_GRO),
        "-t",
        str(STAGE03_CPT),
        "-p",
        str(TOP),
        "-n",
        str(NDX),
        "-o",
        str(STATIC_TPR),
        "-po",
        str(STATIC_MDOUT),
        "-pp",
        str(STATIC_PROCESSED_TOP),
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

    outputs_exist = all(
        path.exists()
        and path.stat().st_size > 0
        for path in (
            STATIC_TPR,
            STATIC_PROCESSED_TOP,
            STATIC_MDOUT,
        )
    )

    restraint_count = (
        count_position_restraints(
            STATIC_PROCESSED_TOP
        )
        if STATIC_PROCESSED_TOP.exists()
        else -1
    )

    return {
        "stage": STAGE03B,
        "kind": "nvt",
        "restraint_k": 1000,
        "grompp_return_code": (
            result.returncode
        ),
        "grompp_pass": (
            result.returncode == 0
            and outputs_exist
        ),
        "position_restraint_entries": (
            restraint_count
        ),
        "expected_position_restraint_entries": (
            1706
        ),
        "restraint_entry_validation": (
            restraint_count == 1706
        ),
        "tpr_path": relative(
            STATIC_TPR
        ),
        "processed_topology_path": (
            relative(
                STATIC_PROCESSED_TOP
            )
        ),
        "grompp_log": relative(
            STATIC_LOG
        ),
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
        fields = list(
            reader.fieldnames or []
        )

    rows = [
        row
        for row in rows
        if row.get("stage") != STAGE03B
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


def update_input_manifest() -> None:
    with INPUT_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    verified = []

    for row in rows:
        path = PROJECT_ROOT / row["path"]

        if not path.exists():
            raise RuntimeError(
                f"Manifest file missing: {path}"
            )

        if sha256(path) != row["sha256"]:
            raise RuntimeError(
                "Existing protocol-input hash "
                f"mismatch: {path}"
            )

        verified.append(row)

    verified = [
        row
        for row in verified
        if row["path"] != relative(
            TARGET_MDP
        )
    ]

    verified.append(
        {
            "path": relative(
                TARGET_MDP
            ),
            "size_bytes": (
                TARGET_MDP.stat().st_size
            ),
            "sha256": sha256(
                TARGET_MDP
            ),
        }
    )

    verified.sort(
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
        writer.writerows(verified)


def main() -> None:
    require_files()
    diagnostic = verify_diagnostic()

    reconcile_stage03_summary()
    rewrite_stage03_manifest()
    create_hold_mdp()

    static_row = run_static_grompp()

    if not static_row["grompp_pass"]:
        raise RuntimeError(
            "Stage03b static grompp failed. "
            f"See {STATIC_LOG}"
        )

    if not static_row[
        "restraint_entry_validation"
    ]:
        raise RuntimeError(
            "Stage03b restraint count failed"
        )

    update_static_csv(
        static_row
    )

    update_input_manifest()

    print(
        "Day021 Stage03b adaptive hold "
        "preparation completed."
    )

    print(
        "Stage03 scientific status: REVISE"
    )

    print(
        "Stage03 numerical instability: NO"
    )

    print(
        "Stage03 rerun required: NO"
    )

    print(
        f"Stage03 HBN atoms above "
        f"0.08/0.10 nm: "
        f"{diagnostic['HBN_atoms_above_0p08_nm']}/"
        f"{diagnostic['HBN_atoms_above_0p10_nm']}"
    )

    print(
        f"Stage03b MDP: "
        f"{relative(TARGET_MDP)}"
    )

    print(
        f"Stage03b static grompp: "
        f"{'PASS' if static_row['grompp_pass'] else 'FAIL'}"
    )

    print(
        "Stage03b position restraints: "
        f"{static_row['position_restraint_entries']}/1706"
    )

    print(
        "Stage04 authorized: NO"
    )


if __name__ == "__main__":
    main()
