#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CONSTRUCTION_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction"
)

XYZ_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_start.xyz"
)

MAP_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

CONSTRUCTION_REPORT = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_CONSTRUCTION_REPORT.json"
)

PRE_QM_REPORT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit/"
    "QM_F06_UPPER_V5B_PRE_QM_AUDIT.json"
)

V4_CONSTRAINT_MAP = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_constraint_design/"
    "QM_F06_UPPER_V4_constraint_map.csv"
)

ORCA = Path(
    "/Users/alejandro/projects/"
    "orca_6_1_1_macosx_intel_openmpi411/orca"
)

PREPARATION_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_orca_input"
)

EXECUTION_ROOT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_executions"
)

EXPECTED_COUNT = 52

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 23,
})

RELEASED_FROM_V4_FIXED = {
    "A:UPPER:14:4",
}

NPROCS = 4
MAXCORE_MB = 2500
SCF_MAXITER = 500
GEOM_MAXITER = 250


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path):
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def read_xyz(path: Path):
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0])
    atoms = []

    for index, line in enumerate(lines[2:2 + count]):
        fields = line.split()

        atoms.append({
            "index": index,
            "element": fields[0],
            "x": float(fields[1]),
            "y": float(fields[2]),
            "z": float(fields[3]),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"XYZ incomplete: expected {count}, "
            f"found {len(atoms)}"
        )

    return atoms


def main() -> None:
    for path in (
        XYZ_PATH,
        MAP_PATH,
        CONSTRUCTION_REPORT,
        PRE_QM_REPORT,
        V4_CONSTRAINT_MAP,
    ):
        require_file(path)

    if not ORCA.is_file() or not os.access(ORCA, os.X_OK):
        raise RuntimeError(
            f"ORCA missing or non-executable: {ORCA}"
        )

    construction = json.loads(
        CONSTRUCTION_REPORT.read_text(
            encoding="utf-8"
        )
    )

    pre_qm = json.loads(
        PRE_QM_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not construction["overall_pass"]:
        raise RuntimeError(
            "Construction report did not pass."
        )

    if not pre_qm["authorization"][
        "constraint_design_authorized"
    ]:
        raise RuntimeError(
            "Constraint design is not authorized."
        )

    atoms = read_xyz(XYZ_PATH)
    provenance_rows = read_csv(MAP_PATH)
    v4_rows = read_csv(V4_CONSTRAINT_MAP)

    if len(atoms) != EXPECTED_COUNT:
        raise RuntimeError(
            f"Unexpected atom count: {len(atoms)}"
        )

    if len(provenance_rows) != EXPECTED_COUNT:
        raise RuntimeError(
            "Provenance-map row count mismatch."
        )

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    if composition != EXPECTED_COMPOSITION:
        raise RuntimeError(
            f"Composition mismatch: {composition}"
        )

    v4_fixed_ids = {
        row["atom_id"]
        for row in v4_rows
        if parse_bool(row["v4_fixed"])
    }

    fixed_ids = (
        v4_fixed_ids
        - RELEASED_FROM_V4_FIXED
    )

    map_ids = {
        row["atom_id"]
        for row in provenance_rows
    }

    missing_fixed = fixed_ids - map_ids

    if missing_fixed:
        raise RuntimeError(
            f"Fixed atoms missing from V5-B: "
            f"{sorted(missing_fixed)}"
        )

    constraint_rows = []

    for atom, row in zip(
        atoms,
        provenance_rows,
        strict=True,
    ):
        index = int(row["index_0based"])

        if index != atom["index"]:
            raise RuntimeError(
                f"Index mismatch at {index}"
            )

        atom_id = row["atom_id"]
        fixed = atom_id in fixed_ids
        mobile = not fixed

        if atom_id == "A:UPPER:14:4" and fixed:
            raise RuntimeError(
                "A:UPPER:14:4 must be mobile in V5-B."
            )

        if parse_bool(row["artificial_cap"]) and fixed:
            raise RuntimeError(
                f"Artificial cap cannot be fixed: {atom_id}"
            )

        if row["atom_role"] in {
            "RESTORED_CANONICAL_R2_ATOM",
            "RESTORED_CANONICAL_PASSIVANT_H",
        } and fixed:
            raise RuntimeError(
                f"Restored atom cannot be fixed: {atom_id}"
            )

        basis = (
            "RETAINED_VALIDATED_CORE"
            if fixed
            else (
                "E2915_FIXED_MOBILE_ASYMMETRY_RELEASE"
                if atom_id == "A:UPPER:14:4"
                else "V5B_REPAIRED_OR_BOUNDARY_MOBILE_REGION"
            )
        )

        constraint_rows.append({
            **row,
            "v5b_fixed": fixed,
            "v5b_mobile": mobile,
            "v5b_mobility_basis": basis,
        })

    fixed_indices = [
        int(row["index_0based"])
        for row in constraint_rows
        if row["v5b_fixed"]
    ]

    mobile_indices = [
        int(row["index_0based"])
        for row in constraint_rows
        if row["v5b_mobile"]
    ]

    if set(fixed_indices) & set(mobile_indices):
        raise RuntimeError(
            "Fixed/mobile overlap detected."
        )

    if (
        set(fixed_indices) | set(mobile_indices)
        != set(range(EXPECTED_COUNT))
    ):
        raise RuntimeError(
            "Fixed/mobile partition is incomplete."
        )

    # Neutral singlet electron-parity gate.
    electron_count = (
        16 * 5
        + 13 * 7
        + 23 * 1
    )

    if electron_count != 194:
        raise RuntimeError(
            f"Unexpected electron count: {electron_count}"
        )

    if electron_count % 2 != 0:
        raise RuntimeError(
            "Neutral singlet has odd electron count."
        )

    PREPARATION_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    constraint_map = (
        PREPARATION_DIR
        / "QM_F06_UPPER_V5B_constraint_map.csv"
    )

    fieldnames = list(constraint_rows[0])

    with constraint_map.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(constraint_rows)

    prepared_xyz = (
        PREPARATION_DIR / "v5b_start.xyz"
    )

    shutil.copy2(
        XYZ_PATH,
        prepared_xyz,
    )

    input_path = (
        PREPARATION_DIR / "v5b.inp"
    )

    lines = [
        (
            "! PBE0 D4 def2-TZVP def2/J "
            "RIJCOSX TightSCF Opt DefGrid3"
        ),
        "",
        "%pal",
        f"  nprocs {NPROCS}",
        "end",
        "",
        f"%maxcore {MAXCORE_MB}",
        "%scf",
        f"  MaxIter {SCF_MAXITER}",
        "end",
        "",
        "%geom",
        f"  MaxIter {GEOM_MAXITER}",
        "  Constraints",
    ]

    for index in fixed_indices:
        lines.append(f"    {{ C {index} C }}")

    lines.extend([
        "  end",
        "end",
        "",
        (
            "# QM_F06 UPPER V5-B selective boundary "
            "expansion and chemical repair"
        ),
        (
            "# Fresh SCF; A:UPPER:14:4 released; "
            "all restored atoms and passivants mobile"
        ),
        "* xyz 0 1",
    ])

    for atom in atoms:
        lines.append(
            f"{atom['element']:2s} "
            f"{atom['x']: .12f} "
            f"{atom['y']: .12f} "
            f"{atom['z']: .12f}"
        )

    lines.append("*")
    lines.append("")

    input_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    input_text = input_path.read_text(
        encoding="utf-8"
    )

    forbidden_reuse = {
        "moread": bool(re.search(
            r"\bmoread\b",
            input_text,
            re.IGNORECASE,
        )),
        "moinp": bool(re.search(
            r"%moinp\b",
            input_text,
            re.IGNORECASE,
        )),
        "gbw_reference": bool(re.search(
            r"\S+\.gbw",
            input_text,
            re.IGNORECASE,
        )),
    }

    gates = {
        "construction_pass": construction["overall_pass"],
        "pre_qm_pass": all(pre_qm["gates"].values()),
        "atom_count": len(atoms) == EXPECTED_COUNT,
        "composition": composition == EXPECTED_COMPOSITION,
        "electron_count": electron_count == 194,
        "electron_parity": electron_count % 2 == 0,
        "fixed_region_nonempty": len(fixed_indices) > 0,
        "mobile_region_nonempty": len(mobile_indices) > 0,
        "A14_4_mobile": (
            next(
                row
                for row in constraint_rows
                if row["atom_id"] == "A:UPPER:14:4"
            )["v5b_mobile"]
        ),
        "all_artificial_caps_mobile": all(
            row["v5b_mobile"]
            for row in constraint_rows
            if parse_bool(row["artificial_cap"])
        ),
        "fresh_scf": not any(
            forbidden_reuse.values()
        ),
        "orca_available": (
            ORCA.is_file()
            and os.access(ORCA, os.X_OK)
        ),
    }

    if not all(gates.values()):
        raise RuntimeError(
            "Execution preflight failed: "
            + json.dumps(gates, indent=2)
        )

    EXECUTION_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    execution_dir = (
        EXECUTION_ROOT
        / f"v5b_{timestamp}"
    )

    execution_dir.mkdir(
        parents=False,
        exist_ok=False,
    )

    execution_input = (
        execution_dir / "v5b.inp"
    )

    execution_xyz = (
        execution_dir / "v5b_start.xyz"
    )

    execution_map = (
        execution_dir
        / "QM_F06_UPPER_V5B_constraint_map.csv"
    )

    shutil.copy2(input_path, execution_input)
    shutil.copy2(prepared_xyz, execution_xyz)
    shutil.copy2(constraint_map, execution_map)
    shutil.copy2(
        PRE_QM_REPORT,
        execution_dir
        / PRE_QM_REPORT.name,
    )
    shutil.copy2(
        CONSTRUCTION_REPORT,
        execution_dir
        / CONSTRUCTION_REPORT.name,
    )

    manifest_path = (
        execution_dir
        / "QM_F06_UPPER_V5B_EXECUTION_MANIFEST.json"
    )

    manifest = {
        "decision": (
            "QM_F06_UPPER_V5B_EXECUTION_PREFLIGHT_PASS"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "overall_pass": True,
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "orca_executable": str(ORCA),
        "input_parameters": {
            "method": (
                "PBE0 D4 def2-TZVP def2/J "
                "RIJCOSX TightSCF Opt DefGrid3"
            ),
            "charge": 0,
            "multiplicity": 1,
            "electron_count": electron_count,
            "nprocs": NPROCS,
            "maxcore_mb_per_process": MAXCORE_MB,
            "scf_maxiter": SCF_MAXITER,
            "geometry_maxiter": GEOM_MAXITER,
            "coordinate_count": EXPECTED_COUNT,
            "composition": dict(
                sorted(composition.items())
            ),
            "fixed_indices": fixed_indices,
            "mobile_indices": mobile_indices,
            "fresh_scf": True,
        },
        "gates": gates,
        "hashes": {
            "input": sha256(execution_input),
            "xyz": sha256(execution_xyz),
            "constraint_map": sha256(
                execution_map
            ),
        },
        "authorization": {
            "orca_execution_authorized": True,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    launch_script = (
        execution_dir
        / "run_v5b_supervised.sh"
    )

    launch_script.write_text(
        f"""#!/bin/bash
set -u

cd "$(dirname "$0")" || exit 1

echo "$$" > v5b.supervisor_pid

"{ORCA}" v5b.inp > v5b.out 2> v5b.stderr &
orca_pid="$!"

echo "$orca_pid" > v5b.orca_pid

wait "$orca_pid"
status="$?"

echo "$status" > v5b.exit_status
exit "$status"
""",
        encoding="utf-8",
    )

    launch_script.chmod(
        launch_script.stat().st_mode
        | stat.S_IXUSR
        | stat.S_IXGRP
        | stat.S_IXOTH
    )

    latest_pointer = (
        EXECUTION_ROOT
        / "LATEST_V5B_EXECUTION.txt"
    )

    latest_pointer.write_text(
        str(execution_dir.relative_to(ROOT))
        + "\n",
        encoding="utf-8",
    )

    # Launch from its own working directory.
    supervisor_log = (
        execution_dir
        / "v5b.supervisor.log"
    )

    with supervisor_log.open(
        "wb"
    ) as log_handle:
        process = subprocess.Popen(
            [str(launch_script)],
            cwd=execution_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )

    (
        execution_dir
        / "v5b.nohup_pid"
    ).write_text(
        f"{process.pid}\n",
        encoding="utf-8",
    )

    print("=" * 88)
    print("QM_F06 UPPER V5-B PREPARATION, PREFLIGHT AND LAUNCH")
    print("=" * 88)

    for name, passed in gates.items():
        print(
            f"{name:38s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Atom count:", EXPECTED_COUNT)
    print("Composition:", dict(composition))
    print("Electron count:", electron_count)
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("Fixed indices:", fixed_indices)
    print()
    print(
        "Decision:",
        manifest["decision"],
    )
    print("Execution directory:", execution_dir)
    print("Manifest:", manifest_path)
    print("Launch script:", launch_script)
    print("Supervisor PID:", process.pid)
    print()
    print("ORCA execution authorized: True")
    print("ORCA execution started: True")


if __name__ == "__main__":
    main()
