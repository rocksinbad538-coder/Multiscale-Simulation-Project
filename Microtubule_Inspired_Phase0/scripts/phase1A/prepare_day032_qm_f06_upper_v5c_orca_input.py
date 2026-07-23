#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SOURCE_INPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_orca_input"
)

SOURCE_INPUT = SOURCE_INPUT_DIR / "v5b.inp"

SOURCE_CONSTRAINT_MAP = SOURCE_INPUT_DIR / (
    "QM_F06_UPPER_V5B_constraint_map.csv"
)

V5C_CONSTRUCTION_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_construction"
)

V5C_SOURCE_XYZ = V5C_CONSTRUCTION_DIR / (
    "QM_F06_UPPER_V5C_start.xyz"
)

V5C_CONSTRUCTION_REPORT = V5C_CONSTRUCTION_DIR / (
    "QM_F06_UPPER_V5C_CONSTRUCTION_REPORT.json"
)

V5C_PRE_QM_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_pre_qm_audit"
)

V5C_PRE_QM_REPORT = V5C_PRE_QM_DIR / (
    "QM_F06_UPPER_V5C_PRE_QM_AUDIT.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_orca_input"
)

OUTPUT_INPUT = OUTPUT_DIR / "v5c.inp"
OUTPUT_XYZ = OUTPUT_DIR / "v5c_start.xyz"

OUTPUT_CONSTRAINT_MAP = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_constraint_map.csv"
)

OUTPUT_REPORT = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_ORCA_INPUT_REPORT.json"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1

REPAIRED_ATOMS = {
    "S:1739",
    "BR4:UPPER:14:1",
    "H4:UPPER:0203:0",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_xyz(path: Path) -> list[dict]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0])
    atoms = []

    for line in lines[2:2 + count]:
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ coordinate line: {line}"
            )

        atoms.append({
            "element": fields[0],
            "x": float(fields[1]),
            "y": float(fields[2]),
            "z": float(fields[3]),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: expected {count}, "
            f"found {len(atoms)}"
        )

    return atoms


def parse_embedded_xyz_block(text: str):
    pattern = re.compile(
        r"(?ms)"
        r"^(?P<prefix>[ \t]*\*[ \t]+xyz[ \t]+"
        r"(?P<charge>[-+]?\d+)[ \t]+"
        r"(?P<multiplicity>\d+)[ \t]*\n)"
        r"(?P<coordinates>.*?)"
        r"^(?P<suffix>[ \t]*\*[ \t]*$)"
    )

    matches = list(pattern.finditer(text))

    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one embedded ORCA "
            f"'* xyz charge multiplicity' block; "
            f"found {len(matches)}."
        )

    match = matches[0]

    coordinate_lines = [
        line
        for line in match.group(
            "coordinates"
        ).splitlines()
        if line.strip()
    ]

    return match, coordinate_lines


def format_coordinates(atoms: list[dict]) -> str:
    lines = []

    for atom in atoms:
        lines.append(
            f"{atom['element']:2s} "
            f"{atom['x']: .12f} "
            f"{atom['y']: .12f} "
            f"{atom['z']: .12f}"
        )

    return "\n".join(lines) + "\n"


def main() -> None:
    required = [
        SOURCE_INPUT,
        SOURCE_CONSTRAINT_MAP,
        V5C_SOURCE_XYZ,
        V5C_CONSTRUCTION_REPORT,
        V5C_PRE_QM_REPORT,
    ]

    for path in required:
        require_file(path)

    construction = json.loads(
        V5C_CONSTRUCTION_REPORT.read_text(
            encoding="utf-8"
        )
    )

    pre_qm = json.loads(
        V5C_PRE_QM_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not construction.get(
        "pre_qm_audit_authorized",
        False,
    ):
        raise RuntimeError(
            "V5-C construction did not authorize "
            "the formal pre-QM audit."
        )

    if not pre_qm.get(
        "authorization",
        {},
    ).get(
        "orca_input_design_authorized",
        False,
    ):
        raise RuntimeError(
            "V5-C pre-QM gate did not authorize "
            "ORCA input design."
        )

    atoms = read_xyz(V5C_SOURCE_XYZ)

    if len(atoms) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_ATOM_COUNT} atoms, "
            f"found {len(atoms)}."
        )

    source_text = SOURCE_INPUT.read_text(
        encoding="utf-8"
    )

    match, source_coordinate_lines = (
        parse_embedded_xyz_block(source_text)
    )

    charge = int(match.group("charge"))
    multiplicity = int(
        match.group("multiplicity")
    )

    if charge != EXPECTED_CHARGE:
        raise RuntimeError(
            f"Unexpected charge: {charge}"
        )

    if multiplicity != EXPECTED_MULTIPLICITY:
        raise RuntimeError(
            f"Unexpected multiplicity: {multiplicity}"
        )

    if len(source_coordinate_lines) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            "Source ORCA coordinate count mismatch: "
            f"{len(source_coordinate_lines)}"
        )

    new_coordinate_text = format_coordinates(
        atoms
    )

    new_block = (
        match.group("prefix")
        + new_coordinate_text
        + match.group("suffix")
    )

    output_text = (
        source_text[:match.start()]
        + new_block
        + source_text[match.end():]
    )

    output_text = output_text.replace(
        "QM_F06 UPPER V5-B selective boundary "
        "expansion and chemical repair",
        "QM_F06 UPPER V5-C baseline-aware "
        "local geometry repair",
    )

    output_text = output_text.replace(
        "# Fresh SCF; A:UPPER:14:4 released; "
        "all restored atoms and passivants mobile",
        "# Fresh SCF; V5-C repaired region mobile; "
        "validated V5-B constraints retained",
    )

    output_match, output_coordinate_lines = (
        parse_embedded_xyz_block(output_text)
    )

    if len(output_coordinate_lines) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            "Generated ORCA coordinate count mismatch: "
            f"{len(output_coordinate_lines)}"
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copy2(
        V5C_SOURCE_XYZ,
        OUTPUT_XYZ,
    )

    OUTPUT_INPUT.write_text(
        output_text,
        encoding="utf-8",
    )

    with SOURCE_CONSTRAINT_MAP.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            "Constraint-map atom-count mismatch: "
            f"{len(rows)}"
        )

    fieldnames = list(rows[0])

    for field in (
        "v5c_fixed",
        "v5c_mobile",
        "v5c_constraint_basis",
    ):
        if field not in fieldnames:
            fieldnames.append(field)

    for row in rows:
        row["v5c_fixed"] = row["v5b_fixed"]
        row["v5c_mobile"] = row["v5b_mobile"]

        if row["atom_id"] in REPAIRED_ATOMS:
            row["v5c_constraint_basis"] = (
                "V5C_REPAIRED_MOBILE_REGION"
            )
        else:
            row["v5c_constraint_basis"] = (
                "INHERITED_VALIDATED_V5B_CONSTRAINT"
            )

    repaired_not_mobile = [
        row["atom_id"]
        for row in rows
        if (
            row["atom_id"] in REPAIRED_ATOMS
            and row["v5c_mobile"].strip().lower()
            != "true"
        )
    ]

    if repaired_not_mobile:
        raise RuntimeError(
            "Repaired atoms are not mobile: "
            f"{repaired_not_mobile}"
        )

    with OUTPUT_CONSTRAINT_MAP.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)

    fixed_indices = sorted(
        int(row["index_0based"])
        for row in rows
        if row["v5c_fixed"].strip().lower()
        == "true"
    )

    mobile_indices = sorted(
        int(row["index_0based"])
        for row in rows
        if row["v5c_mobile"].strip().lower()
        == "true"
    )

    constraint_indices = sorted(
        int(value)
        for value in re.findall(
            r"\{\s*C\s+(\d+)\s+C\s*\}",
            output_text,
        )
    )

    if constraint_indices != fixed_indices:
        raise RuntimeError(
            "Embedded ORCA constraints do not match "
            "the constraint map.\n"
            f"Input constraints: {constraint_indices}\n"
            f"Map fixed atoms: {fixed_indices}"
        )

    report = {
        "decision": (
            "QM_F06_UPPER_V5C_ORCA_INPUT_"
            "PREPARED_EXECUTION_AUDIT_REQUIRED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "coordinate_mode": "EMBEDDED_XYZ",
        "charge": charge,
        "multiplicity": multiplicity,
        "atom_count": len(atoms),
        "fixed_atom_count": len(
            fixed_indices
        ),
        "mobile_atom_count": len(
            mobile_indices
        ),
        "constraint_indices": (
            constraint_indices
        ),
        "repaired_atoms": sorted(
            REPAIRED_ATOMS
        ),
        "input_derivation": (
            "Validated V5-B ORCA protocol retained. "
            "Only the embedded Cartesian coordinates "
            "and descriptive comments were replaced "
            "with the formally audited V5-C geometry."
        ),
        "files": {
            "input": str(
                OUTPUT_INPUT.relative_to(ROOT)
            ),
            "start_xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "constraint_map": str(
                OUTPUT_CONSTRAINT_MAP.relative_to(
                    ROOT
                )
            ),
            "construction_report": str(
                V5C_CONSTRUCTION_REPORT.relative_to(
                    ROOT
                )
            ),
            "pre_qm_report": str(
                V5C_PRE_QM_REPORT.relative_to(ROOT)
            ),
        },
        "sha256": {
            "input": sha256(OUTPUT_INPUT),
            "start_xyz": sha256(OUTPUT_XYZ),
            "constraint_map": sha256(
                OUTPUT_CONSTRAINT_MAP
            ),
            "construction_report": sha256(
                V5C_CONSTRUCTION_REPORT
            ),
            "pre_qm_report": sha256(
                V5C_PRE_QM_REPORT
            ),
        },
        "execution_audit_authorized": True,
        "orca_execution_authorized": False,
        "RESP_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 92)
    print("QM_F06 UPPER V5-C ORCA INPUT PREPARATION")
    print("=" * 92)
    print("Coordinate mode: EMBEDDED_XYZ")
    print("Charge:", charge)
    print("Multiplicity:", multiplicity)
    print("Atom count:", len(atoms))
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("Constraint-map agreement: PASS")
    print("Repaired atoms mobile: PASS")
    print()
    print("Input:", OUTPUT_INPUT)
    print("XYZ:", OUTPUT_XYZ)
    print(
        "Constraint map:",
        OUTPUT_CONSTRAINT_MAP,
    )
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "Decision:",
        report["decision"],
    )
    print(
        "Execution audit authorized: True"
    )
    print(
        "ORCA execution authorized: False"
    )


if __name__ == "__main__":
    main()
