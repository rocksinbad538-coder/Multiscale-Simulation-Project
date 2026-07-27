#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import re


ROOT = Path.cwd()

V7A_DIR = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_formal_construction"
)

PRE_QM_DIR = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_pre_qm_audit"
)

V6B_INPUT_DIR = (
    ROOT
    / "runs/phase1A/day033_qm_f06_upper_v6b_orca_input"
)

V7A_XYZ = (
    V7A_DIR
    / "QM_F06_UPPER_V7A_start.xyz"
)

V7A_MAP = (
    V7A_DIR
    / "QM_F06_UPPER_V7A_atom_role_provenance_map.csv"
)

V7A_CONSTRUCTION_REPORT = (
    V7A_DIR
    / "QM_F06_UPPER_V7A_CONSTRUCTION_REPORT.json"
)

V7A_PRE_QM_REPORT = (
    PRE_QM_DIR
    / "QM_F06_UPPER_V7A_PRE_QM_AUDIT.json"
)

V6B_INPUT = (
    V6B_INPUT_DIR
    / "v6b.inp"
)

V6B_CONSTRAINT_MAP = (
    V6B_INPUT_DIR
    / "QM_F06_UPPER_V6B_constraint_map.csv"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_orca_input"
)

OUTPUT_INPUT = OUTPUT_DIR / "v7a.inp"
OUTPUT_XYZ = OUTPUT_DIR / "v7a_start.xyz"

OUTPUT_CONSTRAINT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_constraint_map.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_ORCA_INPUT_REPORT.json"
)

EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_ATOM_COUNT = 52

NEW_OR_REPAIRED_ATOMS = {
    "A:UPPER:13:1",
    "A:UPPER:14:0",
    "A:UPPER:14:2",
    "HCAPV7:UPPER:A13_1:A11_1",
    "HCAPV7:UPPER:A14_0:A13_M1",
    "H4:UPPER:0045:0",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty source file: {path}"
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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="strict",
    ).splitlines()

    declared = int(lines[0].strip())
    atom_lines = lines[2:2 + declared]

    if len(atom_lines) != declared:
        raise RuntimeError(
            "XYZ atom-count mismatch: "
            f"declared={declared}, parsed={len(atom_lines)}"
        )

    atoms = []

    for index, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ record {index}: {line}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "x_A": float(fields[1]),
            "y_A": float(fields[2]),
            "z_A": float(fields[3]),
        })

    return atoms


def parse_embedded_xyz_block(
    text: str,
) -> tuple[int, int, int, int]:
    matches = list(
        re.finditer(
            r"(?ms)^[ \t]*\*[ \t]+xyz[ \t]+"
            r"(-?\d+)[ \t]+(\d+)[ \t]*\n"
            r"(.*?)"
            r"^[ \t]*\*[ \t]*$",
            text,
        )
    )

    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one embedded ORCA XYZ block; "
            f"found {len(matches)}"
        )

    match = matches[0]

    return (
        int(match.group(1)),
        int(match.group(2)),
        match.start(),
        match.end(),
    )


def rebuild_geom_constraints(
    text: str,
    fixed_indices: list[int],
) -> str:
    """
    Reconstruct the unique ORCA %geom block while preserving
    all non-constraint directives and replacing every prior
    Constraints subsection with exactly one canonical block.
    """

    pattern = re.compile(
        r"(?ms)"
        r"^(?P<indent>[ \t]*)%geom[ \t]*\n"
        r"(?P<body>.*?)"
        r"^(?P=indent)end[ \t]*$"
    )

    matches = list(pattern.finditer(text))

    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one %geom block; "
            f"found {len(matches)}"
        )

    match = matches[0]
    indent = match.group("indent")
    body = match.group("body")

    # Remove every pre-existing Constraints subsection.
    body_without_constraints = re.sub(
        r"(?ms)"
        r"^[ \t]*Constraints[ \t]*\n"
        r".*?"
        r"^[ \t]*end[ \t]*\n?",
        "",
        body,
    )

    preserved_lines = [
        line.rstrip()
        for line in body_without_constraints.splitlines()
        if line.strip()
    ]

    geom_lines = [
        f"{indent}%geom",
        *[
            f"{indent}  {line.strip()}"
            for line in preserved_lines
        ],
        f"{indent}  Constraints",
        *[
            f"{indent}    {{ C {index} C }}"
            for index in fixed_indices
        ],
        f"{indent}  end",
        f"{indent}end",
    ]

    reconstructed = "\n".join(geom_lines)

    return (
        text[:match.start()]
        + reconstructed
        + text[match.end():]
    )



def main() -> None:
    for path in (
        V7A_XYZ,
        V7A_MAP,
        V7A_CONSTRUCTION_REPORT,
        V7A_PRE_QM_REPORT,
        V6B_INPUT,
        V6B_CONSTRAINT_MAP,
    ):
        require_file(path)

    construction = json.loads(
        V7A_CONSTRUCTION_REPORT.read_text(
            encoding="utf-8",
        )
    )

    pre_qm = json.loads(
        V7A_PRE_QM_REPORT.read_text(
            encoding="utf-8",
        )
    )

    if construction.get("decision") != (
        "QM_F06_UPPER_V7A_FORMALLY_CONSTRUCTED_"
        "GLOBAL_PRE_QM_AUDIT_REQUIRED"
    ):
        raise RuntimeError(
            "V7-A construction is not authorized"
        )

    if pre_qm.get("decision") != (
        "QM_F06_UPPER_V7A_PRE_QM_GATE_PASS_"
        "ORCA_INPUT_DESIGN_AUTHORIZED"
    ):
        raise RuntimeError(
            "V7-A pre-QM gate did not authorize input design"
        )

    atoms = read_xyz(V7A_XYZ)
    map_rows = read_csv(V7A_MAP)
    v6b_map_rows = read_csv(
        V6B_CONSTRAINT_MAP
    )

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    if len(atoms) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Unexpected atom count: {len(atoms)}"
        )

    if len(map_rows) != len(atoms):
        raise RuntimeError(
            "V7-A map/XYZ row-count mismatch"
        )

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    if dict(composition) != EXPECTED_COMPOSITION:
        raise RuntimeError(
            f"Unexpected composition: {dict(composition)}"
        )

    for index, (atom, row) in enumerate(
        zip(atoms, map_rows)
    ):
        if int(row["v7a_index_0based"]) != index:
            raise RuntimeError(
                "Nonsequential V7-A map index at "
                f"{index}"
            )

        if atom["element"] != row["element"]:
            raise RuntimeError(
                "V7-A map/XYZ element mismatch at "
                f"{index}: "
                f"{row['atom_id']}"
            )

    v6b_fixed_by_id = {}

    for row in v6b_map_rows:
        atom_id = row["atom_id"]

        raw_fixed = row.get(
            "v6b_fixed",
            "",
        ).strip().lower()

        v6b_fixed_by_id[atom_id] = (
            raw_fixed == "true"
        )

    constraint_records = []
    fixed_indices = []
    mobile_indices = []

    for index, row in enumerate(map_rows):
        atom_id = row["atom_id"]

        inherited_fixed = (
            v6b_fixed_by_id.get(
                atom_id,
                False,
            )
        )

        fixed = (
            inherited_fixed
            and atom_id
            not in NEW_OR_REPAIRED_ATOMS
        )

        mobile = not fixed

        if fixed:
            fixed_indices.append(index)
        else:
            mobile_indices.append(index)

        record = dict(row)
        record.update({
            "v7a_index_0based": index,
            "v7a_inherited_v6b_fixed": (
                inherited_fixed
            ),
            "v7a_fixed": fixed,
            "v7a_mobile": mobile,
            "v7a_constraint_basis": (
                "INHERITED_V6B_FIXED_ATOM"
                if fixed
                else (
                    "V7A_NEW_OR_REPAIRED_MOBILE"
                    if atom_id
                    in NEW_OR_REPAIRED_ATOMS
                    else
                    "INHERITED_MOBILE_REGION"
                )
            ),
        })

        constraint_records.append(record)

    if fixed_indices != sorted(fixed_indices):
        raise RuntimeError(
            "Fixed indices are not sorted"
        )

    if (
        len(fixed_indices)
        + len(mobile_indices)
        != len(atoms)
    ):
        raise RuntimeError(
            "Incomplete fixed/mobile partition"
        )

    modified_mobile = all(
        constraint_records[
            int(row["v7a_index_0based"])
        ]["v7a_mobile"]
        for row in map_rows
        if row["atom_id"]
        in NEW_OR_REPAIRED_ATOMS
    )

    if not modified_mobile:
        raise RuntimeError(
            "At least one V7-A modified atom is fixed"
        )

    template = V6B_INPUT.read_text(
        encoding="utf-8",
        errors="strict",
    )

    charge, multiplicity, start, end = (
        parse_embedded_xyz_block(template)
    )

    if charge != 0 or multiplicity != 1:
        raise RuntimeError(
            "Unexpected V6-B charge/multiplicity"
        )

    protocol = rebuild_geom_constraints(
        template,
        fixed_indices,
    )

    _, _, block_start, block_end = (
        parse_embedded_xyz_block(protocol)
    )

    coordinate_lines = [
        "* xyz 0 1",
        *[
            (
                f"{atom['element']:2s} "
                f"{atom['x_A']: .12f} "
                f"{atom['y_A']: .12f} "
                f"{atom['z_A']: .12f}"
            )
            for atom in atoms
        ],
        "*",
    ]

    new_xyz_block = (
        "\n".join(coordinate_lines)
        + "\n"
    )

    protocol = (
        protocol[:block_start]
        + new_xyz_block
        + protocol[block_end:]
    )

    protocol = re.sub(
        r"(?m)^#.*V6-B.*$",
        "# QM_F06 UPPER V7-A canonical boundary expansion",
        protocol,
    )

    protocol = re.sub(
        r"(?m)^# Fresh SCF.*$",
        (
            "# Fresh SCF; V7-A repaired and newly added "
            "atoms mobile; inherited fixed region preserved"
        ),
        protocol,
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_INPUT.write_text(
        protocol,
        encoding="utf-8",
    )

    OUTPUT_XYZ.write_text(
        V7A_XYZ.read_text(
            encoding="utf-8",
        ),
        encoding="utf-8",
    )

    with OUTPUT_CONSTRAINT_MAP.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                constraint_records[0].keys()
            ),
        )
        writer.writeheader()
        writer.writerows(
            constraint_records
        )

    report = {
        "model": "QM_F06_UPPER_V7A",
        "decision": (
            "QM_F06_UPPER_V7A_ORCA_INPUT_PREPARED_"
            "EXECUTION_AUDIT_REQUIRED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "coordinate_mode": "EMBEDDED_XYZ",
        "protocol_source": str(
            V6B_INPUT.relative_to(ROOT)
        ),
        "protocol": {
            "method": "PBE0",
            "dispersion": "D4",
            "basis": "def2-TZVP",
            "auxiliary_basis": "def2/J",
            "approximation": "RIJCOSX",
            "SCF": "TightSCF",
            "grid": "DefGrid3",
            "optimization": True,
            "charge": 0,
            "multiplicity": 1,
            "nprocs": 4,
            "maxcore_MB_per_process": 2500,
        },
        "structure": {
            "atom_count": len(atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "fixed_atom_count": len(
                fixed_indices
            ),
            "mobile_atom_count": len(
                mobile_indices
            ),
            "fixed_indices_0based": (
                fixed_indices
            ),
            "mobile_indices_0based": (
                mobile_indices
            ),
            "new_or_repaired_atoms_mobile": (
                modified_mobile
            ),
        },
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
                V7A_CONSTRUCTION_REPORT.relative_to(
                    ROOT
                )
            ),
            "pre_qm_report": str(
                V7A_PRE_QM_REPORT.relative_to(ROOT)
            ),
        },
        "sha256": {
            "input": sha256(OUTPUT_INPUT),
            "start_xyz": sha256(OUTPUT_XYZ),
            "constraint_map": sha256(
                OUTPUT_CONSTRAINT_MAP
            ),
            "construction_report": sha256(
                V7A_CONSTRUCTION_REPORT
            ),
            "pre_qm_report": sha256(
                V7A_PRE_QM_REPORT
            ),
        },
        "execution_audit_authorized": True,
        "ORCA_execution_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 100)
    print("QM_F06 UPPER V7-A ORCA INPUT PREPARATION")
    print("=" * 100)
    print("Coordinate mode: EMBEDDED_XYZ")
    print("Charge:", charge)
    print("Multiplicity:", multiplicity)
    print("Atom count:", len(atoms))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("Fixed indices:", fixed_indices)
    print("New/repaired atoms mobile: PASS")

    print()
    print("Input:", OUTPUT_INPUT)
    print("XYZ:", OUTPUT_XYZ)
    print(
        "Constraint map:",
        OUTPUT_CONSTRAINT_MAP,
    )
    print("Report:", OUTPUT_REPORT)

    print()
    print("Decision:", report["decision"])
    print("Execution audit authorized: True")
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
