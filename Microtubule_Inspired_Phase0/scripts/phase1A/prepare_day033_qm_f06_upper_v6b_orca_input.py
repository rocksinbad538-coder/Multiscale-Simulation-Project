#!/usr/bin/env python3
"""
Prepare the formally audited QM_F06 UPPER V6-B ORCA optimization input.

This script:
- requires the V6-B formal pre-QM gate to have passed;
- uses the audited 48-atom V6-B geometry;
- preserves the validated PBE0-D4/def2-TZVP protocol;
- remaps fixed-atom indices after removal of four hydrogens;
- keeps all V6-B-modified atoms mobile;
- writes deterministic input, XYZ, constraint-map and JSON artifacts.

ORCA execution remains blocked until an independent execution audit passes.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V6B_AUDIT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_pre_qm_audit"
)

V6B_XYZ = (
    V6B_AUDIT_DIR
    / "QM_F06_UPPER_V6B_start.xyz"
)

V6B_REPORT = (
    V6B_AUDIT_DIR
    / "QM_F06_UPPER_V6B_PRE_QM_AUDIT.json"
)

V6A_MAP = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6a_topology_closure/"
    "QM_F06_UPPER_V6A_atom_role_provenance_map.csv"
)

V5B_CONSTRAINT_MAP = ROOT / (
    "runs/phase1A/"
    "day031_qm_f06_upper_v5b_orca_input/"
    "QM_F06_UPPER_V5B_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_orca_input"
)

OUTPUT_INPUT = OUTPUT_DIR / "v6b.inp"
OUTPUT_XYZ = OUTPUT_DIR / "v6b_start.xyz"

OUTPUT_CONSTRAINT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_constraint_map.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V6B_ORCA_INPUT_REPORT.json"
)

EXPECTED_DECISION = (
    "QM_F06_UPPER_V6B_PRE_QM_GATE_PASS_"
    "ORCA_INPUT_DESIGN_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 48

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 19,
})

MODIFIED_MOBILE_ATOMS = {
    "S:1739",
    "BR4:UPPER:14:1",
    "HCAPV5B:UPPER:BR4_14_1:BR4_14_2",
}

REMOVED_HYDROGENS = {
    "H4:UPPER:0116:0",
    "H4:UPPER:0117:0",
    "H4:UPPER:0170:0",
    "H4:UPPER:0203:0",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


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
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path) -> list[dict]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(lines[0].strip())

    if len(lines) < atom_count + 2:
        raise RuntimeError(
            f"Incomplete XYZ file: {path}"
        )

    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ record at index {index}"
            )

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    return atoms


def main() -> None:
    for path in (
        V6B_XYZ,
        V6B_REPORT,
        V6A_MAP,
        V5B_CONSTRAINT_MAP,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    audit_report = json.loads(
        V6B_REPORT.read_text(encoding="utf-8")
    )

    if audit_report.get("decision") != EXPECTED_DECISION:
        raise RuntimeError(
            "V6-B formal pre-QM gate has not passed. "
            f"Found decision: {audit_report.get('decision')}"
        )

    atoms = read_xyz(V6B_XYZ)

    if len(atoms) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_ATOM_COUNT} atoms; "
            f"found {len(atoms)}."
        )

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    if composition != EXPECTED_COMPOSITION:
        raise RuntimeError(
            "Unexpected V6-B composition: "
            f"{dict(composition)}"
        )

    v6a_rows = read_csv(V6A_MAP)
    v5b_constraint_rows = read_csv(
        V5B_CONSTRAINT_MAP
    )

    fixed_by_atom = {
        row["atom_id"]: parse_bool(
            row.get("v5b_fixed", "")
        )
        for row in v5b_constraint_rows
    }

    retained_rows = [
        row
        for row in v6a_rows
        if parse_bool(
            row.get("v6a_retained", "false")
        )
    ]

    retained_rows.sort(
        key=lambda row: int(
            row["v6a_index_0based"]
        )
    )

    if len(retained_rows) != len(atoms):
        raise RuntimeError(
            "V6-A retained-map count does not match "
            "the audited V6-B XYZ."
        )

    atom_ids = [
        row["atom_id"]
        for row in retained_rows
    ]

    if any(
        atom_id in REMOVED_HYDROGENS
        for atom_id in atom_ids
    ):
        raise RuntimeError(
            "A formally removed hydrogen remains "
            "in the retained V6-B map."
        )

    constraint_rows = []
    fixed_indices = []
    mobile_indices = []

    for new_index, (
        map_row,
        atom,
    ) in enumerate(zip(retained_rows, atoms)):
        atom_id = map_row["atom_id"]

        expected_element = map_row["element"]

        if atom["element"] != expected_element:
            raise RuntimeError(
                f"Element mismatch at index {new_index}: "
                f"{atom_id}, map={expected_element}, "
                f"XYZ={atom['element']}"
            )

        inherited_fixed = fixed_by_atom.get(
            atom_id,
            False,
        )

        fixed = (
            inherited_fixed
            and atom_id not in MODIFIED_MOBILE_ATOMS
        )

        mobile = not fixed

        if fixed:
            fixed_indices.append(new_index)
        else:
            mobile_indices.append(new_index)

        output_row = dict(map_row)

        output_row.update({
            "v6b_index_0based": new_index,
            "v6b_inherited_fixed": inherited_fixed,
            "v6b_fixed": fixed,
            "v6b_mobile": mobile,
            "v6b_constraint_basis": (
                "V6B_MODIFIED_REGION_MOBILE"
                if atom_id in MODIFIED_MOBILE_ATOMS
                else (
                    "INHERITED_VALIDATED_FIXED_CORE"
                    if fixed
                    else "INHERITED_MOBILE_REGION"
                )
            ),
        })

        constraint_rows.append(output_row)

    for atom_id in MODIFIED_MOBILE_ATOMS:
        if atom_id not in atom_ids:
            raise RuntimeError(
                f"Modified atom missing from V6-B: {atom_id}"
            )

        row = next(
            item
            for item in constraint_rows
            if item["atom_id"] == atom_id
        )

        if parse_bool(row["v6b_fixed"]):
            raise RuntimeError(
                f"Modified atom was incorrectly fixed: {atom_id}"
            )

    OUTPUT_XYZ.write_text(
        V6B_XYZ.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    constraint_lines = "\n".join(
        f"    {{ C {index} C }}"
        for index in fixed_indices
    )

    coordinate_lines = []

    for atom in atoms:
        x_value, y_value, z_value = atom["xyz_A"]

        coordinate_lines.append(
            f"{atom['element']:2s} "
            f"{x_value: .12f} "
            f"{y_value: .12f} "
            f"{z_value: .12f}"
        )

    input_text = f"""! PBE0 D4 def2-TZVP def2/J RIJCOSX TightSCF Opt DefGrid3

%pal
  nprocs 4
end

%maxcore 2500

%scf
  MaxIter 500
end

%geom
  MaxIter 250
  Constraints
{constraint_lines}
  end
end

# QM_F06 UPPER V6-B topology-closed geometric reference
# Fresh SCF; audited embedded XYZ; repaired closure region mobile
* xyz 0 1
{chr(10).join(coordinate_lines)}
*
"""

    OUTPUT_INPUT.write_text(
        input_text,
        encoding="utf-8",
    )

    fieldnames = list(constraint_rows[0].keys())

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
        writer.writerows(constraint_rows)

    report = {
        "decision": (
            "QM_F06_UPPER_V6B_ORCA_INPUT_PREPARED_"
            "EXECUTION_AUDIT_REQUIRED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "model": "QM_F06_UPPER_V6B",
        "coordinate_mode": "EMBEDDED_XYZ",
        "charge": 0,
        "multiplicity": 1,
        "atom_count": len(atoms),
        "composition": dict(
            sorted(composition.items())
        ),
        "fixed_atom_count": len(fixed_indices),
        "mobile_atom_count": len(mobile_indices),
        "fixed_indices_0based": fixed_indices,
        "mobile_indices_0based": mobile_indices,
        "modified_mobile_atoms": sorted(
            MODIFIED_MOBILE_ATOMS
        ),
        "removed_hydrogens": sorted(
            REMOVED_HYDROGENS
        ),
        "protocol": {
            "method": "PBE0",
            "dispersion": "D4",
            "basis": "def2-TZVP",
            "auxiliary_basis": "def2/J",
            "approximation": "RIJCOSX",
            "scf": "TightSCF",
            "grid": "DefGrid3",
            "optimization": True,
            "nprocs": 4,
            "maxcore_MB_per_process": 2500,
            "scf_maxiter": 500,
            "geometry_maxiter": 250,
        },
        "files": {
            "input": str(
                OUTPUT_INPUT.relative_to(ROOT)
            ),
            "start_xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "constraint_map": str(
                OUTPUT_CONSTRAINT_MAP.relative_to(ROOT)
            ),
            "pre_qm_report": str(
                V6B_REPORT.relative_to(ROOT)
            ),
        },
        "sha256": {
            "input": sha256(OUTPUT_INPUT),
            "start_xyz": sha256(OUTPUT_XYZ),
            "constraint_map": sha256(
                OUTPUT_CONSTRAINT_MAP
            ),
            "pre_qm_report": sha256(
                V6B_REPORT
            ),
        },
        "execution_audit_authorized": True,
        "orca_execution_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 96)
    print("QM_F06 UPPER V6-B ORCA INPUT PREPARATION")
    print("=" * 96)
    print("Coordinate mode: EMBEDDED_XYZ")
    print("Charge: 0")
    print("Multiplicity: 1")
    print("Atom count:", len(atoms))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("Modified atoms mobile: PASS")
    print("Removed hydrogens absent: PASS")

    print()
    print("Fixed indices:", fixed_indices)
    print()
    print("Input:", OUTPUT_INPUT)
    print("XYZ:", OUTPUT_XYZ)
    print("Constraint map:", OUTPUT_CONSTRAINT_MAP)
    print("Report:", OUTPUT_REPORT)

    print()
    print("Decision:", report["decision"])
    print("Execution audit authorized: True")
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
