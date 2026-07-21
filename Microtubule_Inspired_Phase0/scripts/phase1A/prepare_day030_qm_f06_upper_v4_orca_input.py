#!/usr/bin/env python3
"""
Prepare the QM_F06 UPPER V4 ORCA optimization input.

The electronic method and runtime controls are inherited exactly from
the validated V3-A2 workflow. Only the geometry, comments, and fixed
indices are replaced with the audited V4 values.

This script prepares the input but does not authorize or execute ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V3A2_INPUT = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow/v3a2.inp"
)

V4_XYZ = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction/"
    "QM_F06_UPPER_V4_start.xyz"
)

CONSTRAINT_MAP = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_constraint_design/"
    "QM_F06_UPPER_V4_constraint_map.csv"
)

CONSTRAINT_REPORT = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_constraint_design/"
    "QM_F06_UPPER_V4_CONSTRAINT_DESIGN.json"
)

PRE_QM_REPORT = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_pre_qm_audit/"
    "QM_F06_UPPER_V4_PRE_QM_AUDIT.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_orca_input"
)

OUTPUT_INPUT = OUTPUT_DIR / "v4.inp"

OUTPUT_XYZ = OUTPUT_DIR / "v4_start.xyz"

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_constraint_map.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_ORCA_INPUT_PREPARATION.json"
)

EXPECTED_KEYWORD_LINE = (
    "! PBE0 D4 def2-TZVP def2/J RIJCOSX "
    "TightSCF Opt DefGrid3"
)

EXPECTED_ATOM_COUNT = 46
EXPECTED_COMPOSITION = Counter({
    "B": 15,
    "N": 11,
    "H": 20,
})

EXPECTED_FIXED_INDICES = list(range(14)) + [16]

EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1
EXPECTED_NPROCS = 4
EXPECTED_MAXCORE_MB = 2500
EXPECTED_SCF_MAXITER = 500
EXPECTED_GEOM_MAXITER = 250


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


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


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
    }


def read_xyz(path: Path):
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    declared = int(lines[0].strip())

    coordinate_lines = lines[2:2 + declared]

    if len(coordinate_lines) != declared:
        raise RuntimeError(
            f"Incomplete XYZ: expected {declared}, "
            f"found {len(coordinate_lines)}"
        )

    atoms = []

    for index, line in enumerate(coordinate_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index}: {line}"
            )

        atoms.append({
            "index": index,
            "element": fields[0],
            "x": float(fields[1]),
            "y": float(fields[2]),
            "z": float(fields[3]),
        })

    return atoms


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for path in (
        V3A2_INPUT,
        V4_XYZ,
        CONSTRAINT_MAP,
        CONSTRAINT_REPORT,
        PRE_QM_REPORT,
    ):
        require_file(path)

    constraint_report = json.loads(
        CONSTRAINT_REPORT.read_text(
            encoding="utf-8"
        )
    )

    expected_constraint_decision = (
        "QM_F06_UPPER_V4_CONSTRAINT_DESIGN_PASS_"
        "ORCA_INPUT_PREPARATION_AUTHORIZED"
    )

    if (
        constraint_report["decision"]
        != expected_constraint_decision
    ):
        raise RuntimeError(
            "Unexpected constraint decision: "
            f"{constraint_report['decision']}"
        )

    if not constraint_report["authorization"][
        "orca_input_preparation_authorized"
    ]:
        raise RuntimeError(
            "ORCA input preparation is not authorized."
        )

    pre_qm_report = json.loads(
        PRE_QM_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not all(
        pre_qm_report["gates"].values()
    ):
        raise RuntimeError(
            "The V4 structural gate is not fully passed."
        )

    template_text = V3A2_INPUT.read_text(
        encoding="utf-8",
        errors="replace",
    )

    keyword_lines = [
        line.strip()
        for line in template_text.splitlines()
        if line.lstrip().startswith("!")
    ]

    if keyword_lines != [EXPECTED_KEYWORD_LINE]:
        raise RuntimeError(
            "Unexpected V3-A2 keyword line: "
            f"{keyword_lines}"
        )

    nprocs_match = re.search(
        r"\bnprocs\s+(\d+)",
        template_text,
        re.IGNORECASE,
    )

    maxcore_match = re.search(
        r"%maxcore\s+(\d+)",
        template_text,
        re.IGNORECASE,
    )

    scf_maxiter_match = re.search(
        r"%scf\b.*?\bMaxIter\s+(\d+).*?\bend\b",
        template_text,
        re.IGNORECASE | re.DOTALL,
    )

    geom_maxiter_match = re.search(
        r"%geom\b.*?\bMaxIter\s+(\d+)",
        template_text,
        re.IGNORECASE | re.DOTALL,
    )

    xyz_header_match = re.search(
        r"^\s*\*\s*xyz\s+([-+]?\d+)\s+(\d+)\s*$",
        template_text,
        re.IGNORECASE | re.MULTILINE,
    )

    if not all((
        nprocs_match,
        maxcore_match,
        scf_maxiter_match,
        geom_maxiter_match,
        xyz_header_match,
    )):
        raise RuntimeError(
            "Could not parse the complete V3-A2 template."
        )

    parsed_template = {
        "nprocs": int(nprocs_match.group(1)),
        "maxcore_mb": int(maxcore_match.group(1)),
        "scf_maxiter": int(
            scf_maxiter_match.group(1)
        ),
        "geom_maxiter": int(
            geom_maxiter_match.group(1)
        ),
        "charge": int(xyz_header_match.group(1)),
        "multiplicity": int(
            xyz_header_match.group(2)
        ),
    }

    expected_template = {
        "nprocs": EXPECTED_NPROCS,
        "maxcore_mb": EXPECTED_MAXCORE_MB,
        "scf_maxiter": EXPECTED_SCF_MAXITER,
        "geom_maxiter": EXPECTED_GEOM_MAXITER,
        "charge": EXPECTED_CHARGE,
        "multiplicity": EXPECTED_MULTIPLICITY,
    }

    if parsed_template != expected_template:
        raise RuntimeError(
            "V3-A2 template parameters differ from "
            f"expected values: {parsed_template}"
        )

    reuse_hits_template = {
        "moread": bool(re.search(
            r"\bmoread\b",
            template_text,
            re.IGNORECASE,
        )),
        "moinp": bool(re.search(
            r"%moinp\b",
            template_text,
            re.IGNORECASE,
        )),
        "gbw_reference": bool(re.search(
            r"\S+\.gbw",
            template_text,
            re.IGNORECASE,
        )),
    }

    if any(reuse_hits_template.values()):
        raise RuntimeError(
            "V3-A2 template unexpectedly reuses SCF data: "
            f"{reuse_hits_template}"
        )

    atoms = read_xyz(V4_XYZ)

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    if len(atoms) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_ATOM_COUNT} atoms; "
            f"found {len(atoms)}"
        )

    if composition != EXPECTED_COMPOSITION:
        raise RuntimeError(
            "Unexpected V4 composition: "
            f"{dict(composition)}"
        )

    constraint_rows = read_csv(CONSTRAINT_MAP)

    if len(constraint_rows) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            "Constraint map atom count mismatch."
        )

    fixed_indices = sorted(
        int(row["index_0based"])
        for row in constraint_rows
        if parse_bool(row["v4_fixed"])
    )

    mobile_indices = sorted(
        int(row["index_0based"])
        for row in constraint_rows
        if parse_bool(row["v4_mobile"])
    )

    if fixed_indices != EXPECTED_FIXED_INDICES:
        raise RuntimeError(
            "Unexpected fixed indices: "
            f"{fixed_indices}"
        )

    if sorted(
        fixed_indices + mobile_indices
    ) != list(range(EXPECTED_ATOM_COUNT)):
        raise RuntimeError(
            "Fixed/mobile indices do not partition "
            "the complete atom set."
        )

    if set(fixed_indices) & set(mobile_indices):
        raise RuntimeError(
            "Fixed and mobile indices overlap."
        )

    constraint_lines = "\n".join(
        f"    {{ C {index} C }}"
        for index in fixed_indices
    )

    coordinate_lines = "\n".join(
        (
            f"{atom['element']:2s} "
            f"{atom['x']: .12f} "
            f"{atom['y']: .12f} "
            f"{atom['z']: .12f}"
        )
        for atom in atoms
    )

    output_text = f"""\
{EXPECTED_KEYWORD_LINE}

%pal
  nprocs {EXPECTED_NPROCS}
end

%maxcore {EXPECTED_MAXCORE_MB}
%scf
  MaxIter {EXPECTED_SCF_MAXITER}
end

%geom
  MaxIter {EXPECTED_GEOM_MAXITER}
  Constraints
{constraint_lines}
  end
end

# QM_F06 UPPER Boundary V4 selected P:1523 expansion
# Fixed: retained validated real V3-A2 core only
# Mobile: all artificial caps, restored atoms, and prior mobile boundary atoms
# Fresh SCF; no MOREAD, %moinp, or GBW reuse
* xyz {EXPECTED_CHARGE} {EXPECTED_MULTIPLICITY}
{coordinate_lines}
*
"""

    OUTPUT_INPUT.write_text(
        output_text,
        encoding="utf-8",
    )

    # Copy the exact audited coordinate and constraint artifacts.
    OUTPUT_XYZ.write_bytes(
        V4_XYZ.read_bytes()
    )

    OUTPUT_MAP.write_bytes(
        CONSTRAINT_MAP.read_bytes()
    )

    prepared_text = OUTPUT_INPUT.read_text(
        encoding="utf-8"
    )

    prepared_reuse_hits = {
        "moread": bool(re.search(
            r"\bmoread\b",
            prepared_text,
            re.IGNORECASE,
        )),
        "moinp": bool(re.search(
            r"%moinp\b",
            prepared_text,
            re.IGNORECASE,
        )),
        "gbw_reference": bool(re.search(
            r"\S+\.gbw",
            prepared_text,
            re.IGNORECASE,
        )),
    }

    # Ignore the explanatory comment when testing MOREAD.
    non_comment_text = "\n".join(
        line
        for line in prepared_text.splitlines()
        if not line.lstrip().startswith("#")
    )

    effective_reuse_hits = {
        "moread": bool(re.search(
            r"\bmoread\b",
            non_comment_text,
            re.IGNORECASE,
        )),
        "moinp": bool(re.search(
            r"%moinp\b",
            non_comment_text,
            re.IGNORECASE,
        )),
        "gbw_reference": bool(re.search(
            r"\S+\.gbw",
            non_comment_text,
            re.IGNORECASE,
        )),
    }

    if any(effective_reuse_hits.values()):
        raise RuntimeError(
            "Prepared V4 input contains forbidden SCF reuse: "
            f"{effective_reuse_hits}"
        )

    electron_count = (
        composition["B"] * 5
        + composition["N"] * 7
        + composition["H"] * 1
        - EXPECTED_CHARGE
    )

    electron_parity_gate = (
        electron_count % 2 == 0
        and EXPECTED_MULTIPLICITY == 1
    )

    gates = {
        "constraint_report_authorization": True,
        "pre_qm_structural_gate": True,
        "template_keyword_line": True,
        "template_parameters": True,
        "template_fresh_scf": True,
        "atom_count": (
            len(atoms) == EXPECTED_ATOM_COUNT
        ),
        "composition": (
            composition == EXPECTED_COMPOSITION
        ),
        "charge": (
            EXPECTED_CHARGE == 0
        ),
        "multiplicity": (
            EXPECTED_MULTIPLICITY == 1
        ),
        "electron_parity": electron_parity_gate,
        "fixed_indices": (
            fixed_indices
            == EXPECTED_FIXED_INDICES
        ),
        "fixed_mobile_partition": (
            len(fixed_indices)
            + len(mobile_indices)
            == EXPECTED_ATOM_COUNT
        ),
        "prepared_fresh_scf": (
            not any(
                effective_reuse_hits.values()
            )
        ),
    }

    overall_pass = all(gates.values())

    decision = (
        "QM_F06_UPPER_V4_ORCA_INPUT_PREPARED_"
        "EXECUTION_PREFLIGHT_REQUIRED"
        if overall_pass
        else
        "QM_F06_UPPER_V4_ORCA_INPUT_PREPARATION_FAIL"
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "template_input": str(
            V3A2_INPUT.relative_to(ROOT)
        ),
        "keyword_line": EXPECTED_KEYWORD_LINE,
        "method_inheritance": (
            "EXACT_V3A2_ELECTRONIC_AND_RUNTIME_TEMPLATE"
        ),
        "parameters": {
            "charge": EXPECTED_CHARGE,
            "multiplicity": EXPECTED_MULTIPLICITY,
            "electron_count": electron_count,
            "nprocs": EXPECTED_NPROCS,
            "maxcore_mb_per_process": (
                EXPECTED_MAXCORE_MB
            ),
            "nominal_total_maxcore_mb": (
                EXPECTED_NPROCS
                * EXPECTED_MAXCORE_MB
            ),
            "scf_maxiter": EXPECTED_SCF_MAXITER,
            "geometry_maxiter": (
                EXPECTED_GEOM_MAXITER
            ),
            "atom_count": len(atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "fixed_indices": fixed_indices,
            "mobile_indices": mobile_indices,
        },
        "template_reuse_hits": (
            reuse_hits_template
        ),
        "prepared_input_reuse_hits_raw": (
            prepared_reuse_hits
        ),
        "prepared_input_reuse_hits_effective": (
            effective_reuse_hits
        ),
        "gates": gates,
        "overall_pass": overall_pass,
        "files": {
            "orca_input": str(
                OUTPUT_INPUT.relative_to(ROOT)
            ),
            "start_xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "constraint_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "template_input": sha256(V3A2_INPUT),
            "source_v4_xyz": sha256(V4_XYZ),
            "source_constraint_map": sha256(
                CONSTRAINT_MAP
            ),
            "constraint_report": sha256(
                CONSTRAINT_REPORT
            ),
            "pre_qm_report": sha256(
                PRE_QM_REPORT
            ),
            "prepared_input": sha256(
                OUTPUT_INPUT
            ),
            "copied_start_xyz": sha256(
                OUTPUT_XYZ
            ),
            "copied_constraint_map": sha256(
                OUTPUT_MAP
            ),
        },
        "authorization": {
            "orca_input_prepared": overall_pass,
            "execution_preflight_authorized": (
                overall_pass
            ),
            "orca_execution_authorized": False,
            "geometric_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 ORCA INPUT PREPARATION")
    print("=" * 78)

    for gate, passed in gates.items():
        print(
            f"{gate:38s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Atom count:", len(atoms))
    print("Composition:", dict(composition))
    print("Electron count:", electron_count)
    print(
        "Charge / multiplicity:",
        EXPECTED_CHARGE,
        EXPECTED_MULTIPLICITY,
    )
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("Fixed indices:", fixed_indices)
    print(
        "Fresh-SCF hits:",
        effective_reuse_hits,
    )
    print()
    print("Decision:", decision)
    print("Input:", OUTPUT_INPUT)
    print("XYZ:", OUTPUT_XYZ)
    print("Constraint map:", OUTPUT_MAP)
    print("Report:", OUTPUT_REPORT)
    print(
        "Execution preflight authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
