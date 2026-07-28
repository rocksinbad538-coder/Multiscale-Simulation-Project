#!/usr/bin/env python3

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import re
import shutil


ROOT = Path(__file__).resolve().parents[2]

SOURCE_EXEC_PARENT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_executions"
)

LATEST_SOURCE_FILE = (
    SOURCE_EXEC_PARENT
    / "LATEST_V7A_EXECUTION.txt"
)

RESTART_GATE_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_step1_restart_gate"
)

RESTART_GATE_REPORT = (
    RESTART_GATE_DIR
    / "QM_F06_UPPER_V7A_STEP1_RESTART_GATE.json"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_orca_input"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_INPUT = OUTPUT_DIR / "v7a_r1.inp"
OUTPUT_XYZ = OUTPUT_DIR / "v7a_r1_start.xyz"
OUTPUT_GBW = OUTPUT_DIR / "v7a_r1_source.gbw"

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_constraint_map.csv"
)

OUTPUT_EDGES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_nominal_edges.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_ORCA_INPUT_REPORT.json"
)


EXPECTED_RESTART_GATE_DECISION = (
    "QM_F06_UPPER_V7A_STEP1_"
    "STRUCTURAL_GATE_PASS_"
    "SCF_RESTART_DESIGN_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1

EXPECTED_FIXED_INDICES = [
    0, 1, 3, 4, 5, 6,
    7, 8, 9, 10, 11, 14,
]


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
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(lines[0].strip())

    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ line {index}: {line}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "x_A": float(fields[1]),
            "y_A": float(fields[2]),
            "z_A": float(fields[3]),
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"Incomplete XYZ file: {path}"
        )

    return atoms


def format_xyz_block(
    atoms: list[dict],
    charge: int,
    multiplicity: int,
) -> str:
    lines = [
        f"* xyz {charge} {multiplicity}"
    ]

    for atom in atoms:
        lines.append(
            f"{atom['element']:2s} "
            f"{atom['x_A']: .12f} "
            f"{atom['y_A']: .12f} "
            f"{atom['z_A']: .12f}"
        )

    lines.append("*")

    return "\n".join(lines)


def replace_exactly_once(
    text: str,
    pattern: str,
    replacement: str,
    label: str,
    flags: int = 0,
) -> str:
    updated, count = re.subn(
        pattern,
        replacement,
        text,
        count=1,
        flags=flags,
    )

    if count != 1:
        raise RuntimeError(
            f"{label}: expected one replacement, "
            f"observed {count}"
        )

    return updated


def main() -> None:
    source_relative = (
        LATEST_SOURCE_FILE
        .read_text(encoding="utf-8")
        .strip()
    )

    source_dir = ROOT / source_relative

    source_input = source_dir / "v7a.inp"
    source_geometry = source_dir / "v7a.xyz"
    source_gbw = source_dir / "v7a.gbw"

    source_map = (
        source_dir
        / "QM_F06_UPPER_V7A_constraint_map.csv"
    )

    source_edges = (
        source_dir
        / "QM_F06_UPPER_V7A_nominal_edges.csv"
    )

    required = (
        source_input,
        source_geometry,
        source_gbw,
        source_map,
        source_edges,
        RESTART_GATE_REPORT,
    )

    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing required restart source: {path}"
            )

    restart_gate = json.loads(
        RESTART_GATE_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if (
        restart_gate.get("decision")
        != EXPECTED_RESTART_GATE_DECISION
    ):
        raise RuntimeError(
            "Restart structural gate is not authorized: "
            f"{restart_gate.get('decision')}"
        )

    atoms = read_xyz(source_geometry)

    if len(atoms) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_ATOM_COUNT} atoms, "
            f"found {len(atoms)}"
        )

    with source_map.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        map_rows = list(csv.DictReader(handle))

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    if len(map_rows) != EXPECTED_ATOM_COUNT:
        raise RuntimeError(
            "Constraint-map atom count mismatch"
        )

    fixed_indices = [
        int(row["v7a_index_0based"])
        for row in map_rows
        if (
            row["v7a_fixed"]
            .strip()
            .lower()
            == "true"
        )
    ]

    mobile_indices = [
        int(row["v7a_index_0based"])
        for row in map_rows
        if (
            row["v7a_mobile"]
            .strip()
            .lower()
            == "true"
        )
    ]

    if fixed_indices != EXPECTED_FIXED_INDICES:
        raise RuntimeError(
            "Unexpected fixed-atom indices: "
            f"{fixed_indices}"
        )

    if sorted(fixed_indices + mobile_indices) != list(
        range(EXPECTED_ATOM_COUNT)
    ):
        raise RuntimeError(
            "Fixed/mobile partition is incomplete"
        )

    for index, row in enumerate(map_rows):
        if atoms[index]["element"] != row["element"]:
            raise RuntimeError(
                "Map/XYZ element mismatch at "
                f"index {index}"
            )

    source_text = source_input.read_text(
        encoding="utf-8",
        errors="strict",
    )

    simple_lines = [
        line
        for line in source_text.splitlines()
        if line.lstrip().startswith("!")
    ]

    if len(simple_lines) != 1:
        raise RuntimeError(
            "Expected exactly one simple-input line"
        )

    original_simple = simple_lines[0]

    required_original_tokens = (
        "PBE0",
        "D4",
        "def2-TZVP",
        "def2/J",
        "RIJCOSX",
        "TightSCF",
        "Opt",
        "DefGrid3",
    )

    if not all(
        token in original_simple
        for token in required_original_tokens
    ):
        raise RuntimeError(
            "Original protocol tokens are incomplete"
        )

    restart_simple = (
        "! PBE0 D4 def2-TZVP def2/J RIJCOSX "
        "TightSCF Opt DefGrid3 TRAH MORead"
    )

    restart_text = replace_exactly_once(
        source_text,
        r"(?m)^[ \t]*!.*$",
        restart_simple,
        "simple_keyword_line",
    )

    restart_scf_block = """%scf
  MaxIter 500
  AutoTRAH false
end

%moinp "v7a_r1_source.gbw\""""

    restart_text = replace_exactly_once(
        restart_text,
        (
            r"(?ims)"
            r"^[ \t]*%scf\b.*?"
            r"^[ \t]*end[ \t]*$"
        ),
        restart_scf_block,
        "scf_block",
    )

    new_xyz_block = format_xyz_block(
        atoms,
        EXPECTED_CHARGE,
        EXPECTED_MULTIPLICITY,
    )

    restart_text = replace_exactly_once(
        restart_text,
        (
            r"(?ims)"
            r"^[ \t]*\*[ \t]*xyz\b.*?"
            r"^[ \t]*\*[ \t]*$"
        ),
        new_xyz_block,
        "embedded_xyz_block",
    )

    if "xyzfile" in restart_text.lower():
        raise RuntimeError(
            "Unexpected xyzfile directive"
        )

    if restart_text.count("%moinp") != 1:
        raise RuntimeError(
            "Expected exactly one %moinp directive"
        )

    if not re.search(
        r"(?i)\bTRAH\b",
        restart_text,
    ):
        raise RuntimeError(
            "TRAH keyword not found"
        )

    if not re.search(
        r"(?i)\bMORead\b",
        restart_text,
    ):
        raise RuntimeError(
            "MORead keyword not found"
        )

    OUTPUT_INPUT.write_text(
        restart_text.rstrip() + "\n",
        encoding="utf-8",
    )

    xyz_lines = [
        str(len(atoms)),
        (
            "QM_F06 UPPER V7-A R1 restart geometry; "
            "validated post-step-1 v7a.xyz"
        ),
    ]

    for atom in atoms:
        xyz_lines.append(
            f"{atom['element']:2s} "
            f"{atom['x_A']: .12f} "
            f"{atom['y_A']: .12f} "
            f"{atom['z_A']: .12f}"
        )

    OUTPUT_XYZ.write_text(
        "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )

    shutil.copy2(
        source_gbw,
        OUTPUT_GBW,
    )

    shutil.copy2(
        source_map,
        OUTPUT_MAP,
    )

    shutil.copy2(
        source_edges,
        OUTPUT_EDGES,
    )

    report = {
        "model": "QM_F06_UPPER_V7A_R1",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "decision": (
            "QM_F06_UPPER_V7A_R1_ORCA_INPUT_"
            "PREPARED_EXECUTION_AUDIT_REQUIRED"
        ),
        "source_execution_directory": str(
            source_dir
        ),
        "source_execution_immutable": True,
        "restart_geometry_source": str(
            source_geometry
        ),
        "orbital_source": str(source_gbw),
        "restart_gate_report": str(
            RESTART_GATE_REPORT
        ),
        "protocol": {
            "method": "PBE0",
            "dispersion": "D4",
            "basis": "def2-TZVP",
            "auxiliary_basis": "def2/J",
            "approximation": "RIJCOSX",
            "grid": "DefGrid3",
            "scf_tolerance": "TightSCF",
            "scf_algorithm": "TRAH",
            "initial_guess": "MORead",
            "moinp": "v7a_r1_source.gbw",
            "geometry_optimization": True,
            "maximum_scf_iterations": 500,
            "maximum_geometry_iterations": 250,
            "charge": EXPECTED_CHARGE,
            "multiplicity": (
                EXPECTED_MULTIPLICITY
            ),
            "nprocs": 4,
            "maxcore_MB_per_process": 2500,
        },
        "system": {
            "atom_count": len(atoms),
            "composition": dict(
                sorted(
                    {
                        element: sum(
                            atom["element"] == element
                            for atom in atoms
                        )
                        for element in {
                            atom["element"]
                            for atom in atoms
                        }
                    }.items()
                )
            ),
            "fixed_indices_0based": (
                fixed_indices
            ),
            "fixed_atom_count": (
                len(fixed_indices)
            ),
            "mobile_atom_count": (
                len(mobile_indices)
            ),
        },
        "files": {
            "input": str(OUTPUT_INPUT),
            "xyz": str(OUTPUT_XYZ),
            "orbital_restart": str(
                OUTPUT_GBW
            ),
            "constraint_map": str(
                OUTPUT_MAP
            ),
            "nominal_edges": str(
                OUTPUT_EDGES
            ),
        },
        "sha256": {
            "source_input": sha256(
                source_input
            ),
            "source_restart_geometry": (
                sha256(source_geometry)
            ),
            "source_gbw": sha256(
                source_gbw
            ),
            "prepared_input": sha256(
                OUTPUT_INPUT
            ),
            "prepared_xyz": sha256(
                OUTPUT_XYZ
            ),
            "copied_gbw": sha256(
                OUTPUT_GBW
            ),
            "constraint_map": sha256(
                OUTPUT_MAP
            ),
            "nominal_edges": sha256(
                OUTPUT_EDGES
            ),
        },
        "authorizations": {
            "execution_audit_authorized": True,
            "ORCA_restart_execution_authorized": (
                False
            ),
            "RESP_authorized": False,
            "force_field_adoption_authorized": (
                False
            ),
            "MD_authorized": False,
        },
    }

    if (
        report["sha256"]["source_gbw"]
        != report["sha256"]["copied_gbw"]
    ):
        raise RuntimeError(
            "Copied GBW hash mismatch"
        )

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 104)
    print(
        "QM_F06 UPPER V7-A R1 "
        "ORCA INPUT PREPARATION"
    )
    print("=" * 104)

    print(
        "Restart gate:",
        restart_gate["decision"],
    )
    print("Atoms:", len(atoms))
    print(
        "Composition:",
        report["system"]["composition"],
    )
    print("Charge:", EXPECTED_CHARGE)
    print(
        "Multiplicity:",
        EXPECTED_MULTIPLICITY,
    )
    print(
        "Fixed atoms:",
        len(fixed_indices),
    )
    print(
        "Mobile atoms:",
        len(mobile_indices),
    )
    print(
        "Fixed indices:",
        fixed_indices,
    )

    print()
    print("SCF strategy: TRAH")
    print("Initial guess: MORead")
    print(
        "Orbital source:",
        OUTPUT_GBW,
    )
    print(
        "GBW copy integrity: PASS"
    )
    print(
        "Coordinate mode: EMBEDDED_XYZ"
    )

    print()
    print("Input:", OUTPUT_INPUT)
    print("XYZ:", OUTPUT_XYZ)
    print("GBW:", OUTPUT_GBW)
    print("Map:", OUTPUT_MAP)
    print("Edges:", OUTPUT_EDGES)
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
        "ORCA restart execution authorized: False"
    )
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
