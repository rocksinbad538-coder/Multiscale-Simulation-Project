#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import csv
import hashlib
import json
import math
import re


ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_orca_input"
)

INPUT_PATH = INPUT_DIR / "v7a_r1.inp"
XYZ_PATH = INPUT_DIR / "v7a_r1_start.xyz"
GBW_PATH = INPUT_DIR / "v7a_r1_source.gbw"

MAP_PATH = (
    INPUT_DIR
    / "QM_F06_UPPER_V7A_R1_constraint_map.csv"
)

EDGES_PATH = (
    INPUT_DIR
    / "QM_F06_UPPER_V7A_R1_nominal_edges.csv"
)

REPORT_PATH = (
    INPUT_DIR
    / "QM_F06_UPPER_V7A_R1_ORCA_INPUT_REPORT.json"
)

EXPECTED_REPORT_DECISION = (
    "QM_F06_UPPER_V7A_R1_ORCA_INPUT_"
    "PREPARED_EXECUTION_AUDIT_REQUIRED"
)

EXPECTED_FIXED_INDICES = [
    0, 1, 3, 4, 5, 6,
    7, 8, 9, 10, 11, 14,
]

EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_ATOM_COUNT = 52
EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1
COORDINATE_TOLERANCE_A = 5.0e-10


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
        errors="strict",
    ).splitlines()

    atom_count = int(lines[0].strip())
    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ record at index {index}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"Incomplete XYZ file: {path}"
        )

    return atoms


def read_embedded_xyz(text: str) -> tuple[int, int, list[dict]]:
    match = re.search(
        r"(?ims)"
        r"^[ \t]*\*[ \t]*xyz[ \t]+"
        r"(-?\d+)[ \t]+(\d+)[ \t]*$"
        r"(.*?)"
        r"^[ \t]*\*[ \t]*$",
        text,
    )

    if match is None:
        raise RuntimeError(
            "Embedded XYZ block not found"
        )

    charge = int(match.group(1))
    multiplicity = int(match.group(2))

    atom_lines = [
        line
        for line in match.group(3).splitlines()
        if line.strip()
    ]

    atoms = []

    for index, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Invalid embedded XYZ record: {line}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    return charge, multiplicity, atoms


def coordinate_difference(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def parse_constraint_indices(text: str) -> list[int]:
    """
    Parse Cartesian constraints from the nested ORCA
    %geom / Constraints blocks.

    This line-oriented parser deliberately avoids a flat
    regular expression because the Constraints block has
    its own `end` nested inside the %geom block.
    """

    lines = text.splitlines()

    in_geom = False
    in_constraints = False
    geom_found = False
    constraints_found = False

    indices: list[int] = []

    geom_pattern = re.compile(
        r"^\s*%geom\b",
        flags=re.IGNORECASE,
    )

    constraints_pattern = re.compile(
        r"^\s*Constraints\b",
        flags=re.IGNORECASE,
    )

    constraint_record_pattern = re.compile(
        r"^\s*\{\s*C\s+(\d+)\s+C\s*\}\s*$",
        flags=re.IGNORECASE,
    )

    for line in lines:
        stripped = line.strip()

        if not in_geom:
            if geom_pattern.match(line):
                in_geom = True
                geom_found = True

            continue

        if not in_constraints:
            if constraints_pattern.match(line):
                in_constraints = True
                constraints_found = True
                continue

            # This is the outer %geom closure only when
            # we are not currently inside Constraints.
            if stripped.lower() == "end":
                break

            continue

        # Inside the Constraints block.
        if stripped.lower() == "end":
            in_constraints = False
            continue

        match = constraint_record_pattern.match(line)

        if match is not None:
            indices.append(
                int(match.group(1))
            )

    if not geom_found:
        raise RuntimeError(
            "%geom block not found"
        )

    if not constraints_found:
        raise RuntimeError(
            "Constraints block not found"
        )

    if not indices:
        raise RuntimeError(
            "Constraints block found, but no Cartesian "
            "constraint records were parsed"
        )

    if len(indices) != len(set(indices)):
        raise RuntimeError(
            "Duplicate Cartesian constraint indices found: "
            f"{indices}"
        )

    return indices


def count_blocks(
    text: str,
    pattern: str,
) -> int:
    return len(
        re.findall(
            pattern,
            text,
            flags=(
                re.IGNORECASE
                | re.MULTILINE
                | re.DOTALL
            ),
        )
    )


def main() -> None:
    required_paths = (
        INPUT_PATH,
        XYZ_PATH,
        GBW_PATH,
        MAP_PATH,
        EDGES_PATH,
        REPORT_PATH,
    )

    for path in required_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing or empty package artifact: {path}"
            )

    text = INPUT_PATH.read_text(
        encoding="utf-8",
        errors="strict",
    )

    report = json.loads(
        REPORT_PATH.read_text(
            encoding="utf-8"
        )
    )

    xyz_atoms = read_xyz(XYZ_PATH)

    charge, multiplicity, embedded_atoms = (
        read_embedded_xyz(text)
    )

    with MAP_PATH.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        map_rows = list(csv.DictReader(handle))

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    with EDGES_PATH.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        edge_rows = list(csv.DictReader(handle))

    fixed_indices_from_map = [
        int(row["v7a_index_0based"])
        for row in map_rows
        if (
            row["v7a_fixed"]
            .strip()
            .lower()
            == "true"
        )
    ]

    mobile_indices_from_map = [
        int(row["v7a_index_0based"])
        for row in map_rows
        if (
            row["v7a_mobile"]
            .strip()
            .lower()
            == "true"
        )
    ]

    constraint_indices = (
        parse_constraint_indices(text)
    )

    composition: dict[str, int] = {}

    for atom in embedded_atoms:
        composition[atom["element"]] = (
            composition.get(atom["element"], 0)
            + 1
        )

    maximum_embedded_xyz_difference = 0.0
    element_order_failures = []

    if len(embedded_atoms) == len(xyz_atoms):
        for index, (
            embedded_atom,
            xyz_atom,
        ) in enumerate(
            zip(embedded_atoms, xyz_atoms)
        ):
            map_element = (
                map_rows[index]["element"]
                if index < len(map_rows)
                else None
            )

            if not (
                embedded_atom["element"]
                == xyz_atom["element"]
                == map_element
            ):
                element_order_failures.append({
                    "index_0based": index,
                    "embedded_element": (
                        embedded_atom["element"]
                    ),
                    "xyz_element": (
                        xyz_atom["element"]
                    ),
                    "map_element": map_element,
                })

            difference = coordinate_difference(
                embedded_atom["xyz_A"],
                xyz_atom["xyz_A"],
            )

            maximum_embedded_xyz_difference = max(
                maximum_embedded_xyz_difference,
                difference,
            )
    else:
        maximum_embedded_xyz_difference = (
            float("inf")
        )

    simple_lines = [
        line.strip()
        for line in text.splitlines()
        if line.lstrip().startswith("!")
    ]

    simple_line = (
        simple_lines[0]
        if len(simple_lines) == 1
        else ""
    )

    required_protocol_tokens = (
        "PBE0",
        "D4",
        "def2-TZVP",
        "def2/J",
        "RIJCOSX",
        "TightSCF",
        "Opt",
        "DefGrid3",
        "TRAH",
        "MORead",
    )

    scf_block_match = re.search(
        r"(?ims)"
        r"^[ \t]*%scf\b"
        r"(.*?)"
        r"^[ \t]*end[ \t]*$",
        text,
    )

    scf_block = (
        scf_block_match.group(1)
        if scf_block_match
        else ""
    )

    moinp_matches = re.findall(
        r'(?im)'
        r'^[ \t]*%moinp[ \t]+'
        r'"([^"]+)"[ \t]*$',
        text,
    )

    pal_match = re.search(
        r"(?ims)"
        r"^[ \t]*%pal\b"
        r"(.*?)"
        r"^[ \t]*end[ \t]*$",
        text,
    )

    pal_block = (
        pal_match.group(1)
        if pal_match
        else ""
    )

    maxcore_match = re.search(
        r"(?im)"
        r"^[ \t]*%maxcore[ \t]+(\d+)[ \t]*$",
        text,
    )

    maximum_iterations_match = re.search(
        r"(?im)"
        r"^[ \t]*MaxIter[ \t]+(\d+)[ \t]*$",
        scf_block,
    )

    gbw_hash = sha256(GBW_PATH)
    input_hash = sha256(INPUT_PATH)
    xyz_hash = sha256(XYZ_PATH)
    map_hash = sha256(MAP_PATH)
    edges_hash = sha256(EDGES_PATH)

    reported_hashes = report.get(
        "sha256",
        {},
    )

    gates = {
        "input_preparation_decision": (
            report.get("decision")
            == EXPECTED_REPORT_DECISION
        ),
        "source_execution_marked_immutable": (
            report.get(
                "source_execution_immutable"
            )
            is True
        ),
        "embedded_xyz_block_present": (
            len(embedded_atoms)
            == EXPECTED_ATOM_COUNT
        ),
        "charge_zero": (
            charge == EXPECTED_CHARGE
        ),
        "singlet_multiplicity": (
            multiplicity
            == EXPECTED_MULTIPLICITY
        ),
        "external_xyz_atom_count": (
            len(xyz_atoms)
            == EXPECTED_ATOM_COUNT
        ),
        "constraint_map_atom_count": (
            len(map_rows)
            == EXPECTED_ATOM_COUNT
        ),
        "nominal_edge_count_57": (
            len(edge_rows) == 57
        ),
        "composition_B17_N14_H21": (
            composition
            == EXPECTED_COMPOSITION
        ),
        "atom_identity_and_order": (
            len(element_order_failures) == 0
        ),
        "embedded_coordinates_match_xyz": (
            maximum_embedded_xyz_difference
            <= COORDINATE_TOLERANCE_A
        ),
        "exactly_one_simple_keyword_line": (
            len(simple_lines) == 1
        ),
        "protocol_PBE0_D4_def2_TZVP": (
            all(
                token in simple_line
                for token
                in required_protocol_tokens
            )
        ),
        "exactly_one_scf_block": (
            count_blocks(
                text,
                (
                    r"^[ \t]*%scf\b.*?"
                    r"^[ \t]*end[ \t]*$"
                ),
            )
            == 1
        ),
        "scf_maxiter_500": (
            maximum_iterations_match
            is not None
            and int(
                maximum_iterations_match.group(1)
            )
            == 500
        ),
        "explicit_TRAH": (
            re.search(
                r"(?i)\bTRAH\b",
                simple_line,
            )
            is not None
        ),
        "explicit_MORead": (
            re.search(
                r"(?i)\bMORead\b",
                simple_line,
            )
            is not None
        ),
        "exactly_one_moinp": (
            moinp_matches
            == ["v7a_r1_source.gbw"]
        ),
        "moinp_file_present": (
            GBW_PATH.is_file()
            and GBW_PATH.stat().st_size > 0
        ),
        "gbw_hash_matches_report": (
            gbw_hash
            == reported_hashes.get(
                "copied_gbw"
            )
        ),
        "input_hash_matches_report": (
            input_hash
            == reported_hashes.get(
                "prepared_input"
            )
        ),
        "xyz_hash_matches_report": (
            xyz_hash
            == reported_hashes.get(
                "prepared_xyz"
            )
        ),
        "constraint_map_hash_matches_report": (
            map_hash
            == reported_hashes.get(
                "constraint_map"
            )
        ),
        "nominal_edges_hash_matches_report": (
            edges_hash
            == reported_hashes.get(
                "nominal_edges"
            )
        ),
        "constraints_match_map": (
            constraint_indices
            == fixed_indices_from_map
            == EXPECTED_FIXED_INDICES
        ),
        "fixed_mobile_partition_complete": (
            sorted(
                fixed_indices_from_map
                + mobile_indices_from_map
            )
            == list(
                range(EXPECTED_ATOM_COUNT)
            )
            and not (
                set(fixed_indices_from_map)
                & set(mobile_indices_from_map)
            )
        ),
        "parallel_configuration_nprocs_4": (
            re.search(
                r"(?im)"
                r"^[ \t]*nprocs[ \t]+4[ \t]*$",
                pal_block,
            )
            is not None
        ),
        "memory_configuration_2500_MB": (
            maxcore_match is not None
            and int(maxcore_match.group(1))
            == 2500
        ),
        "no_xyzfile_directive": (
            "xyzfile"
            not in text.lower()
        ),
        "no_geometry_restart_directive": (
            "restart" not in simple_line.lower()
        ),
    }

    passed = all(gates.values())

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "ORCA EXECUTION AUDIT"
    )
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:62s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Charge:", charge)
    print("Multiplicity:", multiplicity)
    print(
        "Embedded atoms:",
        len(embedded_atoms),
    )
    print("XYZ atoms:", len(xyz_atoms))
    print("Map atoms:", len(map_rows))
    print("Nominal edges:", len(edge_rows))
    print("Composition:", composition)
    print(
        "Fixed atoms:",
        len(fixed_indices_from_map),
    )
    print(
        "Mobile atoms:",
        len(mobile_indices_from_map),
    )
    print(
        "Constraint indices:",
        constraint_indices,
    )
    print(
        "Maximum embedded/XYZ difference A:",
        maximum_embedded_xyz_difference,
    )
    print("MOInp records:", moinp_matches)
    print("GBW size bytes:", GBW_PATH.stat().st_size)

    if element_order_failures:
        print()
        print("Element-order failures:")
        for failure in element_order_failures:
            print("  ", failure)

    print()
    print(
        "Decision:",
        (
            "QM_F06_UPPER_V7A_R1_"
            "ORCA_EXECUTION_GATE_PASS"
            if passed
            else
            "QM_F06_UPPER_V7A_R1_"
            "ORCA_EXECUTION_GATE_FAIL"
        ),
    )

    print(
        "ORCA restart execution authorized:",
        passed,
    )
    print("RESP authorized: False")
    print(
        "Force-field adoption authorized: False"
    )
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
