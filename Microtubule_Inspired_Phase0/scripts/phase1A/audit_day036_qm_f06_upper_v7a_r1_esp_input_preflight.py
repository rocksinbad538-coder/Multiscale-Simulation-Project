#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import re


ROOT = Path(__file__).resolve().parents[2]

DAY036_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_r1_resp_preparation"
)

PREPARATION_REPORT = (
    DAY036_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_PREPARATION.json"
)

EXECUTION_MANIFEST = (
    DAY036_DIR
    / "QM_F06_UPPER_V7A_R1_RESP_EXECUTION_MANIFEST.json"
)

ESP_INPUT = (
    DAY036_DIR
    / "QM_F06_UPPER_V7A_R1_ESP_ORCA.inp"
)

ADOPTION_REPORT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_coordinate_adoption/"
    / "QM_F06_UPPER_V7A_COORDINATE_ADOPTION.json"
)

OUTPUT_REPORT = (
    DAY036_DIR
    / "QM_F06_UPPER_V7A_R1_ESP_INPUT_PREFLIGHT.json"
)

EXPECTED_PREPARATION_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "RESP_PREPARATION_PASS_"
    "ESP_INPUT_PREFLIGHT_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1

COORDINATE_TOLERANCE_A = 5.0e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(
            f"Required JSON is missing: {path}"
        )

    return json.loads(
        path.read_text(encoding="utf-8")
    )


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    atom_count = int(lines[0].strip())
    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ record {index}: {path}"
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
            f"XYZ atom-count mismatch: {path}"
        )

    return atoms


def read_orca_xyz_block(
    path: Path,
) -> tuple[int, int, list[dict]]:
    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    header_index = None
    charge = None
    multiplicity = None

    pattern = re.compile(
        r"^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$",
        re.IGNORECASE,
    )

    for index, line in enumerate(lines):
        match = pattern.match(line)

        if match:
            header_index = index
            charge = int(match.group(1))
            multiplicity = int(match.group(2))
            break

    if header_index is None:
        raise RuntimeError(
            "No '* xyz charge multiplicity' block found."
        )

    atoms = []

    for line in lines[header_index + 1:]:
        stripped = line.strip()

        if stripped == "*":
            break

        if not stripped:
            continue

        fields = stripped.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Invalid ORCA coordinate row: {line}"
            )

        atoms.append({
            "index_0based": len(atoms),
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    return charge, multiplicity, atoms


def maximum_coordinate_difference(
    first: list[dict],
    second: list[dict],
) -> float:
    if len(first) != len(second):
        return float("inf")

    maximum = 0.0

    for atom_a, atom_b in zip(first, second):
        if atom_a["element"] != atom_b["element"]:
            return float("inf")

        difference = sum(
            (a - b) ** 2
            for a, b in zip(
                atom_a["xyz_A"],
                atom_b["xyz_A"],
            )
        ) ** 0.5

        maximum = max(
            maximum,
            difference,
        )

    return maximum


def main() -> None:
    preparation = read_json(
        PREPARATION_REPORT
    )

    manifest = read_json(
        EXECUTION_MANIFEST
    )

    adoption = read_json(
        ADOPTION_REPORT
    )

    adopted_xyz = Path(
        adoption["outputs"]["adopted_final_XYZ"]
    )

    if not adopted_xyz.is_file():
        raise RuntimeError(
            f"Adopted XYZ is missing: {adopted_xyz}"
        )

    orca_executable = Path(
        manifest["candidate_ORCA_executable"]
    )

    esp_text = ESP_INPUT.read_text(
        encoding="utf-8"
    )

    charge, multiplicity, esp_atoms = (
        read_orca_xyz_block(ESP_INPUT)
    )

    adopted_atoms = read_xyz(
        adopted_xyz
    )

    maximum_difference = (
        maximum_coordinate_difference(
            adopted_atoms,
            esp_atoms,
        )
    )

    composition = Counter(
        atom["element"]
        for atom in esp_atoms
    )

    gates = {
        "preparation_decision_matches": (
            preparation.get("decision")
            == EXPECTED_PREPARATION_DECISION
        ),
        "ESP_input_preflight_authorized_upstream": (
            preparation.get(
                "authorizations",
                {},
            ).get(
                "ESP_input_preflight_authorized"
            )
            is True
        ),
        "ESP_execution_blocked_upstream": (
            preparation.get(
                "authorizations",
                {},
            ).get(
                "ESP_execution_authorized"
            )
            is False
        ),
        "manifest_execution_blocked": (
            manifest.get(
                "execution_authorized"
            )
            is False
        ),
        "ORCA_executable_exists": (
            orca_executable.is_file()
        ),
        "ORCA_executable_is_executable": (
            orca_executable.is_file()
            and os.access(
                orca_executable,
                os.X_OK,
            )
        ),
        "ESP_input_exists": (
            ESP_INPUT.is_file()
            and ESP_INPUT.stat().st_size > 0
        ),
        "ESP_input_hash_matches_manifest": (
            sha256(ESP_INPUT)
            == manifest.get("input_sha256")
        ),
        "adopted_geometry_hash_matches_manifest": (
            sha256(adopted_xyz)
            == manifest.get("geometry_sha256")
        ),
        "method_PBE0_D4": (
            "PBE0 D4" in esp_text
        ),
        "basis_def2_TZVP": (
            "def2-TZVP" in esp_text
        ),
        "RIJCOSX_present": (
            "RIJCOSX" in esp_text
        ),
        "TightSCF_present": (
            "TightSCF" in esp_text
        ),
        "DefGrid3_present": (
            "DefGrid3" in esp_text
        ),
        "CHELPG_keyword_present": (
            "CHELPG" in esp_text
        ),
        "CHELPG_block_present": (
            "%chelpg" in esp_text.lower()
        ),
        "net_charge_zero": (
            charge == EXPECTED_CHARGE
        ),
        "multiplicity_one": (
            multiplicity
            == EXPECTED_MULTIPLICITY
        ),
        "ESP_input_atom_count_52": (
            len(esp_atoms)
            == EXPECTED_ATOM_COUNT
        ),
        "ESP_input_composition_B17_N14_H21": (
            dict(composition)
            == EXPECTED_COMPOSITION
        ),
        "ESP_atom_order_matches_adopted_geometry": (
            len(esp_atoms)
            == len(adopted_atoms)
            and all(
                atom_a["element"]
                == atom_b["element"]
                for atom_a, atom_b in zip(
                    esp_atoms,
                    adopted_atoms,
                )
            )
        ),
        "ESP_coordinates_match_adopted_geometry": (
            maximum_difference
            <= COORDINATE_TOLERANCE_A
        ),
        "input_isolated_in_day036_directory": (
            ESP_INPUT.parent.resolve()
            == DAY036_DIR.resolve()
        ),
    }

    passed = all(gates.values())

    decision = (
        "QM_F06_UPPER_V7A_R1_"
        "ESP_INPUT_PREFLIGHT_PASS_"
        "ESP_EXECUTION_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V7A_R1_"
        "ESP_INPUT_PREFLIGHT_FAIL_"
        "ESP_EXECUTION_BLOCKED"
    )

    report = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "model": "QM_F06_UPPER_V7A_R1",
        "stage": "DAY036_ESP_INPUT_PREFLIGHT",
        "decision": decision,
        "gates": gates,
        "summary": {
            "ORCA_executable": str(
                orca_executable
            ),
            "ESP_input": str(
                ESP_INPUT
            ),
            "adopted_geometry": str(
                adopted_xyz
            ),
            "net_charge": charge,
            "multiplicity": multiplicity,
            "atom_count": len(esp_atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "maximum_adopted_vs_ESP_coordinate_difference_A": (
                maximum_difference
            ),
            "coordinate_tolerance_A": (
                COORDINATE_TOLERANCE_A
            ),
        },
        "authorizations": {
            "ESP_execution_authorized": passed,
            "RESP_input_generation_authorized": False,
            "RESP_execution_authorized": False,
            "RESP_validation_authorized": False,
            "charge_adoption_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 104)
    print("QM_F06 UPPER V7-A R1 ESP INPUT PREFLIGHT")
    print("=" * 104)

    for name, value in gates.items():
        print(
            f"{name:60s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(
        "Maximum adopted/input coordinate difference A:",
        maximum_difference,
    )
    print("Decision:", decision)
    print("Report:", OUTPUT_REPORT)
    print(
        "ESP execution authorized:",
        passed,
    )
    print("RESP execution authorized: False")
    print("Force-field adoption authorized: False")
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
