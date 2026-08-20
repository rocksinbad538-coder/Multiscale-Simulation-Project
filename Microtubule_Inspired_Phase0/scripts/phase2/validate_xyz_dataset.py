#!/usr/bin/env python3

from pathlib import Path
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[2]

XYZDIR = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
    / "tddft_xyz"
)

EXPECTED_NATOMS = 37
ALLOWED_ELEMENTS = {"B", "N", "H"}

xyz_files = sorted(
    XYZDIR.glob("*.xyz")
)

if not xyz_files:
    raise RuntimeError(
        f"No XYZ structures found in {XYZDIR}"
    )

summary = []
reference_composition = None
all_valid = True

for xyz in xyz_files:

    lines = xyz.read_text().splitlines()

    valid = True
    errors = []

    if len(lines) < 2:

        valid = False
        errors.append("file_too_short")

        natoms = -1
        atom_lines = []

    else:

        natoms = int(lines[0])
        atom_lines = lines[2:]

    composition = {}
    parsed_atoms = 0

    for line in atom_lines:

        if not line.strip():
            continue

        fields = line.split()

        if len(fields) != 4:

            valid = False
            errors.append(
                "invalid_xyz_record"
            )
            continue

        element = fields[0]

        if element not in ALLOWED_ELEMENTS:

            valid = False
            errors.append(
                f"invalid_element_{element}"
            )

        try:
            x, y, z = map(
                float,
                fields[1:4]
            )
        except ValueError:

            valid = False
            errors.append(
                "non_numeric_coordinate"
            )
            continue

        if not all(
            math.isfinite(v)
            for v in (x, y, z)
        ):

            valid = False
            errors.append(
                "non_finite_coordinate"
            )

        composition[element] = (
            composition.get(element, 0) + 1
        )

        parsed_atoms += 1

    if natoms != EXPECTED_NATOMS:
        valid = False
        errors.append(
            f"natoms_{natoms}_expected_{EXPECTED_NATOMS}"
        )

    if parsed_atoms != natoms:
        valid = False
        errors.append(
            f"parsed_{parsed_atoms}_header_{natoms}"
        )

    if sum(composition.values()) != EXPECTED_NATOMS:
        valid = False
        errors.append(
            "composition_count_mismatch"
        )

    if reference_composition is None:

        reference_composition = dict(
            sorted(composition.items())
        )

    elif dict(sorted(composition.items())) != reference_composition:

        valid = False
        errors.append(
            "composition_not_constant_across_dataset"
        )

    all_valid &= valid

    summary.append({
        "file": xyz.name,
        "natoms": natoms,
        "parsed_atoms": parsed_atoms,
        "composition": composition,
        "valid": valid,
        "errors": sorted(set(errors)),
    })


outfile = (
    XYZDIR
    / "XYZ_VALIDATION.json"
)

outfile.write_text(
    json.dumps(
        {
            "expected_natoms": EXPECTED_NATOMS,
            "reference_composition":
                reference_composition,
            "structures": summary,
            "dataset_valid": bool(all_valid),
        },
        indent=2,
    )
    + "\n"
)

print("="*90)
print("PHASE5-D51")
print("XYZ DATASET VALIDATION")
print("="*90)

for item in summary:

    print(
        f"{item['file']:35s} "
        f"atoms={item['natoms']:2d} "
        f"parsed={item['parsed_atoms']:2d} "
        f"{'PASS' if item['valid'] else 'FAIL'}"
    )

print()
print(
    "REFERENCE_COMPOSITION="
    f"{reference_composition}"
)

print(
    "XYZ_DATASET="
    f"{'PASS' if all_valid else 'FAIL'}"
)

print(outfile)

sys.exit(0 if all_valid else 1)
