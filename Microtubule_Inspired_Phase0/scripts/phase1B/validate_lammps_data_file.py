#!/usr/bin/env python3

from pathlib import Path
import json
import datetime

ROOT = Path(__file__).resolve().parents[2]

DATA = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_export"
    / "data.lammps"
)

MODEL = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_model"
    / "PHASE1B_LAMMPS_MODEL.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_data_validation"
)

RUN.mkdir(parents=True, exist_ok=True)

assert DATA.exists()
assert MODEL.exists()

model = json.load(open(MODEL))

text = DATA.read_text()

required_sections = [
    "Masses",
    "Bond Coeffs",
    "Angle Coeffs",
    "Improper Coeffs",
    "Dihedral Coeffs",
    "Atoms",
    "Bonds",
    "Angles",
    "Impropers",
    "Dihedrals"
]

section_status = {
    s: (s in text)
    for s in required_sections
}

report = {

    "timestamp":
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z"),

    "atom_count":
        len(model["system"]["atoms"]),

    "bond_count":
        len(model["system"]["bonds"]),

    "angle_count":
        len(model["system"]["angles"]),

    "improper_count":
        len(model["system"]["impropers"]),

    "dihedral_count":
        len(model["system"]["dihedrals"]),

    "sections":
        section_status,

    "all_sections_present":
        all(section_status.values()),

    "decision":
        "PASS" if all(section_status.values()) else "FAIL"
}

outfile = RUN / "LAMMPS_DATA_FILE_VALIDATION.json"

outfile.write_text(
    json.dumps(report, indent=2)
)

print("="*100)
print("DAY041 / PHASE1B-A24")
print("LAMMPS DATA FILE VALIDATION")
print("="*100)
print()

for s,v in section_status.items():
    print(f"{s:18s} {v}")

print()

print("[OUTPUT]")
print(outfile)

print()

print("[DECISION]")
print(report["decision"])
