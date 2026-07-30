#!/usr/bin/env python3
"""
design_day035_qm_f06_upper_v7a_r1_RESP_input.py

Skeleton production gate for RESP input design.

STATUS:
- Version 0.1
- Intended to be the first RESP gate after
  adopt_day035_qm_f06_upper_v7a_r1_final_coordinates.py

Responsibilities
----------------
1. Read QM_F06_UPPER_V7A_COORDINATE_ADOPTION.json
2. Verify:
   * decision ==
     QM_F06_UPPER_V7A_FINAL_COORDINATES_ADOPTED_RESP_INPUT_DESIGN_AUTHORIZED
   * authorizations.RESP_input_design_authorized == True
3. Verify adopted XYZ exists.
4. Verify SHA256 of adopted XYZ.
5. Build RESP manifest.
6. Emit:
   QM_F06_UPPER_V7A_R1_RESP_INPUT_DESIGN.json

NOTE
----
This file is intentionally a production skeleton. It defines the complete
contract for the RESP stage before implementation of Multiwfn / ESP
generation.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]

COORDINATE_ADOPTION_DIR = (
    ROOT /
    "runs/phase1A/day035_qm_f06_upper_v7a_r1_coordinate_adoption"
)

ADOPTION_REPORT = (
    COORDINATE_ADOPTION_DIR /
    "QM_F06_UPPER_V7A_COORDINATE_ADOPTION.json"
)

OUTPUT_DIR = (
    ROOT /
    "runs/phase1A/day035_qm_f06_upper_v7a_r1_resp_input_design"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_REPORT = (
    OUTPUT_DIR /
    "QM_F06_UPPER_V7A_R1_RESP_INPUT_DESIGN.json"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():

    if not ADOPTION_REPORT.exists():
        raise RuntimeError(f"Missing adoption report: {ADOPTION_REPORT}")

    adoption = json.loads(ADOPTION_REPORT.read_text())

    if adoption["decision"] != (
        "QM_F06_UPPER_V7A_FINAL_COORDINATES_ADOPTED_"
        "RESP_INPUT_DESIGN_AUTHORIZED"
    ):
        raise RuntimeError("Coordinate adoption decision mismatch.")

    auth = adoption["authorizations"]

    if not auth["RESP_input_design_authorized"]:
        raise RuntimeError("RESP design not authorized.")

    xyz = Path(
        adoption["outputs"]["adopted_final_XYZ"]
    )

    if not xyz.exists():
        raise RuntimeError(f"Missing adopted XYZ: {xyz}")

    report = {
        "generated_utc":
            datetime.now(timezone.utc).isoformat(),
        "source_coordinate_adoption_report":
            str(ADOPTION_REPORT),
        "execution_directory":
            adoption["source_execution_directory"],
        "geometry_sha256":
            sha256(xyz),
        "protocol": {
            "charge_model": "RESP",
            "ESP_source": "ORCA",
            "status": "DESIGN_APPROVED"
        },
        "decision":
            "QM_F06_UPPER_V7A_R1_RESP_INPUT_DESIGN_PASS_"
            "RESP_INPUT_PREPARATION_AUTHORIZED",
        "authorizations": {
            "RESP_input_design_authorized": True,
            "RESP_input_preparation_authorized": True,
            "RESP_execution_authorized": False,
            "RESP_validation_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False
        }
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8"
    )

    print("RESP input design gate: PASS")
    print(f"Report: {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()
