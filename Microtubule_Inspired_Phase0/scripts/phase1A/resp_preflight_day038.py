#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json

from resp_common import (
    read_orca_vpot,
    read_amber_esp,
)

ROOT = Path(__file__).resolve().parents[2]

VPOT = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "esp_upper_v7a_r1_20260731T174832Z"
    / "esp_upper_v7a_r1.vpot"
)

ESP = (
    ROOT
    / "runs/phase1A/day038_resp_generation"
    / "candidate_from_orca_vpot.esp"
)

REPORT = (
    ROOT
    / "runs/phase1A/day038_resp_generation"
    / "DAY038_RESP_PREFLIGHT.json"
)

AUTHORIZED_VPOT_SHA = (
    "73df47796c29d0b2a88a03d83efcb723d4f6e583d2e1789e1ea1a43d82c1064d"
)


def sha256(path: Path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        while True:

            block = f.read(1024 * 1024)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


print("=" * 100)
print("DAY038 / D038-C3")
print("RESP PREFLIGHT")
print("=" * 100)

vpot = read_orca_vpot(VPOT)
esp = read_amber_esp(ESP)

report = {

    "decision":
        "D038_RESP_PREFLIGHT_PASS",

    "source_execution":
        "day036_qm_f06_upper_v7a_r1_esp_executions",

    "vpot":{

        "path":str(VPOT),
        "sha256":sha256(VPOT),
        "authorized_sha256":AUTHORIZED_VPOT_SHA,
        "sha_match":
            sha256(VPOT)==AUTHORIZED_VPOT_SHA,

        "atom_count":vpot.atom_count,
        "grid_points":vpot.grid_point_count,
    },

    "amber_esp":{

        "path":str(ESP),
        "sha256":sha256(ESP),
        "atom_count":esp["natoms"],
        "grid_points":esp["npoints"],
    },

    "gates":{

        "authorized_vpot":(
            sha256(VPOT)==AUTHORIZED_VPOT_SHA
        ),

        "atom_count_52":
            vpot.atom_count==52,

        "grid_point_count_24835":
            vpot.grid_point_count==24835,

        "amber_matches_vpot":

            (
                esp["natoms"]==vpot.atom_count
            )

            and

            (
                esp["npoints"]==vpot.grid_point_count
            ),
    }
}

REPORT.write_text(
    json.dumps(report,indent=2)
)

print()

print(json.dumps(report,indent=2))

print()

print("Report written to")

print(REPORT)
