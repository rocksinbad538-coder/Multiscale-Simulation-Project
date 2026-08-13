#!/usr/bin/env python3

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]

protocol = {
    "current_outputs": [
        "coordinates",
        "box_dimensions",
        "potential_energy",
        "temperature",
        "pressure",
        "Rg",
        "RMSD",
        "RMSF",
        "shape_descriptors"
    ],

    "candidate_future_outputs": [
        {
            "observable": "velocities",
            "required_for": [
                "VACF",
                "diffusion",
                "kinetic analyses"
            ],
            "status": "REVIEW_REQUIRED"
        },
        {
            "observable": "forces",
            "required_for": [
                "mechanical response",
                "local force analysis"
            ],
            "status": "REVIEW_REQUIRED"
        },
        {
            "observable": "per_atom_energy",
            "required_for": [
                "energy localization"
            ],
            "status": "REVIEW_REQUIRED"
        },
        {
            "observable": "stress_tensor",
            "required_for": [
                "elastic response"
            ],
            "status": "REVIEW_REQUIRED"
        },
        {
            "observable": "dipole_related_quantities",
            "required_for": [
                "confined water",
                "dielectric response",
                "excitonic model"
            ],
            "status": "BLOCKED_UNTIL_MODEL_REVIEW"
        }
    ],

    "decision":
        "DO_NOT_REPEAT_MD_UNTIL_FINAL_PROTOCOL_IS_APPROVED"
}

out = ROOT / "runs" / "phase2" / "FINAL_MD_PROTOCOL_REVIEW.json"

out.write_text(
    json.dumps(protocol, indent=2)
)

print("="*90)
print("DAY046 / PHASE2-A31")
print("FINAL MD PROTOCOL REVIEW")
print("="*90)
print(out)
