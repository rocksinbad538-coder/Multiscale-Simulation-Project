#!/usr/bin/env python3
"""
Pair-specific steric audit for capped QM_F06 fragments.

The previous universal 0.70 Å threshold only detected near-overlapping atoms.
This audit evaluates each nonbonded pair relative to the sum of elemental
van der Waals radii.

Classification:
    ratio = distance / (vdW_radius_i + vdW_radius_j)

    ratio < 0.70  -> SEVERE_CLASH
    ratio < 0.80  -> STRONG_COMPRESSION
    ratio < 0.90  -> CLOSE_CONTACT
    otherwise     -> ACCEPTABLE

No coordinates are modified.
No QM calculation is executed or authorized.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")

# Bondi-type van der Waals radii, in angstrom.
VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def coordinates(
    atom: dict[str, str],
) -> tuple[float, float, float]:
    return (
        float(atom["x_angstrom"]),
        float(atom["y_angstrom"]),
        float(atom["z_angstrom"]),
    )


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second, strict=True)
        )
    )


def canonical_edge(
    atom_1: str,
    atom_2: str,
) -> tuple[str, str]:
    return tuple(sorted((atom_1, atom_2)))


def classify_ratio(ratio: float) -> str:
    if ratio < 0.70:
        return "SEVERE_CLASH"

    if ratio < 0.80:
        return "STRONG_COMPRESSION"

    if ratio < 0.90:
        return "CLOSE_CONTACT"

    return "ACCEPTABLE"


def main() -> None:
    combined_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []
    fragment_passes: list[bool] = []

    for end in ENDS:
        label = f"QM_F06_{end}_CAPPED"

        atoms = read_csv(
            F06_DIR / f"{label}_atoms.csv"
        )

        internal_edges = read_csv(
            F06_DIR / f"QM_F06_{end}_internal_edges.csv"
        )

        boundary_audit = read_csv(
            F06_DIR / f"QM_F06_{end}_boundary_edge_audit.csv"
        )

        caps = read_csv(
            F06_DIR / f"{label}_caps.csv"
        )

        atom_lookup = {
            atom["atom_id"]: atom
            for atom in atoms
        }

        bonded_pairs: set[tuple[str, str]] = set()

        for edge in internal_edges:
            bonded_pairs.add(
                canonical_edge(
                    edge["source_node"],
                    edge["target_node"],
                )
            )

        for edge in boundary_audit:
            if (
                edge["preliminary_action"]
                == "INCLUDE_EXISTING_HYDROGEN"
            ):
                bonded_pairs.add(
                    canonical_edge(
                        edge["inside_node"],
                        edge["outside_node"],
                    )
                )

        for cap in caps:
            bonded_pairs.add(
                canonical_edge(
                    cap["parent_inside_node"],
                    cap["cap_id"],
                )
            )

        atom_ids = sorted(atom_lookup)
        audit_rows: list[dict[str, Any]] = []

        for index, atom_1_id in enumerate(atom_ids):
            for atom_2_id in atom_ids[index + 1:]:
                pair = canonical_edge(
                    atom_1_id,
                    atom_2_id,
                )

                if pair in bonded_pairs:
                    continue

                atom_1 = atom_lookup[atom_1_id]
                atom_2 = atom_lookup[atom_2_id]

                element_1 = atom_1["element"]
                element_2 = atom_2["element"]

                measured = distance(
                    coordinates(atom_1),
                    coordinates(atom_2),
                )

                vdw_sum = (
                    VDW_RADII[element_1]
                    + VDW_RADII[element_2]
                )

                ratio = measured / vdw_sum
                classification = classify_ratio(ratio)

                involves_artificial_cap = (
                    str(atom_1["artificial_cap"]).lower()
                    == "true"
                    or str(atom_2["artificial_cap"]).lower()
                    == "true"
                )

                involves_existing_h = (
                    atom_1["atom_role"]
                    == "EXISTING_R2_HYDROGEN_ADDED"
                    or atom_2["atom_role"]
                    == "EXISTING_R2_HYDROGEN_ADDED"
                )

                # Preserve all contacts below the vdW sum.
                if ratio < 1.0:
                    row = {
                        "fragment": label,
                        "atom_1": atom_1_id,
                        "element_1": element_1,
                        "role_1": atom_1["atom_role"],
                        "atom_2": atom_2_id,
                        "element_2": element_2,
                        "role_2": atom_2["atom_role"],
                        "distance_angstrom": f"{measured:.10f}",
                        "vdw_sum_angstrom": f"{vdw_sum:.10f}",
                        "distance_over_vdw_sum": f"{ratio:.10f}",
                        "classification": classification,
                        "involves_artificial_cap": (
                            involves_artificial_cap
                        ),
                        "involves_existing_hydrogen": (
                            involves_existing_h
                        ),
                        "requires_geometry_repair": (
                            classification
                            in {
                                "SEVERE_CLASH",
                                "STRONG_COMPRESSION",
                            }
                        ),
                    }

                    audit_rows.append(row)
                    combined_rows.append(row)

        class_counts = Counter(
            row["classification"]
            for row in audit_rows
        )

        severe_or_strong = sum(
            class_counts[classification]
            for classification in (
                "SEVERE_CLASH",
                "STRONG_COMPRESSION",
            )
        )

        fragment_pass = severe_or_strong == 0
        fragment_passes.append(fragment_pass)

        write_csv(
            F06_DIR
            / f"{label}_pair_specific_steric_audit.csv",
            audit_rows,
        )

        minimum_row = min(
            audit_rows,
            key=lambda row: float(
                row["distance_over_vdw_sum"]
            ),
        )

        report_sections.extend(
            [
                f"## {label}",
                "",
                (
                    f"- Nonbonded contacts below vdW sum: "
                    f"**{len(audit_rows)}**"
                ),
                (
                    f"- Classification counts: "
                    f"`{dict(sorted(class_counts.items()))}`"
                ),
                (
                    f"- Most compressed contact: "
                    f"`{minimum_row['atom_1']} — "
                    f"{minimum_row['atom_2']}`"
                ),
                (
                    f"- Distance: "
                    f"**{float(minimum_row['distance_angstrom']):.6f} Å**"
                ),
                (
                    f"- Distance/vdW-sum ratio: "
                    f"**{float(minimum_row['distance_over_vdw_sum']):.4f}**"
                ),
                (
                    f"- Pair-specific steric gate: "
                    f"**{'PASS' if fragment_pass else 'FAIL'}**"
                ),
                "",
            ]
        )

    overall_pass = all(fragment_passes)

    decision = (
        "QM_F06_CAPPED_FRAGMENTS_PASS_PAIR_SPECIFIC_STERIC_GATE"
        if overall_pass
        else
        "QM_F06_CAPPED_FRAGMENTS_REQUIRE_CAP_GEOMETRY_REPAIR"
    )

    write_csv(
        F06_DIR
        / "QM_F06_pair_specific_steric_audit_combined.csv",
        combined_rows,
    )

    (
        F06_DIR
        / "QM_F06_PAIR_SPECIFIC_STERIC_AUDIT.md"
    ).write_text(
        "\n".join(
            [
                "# QM_F06 Pair-Specific Steric Audit — Day026",
                "",
                "## Rationale",
                "",
                (
                    "The previous 0.70 Å universal threshold only "
                    "excluded near-coincident atoms. The present audit "
                    "uses interatomic distances normalized by the sum "
                    "of elemental van der Waals radii."
                ),
                "",
                f"## Decision: **{decision}**",
                "",
                *report_sections,
                "## Interpretation",
                "",
                (
                    "SEVERE_CLASH and STRONG_COMPRESSION contacts must "
                    "be corrected before electronic-structure input "
                    "preparation. CLOSE_CONTACT geometries may be retained "
                    "as initial structures but must be monitored during "
                    "optimization."
                ),
                "",
                "## Authorization state",
                "",
                "- Pair-specific steric audit: **COMPLETED**",
                (
                    "- QM input preparation: "
                    f"**{'AUTHORIZED' if overall_pass else 'NOT AUTHORIZED'}**"
                ),
                "- QM calculation executed: **NO**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "pair_specific_steric_gate_pass": overall_pass,
        "qm_input_preparation_authorized": overall_pass,
        "qm_calculation_executed": False,
        "required_next_step": (
            "PREPARE_QM_INPUTS"
            if overall_pass
            else "REPAIR_ARTIFICIAL_CAP_GEOMETRIES"
        ),
    }

    (
        F06_DIR
        / "QM_F06_pair_specific_steric_audit_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Pair-specific QM_F06 steric audit completed.")
    print(f"Decision: {decision}")
    print(
        "QM input preparation authorized:",
        overall_pass,
    )


if __name__ == "__main__":
    main()
