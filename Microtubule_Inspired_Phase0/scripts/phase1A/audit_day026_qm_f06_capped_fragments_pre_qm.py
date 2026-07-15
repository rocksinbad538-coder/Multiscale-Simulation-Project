#!/usr/bin/env python3
"""
Pre-QM audit of capped QM_F06 LOWER and UPPER fragments.

Checks:
- atom and elemental counts;
- reconstructed fragment connectivity;
- expected graph degree (B/N=3, H=1);
- bonded-distance ranges;
- minimum nonbonded interatomic distance and steric clashes;
- valence-electron parity for candidate charge states;
- provisional charge/multiplicity recommendation.

This script does not execute or authorize a QM calculation.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")

VALENCE_ELECTRONS = {
    "H": 1,
    "B": 3,
    "N": 5,
}

EXPECTED_DEGREE = {
    "H": 1,
    "B": 3,
    "N": 3,
}

BONDED_DISTANCE_RANGES = {
    frozenset(("B", "N")): (1.20, 1.75),
    frozenset(("B", "H")): (0.90, 1.40),
    frozenset(("N", "H")): (0.80, 1.25),
}

# Conservative universal nonbonded clash gate.
NONBONDED_CLASH_THRESHOLD_ANGSTROM = 0.70


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

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


def xyz(row: dict[str, str]) -> tuple[float, float, float]:
    return (
        float(row["x_angstrom"]),
        float(row["y_angstrom"]),
        float(row["z_angstrom"]),
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
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def main() -> None:
    all_summary_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []

    all_fragments_pass = True

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
            row["atom_id"]: row
            for row in atoms
        }

        edges: set[tuple[str, str]] = set()
        edge_source: dict[tuple[str, str], str] = {}

        # Original internal graph edges.
        for row in internal_edges:
            edge = canonical_edge(
                row["source_node"],
                row["target_node"],
            )
            edges.add(edge)
            edge_source[edge] = "ORIGINAL_INTERNAL_EDGE"

        # Existing H atoms restored from the original R2 graph.
        for row in boundary_audit:
            if row["preliminary_action"] != "INCLUDE_EXISTING_HYDROGEN":
                continue

            edge = canonical_edge(
                row["inside_node"],
                row["outside_node"],
            )
            edges.add(edge)
            edge_source[edge] = "RESTORED_EXISTING_H_EDGE"

        # Artificial capping bonds.
        for row in caps:
            edge = canonical_edge(
                row["parent_inside_node"],
                row["cap_id"],
            )
            edges.add(edge)
            edge_source[edge] = "ARTIFICIAL_CAP_BOND"

        unknown_edge_atoms = sorted(
            {
                atom_id
                for edge in edges
                for atom_id in edge
                if atom_id not in atom_lookup
            }
        )

        if unknown_edge_atoms:
            raise RuntimeError(
                f"{label}: connectivity references unknown atoms: "
                f"{unknown_edge_atoms}"
            )

        adjacency: dict[str, set[str]] = defaultdict(set)

        for first, second in edges:
            adjacency[first].add(second)
            adjacency[second].add(first)

        valence_rows: list[dict[str, Any]] = []
        valence_failures = 0

        for atom_id in sorted(atom_lookup):
            atom = atom_lookup[atom_id]
            element = atom["element"]
            observed_degree = len(adjacency[atom_id])
            expected_degree = EXPECTED_DEGREE[element]
            passed = observed_degree == expected_degree

            if not passed:
                valence_failures += 1

            valence_rows.append(
                {
                    "fragment": label,
                    "atom_id": atom_id,
                    "element": element,
                    "atom_role": atom["atom_role"],
                    "node_type": atom["node_type"],
                    "observed_graph_degree": observed_degree,
                    "expected_graph_degree": expected_degree,
                    "degree_gate_pass": passed,
                    "neighbors": "|".join(
                        sorted(adjacency[atom_id])
                    ),
                }
            )

        bonded_rows: list[dict[str, Any]] = []
        bonded_failures = 0

        for first, second in sorted(edges):
            atom_a = atom_lookup[first]
            atom_b = atom_lookup[second]

            pair = frozenset(
                (atom_a["element"], atom_b["element"])
            )
            measured = distance(
                xyz(atom_a),
                xyz(atom_b),
            )

            limits = BONDED_DISTANCE_RANGES.get(pair)

            if limits is None:
                passed = False
                minimum = ""
                maximum = ""
            else:
                minimum, maximum = limits
                passed = minimum <= measured <= maximum

            if not passed:
                bonded_failures += 1

            bonded_rows.append(
                {
                    "fragment": label,
                    "atom_1": first,
                    "element_1": atom_a["element"],
                    "atom_2": second,
                    "element_2": atom_b["element"],
                    "edge_source": edge_source[(first, second)],
                    "distance_angstrom": f"{measured:.10f}",
                    "allowed_min_angstrom": minimum,
                    "allowed_max_angstrom": maximum,
                    "bond_distance_gate_pass": passed,
                }
            )

        atom_ids = sorted(atom_lookup)
        nonbonded_rows: list[dict[str, Any]] = []

        minimum_nonbonded_distance = float("inf")
        minimum_nonbonded_pair = ("", "")
        nonbonded_clashes = 0

        for i, first in enumerate(atom_ids):
            for second in atom_ids[i + 1:]:
                edge = canonical_edge(first, second)

                if edge in edges:
                    continue

                measured = distance(
                    xyz(atom_lookup[first]),
                    xyz(atom_lookup[second]),
                )

                if measured < minimum_nonbonded_distance:
                    minimum_nonbonded_distance = measured
                    minimum_nonbonded_pair = (first, second)

                clash = (
                    measured
                    < NONBONDED_CLASH_THRESHOLD_ANGSTROM
                )

                if clash:
                    nonbonded_clashes += 1

                # Preserve only potentially relevant close contacts.
                if measured < 1.50:
                    nonbonded_rows.append(
                        {
                            "fragment": label,
                            "atom_1": first,
                            "element_1": atom_lookup[first]["element"],
                            "atom_2": second,
                            "element_2": atom_lookup[second]["element"],
                            "distance_angstrom": f"{measured:.10f}",
                            "clash_threshold_angstrom": (
                                NONBONDED_CLASH_THRESHOLD_ANGSTROM
                            ),
                            "steric_clash": clash,
                        }
                    )

        element_counts = Counter(
            row["element"]
            for row in atoms
        )

        formula = "".join(
            f"{element}{element_counts[element]}"
            for element in ("B", "N", "H")
            if element_counts[element]
        )

        neutral_valence_electrons = sum(
            VALENCE_ELECTRONS[row["element"]]
            for row in atoms
        )

        charge_rows: list[dict[str, Any]] = []

        for charge in range(-2, 3):
            electron_count = (
                neutral_valence_electrons - charge
            )
            parity = (
                "EVEN"
                if electron_count % 2 == 0
                else "ODD"
            )

            minimum_multiplicity = (
                1
                if parity == "EVEN"
                else 2
            )

            charge_rows.append(
                {
                    "fragment": label,
                    "candidate_charge": charge,
                    "valence_electron_count": electron_count,
                    "electron_parity": parity,
                    "minimum_spin_multiplicity": (
                        minimum_multiplicity
                    ),
                    "recommended_initial_state": (
                        charge == 0
                        and minimum_multiplicity == 1
                    ),
                }
            )

        provisional_charge = 0
        provisional_multiplicity = 1

        atom_count_pass = len(atoms) == 22
        formula_pass = (
            element_counts["B"] == 5
            and element_counts["N"] == 5
            and element_counts["H"] == 12
        )
        electron_parity_pass = (
            neutral_valence_electrons % 2 == 0
        )

        fragment_pass = all(
            (
                atom_count_pass,
                formula_pass,
                valence_failures == 0,
                bonded_failures == 0,
                nonbonded_clashes == 0,
                electron_parity_pass,
            )
        )

        all_fragments_pass = (
            all_fragments_pass and fragment_pass
        )

        write_csv(
            F06_DIR / f"{label}_valence_audit.csv",
            valence_rows,
        )
        write_csv(
            F06_DIR / f"{label}_bond_distance_audit.csv",
            bonded_rows,
        )

        # Always write a nonbonded output, even if no close contacts.
        if nonbonded_rows:
            write_csv(
                F06_DIR
                / f"{label}_close_nonbonded_contacts.csv",
                nonbonded_rows,
            )
        else:
            (
                F06_DIR
                / f"{label}_close_nonbonded_contacts.csv"
            ).write_text(
                (
                    "fragment,atom_1,element_1,atom_2,"
                    "element_2,distance_angstrom,"
                    "clash_threshold_angstrom,steric_clash\n"
                ),
                encoding="utf-8",
            )

        write_csv(
            F06_DIR
            / f"{label}_charge_multiplicity_candidates.csv",
            charge_rows,
        )

        summary = {
            "fragment": label,
            "atom_count": len(atoms),
            "formula": formula,
            "B_atoms": element_counts["B"],
            "N_atoms": element_counts["N"],
            "H_atoms": element_counts["H"],
            "connectivity_edges": len(edges),
            "valence_failures": valence_failures,
            "bond_distance_failures": bonded_failures,
            "minimum_nonbonded_distance_angstrom": (
                f"{minimum_nonbonded_distance:.10f}"
            ),
            "minimum_nonbonded_pair": (
                "|".join(minimum_nonbonded_pair)
            ),
            "nonbonded_clashes": nonbonded_clashes,
            "neutral_valence_electrons": (
                neutral_valence_electrons
            ),
            "provisional_charge": provisional_charge,
            "provisional_multiplicity": (
                provisional_multiplicity
            ),
            "pre_qm_gate_pass": fragment_pass,
        }

        all_summary_rows.append(summary)

        report_sections.extend(
            [
                f"## {label}",
                "",
                f"- Formula: **{formula}**",
                f"- Atoms: **{len(atoms)}**",
                f"- Reconstructed bonds: **{len(edges)}**",
                f"- Degree/valence failures: **{valence_failures}**",
                (
                    f"- Bond-distance failures: "
                    f"**{bonded_failures}**"
                ),
                (
                    f"- Minimum nonbonded distance: "
                    f"**{minimum_nonbonded_distance:.6f} Å** "
                    f"(`{' — '.join(minimum_nonbonded_pair)}`)"
                ),
                f"- Steric clashes: **{nonbonded_clashes}**",
                (
                    f"- Neutral valence-electron count: "
                    f"**{neutral_valence_electrons}**"
                ),
                (
                    "- Provisional electronic state: "
                    f"**charge {provisional_charge}, "
                    f"multiplicity {provisional_multiplicity}**"
                ),
                (
                    f"- Pre-QM gate: "
                    f"**{'PASS' if fragment_pass else 'FAIL'}**"
                ),
                "",
            ]
        )

    summary_path = F06_DIR / (
        "QM_F06_capped_fragment_pre_qm_summary.csv"
    )
    write_csv(summary_path, all_summary_rows)

    decision = (
        "QM_F06_CAPPED_FRAGMENTS_PASS_PRE_QM_GATE"
        if all_fragments_pass
        else "QM_F06_CAPPED_FRAGMENTS_FAIL_PRE_QM_GATE"
    )

    report_path = F06_DIR / (
        "QM_F06_CAPPED_FRAGMENT_PRE_QM_AUDIT.md"
    )
    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 Capped-Fragment Pre-QM Audit — Day026",
                "",
                "## Scope",
                "",
                (
                    "The chemically capped LOWER and UPPER bridge "
                    "fragments were audited for atom counts, graph "
                    "valence, bonded distances, nonbonded clashes and "
                    "electronic-state parity."
                ),
                "",
                f"## Decision: **{decision}**",
                "",
                *report_sections,
                "## Charge and multiplicity interpretation",
                "",
                (
                    "Both fragments have formula B5N5H12 and 52 neutral "
                    "valence electrons. A neutral closed-shell singlet "
                    "(charge 0, multiplicity 1) is therefore the provisional "
                    "initial electronic state."
                ),
                "",
                (
                    "This is an electron-count and stoichiometric assignment, "
                    "not yet an electronic-structure validation. The eventual "
                    "QM workflow must confirm SCF stability and absence of a "
                    "lower-energy open-shell solution."
                ),
                "",
                "## Authorization state",
                "",
                "- Geometry construction: **COMPLETED**",
                "- Pre-QM audit: **COMPLETED**",
                "- QM input preparation: **PENDING THIS GATE**",
                "- QM calculation executed: **NO**",
                "",
            ]
        ),
        encoding="utf-8",
    )

    json_path = F06_DIR / (
        "QM_F06_capped_fragment_pre_qm_summary.json"
    )
    json_path.write_text(
        json.dumps(
            {
                "decision": decision,
                "all_fragments_pass": all_fragments_pass,
                "provisional_charge": 0,
                "provisional_multiplicity": 1,
                "qm_calculation_executed": False,
                "required_next_step": (
                    "PREPARE_QM_INPUTS"
                    if all_fragments_pass
                    else "CORRECT_FAILED_PRE_QM_GATES"
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 capped-fragment pre-QM audit completed.")
    print(f"Decision: {decision}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
