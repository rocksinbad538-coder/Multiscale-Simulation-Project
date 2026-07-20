#!/usr/bin/env python3
"""
Prepare QM_F06 UPPER V3-A2 after identifying a fixed geminal cap artifact.

Correction:
- retain Boundary V3 topology;
- use the final available V3-A restart4 geometry as starting structure;
- keep P:1641 / local index 11 fixed;
- release HCAP:UPPER:01 / local index 14;
- release HCAP:UPPER:04 / local index 17;
- reconstruct both N-H vectors as a symmetric trigonal-planar pair;
- preserve all other V3-A constraints;
- use a fresh SCF guess;
- do not execute ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V3_ROOT = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3"
)

WORKFLOW_ROOT = V3_ROOT / "orca_v3_workflow"
MAP_PATH = WORKFLOW_ROOT / "v3a_atom_role_constraint_map.csv"

RESTART4_EXECUTION = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_upper_boundary_v3a_restart4_executions/"
    "restart4_20260720T111123"
)

SOURCE_XYZ = RESTART4_EXECUTION / "restart4.xyz"
SOURCE_INPUT = WORKFLOW_ROOT / "restart4/v3a_restart4.inp"

OUTPUT_ROOT = V3_ROOT / "orca_v3a2_workflow"
OUTPUT_XYZ = OUTPUT_ROOT / "v3a2_start.xyz"
OUTPUT_INPUT = OUTPUT_ROOT / "v3a2.inp"
OUTPUT_MAP = OUTPUT_ROOT / "v3a2_atom_role_constraint_map.csv"
OUTPUT_SUMMARY = OUTPUT_ROOT / "v3a2_preparation_summary.json"

CENTER_ID = "P:1641"
CAP_IDS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

NH_DISTANCE_A = 1.010
TARGET_HNH_ANGLE_DEG = 120.0


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(vector, scalar):
    return tuple(value * scalar for value in vector)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(vector):
    return math.sqrt(dot(vector, vector))


def normalize(vector):
    length = norm(vector)

    if length < 1.0e-12:
        raise RuntimeError("Cannot normalize zero-length vector.")

    return scale(vector, 1.0 / length)


def distance(a, b):
    return norm(subtract(a, b))


def angle_deg(a, vertex, c):
    u = normalize(subtract(a, vertex))
    v = normalize(subtract(c, vertex))

    cosine = max(-1.0, min(1.0, dot(u, v)))
    return math.degrees(math.acos(cosine))


def read_xyz(path: Path):
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(lines[0].strip())
    atoms = []

    for index, line in enumerate(lines[2:2 + atom_count]):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ atom line {index}: {line}"
            )

        atoms.append(
            {
                "index": index,
                "element": fields[0],
                "xyz": (
                    float(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                ),
            }
        )

    if len(atoms) != atom_count:
        raise RuntimeError("Incomplete XYZ geometry.")

    return atoms


def write_xyz(path: Path, comment: str, atoms) -> None:
    lines = [
        str(len(atoms)),
        comment,
    ]

    for atom in atoms:
        x, y, z = atom["xyz"]

        lines.append(
            f"{atom['element']:2s} "
            f"{x: .10f} "
            f"{y: .10f} "
            f"{z: .10f}"
        )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def read_map(path: Path):
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if "index_0based" not in fieldnames:
        raise RuntimeError(
            "Expected index_0based column is absent."
        )

    return fieldnames, rows


def nearest_heavy_neighbor(
    atoms,
    center_index: int,
    excluded_indices: set[int],
) -> int:
    """
    Identify the unique chemically compatible h-BN neighbor.

    Geometric proximity alone is insufficient because the compact
    fragment contains close nonbonded heavy atoms. For an N center,
    only B is a chemically compatible heavy neighbor; for a B center,
    only N is compatible.
    """
    center = atoms[center_index]
    center_element = center["element"]

    expected_neighbor_element = {
        "N": "B",
        "B": "N",
    }.get(center_element)

    if expected_neighbor_element is None:
        raise RuntimeError(
            "Unsupported center element for h-BN neighbor search: "
            f"{center_element}"
        )

    candidates = []

    for atom in atoms:
        if atom["index"] == center_index:
            continue

        if atom["index"] in excluded_indices:
            continue

        if atom["element"] != expected_neighbor_element:
            continue

        value = distance(
            center["xyz"],
            atom["xyz"],
        )

        if 1.20 <= value <= 1.85:
            candidates.append(
                (value, atom["index"])
            )

    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one chemically compatible h-BN "
            f"neighbor for center {center_index} "
            f"({center_element}); found {candidates}"
        )

    return candidates[0][1]


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    require_file(SOURCE_XYZ)
    require_file(SOURCE_INPUT)
    require_file(MAP_PATH)

    atoms = read_xyz(SOURCE_XYZ)
    fieldnames, map_rows = read_map(MAP_PATH)

    by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in map_rows
    }

    center_index = by_id[CENTER_ID]
    cap_indices = sorted(
        by_id[cap_id]
        for cap_id in CAP_IDS
    )

    if center_index != 11:
        raise RuntimeError(
            f"Unexpected center index: {center_index}"
        )

    if cap_indices != [14, 17]:
        raise RuntimeError(
            f"Unexpected cap indices: {cap_indices}"
        )

    if atoms[center_index]["element"] != "N":
        raise RuntimeError("Center atom is not nitrogen.")

    if any(
        atoms[index]["element"] != "H"
        for index in cap_indices
    ):
        raise RuntimeError("Target caps are not hydrogens.")

    heavy_neighbor_index = nearest_heavy_neighbor(
        atoms,
        center_index,
        set(cap_indices),
    )

    center_xyz = atoms[center_index]["xyz"]
    heavy_xyz = atoms[heavy_neighbor_index]["xyz"]

    # Direction opposite to the unique N-B bond.
    outward_axis = normalize(
        subtract(center_xyz, heavy_xyz)
    )

    old_h1_vector = normalize(
        subtract(
            atoms[cap_indices[0]]["xyz"],
            center_xyz,
        )
    )

    old_h2_vector = normalize(
        subtract(
            atoms[cap_indices[1]]["xyz"],
            center_xyz,
        )
    )

    # Use the original cap geometry only to define the local plane.
    plane_normal = cross(
        old_h1_vector,
        old_h2_vector,
    )

    if norm(plane_normal) < 1.0e-8:
        # Fallback: choose a stable perpendicular direction.
        reference = (
            (1.0, 0.0, 0.0)
            if abs(outward_axis[0]) < 0.8
            else (0.0, 1.0, 0.0)
        )
        plane_normal = cross(
            outward_axis,
            reference,
        )

    plane_normal = normalize(plane_normal)

    in_plane_perpendicular = normalize(
        cross(
            plane_normal,
            outward_axis,
        )
    )

    half_angle_rad = math.radians(
        TARGET_HNH_ANGLE_DEG / 2.0
    )

    direction_1 = normalize(
        add(
            scale(
                outward_axis,
                math.cos(half_angle_rad),
            ),
            scale(
                in_plane_perpendicular,
                math.sin(half_angle_rad),
            ),
        )
    )

    direction_2 = normalize(
        add(
            scale(
                outward_axis,
                math.cos(half_angle_rad),
            ),
            scale(
                in_plane_perpendicular,
                -math.sin(half_angle_rad),
            ),
        )
    )

    atoms[cap_indices[0]]["xyz"] = add(
        center_xyz,
        scale(direction_1, NH_DISTANCE_A),
    )

    atoms[cap_indices[1]]["xyz"] = add(
        center_xyz,
        scale(direction_2, NH_DISTANCE_A),
    )

    h1_xyz = atoms[cap_indices[0]]["xyz"]
    h2_xyz = atoms[cap_indices[1]]["xyz"]

    repaired_metrics = {
        "N_H1_A": distance(center_xyz, h1_xyz),
        "N_H2_A": distance(center_xyz, h2_xyz),
        "H1_H2_A": distance(h1_xyz, h2_xyz),
        "H_N_H_deg": angle_deg(
            h1_xyz,
            center_xyz,
            h2_xyz,
        ),
    }

    write_xyz(
        OUTPUT_XYZ,
        (
            "QM_F06_UPPER_V3A2_START; "
            "local geminal cap repair; "
            "P:1641 fixed; HCAP:UPPER:01 and "
            "HCAP:UPPER:04 mobile; no QM executed"
        ),
        atoms,
    )

    # Preserve all previous fixed atoms except the two repaired caps.
    source_input_text = SOURCE_INPUT.read_text(
        encoding="utf-8"
    )

    previous_fixed = {
        int(value)
        for value in re.findall(
            r"\{\s*C\s+(\d+)\s+C\s*\}",
            source_input_text,
        )
    }

    new_fixed = sorted(
        previous_fixed - set(cap_indices)
    )

    if center_index not in new_fixed:
        raise RuntimeError(
            "The real nitrogen center must remain fixed in V3-A2."
        )

    if any(index in new_fixed for index in cap_indices):
        raise RuntimeError(
            "Target caps were not released."
        )

    xyz_block = "\n".join(
        (
            f"{atom['element']:2s} "
            f"{atom['xyz'][0]: .10f} "
            f"{atom['xyz'][1]: .10f} "
            f"{atom['xyz'][2]: .10f}"
        )
        for atom in atoms
    )

    input_text = re.sub(
        r"(?ms)^\s*\*\s+xyz\s+-?\d+\s+\d+\s*$"
        r".*?"
        r"^\s*\*\s*$",
        "* xyz 0 1\n"
        + xyz_block
        + "\n*",
        source_input_text,
        count=1,
    )

    constraints_block = "\n".join(
        f"    {{ C {index} C }}"
        for index in new_fixed
    )

    input_text, substitutions = re.subn(
        r"(?ms)(%geom.*?Constraints\s*\n)"
        r".*?"
        r"(\s*end\s*\n\s*end)",
        (
            r"\1"
            + constraints_block
            + "\n"
            + r"\2"
        ),
        input_text,
        count=1,
    )

    if substitutions != 1:
        raise RuntimeError(
            "Could not replace the ORCA constraint block."
        )

    input_text = re.sub(
        r"(?im)^\s*nprocs\s+\d+\s*$",
        "  nprocs 4",
        input_text,
        count=1,
    )

    input_text = re.sub(
        r"(?im)^%maxcore\s+\d+\s*$",
        "%maxcore 2500",
        input_text,
        count=1,
    )

    input_text = input_text.replace(
        "# QM_F06 UPPER Boundary V3-A clean restart",
        "# QM_F06 UPPER Boundary V3-A2 local geminal-cap repair",
    )

    OUTPUT_INPUT.write_text(
        input_text,
        encoding="utf-8",
    )

    updated_rows = []

    for row in map_rows:
        updated = dict(row)
        index = int(row["index_0based"])

        updated["v3a2_fixed"] = str(
            index in new_fixed
        )
        updated["v3a2_mobile"] = str(
            index not in new_fixed
        )

        if index in cap_indices:
            updated["v3a2_mobility_basis"] = (
                "REPAIRED_GEMINAL_ARTIFICIAL_CAP"
            )
        elif index == center_index:
            updated["v3a2_mobility_basis"] = (
                "FIXED_REAL_CENTER_LOCAL_CAP_REPAIR"
            )
        else:
            updated["v3a2_mobility_basis"] = (
                row["mobility_basis"]
            )

        updated_rows.append(updated)

    output_fields = [
        *fieldnames,
        "v3a2_fixed",
        "v3a2_mobile",
        "v3a2_mobility_basis",
    ]

    with OUTPUT_MAP.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_fields,
        )
        writer.writeheader()
        writer.writerows(updated_rows)

    summary = {
        "decision": (
            "QM_F06_UPPER_V3A2_LOCAL_CAP_REPAIR_PREPARED_"
            "PRE_QM_AUDIT_REQUIRED"
        ),
        "source_geometry": str(
            SOURCE_XYZ.relative_to(ROOT)
        ),
        "source_input": str(
            SOURCE_INPUT.relative_to(ROOT)
        ),
        "center": {
            "index_0based": center_index,
            "atom_id": CENTER_ID,
            "element": "N",
            "fixed": True,
        },
        "heavy_neighbor": {
            "index_0based": heavy_neighbor_index,
            "element": atoms[
                heavy_neighbor_index
            ]["element"],
        },
        "repaired_caps": [
            {
                "index_0based": index,
                "atom_id": next(
                    atom_id
                    for atom_id, mapped_index
                    in by_id.items()
                    if mapped_index == index
                ),
                "fixed_before": True,
                "fixed_v3a2": False,
            }
            for index in cap_indices
        ],
        "previous_fixed_count": len(previous_fixed),
        "v3a2_fixed_count": len(new_fixed),
        "v3a2_mobile_count": (
            len(atoms) - len(new_fixed)
        ),
        "target_NH_distance_A": NH_DISTANCE_A,
        "target_HNH_angle_deg": (
            TARGET_HNH_ANGLE_DEG
        ),
        "repaired_metrics": repaired_metrics,
        "fresh_scf_guess": True,
        "gbw_reused": False,
        "nprocs": 4,
        "maxcore_mb_per_process": 2500,
        "files": {
            "xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "input": str(
                OUTPUT_INPUT.relative_to(ROOT)
            ),
            "constraint_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
        },
        "sha256": {
            "xyz": sha256(OUTPUT_XYZ),
            "input": sha256(OUTPUT_INPUT),
            "constraint_map": sha256(
                OUTPUT_MAP
            ),
        },
        "authorization": {
            "pre_qm_audit_authorized": True,
            "orca_execution_authorized": False,
            "geometry_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_SUMMARY.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V3-A2 PREPARATION")
    print("=" * 78)
    print("Center:", CENTER_ID, center_index)
    print(
        "Heavy neighbor index:",
        heavy_neighbor_index,
    )
    print("Released caps:", cap_indices)
    print("Previous fixed atoms:", len(previous_fixed))
    print("V3-A2 fixed atoms:", len(new_fixed))
    print(
        "V3-A2 mobile atoms:",
        len(atoms) - len(new_fixed),
    )
    print(
        "N-H distances:",
        repaired_metrics["N_H1_A"],
        repaired_metrics["N_H2_A"],
    )
    print(
        "H-H distance:",
        repaired_metrics["H1_H2_A"],
    )
    print(
        "H-N-H angle:",
        repaired_metrics["H_N_H_deg"],
    )
    print("XYZ:", OUTPUT_XYZ)
    print("Input:", OUTPUT_INPUT)
    print("Summary:", OUTPUT_SUMMARY)
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
