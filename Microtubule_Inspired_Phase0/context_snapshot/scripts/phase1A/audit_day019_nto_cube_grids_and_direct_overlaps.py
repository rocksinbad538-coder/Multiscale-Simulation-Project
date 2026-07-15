#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CUBE_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/day019_nto_cube_generation"
)

SELECTION_MANIFEST = (
    CUBE_ROOT / "CUBE_SELECTION_MANIFEST_DAY019.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_nto_cube_overlap_analysis"
)

CUBE_AUDIT_CSV = OUTPUT_ROOT / "cube_integrity_and_norm_audit.csv"
DIRECT_OVERLAPS_CSV = OUTPUT_ROOT / "vacuum_embedding_direct_overlaps.csv"
PAIR_SUMMARY_CSV = OUTPUT_ROOT / "vacuum_embedding_pair_similarity.csv"
REPORT_MD = OUTPUT_ROOT / "DAY019_CUBE_GRID_AND_DIRECT_OVERLAP_AUDIT.md"

GRID_TOL = 1.0e-10
ATOM_TOL = 1.0e-8
MIN_NORM = 1.0e-12


@dataclass(frozen=True)
class CubeMetadata:
    natoms: int
    natoms_raw: int
    origin: np.ndarray
    shape: tuple[int, int, int]
    vectors: np.ndarray
    atomic_numbers: np.ndarray
    atom_charges: np.ndarray
    atom_coordinates: np.ndarray
    n_datasets: int
    dataset_ids: tuple[int, ...]
    voxel_volume: float
    n_values: int


@dataclass
class CubeData:
    metadata: CubeMetadata
    values: np.ndarray


def log(message: str = "") -> None:
    print(message, flush=True)


def read_selection_manifest() -> list[dict[str, str]]:
    if not SELECTION_MANIFEST.is_file():
        raise SystemExit(
            f"Missing cube-selection manifest: {SELECTION_MANIFEST}"
        )

    with SELECTION_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 48:
        raise SystemExit(
            f"Expected 48 cube selections, found {len(rows)}."
        )

    cube_paths = [row["cube_file"] for row in rows]

    if len(set(cube_paths)) != 48:
        raise SystemExit(
            "Cube-selection manifest does not contain 48 unique cube files."
        )

    return rows


def parse_int_like(value: str) -> int:
    return int(float(value))


def read_cube(path: Path, load_values: bool = True) -> CubeData:
    """Read a Gaussian/ORCA cube file, including MO dataset-ID records.

    For molecular-orbital cubes, a negative NATOMS value indicates that one
    or more dataset identifiers are stored immediately after the atom records.
    Those identifiers are header metadata and must not be counted as volumetric
    values. The selected ORCA NTO cubes are expected to contain one dataset.
    """
    if not path.is_file():
        raise RuntimeError(f"Missing cube file: {path}")

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        comment_1 = handle.readline()
        comment_2 = handle.readline()

        if not comment_1 or not comment_2:
            raise RuntimeError(f"Truncated cube header: {path}")

        line = handle.readline().split()
        if len(line) < 4:
            raise RuntimeError(f"Invalid atom/origin line: {path}")

        natoms_raw = parse_int_like(line[0])
        natoms = abs(natoms_raw)
        origin = np.array(
            [float(line[1]), float(line[2]), float(line[3])],
            dtype=np.float64,
        )

        # Some cube variants include NVAL as a fifth field when NATOMS is
        # positive. ORCA MO cubes use negative NATOMS plus a DSET_IDS record.
        nval_header = (
            parse_int_like(line[4])
            if len(line) >= 5
            else None
        )

        shape: list[int] = []
        vectors: list[list[float]] = []

        for axis in range(3):
            line = handle.readline().split()
            if len(line) < 4:
                raise RuntimeError(
                    f"Invalid grid line {axis + 1}: {path}"
                )

            count = abs(parse_int_like(line[0]))
            vector = [
                float(line[1]),
                float(line[2]),
                float(line[3]),
            ]

            shape.append(count)
            vectors.append(vector)

        atomic_numbers: list[int] = []
        atom_charges: list[float] = []
        atom_coordinates: list[list[float]] = []

        for atom_index in range(natoms):
            line = handle.readline().split()
            if len(line) < 5:
                raise RuntimeError(
                    f"Invalid atom line {atom_index + 1}: {path}"
                )

            atomic_numbers.append(parse_int_like(line[0]))
            atom_charges.append(float(line[1]))
            atom_coordinates.append(
                [float(line[2]), float(line[3]), float(line[4])]
            )

        dataset_ids: list[int] = []

        if natoms_raw < 0:
            dataset_line = handle.readline()

            if not dataset_line:
                raise RuntimeError(
                    f"Negative NATOMS but missing dataset-ID record: {path}"
                )

            dataset_tokens = dataset_line.split()

            if not dataset_tokens:
                raise RuntimeError(
                    f"Empty dataset-ID record: {path}"
                )

            n_datasets = parse_int_like(dataset_tokens[0])

            if n_datasets <= 0:
                raise RuntimeError(
                    f"Invalid dataset count {n_datasets}: {path}"
                )

            dataset_ids.extend(
                parse_int_like(token)
                for token in dataset_tokens[1:]
            )

            while len(dataset_ids) < n_datasets:
                continuation = handle.readline()

                if not continuation:
                    raise RuntimeError(
                        f"Truncated dataset-ID record: {path}"
                    )

                dataset_ids.extend(
                    parse_int_like(token)
                    for token in continuation.split()
                )

            if len(dataset_ids) != n_datasets:
                raise RuntimeError(
                    f"Dataset-ID count mismatch for {path}: "
                    f"expected {n_datasets}, found {len(dataset_ids)}"
                )
        else:
            n_datasets = nval_header if nval_header is not None else 1

            if n_datasets <= 0:
                raise RuntimeError(
                    f"Invalid NVAL/dataset count {n_datasets}: {path}"
                )

        # Each selected file was generated for one individual NTO orbital.
        # Refuse multi-dataset files rather than silently mixing orbitals.
        if n_datasets != 1:
            raise RuntimeError(
                f"Expected one orbital dataset in {path}, "
                f"found {n_datasets} with IDs {dataset_ids}"
            )

        expected_values = int(np.prod(shape)) * n_datasets

        if load_values:
            values = np.fromiter(
                (
                    float(token)
                    for line in handle
                    for token in line.split()
                ),
                dtype=np.float64,
            )

            if values.size != expected_values:
                raise RuntimeError(
                    f"Cube data-count mismatch for {path}: "
                    f"expected {expected_values}, found {values.size}; "
                    f"NATOMS={natoms_raw}, n_datasets={n_datasets}, "
                    f"dataset_ids={dataset_ids}"
                )
        else:
            values = np.empty(0, dtype=np.float64)

    vectors_array = np.array(vectors, dtype=np.float64)
    voxel_volume = abs(float(np.linalg.det(vectors_array)))

    if not math.isfinite(voxel_volume) or voxel_volume <= 0.0:
        raise RuntimeError(
            f"Invalid voxel volume {voxel_volume} for {path}"
        )

    metadata = CubeMetadata(
        natoms=natoms,
        natoms_raw=natoms_raw,
        origin=origin,
        shape=(shape[0], shape[1], shape[2]),
        vectors=vectors_array,
        atomic_numbers=np.array(atomic_numbers, dtype=np.int64),
        atom_charges=np.array(atom_charges, dtype=np.float64),
        atom_coordinates=np.array(
            atom_coordinates,
            dtype=np.float64,
        ),
        n_datasets=n_datasets,
        dataset_ids=tuple(dataset_ids),
        voxel_volume=voxel_volume,
        n_values=expected_values,
    )

    return CubeData(
        metadata=metadata,
        values=values,
    )


def cube_norm(cube: CubeData) -> float:
    integral = float(
        np.dot(cube.values, cube.values)
        * cube.metadata.voxel_volume
    )

    if not math.isfinite(integral) or integral <= MIN_NORM:
        raise RuntimeError(
            f"Invalid orbital norm-squared integral: {integral}"
        )

    return math.sqrt(integral)


def max_abs_difference(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    if first.shape != second.shape:
        return float("inf")

    if first.size == 0:
        return 0.0

    return float(np.max(np.abs(first - second)))


def compare_metadata(
    reference: CubeMetadata,
    moving: CubeMetadata,
) -> dict[str, object]:
    same_shape = reference.shape == moving.shape
    same_natoms = reference.natoms == moving.natoms

    origin_max_abs_diff = max_abs_difference(
        reference.origin,
        moving.origin,
    )

    vectors_max_abs_diff = max_abs_difference(
        reference.vectors,
        moving.vectors,
    )

    atomic_numbers_equal = bool(
        np.array_equal(
            reference.atomic_numbers,
            moving.atomic_numbers,
        )
    )

    atom_charges_max_abs_diff = max_abs_difference(
        reference.atom_charges,
        moving.atom_charges,
    )

    atom_coordinates_max_abs_diff = max_abs_difference(
        reference.atom_coordinates,
        moving.atom_coordinates,
    )

    same_grid = (
        same_shape
        and origin_max_abs_diff <= GRID_TOL
        and vectors_max_abs_diff <= GRID_TOL
    )

    same_geometry = (
        same_natoms
        and atomic_numbers_equal
        and atom_charges_max_abs_diff <= ATOM_TOL
        and atom_coordinates_max_abs_diff <= ATOM_TOL
    )

    return {
        "same_shape": same_shape,
        "same_natoms": same_natoms,
        "origin_max_abs_diff": origin_max_abs_diff,
        "vectors_max_abs_diff": vectors_max_abs_diff,
        "atomic_numbers_equal": atomic_numbers_equal,
        "atom_charges_max_abs_diff": atom_charges_max_abs_diff,
        "atom_coordinates_max_abs_diff": atom_coordinates_max_abs_diff,
        "same_grid": same_grid,
        "same_geometry": same_geometry,
    }


def normalized_overlap(
    reference: CubeData,
    moving: CubeData,
) -> dict[str, float]:
    metadata_comparison = compare_metadata(
        reference.metadata,
        moving.metadata,
    )

    if not bool(metadata_comparison["same_grid"]):
        raise RuntimeError(
            "Direct overlap requested for nonidentical cube grids."
        )

    if reference.values.shape != moving.values.shape:
        raise RuntimeError(
            "Direct overlap requested for arrays with different sizes."
        )

    voxel_volume = reference.metadata.voxel_volume

    norm_reference = cube_norm(reference)
    norm_moving = cube_norm(moving)

    signed = float(
        np.dot(reference.values, moving.values)
        * voxel_volume
        / (norm_reference * norm_moving)
    )

    signed = max(-1.0, min(1.0, signed))
    absolute = abs(signed)

    optimal_sign = 1.0 if signed >= 0.0 else -1.0

    reference_normalized = reference.values / norm_reference
    moving_normalized = (
        optimal_sign * moving.values / norm_moving
    )

    rms_l2_difference = math.sqrt(
        float(
            np.dot(
                reference_normalized - moving_normalized,
                reference_normalized - moving_normalized,
            )
            * voxel_volume
        )
    )

    return {
        "reference_norm": norm_reference,
        "moving_norm": norm_moving,
        "signed_overlap": signed,
        "absolute_overlap": absolute,
        "optimal_global_sign": optimal_sign,
        "l2_difference_after_sign_alignment": rms_l2_difference,
    }


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def bool_from_csv(value: str) -> bool:
    return value.strip().lower() == "true"


def main() -> None:
    selection = read_selection_manifest()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    log("Day019 cube-grid and direct-overlap audit")
    log(f"Cube selections: {len(selection)}")

    cube_audit_rows: list[dict[str, object]] = []
    metadata_cache: dict[str, CubeMetadata] = {}

    for index, row in enumerate(selection, start=1):
        cube_path = (
            PROJECT_ROOT / row["cube_file"]
        ).resolve()

        cube = read_cube(cube_path, load_values=True)
        norm = cube_norm(cube)

        metadata_cache[row["cube_file"]] = cube.metadata

        cube_audit_rows.append(
            {
                "calculation_type": row["calculation_type"],
                "frame": int(row["frame"]),
                "cluster": row["cluster"],
                "job": row["job"],
                "root": int(row["root"]),
                "is_tracked_root": bool_from_csv(
                    row["is_tracked_root"]
                ),
                "pair_rank": int(row["pair_rank"]),
                "pair_label": row["pair_label"],
                "orbital_role": row["orbital_role"],
                "orbital_index": int(row["orbital_index"]),
                "orbital_spin": row["orbital_spin"],
                "natoms": cube.metadata.natoms,
                "natoms_raw": cube.metadata.natoms_raw,
                "n_datasets": cube.metadata.n_datasets,
                "dataset_ids": ";".join(
                    str(value)
                    for value in cube.metadata.dataset_ids
                ),
                "nx": cube.metadata.shape[0],
                "ny": cube.metadata.shape[1],
                "nz": cube.metadata.shape[2],
                "n_values": cube.metadata.n_values,
                "voxel_volume_bohr3": cube.metadata.voxel_volume,
                "orbital_norm": norm,
                "cube_size_bytes": cube_path.stat().st_size,
                "cube_file": row["cube_file"],
            }
        )

        log(
            f"[{index:02d}/48] parsed "
            f"{cube_path.name} "
            f"norm={norm:.8f}"
        )

    cube_audit_fieldnames = [
        "calculation_type",
        "frame",
        "cluster",
        "job",
        "root",
        "is_tracked_root",
        "pair_rank",
        "pair_label",
        "orbital_role",
        "orbital_index",
        "orbital_spin",
        "natoms",
        "natoms_raw",
        "n_datasets",
        "dataset_ids",
        "nx",
        "ny",
        "nz",
        "n_values",
        "voxel_volume_bohr3",
        "orbital_norm",
        "cube_size_bytes",
        "cube_file",
    ]

    write_csv(
        CUBE_AUDIT_CSV,
        cube_audit_rows,
        cube_audit_fieldnames,
    )

    vacuum_rows: dict[
        tuple[str, int, int, str],
        dict[str, str],
    ] = {}

    embedded_rows: list[dict[str, str]] = []

    for row in selection:
        key = (
            row["cluster"],
            int(row["root"]),
            int(row["orbital_index"]),
            row["orbital_spin"],
        )

        if row["calculation_type"] == "vacuum_reference":
            vacuum_rows[key] = row
        elif row["calculation_type"] == "embedded":
            embedded_rows.append(row)

    direct_rows: list[dict[str, object]] = []

    for index, embedded_row in enumerate(
        embedded_rows,
        start=1,
    ):
        key = (
            embedded_row["cluster"],
            int(embedded_row["root"]),
            int(embedded_row["orbital_index"]),
            embedded_row["orbital_spin"],
        )

        if key not in vacuum_rows:
            raise RuntimeError(
                "Missing matching vacuum cube for embedded selection: "
                f"{key}"
            )

        vacuum_row = vacuum_rows[key]

        vacuum_path = (
            PROJECT_ROOT / vacuum_row["cube_file"]
        ).resolve()

        embedded_path = (
            PROJECT_ROOT / embedded_row["cube_file"]
        ).resolve()

        vacuum_cube = read_cube(vacuum_path, load_values=True)
        embedded_cube = read_cube(
            embedded_path,
            load_values=True,
        )

        metadata_comparison = compare_metadata(
            vacuum_cube.metadata,
            embedded_cube.metadata,
        )

        if not bool(metadata_comparison["same_grid"]):
            raise RuntimeError(
                "Frozen-geometry vacuum/embedded grids do not match: "
                f"{vacuum_path.name} vs {embedded_path.name}"
            )

        if not bool(metadata_comparison["same_geometry"]):
            raise RuntimeError(
                "Frozen-geometry vacuum/embedded atom records do not match: "
                f"{vacuum_path.name} vs {embedded_path.name}"
            )

        overlap = normalized_overlap(
            vacuum_cube,
            embedded_cube,
        )

        direct_rows.append(
            {
                "frame": int(embedded_row["frame"]),
                "cluster": embedded_row["cluster"],
                "root": int(embedded_row["root"]),
                "is_tracked_root": bool_from_csv(
                    embedded_row["is_tracked_root"]
                ),
                "pair_rank": int(embedded_row["pair_rank"]),
                "pair_label": embedded_row["pair_label"],
                "orbital_role": embedded_row["orbital_role"],
                "orbital_index": int(
                    embedded_row["orbital_index"]
                ),
                "orbital_spin": embedded_row["orbital_spin"],
                **metadata_comparison,
                **overlap,
                "vacuum_cube": vacuum_row["cube_file"],
                "embedded_cube": embedded_row["cube_file"],
            }
        )

        log(
            f"[direct {index:02d}/{len(embedded_rows)}] "
            f"frame{int(embedded_row['frame']):03d} "
            f"{embedded_row['cluster']} "
            f"S{embedded_row['root']} "
            f"MO{embedded_row['orbital_index']}"
            f"{embedded_row['orbital_spin']} "
            f"|S|={overlap['absolute_overlap']:.8f}"
        )

    if len(direct_rows) != 24:
        raise RuntimeError(
            f"Expected 24 direct vacuum/embedding overlaps, "
            f"found {len(direct_rows)}."
        )

    direct_fieldnames = [
        "frame",
        "cluster",
        "root",
        "is_tracked_root",
        "pair_rank",
        "pair_label",
        "orbital_role",
        "orbital_index",
        "orbital_spin",
        "same_shape",
        "same_natoms",
        "origin_max_abs_diff",
        "vectors_max_abs_diff",
        "atomic_numbers_equal",
        "atom_charges_max_abs_diff",
        "atom_coordinates_max_abs_diff",
        "same_grid",
        "same_geometry",
        "reference_norm",
        "moving_norm",
        "signed_overlap",
        "absolute_overlap",
        "optimal_global_sign",
        "l2_difference_after_sign_alignment",
        "vacuum_cube",
        "embedded_cube",
    ]

    write_csv(
        DIRECT_OVERLAPS_CSV,
        direct_rows,
        direct_fieldnames,
    )

    grouped: dict[
        tuple[int, str, int, bool, int, str],
        dict[str, dict[str, object]],
    ] = {}

    for row in direct_rows:
        key = (
            int(row["frame"]),
            str(row["cluster"]),
            int(row["root"]),
            bool(row["is_tracked_root"]),
            int(row["pair_rank"]),
            str(row["pair_label"]),
        )

        grouped.setdefault(key, {})[
            str(row["orbital_role"])
        ] = row

    pair_rows: list[dict[str, object]] = []

    for key, roles in sorted(grouped.items()):
        if set(roles) != {"hole", "particle"}:
            raise RuntimeError(
                f"Incomplete hole/particle overlap pair for {key}: "
                f"{sorted(roles)}"
            )

        hole_overlap = float(
            roles["hole"]["absolute_overlap"]
        )
        particle_overlap = float(
            roles["particle"]["absolute_overlap"]
        )

        pair_rows.append(
            {
                "frame": key[0],
                "cluster": key[1],
                "root": key[2],
                "is_tracked_root": key[3],
                "pair_rank": key[4],
                "pair_label": key[5],
                "hole_absolute_overlap": hole_overlap,
                "particle_absolute_overlap": particle_overlap,
                "pair_geometric_mean_overlap": math.sqrt(
                    hole_overlap * particle_overlap
                ),
                "pair_minimum_overlap": min(
                    hole_overlap,
                    particle_overlap,
                ),
            }
        )

    if len(pair_rows) != 12:
        raise RuntimeError(
            f"Expected 12 direct NTO-pair comparisons, "
            f"found {len(pair_rows)}."
        )

    pair_fieldnames = [
        "frame",
        "cluster",
        "root",
        "is_tracked_root",
        "pair_rank",
        "pair_label",
        "hole_absolute_overlap",
        "particle_absolute_overlap",
        "pair_geometric_mean_overlap",
        "pair_minimum_overlap",
    ]

    write_csv(
        PAIR_SUMMARY_CSV,
        pair_rows,
        pair_fieldnames,
    )

    norms = [
        float(row["orbital_norm"])
        for row in cube_audit_rows
    ]

    absolute_overlaps = [
        float(row["absolute_overlap"])
        for row in direct_rows
    ]

    tracked_pair_overlaps = [
        float(row["pair_geometric_mean_overlap"])
        for row in pair_rows
        if bool(row["is_tracked_root"])
    ]

    alternate_pair_overlaps = [
        float(row["pair_geometric_mean_overlap"])
        for row in pair_rows
        if not bool(row["is_tracked_root"])
    ]

    grids_ok = all(
        bool(row["same_grid"])
        for row in direct_rows
    )

    geometries_ok = all(
        bool(row["same_geometry"])
        for row in direct_rows
    )

    finite_norms = all(
        math.isfinite(value) and value > 0.0
        for value in norms
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 cube-grid and direct vacuum/embedding overlap audit\n\n"
        )

        handle.write("## Scope\n\n")
        handle.write(
            "- Parsed and numerically integrated all 48 selected NTO cubes, "
            "including ORCA molecular-orbital dataset-ID records.\n"
        )
        handle.write(
            "- Verified the cube-grid and frozen-geometry identity for "
            "the 24 vacuum/embedding orbital comparisons available for "
            "PYR2 and PYR5.\n"
        )
        handle.write(
            "- Computed sign-invariant normalized orbital overlaps without "
            "interpolation only where grids and atom records were identical.\n\n"
        )

        handle.write("## Integrity results\n\n")
        handle.write(
            f"- Cubes parsed successfully: {len(cube_audit_rows)}/48\n"
        )
        handle.write(
            f"- Finite positive orbital norms: "
            f"{sum(math.isfinite(v) and v > 0.0 for v in norms)}/48\n"
        )
        handle.write(
            f"- Orbital-norm range: {min(norms):.10f} to "
            f"{max(norms):.10f}\n"
        )
        handle.write(
            f"- Direct vacuum/embedding comparisons: "
            f"{len(direct_rows)}/24\n"
        )
        handle.write(
            f"- Identical grids: "
            f"{sum(bool(row['same_grid']) for row in direct_rows)}/24\n"
        )
        handle.write(
            f"- Identical frozen atom records: "
            f"{sum(bool(row['same_geometry']) for row in direct_rows)}/24\n\n"
        )

        handle.write("## Direct orbital-shape results\n\n")
        handle.write(
            f"- Absolute-overlap range across 24 orbitals: "
            f"{min(absolute_overlaps):.8f} to "
            f"{max(absolute_overlaps):.8f}\n"
        )
        handle.write(
            f"- Tracked-pair geometric-mean overlap range: "
            f"{min(tracked_pair_overlaps):.8f} to "
            f"{max(tracked_pair_overlaps):.8f}\n"
        )
        handle.write(
            f"- Alternate-pair geometric-mean overlap range: "
            f"{min(alternate_pair_overlaps):.8f} to "
            f"{max(alternate_pair_overlaps):.8f}\n\n"
        )

        handle.write("## Pair-resolved comparison\n\n")
        handle.write(
            "| Frame | Site | Root | Tracked | Pair | Hole | Particle | "
            "Geometric mean | Minimum |\n"
        )
        handle.write(
            "|---:|---|---:|---:|---|---:|---:|---:|---:|\n"
        )

        for row in pair_rows:
            handle.write(
                f"| {row['frame']} "
                f"| {row['cluster']} "
                f"| S{row['root']} "
                f"| {row['is_tracked_root']} "
                f"| `{row['pair_label']}` "
                f"| {float(row['hole_absolute_overlap']):.8f} "
                f"| {float(row['particle_absolute_overlap']):.8f} "
                f"| {float(row['pair_geometric_mean_overlap']):.8f} "
                f"| {float(row['pair_minimum_overlap']):.8f} |\n"
            )

        handle.write("\n## Acceptance\n\n")

        if (
            len(cube_audit_rows) == 48
            and len(direct_rows) == 24
            and len(pair_rows) == 12
            and grids_ok
            and geometries_ok
            and finite_norms
        ):
            handle.write(
                "**Day019 cube-grid and direct-overlap audit: PASS.**\n\n"
            )
        else:
            handle.write(
                "**Day019 cube-grid and direct-overlap audit: FAIL.**\n\n"
            )

        handle.write(
            "Cross-site comparisons are deliberately excluded here because "
            "PYR2-PYR5 occupy different Cartesian frames. Those comparisons "
            "require atom-based rigid alignment and field interpolation; "
            "direct voxel-wise overlap would be invalid.\n"
        )

    log("")
    log("Day019 cube-grid and direct-overlap audit completed.")
    log(f"Cubes parsed: {len(cube_audit_rows)}/48")
    log(f"Direct overlaps: {len(direct_rows)}/24")
    log(f"Pair comparisons: {len(pair_rows)}/12")
    log(
        "Identical grids and frozen geometries: "
        f"{sum(bool(row['same_grid']) and bool(row['same_geometry']) for row in direct_rows)}/24"
    )
    log(
        f"Orbital norm range: "
        f"{min(norms):.10f}-{max(norms):.10f}"
    )
    log(
        f"Direct |overlap| range: "
        f"{min(absolute_overlaps):.8f}-{max(absolute_overlaps):.8f}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
