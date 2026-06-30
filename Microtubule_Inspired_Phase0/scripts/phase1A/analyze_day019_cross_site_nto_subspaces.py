#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CUBE_GENERATION_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_nto_cube_generation"
)

SELECTION_MANIFEST = (
    CUBE_GENERATION_ROOT / "CUBE_SELECTION_MANIFEST_DAY019.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_nto_cross_site_overlap_analysis"
)

DIRECTIONAL_CSV = OUTPUT_ROOT / "cross_site_directional_metrics.csv"
SYMMETRIC_CSV = OUTPUT_ROOT / "cross_site_symmetric_metrics.csv"
ALIGNMENT_CSV = OUTPUT_ROOT / "cross_site_kabsch_alignment.csv"
REPORT_MD = OUTPUT_ROOT / "DAY019_CROSS_SITE_NTO_SUBSPACE_ANALYSIS.md"

BOHR_TO_ANGSTROM = 0.529177210903
MIN_NORM = 1.0e-12
INTERPOLATION_CHUNK = 100_000


@dataclass(frozen=True)
class CubeMetadata:
    natoms_raw: int
    natoms: int
    origin: np.ndarray
    shape: tuple[int, int, int]
    vectors: np.ndarray
    atomic_numbers: np.ndarray
    atom_charges: np.ndarray
    atom_coordinates: np.ndarray
    voxel_volume: float
    n_datasets: int
    dataset_ids: tuple[int, ...]


@dataclass
class CubeData:
    metadata: CubeMetadata
    values: np.ndarray


def log(message: str = "") -> None:
    print(message, flush=True)


def bool_from_csv(value: str) -> bool:
    return value.strip().lower() == "true"


def read_selection_manifest() -> list[dict[str, str]]:
    if not SELECTION_MANIFEST.is_file():
        raise SystemExit(
            f"Missing selection manifest: {SELECTION_MANIFEST}"
        )

    with SELECTION_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    vacuum_rows = [
        row
        for row in rows
        if row["calculation_type"] == "vacuum_reference"
    ]

    if len(vacuum_rows) != 24:
        raise SystemExit(
            f"Expected 24 vacuum-reference cube rows, "
            f"found {len(vacuum_rows)}."
        )

    return vacuum_rows


def parse_int_like(value: str) -> int:
    return int(float(value))


def read_cube(path: Path) -> CubeData:
    if not path.is_file():
        raise RuntimeError(f"Missing cube file: {path}")

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        if not handle.readline() or not handle.readline():
            raise RuntimeError(f"Truncated cube comments: {path}")

        origin_line = handle.readline().split()
        if len(origin_line) < 4:
            raise RuntimeError(f"Invalid cube origin line: {path}")

        natoms_raw = parse_int_like(origin_line[0])
        natoms = abs(natoms_raw)

        origin = np.array(
            [
                float(origin_line[1]),
                float(origin_line[2]),
                float(origin_line[3]),
            ],
            dtype=np.float64,
        )

        shape: list[int] = []
        vectors: list[list[float]] = []

        for _ in range(3):
            line = handle.readline().split()
            if len(line) < 4:
                raise RuntimeError(f"Invalid cube grid line: {path}")

            shape.append(abs(parse_int_like(line[0])))
            vectors.append(
                [float(line[1]), float(line[2]), float(line[3])]
            )

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

        n_datasets = 1
        dataset_ids: list[int] = []

        if natoms_raw < 0:
            tokens: list[str] = []

            while len(tokens) < 1:
                line = handle.readline()
                if not line:
                    raise RuntimeError(
                        f"Missing DSET_IDS record: {path}"
                    )
                tokens.extend(line.split())

            n_datasets = parse_int_like(tokens[0])
            required = 1 + n_datasets

            while len(tokens) < required:
                line = handle.readline()
                if not line:
                    raise RuntimeError(
                        f"Truncated DSET_IDS record: {path}"
                    )
                tokens.extend(line.split())

            dataset_ids = [
                parse_int_like(token)
                for token in tokens[1:required]
            ]

            leftover_tokens = tokens[required:]
        else:
            leftover_tokens = []

        expected_values = int(np.prod(shape)) * n_datasets

        data_tokens = (
            leftover_tokens
            + [
                token
                for line in handle
                for token in line.split()
            ]
        )

        values = np.fromiter(
            (float(token) for token in data_tokens),
            dtype=np.float64,
        )

    if values.size != expected_values:
        raise RuntimeError(
            f"Cube data-count mismatch for {path}: "
            f"expected {expected_values}, found {values.size}"
        )

    if n_datasets != 1:
        raise RuntimeError(
            f"Expected one orbital dataset in {path}, "
            f"found {n_datasets}"
        )

    vectors_array = np.array(vectors, dtype=np.float64)
    voxel_volume = abs(float(np.linalg.det(vectors_array)))

    if voxel_volume <= 0.0 or not math.isfinite(voxel_volume):
        raise RuntimeError(
            f"Invalid voxel volume in {path}: {voxel_volume}"
        )

    values_3d = values.reshape(
        (shape[0], shape[1], shape[2]),
        order="C",
    )

    metadata = CubeMetadata(
        natoms_raw=natoms_raw,
        natoms=natoms,
        origin=origin,
        shape=(shape[0], shape[1], shape[2]),
        vectors=vectors_array,
        atomic_numbers=np.array(
            atomic_numbers,
            dtype=np.int64,
        ),
        atom_charges=np.array(
            atom_charges,
            dtype=np.float64,
        ),
        atom_coordinates=np.array(
            atom_coordinates,
            dtype=np.float64,
        ),
        voxel_volume=voxel_volume,
        n_datasets=n_datasets,
        dataset_ids=tuple(dataset_ids),
    )

    return CubeData(
        metadata=metadata,
        values=values_3d,
    )


def cube_norm(values: np.ndarray, voxel_volume: float) -> float:
    norm_squared = float(
        np.dot(values.ravel(), values.ravel())
        * voxel_volume
    )

    if norm_squared <= MIN_NORM or not math.isfinite(norm_squared):
        raise RuntimeError(
            f"Invalid field norm-squared: {norm_squared}"
        )

    return math.sqrt(norm_squared)


def kabsch_row_transform(
    moving: np.ndarray,
    reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    if moving.shape != reference.shape:
        raise RuntimeError(
            f"Kabsch shape mismatch: {moving.shape} vs {reference.shape}"
        )

    moving_center = moving.mean(axis=0)
    reference_center = reference.mean(axis=0)

    moving_centered = moving - moving_center
    reference_centered = reference - reference_center

    covariance = moving_centered.T @ reference_centered
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)

    rotation = u_matrix @ vt_matrix

    if np.linalg.det(rotation) < 0.0:
        u_matrix[:, -1] *= -1.0
        rotation = u_matrix @ vt_matrix

    translation = reference_center - moving_center @ rotation

    aligned = moving @ rotation + translation
    differences = aligned - reference

    rmsd = math.sqrt(
        float(np.mean(np.sum(differences * differences, axis=1)))
    )

    max_displacement = math.sqrt(
        float(np.max(np.sum(differences * differences, axis=1)))
    )

    return rotation, translation, rmsd, max_displacement


def reference_grid_points(metadata: CubeMetadata) -> np.ndarray:
    nx, ny, nz = metadata.shape

    indices = np.indices(
        (nx, ny, nz),
        dtype=np.float64,
    ).reshape(3, -1).T

    return metadata.origin + indices @ metadata.vectors


def trilinear_resample_moving_to_reference(
    moving: CubeData,
    reference_metadata: CubeMetadata,
    reference_points: np.ndarray,
    rotation_moving_to_reference: np.ndarray,
    translation_moving_to_reference: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    total_points = reference_points.shape[0]

    output = np.zeros(total_points, dtype=np.float64)
    valid_points = 0

    inverse_moving_vectors = np.linalg.inv(
        moving.metadata.vectors
    )

    nx, ny, nz = moving.metadata.shape
    moving_values = moving.values

    for start in range(0, total_points, INTERPOLATION_CHUNK):
        stop = min(start + INTERPOLATION_CHUNK, total_points)

        reference_chunk = reference_points[start:stop]

        moving_cartesian = (
            reference_chunk - translation_moving_to_reference
        ) @ rotation_moving_to_reference.T

        fractional = (
            moving_cartesian - moving.metadata.origin
        ) @ inverse_moving_vectors

        valid = np.all(
            (fractional >= 0.0)
            & (
                fractional
                <= np.array([nx - 1, ny - 1, nz - 1])
            ),
            axis=1,
        )

        valid_points += int(np.count_nonzero(valid))

        if not np.any(valid):
            continue

        frac_valid = fractional[valid]

        lower = np.floor(frac_valid).astype(np.int64)
        upper = np.minimum(
            lower + 1,
            np.array([nx - 1, ny - 1, nz - 1]),
        )

        delta = frac_valid - lower

        x0, y0, z0 = lower.T
        x1, y1, z1 = upper.T
        dx, dy, dz = delta.T

        c000 = moving_values[x0, y0, z0]
        c001 = moving_values[x0, y0, z1]
        c010 = moving_values[x0, y1, z0]
        c011 = moving_values[x0, y1, z1]
        c100 = moving_values[x1, y0, z0]
        c101 = moving_values[x1, y0, z1]
        c110 = moving_values[x1, y1, z0]
        c111 = moving_values[x1, y1, z1]

        c00 = c000 * (1.0 - dz) + c001 * dz
        c01 = c010 * (1.0 - dz) + c011 * dz
        c10 = c100 * (1.0 - dz) + c101 * dz
        c11 = c110 * (1.0 - dz) + c111 * dz

        c0 = c00 * (1.0 - dy) + c01 * dy
        c1 = c10 * (1.0 - dy) + c11 * dy

        interpolated = c0 * (1.0 - dx) + c1 * dx

        chunk_output = output[start:stop]
        chunk_output[valid] = interpolated
        output[start:stop] = chunk_output

    captured_norm = cube_norm(
        output,
        reference_metadata.voxel_volume,
    )

    coverage_fraction = valid_points / total_points

    return (
        output.reshape(reference_metadata.shape, order="C"),
        captured_norm,
        coverage_fraction,
    )


def weighted_orthonormal_basis(
    vectors: list[np.ndarray],
    voxel_volume: float,
) -> np.ndarray:
    matrix = np.column_stack(
        [vector.ravel() for vector in vectors]
    )

    weighted = matrix * math.sqrt(voxel_volume)
    q_matrix, r_matrix = np.linalg.qr(weighted)

    diagonal = np.abs(np.diag(r_matrix))

    if np.any(diagonal <= 1.0e-10):
        raise RuntimeError(
            f"Linearly dependent orbital subspace; QR diagonal={diagonal}"
        )

    return q_matrix


def subspace_principal_cosines(
    reference_vectors: list[np.ndarray],
    moving_vectors_on_reference_grid: list[np.ndarray],
    voxel_volume: float,
) -> np.ndarray:
    if len(reference_vectors) != len(
        moving_vectors_on_reference_grid
    ):
        raise RuntimeError(
            "Subspace dimensions do not match."
        )

    reference_basis = weighted_orthonormal_basis(
        reference_vectors,
        voxel_volume,
    )

    moving_basis = weighted_orthonormal_basis(
        moving_vectors_on_reference_grid,
        voxel_volume,
    )

    overlap_matrix = reference_basis.T @ moving_basis
    singular_values = np.linalg.svd(
        overlap_matrix,
        compute_uv=False,
    )

    return np.clip(singular_values, 0.0, 1.0)


def rms_principal_similarity(
    singular_values: np.ndarray,
) -> float:
    return math.sqrt(
        float(np.mean(singular_values * singular_values))
    )


def prepare_site_records(
    rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    by_site: dict[str, list[dict[str, str]]] = {}

    for row in rows:
        by_site.setdefault(row["cluster"], []).append(row)

    expected_sites = {"PYR2", "PYR3", "PYR4", "PYR5"}

    if set(by_site) != expected_sites:
        raise RuntimeError(
            f"Unexpected vacuum sites: {sorted(by_site)}"
        )

    records: dict[str, dict[str, object]] = {}

    for site, site_rows in by_site.items():
        jobs = {row["job"] for row in site_rows}

        if len(jobs) != 1:
            raise RuntimeError(
                f"Expected one vacuum job for {site}, found {jobs}"
            )

        tracked_roots = {
            int(row["root"])
            for row in site_rows
            if bool_from_csv(row["is_tracked_root"])
        }

        alternate_roots = {
            int(row["root"])
            for row in site_rows
            if not bool_from_csv(row["is_tracked_root"])
        }

        if len(tracked_roots) != 1 or len(alternate_roots) != 1:
            raise RuntimeError(
                f"Invalid tracked/alternate roots for {site}"
            )

        tracked_root = next(iter(tracked_roots))
        alternate_root = next(iter(alternate_roots))

        def select(
            root: int,
            role: str,
        ) -> list[dict[str, str]]:
            selected = [
                row
                for row in site_rows
                if int(row["root"]) == root
                and row["orbital_role"] == role
            ]

            selected.sort(
                key=lambda row: int(row["pair_rank"])
            )

            return selected

        tracked_hole = select(tracked_root, "hole")
        tracked_particle = select(tracked_root, "particle")
        alternate_hole = select(alternate_root, "hole")
        alternate_particle = select(alternate_root, "particle")

        if (
            len(tracked_hole) != 1
            or len(tracked_particle) != 1
            or len(alternate_hole) != 2
            or len(alternate_particle) != 2
        ):
            raise RuntimeError(
                f"Unexpected selected orbital counts for {site}"
            )

        all_paths = [
            row["cube_file"]
            for row in site_rows
        ]

        records[site] = {
            "job": next(iter(jobs)),
            "tracked_root": tracked_root,
            "alternate_root": alternate_root,
            "tracked_hole": tracked_hole,
            "tracked_particle": tracked_particle,
            "alternate_hole": alternate_hole,
            "alternate_particle": alternate_particle,
            "all_cube_paths": all_paths,
        }

    return records


def load_cube_cache(
    rows: list[dict[str, str]],
) -> dict[str, CubeData]:
    cache: dict[str, CubeData] = {}

    unique_paths = sorted(
        {row["cube_file"] for row in rows}
    )

    for index, relative_path in enumerate(
        unique_paths,
        start=1,
    ):
        path = (PROJECT_ROOT / relative_path).resolve()
        cache[relative_path] = read_cube(path)

        log(
            f"[load {index:02d}/{len(unique_paths)}] "
            f"{path.name}"
        )

    return cache


def validate_site_internal_metadata(
    records: dict[str, dict[str, object]],
    cache: dict[str, CubeData],
) -> None:
    for site, record in records.items():
        paths = list(record["all_cube_paths"])
        reference = cache[paths[0]].metadata

        for relative_path in paths[1:]:
            metadata = cache[relative_path].metadata

            if metadata.shape != reference.shape:
                raise RuntimeError(
                    f"Inconsistent cube shape within {site}"
                )

            if not np.allclose(
                metadata.origin,
                reference.origin,
                atol=1.0e-10,
                rtol=0.0,
            ):
                raise RuntimeError(
                    f"Inconsistent cube origin within {site}"
                )

            if not np.allclose(
                metadata.vectors,
                reference.vectors,
                atol=1.0e-10,
                rtol=0.0,
            ):
                raise RuntimeError(
                    f"Inconsistent cube vectors within {site}"
                )

            if not np.array_equal(
                metadata.atomic_numbers,
                reference.atomic_numbers,
            ):
                raise RuntimeError(
                    f"Inconsistent atom identities within {site}"
                )

            if not np.allclose(
                metadata.atom_coordinates,
                reference.atom_coordinates,
                atol=1.0e-10,
                rtol=0.0,
            ):
                raise RuntimeError(
                    f"Inconsistent atom coordinates within {site}"
                )


def rows_to_vectors(
    rows: list[dict[str, str]],
    cache: dict[str, CubeData],
) -> list[np.ndarray]:
    return [
        cache[row["cube_file"]].values
        for row in rows
    ]


def directional_compare(
    reference_site: str,
    moving_site: str,
    records: dict[str, dict[str, object]],
    cache: dict[str, CubeData],
) -> tuple[dict[str, object], dict[str, object]]:
    reference_record = records[reference_site]
    moving_record = records[moving_site]

    reference_any = cache[
        reference_record["all_cube_paths"][0]
    ]
    moving_any = cache[
        moving_record["all_cube_paths"][0]
    ]

    if not np.array_equal(
        reference_any.metadata.atomic_numbers,
        moving_any.metadata.atomic_numbers,
    ):
        raise RuntimeError(
            f"Atom identity mismatch: {reference_site} vs {moving_site}"
        )

    heavy_mask = reference_any.metadata.atomic_numbers == 6

    if int(np.count_nonzero(heavy_mask)) != 16:
        raise RuntimeError(
            f"Expected 16 heavy carbon atoms for {reference_site}"
        )

    rotation, translation, heavy_rmsd_bohr, heavy_max_bohr = (
        kabsch_row_transform(
            moving_any.metadata.atom_coordinates[heavy_mask],
            reference_any.metadata.atom_coordinates[heavy_mask],
        )
    )

    aligned_all_atoms = (
        moving_any.metadata.atom_coordinates @ rotation
        + translation
    )

    all_atom_diff = (
        aligned_all_atoms
        - reference_any.metadata.atom_coordinates
    )

    all_atom_rmsd_bohr = math.sqrt(
        float(
            np.mean(
                np.sum(all_atom_diff * all_atom_diff, axis=1)
            )
        )
    )

    all_atom_max_bohr = math.sqrt(
        float(
            np.max(
                np.sum(all_atom_diff * all_atom_diff, axis=1)
            )
        )
    )

    groups = [
        "tracked_hole",
        "tracked_particle",
        "alternate_hole",
        "alternate_particle",
    ]

    reference_points = reference_grid_points(
        reference_any.metadata
    )

    resampled: dict[str, list[np.ndarray]] = {}
    captured_norms: list[float] = []
    coverage_fractions: list[float] = []

    for group in groups:
        moving_rows = moving_record[group]
        resampled[group] = []

        for row in moving_rows:
            moving_cube = cache[row["cube_file"]]

            (
                interpolated,
                captured_norm,
                coverage_fraction,
            ) = trilinear_resample_moving_to_reference(
                moving=moving_cube,
                reference_metadata=reference_any.metadata,
                reference_points=reference_points,
                rotation_moving_to_reference=rotation,
                translation_moving_to_reference=translation,
            )

            resampled[group].append(interpolated)
            captured_norms.append(captured_norm)
            coverage_fractions.append(coverage_fraction)

            log(
                f"    {moving_site}->{reference_site} "
                f"S{row['root']} "
                f"MO{row['orbital_index']}{row['orbital_spin']} "
                f"captured_norm={captured_norm:.8f} "
                f"coverage={coverage_fraction:.6f}"
            )

    voxel_volume = reference_any.metadata.voxel_volume

    tracked_hole_sv = subspace_principal_cosines(
        rows_to_vectors(
            reference_record["tracked_hole"],
            cache,
        ),
        resampled["tracked_hole"],
        voxel_volume,
    )

    tracked_particle_sv = subspace_principal_cosines(
        rows_to_vectors(
            reference_record["tracked_particle"],
            cache,
        ),
        resampled["tracked_particle"],
        voxel_volume,
    )

    alternate_hole_sv = subspace_principal_cosines(
        rows_to_vectors(
            reference_record["alternate_hole"],
            cache,
        ),
        resampled["alternate_hole"],
        voxel_volume,
    )

    alternate_particle_sv = subspace_principal_cosines(
        rows_to_vectors(
            reference_record["alternate_particle"],
            cache,
        ),
        resampled["alternate_particle"],
        voxel_volume,
    )

    tracked_hole_similarity = float(tracked_hole_sv[0])
    tracked_particle_similarity = float(
        tracked_particle_sv[0]
    )

    tracked_pair_similarity = math.sqrt(
        tracked_hole_similarity
        * tracked_particle_similarity
    )

    alternate_hole_rms = rms_principal_similarity(
        alternate_hole_sv
    )

    alternate_particle_rms = rms_principal_similarity(
        alternate_particle_sv
    )

    alternate_transition_similarity = math.sqrt(
        alternate_hole_rms
        * alternate_particle_rms
    )

    metrics = {
        "reference_site": reference_site,
        "moving_site": moving_site,
        "reference_tracked_root": reference_record["tracked_root"],
        "moving_tracked_root": moving_record["tracked_root"],
        "reference_alternate_root": reference_record["alternate_root"],
        "moving_alternate_root": moving_record["alternate_root"],
        "heavy_atom_rmsd_A": heavy_rmsd_bohr * BOHR_TO_ANGSTROM,
        "heavy_atom_max_displacement_A": (
            heavy_max_bohr * BOHR_TO_ANGSTROM
        ),
        "all_atom_rmsd_A": all_atom_rmsd_bohr * BOHR_TO_ANGSTROM,
        "all_atom_max_displacement_A": (
            all_atom_max_bohr * BOHR_TO_ANGSTROM
        ),
        "minimum_captured_norm": min(captured_norms),
        "maximum_captured_norm": max(captured_norms),
        "minimum_grid_coverage_fraction": min(
            coverage_fractions
        ),
        "tracked_hole_similarity": tracked_hole_similarity,
        "tracked_particle_similarity": tracked_particle_similarity,
        "tracked_pair_geometric_mean_similarity": (
            tracked_pair_similarity
        ),
        "alternate_hole_principal_cosine_1": float(
            alternate_hole_sv[0]
        ),
        "alternate_hole_principal_cosine_2": float(
            alternate_hole_sv[1]
        ),
        "alternate_hole_subspace_rms_similarity": (
            alternate_hole_rms
        ),
        "alternate_particle_principal_cosine_1": float(
            alternate_particle_sv[0]
        ),
        "alternate_particle_principal_cosine_2": float(
            alternate_particle_sv[1]
        ),
        "alternate_particle_subspace_rms_similarity": (
            alternate_particle_rms
        ),
        "alternate_transition_subspace_geometric_mean_similarity": (
            alternate_transition_similarity
        ),
        "alternate_minimum_principal_cosine": min(
            float(alternate_hole_sv[-1]),
            float(alternate_particle_sv[-1]),
        ),
    }

    alignment = {
        "reference_site": reference_site,
        "moving_site": moving_site,
        "heavy_atom_rmsd_A": metrics["heavy_atom_rmsd_A"],
        "heavy_atom_max_displacement_A": metrics[
            "heavy_atom_max_displacement_A"
        ],
        "all_atom_rmsd_A": metrics["all_atom_rmsd_A"],
        "all_atom_max_displacement_A": metrics[
            "all_atom_max_displacement_A"
        ],
        "rotation_00": rotation[0, 0],
        "rotation_01": rotation[0, 1],
        "rotation_02": rotation[0, 2],
        "rotation_10": rotation[1, 0],
        "rotation_11": rotation[1, 1],
        "rotation_12": rotation[1, 2],
        "rotation_20": rotation[2, 0],
        "rotation_21": rotation[2, 1],
        "rotation_22": rotation[2, 2],
        "translation_x_bohr": translation[0],
        "translation_y_bohr": translation[1],
        "translation_z_bohr": translation[2],
    }

    return metrics, alignment


def average_directional_rows(
    first: dict[str, object],
    second: dict[str, object],
) -> dict[str, object]:
    site_a = str(first["reference_site"])
    site_b = str(first["moving_site"])

    if {
        site_a,
        site_b,
    } != {
        str(second["reference_site"]),
        str(second["moving_site"]),
    }:
        raise RuntimeError(
            "Directional rows do not describe the same site pair."
        )

    numeric_fields = [
        "heavy_atom_rmsd_A",
        "heavy_atom_max_displacement_A",
        "all_atom_rmsd_A",
        "all_atom_max_displacement_A",
        "minimum_captured_norm",
        "maximum_captured_norm",
        "minimum_grid_coverage_fraction",
        "tracked_hole_similarity",
        "tracked_particle_similarity",
        "tracked_pair_geometric_mean_similarity",
        "alternate_hole_principal_cosine_1",
        "alternate_hole_principal_cosine_2",
        "alternate_hole_subspace_rms_similarity",
        "alternate_particle_principal_cosine_1",
        "alternate_particle_principal_cosine_2",
        "alternate_particle_subspace_rms_similarity",
        "alternate_transition_subspace_geometric_mean_similarity",
        "alternate_minimum_principal_cosine",
    ]

    result: dict[str, object] = {
        "site_a": min(site_a, site_b),
        "site_b": max(site_a, site_b),
    }

    for field in numeric_fields:
        first_value = float(first[field])
        second_value = float(second[field])

        result[f"{field}_mean"] = (
            first_value + second_value
        ) / 2.0

        result[f"{field}_directional_difference"] = abs(
            first_value - second_value
        )

        result[f"{field}_minimum"] = min(
            first_value,
            second_value,
        )

    return result


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    fieldnames = list(rows[0].keys())

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


def main() -> None:
    vacuum_rows = read_selection_manifest()
    records = prepare_site_records(vacuum_rows)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    log("Day019 cross-site NTO aligned-subspace analysis")
    log("Vacuum sites: PYR2, PYR3, PYR4, PYR5")
    log("Loading 24 selected vacuum-reference cubes.")

    cache = load_cube_cache(vacuum_rows)
    validate_site_internal_metadata(records, cache)

    directional_rows: list[dict[str, object]] = []
    alignment_rows: list[dict[str, object]] = []
    symmetric_rows: list[dict[str, object]] = []

    site_pairs = list(
        combinations(
            ["PYR2", "PYR3", "PYR4", "PYR5"],
            2,
        )
    )

    for pair_index, (site_a, site_b) in enumerate(
        site_pairs,
        start=1,
    ):
        log("")
        log("=" * 88)
        log(
            f"[pair {pair_index}/6] "
            f"{site_a} <-> {site_b}"
        )
        log("=" * 88)

        first, first_alignment = directional_compare(
            reference_site=site_a,
            moving_site=site_b,
            records=records,
            cache=cache,
        )

        second, second_alignment = directional_compare(
            reference_site=site_b,
            moving_site=site_a,
            records=records,
            cache=cache,
        )

        directional_rows.extend([first, second])
        alignment_rows.extend(
            [first_alignment, second_alignment]
        )

        symmetric = average_directional_rows(
            first,
            second,
        )

        symmetric_rows.append(symmetric)

        log(
            f"  symmetric tracked-pair similarity="
            f"{float(symmetric['tracked_pair_geometric_mean_similarity_mean']):.8f}"
        )
        log(
            f"  symmetric alternate-subspace similarity="
            f"{float(symmetric['alternate_transition_subspace_geometric_mean_similarity_mean']):.8f}"
        )

    write_csv(DIRECTIONAL_CSV, directional_rows)
    write_csv(SYMMETRIC_CSV, symmetric_rows)
    write_csv(ALIGNMENT_CSV, alignment_rows)

    tracked_values = [
        float(
            row[
                "tracked_pair_geometric_mean_similarity_mean"
            ]
        )
        for row in symmetric_rows
    ]

    alternate_values = [
        float(
            row[
                "alternate_transition_subspace_geometric_mean_similarity_mean"
            ]
        )
        for row in symmetric_rows
    ]

    minimum_captured_norm = min(
        float(row["minimum_captured_norm_minimum"])
        for row in symmetric_rows
    )

    minimum_coverage = min(
        float(
            row[
                "minimum_grid_coverage_fraction_minimum"
            ]
        )
        for row in symmetric_rows
    )

    maximum_directional_asymmetry = max(
        max(
            float(
                row[
                    "tracked_pair_geometric_mean_similarity_directional_difference"
                ]
            ),
            float(
                row[
                    "alternate_transition_subspace_geometric_mean_similarity_directional_difference"
                ]
            ),
        )
        for row in symmetric_rows
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 cross-site NTO aligned-subspace analysis\n\n"
        )

        handle.write("## Method\n\n")
        handle.write(
            "- The four vacuum-reference chromophores were compared "
            "pairwise.\n"
        )
        handle.write(
            "- Moving geometries were aligned to reference geometries by "
            "a proper-rotation Kabsch fit over the 16 carbon atoms.\n"
        )
        handle.write(
            "- Moving orbital fields were evaluated on the reference cube "
            "grid by trilinear interpolation.\n"
        )
        handle.write(
            "- Both directions were calculated for every site pair and "
            "averaged to expose interpolation asymmetry.\n"
        )
        handle.write(
            "- The tracked bright state was compared as a one-dimensional "
            "hole/particle pair.\n"
        )
        handle.write(
            "- The alternate low state was compared as separate "
            "two-dimensional hole and particle subspaces. This is invariant "
            "to arbitrary rotations between near-degenerate NTO pairs.\n\n"
        )

        handle.write("## Numerical controls\n\n")
        handle.write(
            f"- Pairwise site comparisons: {len(symmetric_rows)}/6\n"
        )
        handle.write(
            f"- Directional calculations: {len(directional_rows)}/12\n"
        )
        handle.write(
            f"- Minimum interpolated-orbital captured norm: "
            f"{minimum_captured_norm:.8f}\n"
        )
        handle.write(
            f"- Minimum reference-grid coverage fraction: "
            f"{minimum_coverage:.8f}\n"
        )
        handle.write(
            f"- Maximum directional similarity asymmetry: "
            f"{maximum_directional_asymmetry:.8f}\n\n"
        )

        handle.write("## Symmetric cross-site results\n\n")
        handle.write(
            "| Site A | Site B | Heavy RMSD (Ã) | Tracked hole | "
            "Tracked particle | Tracked pair | Alternate hole subspace | "
            "Alternate particle subspace | Alternate transition subspace | "
            "Minimum principal cosine |\n"
        )
        handle.write(
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n"
        )

        for row in symmetric_rows:
            handle.write(
                f"| {row['site_a']} "
                f"| {row['site_b']} "
                f"| {float(row['heavy_atom_rmsd_A_mean']):.6f} "
                f"| {float(row['tracked_hole_similarity_mean']):.8f} "
                f"| {float(row['tracked_particle_similarity_mean']):.8f} "
                f"| {float(row['tracked_pair_geometric_mean_similarity_mean']):.8f} "
                f"| {float(row['alternate_hole_subspace_rms_similarity_mean']):.8f} "
                f"| {float(row['alternate_particle_subspace_rms_similarity_mean']):.8f} "
                f"| {float(row['alternate_transition_subspace_geometric_mean_similarity_mean']):.8f} "
                f"| {float(row['alternate_minimum_principal_cosine_minimum']):.8f} |\n"
            )

        handle.write("\n## Aggregate ranges\n\n")
        handle.write(
            f"- Tracked-pair cross-site similarity range: "
            f"{min(tracked_values):.8f} to {max(tracked_values):.8f}\n"
        )
        handle.write(
            f"- Alternate-transition-subspace similarity range: "
            f"{min(alternate_values):.8f} to "
            f"{max(alternate_values):.8f}\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "This analysis establishes spatial similarity after rigid "
            "alignment. It does not by itself provide diabatic state phases "
            "or interstate couplings. The tracked-root and alternate-root "
            "spaces should only be frozen into a Hamiltonian after these "
            "overlap results are reviewed together with the S1-S2 energy "
            "splittings, oscillator strengths, and subsequent coupling "
            "calculations.\n"
        )

    log("")
    log("Day019 cross-site NTO analysis completed.")
    log(f"Site pairs analyzed: {len(symmetric_rows)}/6")
    log(f"Directional calculations: {len(directional_rows)}/12")
    log(
        f"Tracked-pair similarity range: "
        f"{min(tracked_values):.8f}-{max(tracked_values):.8f}"
    )
    log(
        f"Alternate-subspace similarity range: "
        f"{min(alternate_values):.8f}-{max(alternate_values):.8f}"
    )
    log(
        f"Minimum captured norm: {minimum_captured_norm:.8f}"
    )
    log(
        f"Maximum directional asymmetry: "
        f"{maximum_directional_asymmetry:.8f}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
