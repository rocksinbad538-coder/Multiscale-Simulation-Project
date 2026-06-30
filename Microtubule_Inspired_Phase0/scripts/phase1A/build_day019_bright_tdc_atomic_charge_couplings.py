#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from itertools import combinations
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

TDENS_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_transition_density_N80"
)

TDENS_METRICS = (
    TDENS_ROOT
    / "BRIGHT_TRANSITION_DENSITY_N80_METRICS_DAY019.csv"
)

POINT_COUPLINGS = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_point_transition_dipole_couplings/"
    "point_dipole_couplings_long.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_tdc_atomic_charge_couplings"
)

CHARGES_CSV = OUTPUT_ROOT / "transition_density_atomic_charges.csv"
SITE_VALIDATION_CSV = OUTPUT_ROOT / "site_charge_validation.csv"
COUPLINGS_CSV = OUTPUT_ROOT / "tdc_atomic_charge_couplings_frame000.csv"
POINT_MATRIX_CSV = OUTPUT_ROOT / "point_dipole_bright_matrix_meV.csv"
TDCAC_MATRIX_CSV = OUTPUT_ROOT / "tdc_atomic_charge_matrix_meV.csv"
DELTA_MATRIX_CSV = OUTPUT_ROOT / "tdcac_minus_point_matrix_meV.csv"
REPORT_MD = OUTPUT_ROOT / "TDC_ATOMIC_CHARGE_COUPLING_BENCHMARK_DAY019.md"

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
SITE_INDEX = {site: index for index, site in enumerate(SITES)}
ATOMIC_SYMBOL = {1: "H", 6: "C"}

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_MEV = HARTREE_TO_EV * 1000.0
HARTREE_TO_CM1 = 219474.6313705

CHUNK_SIZE = 20000
PARTITION_CHARGE_TOL = 1.0e-10
CONSTRAINT_TOL = 1.0e-10
POINT_REPRODUCTION_TOL_MEV = 1.0e-5


def log(message: str = "") -> None:
    print(message, flush=True)


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    if fieldnames is None:
        fieldnames = list(rows[0].keys())

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_matrix_csv(
    path: Path,
    labels: tuple[str, ...],
    matrix: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["site", *labels])

        for label, row in zip(labels, matrix):
            writer.writerow(
                [label, *[f"{float(value):.12g}" for value in row]]
            )


def parse_cube(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise SystemExit(f"Missing cube file: {path}")

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        comment_1 = handle.readline().rstrip("\n")
        comment_2 = handle.readline().rstrip("\n")

        tokens = handle.readline().split()

        if len(tokens) < 4:
            raise RuntimeError(f"Invalid cube origin line: {path}")

        natoms_raw = int(tokens[0])
        natoms = abs(natoms_raw)
        origin = np.array(
            [float(tokens[1]), float(tokens[2]), float(tokens[3])],
            dtype=np.float64,
        )

        counts: list[int] = []
        axes: list[np.ndarray] = []

        for _ in range(3):
            tokens = handle.readline().split()

            if len(tokens) < 4:
                raise RuntimeError(f"Invalid cube axis line: {path}")

            counts.append(abs(int(tokens[0])))
            axes.append(
                np.array(
                    [float(tokens[1]), float(tokens[2]), float(tokens[3])],
                    dtype=np.float64,
                )
            )

        atom_records: list[tuple[int, float, float, float, float]] = []

        for _ in range(natoms):
            tokens = handle.readline().split()

            if len(tokens) < 5:
                raise RuntimeError(f"Invalid cube atom line: {path}")

            atom_records.append(
                (
                    int(float(tokens[0])),
                    float(tokens[1]),
                    float(tokens[2]),
                    float(tokens[3]),
                    float(tokens[4]),
                )
            )

        if natoms_raw < 0:
            handle.readline()

        values: list[float] = []

        for line in handle:
            values.extend(float(token) for token in line.split())

    nx, ny, nz = counts
    expected_values = nx * ny * nz

    if len(values) != expected_values:
        raise RuntimeError(
            f"{path}: expected {expected_values} values, "
            f"found {len(values)}."
        )

    axes_array = np.vstack(axes)
    voxel_volume = abs(float(np.linalg.det(axes_array)))

    if voxel_volume <= 0.0:
        raise RuntimeError(f"Nonpositive voxel volume: {path}")

    return {
        "comment_1": comment_1,
        "comment_2": comment_2,
        "natoms": natoms,
        "origin": origin,
        "counts": np.asarray(counts, dtype=np.int64),
        "axes": axes_array,
        "atom_records": atom_records,
        "density": np.asarray(values, dtype=np.float64),
        "voxel_volume": voxel_volume,
    }


def read_density_manifest() -> dict[str, dict[str, object]]:
    if not TDENS_METRICS.is_file():
        raise SystemExit(f"Missing transition-density metrics: {TDENS_METRICS}")

    data: dict[str, dict[str, object]] = {}

    with TDENS_METRICS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            site = row["site"]

            if site not in SITES:
                continue

            if not as_bool(row["overall_site_pass"]):
                raise RuntimeError(
                    f"Transition-density validation did not pass for {site}."
                )

            vector = np.array(
                [
                    float(row["ORCA_mu_x_au"]),
                    float(row["ORCA_mu_y_au"]),
                    float(row["ORCA_mu_z_au"]),
                ],
                dtype=np.float64,
            )

            data[site] = {
                "root": int(row["bright_root"]),
                "cube_path": (
                    PROJECT_ROOT / row["normalized_cube"]
                ).resolve(),
                "orca_mu_au": vector,
                "orca_mu_magnitude_au": float(
                    row["ORCA_mu_magnitude_au"]
                ),
                "normalized_cube_mu_magnitude_au": float(
                    row["sqrt2_normalized_mu_magnitude_au"]
                ),
                "density_scale": int(row["gauge_sign"]) * math.sqrt(2.0),
            }

    if set(data) != set(SITES):
        raise RuntimeError(
            f"Expected normalized cubes for {SITES}, found {sorted(data)}."
        )

    return data


def read_frame000_point_couplings() -> dict[tuple[str, str], float]:
    if not POINT_COUPLINGS.is_file():
        raise SystemExit(f"Missing point-dipole couplings: {POINT_COUPLINGS}")

    data: dict[tuple[str, str], float] = {}

    with POINT_COUPLINGS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["frame"]) != 0:
                continue

            if (
                row["family_a"] != "bright_like"
                or row["family_b"] != "bright_like"
            ):
                continue

            key = (row["site_a"], row["site_b"])
            data[key] = float(row["J_meV"])

    expected = set(combinations(SITES, 2))

    if set(data) != expected:
        raise RuntimeError(
            f"Expected six frame000 bright-bright point couplings; "
            f"found {sorted(data)}."
        )

    return data


def voxel_coordinates(
    start: int,
    stop: int,
    counts: np.ndarray,
    origin: np.ndarray,
    axes: np.ndarray,
) -> np.ndarray:
    nx, ny, nz = [int(value) for value in counts]
    flat = np.arange(start, stop, dtype=np.int64)

    i = flat // (ny * nz)
    remainder = flat % (ny * nz)
    j = remainder // nz
    k = remainder % nz

    coordinates = (
        origin[None, :]
        + i[:, None] * axes[0][None, :]
        + j[:, None] * axes[1][None, :]
        + k[:, None] * axes[2][None, :]
    )

    return coordinates


def partition_to_nearest_atoms(
    cube: dict[str, object],
) -> tuple[np.ndarray, float]:
    atom_records = list(cube["atom_records"])
    atom_coordinates = np.array(
        [[record[2], record[3], record[4]] for record in atom_records],
        dtype=np.float64,
    )

    density = np.asarray(cube["density"], dtype=np.float64)
    counts = np.asarray(cube["counts"], dtype=np.int64)
    origin = np.asarray(cube["origin"], dtype=np.float64)
    axes = np.asarray(cube["axes"], dtype=np.float64)
    voxel_volume = float(cube["voxel_volume"])

    charges = np.zeros(len(atom_records), dtype=np.float64)

    for start in range(0, density.size, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, density.size)
        coordinates = voxel_coordinates(
            start=start,
            stop=stop,
            counts=counts,
            origin=origin,
            axes=axes,
        )

        differences = (
            coordinates[:, None, :]
            - atom_coordinates[None, :, :]
        )
        distances_squared = np.einsum(
            "caj,caj->ca",
            differences,
            differences,
            optimize=True,
        )
        nearest = np.argmin(distances_squared, axis=1)
        voxel_charges = density[start:stop] * voxel_volume

        charges += np.bincount(
            nearest,
            weights=voxel_charges,
            minlength=len(atom_records),
        )

    cube_integral = float(np.sum(density) * voxel_volume)
    return charges, cube_integral


def constrained_charge_correction(
    raw_charges: np.ndarray,
    atom_coordinates: np.ndarray,
    carbon_centroid: np.ndarray,
    target_mu: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    relative_coordinates = atom_coordinates - carbon_centroid[None, :]

    constraints = np.vstack(
        [
            np.ones(raw_charges.size, dtype=np.float64),
            relative_coordinates[:, 0],
            relative_coordinates[:, 1],
            relative_coordinates[:, 2],
        ]
    )

    target = np.concatenate(
        [
            np.array([0.0], dtype=np.float64),
            np.asarray(target_mu, dtype=np.float64),
        ]
    )

    residual = target - constraints @ raw_charges
    gram = constraints @ constraints.T

    correction = constraints.T @ np.linalg.solve(gram, residual)
    corrected = raw_charges + correction

    final_residual = constraints @ corrected - target

    return corrected, final_residual


def dipole_from_charges(
    charges: np.ndarray,
    atom_coordinates: np.ndarray,
    origin: np.ndarray,
) -> np.ndarray:
    relative = atom_coordinates - origin[None, :]
    return np.sum(charges[:, None] * relative, axis=0)


def point_dipole_coupling_meV(
    mu_a: np.ndarray,
    mu_b: np.ndarray,
    center_a: np.ndarray,
    center_b: np.ndarray,
) -> float:
    displacement = center_b - center_a
    distance = float(np.linalg.norm(displacement))
    rhat = displacement / distance

    coupling_hartree = (
        float(np.dot(mu_a, mu_b))
        - 3.0
        * float(np.dot(mu_a, rhat))
        * float(np.dot(mu_b, rhat))
    ) / (distance**3)

    return coupling_hartree * HARTREE_TO_MEV


def charge_coupling_hartree(
    charges_a: np.ndarray,
    coordinates_a: np.ndarray,
    charges_b: np.ndarray,
    coordinates_b: np.ndarray,
) -> float:
    displacement = (
        coordinates_a[:, None, :]
        - coordinates_b[None, :, :]
    )
    distances = np.linalg.norm(displacement, axis=2)

    if np.any(distances <= 0.0):
        raise RuntimeError("Zero interatomic distance in charge coupling.")

    return float(
        np.sum(
            charges_a[:, None]
            * charges_b[None, :]
            / distances
        )
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    density_manifest = read_density_manifest()
    saved_point = read_frame000_point_couplings()

    site_data: dict[str, dict[str, object]] = {}
    charge_rows: list[dict[str, object]] = []
    validation_rows: list[dict[str, object]] = []

    log("Day019 bright transition-density atomic-charge benchmark")
    log("Partition: nearest-atom Voronoi on normalized N80 cubes")
    log("Correction: minimum-norm atomic-charge adjustment")
    log("Constraints: exact Q=0 and exact ORCA transition dipole")

    for site in SITES:
        metadata = density_manifest[site]
        cube = parse_cube(Path(metadata["cube_path"]))

        atom_records = list(cube["atom_records"])
        atomic_numbers = np.array(
            [record[0] for record in atom_records],
            dtype=np.int64,
        )
        atom_coordinates = np.array(
            [[record[2], record[3], record[4]] for record in atom_records],
            dtype=np.float64,
        )

        carbon_mask = atomic_numbers == 6

        if int(np.count_nonzero(carbon_mask)) != 16:
            raise RuntimeError(
                f"{site}: expected 16 carbon atoms, found "
                f"{int(np.count_nonzero(carbon_mask))}."
            )

        carbon_centroid = np.mean(
            atom_coordinates[carbon_mask],
            axis=0,
        )

        raw_charges, cube_integral = partition_to_nearest_atoms(cube)
        partition_residual = float(np.sum(raw_charges) - cube_integral)

        if abs(partition_residual) > PARTITION_CHARGE_TOL:
            raise RuntimeError(
                f"{site}: Voronoi charge partition does not recover the "
                f"cube integral. Residual={partition_residual:.3e}"
            )

        target_mu = np.asarray(
            metadata["orca_mu_au"],
            dtype=np.float64,
        )

        raw_mu = dipole_from_charges(
            charges=raw_charges,
            atom_coordinates=atom_coordinates,
            origin=carbon_centroid,
        )

        corrected_charges, final_residual = constrained_charge_correction(
            raw_charges=raw_charges,
            atom_coordinates=atom_coordinates,
            carbon_centroid=carbon_centroid,
            target_mu=target_mu,
        )

        corrected_mu = dipole_from_charges(
            charges=corrected_charges,
            atom_coordinates=atom_coordinates,
            origin=carbon_centroid,
        )

        correction = corrected_charges - raw_charges

        maximum_constraint_residual = float(
            np.max(np.abs(final_residual))
        )

        if maximum_constraint_residual > CONSTRAINT_TOL:
            raise RuntimeError(
                f"{site}: corrected charge constraints failed. "
                f"Max residual={maximum_constraint_residual:.3e}"
            )

        raw_mu_error = float(
            np.linalg.norm(raw_mu - target_mu)
            / np.linalg.norm(target_mu)
        )
        corrected_mu_error = float(
            np.linalg.norm(corrected_mu - target_mu)
            / np.linalg.norm(target_mu)
        )

        site_data[site] = {
            "root": int(metadata["root"]),
            "atomic_numbers": atomic_numbers,
            "atom_coordinates_bohr": atom_coordinates,
            "carbon_centroid_bohr": carbon_centroid,
            "raw_charges": raw_charges,
            "corrected_charges": corrected_charges,
            "target_mu_au": target_mu,
            "raw_mu_au": raw_mu,
            "corrected_mu_au": corrected_mu,
        }

        for atom_index, (
            atomic_number,
            coordinate,
            raw_charge,
            corrected_charge,
        ) in enumerate(
            zip(
                atomic_numbers,
                atom_coordinates,
                raw_charges,
                corrected_charges,
            ),
            start=1,
        ):
            charge_rows.append(
                {
                    "frame": 0,
                    "site": site,
                    "bright_root": int(metadata["root"]),
                    "atom_index": atom_index,
                    "atomic_number": int(atomic_number),
                    "element": ATOMIC_SYMBOL.get(
                        int(atomic_number),
                        f"Z{int(atomic_number)}",
                    ),
                    "x_bohr": coordinate[0],
                    "y_bohr": coordinate[1],
                    "z_bohr": coordinate[2],
                    "raw_voronoi_transition_charge_e": raw_charge,
                    "charge_correction_e": (
                        corrected_charge - raw_charge
                    ),
                    "corrected_transition_charge_e": corrected_charge,
                }
            )

        validation_rows.append(
            {
                "frame": 0,
                "site": site,
                "bright_root": int(metadata["root"]),
                "n_atoms": len(atom_records),
                "cube_net_transition_charge_e": cube_integral,
                "partitioned_net_charge_e": float(
                    np.sum(raw_charges)
                ),
                "partition_recovery_residual_e": partition_residual,
                "raw_atomic_charge_net_e": float(
                    np.sum(raw_charges)
                ),
                "corrected_atomic_charge_net_e": float(
                    np.sum(corrected_charges)
                ),
                "target_mu_x_au": target_mu[0],
                "target_mu_y_au": target_mu[1],
                "target_mu_z_au": target_mu[2],
                "raw_atomic_mu_x_au": raw_mu[0],
                "raw_atomic_mu_y_au": raw_mu[1],
                "raw_atomic_mu_z_au": raw_mu[2],
                "corrected_atomic_mu_x_au": corrected_mu[0],
                "corrected_atomic_mu_y_au": corrected_mu[1],
                "corrected_atomic_mu_z_au": corrected_mu[2],
                "raw_atomic_dipole_relative_error": raw_mu_error,
                "corrected_atomic_dipole_relative_error": (
                    corrected_mu_error
                ),
                "charge_correction_rms_e": float(
                    np.sqrt(np.mean(correction * correction))
                ),
                "charge_correction_max_abs_e": float(
                    np.max(np.abs(correction))
                ),
                "raw_charge_max_abs_e": float(
                    np.max(np.abs(raw_charges))
                ),
                "corrected_charge_max_abs_e": float(
                    np.max(np.abs(corrected_charges))
                ),
                "maximum_constraint_residual": (
                    maximum_constraint_residual
                ),
                "numerical_validation_pass": (
                    abs(partition_residual)
                    <= PARTITION_CHARGE_TOL
                    and maximum_constraint_residual
                    <= CONSTRAINT_TOL
                ),
            }
        )

        log(
            f"[{site}] S{int(metadata['root'])}: "
            f"Qcube={cube_integral:+.3e}, "
            f"raw dipole error={raw_mu_error:.3%}, "
            f"corrected dipole error={corrected_mu_error:.3e}, "
            f"RMS dq={float(np.sqrt(np.mean(correction * correction))):.3e} e"
        )

    write_csv(CHARGES_CSV, charge_rows)
    write_csv(SITE_VALIDATION_CSV, validation_rows)

    coupling_rows: list[dict[str, object]] = []
    point_matrix = np.zeros((4, 4), dtype=np.float64)
    tdcac_matrix = np.zeros((4, 4), dtype=np.float64)
    delta_matrix = np.zeros((4, 4), dtype=np.float64)

    for site_a, site_b in combinations(SITES, 2):
        data_a = site_data[site_a]
        data_b = site_data[site_b]

        raw_hartree = charge_coupling_hartree(
            charges_a=np.asarray(
                data_a["raw_charges"],
                dtype=np.float64,
            ),
            coordinates_a=np.asarray(
                data_a["atom_coordinates_bohr"],
                dtype=np.float64,
            ),
            charges_b=np.asarray(
                data_b["raw_charges"],
                dtype=np.float64,
            ),
            coordinates_b=np.asarray(
                data_b["atom_coordinates_bohr"],
                dtype=np.float64,
            ),
        )

        corrected_hartree = charge_coupling_hartree(
            charges_a=np.asarray(
                data_a["corrected_charges"],
                dtype=np.float64,
            ),
            coordinates_a=np.asarray(
                data_a["atom_coordinates_bohr"],
                dtype=np.float64,
            ),
            charges_b=np.asarray(
                data_b["corrected_charges"],
                dtype=np.float64,
            ),
            coordinates_b=np.asarray(
                data_b["atom_coordinates_bohr"],
                dtype=np.float64,
            ),
        )

        point_recomputed_meV = point_dipole_coupling_meV(
            mu_a=np.asarray(data_a["target_mu_au"], dtype=np.float64),
            mu_b=np.asarray(data_b["target_mu_au"], dtype=np.float64),
            center_a=np.asarray(
                data_a["carbon_centroid_bohr"],
                dtype=np.float64,
            ),
            center_b=np.asarray(
                data_b["carbon_centroid_bohr"],
                dtype=np.float64,
            ),
        )

        point_saved_meV = saved_point[(site_a, site_b)]
        point_reproduction_error = (
            point_recomputed_meV - point_saved_meV
        )

        if abs(point_reproduction_error) > POINT_REPRODUCTION_TOL_MEV:
            raise RuntimeError(
                f"{site_a}-{site_b}: independently recomputed point "
                f"coupling differs from saved baseline by "
                f"{point_reproduction_error:.6e} meV."
            )

        raw_meV = raw_hartree * HARTREE_TO_MEV
        corrected_meV = corrected_hartree * HARTREE_TO_MEV
        delta_meV = corrected_meV - point_saved_meV

        ratio_to_point = (
            corrected_meV / point_saved_meV
            if abs(point_saved_meV) > 1.0e-15
            else float("nan")
        )
        relative_deviation = (
            delta_meV / point_saved_meV
            if abs(point_saved_meV) > 1.0e-15
            else float("nan")
        )

        center_a = np.asarray(
            data_a["carbon_centroid_bohr"],
            dtype=np.float64,
        )
        center_b = np.asarray(
            data_b["carbon_centroid_bohr"],
            dtype=np.float64,
        )
        centroid_distance_bohr = float(
            np.linalg.norm(center_b - center_a)
        )

        minimum_atom_distance_bohr = float(
            np.min(
                np.linalg.norm(
                    np.asarray(
                        data_a["atom_coordinates_bohr"],
                        dtype=np.float64,
                    )[:, None, :]
                    - np.asarray(
                        data_b["atom_coordinates_bohr"],
                        dtype=np.float64,
                    )[None, :, :],
                    axis=2,
                )
            )
        )

        coupling_rows.append(
            {
                "frame": 0,
                "site_a": site_a,
                "bright_root_a": int(data_a["root"]),
                "site_b": site_b,
                "bright_root_b": int(data_b["root"]),
                "centroid_distance_bohr": centroid_distance_bohr,
                "minimum_atom_distance_bohr": (
                    minimum_atom_distance_bohr
                ),
                "raw_voronoi_J_hartree": raw_hartree,
                "raw_voronoi_J_meV": raw_meV,
                "corrected_TDCAC_J_hartree": corrected_hartree,
                "corrected_TDCAC_J_eV": (
                    corrected_hartree * HARTREE_TO_EV
                ),
                "corrected_TDCAC_J_meV": corrected_meV,
                "corrected_TDCAC_J_cm-1": (
                    corrected_hartree * HARTREE_TO_CM1
                ),
                "saved_point_dipole_J_meV": point_saved_meV,
                "recomputed_point_dipole_J_meV": (
                    point_recomputed_meV
                ),
                "point_reproduction_error_meV": (
                    point_reproduction_error
                ),
                "TDCAC_minus_point_meV": delta_meV,
                "TDCAC_over_point_ratio": ratio_to_point,
                "relative_deviation_from_point": relative_deviation,
                "absolute_relative_deviation_from_point": abs(
                    relative_deviation
                ),
            }
        )

        index_a = SITE_INDEX[site_a]
        index_b = SITE_INDEX[site_b]

        point_matrix[index_a, index_b] = point_saved_meV
        point_matrix[index_b, index_a] = point_saved_meV

        tdcac_matrix[index_a, index_b] = corrected_meV
        tdcac_matrix[index_b, index_a] = corrected_meV

        delta_matrix[index_a, index_b] = delta_meV
        delta_matrix[index_b, index_a] = delta_meV

        log(
            f"[{site_a}-{site_b}] "
            f"point={point_saved_meV:+.6f} meV, "
            f"TDC-AC={corrected_meV:+.6f} meV, "
            f"delta={delta_meV:+.6f} meV, "
            f"ratio={ratio_to_point:+.6f}"
        )

    write_csv(COUPLINGS_CSV, coupling_rows)
    write_matrix_csv(POINT_MATRIX_CSV, SITES, point_matrix)
    write_matrix_csv(TDCAC_MATRIX_CSV, SITES, tdcac_matrix)
    write_matrix_csv(DELTA_MATRIX_CSV, SITES, delta_matrix)

    nearest_rows = sorted(
        coupling_rows,
        key=lambda row: float(row["centroid_distance_bohr"]),
    )[:3]
    distant_rows = sorted(
        coupling_rows,
        key=lambda row: float(row["centroid_distance_bohr"]),
    )[3:]

    maximum_point_reproduction_error = max(
        abs(float(row["point_reproduction_error_meV"]))
        for row in coupling_rows
    )
    maximum_abs_delta = max(
        abs(float(row["TDCAC_minus_point_meV"]))
        for row in coupling_rows
    )
    maximum_abs_relative_deviation = max(
        float(row["absolute_relative_deviation_from_point"])
        for row in coupling_rows
    )
    maximum_rms_correction = max(
        float(row["charge_correction_rms_e"])
        for row in validation_rows
    )
    maximum_charge_correction = max(
        float(row["charge_correction_max_abs_e"])
        for row in validation_rows
    )

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 bright TDC-derived atomic-charge coupling benchmark\n\n"
        )

        handle.write("## Method\n\n")
        handle.write(
            "- Input: four phase-aligned, sqrt(2)-normalized N80 bright "
            "transition-density cubes.\n"
        )
        handle.write(
            "- Each voxel charge was assigned to its nearest atom "
            "(atom-centered Voronoi partition).\n"
        )
        handle.write(
            "- A minimum-norm correction imposed exactly zero total "
            "transition charge and the independently printed ORCA "
            "transition dipole.\n"
        )
        handle.write(
            "- Coulomb couplings were evaluated between the resulting "
            "26 atom-centered transition charges on each monomer.\n"
        )
        handle.write(
            "- This is a transition-density-derived atomic-charge "
            "(TDC-AC) model, not the exact double integral over the two "
            "continuous three-dimensional densities.\n\n"
        )

        handle.write("## Numerical validation\n\n")
        handle.write(
            f"- Sites partitioned and constrained: "
            f"{sum(bool(row['numerical_validation_pass']) for row in validation_rows)}/4\n"
        )
        handle.write(
            f"- Maximum independent point-dipole reproduction error: "
            f"{maximum_point_reproduction_error:.6e} meV\n"
        )
        handle.write(
            f"- Maximum RMS atomic-charge correction: "
            f"{maximum_rms_correction:.6e} e\n"
        )
        handle.write(
            f"- Maximum absolute atomic-charge correction: "
            f"{maximum_charge_correction:.6e} e\n\n"
        )

        handle.write("## Site charge validation\n\n")
        handle.write(
            "| Site | Root | Cube Q | Raw dipole error | "
            "Corrected dipole error | RMS dq (e) | Max |dq| (e) |\n"
        )
        handle.write("|---|---:|---:|---:|---:|---:|---:|\n")

        for row in validation_rows:
            handle.write(
                f"| {row['site']} "
                f"| S{row['bright_root']} "
                f"| {float(row['cube_net_transition_charge_e']):+.3e} "
                f"| {float(row['raw_atomic_dipole_relative_error']):.4%} "
                f"| {float(row['corrected_atomic_dipole_relative_error']):.3e} "
                f"| {float(row['charge_correction_rms_e']):.3e} "
                f"| {float(row['charge_correction_max_abs_e']):.3e} |\n"
            )

        handle.write("\n## Coupling comparison\n\n")
        handle.write(
            "| Pair | Point dipole (meV) | Raw Voronoi (meV) | "
            "Corrected TDC-AC (meV) | Delta (meV) | "
            "TDC-AC / point |\n"
        )
        handle.write("|---|---:|---:|---:|---:|---:|\n")

        for row in coupling_rows:
            handle.write(
                f"| {row['site_a']}-{row['site_b']} "
                f"| {float(row['saved_point_dipole_J_meV']):+.6f} "
                f"| {float(row['raw_voronoi_J_meV']):+.6f} "
                f"| {float(row['corrected_TDCAC_J_meV']):+.6f} "
                f"| {float(row['TDCAC_minus_point_meV']):+.6f} "
                f"| {float(row['TDCAC_over_point_ratio']):+.6f} |\n"
            )

        handle.write("\n## Aggregate deviations\n\n")
        handle.write(
            f"- Maximum |TDC-AC - point|: "
            f"{maximum_abs_delta:.6f} meV\n"
        )
        handle.write(
            f"- Maximum absolute relative deviation: "
            f"{maximum_abs_relative_deviation:.4%}\n"
        )
        handle.write(
            "- Mean absolute relative deviation, three nearest pairs: "
            f"{np.mean([float(row['absolute_relative_deviation_from_point']) for row in nearest_rows]):.4%}\n"
        )
        handle.write(
            "- Mean absolute relative deviation, three distant pairs: "
            f"{np.mean([float(row['absolute_relative_deviation_from_point']) for row in distant_rows]):.4%}\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "The corrected TDC-AC model reproduces the exact transition "
            "charge and dipole by construction, so differences from the "
            "point-dipole result originate from the atom-resolved finite "
            "extent of the transition density. Agreement at the more distant "
            "pairs supports the far-field limit. Deviations at neighboring "
            "pairs quantify finite-size corrections within this atom-centered "
            "discretization. Because the Voronoi partition is not unique, "
            "these values should be treated as a finite-size benchmark rather "
            "than as an exact continuous-density Coulomb integral. No "
            "dielectric screening is included.\n"
        )

    log("")
    log("Day019 bright TDC-derived atomic-charge benchmark completed.")
    log("Validated sites: 4/4")
    log("Coupling pairs: 6/6")
    log(
        f"Maximum point reproduction error: "
        f"{maximum_point_reproduction_error:.3e} meV"
    )
    log(
        f"Maximum |TDC-AC - point|: "
        f"{maximum_abs_delta:.6f} meV"
    )
    log(
        f"Maximum relative deviation: "
        f"{maximum_abs_relative_deviation:.4%}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
