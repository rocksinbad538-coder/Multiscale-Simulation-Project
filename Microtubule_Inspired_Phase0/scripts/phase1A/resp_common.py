#!/usr/bin/env python3
"""
resp_common.py

Shared utilities for the RESP pipeline.

Version: 0.2 (ORCA ESP parsing foundation)

Purpose
-------
Centralize all reusable validation logic for the RESP pipeline so that
all gates share identical provenance, hashing, JSON handling and report
generation.

Current utilities
-----------------
- sha256()
- load_json()
- save_json()
- require_file()
- require_authorization()
- require_decision()
- validate_execution_binding()
- utc_now()
- read_orca_vpot()
- validate_orca_vpot()
- read_orca_pc_chelpg()
- validate_orca_pc_chelpg()
- coulomb_potential_at_point()
- reconstruct_coulomb_potential()
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Iterable, Sequence


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def require_file(path: Path) -> None:
    if (not path.exists()) or path.stat().st_size == 0:
        raise RuntimeError(f"Required file missing: {path}")


def load_json(path: Path) -> dict:
    require_file(path)
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    require_file(path)

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def require_decision(report: dict, expected: str) -> None:
    decision = report.get("decision")

    if decision != expected:
        raise RuntimeError(
            f"Decision mismatch.\n"
            f"Expected: {expected}\n"
            f"Observed: {decision}"
        )


def require_authorization(
    report: dict,
    authorization_name: str,
) -> None:

    authorizations = report.get("authorizations", {})

    if authorizations.get(authorization_name) is not True:
        raise RuntimeError(
            f"Authorization '{authorization_name}' not granted."
        )


def validate_execution_binding(
    previous_report: dict,
    execution_directory: str,
) -> None:

    observed = previous_report.get(
        "source_execution_directory",
        previous_report.get("execution_directory"),
    )

    if observed != execution_directory:
        raise RuntimeError(
            "Execution binding failed.\n"
            f"Expected: {execution_directory}\n"
            f"Observed: {observed}"
        )


def build_report(
    *,
    decision: str,
    execution_directory: str,
    authorizations: dict,
    sha256_items: dict | None = None,
    summary: dict | None = None,
    metadata: dict | None = None,
) -> dict:

    report = {
        "generated_utc": utc_now(),
        "decision": decision,
        "execution_directory": execution_directory,
        "authorizations": authorizations,
    }

    if sha256_items:
        report["sha256"] = sha256_items

    if summary:
        report["summary"] = summary

    if metadata:
        report["metadata"] = metadata

    return report



BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM


@dataclass(frozen=True)
class OrcaVpotData:
    """
    Parsed ORCA CHELPG electrostatic-potential dataset.

    Confirmed Day038 schema
    -----------------------
    Header:
        atom_count grid_point_count

    Atomic block:
        x_bohr y_bohr z_bohr

    Grid block:
        potential_au x_bohr y_bohr z_bohr
    """

    source_path: Path
    source_sha256: str
    atom_count: int
    grid_point_count: int
    atom_coordinates_bohr: tuple[
        tuple[float, float, float],
        ...,
    ]
    grid_potential_au: tuple[float, ...]
    grid_coordinates_bohr: tuple[
        tuple[float, float, float],
        ...,
    ]


@dataclass(frozen=True)
class OrcaChelpgPointCharges:
    """
    Parsed ORCA `.pc_chelpg` point-charge artifact.

    Confirmed Day038 schema
    -----------------------
    First line:
        atom_count

    Second line:
        comment

    Charge rows:
        Q charge_e x_angstrom y_angstrom z_angstrom
    """

    source_path: Path
    source_sha256: str
    atom_count: int
    comment: str
    charges_e: tuple[float, ...]
    coordinates_angstrom: tuple[
        tuple[float, float, float],
        ...,
    ]
    coordinates_bohr: tuple[
        tuple[float, float, float],
        ...,
    ]


def parse_scientific_float(token: str) -> float:
    """
    Parse standard or Fortran-style scientific notation.
    """

    value = float(
        token.replace("D", "E").replace("d", "e")
    )

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite numerical value encountered: {token!r}"
        )

    return value


def convert_bohr_to_angstrom(
    coordinates_bohr: Sequence[float],
) -> tuple[float, ...]:
    """
    Convert one coordinate vector from bohr to angstrom.
    """

    return tuple(
        float(value) * BOHR_TO_ANGSTROM
        for value in coordinates_bohr
    )


def convert_angstrom_to_bohr(
    coordinates_angstrom: Sequence[float],
) -> tuple[float, ...]:
    """
    Convert one coordinate vector from angstrom to bohr.
    """

    return tuple(
        float(value) * ANGSTROM_TO_BOHR
        for value in coordinates_angstrom
    )


def read_orca_vpot(
    path: Path,
    *,
    expected_sha256: str | None = None,
) -> OrcaVpotData:
    """
    Parse an ORCA CHELPG `.vpot` file.

    The function validates:

    - source existence and non-zero size;
    - optional source SHA256;
    - two-integer header;
    - exact line-count partition;
    - three numerical fields per atom row;
    - four numerical fields per grid row;
    - finite coordinates and potential values.

    No unit conversion is applied. Native `.vpot` coordinates are
    returned in bohr and the potential column in atomic units.
    """

    require_file(path)

    observed_sha256 = sha256(path)

    if (
        expected_sha256 is not None
        and observed_sha256 != expected_sha256
    ):
        raise RuntimeError(
            "ORCA VPOT SHA256 mismatch.\n"
            f"Expected: {expected_sha256}\n"
            f"Observed: {observed_sha256}\n"
            f"Path: {path}"
        )

    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    if not lines:
        raise RuntimeError(
            f"ORCA VPOT file is empty: {path}"
        )

    header = lines[0].split()

    if len(header) != 2:
        raise RuntimeError(
            "Invalid ORCA VPOT header.\n"
            "Expected two integer fields.\n"
            f"Observed: {lines[0]!r}"
        )

    try:
        atom_count = int(header[0])
        grid_point_count = int(header[1])
    except ValueError as exc:
        raise RuntimeError(
            f"Non-integer ORCA VPOT header: {header!r}"
        ) from exc

    if atom_count <= 0:
        raise RuntimeError(
            f"Invalid ORCA VPOT atom count: {atom_count}"
        )

    if grid_point_count <= 0:
        raise RuntimeError(
            "Invalid ORCA VPOT grid-point count: "
            f"{grid_point_count}"
        )

    expected_line_count = (
        1 + atom_count + grid_point_count
    )

    if len(lines) != expected_line_count:
        raise RuntimeError(
            "ORCA VPOT line-count mismatch.\n"
            f"Expected: {expected_line_count}\n"
            f"Observed: {len(lines)}"
        )

    atom_lines = lines[1:1 + atom_count]
    grid_lines = lines[1 + atom_count:]

    atom_coordinates: list[
        tuple[float, float, float]
    ] = []

    for relative_index, line in enumerate(
        atom_lines,
        start=1,
    ):
        tokens = line.split()

        if len(tokens) != 3:
            raise RuntimeError(
                "Malformed ORCA VPOT atom row.\n"
                f"Relative atom row: {relative_index}\n"
                f"Expected fields: 3\n"
                f"Observed fields: {len(tokens)}\n"
                f"Raw line: {line!r}"
            )

        atom_coordinates.append(
            tuple(
                parse_scientific_float(token)
                for token in tokens
            )
        )

    grid_potential: list[float] = []
    grid_coordinates: list[
        tuple[float, float, float]
    ] = []

    for relative_index, line in enumerate(
        grid_lines,
        start=1,
    ):
        tokens = line.split()

        if len(tokens) != 4:
            raise RuntimeError(
                "Malformed ORCA VPOT grid row.\n"
                f"Relative grid row: {relative_index}\n"
                f"Expected fields: 4\n"
                f"Observed fields: {len(tokens)}\n"
                f"Raw line: {line!r}"
            )

        values = tuple(
            parse_scientific_float(token)
            for token in tokens
        )

        potential, x, y, z = values

        grid_potential.append(potential)
        grid_coordinates.append((x, y, z))

    return OrcaVpotData(
        source_path=path.resolve(),
        source_sha256=observed_sha256,
        atom_count=atom_count,
        grid_point_count=grid_point_count,
        atom_coordinates_bohr=tuple(atom_coordinates),
        grid_potential_au=tuple(grid_potential),
        grid_coordinates_bohr=tuple(grid_coordinates),
    )


def validate_orca_vpot(
    dataset: OrcaVpotData,
    *,
    expected_atom_count: int | None = None,
    expected_grid_point_count: int | None = None,
) -> dict:
    """
    Validate an already parsed ORCA VPOT dataset.

    Returns a machine-readable validation summary. Raises an exception
    if any required invariant fails.
    """

    if dataset.atom_count != len(
        dataset.atom_coordinates_bohr
    ):
        raise RuntimeError(
            "ORCA VPOT atom-count invariant failed"
        )

    if dataset.grid_point_count != len(
        dataset.grid_potential_au
    ):
        raise RuntimeError(
            "ORCA VPOT potential-count invariant failed"
        )

    if dataset.grid_point_count != len(
        dataset.grid_coordinates_bohr
    ):
        raise RuntimeError(
            "ORCA VPOT grid-coordinate-count invariant failed"
        )

    if (
        expected_atom_count is not None
        and dataset.atom_count != expected_atom_count
    ):
        raise RuntimeError(
            "Unexpected ORCA VPOT atom count.\n"
            f"Expected: {expected_atom_count}\n"
            f"Observed: {dataset.atom_count}"
        )

    if (
        expected_grid_point_count is not None
        and dataset.grid_point_count
        != expected_grid_point_count
    ):
        raise RuntimeError(
            "Unexpected ORCA VPOT grid-point count.\n"
            f"Expected: {expected_grid_point_count}\n"
            f"Observed: {dataset.grid_point_count}"
        )

    atom_values = [
        value
        for row in dataset.atom_coordinates_bohr
        for value in row
    ]

    grid_coordinate_values = [
        value
        for row in dataset.grid_coordinates_bohr
        for value in row
    ]

    all_values = (
        atom_values
        + grid_coordinate_values
        + list(dataset.grid_potential_au)
    )

    if not all(math.isfinite(value) for value in all_values):
        raise RuntimeError(
            "ORCA VPOT contains non-finite values"
        )

    return {
        "source_path": str(dataset.source_path),
        "source_sha256": dataset.source_sha256,
        "atom_count": dataset.atom_count,
        "grid_point_count": dataset.grid_point_count,
        "coordinate_unit": "bohr",
        "potential_unit": "atomic_unit_Eh_per_e",
        "atom_coordinate_min_bohr": min(atom_values),
        "atom_coordinate_max_bohr": max(atom_values),
        "grid_coordinate_min_bohr": min(
            grid_coordinate_values
        ),
        "grid_coordinate_max_bohr": max(
            grid_coordinate_values
        ),
        "potential_min_au": min(
            dataset.grid_potential_au
        ),
        "potential_max_au": max(
            dataset.grid_potential_au
        ),
        "potential_mean_au": (
            sum(dataset.grid_potential_au)
            / dataset.grid_point_count
        ),
        "finite_value_gate": True,
        "count_invariant_gate": True,
    }


def read_orca_pc_chelpg(
    path: Path,
    *,
    expected_sha256: str | None = None,
) -> OrcaChelpgPointCharges:
    """
    Parse an ORCA `.pc_chelpg` XYZ-style point-charge file.
    """

    require_file(path)

    observed_sha256 = sha256(path)

    if (
        expected_sha256 is not None
        and observed_sha256 != expected_sha256
    ):
        raise RuntimeError(
            "ORCA .pc_chelpg SHA256 mismatch.\n"
            f"Expected: {expected_sha256}\n"
            f"Observed: {observed_sha256}\n"
            f"Path: {path}"
        )

    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"Incomplete ORCA .pc_chelpg file: {path}"
        )

    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise RuntimeError(
            "The first .pc_chelpg line is not an integer "
            "atom count"
        ) from exc

    comment = lines[1]

    charges: list[float] = []
    coordinates_angstrom: list[
        tuple[float, float, float]
    ] = []

    for line_number, line in enumerate(
        lines[2:],
        start=3,
    ):
        if not line.strip():
            continue

        tokens = line.split()

        if len(tokens) != 5:
            raise RuntimeError(
                "Malformed ORCA .pc_chelpg row.\n"
                f"Line: {line_number}\n"
                f"Expected fields: 5\n"
                f"Observed fields: {len(tokens)}\n"
                f"Raw line: {line!r}"
            )

        if tokens[0].upper() != "Q":
            raise RuntimeError(
                "Unexpected ORCA .pc_chelpg row marker.\n"
                f"Line: {line_number}\n"
                f"Expected: Q\n"
                f"Observed: {tokens[0]!r}"
            )

        charge = parse_scientific_float(tokens[1])

        coordinates = tuple(
            parse_scientific_float(token)
            for token in tokens[2:5]
        )

        charges.append(charge)
        coordinates_angstrom.append(coordinates)

    if len(charges) != atom_count:
        raise RuntimeError(
            "ORCA .pc_chelpg count mismatch.\n"
            f"Declared: {atom_count}\n"
            f"Parsed: {len(charges)}"
        )

    coordinates_bohr = tuple(
        convert_angstrom_to_bohr(row)
        for row in coordinates_angstrom
    )

    return OrcaChelpgPointCharges(
        source_path=path.resolve(),
        source_sha256=observed_sha256,
        atom_count=atom_count,
        comment=comment,
        charges_e=tuple(charges),
        coordinates_angstrom=tuple(
            coordinates_angstrom
        ),
        coordinates_bohr=coordinates_bohr,
    )


def validate_orca_pc_chelpg(
    dataset: OrcaChelpgPointCharges,
    *,
    expected_atom_count: int | None = None,
    expected_total_charge_e: float | None = None,
    total_charge_tolerance_e: float = 1.0e-5,
) -> dict:
    """
    Validate an already parsed ORCA `.pc_chelpg` dataset.
    """

    if dataset.atom_count != len(dataset.charges_e):
        raise RuntimeError(
            "ORCA .pc_chelpg charge-count invariant failed"
        )

    if dataset.atom_count != len(
        dataset.coordinates_angstrom
    ):
        raise RuntimeError(
            "ORCA .pc_chelpg angstrom-coordinate-count "
            "invariant failed"
        )

    if dataset.atom_count != len(
        dataset.coordinates_bohr
    ):
        raise RuntimeError(
            "ORCA .pc_chelpg bohr-coordinate-count "
            "invariant failed"
        )

    if (
        expected_atom_count is not None
        and dataset.atom_count != expected_atom_count
    ):
        raise RuntimeError(
            "Unexpected ORCA .pc_chelpg atom count.\n"
            f"Expected: {expected_atom_count}\n"
            f"Observed: {dataset.atom_count}"
        )

    total_charge = sum(dataset.charges_e)

    if (
        expected_total_charge_e is not None
        and abs(total_charge - expected_total_charge_e)
        > total_charge_tolerance_e
    ):
        raise RuntimeError(
            "ORCA .pc_chelpg total-charge gate failed.\n"
            f"Expected: {expected_total_charge_e}\n"
            f"Observed: {total_charge}\n"
            f"Tolerance: {total_charge_tolerance_e}"
        )

    all_values = (
        list(dataset.charges_e)
        + [
            value
            for row in dataset.coordinates_angstrom
            for value in row
        ]
        + [
            value
            for row in dataset.coordinates_bohr
            for value in row
        ]
    )

    if not all(math.isfinite(value) for value in all_values):
        raise RuntimeError(
            "ORCA .pc_chelpg contains non-finite values"
        )

    return {
        "source_path": str(dataset.source_path),
        "source_sha256": dataset.source_sha256,
        "atom_count": dataset.atom_count,
        "charge_unit": "elementary_charge",
        "coordinate_unit_native": "angstrom",
        "total_charge_e": total_charge,
        "minimum_charge_e": min(dataset.charges_e),
        "maximum_charge_e": max(dataset.charges_e),
        "finite_value_gate": True,
        "count_invariant_gate": True,
        "total_charge_gate": True,
    }


def coordinate_difference_metrics(
    reference: Sequence[Sequence[float]],
    candidate: Sequence[Sequence[float]],
) -> dict:
    """
    Calculate component-wise and vector coordinate differences.

    Both inputs must use the same units and preserve atom ordering.
    """

    if len(reference) != len(candidate):
        raise RuntimeError(
            "Coordinate sequence lengths differ.\n"
            f"Reference: {len(reference)}\n"
            f"Candidate: {len(candidate)}"
        )

    differences = [
        tuple(
            candidate_value - reference_value
            for reference_value, candidate_value
            in zip(reference_row, candidate_row)
        )
        for reference_row, candidate_row
        in zip(reference, candidate)
    ]

    if any(len(row) != 3 for row in differences):
        raise RuntimeError(
            "Coordinate rows must contain three components"
        )

    absolute_components = [
        abs(value)
        for row in differences
        for value in row
    ]

    displacement_norms = [
        math.sqrt(
            sum(value * value for value in row)
        )
        for row in differences
    ]

    squared_components = [
        value * value
        for row in differences
        for value in row
    ]

    return {
        "maximum_absolute_component": max(
            absolute_components
        ),
        "mean_absolute_component": (
            sum(absolute_components)
            / len(absolute_components)
        ),
        "maximum_vector_displacement": max(
            displacement_norms
        ),
        "mean_vector_displacement": (
            sum(displacement_norms)
            / len(displacement_norms)
        ),
        "component_rmsd": math.sqrt(
            sum(squared_components)
            / len(squared_components)
        ),
        "vector_rmsd": math.sqrt(
            sum(value * value for value in displacement_norms)
            / len(displacement_norms)
        ),
    }


def coulomb_potential_at_point(
    point_bohr: Sequence[float],
    charges_e: Sequence[float],
    charge_coordinates_bohr: Sequence[
        Sequence[float]
    ],
    *,
    minimum_distance_bohr: float = 1.0e-12,
) -> tuple[float, float, int]:
    """
    Reconstruct the Coulomb potential at one point.

    With charges in elementary-charge units and distances in bohr,
    the returned potential is in atomic units (Eh/e).
    """

    if len(point_bohr) != 3:
        raise RuntimeError(
            "A Coulomb evaluation point must have "
            "three coordinates"
        )

    if len(charges_e) != len(charge_coordinates_bohr):
        raise RuntimeError(
            "Charge and coordinate counts differ"
        )

    px, py, pz = (
        float(point_bohr[0]),
        float(point_bohr[1]),
        float(point_bohr[2]),
    )

    potential = 0.0
    nearest_distance = math.inf
    nearest_index = -1

    for atom_index, (
        charge,
        coordinates,
    ) in enumerate(
        zip(charges_e, charge_coordinates_bohr)
    ):
        if len(coordinates) != 3:
            raise RuntimeError(
                "Charge-coordinate rows must contain "
                "three values"
            )

        x, y, z = (
            float(coordinates[0]),
            float(coordinates[1]),
            float(coordinates[2]),
        )

        dx = px - x
        dy = py - y
        dz = pz - z

        distance = math.sqrt(
            dx * dx + dy * dy + dz * dz
        )

        if distance < nearest_distance:
            nearest_distance = distance
            nearest_index = atom_index

        if distance <= minimum_distance_bohr:
            raise RuntimeError(
                "ESP grid point is too close to a point charge.\n"
                f"Atom index: {atom_index}\n"
                f"Distance (bohr): {distance}"
            )

        potential += float(charge) / distance

    if not math.isfinite(potential):
        raise RuntimeError(
            "Non-finite reconstructed Coulomb potential"
        )

    return potential, nearest_distance, nearest_index


def reconstruct_coulomb_potential(
    grid_coordinates_bohr: Iterable[
        Sequence[float]
    ],
    charges_e: Sequence[float],
    charge_coordinates_bohr: Sequence[
        Sequence[float]
    ],
) -> tuple[
    tuple[float, ...],
    tuple[float, ...],
    tuple[int, ...],
]:
    """
    Reconstruct Coulomb potentials for an ESP grid.

    Returns:

    - potential values in atomic units;
    - nearest-charge distances in bohr;
    - nearest-charge indices.
    """

    potential_values: list[float] = []
    nearest_distances: list[float] = []
    nearest_indices: list[int] = []

    for point in grid_coordinates_bohr:
        potential, distance, index = (
            coulomb_potential_at_point(
                point,
                charges_e,
                charge_coordinates_bohr,
            )
        )

        potential_values.append(potential)
        nearest_distances.append(distance)
        nearest_indices.append(index)

    return (
        tuple(potential_values),
        tuple(nearest_distances),
        tuple(nearest_indices),
    )


def print_gate_banner(title: str) -> None:
    print("=" * 100)
    print(title)
    print("=" * 100)

# ======================================================================================
# ORCA VPOT -> Amber ESP writer
# Day038 / D038-C1
# ======================================================================================

from pathlib import Path


def write_amber_esp(
    output_path: Path,
    atom_xyz_bohr,
    esp_xyz_bohr,
    esp_values_au,
) -> None:
    """
    Write an Amber-compatible ESP file from an ORCA VPOT dataset.

    Parameters
    ----------
    output_path
        Destination file.

    atom_xyz_bohr
        ndarray (Natoms,3)

    esp_xyz_bohr
        ndarray (Npoints,3)

    esp_values_au
        ndarray (Npoints,)
    """

    natoms = len(atom_xyz_bohr)
    npoints = len(esp_xyz_bohr)

    with output_path.open("w") as f:

        f.write(f"{natoms:5d}{npoints:5d}\n")

        for xyz in atom_xyz_bohr:
            f.write(
                " " * 17
                + f"{xyz[0]:16.7E}"
                + f"{xyz[1]:16.7E}"
                + f"{xyz[2]:16.7E}\n"
            )

        for value, xyz in zip(esp_values_au, esp_xyz_bohr):
            f.write(
                " "
                + f"{value:16.7E}"
                + f"{xyz[0]:16.7E}"
                + f"{xyz[1]:16.7E}"
                + f"{xyz[2]:16.7E}\n"
            )


def read_amber_esp(path: Path):
    """
    Read an Amber RESP fixed-width electrostatic-potential file.

    Confirmed Amber RESP records
    ----------------------------
    Header:
        (2I5)

    Atomic coordinates:
        (17X,3E16.7)

    ESP-grid records:
        (1X,4E16.7)

    Returns
    -------
    dict
        Counts and NumPy arrays in native RESP units:
        coordinates in bohr and potential in atomic units.
    """

    require_file(path)

    import numpy as np

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:

        header_raw = handle.readline().rstrip("\n")

        if len(header_raw) < 10:
            raise RuntimeError(
                "Amber ESP header is shorter than the "
                "required (2I5) record"
            )

        try:
            natoms = int(header_raw[0:5])
            npoints = int(header_raw[5:10])
        except ValueError as exc:
            raise RuntimeError(
                "Amber ESP header does not satisfy (2I5).\n"
                f"Observed: {header_raw!r}"
            ) from exc

        if natoms <= 0 or npoints <= 0:
            raise RuntimeError(
                "Amber ESP contains invalid counts.\n"
                f"natoms={natoms}, npoints={npoints}"
            )

        atom_rows = []

        for atom_index in range(natoms):
            line = handle.readline().rstrip("\n")

            if len(line) < 65:
                raise RuntimeError(
                    "Amber ESP atomic record is shorter "
                    "than (17X,3E16.7).\n"
                    f"Atom index: {atom_index}\n"
                    f"Observed length: {len(line)}\n"
                    f"Raw line: {line!r}"
                )

            try:
                row = [
                    parse_scientific_float(line[17:33]),
                    parse_scientific_float(line[33:49]),
                    parse_scientific_float(line[49:65]),
                ]
            except ValueError as exc:
                raise RuntimeError(
                    "Failed to parse Amber ESP atomic "
                    f"record {atom_index}"
                ) from exc

            atom_rows.append(row)

        grid_rows = []

        for point_index in range(npoints):
            line = handle.readline().rstrip("\n")

            if len(line) < 65:
                raise RuntimeError(
                    "Amber ESP grid record is shorter "
                    "than (1X,4E16.7).\n"
                    f"Grid index: {point_index}\n"
                    f"Observed length: {len(line)}\n"
                    f"Raw line: {line!r}"
                )

            try:
                row = [
                    parse_scientific_float(line[1:17]),
                    parse_scientific_float(line[17:33]),
                    parse_scientific_float(line[33:49]),
                    parse_scientific_float(line[49:65]),
                ]
            except ValueError as exc:
                raise RuntimeError(
                    "Failed to parse Amber ESP grid "
                    f"record {point_index}"
                ) from exc

            grid_rows.append(row)

        trailing_nonempty = [
            line
            for line in handle
            if line.strip()
        ]

        if trailing_nonempty:
            raise RuntimeError(
                "Amber ESP contains unexpected non-empty "
                "records after the declared dataset"
            )

    atoms = np.asarray(
        atom_rows,
        dtype=float,
    )

    grid = np.asarray(
        grid_rows,
        dtype=float,
    )

    if atoms.shape != (natoms, 3):
        raise RuntimeError(
            "Amber ESP atomic-array shape mismatch.\n"
            f"Expected: {(natoms, 3)}\n"
            f"Observed: {atoms.shape}"
        )

    if grid.shape != (npoints, 4):
        raise RuntimeError(
            "Amber ESP grid-array shape mismatch.\n"
            f"Expected: {(npoints, 4)}\n"
            f"Observed: {grid.shape}"
        )

    return {
        "natoms": natoms,
        "npoints": npoints,
        "atom_xyz_bohr": atoms,
        "esp_values_au": grid[:, 0],
        "esp_xyz_bohr": grid[:, 1:4],
        "format": {
            "header": "2I5",
            "atom_records": "17X,3E16.7",
            "grid_records": "1X,4E16.7",
        },
    }
