#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

SOURCE_SELECTION = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_full_two_state_manifold_analysis/"
    "SOURCE_OUTPUT_SELECTION_DAY019.csv"
)

FRAME_SUMMARY = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_full_two_state_manifold_analysis/"
    "two_state_frame_summary.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_dipole_geometry_audit"
)

DIPOLES_CSV = OUTPUT_ROOT / "transition_dipole_observations.csv"
GEOMETRY_FRAMES_CSV = OUTPUT_ROOT / "geometry_frame_audit.csv"
GEOMETRY_CONSISTENCY_CSV = OUTPUT_ROOT / "site_geometry_consistency.csv"
SITE_PAIRS_CSV = OUTPUT_ROOT / "site_pair_geometry.csv"
REPORT_MD = OUTPUT_ROOT / "TRANSITION_DIPOLE_GEOMETRY_AUDIT_DAY019.md"

NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
SCF_MARKER = "SCF CONVERGED"
TDDFT_MARKER = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"

STATE_ENERGY_RE = re.compile(
    r"STATE\s+(?P<root>\d+)\s*:\s*"
    r"E\s*=\s*(?P<energy_au>[-+0-9.Ee]+)\s+au\s+"
    r"(?P<energy_eV>[-+0-9.Ee]+)\s+eV",
    re.IGNORECASE,
)

TARGET_STATE_RE = re.compile(r"^(?P<root>\d+)-")

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
FRAMES = tuple(range(21))
FRAME_SPACING_PS = 5.0

EXPECTED_BRIGHT_ROOT = {
    "PYR2": 2,
    "PYR3": 2,
    "PYR4": 2,
    "PYR5": 1,
}

CARBON_SYMBOL = "C"
# ORCA prints D2 and Cartesian dipole components with finite decimal
# precision. A 5e-5 absolute tolerance accommodates the maximum
# reconstruction error expected from those rounded table values while
# remaining negligible relative to bright-state D2 values (~4-7 au^2).
D2_ABS_TOL = 5.0e-5
FOSC_ABS_TOL = 2.0e-5
GEOMETRY_TOL_A = 1.0e-8


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


def read_selected_outputs() -> dict[tuple[int, str], Path]:
    if not SOURCE_SELECTION.is_file():
        raise SystemExit(
            f"Missing selected-output table: {SOURCE_SELECTION}"
        )

    selected: dict[tuple[int, str], Path] = {}

    with SOURCE_SELECTION.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        for row in csv.DictReader(handle):
            if not as_bool(row["selected"]):
                continue

            key = (int(row["frame"]), row["site"])
            path = (PROJECT_ROOT / row["path"]).resolve()

            if key in selected:
                raise RuntimeError(
                    f"Duplicate selected output for {key}: "
                    f"{selected[key]} and {path}"
                )

            selected[key] = path

    expected = {
        (frame, site)
        for frame in FRAMES
        for site in SITES
    }

    if set(selected) != expected:
        missing = sorted(expected - set(selected))
        extra = sorted(set(selected) - expected)
        raise RuntimeError(
            f"Selected-output coverage mismatch. "
            f"Missing={missing}; extra={extra}"
        )

    return selected


def read_frame_mapping() -> dict[tuple[int, str], dict[str, str]]:
    if not FRAME_SUMMARY.is_file():
        raise SystemExit(f"Missing frame summary: {FRAME_SUMMARY}")

    mapping: dict[tuple[int, str], dict[str, str]] = {}

    with FRAME_SUMMARY.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        for row in csv.DictReader(handle):
            key = (int(row["frame"]), row["site"])

            if key in mapping:
                raise RuntimeError(f"Duplicate frame mapping: {key}")

            mapping[key] = row

    if len(mapping) != 84:
        raise RuntimeError(
            f"Expected 84 frame mappings, found {len(mapping)}."
        )

    return mapping


def float_tokens_after_target(line: str) -> tuple[int, list[float]] | None:
    tokens = line.split()

    if "->" not in tokens:
        return None

    arrow_index = tokens.index("->")

    if arrow_index + 1 >= len(tokens):
        return None

    target = tokens[arrow_index + 1]
    match = TARGET_STATE_RE.match(target)

    if not match:
        return None

    root = int(match.group("root"))
    values: list[float] = []

    for token in tokens[arrow_index + 2 :]:
        try:
            values.append(float(token))
        except ValueError:
            break

    if len(values) < 8:
        return None

    return root, values[:8]


def parse_transition_dipoles(
    text: str,
) -> dict[int, dict[str, float]]:
    lines = text.splitlines()
    in_length_table = False
    parsed: dict[int, dict[str, float]] = {}

    for line in lines:
        if (
            "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC "
            "DIPOLE MOMENTS" in line
        ):
            in_length_table = True
            continue

        if (
            in_length_table
            and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY"
            in line
        ):
            break

        if not in_length_table:
            continue

        result = float_tokens_after_target(line)

        if result is None:
            continue

        root, values = result

        if root not in (1, 2):
            continue

        (
            energy_eV,
            energy_cm1,
            wavelength_nm,
            fosc,
            d2,
            dx,
            dy,
            dz,
        ) = values

        parsed[root] = {
            "table_energy_eV": energy_eV,
            "energy_cm1": energy_cm1,
            "wavelength_nm": wavelength_nm,
            "fosc": fosc,
            "D2_au2": d2,
            "DX_au": dx,
            "DY_au": dy,
            "DZ_au": dz,
        }

    if set(parsed) != {1, 2}:
        raise RuntimeError(
            f"Expected electric-dipole rows for S1/S2; "
            f"found roots {sorted(parsed)}"
        )

    return parsed


def parse_state_energies(
    text: str,
) -> dict[int, dict[str, float]]:
    parsed: dict[int, dict[str, float]] = {}

    for match in STATE_ENERGY_RE.finditer(text):
        root = int(match.group("root"))

        if root not in (1, 2):
            continue

        parsed[root] = {
            "energy_au": float(match.group("energy_au")),
            "state_energy_eV": float(match.group("energy_eV")),
        }

    if set(parsed) != {1, 2}:
        raise RuntimeError(
            f"Expected STATE energies for S1/S2; "
            f"found roots {sorted(parsed)}"
        )

    return parsed


def is_float(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def parse_coordinate_line(
    line: str,
) -> tuple[str, float, float, float] | None:
    tokens = line.split()

    if len(tokens) >= 4:
        if (
            tokens[0].isalpha()
            and is_float(tokens[1])
            and is_float(tokens[2])
            and is_float(tokens[3])
        ):
            return (
                tokens[0],
                float(tokens[1]),
                float(tokens[2]),
                float(tokens[3]),
            )

    if len(tokens) >= 5:
        if (
            tokens[0].lstrip("+-").isdigit()
            and tokens[1].isalpha()
            and is_float(tokens[2])
            and is_float(tokens[3])
            and is_float(tokens[4])
        ):
            return (
                tokens[1],
                float(tokens[2]),
                float(tokens[3]),
                float(tokens[4]),
            )

    return None


def parse_coordinate_blocks(
    text: str,
) -> list[tuple[list[str], np.ndarray]]:
    lines = text.splitlines()
    blocks: list[tuple[list[str], np.ndarray]] = []

    for index, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" not in line:
            continue

        symbols: list[str] = []
        coordinates: list[list[float]] = []
        started = False

        for candidate in lines[index + 1 :]:
            parsed = parse_coordinate_line(candidate)

            if parsed is None:
                if started:
                    break
                continue

            started = True
            symbol, x, y, z = parsed
            symbols.append(symbol)
            coordinates.append([x, y, z])

        if symbols:
            blocks.append(
                (
                    symbols,
                    np.array(coordinates, dtype=np.float64),
                )
            )

    return blocks


def choose_qm_geometry(
    text: str,
) -> tuple[list[str], np.ndarray]:
    blocks = parse_coordinate_blocks(text)

    if not blocks:
        raise RuntimeError(
            "No CARTESIAN COORDINATES (ANGSTROEM) block found."
        )

    max_atoms = max(len(symbols) for symbols, _ in blocks)
    largest = [
        block
        for block in blocks
        if len(block[0]) == max_atoms
    ]

    symbols, coordinates = largest[-1]

    if max_atoms != 26:
        raise RuntimeError(
            f"Expected 26 QM atoms, found {max_atoms}."
        )

    if symbols.count(CARBON_SYMBOL) != 16:
        raise RuntimeError(
            f"Expected 16 carbon atoms, found "
            f"{symbols.count(CARBON_SYMBOL)}."
        )

    return symbols, coordinates


def minimum_interatomic_distance(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    differences = first[:, None, :] - second[None, :, :]
    distances_squared = np.sum(
        differences * differences,
        axis=2,
    )
    return math.sqrt(float(np.min(distances_squared)))


def main() -> None:
    selected = read_selected_outputs()
    frame_mapping = read_frame_mapping()

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    log("Day019 transition-dipole and geometry audit")
    log("Selected embedded outputs: 84/84")

    dipole_rows: list[dict[str, object]] = []
    geometry_rows: list[dict[str, object]] = []
    geometries: dict[
        tuple[int, str],
        tuple[list[str], np.ndarray],
    ] = {}

    maximum_d2_error = 0.0
    maximum_fosc_error = 0.0
    maximum_table_energy_error = 0.0

    for index, ((frame, site), path) in enumerate(
        sorted(selected.items()),
        start=1,
    ):
        text = path.read_text(errors="ignore")

        for marker in (
            NORMAL_MARKER,
            SCF_MARKER,
            TDDFT_MARKER,
        ):
            if marker not in text:
                raise RuntimeError(
                    f"Missing marker {marker!r}: {path}"
                )

        dipoles = parse_transition_dipoles(text)
        energies = parse_state_energies(text)
        symbols, coordinates = choose_qm_geometry(text)

        geometries[(frame, site)] = (
            symbols,
            coordinates,
        )

        symbol_array = np.array(symbols)
        carbon_mask = symbol_array == CARBON_SYMBOL
        carbon_coordinates = coordinates[carbon_mask]

        carbon_centroid = np.mean(
            carbon_coordinates,
            axis=0,
        )
        all_atom_centroid = np.mean(
            coordinates,
            axis=0,
        )

        geometry_rows.append(
            {
                "frame": frame,
                "time_ps": frame * FRAME_SPACING_PS,
                "site": site,
                "n_atoms": len(symbols),
                "n_carbon_atoms": int(np.count_nonzero(carbon_mask)),
                "carbon_centroid_x_A": carbon_centroid[0],
                "carbon_centroid_y_A": carbon_centroid[1],
                "carbon_centroid_z_A": carbon_centroid[2],
                "all_atom_centroid_x_A": all_atom_centroid[0],
                "all_atom_centroid_y_A": all_atom_centroid[1],
                "all_atom_centroid_z_A": all_atom_centroid[2],
                "source_output": str(path.relative_to(PROJECT_ROOT)),
            }
        )

        mapping = frame_mapping[(frame, site)]
        bright_root = int(mapping["bright_root"])

        for root in (1, 2):
            dipole = dipoles[root]
            state_energy = energies[root]

            dx = float(dipole["DX_au"])
            dy = float(dipole["DY_au"])
            dz = float(dipole["DZ_au"])
            d2 = float(dipole["D2_au2"])

            reconstructed_d2 = dx * dx + dy * dy + dz * dz
            d2_error = abs(d2 - reconstructed_d2)

            reconstructed_fosc = (
                (2.0 / 3.0)
                * float(state_energy["energy_au"])
                * d2
            )
            fosc_error = abs(
                float(dipole["fosc"])
                - reconstructed_fosc
            )

            table_energy_error = abs(
                float(dipole["table_energy_eV"])
                - float(state_energy["state_energy_eV"])
            )

            maximum_d2_error = max(
                maximum_d2_error,
                d2_error,
            )
            maximum_fosc_error = max(
                maximum_fosc_error,
                fosc_error,
            )
            maximum_table_energy_error = max(
                maximum_table_energy_error,
                table_energy_error,
            )

            family = (
                "bright_like"
                if root == bright_root
                else "alternate_like"
            )

            dipole_rows.append(
                {
                    "frame": frame,
                    "time_ps": frame * FRAME_SPACING_PS,
                    "site": site,
                    "root": root,
                    "family": family,
                    "is_bright_root": root == bright_root,
                    "energy_au": state_energy["energy_au"],
                    "state_energy_eV": state_energy["state_energy_eV"],
                    "table_energy_eV": dipole["table_energy_eV"],
                    "energy_cm1": dipole["energy_cm1"],
                    "wavelength_nm": dipole["wavelength_nm"],
                    "fosc": dipole["fosc"],
                    "D2_au2": d2,
                    "DX_au": dx,
                    "DY_au": dy,
                    "DZ_au": dz,
                    "dipole_magnitude_au": math.sqrt(max(d2, 0.0)),
                    "reconstructed_D2_au2": reconstructed_d2,
                    "D2_absolute_error": d2_error,
                    "reconstructed_fosc": reconstructed_fosc,
                    "fosc_absolute_error": fosc_error,
                    "table_state_energy_difference_eV": (
                        table_energy_error
                    ),
                    "carbon_centroid_x_A": carbon_centroid[0],
                    "carbon_centroid_y_A": carbon_centroid[1],
                    "carbon_centroid_z_A": carbon_centroid[2],
                    "source_output": str(path.relative_to(PROJECT_ROOT)),
                }
            )

        log(
            f"[{index:02d}/84] frame{frame:03d} {site} "
            f"S1|mu|={math.sqrt(dipoles[1]['D2_au2']):.6f} au "
            f"S2|mu|={math.sqrt(dipoles[2]['D2_au2']):.6f} au"
        )

    if len(dipole_rows) != 168:
        raise RuntimeError(
            f"Expected 168 dipole rows, found {len(dipole_rows)}."
        )

    if len(geometry_rows) != 84:
        raise RuntimeError(
            f"Expected 84 geometry rows, found {len(geometry_rows)}."
        )

    log("")
    log(
        f"Maximum D2 reconstruction error: "
        f"{maximum_d2_error:.8e} "
        f"(tolerance {D2_ABS_TOL:.8e})"
    )
    log(
        f"Maximum oscillator-strength reconstruction error: "
        f"{maximum_fosc_error:.8e} "
        f"(tolerance {FOSC_ABS_TOL:.8e})"
    )
    log(
        f"Maximum state/table energy discrepancy: "
        f"{maximum_table_energy_error:.8e} eV"
    )

    # Write the fully parsed observations before the final numerical
    # acceptance checks. This preserves diagnostics if a later tolerance
    # check fails for any reason.
    write_csv(DIPOLES_CSV, dipole_rows)
    write_csv(GEOMETRY_FRAMES_CSV, geometry_rows)

    if maximum_d2_error > D2_ABS_TOL:
        raise RuntimeError(
            f"D2 reconstruction failed after rounding-aware validation: "
            f"max error {maximum_d2_error:.8e} > "
            f"{D2_ABS_TOL:.8e}"
        )

    if maximum_fosc_error > FOSC_ABS_TOL:
        raise RuntimeError(
            f"Oscillator-strength reconstruction failed: max error "
            f"{maximum_fosc_error:.8e} > {FOSC_ABS_TOL:.8e}"
        )

    consistency_rows: list[dict[str, object]] = []

    for site in SITES:
        reference_symbols, reference_coordinates = geometries[(0, site)]
        carbon_mask = np.array(reference_symbols) == CARBON_SYMBOL
        reference_carbon_centroid = np.mean(
            reference_coordinates[carbon_mask],
            axis=0,
        )

        maximum_atom_displacement = 0.0
        maximum_carbon_centroid_shift = 0.0

        for frame in FRAMES:
            symbols, coordinates = geometries[(frame, site)]

            if symbols != reference_symbols:
                raise RuntimeError(
                    f"Atom-order mismatch for frame{frame:03d} {site}."
                )

            displacements = np.linalg.norm(
                coordinates - reference_coordinates,
                axis=1,
            )

            maximum_atom_displacement = max(
                maximum_atom_displacement,
                float(np.max(displacements)),
            )

            carbon_centroid = np.mean(
                coordinates[carbon_mask],
                axis=0,
            )

            maximum_carbon_centroid_shift = max(
                maximum_carbon_centroid_shift,
                float(
                    np.linalg.norm(
                        carbon_centroid
                        - reference_carbon_centroid
                    )
                ),
            )

        consistency_rows.append(
            {
                "site": site,
                "n_frames": 21,
                "n_atoms": len(reference_symbols),
                "n_carbon_atoms": int(
                    np.count_nonzero(carbon_mask)
                ),
                "maximum_atom_displacement_A": (
                    maximum_atom_displacement
                ),
                "maximum_carbon_centroid_shift_A": (
                    maximum_carbon_centroid_shift
                ),
                "frozen_geometry_pass": (
                    maximum_atom_displacement
                    <= GEOMETRY_TOL_A
                ),
            }
        )

    write_csv(
        GEOMETRY_CONSISTENCY_CSV,
        consistency_rows,
    )

    pair_rows: list[dict[str, object]] = []

    for site_a, site_b in combinations(SITES, 2):
        symbols_a, coordinates_a = geometries[(0, site_a)]
        symbols_b, coordinates_b = geometries[(0, site_b)]

        carbon_a = coordinates_a[
            np.array(symbols_a) == CARBON_SYMBOL
        ]
        carbon_b = coordinates_b[
            np.array(symbols_b) == CARBON_SYMBOL
        ]

        centroid_a = np.mean(carbon_a, axis=0)
        centroid_b = np.mean(carbon_b, axis=0)

        displacement = centroid_b - centroid_a
        distance = float(np.linalg.norm(displacement))

        if distance <= 0.0:
            raise RuntimeError(
                f"Nonpositive centroid distance: {site_a}, {site_b}"
            )

        pair_rows.append(
            {
                "site_a": site_a,
                "site_b": site_b,
                "carbon_centroid_distance_A": distance,
                "Rhat_x": displacement[0] / distance,
                "Rhat_y": displacement[1] / distance,
                "Rhat_z": displacement[2] / distance,
                "minimum_all_atom_distance_A": (
                    minimum_interatomic_distance(
                        coordinates_a,
                        coordinates_b,
                    )
                ),
                "minimum_carbon_distance_A": (
                    minimum_interatomic_distance(
                        carbon_a,
                        carbon_b,
                    )
                ),
            }
        )

    write_csv(SITE_PAIRS_CSV, pair_rows)

    frozen_pass_count = sum(
        bool(row["frozen_geometry_pass"])
        for row in consistency_rows
    )

    minimum_centroid_distance = min(
        float(row["carbon_centroid_distance_A"])
        for row in pair_rows
    )

    minimum_atomic_distance = min(
        float(row["minimum_all_atom_distance_A"])
        for row in pair_rows
    )

    bright_magnitudes = [
        float(row["dipole_magnitude_au"])
        for row in dipole_rows
        if bool(row["is_bright_root"])
    ]

    alternate_magnitudes = [
        float(row["dipole_magnitude_au"])
        for row in dipole_rows
        if not bool(row["is_bright_root"])
    ]

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 transition-dipole and geometry audit\n\n"
        )

        handle.write("## Validation\n\n")
        handle.write("- Embedded ORCA outputs audited: 84/84\n")
        handle.write("- S1/S2 transition-dipole records: 168/168\n")
        handle.write("- QM geometry records: 84/84\n")
        handle.write(
            f"- Frozen site geometries: {frozen_pass_count}/4 PASS\n"
        )
        handle.write(
            f"- Maximum D2 reconstruction error: "
            f"{maximum_d2_error:.8e} auÂ² "
            f"(tolerance {D2_ABS_TOL:.8e})\n"
        )
        handle.write(
            f"- Maximum oscillator-strength reconstruction error: "
            f"{maximum_fosc_error:.8e} "
            f"(tolerance {FOSC_ABS_TOL:.8e})\n"
        )
        handle.write(
            f"- Maximum state/table energy discrepancy: "
            f"{maximum_table_energy_error:.8e} eV\n\n"
        )

        handle.write("## Transition-dipole magnitudes\n\n")
        handle.write(
            f"- Bright-like |mu| range: "
            f"{min(bright_magnitudes):.6f} to "
            f"{max(bright_magnitudes):.6f} au\n"
        )
        handle.write(
            f"- Alternate-like |mu| range: "
            f"{min(alternate_magnitudes):.6f} to "
            f"{max(alternate_magnitudes):.6f} au\n\n"
        )

        handle.write("## Site-pair geometry\n\n")
        handle.write(
            "| Site A | Site B | Carbon-centroid distance (Ã) | "
            "Minimum atom distance (Ã) | Minimum C-C distance (Ã) |\n"
        )
        handle.write("|---|---|---:|---:|---:|\n")

        for row in pair_rows:
            handle.write(
                f"| {row['site_a']} "
                f"| {row['site_b']} "
                f"| {float(row['carbon_centroid_distance_A']):.6f} "
                f"| {float(row['minimum_all_atom_distance_A']):.6f} "
                f"| {float(row['minimum_carbon_distance_A']):.6f} |\n"
            )

        handle.write("\n## Coupling-model boundary\n\n")
        handle.write(
            "The audited electric transition dipoles and fixed site "
            "geometry are sufficient to construct a point transition-dipole "
            "coupling baseline for all 24 intersite state pairs per frame. "
            "That baseline is not automatically a final coupling model: its "
            "validity must be judged against chromophore size and separation, "
            "and it should be superseded by a transition-charge or transition-"
            "density treatment when short-range effects are material.\n"
        )

    log("")
    log("Day019 transition-dipole and geometry audit completed.")
    log("Embedded outputs: 84/84")
    log("Transition-dipole records: 168/168")
    log(f"Frozen geometries: {frozen_pass_count}/4")
    log(
        f"Maximum D2 error: {maximum_d2_error:.3e}"
    )
    log(
        f"Maximum fosc error: {maximum_fosc_error:.3e}"
    )
    log(
        f"Minimum carbon-centroid distance: "
        f"{minimum_centroid_distance:.6f} A"
    )
    log(
        f"Minimum interatomic distance: "
        f"{minimum_atomic_distance:.6f} A"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
