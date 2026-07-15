#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path
from statistics import mean, pstdev

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

AUDIT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_dipole_geometry_audit"
)

DIPOLES_CSV = AUDIT_ROOT / "transition_dipole_observations.csv"
SITE_PAIRS_CSV = AUDIT_ROOT / "site_pair_geometry.csv"

DIAGONAL_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_eight_state_diagonal_model"
)

DIAGONAL_LONG_CSV = DIAGONAL_ROOT / "eight_state_diagonal_long.csv"
BASIS_TXT = DIAGONAL_ROOT / "EIGHT_STATE_BASIS_ORDER_DAY019.txt"

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_point_transition_dipole_couplings"
)

ALIGNED_DIPOLES_CSV = OUTPUT_ROOT / "phase_aligned_transition_dipoles.csv"
GAUGE_AUDIT_CSV = OUTPUT_ROOT / "transition_dipole_gauge_audit.csv"
COUPLINGS_LONG_CSV = OUTPUT_ROOT / "point_dipole_couplings_long.csv"
COUPLING_STATS_CSV = OUTPUT_ROOT / "point_dipole_coupling_statistics.csv"
MEAN_MATRIX_CSV = OUTPUT_ROOT / "mean_coupling_matrix_meV.csv"
SD_MATRIX_CSV = OUTPUT_ROOT / "sd_coupling_matrix_meV.csv"
MAXABS_MATRIX_CSV = OUTPUT_ROOT / "maximum_absolute_coupling_matrix_meV.csv"
HAMILTONIAN_DIR = OUTPUT_ROOT / "hamiltonian_snapshots"
REPORT_MD = OUTPUT_ROOT / "POINT_TRANSITION_DIPOLE_BASELINE_DAY019.md"

BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_MEV = HARTREE_TO_EV * 1000.0
HARTREE_TO_CM1 = 219474.6313705

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
FAMILIES = ("alternate_like", "bright_like")
FRAMES = tuple(range(21))
FRAME_SPACING_PS = 5.0

MIN_VECTOR_NORM = 1.0e-12
GAUGE_WARNING_COSINE = 0.80

EXPECTED_BASIS = (
    ("PYR2_alternate", "PYR2", "alternate_like"),
    ("PYR2_bright", "PYR2", "bright_like"),
    ("PYR3_alternate", "PYR3", "alternate_like"),
    ("PYR3_bright", "PYR3", "bright_like"),
    ("PYR4_alternate", "PYR4", "alternate_like"),
    ("PYR4_bright", "PYR4", "bright_like"),
    ("PYR5_alternate", "PYR5", "alternate_like"),
    ("PYR5_bright", "PYR5", "bright_like"),
)


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build an unscreened/screened point-transition-dipole coupling "
            "baseline for the Day019 character-indexed eight-state manifold."
        )
    )
    parser.add_argument(
        "--relative-permittivity",
        type=float,
        default=1.0,
        help=(
            "Scalar screening divisor applied to point-dipole couplings. "
            "Use 1.0 for the unscreened baseline."
        ),
    )
    parser.add_argument(
        "--energy-zero",
        choices=("absolute", "global_minimum", "per_frame_minimum"),
        default="global_minimum",
        help="Energy-zero convention for exported full Hamiltonians.",
    )
    return parser.parse_args()


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
    labels: list[str],
    matrix: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["state", *labels])

        for label, values in zip(labels, matrix):
            writer.writerow(
                [label, *[f"{float(value):.12g}" for value in values]]
            )


def read_basis() -> list[dict[str, object]]:
    if not BASIS_TXT.is_file():
        raise SystemExit(f"Missing basis file: {BASIS_TXT}")

    rows: list[dict[str, object]] = []

    for line in BASIS_TXT.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()

        if not stripped or stripped.startswith("#"):
            continue

        tokens = stripped.split()

        if len(tokens) != 5:
            raise RuntimeError(f"Invalid basis line: {line!r}")

        index = int(tokens[0])
        label = tokens[1]
        site = tokens[2]
        family = tokens[3]
        root_token = tokens[4]

        rows.append(
            {
                "state_index": index,
                "state_label": label,
                "site": site,
                "family": family,
                "root": int(root_token.removeprefix("S")),
            }
        )

    rows.sort(key=lambda row: int(row["state_index"]))

    observed = tuple(
        (
            str(row["state_label"]),
            str(row["site"]),
            str(row["family"]),
        )
        for row in rows
    )

    if observed != EXPECTED_BASIS:
        raise RuntimeError(
            f"Unexpected eight-state basis:\n"
            f"observed={observed}\nexpected={EXPECTED_BASIS}"
        )

    return rows


def read_dipoles() -> list[dict[str, object]]:
    if not DIPOLES_CSV.is_file():
        raise SystemExit(f"Missing dipole table: {DIPOLES_CSV}")

    rows: list[dict[str, object]] = []

    with DIPOLES_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "frame": int(raw["frame"]),
                    "time_ps": float(raw["time_ps"]),
                    "site": raw["site"],
                    "root": int(raw["root"]),
                    "family": raw["family"],
                    "is_bright_root": as_bool(raw["is_bright_root"]),
                    "energy_eV": float(raw["state_energy_eV"]),
                    "fosc": float(raw["fosc"]),
                    "D2_au2": float(raw["D2_au2"]),
                    "DX_raw_au": float(raw["DX_au"]),
                    "DY_raw_au": float(raw["DY_au"]),
                    "DZ_raw_au": float(raw["DZ_au"]),
                    "dipole_magnitude_au": float(
                        raw["dipole_magnitude_au"]
                    ),
                }
            )

    if len(rows) != 168:
        raise RuntimeError(
            f"Expected 168 transition-dipole rows, found {len(rows)}."
        )

    expected = {
        (frame, site, family)
        for frame in FRAMES
        for site in SITES
        for family in FAMILIES
    }

    observed = {
        (
            int(row["frame"]),
            str(row["site"]),
            str(row["family"]),
        )
        for row in rows
    }

    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise RuntimeError(
            f"Dipole coverage mismatch. Missing={missing}; extra={extra}"
        )

    return rows


def read_site_pairs() -> dict[tuple[str, str], dict[str, object]]:
    if not SITE_PAIRS_CSV.is_file():
        raise SystemExit(
            f"Missing site-pair geometry: {SITE_PAIRS_CSV}"
        )

    pairs: dict[tuple[str, str], dict[str, object]] = {}

    with SITE_PAIRS_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        for raw in csv.DictReader(handle):
            site_a = raw["site_a"]
            site_b = raw["site_b"]

            key = (site_a, site_b)

            pairs[key] = {
                "site_a": site_a,
                "site_b": site_b,
                "distance_A": float(
                    raw["carbon_centroid_distance_A"]
                ),
                "rhat": np.array(
                    [
                        float(raw["Rhat_x"]),
                        float(raw["Rhat_y"]),
                        float(raw["Rhat_z"]),
                    ],
                    dtype=np.float64,
                ),
                "minimum_all_atom_distance_A": float(
                    raw["minimum_all_atom_distance_A"]
                ),
                "minimum_carbon_distance_A": float(
                    raw["minimum_carbon_distance_A"]
                ),
            }

    expected_pairs = set(combinations(SITES, 2))

    if set(pairs) != expected_pairs:
        raise RuntimeError(
            f"Site-pair coverage mismatch: "
            f"observed={sorted(pairs)}, expected={sorted(expected_pairs)}"
        )

    for key, row in pairs.items():
        rhat = np.asarray(row["rhat"], dtype=np.float64)
        norm = float(np.linalg.norm(rhat))

        if not math.isclose(norm, 1.0, abs_tol=1.0e-10):
            raise RuntimeError(
                f"Rhat is not normalized for {key}: norm={norm}"
            )

    return pairs


def read_diagonal_energies(
    basis: list[dict[str, object]],
) -> np.ndarray:
    if not DIAGONAL_LONG_CSV.is_file():
        raise SystemExit(
            f"Missing diagonal model table: {DIAGONAL_LONG_CSV}"
        )

    energies = np.full((21, 8), np.nan, dtype=np.float64)
    labels = [str(row["state_label"]) for row in basis]
    label_to_index = {
        label: index for index, label in enumerate(labels)
    }

    with DIAGONAL_LONG_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 168:
        raise RuntimeError(
            f"Expected 168 diagonal rows, found {len(rows)}."
        )

    for row in rows:
        frame = int(row["frame"])
        label = row["state_label"]

        if label not in label_to_index:
            raise RuntimeError(f"Unknown state label: {label}")

        index = label_to_index[label]

        if math.isfinite(energies[frame, index]):
            raise RuntimeError(
                f"Duplicate diagonal energy for frame={frame}, label={label}"
            )

        energies[frame, index] = float(row["energy_eV"])

    if not np.all(np.isfinite(energies)):
        raise RuntimeError("Incomplete diagonal-energy matrix.")

    return energies


def align_dipole_gauge(
    rows: list[dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
]:
    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in rows:
        grouped[
            (str(row["site"]), str(row["family"]))
        ].append(row)

    aligned_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []

    for key in sorted(grouped):
        group = sorted(
            grouped[key],
            key=lambda row: int(row["frame"]),
        )

        if len(group) != 21:
            raise RuntimeError(
                f"Expected 21 dipoles for {key}, found {len(group)}."
            )

        reference_raw = np.array(
            [
                float(group[0]["DX_raw_au"]),
                float(group[0]["DY_raw_au"]),
                float(group[0]["DZ_raw_au"]),
            ],
            dtype=np.float64,
        )

        reference_norm = float(np.linalg.norm(reference_raw))

        if reference_norm <= MIN_VECTOR_NORM:
            raise RuntimeError(
                f"Reference transition dipole is zero for {key}."
            )

        n_flips = 0
        minimum_abs_raw_cosine = 1.0
        minimum_aligned_cosine = 1.0

        for row in group:
            raw = np.array(
                [
                    float(row["DX_raw_au"]),
                    float(row["DY_raw_au"]),
                    float(row["DZ_raw_au"]),
                ],
                dtype=np.float64,
            )

            raw_norm = float(np.linalg.norm(raw))

            if raw_norm <= MIN_VECTOR_NORM:
                raise RuntimeError(
                    f"Zero transition dipole for "
                    f"frame={row['frame']}, key={key}"
                )

            raw_cosine = float(
                np.dot(raw, reference_raw)
                / (raw_norm * reference_norm)
            )

            sign_factor = 1.0 if raw_cosine >= 0.0 else -1.0
            aligned = sign_factor * raw
            aligned_cosine = abs(raw_cosine)

            if sign_factor < 0.0:
                n_flips += 1

            minimum_abs_raw_cosine = min(
                minimum_abs_raw_cosine,
                abs(raw_cosine),
            )
            minimum_aligned_cosine = min(
                minimum_aligned_cosine,
                aligned_cosine,
            )

            aligned_rows.append(
                {
                    **row,
                    "reference_frame": 0,
                    "gauge_sign_factor": int(sign_factor),
                    "raw_cosine_to_reference": raw_cosine,
                    "aligned_cosine_to_reference": aligned_cosine,
                    "DX_aligned_au": aligned[0],
                    "DY_aligned_au": aligned[1],
                    "DZ_aligned_au": aligned[2],
                }
            )

        audit_rows.append(
            {
                "site": key[0],
                "family": key[1],
                "reference_frame": 0,
                "n_frames": 21,
                "n_sign_flips": n_flips,
                "minimum_absolute_raw_cosine_to_reference": (
                    minimum_abs_raw_cosine
                ),
                "minimum_aligned_cosine_to_reference": (
                    minimum_aligned_cosine
                ),
                "gauge_stability_pass": (
                    minimum_aligned_cosine
                    >= GAUGE_WARNING_COSINE
                ),
            }
        )

    aligned_rows.sort(
        key=lambda row: (
            int(row["frame"]),
            str(row["site"]),
            str(row["family"]),
        )
    )

    return aligned_rows, audit_rows


def coupling_au(
    mu_a: np.ndarray,
    mu_b: np.ndarray,
    rhat: np.ndarray,
    distance_bohr: float,
    relative_permittivity: float,
) -> tuple[float, float]:
    norm_a = float(np.linalg.norm(mu_a))
    norm_b = float(np.linalg.norm(mu_b))

    if norm_a <= MIN_VECTOR_NORM or norm_b <= MIN_VECTOR_NORM:
        raise RuntimeError("Zero transition dipole in coupling evaluation.")

    orientation_factor = float(
        np.dot(mu_a / norm_a, mu_b / norm_b)
        - 3.0
        * np.dot(mu_a / norm_a, rhat)
        * np.dot(mu_b / norm_b, rhat)
    )

    numerator = float(
        np.dot(mu_a, mu_b)
        - 3.0
        * np.dot(mu_a, rhat)
        * np.dot(mu_b, rhat)
    )

    value = (
        numerator
        / (distance_bohr**3)
        / relative_permittivity
    )

    return value, orientation_factor


def main() -> None:
    args = parse_args()

    if args.relative_permittivity <= 0.0:
        raise SystemExit(
            "--relative-permittivity must be strictly positive."
        )

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    HAMILTONIAN_DIR.mkdir(parents=True, exist_ok=True)

    basis = read_basis()
    labels = [str(row["state_label"]) for row in basis]
    state_index = {
        (
            str(row["site"]),
            str(row["family"]),
        ): int(row["state_index"])
        for row in basis
    }

    dipole_rows = read_dipoles()
    site_pairs = read_site_pairs()
    diagonal_energies = read_diagonal_energies(basis)

    aligned_rows, gauge_audit = align_dipole_gauge(
        dipole_rows
    )

    write_csv(ALIGNED_DIPOLES_CSV, aligned_rows)
    write_csv(GAUGE_AUDIT_CSV, gauge_audit)

    aligned_lookup = {
        (
            int(row["frame"]),
            str(row["site"]),
            str(row["family"]),
        ): np.array(
            [
                float(row["DX_aligned_au"]),
                float(row["DY_aligned_au"]),
                float(row["DZ_aligned_au"]),
            ],
            dtype=np.float64,
        )
        for row in aligned_rows
    }

    coupling_rows: list[dict[str, object]] = []
    coupling_matrices_meV = np.zeros(
        (21, 8, 8),
        dtype=np.float64,
    )

    log("Day019 point-transition-dipole coupling baseline")
    log(
        f"Relative permittivity divisor: "
        f"{args.relative_permittivity:.6f}"
    )
    log("Gauge: each site/family aligned to its frame000 dipole")

    for frame in FRAMES:
        frame_values: list[float] = []

        for site_a, site_b in combinations(SITES, 2):
            geometry = site_pairs[(site_a, site_b)]
            rhat = np.asarray(geometry["rhat"], dtype=np.float64)
            distance_A = float(geometry["distance_A"])
            distance_bohr = distance_A / BOHR_TO_ANGSTROM

            for family_a, family_b in product(
                FAMILIES,
                FAMILIES,
            ):
                mu_a = aligned_lookup[
                    (frame, site_a, family_a)
                ]
                mu_b = aligned_lookup[
                    (frame, site_b, family_b)
                ]

                j_au, kappa = coupling_au(
                    mu_a=mu_a,
                    mu_b=mu_b,
                    rhat=rhat,
                    distance_bohr=distance_bohr,
                    relative_permittivity=args.relative_permittivity,
                )

                j_eV = j_au * HARTREE_TO_EV
                j_meV = j_au * HARTREE_TO_MEV
                j_cm1 = j_au * HARTREE_TO_CM1

                index_a = state_index[(site_a, family_a)]
                index_b = state_index[(site_b, family_b)]

                coupling_matrices_meV[
                    frame,
                    index_a,
                    index_b,
                ] = j_meV
                coupling_matrices_meV[
                    frame,
                    index_b,
                    index_a,
                ] = j_meV

                frame_values.append(j_meV)

                coupling_rows.append(
                    {
                        "frame": frame,
                        "time_ps": frame * FRAME_SPACING_PS,
                        "site_a": site_a,
                        "family_a": family_a,
                        "state_index_a": index_a,
                        "state_label_a": labels[index_a],
                        "site_b": site_b,
                        "family_b": family_b,
                        "state_index_b": index_b,
                        "state_label_b": labels[index_b],
                        "distance_A": distance_A,
                        "distance_bohr": distance_bohr,
                        "minimum_all_atom_distance_A": geometry[
                            "minimum_all_atom_distance_A"
                        ],
                        "relative_permittivity": (
                            args.relative_permittivity
                        ),
                        "orientation_factor_kappa": kappa,
                        "mu_a_magnitude_au": float(
                            np.linalg.norm(mu_a)
                        ),
                        "mu_b_magnitude_au": float(
                            np.linalg.norm(mu_b)
                        ),
                        "J_hartree": j_au,
                        "J_eV": j_eV,
                        "J_meV": j_meV,
                        "J_cm-1": j_cm1,
                        "absolute_J_meV": abs(j_meV),
                    }
                )

        log(
            f"[frame {frame:03d}/020] "
            f"couplings=24 "
            f"range={min(frame_values):+.6f} to "
            f"{max(frame_values):+.6f} meV "
            f"max|J|={max(abs(v) for v in frame_values):.6f} meV"
        )

    if len(coupling_rows) != 504:
        raise RuntimeError(
            f"Expected 504 couplings, found {len(coupling_rows)}."
        )

    for frame in FRAMES:
        matrix = coupling_matrices_meV[frame]

        if not np.allclose(
            matrix,
            matrix.T,
            atol=1.0e-14,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"Coupling matrix is not symmetric at frame {frame}."
            )

        if not np.allclose(
            np.diag(matrix),
            0.0,
            atol=1.0e-14,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"Coupling-matrix diagonal is not zero at frame {frame}."
            )

    write_csv(COUPLINGS_LONG_CSV, coupling_rows)

    grouped: dict[
        tuple[str, str, str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in coupling_rows:
        grouped[
            (
                str(row["site_a"]),
                str(row["family_a"]),
                str(row["site_b"]),
                str(row["family_b"]),
            )
        ].append(row)

    statistics_rows: list[dict[str, object]] = []

    for key in sorted(grouped):
        group = grouped[key]
        values = [float(row["J_meV"]) for row in group]
        absolute_values = [abs(value) for value in values]
        kappas = [
            float(row["orientation_factor_kappa"])
            for row in group
        ]

        statistics_rows.append(
            {
                "site_a": key[0],
                "family_a": key[1],
                "site_b": key[2],
                "family_b": key[3],
                "n_frames": len(group),
                "distance_A": float(group[0]["distance_A"]),
                "relative_permittivity": (
                    args.relative_permittivity
                ),
                "mean_J_meV": mean(values),
                "sd_J_meV": pstdev(values),
                "minimum_J_meV": min(values),
                "maximum_J_meV": max(values),
                "mean_absolute_J_meV": mean(absolute_values),
                "maximum_absolute_J_meV": max(absolute_values),
                "mean_kappa": mean(kappas),
                "minimum_kappa": min(kappas),
                "maximum_kappa": max(kappas),
            }
        )

    write_csv(COUPLING_STATS_CSV, statistics_rows)

    mean_matrix = np.mean(coupling_matrices_meV, axis=0)
    sd_matrix = np.std(coupling_matrices_meV, axis=0)
    maxabs_matrix = np.max(
        np.abs(coupling_matrices_meV),
        axis=0,
    )

    write_matrix_csv(MEAN_MATRIX_CSV, labels, mean_matrix)
    write_matrix_csv(SD_MATRIX_CSV, labels, sd_matrix)
    write_matrix_csv(MAXABS_MATRIX_CSV, labels, maxabs_matrix)

    if args.energy_zero == "absolute":
        diagonal_for_export = diagonal_energies.copy()
        zero_description = "absolute TDDFT excitation energies"
    elif args.energy_zero == "global_minimum":
        global_minimum = float(np.min(diagonal_energies))
        diagonal_for_export = (
            diagonal_energies - global_minimum
        )
        zero_description = (
            f"global minimum shifted to zero "
            f"(E0={global_minimum:.9f} eV)"
        )
    else:
        frame_minimum = np.min(
            diagonal_energies,
            axis=1,
            keepdims=True,
        )
        diagonal_for_export = (
            diagonal_energies - frame_minimum
        )
        zero_description = (
            "minimum state energy shifted to zero in each frame"
        )

    for frame in FRAMES:
        hamiltonian_eV = np.diag(
            diagonal_for_export[frame]
        )

        hamiltonian_eV += (
            coupling_matrices_meV[frame] / 1000.0
        )

        output = (
            HAMILTONIAN_DIR
            / f"H_point_dipole_frame{frame:03d}.dat"
        )

        with output.open("w", encoding="utf-8") as handle:
            handle.write(
                f"# Day019 eight-state point-transition-dipole baseline, "
                f"frame={frame:03d}, "
                f"time_ps={frame * FRAME_SPACING_PS:.6f}\n"
            )
            handle.write(
                f"# relative_permittivity="
                f"{args.relative_permittivity:.9f}\n"
            )
            handle.write(
                f"# energy_zero: {zero_description}\n"
            )
            handle.write(
                "# same-site alternate/bright electronic coupling is zero "
                "in the local adiabatic basis; derivative couplings are "
                "not included\n"
            )
            handle.write(
                "# basis: " + " ".join(labels) + "\n"
            )
            np.savetxt(handle, hamiltonian_eV, fmt="%.12f")

    all_values = np.array(
        [float(row["J_meV"]) for row in coupling_rows],
        dtype=np.float64,
    )

    all_abs = np.abs(all_values)
    maximum_abs_j = float(np.max(all_abs))
    mean_abs_j = float(np.mean(all_abs))
    sd_all_j = float(np.std(all_values))

    bright_bright = [
        abs(float(row["J_meV"]))
        for row in coupling_rows
        if row["family_a"] == "bright_like"
        and row["family_b"] == "bright_like"
    ]

    mixed = [
        abs(float(row["J_meV"]))
        for row in coupling_rows
        if (
            row["family_a"] != row["family_b"]
        )
    ]

    alternate_alternate = [
        abs(float(row["J_meV"]))
        for row in coupling_rows
        if row["family_a"] == "alternate_like"
        and row["family_b"] == "alternate_like"
    ]

    minimum_gauge_cosine = min(
        float(
            row["minimum_aligned_cosine_to_reference"]
        )
        for row in gauge_audit
    )

    gauge_pass_count = sum(
        bool(row["gauge_stability_pass"])
        for row in gauge_audit
    )

    minimum_site_distance = min(
        float(row["distance_A"])
        for row in site_pairs.values()
    )

    minimum_atomic_distance = min(
        float(row["minimum_all_atom_distance_A"])
        for row in site_pairs.values()
    )

    minimum_local_gap_meV = 53.0
    maximum_diagonal_sd_meV = 16.034118950586075

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 point-transition-dipole coupling baseline\n\n"
        )

        handle.write("## Scope\n\n")
        handle.write(
            "- Character-indexed eight-state local basis.\n"
        )
        handle.write(
            "- Six intersite pairs Ã four state-family combinations Ã "
            "21 frames = 504 couplings.\n"
        )
        handle.write(
            "- Local-state transition-dipole signs were gauge-aligned to "
            "the corresponding frame000 vector before evaluating signed "
            "couplings.\n"
        )
        handle.write(
            "- Same-site alternate/bright electronic couplings are zero "
            "in the local adiabatic basis used here. Time-derivative "
            "nonadiabatic couplings are not included.\n\n"
        )

        handle.write("## Numerical controls\n\n")
        handle.write("- Transition-dipole observations: 168/168\n")
        handle.write("- Intersite coupling observations: 504/504\n")
        handle.write("- Hamiltonian snapshots: 21/21\n")
        handle.write(
            f"- Relative-permittivity divisor: "
            f"{args.relative_permittivity:.6f}\n"
        )
        handle.write(
            f"- Gauge-stability groups passing cosine >= "
            f"{GAUGE_WARNING_COSINE:.2f}: {gauge_pass_count}/8\n"
        )
        handle.write(
            f"- Minimum aligned cosine to frame000: "
            f"{minimum_gauge_cosine:.6f}\n"
        )
        handle.write(
            f"- Minimum site-centroid distance: "
            f"{minimum_site_distance:.6f} Ã\n"
        )
        handle.write(
            f"- Minimum interatomic distance: "
            f"{minimum_atomic_distance:.6f} Ã\n\n"
        )

        handle.write("## Coupling magnitudes\n\n")
        handle.write(
            f"- Overall mean |J|: {mean_abs_j:.6f} meV\n"
        )
        handle.write(
            f"- Overall maximum |J|: {maximum_abs_j:.6f} meV\n"
        )
        handle.write(
            f"- Overall signed-J SD: {sd_all_j:.6f} meV\n"
        )
        handle.write(
            f"- Bright-bright mean |J|: "
            f"{mean(bright_bright):.6f} meV\n"
        )
        handle.write(
            f"- Bright-bright maximum |J|: "
            f"{max(bright_bright):.6f} meV\n"
        )
        handle.write(
            f"- Mixed-family mean |J|: "
            f"{mean(mixed):.6f} meV\n"
        )
        handle.write(
            f"- Alternate-alternate mean |J|: "
            f"{mean(alternate_alternate):.6f} meV\n\n"
        )

        handle.write("## Scale comparison\n\n")
        handle.write(
            f"- Minimum local S1-S2 gap: "
            f"{minimum_local_gap_meV:.3f} meV\n"
        )
        handle.write(
            f"- Maximum diagonal energy SD: "
            f"{maximum_diagonal_sd_meV:.3f} meV\n"
        )
        handle.write(
            f"- max|J| / minimum local gap: "
            f"{maximum_abs_j / minimum_local_gap_meV:.6f}\n"
        )
        handle.write(
            f"- max|J| / maximum diagonal SD: "
            f"{maximum_abs_j / maximum_diagonal_sd_meV:.6f}\n\n"
        )

        handle.write("## Pair-resolved statistics\n\n")
        handle.write(
            "| Site A | Family A | Site B | Family B | Distance (Ã) | "
            "Mean J (meV) | SD (meV) | Mean |J| (meV) | "
            "Max |J| (meV) |\n"
        )
        handle.write(
            "|---|---|---|---|---:|---:|---:|---:|---:|\n"
        )

        for row in statistics_rows:
            handle.write(
                f"| {row['site_a']} "
                f"| {row['family_a']} "
                f"| {row['site_b']} "
                f"| {row['family_b']} "
                f"| {float(row['distance_A']):.6f} "
                f"| {float(row['mean_J_meV']):.6f} "
                f"| {float(row['sd_J_meV']):.6f} "
                f"| {float(row['mean_absolute_J_meV']):.6f} "
                f"| {float(row['maximum_absolute_J_meV']):.6f} |\n"
            )

        handle.write("\n## Interpretation boundary\n\n")
        handle.write(
            "These are point-transition-dipole couplings and constitute a "
            "controlled baseline, not a transition-density benchmark. The "
            "minimum interatomic separation exceeds 20 Ã, which supports "
            "using the dipolar term as a first approximation; nevertheless, "
            "finite-size multipolar and dielectric-screening effects remain "
            "unquantified. Coupling signs are reported in the explicit "
            "frame000 local-state gauge and are not independently observable "
            "under arbitrary local basis-phase changes.\n"
        )

    log("")
    log("Day019 point-transition-dipole baseline completed.")
    log("Dipole observations: 168/168")
    log("Coupling observations: 504/504")
    log("Hamiltonian snapshots: 21/21")
    log(
        f"Gauge-stability groups: {gauge_pass_count}/8 "
        f"(minimum cosine {minimum_gauge_cosine:.6f})"
    )
    log(
        f"Overall |J| range: "
        f"{float(np.min(all_abs)):.6f}-"
        f"{maximum_abs_j:.6f} meV"
    )
    log(
        f"Bright-bright maximum |J|: "
        f"{max(bright_bright):.6f} meV"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
