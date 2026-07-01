#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HAMILTONIAN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians/"
    "hamiltonian_snapshots_bright4"
)

DYNAMICS_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_haken_strobl_high_dephasing"
)

NPZ_PATH = DYNAMICS_ROOT / "dephasing_population_trajectories.npz"

OUTPUT_ROOT = DYNAMICS_ROOT / "zeno_scaling_audit"
PAIR_RATES_CSV = OUTPUT_ROOT / "pair_rate_summary.csv"
SCALING_CSV = OUTPUT_ROOT / "pair_scaling_summary.csv"
CLASSICAL_CSV = OUTPUT_ROOT / "classical_reduction_validation.csv"
REPORT_MD = OUTPUT_ROOT / "ZENO_SCALING_AUDIT_DAY020.md"

HBAR_EV_PS = 6.582119569e-4
SELECTED_GAMMAS = (20.0, 50.0, 100.0, 200.0)
EXPECTED_FRAMES = 21
N_STATES = 4

HIGH_INDICES = (0, 1, 2)
PYR5_INDEX = 3


def log(message: str = "") -> None:
    print(message, flush=True)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def classical_generator(
    hamiltonian_ev: np.ndarray,
    gamma_ps_inv: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    energies = np.diag(hamiltonian_ev)
    rates = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    detunings = np.zeros((N_STATES, N_STATES), dtype=np.float64)

    for i in range(N_STATES):
        for j in range(i + 1, N_STATES):
            delta_omega = abs(energies[i] - energies[j]) / HBAR_EV_PS
            coupling_frequency = abs(hamiltonian_ev[i, j]) / HBAR_EV_PS

            rate = (
                2.0
                * coupling_frequency**2
                * gamma_ps_inv
                / (gamma_ps_inv**2 + delta_omega**2)
            )

            rates[i, j] = rate
            rates[j, i] = rate
            detunings[i, j] = delta_omega
            detunings[j, i] = delta_omega

    generator = rates.copy()

    for source in range(N_STATES):
        generator[source, source] = -np.sum(rates[:, source])

    return generator, rates, detunings


def propagate_classical(
    generator: np.ndarray,
    times_ps: np.ndarray,
    initial_index: int,
) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(generator)

    initial = np.zeros(N_STATES, dtype=np.float64)
    initial[initial_index] = 1.0

    coefficients = eigenvectors.T @ initial
    exponentials = np.exp(np.outer(eigenvalues, times_ps))

    populations = (
        (eigenvectors * coefficients[np.newaxis, :])
        @ exponentials
    ).T

    return populations


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    hamiltonian_files = sorted(
        HAMILTONIAN_ROOT.glob("H_bright4_tdcac_frame*.dat")
    )

    if len(hamiltonian_files) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} Hamiltonians, "
            f"found {len(hamiltonian_files)}"
        )

    hamiltonians = np.stack(
        [np.loadtxt(path, comments="#") for path in hamiltonian_files]
    )

    if hamiltonians.shape != (
        EXPECTED_FRAMES,
        N_STATES,
        N_STATES,
    ):
        raise RuntimeError(
            f"Unexpected Hamiltonian shape: {hamiltonians.shape}"
        )

    symmetry_error = float(
        np.max(
            np.abs(
                hamiltonians
                - np.swapaxes(hamiltonians, 1, 2)
            )
        )
    )

    if symmetry_error > 1.0e-12:
        raise RuntimeError(
            f"Hamiltonian symmetry failure: {symmetry_error:.3e}"
        )

    with np.load(NPZ_PATH) as data:
        gammas = data["gamma_phi_ps_inv"].astype(np.float64)
        times_ps = data["times_ps"].astype(np.float64)
        site_labels = data["site_labels"].astype(str)
        populations = data["populations"].astype(np.float64)

    expected_shape = (
        gammas.size,
        N_STATES,
        EXPECTED_FRAMES,
        times_ps.size,
        N_STATES,
    )

    if populations.shape != expected_shape:
        raise RuntimeError(
            f"Unexpected population shape: {populations.shape}; "
            f"expected {expected_shape}"
        )

    gamma_indices: dict[float, int] = {}

    for gamma in SELECTED_GAMMAS:
        matches = np.where(np.isclose(gammas, gamma))[0]

        if matches.size != 1:
            raise RuntimeError(
                f"Could not uniquely locate gamma={gamma}"
            )

        gamma_indices[gamma] = int(matches[0])

    pairs = [
        (i, j)
        for i in range(N_STATES)
        for j in range(i + 1, N_STATES)
    ]

    pair_rate_rows: list[dict[str, object]] = []
    pair_scaling_rows: list[dict[str, object]] = []
    classical_rows: list[dict[str, object]] = []

    pair_data: dict[
        tuple[int, int],
        dict[float, list[float]],
    ] = {
        pair: {gamma: [] for gamma in SELECTED_GAMMAS}
        for pair in pairs
    }

    pair_detunings: dict[
        tuple[int, int],
        list[float],
    ] = {
        pair: []
        for pair in pairs
    }

    maximum_classical_norm_error = 0.0

    for gamma in SELECTED_GAMMAS:
        gamma_index = gamma_indices[gamma]
        all_errors: list[float] = []
        high_errors: list[float] = []
        pyr5_errors: list[float] = []
        maximum_errors: list[float] = []

        start_time_ps = max(5.0 / gamma, 0.05)
        time_mask = times_ps >= start_time_ps

        for frame_index, hamiltonian in enumerate(hamiltonians):
            generator, rates, detunings = classical_generator(
                hamiltonian,
                gamma,
            )

            if np.min(rates) < -1.0e-15:
                raise RuntimeError("Negative classical rate detected")

            for i, j in pairs:
                pair_data[(i, j)][gamma].append(
                    float(rates[i, j])
                )

                if gamma == SELECTED_GAMMAS[0]:
                    pair_detunings[(i, j)].append(
                        float(detunings[i, j])
                    )

            for initial_index in range(N_STATES):
                classical_population = propagate_classical(
                    generator,
                    times_ps,
                    initial_index,
                )

                norm_error = float(
                    np.max(
                        np.abs(
                            np.sum(
                                classical_population,
                                axis=1,
                            )
                            - 1.0
                        )
                    )
                )

                maximum_classical_norm_error = max(
                    maximum_classical_norm_error,
                    norm_error,
                )

                quantum_population = populations[
                    gamma_index,
                    initial_index,
                    frame_index,
                    :,
                    :,
                ]

                difference = (
                    quantum_population[time_mask]
                    - classical_population[time_mask]
                )

                all_errors.extend(
                    np.ravel(difference).tolist()
                )

                if initial_index in HIGH_INDICES:
                    high_difference = difference[
                        :,
                        list(HIGH_INDICES),
                    ]

                    high_errors.extend(
                        np.ravel(high_difference).tolist()
                    )

                    pyr5_errors.extend(
                        difference[:, PYR5_INDEX].tolist()
                    )

                maximum_errors.append(
                    float(np.max(np.abs(difference)))
                )

        all_errors_array = np.asarray(
            all_errors,
            dtype=np.float64,
        )
        high_errors_array = np.asarray(
            high_errors,
            dtype=np.float64,
        )
        pyr5_errors_array = np.asarray(
            pyr5_errors,
            dtype=np.float64,
        )

        classical_rows.append(
            {
                "gamma_phi_ps_inv": gamma,
                "comparison_start_time_ps": start_time_ps,
                "population_rmse_all": float(
                    np.sqrt(
                        np.mean(all_errors_array**2)
                    )
                ),
                "population_rmse_high_manifold": float(
                    np.sqrt(
                        np.mean(high_errors_array**2)
                    )
                ),
                "population_rmse_PYR5_from_high_states": float(
                    np.sqrt(
                        np.mean(pyr5_errors_array**2)
                    )
                ),
                "maximum_absolute_population_error": max(
                    maximum_errors
                ),
            }
        )

    for i, j in pairs:
        category = (
            "high_high"
            if i in HIGH_INDICES and j in HIGH_INDICES
            else "includes_PYR5"
        )

        detuning_values = np.asarray(
            pair_detunings[(i, j)],
            dtype=np.float64,
        )

        for gamma in SELECTED_GAMMAS:
            rate_values = np.asarray(
                pair_data[(i, j)][gamma],
                dtype=np.float64,
            )

            pair_rate_rows.append(
                {
                    "site_i": site_labels[i],
                    "site_j": site_labels[j],
                    "category": category,
                    "gamma_phi_ps_inv": gamma,
                    "mean_rate_ps_inv": float(
                        np.mean(rate_values)
                    ),
                    "sd_rate_ps_inv": float(
                        np.std(rate_values)
                    ),
                    "minimum_rate_ps_inv": float(
                        np.min(rate_values)
                    ),
                    "maximum_rate_ps_inv": float(
                        np.max(rate_values)
                    ),
                }
            )

        mean_rates = {
            gamma: float(
                np.mean(
                    np.asarray(
                        pair_data[(i, j)][gamma],
                        dtype=np.float64,
                    )
                )
            )
            for gamma in SELECTED_GAMMAS
        }

        tested_peak_gamma = max(
            SELECTED_GAMMAS,
            key=lambda gamma: mean_rates[gamma],
        )

        k100 = mean_rates[100.0]
        k200 = mean_rates[200.0]

        pair_scaling_rows.append(
            {
                "site_i": site_labels[i],
                "site_j": site_labels[j],
                "category": category,
                "mean_predicted_turnover_gamma_ps_inv": float(
                    np.mean(detuning_values)
                ),
                "minimum_predicted_turnover_gamma_ps_inv": float(
                    np.min(detuning_values)
                ),
                "maximum_predicted_turnover_gamma_ps_inv": float(
                    np.max(detuning_values)
                ),
                "tested_gamma_maximizing_mean_rate_ps_inv": (
                    tested_peak_gamma
                ),
                "mean_rate_at_20_ps_inv": mean_rates[20.0],
                "mean_rate_at_50_ps_inv": mean_rates[50.0],
                "mean_rate_at_100_ps_inv": k100,
                "mean_rate_at_200_ps_inv": k200,
                "k100_over_k200": (
                    k100 / k200
                    if k200 > 0.0
                    else float("nan")
                ),
                "gamma_k_100": 100.0 * k100,
                "gamma_k_200": 200.0 * k200,
                "fraction_frames_turnover_below_200": float(
                    np.mean(detuning_values < 200.0)
                ),
            }
        )

    write_csv(PAIR_RATES_CSV, pair_rate_rows)
    write_csv(SCALING_CSV, pair_scaling_rows)
    write_csv(CLASSICAL_CSV, classical_rows)

    high_high_rows = [
        row
        for row in pair_scaling_rows
        if row["category"] == "high_high"
    ]

    pyr5_rows = [
        row
        for row in pair_scaling_rows
        if row["category"] == "includes_PYR5"
    ]

    high_turnover_mean = float(
        np.mean(
            [
                float(
                    row[
                        "mean_predicted_turnover_gamma_ps_inv"
                    ]
                )
                for row in high_high_rows
            ]
        )
    )

    pyr5_turnover_mean = float(
        np.mean(
            [
                float(
                    row[
                        "mean_predicted_turnover_gamma_ps_inv"
                    ]
                )
                for row in pyr5_rows
            ]
        )
    )

    high_k_ratio_mean = float(
        np.mean(
            [
                float(row["k100_over_k200"])
                for row in high_high_rows
            ]
        )
    )

    pyr5_k_ratio_mean = float(
        np.mean(
            [
                float(row["k100_over_k200"])
                for row in pyr5_rows
            ]
        )
    )

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day020 High-Dephasing and Zeno-Scaling Audit\n\n"
        )

        handle.write("## Model\n\n")
        handle.write(
            "The effective Haken-Strobl pair rate was evaluated as\n\n"
        )
        handle.write(
            "\\[\n"
            "k_{ij} = "
            "\\frac{2(J_{ij}/\\hbar)^2\\gamma_\\phi}"
            "{\\gamma_\\phi^2 + "
            "(\\Delta E_{ij}/\\hbar)^2}.\n"
            "\\]\n\n"
        )

        handle.write(
            "The rate is maximal when "
            "`gamma_phi = |Delta E|/hbar`.\n\n"
        )

        handle.write("## Validation\n\n")
        handle.write(
            f"- Hamiltonians: {len(hamiltonian_files)}/"
            f"{EXPECTED_FRAMES}.\n"
        )
        handle.write(
            f"- Maximum Hamiltonian symmetry error: "
            f"{symmetry_error:.3e} eV.\n"
        )
        handle.write(
            f"- Maximum classical probability-norm error: "
            f"{maximum_classical_norm_error:.3e}.\n"
        )
        handle.write(
            "- Selected gamma values: "
            + ", ".join(
                f"{value:g}" for value in SELECTED_GAMMAS
            )
            + " ps^-1.\n\n"
        )

        handle.write("## Characteristic turnover scales\n\n")
        handle.write(
            f"- Mean PYR2-PYR4 predicted turnover scale: "
            f"{high_turnover_mean:.3f} ps^-1.\n"
        )
        handle.write(
            f"- Mean turnover scale for pairs involving PYR5: "
            f"{pyr5_turnover_mean:.3f} ps^-1.\n"
        )
        handle.write(
            f"- Mean high-manifold k(100)/k(200): "
            f"{high_k_ratio_mean:.6f}.\n"
        )
        handle.write(
            f"- Mean PYR5-pair k(100)/k(200): "
            f"{pyr5_k_ratio_mean:.6f}.\n\n"
        )

        handle.write("## Interpretation\n\n")
        handle.write(
            "The PYR2-PYR4 pairs have detuning frequencies "
            "within or below the tested 20-50 ps^-1 turnover "
            "region. Their effective rates therefore decrease "
            "at sufficiently large gamma_phi, consistent with "
            "high-dephasing or Zeno suppression.\n\n"
        )
        handle.write(
            "Pairs involving PYR5 have substantially larger "
            "detuning frequencies because PYR5 lies roughly "
            "300 meV below the high-energy manifold. The tested "
            "range up to 200 ps^-1 does not necessarily exceed "
            "their pair-specific turnover scales. Their rates "
            "and PYR5 populations can therefore continue to "
            "increase while PYR2-PYR4 transfer is already being "
            "suppressed.\n\n"
        )
        handle.write(
            "The PYR5 increase remains a phenomenological "
            "pure-dephasing result. It is not thermal downhill "
            "relaxation because the model has no detailed "
            "balance or population-relaxation bath.\n\n"
        )

        handle.write("## Classical-reduction comparison\n\n")
        handle.write(
            "| gamma_phi (ps^-1) | Start time (ps) | "
            "RMSE all | RMSE PYR2-PYR4 | RMSE PYR5 |\n"
        )
        handle.write("|---:|---:|---:|---:|---:|\n")

        for row in classical_rows:
            handle.write(
                f"| {float(row['gamma_phi_ps_inv']):g} "
                f"| {float(row['comparison_start_time_ps']):.4f} "
                f"| {float(row['population_rmse_all']):.6e} "
                f"| {float(row['population_rmse_high_manifold']):.6e} "
                f"| {float(row['population_rmse_PYR5_from_high_states']):.6e} "
                "|\n"
            )

    log("Day020 high-dephasing Zeno-scaling audit completed.")
    log(f"Hamiltonians: {len(hamiltonian_files)}/{EXPECTED_FRAMES}")
    log(
        "Maximum Hamiltonian symmetry error: "
        f"{symmetry_error:.3e} eV"
    )
    log(
        "Maximum classical norm error: "
        f"{maximum_classical_norm_error:.3e}"
    )
    log(
        "Mean predicted PYR2-PYR4 turnover scale: "
        f"{high_turnover_mean:.3f} ps^-1"
    )
    log(
        "Mean predicted PYR5-pair turnover scale: "
        f"{pyr5_turnover_mean:.3f} ps^-1"
    )
    log(
        "Mean high-manifold k(100)/k(200): "
        f"{high_k_ratio_mean:.6f}"
    )
    log(
        "Mean PYR5-pair k(100)/k(200): "
        f"{pyr5_k_ratio_mean:.6f}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
