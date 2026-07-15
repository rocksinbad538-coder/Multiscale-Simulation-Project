from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_dynamic_disorder_lindblad")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514  # eV fs
KB = 8.617333262145e-5
T = 300.0
kBT = KB * T

labels = ["PYR2", "PYR3", "PYR4", "PYR5"]

H0 = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
], dtype=float)

# Dissipator basis fixed from static Hamiltonian.
E, U = np.linalg.eigh(H0)

k_down_ps = 1.0
gamma_phi_ps = 1.0
k_down_fs = k_down_ps / 1000.0
gamma_phi_fs = gamma_phi_ps / 1000.0

rng = np.random.default_rng(20260624)

t_end_fs = 5000.0
dt_noise_fs = 5.0
time_grid_fs = np.arange(0.0, t_end_fs + dt_noise_fs, dt_noise_fs)

n_traj = 30

cases = [
    {
        "case": "static",
        "sigma_E_meV": 0.0,
        "sigma_J_meV": 0.0,
        "tau_c_ps": 0.10,
    },
    {
        "case": "site_disorder_5meV",
        "sigma_E_meV": 5.0,
        "sigma_J_meV": 0.0,
        "tau_c_ps": 0.10,
    },
    {
        "case": "site_disorder_10meV",
        "sigma_E_meV": 10.0,
        "sigma_J_meV": 0.0,
        "tau_c_ps": 0.10,
    },
    {
        "case": "coupling_disorder_2meV",
        "sigma_E_meV": 0.0,
        "sigma_J_meV": 2.0,
        "tau_c_ps": 0.10,
    },
    {
        "case": "site10_coupling2meV",
        "sigma_E_meV": 10.0,
        "sigma_J_meV": 2.0,
        "tau_c_ps": 0.10,
    },
]

def ou_noise(n_steps, sigma_eV, tau_c_fs):
    if sigma_eV == 0.0:
        return np.zeros(n_steps)

    x = np.zeros(n_steps)
    a = np.exp(-dt_noise_fs / tau_c_fs)
    b = sigma_eV * np.sqrt(1.0 - a*a)

    for i in range(1, n_steps):
        x[i] = a * x[i-1] + b * rng.normal()

    return x

def build_noise(case):
    n = len(time_grid_fs)
    tau_c_fs = case["tau_c_ps"] * 1000.0
    sigma_E = case["sigma_E_meV"] / 1000.0
    sigma_J = case["sigma_J_meV"] / 1000.0

    noise_E = np.vstack([
        ou_noise(n, sigma_E, tau_c_fs) for _ in range(4)
    ]).T

    noise_J = np.vstack([
        ou_noise(n, sigma_J, tau_c_fs) for _ in range(3)
    ]).T

    return noise_E, noise_J

def H_at_time(t_fs, noise_E, noise_J):
    idx = int(np.clip(round(t_fs / dt_noise_fs), 0, len(time_grid_fs)-1))
    H = H0.copy()

    for i in range(4):
        H[i, i] += noise_E[idx, i]

    H[0, 1] += noise_J[idx, 0]
    H[1, 0] += noise_J[idx, 0]

    H[1, 2] += noise_J[idx, 1]
    H[2, 1] += noise_J[idx, 1]

    H[2, 3] += noise_J[idx, 2]
    H[3, 2] += noise_J[idx, 2]

    return H

def dissipator_exciton_basis(rho_site):
    rho_exc = U.conj().T @ rho_site @ U
    drho_exc = np.zeros_like(rho_exc, dtype=complex)

    for high in range(4):
        for low in range(4):
            if E[high] <= E[low]:
                continue

            dE = E[high] - E[low]
            k_up_fs = k_down_fs * np.exp(-dE / kBT)

            Ld = np.zeros((4,4), dtype=complex)
            Ld[low, high] = np.sqrt(k_down_fs)
            drho_exc += Ld @ rho_exc @ Ld.conj().T
            drho_exc -= 0.5 * (Ld.conj().T @ Ld @ rho_exc + rho_exc @ Ld.conj().T @ Ld)

            Lu = np.zeros((4,4), dtype=complex)
            Lu[high, low] = np.sqrt(k_up_fs)
            drho_exc += Lu @ rho_exc @ Lu.conj().T
            drho_exc -= 0.5 * (Lu.conj().T @ Lu @ rho_exc + rho_exc @ Lu.conj().T @ Lu)

    for i in range(4):
        L = np.zeros((4,4), dtype=complex)
        L[i, i] = np.sqrt(gamma_phi_fs)
        drho_exc += L @ rho_exc @ L.conj().T
        drho_exc -= 0.5 * (L.conj().T @ L @ rho_exc + rho_exc @ L.conj().T @ L)

    return U @ drho_exc @ U.conj().T

def rhs(t, y, noise_E, noise_J):
    rho = y.reshape((4,4)).astype(complex)
    H = H_at_time(t, noise_E, noise_J)

    drho = -1j / HBAR * (H @ rho - rho @ H)
    drho += dissipator_exciton_basis(rho)

    return drho.reshape(-1)

def first_crossing(time_fs, vals, thr):
    idx = np.where(vals >= thr)[0]
    if len(idx) == 0:
        return np.nan
    return float(time_fs[idx[0]])

all_summary = []

rho0 = np.zeros((4,4), dtype=complex)
rho0[3,3] = 1.0

t_eval = np.linspace(0, t_end_fs, 501)

for case in cases:
    traj_rows = []
    pyr2_curves = []

    for tr in range(n_traj):
        noise_E, noise_J = build_noise(case)

        sol = solve_ivp(
            rhs,
            (0, t_end_fs),
            rho0.reshape(-1),
            t_eval=t_eval,
            args=(noise_E, noise_J),
            rtol=1e-7,
            atol=1e-9,
        )

        pops = []
        min_eig = np.inf
        max_trace_error = 0.0

        for k, t in enumerate(sol.t):
            rho = sol.y[:, k].reshape((4,4))
            rho_h = (rho + rho.conj().T) / 2.0
            eigs = np.linalg.eigvalsh(rho_h)
            min_eig = min(min_eig, float(np.min(eigs)))
            max_trace_error = max(max_trace_error, abs(np.trace(rho) - 1.0))
            pops.append(np.real(np.diag(rho_h)))

        pops = np.array(pops)
        pyr2 = pops[:,0]
        pyr2_curves.append(pyr2)

        traj_rows.append({
            "case": case["case"],
            "trajectory": tr,
            "solver_success": bool(sol.success),
            "max_trace_error": max_trace_error,
            "min_density_eigenvalue": min_eig,
            "max_PYR2": float(np.max(pyr2)),
            "final_PYR2": float(pyr2[-1]),
            "t_PYR2_10pct_fs": first_crossing(sol.t, pyr2, 0.10),
            "t_PYR2_20pct_fs": first_crossing(sol.t, pyr2, 0.20),
            "t_PYR2_30pct_fs": first_crossing(sol.t, pyr2, 0.30),
        })

    traj_df = pd.DataFrame(traj_rows)
    traj_df.to_csv(OUT / f"{case['case']}_trajectory_summary.csv", index=False)

    pyr2_curves = np.array(pyr2_curves)
    mean = pyr2_curves.mean(axis=0)
    lo = np.percentile(pyr2_curves, 10, axis=0)
    hi = np.percentile(pyr2_curves, 90, axis=0)

    curve_df = pd.DataFrame({
        "time_fs": t_eval,
        "PYR2_mean": mean,
        "PYR2_p10": lo,
        "PYR2_p90": hi,
    })
    curve_df.to_csv(OUT / f"{case['case']}_PYR2_ensemble_curve.csv", index=False)

    plt.figure(figsize=(6,4))
    plt.plot(t_eval/1000.0, mean, label="mean PYR2")
    plt.fill_between(t_eval/1000.0, lo, hi, alpha=0.25, label="10-90%")
    plt.xlabel("time (ps)")
    plt.ylabel("PYR2 population")
    plt.title(case["case"])
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"{case['case']}_PYR2_ensemble.png", dpi=300)
    plt.close()

    all_summary.append({
        "case": case["case"],
        "sigma_E_meV": case["sigma_E_meV"],
        "sigma_J_meV": case["sigma_J_meV"],
        "tau_c_ps": case["tau_c_ps"],
        "n_traj": n_traj,
        "mean_max_PYR2": traj_df["max_PYR2"].mean(),
        "std_max_PYR2": traj_df["max_PYR2"].std(),
        "mean_final_PYR2": traj_df["final_PYR2"].mean(),
        "std_final_PYR2": traj_df["final_PYR2"].std(),
        "mean_t_PYR2_10pct_fs": traj_df["t_PYR2_10pct_fs"].mean(),
        "mean_t_PYR2_20pct_fs": traj_df["t_PYR2_20pct_fs"].mean(),
        "mean_t_PYR2_30pct_fs": traj_df["t_PYR2_30pct_fs"].mean(),
        "min_density_eigenvalue_min": traj_df["min_density_eigenvalue"].min(),
        "max_trace_error_max": traj_df["max_trace_error"].max(),
    })

summary = pd.DataFrame(all_summary)
summary.to_csv(OUT / "dynamic_disorder_lindblad_summary.csv", index=False)

md = """# Day015 Dynamic-Disorder Lindblad Scan

## Purpose

Evaluate whether dynamic fluctuations in site energies and couplings modify PYR2 access in the hydrated four-site excitonic model.

## Model

The static hydrated Hamiltonian is perturbed by Ornstein-Uhlenbeck stochastic fluctuations in either site energies, couplings, or both. Lindblad relaxation is kept fixed using the previously validated exciton-basis thermal model.

## Caveat

This is not yet a production MD-derived disorder model. It is a controlled synthetic-disorder stress test used to evaluate model robustness and expected sensitivity.
"""
(OUT / "DYNAMIC_DISORDER_LINDBLAD_DAY015.md").write_text(md)

print(summary.to_string(index=False))
print("Wrote:", OUT)
