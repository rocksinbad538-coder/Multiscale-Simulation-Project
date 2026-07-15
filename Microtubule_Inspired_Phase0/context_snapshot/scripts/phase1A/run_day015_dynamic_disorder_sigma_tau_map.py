from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_dynamic_disorder_sigma_tau_map")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514
KB = 8.617333262145e-5
T = 300.0
kBT = KB*T

H0 = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
], float)

E, U = np.linalg.eigh(H0)

k_down_ps = 1.0
gamma_phi_ps = 1.0
k_down_fs = k_down_ps / 1000.0
gamma_phi_fs = gamma_phi_ps / 1000.0

rng = np.random.default_rng(20260624)

sigma_E_values = [0, 2, 5, 10, 20]
tau_c_values = [0.02, 0.05, 0.10, 0.20, 0.50, 1.00]

n_traj = 8
dt_fs = 5.0
t_end_fs = 3000.0
time = np.arange(0, t_end_fs + dt_fs, dt_fs)

rho0 = np.zeros((4,4), complex)
rho0[3,3] = 1.0

def dissipator(rho):
    rho_e = U.conj().T @ rho @ U
    dr = np.zeros_like(rho_e)

    for high in range(4):
        for low in range(4):
            if E[high] <= E[low]:
                continue
            dE = E[high] - E[low]
            kup = k_down_fs * np.exp(-dE/kBT)

            Ld = np.zeros((4,4), complex)
            Ld[low, high] = np.sqrt(k_down_fs)
            dr += Ld @ rho_e @ Ld.conj().T
            dr -= 0.5*(Ld.conj().T@Ld@rho_e + rho_e@Ld.conj().T@Ld)

            Lu = np.zeros((4,4), complex)
            Lu[high, low] = np.sqrt(kup)
            dr += Lu @ rho_e @ Lu.conj().T
            dr -= 0.5*(Lu.conj().T@Lu@rho_e + rho_e@Lu.conj().T@Lu)

    for i in range(4):
        L = np.zeros((4,4), complex)
        L[i,i] = np.sqrt(gamma_phi_fs)
        dr += L @ rho_e @ L.conj().T
        dr -= 0.5*(L.conj().T@L@rho_e + rho_e@L.conj().T@L)

    return U @ dr @ U.conj().T

def ou_series(sigma_eV, tau_fs):
    x = np.zeros((len(time), 4))
    if sigma_eV == 0:
        return x
    a = np.exp(-dt_fs/tau_fs)
    b = sigma_eV*np.sqrt(1-a*a)
    for k in range(1, len(time)):
        x[k] = a*x[k-1] + b*rng.normal(size=4)
    return x

def rhs(rho, H):
    return -1j/HBAR*(H@rho - rho@H) + dissipator(rho)

def rk4_step(rho, H):
    k1 = rhs(rho, H)
    k2 = rhs(rho + 0.5*dt_fs*k1, H)
    k3 = rhs(rho + 0.5*dt_fs*k2, H)
    k4 = rhs(rho + dt_fs*k3, H)
    return rho + dt_fs*(k1 + 2*k2 + 2*k3 + k4)/6.0

def first_cross(vals, thr):
    idx = np.where(vals >= thr)[0]
    return np.nan if len(idx) == 0 else float(time[idx[0]])

rows = []

for sigma in sigma_E_values:
    for tau in tau_c_values:
        vals = []

        for tr in range(n_traj):
            noise = ou_series(sigma/1000.0, tau*1000.0)
            rho = rho0.copy()
            pyr2 = []

            min_eig = 1.0
            max_trace_err = 0.0

            for k, t in enumerate(time):
                rh = (rho + rho.conj().T)/2
                pyr2.append(np.real(rh[0,0]))
                min_eig = min(min_eig, float(np.linalg.eigvalsh(rh).min()))
                max_trace_err = max(max_trace_err, abs(np.trace(rho)-1.0))

                H = H0.copy()
                for i in range(4):
                    H[i,i] += noise[k,i]

                rho = rk4_step(rho, H)

            pyr2 = np.array(pyr2)

            vals.append({
                "max_PYR2": pyr2.max(),
                "final_PYR2": pyr2[-1],
                "t10_fs": first_cross(pyr2, 0.10),
                "t20_fs": first_cross(pyr2, 0.20),
                "t30_fs": first_cross(pyr2, 0.30),
                "min_eig": min_eig,
                "max_trace_err": max_trace_err,
            })

        df = pd.DataFrame(vals)

        rows.append({
            "sigma_E_meV": sigma,
            "tau_c_ps": tau,
            "n_traj": n_traj,
            "mean_max_PYR2": df["max_PYR2"].mean(),
            "std_max_PYR2": df["max_PYR2"].std(),
            "mean_final_PYR2": df["final_PYR2"].mean(),
            "mean_t10_fs": df["t10_fs"].mean(),
            "mean_t20_fs": df["t20_fs"].mean(),
            "mean_t30_fs": df["t30_fs"].mean(),
            "min_eig_min": df["min_eig"].min(),
            "max_trace_err": df["max_trace_err"].max(),
        })

        print(f"done sigma={sigma} meV tau={tau} ps")

summary = pd.DataFrame(rows)
summary.to_csv(OUT/"sigma_tau_dynamic_disorder_summary.csv", index=False)

for metric in ["mean_max_PYR2", "mean_t30_fs"]:
    pivot = summary.pivot(index="sigma_E_meV", columns="tau_c_ps", values=metric)
    plt.figure(figsize=(7,4))
    plt.imshow(pivot.values, aspect="auto", origin="lower")
    plt.xticks(range(len(pivot.columns)), pivot.columns)
    plt.yticks(range(len(pivot.index)), pivot.index)
    plt.xlabel("tau_c (ps)")
    plt.ylabel("sigma_E (meV)")
    plt.colorbar(label=metric)
    plt.tight_layout()
    plt.savefig(OUT/f"heatmap_{metric}.png", dpi=300)
    plt.close()

md = """# Day015 Dynamic-Disorder Sigma-Tau Map

## Purpose

Map how stochastic site-energy disorder amplitude and correlation time affect access to PYR2.

## Model

- Static hydrated four-site Hamiltonian.
- Site-energy Ornstein-Uhlenbeck fluctuations.
- Fixed Lindblad exciton-relaxation bath.
- Initial excitation on PYR5.

## Caveat

This is a quick exploratory map using synthetic disorder. It is intended to locate enhancement/suppression regimes before replacing disorder parameters with MD-derived fluctuations.
"""
(OUT/"SIGMA_TAU_DYNAMIC_DISORDER_MAP_DAY015.md").write_text(md)

print(summary.to_string(index=False))
print("Wrote:", OUT)
