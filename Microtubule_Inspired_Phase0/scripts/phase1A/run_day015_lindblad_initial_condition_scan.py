from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_lindblad_initial_conditions")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514
KB = 8.617333262145e-5
T = 300.0
kBT = KB * T

labels_site = ["PYR2", "PYR3", "PYR4", "PYR5"]
labels_exc = ["X1", "X2", "X3", "X4"]

H_site = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
], dtype=float)

E, U = np.linalg.eigh(H_site)

k_down_ps = 1.0
gamma_phi_ps = 1.0
k_down_fs = k_down_ps / 1000.0
gamma_phi_fs = gamma_phi_ps / 1000.0

def rhs(t, y):
    rho_site = y.reshape((4, 4)).astype(complex)

    drho = -1j / HBAR * (H_site @ rho_site - rho_site @ H_site)

    rho_exc = U.conj().T @ rho_site @ U
    drho_exc = np.zeros_like(rho_exc)

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

    drho += U @ drho_exc @ U.conj().T
    return drho.reshape(-1)

def first_crossing(df, site, thr):
    hit = df[df[site] >= thr]
    if len(hit) == 0:
        return np.nan
    return float(hit.iloc[0]["time_fs"])

summary = []
t_eval = np.linspace(0, 10000, 1001)

for init_idx, init_site in enumerate(labels_site):
    rho0 = np.zeros((4,4), dtype=complex)
    rho0[init_idx, init_idx] = 1.0

    sol = solve_ivp(
        rhs,
        (0, 10000),
        rho0.reshape(-1),
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-11,
    )

    records = []
    min_eig = np.inf
    max_trace_error = 0.0

    for k, t in enumerate(sol.t):
        rho = sol.y[:, k].reshape((4,4))
        rho_h = (rho + rho.conj().T) / 2
        eigs = np.linalg.eigvalsh(rho_h)
        min_eig = min(min_eig, float(np.min(eigs)))
        max_trace_error = max(max_trace_error, abs(np.trace(rho) - 1.0))
        pops_site = np.real(np.diag(rho_h))
        rho_exc = U.conj().T @ rho_h @ U
        pops_exc = np.real(np.diag(rho_exc))

        records.append({
            "time_fs": t,
            **{labels_site[i]: pops_site[i] for i in range(4)},
            **{labels_exc[i]: pops_exc[i] for i in range(4)},
        })

    df = pd.DataFrame(records)
    tag = init_site
    df.to_csv(OUT / f"lindblad_init_{tag}.csv", index=False)

    plt.figure(figsize=(6,4))
    for lab in labels_site:
        plt.plot(df["time_fs"]/1000.0, df[lab], label=lab)
    plt.xlabel("time (ps)")
    plt.ylabel("site population")
    plt.title(f"Initial excitation: {init_site}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"lindblad_site_populations_init_{tag}.png", dpi=300)
    plt.close()

    summary.append({
        "initial_site": init_site,
        "solver_success": bool(sol.success),
        "max_trace_error": max_trace_error,
        "min_density_eigenvalue": min_eig,
        "final_PYR2": float(df.iloc[-1]["PYR2"]),
        "final_PYR3": float(df.iloc[-1]["PYR3"]),
        "final_PYR4": float(df.iloc[-1]["PYR4"]),
        "final_PYR5": float(df.iloc[-1]["PYR5"]),
        "max_PYR2": float(df["PYR2"].max()),
        "t_PYR2_10pct_fs": first_crossing(df, "PYR2", 0.10),
        "t_PYR2_20pct_fs": first_crossing(df, "PYR2", 0.20),
        "t_PYR2_30pct_fs": first_crossing(df, "PYR2", 0.30),
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / "lindblad_initial_condition_scan_summary.csv", index=False)

print(summary_df.to_string(index=False))
print("Wrote:", OUT)
