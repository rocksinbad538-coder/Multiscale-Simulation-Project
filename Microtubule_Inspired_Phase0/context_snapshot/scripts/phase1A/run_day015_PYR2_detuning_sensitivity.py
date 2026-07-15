from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_PYR2_detuning_sensitivity")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514  # eV fs
labels = ["PYR2", "PYR3", "PYR4", "PYR5"]

E_ref = 3.77750
J23 = 0.00237
J34 = 0.01700
J45 = 0.03300

gamma_ps = 10.0
gamma_fs = gamma_ps / 1000.0

# E_PYR2 = E_ref - detuning
detuning_meV_values = [0, 5, 10, 15, 20, 22.5, 30, 40]


def build_H(detuning_meV):
    E2 = E_ref - detuning_meV / 1000.0
    H = np.array([
        [E2, J23, 0.0, 0.0],
        [J23, E_ref, J34, 0.0],
        [0.0, J34, E_ref, J45],
        [0.0, 0.0, J45, E_ref],
    ], dtype=float)
    return H


def rhs(t, y, H):
    rho = y.reshape((4, 4)).astype(complex)
    drho = -1j / HBAR * (H @ rho - rho @ H)

    for i in range(4):
        for j in range(4):
            if i != j:
                drho[i, j] -= gamma_fs * rho[i, j]

    return drho.reshape(-1)


def first_crossing(df, site, thr):
    hit = df[df[site] >= thr]
    if len(hit) == 0:
        return np.nan
    return float(hit.iloc[0]["time_fs"])


summary = []

for detuning_meV in detuning_meV_values:
    H = build_H(detuning_meV)

    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[3, 3] = 1.0

    t_eval = np.linspace(0, 10000, 1001)

    sol = solve_ivp(
        rhs,
        (0, 10000),
        rho0.reshape(-1),
        t_eval=t_eval,
        args=(H,),
        rtol=1e-9,
        atol=1e-11,
    )

    pops = []
    min_eig = np.inf
    max_trace_error = 0.0

    for k in range(len(sol.t)):
        rho = sol.y[:, k].reshape((4, 4))
        rho_h = (rho + rho.conj().T) / 2.0
        eig = np.linalg.eigvalsh(rho_h)
        min_eig = min(min_eig, float(np.min(eig)))
        max_trace_error = max(max_trace_error, abs(np.trace(rho) - 1.0))
        pops.append(np.real(np.diag(rho)))

    pops = np.array(pops)

    df = pd.DataFrame({
        "time_fs": sol.t,
        "PYR2": pops[:, 0],
        "PYR3": pops[:, 1],
        "PYR4": pops[:, 2],
        "PYR5": pops[:, 3],
    })

    tag = f"detuning_{detuning_meV:g}meV".replace(".", "p")
    df.to_csv(OUT / f"haken_strobl_gamma10ps_{tag}.csv", index=False)

    summary.append({
        "detuning_meV": detuning_meV,
        "E_PYR2_eV": E_ref - detuning_meV / 1000.0,
        "gamma_ps": gamma_ps,
        "solver_success": bool(sol.success),
        "max_trace_error": max_trace_error,
        "min_density_eigenvalue": min_eig,
        "max_PYR2": float(df["PYR2"].max()),
        "time_max_PYR2_fs": float(df.loc[df["PYR2"].idxmax(), "time_fs"]),
        "final_PYR2": float(df.iloc[-1]["PYR2"]),
        "t_PYR2_1pct_fs": first_crossing(df, "PYR2", 0.01),
        "t_PYR2_5pct_fs": first_crossing(df, "PYR2", 0.05),
        "t_PYR2_10pct_fs": first_crossing(df, "PYR2", 0.10),
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / "PYR2_detuning_sensitivity_summary.csv", index=False)

plt.figure(figsize=(6,4))
plt.plot(summary_df["detuning_meV"], summary_df["max_PYR2"], "o-")
plt.xlabel("PYR2 detuning below PYR3/PYR4/PYR5 (meV)")
plt.ylabel("Maximum PYR2 population")
plt.tight_layout()
plt.savefig(OUT / "PYR2_detuning_sensitivity_max_PYR2.png", dpi=300)
plt.close()

plt.figure(figsize=(6,4))
for col in ["t_PYR2_1pct_fs", "t_PYR2_5pct_fs", "t_PYR2_10pct_fs"]:
    plt.plot(summary_df["detuning_meV"], summary_df[col], "o-", label=col)
plt.xlabel("PYR2 detuning below PYR3/PYR4/PYR5 (meV)")
plt.ylabel("First arrival time (fs)")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "PYR2_detuning_sensitivity_arrival_times.png", dpi=300)
plt.close()

print(summary_df.to_string(index=False))
print("Wrote:", OUT)
