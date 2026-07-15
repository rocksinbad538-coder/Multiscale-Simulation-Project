from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_lindblad_relaxation")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514  # eV fs
KB = 8.617333262145e-5  # eV/K
T = 300.0
kBT = KB * T

labels_site = ["PYR2", "PYR3", "PYR4", "PYR5"]

H_site = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
], dtype=float)

E, U = np.linalg.eigh(H_site)
# columns of U are exciton eigenvectors in site basis

labels_exc = [f"X{i+1}" for i in range(4)]

# relaxation rates in ps^-1.
# test a physically modest range.
k_down_values_ps = [0.1, 1.0, 5.0, 10.0]

# optional pure dephasing in exciton basis, ps^-1.
gamma_phi_ps = 1.0


def lindblad_rhs(t, y, k_down_ps):
    rho_site = y.reshape((4, 4)).astype(complex)

    # Hamiltonian unitary evolution in site basis
    drho = -1j / HBAR * (H_site @ rho_site - rho_site @ H_site)

    # transform rho to exciton basis for dissipators
    rho_exc = U.conj().T @ rho_site @ U
    drho_exc = np.zeros_like(rho_exc, dtype=complex)

    # population relaxation between exciton eigenstates
    for high in range(4):
        for low in range(4):
            if E[high] <= E[low]:
                continue

            dE = E[high] - E[low]

            k_down_fs = k_down_ps / 1000.0
            k_up_fs = k_down_fs * np.exp(-dE / kBT)

            # high -> low
            L = np.zeros((4, 4), dtype=complex)
            L[low, high] = np.sqrt(k_down_fs)
            drho_exc += L @ rho_exc @ L.conj().T
            drho_exc -= 0.5 * (L.conj().T @ L @ rho_exc + rho_exc @ L.conj().T @ L)

            # low -> high
            Lup = np.zeros((4, 4), dtype=complex)
            Lup[high, low] = np.sqrt(k_up_fs)
            drho_exc += Lup @ rho_exc @ Lup.conj().T
            drho_exc -= 0.5 * (Lup.conj().T @ Lup @ rho_exc + rho_exc @ Lup.conj().T @ Lup)

    # pure dephasing in exciton basis
    gamma_phi_fs = gamma_phi_ps / 1000.0
    for i in range(4):
        L = np.zeros((4, 4), dtype=complex)
        L[i, i] = np.sqrt(gamma_phi_fs)
        drho_exc += L @ rho_exc @ L.conj().T
        drho_exc -= 0.5 * (L.conj().T @ L @ rho_exc + rho_exc @ L.conj().T @ L)

    # transform dissipator back to site basis and add to unitary part
    drho += U @ drho_exc @ U.conj().T

    return drho.reshape(-1)


def first_crossing(df, site, thr):
    hit = df[df[site] >= thr]
    if len(hit) == 0:
        return np.nan
    return float(hit.iloc[0]["time_fs"])


summary = []

# initial site excitation on PYR5
rho0 = np.zeros((4, 4), dtype=complex)
rho0[3, 3] = 1.0

t_eval = np.linspace(0, 10000, 1001)

for k_down_ps in k_down_values_ps:
    sol = solve_ivp(
        lindblad_rhs,
        (0, 10000),
        rho0.reshape(-1),
        t_eval=t_eval,
        args=(k_down_ps,),
        rtol=1e-9,
        atol=1e-11,
    )

    records = []
    min_eig = np.inf
    max_trace_error = 0.0
    max_herm_error = 0.0

    for idx, t in enumerate(sol.t):
        rho = sol.y[:, idx].reshape((4, 4))
        rho_h = (rho + rho.conj().T) / 2
        eigs = np.linalg.eigvalsh(rho_h)

        trace = np.trace(rho)
        max_trace_error = max(max_trace_error, abs(trace - 1.0))
        max_herm_error = max(max_herm_error, np.max(np.abs(rho - rho.conj().T)))
        min_eig = min(min_eig, float(np.min(eigs)))

        pops_site = np.real(np.diag(rho_h))
        rho_exc = U.conj().T @ rho_h @ U
        pops_exc = np.real(np.diag(rho_exc))

        records.append({
            "time_fs": t,
            **{labels_site[i]: pops_site[i] for i in range(4)},
            **{labels_exc[i]: pops_exc[i] for i in range(4)},
            "purity": float(np.real(np.trace(rho_h @ rho_h))),
            "trace_real": float(np.real(trace)),
            "trace_imag": float(np.imag(trace)),
        })

    df = pd.DataFrame(records)
    tag = f"kdown_{k_down_ps:g}ps".replace(".", "p")
    df.to_csv(OUT / f"lindblad_relaxation_{tag}.csv", index=False)

    plt.figure(figsize=(6,4))
    for lab in labels_site:
        plt.plot(df["time_fs"] / 1000.0, df[lab], label=lab)
    plt.xlabel("time (ps)")
    plt.ylabel("site population")
    plt.title(f"Exciton relaxation, k_down={k_down_ps:g} ps^-1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"lindblad_site_populations_{tag}.png", dpi=300)
    plt.close()

    plt.figure(figsize=(6,4))
    for lab in labels_exc:
        plt.plot(df["time_fs"] / 1000.0, df[lab], label=lab)
    plt.xlabel("time (ps)")
    plt.ylabel("exciton population")
    plt.title(f"Exciton-basis populations, k_down={k_down_ps:g} ps^-1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / f"lindblad_exciton_populations_{tag}.png", dpi=300)
    plt.close()

    summary.append({
        "k_down_ps": k_down_ps,
        "gamma_phi_ps": gamma_phi_ps,
        "solver_success": bool(sol.success),
        "max_trace_error": max_trace_error,
        "max_hermiticity_error": max_herm_error,
        "min_density_eigenvalue": min_eig,
        "final_PYR2": float(df.iloc[-1]["PYR2"]),
        "final_PYR3": float(df.iloc[-1]["PYR3"]),
        "final_PYR4": float(df.iloc[-1]["PYR4"]),
        "final_PYR5": float(df.iloc[-1]["PYR5"]),
        "max_PYR2": float(df["PYR2"].max()),
        "time_max_PYR2_fs": float(df.loc[df["PYR2"].idxmax(), "time_fs"]),
        "t_PYR2_1pct_fs": first_crossing(df, "PYR2", 0.01),
        "t_PYR2_5pct_fs": first_crossing(df, "PYR2", 0.05),
        "t_PYR2_10pct_fs": first_crossing(df, "PYR2", 0.10),
        "final_purity": float(df.iloc[-1]["purity"]),
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / "lindblad_relaxation_summary.csv", index=False)

pd.DataFrame(H_site, index=labels_site, columns=labels_site).to_csv(OUT / "lindblad_input_hamiltonian_site_basis.csv")
pd.DataFrame({"exciton": labels_exc, "energy_eV": E}).to_csv(OUT / "lindblad_exciton_energies.csv", index=False)
pd.DataFrame(U, index=labels_site, columns=labels_exc).to_csv(OUT / "lindblad_exciton_eigenvectors.csv")

print(summary_df.to_string(index=False))
print("Exciton energies:")
for lab, e in zip(labels_exc, E):
    print(lab, f"{e:.6f} eV")
print("Wrote:", OUT)
