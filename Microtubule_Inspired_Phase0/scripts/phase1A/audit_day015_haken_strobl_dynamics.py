from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

OUT = Path("runs/phase1A/day015_exciton_model")
OUT.mkdir(parents=True, exist_ok=True)

HBAR = 0.6582119514  # eV fs

H = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
], dtype=float)

labels = ["PYR2", "PYR3", "PYR4", "PYR5"]
gammas_ps = [0, 1, 5, 10, 50]


def rhs(t, y, gamma_fs):
    n = H.shape[0]
    rho = y.reshape((n, n)).astype(complex)

    drho = -1j / HBAR * (H @ rho - rho @ H)

    # Site-basis pure dephasing.
    if gamma_fs > 0:
        for i in range(n):
            for j in range(n):
                if i != j:
                    drho[i, j] -= gamma_fs * rho[i, j]

    return drho.reshape(-1)


def run_case(gamma_ps):
    gamma_fs = gamma_ps / 1000.0

    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[3, 3] = 1.0  # initial excitation on PYR5

    t_eval = np.linspace(0, 10000, 1001)

    sol = solve_ivp(
        rhs,
        (0, 10000),
        rho0.reshape(-1),
        t_eval=t_eval,
        args=(gamma_fs,),
        rtol=1e-9,
        atol=1e-11,
    )

    records = []
    min_eig_global = +np.inf
    max_trace_error = 0.0
    max_hermiticity_error = 0.0
    max_population_error = 0.0

    for k, t in enumerate(sol.t):
        rho = sol.y[:, k].reshape((4, 4))

        trace = np.trace(rho)
        herm_err = np.max(np.abs(rho - rho.conj().T))
        evals = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
        min_eig = np.min(evals)

        pops = np.real(np.diag(rho))
        pop_sum = np.sum(pops)
        pop_err = abs(pop_sum - 1.0)

        coherence_l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
        purity = np.real(np.trace(rho @ rho))

        min_eig_global = min(min_eig_global, min_eig)
        max_trace_error = max(max_trace_error, abs(trace - 1.0))
        max_hermiticity_error = max(max_hermiticity_error, herm_err)
        max_population_error = max(max_population_error, pop_err)

        records.append({
            "time_fs": t,
            "trace_real": np.real(trace),
            "trace_imag": np.imag(trace),
            "min_density_eigenvalue": min_eig,
            "hermiticity_error": herm_err,
            "population_sum": pop_sum,
            "population_error": pop_err,
            "coherence_l1": coherence_l1,
            "purity": purity,
            **{labels[i]: pops[i] for i in range(4)}
        })

    df = pd.DataFrame(records)
    df.to_csv(OUT / f"haken_strobl_audit_gamma{gamma_ps}ps.csv", index=False)

    def first_crossing(site, thr):
        hit = df[df[site] >= thr]
        if len(hit) == 0:
            return np.nan
        return float(hit.iloc[0]["time_fs"])

    summary = {
        "gamma_ps": gamma_ps,
        "solver_success": bool(sol.success),
        "max_trace_error": max_trace_error,
        "max_hermiticity_error": max_hermiticity_error,
        "max_population_error": max_population_error,
        "min_density_eigenvalue": min_eig_global,
        "final_PYR2": float(df.iloc[-1]["PYR2"]),
        "final_PYR3": float(df.iloc[-1]["PYR3"]),
        "final_PYR4": float(df.iloc[-1]["PYR4"]),
        "final_PYR5": float(df.iloc[-1]["PYR5"]),
        "final_coherence_l1": float(df.iloc[-1]["coherence_l1"]),
        "final_purity": float(df.iloc[-1]["purity"]),
        "t_PYR2_1pct_fs": first_crossing("PYR2", 0.01),
        "t_PYR2_5pct_fs": first_crossing("PYR2", 0.05),
        "t_PYR2_10pct_fs": first_crossing("PYR2", 0.10),
        "max_PYR2": float(df["PYR2"].max()),
        "time_max_PYR2_fs": float(df.loc[df["PYR2"].idxmax(), "time_fs"]),
    }

    return summary


summaries = [run_case(g) for g in gammas_ps]
summary_df = pd.DataFrame(summaries)
summary_df.to_csv(OUT / "haken_strobl_audit_summary.csv", index=False)

md = ["# Day015 Haken-Strobl Dynamics Audit", ""]
md.append("## Purpose")
md.append("")
md.append("Audit trace conservation, Hermiticity, positivity, population conservation, coherence decay, and early PYR2 arrival times for the hydrated four-site Hamiltonian.")
md.append("")
md.append("## Interpretation rule")
md.append("")
md.append("The long-time uniform population reached by the pure-dephasing Haken-Strobl model should be treated as a model artifact, not as thermodynamic relaxation. Early-time transfer metrics are the physically useful diagnostic at this stage.")
md.append("")
md.append("## Summary table")
md.append("")
def df_to_markdown(df):
    cols = list(df.columns)
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                vals.append(f"{v:.6g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

md.append(df_to_markdown(summary_df))
md.append("")
md.append("## Files")
md.append("")
md.append("- `haken_strobl_audit_summary.csv`")
md.append("- `haken_strobl_audit_gamma0ps.csv`")
md.append("- `haken_strobl_audit_gamma1ps.csv`")
md.append("- `haken_strobl_audit_gamma5ps.csv`")
md.append("- `haken_strobl_audit_gamma10ps.csv`")
md.append("- `haken_strobl_audit_gamma50ps.csv`")

(OUT / "HAKEN_STROBL_AUDIT_DAY015.md").write_text("\n".join(md))

print(summary_df.to_string(index=False))
print("Wrote:", OUT / "haken_strobl_audit_summary.csv")
print("Wrote:", OUT / "HAKEN_STROBL_AUDIT_DAY015.md")
