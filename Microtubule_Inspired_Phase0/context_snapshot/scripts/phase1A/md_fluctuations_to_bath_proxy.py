from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HBAR_eV_ps = 0.0006582119514  # eV ps


def autocorr(x):
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    n = len(x)
    f = np.fft.fft(x, n=2*n)
    ac = np.fft.ifft(f * np.conjugate(f))[:n].real
    norm = np.arange(n, 0, -1)
    ac = ac / norm
    if ac[0] != 0:
        ac_norm = ac / ac[0]
    else:
        ac_norm = ac
    return ac, ac_norm


def correlation_time(time_ps, ac_norm):
    idx = len(ac_norm)

    for i in range(1, len(ac_norm)):
        if ac_norm[i] <= 0:
            idx = i
            break

    return np.trapezoid(
        ac_norm[:idx],
        x=time_ps[:idx]
    )


def spectral_density_proxy(time_ps, ac):
    dt = time_ps[1] - time_ps[0]
    freq = np.fft.rfftfreq(len(ac), d=dt)  # ps^-1
    spec = np.real(np.fft.rfft(ac))
    spec = np.maximum(spec, 0.0)
    return freq, spec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CSV with time_ps and fluctuation variables")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--columns", nargs="+", required=True)
    args = ap.parse_args()

    inp = Path(args.input)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inp)
    if "time_ps" not in df.columns:
        raise ValueError("Input CSV must contain time_ps column.")

    time = df["time_ps"].values
    if len(time) < 3:
        raise ValueError("Need at least 3 time points.")

    dt = np.mean(np.diff(time))

    rows = []

    for col in args.columns:
        if col not in df.columns:
            raise ValueError(f"Column {col} not found in input.")

        x = df[col].values
        dx = x - np.mean(x)

        ac, ac_norm = autocorr(x)
        tau_ps = correlation_time(time - time[0], ac_norm)

        variance = float(np.var(dx))
        std = float(np.std(dx))

        # Very simple motional-narrowing dephasing proxy:
        # gamma_phi ~ variance * tau_c / hbar^2
        gamma_phi_ps = variance * tau_ps / (HBAR_eV_ps ** 2)

        freq, spec = spectral_density_proxy(time - time[0], ac)

        ac_df = pd.DataFrame({
            "lag_ps": time - time[0],
            "autocorrelation": ac,
            "autocorrelation_normalized": ac_norm,
        })
        ac_df.to_csv(out / f"{col}_autocorrelation.csv", index=False)

        spec_df = pd.DataFrame({
            "frequency_ps^-1": freq,
            "spectral_density_proxy": spec,
        })
        spec_df.to_csv(out / f"{col}_spectral_density_proxy.csv", index=False)

        plt.figure(figsize=(6,4))
        plt.plot(time - time[0], ac_norm)
        plt.xlabel("lag time (ps)")
        plt.ylabel("normalized autocorrelation")
        plt.tight_layout()
        plt.savefig(out / f"{col}_autocorrelation.png", dpi=300)
        plt.close()

        plt.figure(figsize=(6,4))
        plt.plot(freq, spec)
        plt.xlabel("frequency (ps^-1)")
        plt.ylabel("spectral density proxy")
        plt.tight_layout()
        plt.savefig(out / f"{col}_spectral_density_proxy.png", dpi=300)
        plt.close()

        rows.append({
            "quantity": col,
            "mean_eV": float(np.mean(x)),
            "std_eV": std,
            "variance_eV2": variance,
            "correlation_time_ps": tau_ps,
            "gamma_phi_proxy_ps^-1": gamma_phi_ps,
            "dt_ps": dt,
            "n_samples": len(time),
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "bath_proxy_summary.csv", index=False)

    md = """# MD Fluctuation to Bath-Parameter Proxy

## Purpose

Convert MD-derived site-energy or coupling time series into first-pass bath parameters for open-system exciton dynamics.

## Required input format

CSV with:

- `time_ps`
- one or more fluctuation columns, for example `PYR2`, `PYR3`, `PYR4`, `PYR5`, `J23`, `J34`, `J45`

## Outputs

For each selected column:

- fluctuation autocorrelation
- normalized autocorrelation
- correlation time
- spectral-density proxy
- dephasing-rate proxy

## Caveat

The dephasing estimate is a motional-narrowing proxy, not a final quantum bath model. It should be used as a diagnostic bridge between MD and Lindblad/Haken-Strobl parameters.
"""
    (out / "MD_TO_BATH_PROXY_README.md").write_text(md)

    print(summary.to_string(index=False))
    print("Wrote:", out)


if __name__ == "__main__":
    main()
