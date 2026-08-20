from __future__ import annotations

import numpy as np
import pandas as pd


def block_statistics(df, columns, nblocks=5):

    n = len(df)

    edges = np.linspace(
        0,
        n,
        nblocks + 1,
        dtype=int,
    )

    rows = []

    for i in range(nblocks):

        part = df.iloc[
            edges[i]:edges[i + 1]
        ]

        row = {
            "Block": i + 1,
            "Frames": len(part),
        }

        for c in columns:

            row[f"{c}_mean"] = part[c].mean()
            row[f"{c}_std"] = part[c].std()

        rows.append(row)

    return pd.DataFrame(rows)


def running_average(x):

    x = np.asarray(x)

    return np.cumsum(x) / np.arange(
        1,
        len(x) + 1,
    )


def running_std(x):

    x = np.asarray(x)

    out = np.zeros_like(
        x,
        dtype=float,
    )

    for i in range(len(x)):

        out[i] = np.std(
            x[: i + 1]
        )

    return out


def autocorrelation(x, max_lag=None):
    """
    Normalized autocorrelation function.

    Returns rho(k) for k = 0 ... max_lag.
    """
    x = np.asarray(x, dtype=float)

    if x.ndim != 1:
        raise ValueError("autocorrelation expects a 1D array")

    n = len(x)

    if n < 2:
        raise ValueError("At least two samples are required")

    x = x - np.mean(x)

    variance = np.dot(x, x) / n

    if variance <= 0.0:
        return np.ones(1, dtype=float)

    if max_lag is None:
        max_lag = min(n - 1, n // 2)

    max_lag = int(min(max_lag, n - 1))

    rho = np.empty(max_lag + 1, dtype=float)
    rho[0] = 1.0

    for lag in range(1, max_lag + 1):
        covariance = np.dot(
            x[:-lag],
            x[lag:]
        ) / (n - lag)

        rho[lag] = covariance / variance

    return rho


def integrated_autocorrelation_time(x, max_lag=None):
    """
    Integrated autocorrelation time estimated using Geyer's
    initial-positive-sequence (IPS) truncation.

    For normalized autocorrelation rho(k), adjacent positive-lag
    terms are grouped as

        Gamma_m = rho(2m-1) + rho(2m)

    and accumulation stops at the first non-positive Gamma_m.

    tau_int = 0.5 + sum_m Gamma_m

    The result is expressed in units of samples.
    """
    rho = autocorrelation(
        x,
        max_lag=max_lag
    )

    tau = 0.5

    # Positive lags are grouped into adjacent pairs:
    # (rho_1 + rho_2), (rho_3 + rho_4), ...
    k = 1

    while k < len(rho):

        if k + 1 < len(rho):
            pair_sum = rho[k] + rho[k + 1]
        else:
            pair_sum = rho[k]

        if not np.isfinite(pair_sum):
            break

        if pair_sum <= 0.0:
            break

        tau += pair_sum

        k += 2

    return max(float(tau), 0.5)


def effective_sample_size(x, max_lag=None):
    """
    Effective number of statistically independent samples.

    N_eff = N / (2 tau_int)
    """
    x = np.asarray(x, dtype=float)

    tau = integrated_autocorrelation_time(
        x,
        max_lag=max_lag
    )

    neff = len(x) / (2.0 * tau)

    return float(min(len(x), max(1.0, neff)))


def statistical_inefficiency(x, max_lag=None):
    """
    Statistical inefficiency:
        g = 2 tau_int
    """
    tau = integrated_autocorrelation_time(
        x,
        max_lag=max_lag
    )

    return float(2.0 * tau)
