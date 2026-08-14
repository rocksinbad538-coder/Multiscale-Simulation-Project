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
