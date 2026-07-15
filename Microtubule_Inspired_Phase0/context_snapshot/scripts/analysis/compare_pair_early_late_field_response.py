#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


INPUTS = {
    "original": Path(
        "results/contacts_by_pair/"
        "day007_fieldZ_0_to_20k_contact_by_frame_and_pair.csv"
    ),
    "replica01": Path(
        "results/contacts_by_pair/"
        "day008_fieldZ_replica01_0_to_20k_contact_by_frame_and_pair.csv"
    ),
}

OUTPUT_DIR = Path("results/contacts_by_pair")
FIGURE_DIR = Path("figures/day009")

EARLY_MIN = 0
EARLY_MAX = 5000
LATE_MIN = 15000
LATE_MAX = 20000


def summarize_condition(name: str, path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    required = {
        "pair_index",
        "positive_atom_id",
        "negative_atom_id",
        "timestep",
        "contact_count",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    early = df[
        df["timestep"].between(EARLY_MIN, EARLY_MAX, inclusive="both")
    ]
    late = df[
        df["timestep"].between(LATE_MIN, LATE_MAX, inclusive="both")
    ]

    group_columns = [
        "pair_index",
        "positive_atom_id",
        "negative_atom_id",
    ]

    early_summary = (
        early.groupby(group_columns, as_index=False)
        .agg(
            early_mean=("contact_count", "mean"),
            early_std=("contact_count", "std"),
            early_n_frames=("contact_count", "size"),
        )
    )

    late_summary = (
        late.groupby(group_columns, as_index=False)
        .agg(
            late_mean=("contact_count", "mean"),
            late_std=("contact_count", "std"),
            late_n_frames=("contact_count", "size"),
        )
    )

    result = early_summary.merge(
        late_summary,
        on=group_columns,
        how="outer",
        validate="one_to_one",
    )

    result["condition"] = name
    result["delta_late_minus_early"] = (
        result["late_mean"] - result["early_mean"]
    )

    result["relative_change"] = np.where(
        result["early_mean"] > 0,
        result["delta_late_minus_early"] / result["early_mean"],
        np.nan,
    )

    return result


def classify_sign(value: float, tolerance: float = 0.20) -> str:
    """
    Changes smaller than 0.20 water molecules per frame are treated
    as approximately unchanged rather than assigning a fragile sign.
    """
    if value > tolerance:
        return "increase"
    if value < -tolerance:
        return "decrease"
    return "unchanged"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    original = summarize_condition("original", INPUTS["original"])
    replica = summarize_condition("replica01", INPUTS["replica01"])

    original.to_csv(
        OUTPUT_DIR / "day009_original_pair_early_late_summary.csv",
        index=False,
    )
    replica.to_csv(
        OUTPUT_DIR / "day009_replica01_pair_early_late_summary.csv",
        index=False,
    )

    keys = [
        "pair_index",
        "positive_atom_id",
        "negative_atom_id",
    ]

    comparison = original[keys + [
        "early_mean",
        "late_mean",
        "delta_late_minus_early",
    ]].rename(
        columns={
            "early_mean": "original_early_mean",
            "late_mean": "original_late_mean",
            "delta_late_minus_early": "original_delta",
        }
    )

    comparison = comparison.merge(
        replica[keys + [
            "early_mean",
            "late_mean",
            "delta_late_minus_early",
        ]].rename(
            columns={
                "early_mean": "replica_early_mean",
                "late_mean": "replica_late_mean",
                "delta_late_minus_early": "replica_delta",
            }
        ),
        on=keys,
        how="outer",
        validate="one_to_one",
    )

    comparison["original_response"] = comparison[
        "original_delta"
    ].map(classify_sign)

    comparison["replica_response"] = comparison[
        "replica_delta"
    ].map(classify_sign)

    comparison["same_response_class"] = (
        comparison["original_response"]
        == comparison["replica_response"]
    )

    comparison["both_active"] = (
        comparison[
            [
                "original_early_mean",
                "original_late_mean",
                "replica_early_mean",
                "replica_late_mean",
            ]
        ].max(axis=1)
        >= 0.5
    )

    output_csv = (
        OUTPUT_DIR
        / "day009_pair_early_late_field_response_comparison.csv"
    )
    comparison.to_csv(output_csv, index=False)

    print("\nEarly window: steps 0–5000")
    print("Late window: steps 15000–20000")
    print("Response tolerance: ±0.20 water molecules per frame\n")

    print(comparison.to_string(index=False))

    active = comparison[comparison["both_active"]]

    print(
        "\nSame response class, all pairs: "
        f"{comparison['same_response_class'].sum()}/{len(comparison)}"
    )
    print(
        "Same response class, active pairs: "
        f"{active['same_response_class'].sum()}/{len(active)}"
    )
    print(f"\nWrote {output_csv}")

    x = comparison["pair_index"].to_numpy()
    width = 0.36

    fig, axis = plt.subplots(figsize=(11, 6))

    axis.bar(
        x - width / 2,
        comparison["original_delta"],
        width=width,
        label="Day 007 original",
    )
    axis.bar(
        x + width / 2,
        comparison["replica_delta"],
        width=width,
        label="Day 008 replica01",
    )

    axis.axhline(0.0, linewidth=1)
    axis.axhline(0.20, linestyle="--", linewidth=1)
    axis.axhline(-0.20, linestyle="--", linewidth=1)

    axis.set_xticks(x)
    axis.set_xlabel("Chromophore dipole-pair index")
    axis.set_ylabel(
        "Late minus early occupancy\n"
        "(water oxygens within 6 Å)"
    )
    axis.set_title(
        "Pair-resolved temporal field response\n"
        "late fieldZ occupancy minus early fieldZ occupancy"
    )
    axis.legend()
    fig.tight_layout()

    output_png = (
        FIGURE_DIR
        / "day009_pair_early_late_field_response.png"
    )
    fig.savefig(output_png, dpi=200)
    plt.close(fig)

    print(f"Wrote {output_png}")


if __name__ == "__main__":
    main()
