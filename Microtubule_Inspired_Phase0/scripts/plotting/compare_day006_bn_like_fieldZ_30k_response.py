#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def load_confinement(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing confinement file: {path}")
    df = pd.read_csv(path)
    print(f"Loaded confinement: {path}")
    print("Columns:", list(df.columns))
    return df


def load_dipole(path):
    """
    Load dipole-orientation summary/records robustly.

    Expected structure in current workflow:
    label,timestep,n_water,cos_z_mean,cos_z_std,abs_cos_z_mean,S_z,...
    but column names may vary depending on script version.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing dipole file: {path}")

    df = pd.read_csv(path)
    print(f"Loaded dipoles: {path}")
    print("Columns:", list(df.columns))

    # Ensure timestep column.
    if "timestep" not in df.columns:
        # fallback: second column is timestep
        df["timestep"] = df.iloc[:, 1]

    # Robust cos_z_mean.
    if "cos_z_mean" not in df.columns:
        possible = [
            "mean_cos_theta_z",
            "cos_theta_z_mean",
            "mean_cos_z",
            "cosz_mean",
        ]
        found = None
        for col in possible:
            if col in df.columns:
                found = col
                break
        if found is not None:
            df["cos_z_mean"] = df[found]
        else:
            # fallback: fourth column, based on current CSV row structure
            df["cos_z_mean"] = df.iloc[:, 3]

    # Robust S_z.
    if "S_z" not in df.columns:
        possible = [
            "S_z_mean",
            "Sz",
            "S",
            "order_parameter_z",
        ]
        found = None
        for col in possible:
            if col in df.columns:
                found = col
                break
        if found is not None:
            df["S_z"] = df[found]
        else:
            # fallback: seventh column, based on current CSV row structure
            df["S_z"] = df.iloc[:, 6]

    return df


def main():
    outdir = Path("figures/day006")
    outdir.mkdir(parents=True, exist_ok=True)

    # Exact files present in your repository.
    conf_10k_path = "results/confinement/nvt300_fieldZ_scaffold_water_confinement_summary.csv"
    dip_10k_path = "results/dipoles/nvt300_fieldZ_scaffold_water_water_dipole_orientation_summary.csv"

    conf_ext_path = "results/confinement/bn_like_scaffold_water_30000w_nvt300_fieldZ_contained_extend_10k_to_30k_confinement_summary.csv"
    dip_ext_path = "results/dipoles/nvt300_fieldZ_contained_bn_like_extend_10k_to_30k_water_dipole_orientation_summary.csv"

    conf_10k = load_confinement(conf_10k_path)
    dip_10k = load_dipole(dip_10k_path)
    conf_ext = load_confinement(conf_ext_path)
    dip_ext = load_dipole(dip_ext_path)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Panel A: maximum radial water oxygen position.
    axes[0, 0].plot(
        conf_10k["timestep"],
        conf_10k["water_r_max_A"],
        label="BN-like fieldZ 0-10k",
    )
    axes[0, 0].plot(
        conf_ext["timestep"],
        conf_ext["water_r_max_A"],
        label="BN-like fieldZ 10k-30k",
    )
    axes[0, 0].axhline(70.0, linestyle="--", linewidth=1, label="nominal lumen radius")
    axes[0, 0].set_title("Maximum radial water oxygen position")
    axes[0, 0].set_xlabel("fieldZ step")
    axes[0, 0].set_ylabel("Max r, Angstrom")
    axes[0, 0].legend(fontsize=8)

    # Panel B: fraction inside lumen segment.
    axes[0, 1].plot(
        conf_10k["timestep"],
        conf_10k["fraction_inside_lumen_segment"],
        label="BN-like fieldZ 0-10k",
    )
    axes[0, 1].plot(
        conf_ext["timestep"],
        conf_ext["fraction_inside_lumen_segment"],
        label="BN-like fieldZ 10k-30k",
    )
    axes[0, 1].set_title("Fraction inside nominal lumen segment")
    axes[0, 1].set_xlabel("fieldZ step")
    axes[0, 1].set_ylabel("Fraction inside")
    axes[0, 1].legend(fontsize=8)

    # Panel C: mean cos(theta_z).
    axes[1, 0].plot(
        dip_10k["timestep"],
        dip_10k["cos_z_mean"],
        label="BN-like fieldZ 0-10k",
    )
    axes[1, 0].plot(
        dip_ext["timestep"],
        dip_ext["cos_z_mean"],
        label="BN-like fieldZ 10k-30k",
    )
    axes[1, 0].set_title("Axial water-dipole projection")
    axes[1, 0].set_xlabel("fieldZ step")
    axes[1, 0].set_ylabel("mean cos(theta_z)")
    axes[1, 0].legend(fontsize=8)

    # Panel D: S_z.
    axes[1, 1].plot(
        dip_10k["timestep"],
        dip_10k["S_z"],
        label="BN-like fieldZ 0-10k",
    )
    axes[1, 1].plot(
        dip_ext["timestep"],
        dip_ext["S_z"],
        label="BN-like fieldZ 10k-30k",
    )
    axes[1, 1].set_title("Axial dipolar order parameter")
    axes[1, 1].set_xlabel("fieldZ step")
    axes[1, 1].set_ylabel("S_z")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Day 006 BN-like polar scaffold-water fieldZ extension to 30k", fontsize=15)
    fig.tight_layout()

    out = outdir / "day006_bn_like_fieldZ_30k_response.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
