#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


def last_row(path):
    df = pd.read_csv(path)
    return df.iloc[-1].to_dict()


def safe(row, key):
    return row.get(key, None)


def build_record(system, condition, conf_path, dip_path):
    conf = last_row(conf_path)
    dip = last_row(dip_path)

    return {
        "system": system,
        "condition": condition,
        "timestep": safe(conf, "timestep"),
        "n_water": safe(conf, "n_water"),
        "final_water_r_mean_A": safe(conf, "water_r_mean_A"),
        "final_water_r_max_A": safe(conf, "water_r_max_A"),
        "final_water_abs_z_mean_A": safe(conf, "water_abs_z_mean_A"),
        "final_fraction_inside_lumen_segment": safe(conf, "fraction_inside_lumen_segment"),
        "final_fraction_radial_outside_lumen": safe(conf, "fraction_radial_outside_lumen"),
        "final_fraction_axial_outside_segment": safe(conf, "fraction_axial_outside_segment"),
        "final_fraction_outside_outer_radius": safe(conf, "fraction_outside_outer_radius"),
        "final_mean_cos_theta_z": safe(dip, "cos_theta_z_mean"),
        "final_abs_cos_theta_z_mean": safe(dip, "abs_cos_theta_z_mean"),
        "final_S_z": safe(dip, "S_z_mean"),
        "final_mean_cos_theta_radial": safe(dip, "cos_theta_radial_mean"),
        "final_abs_cos_theta_radial_mean": safe(dip, "abs_cos_theta_radial_mean"),
    }


def main():
    records = [
        build_record(
            "BN-like polar scaffold-water",
            "no_field_20k",
            "results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv",
            "results/dipoles/nvt300_hold_extended_contained_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "BN-neutralized scaffold-water",
            "no_field_20k",
            "results/confinement/bn_neutralized_scaffold_water_30000w_nvt300_contained_20k_confinement_summary.csv",
            "results/dipoles/nvt300_contained_bn_neutralized_20k_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "carbon-like neutral scaffold-water",
            "no_field_20k",
            "results/confinement/carbon_like_scaffold_water_30000w_nvt300_contained_20k_confinement_summary.csv",
            "results/dipoles/nvt300_contained_carbon_like_20k_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "water-only contained",
            "no_field_20k",
            "results/confinement/water_only_30000w_nvt300_contained_confinement_summary.csv",
            "results/dipoles/nvt300_water_only_contained_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "BN-like polar scaffold-water",
            "fieldZ_10k",
            "results/confinement/nvt300_fieldZ_scaffold_water_confinement_summary.csv",
            "results/dipoles/nvt300_fieldZ_scaffold_water_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "BN-neutralized scaffold-water",
            "fieldZ_10k",
            "results/confinement/bn_neutralized_scaffold_water_30000w_fieldZ_confinement_summary.csv",
            "results/dipoles/nvt300_fieldZ_bn_neutralized_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "carbon-like neutral scaffold-water",
            "fieldZ_10k",
            "results/confinement/carbon_like_scaffold_water_30000w_fieldZ_confinement_summary.csv",
            "results/dipoles/nvt300_fieldZ_carbon_like_water_dipole_orientation_summary.csv",
        ),
        build_record(
            "water-only contained",
            "fieldZ_10k",
            "results/confinement/nvt300_fieldZ_water_only_confinement_summary.csv",
            "results/dipoles/nvt300_fieldZ_water_only_water_dipole_orientation_summary.csv",
        ),
    ]

    out = Path("results/phase0/day005_charge_polarity_screening_summary.csv")
    out.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(records)
    df.to_csv(out, index=False)

    print(f"Wrote {out}")
    print(df[[
        "system",
        "condition",
        "final_water_r_max_A",
        "final_mean_cos_theta_z",
        "final_S_z",
        "final_fraction_inside_lumen_segment",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
