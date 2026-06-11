#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


def last_row(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return df.iloc[-1].copy()


def main():
    outdir = Path("results/phase0")
    outdir.mkdir(parents=True, exist_ok=True)

    cases = [
        {
            "case": "hybrid_field_free_5k",
            "confinement": "results/confinement/bn_like_chromophore_12dipoles_carved_nvt300_contained_5k_confinement_summary.csv",
            "dipole": "results/dipoles/nvt300_contained_bn_like_chromophore_12dipoles_carved_5k_water_dipole_orientation_summary.csv",
            "thermo": "results/thermo/bn_like_chromophore_12dipoles_carved_nvt300_contained_5k_thermo_records.csv",
        },
        {
            "case": "hybrid_field_free_20k",
            "confinement": "results/confinement/bn_like_chromophore_12dipoles_carved_nvt300_contained_extend_5k_to_20k_confinement_summary.csv",
            "dipole": "results/dipoles/nvt300_contained_bn_like_chromophore_12dipoles_carved_extend_5k_to_20k_water_dipole_orientation_summary.csv",
            "thermo": "results/thermo/bn_like_chromophore_12dipoles_carved_nvt300_contained_extend_5k_to_20k_thermo_records.csv",
        },
        {
            "case": "hybrid_fieldZ_10k",
            "confinement": "results/confinement/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_10k_confinement_summary.csv",
            "dipole": "results/dipoles/nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_10k_water_dipole_orientation_summary.csv",
            "thermo": "results/thermo/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_10k_thermo_records.csv",
        },
        {
            "case": "hybrid_fieldZ_20k",
            "confinement": "results/confinement/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_extend_10k_to_20k_confinement_summary.csv",
            "dipole": "results/dipoles/nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_extend_10k_to_20k_water_dipole_orientation_summary.csv",
            "thermo": "results/thermo/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_extend_10k_to_20k_thermo_records.csv",
        },
        {
            "case": "hybrid_fieldZ_30k",
            "confinement": "results/confinement/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_extend_20k_to_30k_confinement_summary.csv",
            "dipole": "results/dipoles/nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_extend_20k_to_30k_water_dipole_orientation_summary.csv",
            "thermo": "results/thermo/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_extend_20k_to_30k_thermo_records.csv",
        },
    ]

    rows = []

    for c in cases:
        conf = last_row(c["confinement"])
        dip = last_row(c["dipole"])
        thermo = last_row(c["thermo"])

        rows.append({
            "case": c["case"],
            "step": int(conf["timestep"]),
            "n_water": int(conf["n_water_oxygen"]),
            "water_r_mean_A": float(conf["water_r_mean_A"]),
            "water_r_max_A": float(conf["water_r_max_A"]),
            "water_abs_z_mean_A": float(conf["water_abs_z_mean_A"]),
            "fraction_inside_lumen_segment": float(conf["fraction_inside_lumen_segment"]),
            "fraction_radial_outside_lumen": float(conf["fraction_radial_outside_lumen"]),
            "fraction_axial_outside_segment": float(conf["fraction_axial_outside_segment"]),
            "fraction_outside_outer_radius": float(conf["fraction_radial_outside_outer_radius"]),
            "cos_theta_z_mean": float(dip["cos_theta_z_mean"]),
            "abs_cos_theta_z_mean": float(dip["abs_cos_theta_z_mean"]),
            "S_z_mean": float(dip["S_z_mean"]),
            "c_twater_K": float(thermo["c_twater"]),
            "PotEng_kcal_mol": float(thermo["PotEng"]),
            "TotEng_kcal_mol": float(thermo["TotEng"]),
            "Press": float(thermo["Press"]),
            "water_msd_A2": float(thermo["c_msdwater[4]"]),
        })

    df = pd.DataFrame(rows)
    out = outdir / "day007_lead_hybrid_candidate_summary.csv"
    df.to_csv(out, index=False)

    print(df.to_string(index=False))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
