#!/usr/bin/env python3
"""
Generate a LAMMPS read-test input for BN-like scaffold + confined water systems.

This avoids manually rewriting LAMMPS input files for each water-count scaling test.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def write_input(n_water: int, output_path: Path) -> None:
    data_file = (
        f"../../../systems/inorganic/bn_like_scaffold_water/lammps/"
        f"bn_like_scaffold_water_24OD_14ID_20L_{n_water}w.data"
    )

    checked_file = (
        f"../../../systems/inorganic/bn_like_scaffold_water/outputs/"
        f"bn_like_scaffold_water_{n_water}w_checked.data"
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    text = f"""# Phase 0 — BN-like scaffold + confined water LAMMPS read test
# Water molecules: {n_water}
# Purpose: verify that the generated scaffold-water LAMMPS data file can be read correctly.

units           real
atom_style      full
boundary        f f f

read_data       {data_file}

group           scaffold type 1 2
group           water type 3 4

mass            1 10.811
mass            2 14.007
mass            3 15.9994
mass            4 1.008

pair_style      lj/cut/coul/cut 10.0

# Placeholder nonbonded parameters for read/initialization only.
# Not final production force field.
pair_coeff      1 1 0.0500 3.40
pair_coeff      2 2 0.0500 3.40
pair_coeff      1 2 0.0500 3.40
pair_coeff      3 3 0.1521 3.1507
pair_coeff      4 4 0.0000 0.0000
pair_coeff      1 3 0.0872 3.275
pair_coeff      2 3 0.0872 3.275
pair_coeff      1 4 0.0000 0.0000
pair_coeff      2 4 0.0000 0.0000
pair_coeff      3 4 0.0000 0.0000

bond_style      harmonic
angle_style     harmonic

# Placeholder strong water geometry terms.
# Production may use SHAKE instead.
bond_coeff      1 450.0 0.9572
angle_coeff     1 55.0 104.52

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

thermo          1
thermo_style    custom step atoms temp pe ke etotal press

run             0

write_data      {checked_file}
"""

    output_path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-water", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    write_input(args.n_water, args.output)
    print(f"Wrote LAMMPS read-test input: {args.output}")


if __name__ == "__main__":
    main()
