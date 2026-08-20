#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs"
    / "phase2"
    / "day044_md_protocol"
    / "MD_PROTOCOL.json"
)

OUT = (
    ROOT
    / "runs"
    / "phase2"
    / "day044_md_protocol"
)

cfg = json.loads(PROTOCOL.read_text())

# ------------------------------------------------------------------
# Automatic conversion of physical simulation times into MD steps
# ------------------------------------------------------------------

HEATING_STEPS = int(
    cfg["heating_ps"] * 1000.0 / cfg["timestep_fs"]
)

EQUILIBRATION_STEPS = int(
    cfg["equilibration_ps"] * 1000.0 / cfg["timestep_fs"]
)

PRODUCTION_STEPS = int(
    cfg["production_ns"] * 1_000_000.0 / cfg["timestep_fs"]
)



def ps_to_steps(ps):
    return int(ps * 1000.0 / cfg["timestep_fs"])

def ns_to_steps(ns):
    return int(ns * 1_000_000.0 / cfg["timestep_fs"])


DATA = (
    "../../phase1B/day041_lammps_export/data.lammps"
)


def header(datafile):

    return f"""units {cfg["units"]}
atom_style {cfg["atom_style"]}

boundary p p p

pair_style {cfg["pair_style"]} {cfg["pair_cutoff_A"]}

pair_modify mix {cfg["pair_mixing_rule"]}

kspace_style {cfg["kspace_style"]} {cfg["kspace_accuracy"]}

bond_style {cfg["bond_style"]}
angle_style {cfg["angle_style"]}
dihedral_style {cfg["dihedral_style"]} nocoeff
improper_style {cfg["improper_style"]}

special_bonds lj/coul 0.0 0.0 0.0

read_data {datafile}

neighbor {cfg["neighbor_skin_A"]} {cfg["neighbor_style"]}
neigh_modify every 1 delay 0 check yes

timestep {cfg["timestep_fs"]}

thermo_style custom step temp pe ke etotal ebond eangle edihed eimp evdwl press

"""

def write(name, datafile, body):

    path = OUT / name

    path.write_text(
        header(datafile) + body
    )

    print(path)


write(

    "in.minimize",
    "../../phase1B/day041_lammps_export/data.lammps",

    """
log minimize.log

reset_timestep 0

thermo 100

min_style cg

minimize 1.0e-8 1.0e-8 5000 10000

write_data minimized.data
"""

)

write(

    "in.heating",
    "minimized.data",

    f"""
log heating.log\n\nreset_timestep 0\n\nvelocity all create {cfg["temperature_initial_K"]} 12345 mom yes rot yes dist gaussian

fix 1 all nvt temp {cfg["temperature_initial_K"]} {cfg["temperature_final_K"]} 100.0

thermo {cfg["thermo_every"]}

dump 1 all custom {cfg["dump_every"]} heating.xyz id type x y z\n\ndump_modify 1 sort id

run {ps_to_steps(cfg["heating_ps"])}

write_data heated.data

unfix 1
"""
)

write(

    "in.nvt",
    "heated.data",

    f"""
log nvt.log\n\nreset_timestep 0\n\nfix 1 all nvt temp 300.0 300.0 100.0

thermo {cfg["thermo_every"]}

dump 1 all custom {cfg["dump_every"]} nvt.xyz id type x y z\n\ndump_modify 1 sort id

run {ps_to_steps(cfg["equilibration_ps"])}

write_data equilibrated.data

unfix 1
"""
)

write(

    "in.production",
    "equilibrated.data",

    f"""
log production.log

reset_timestep 0

fix 1 all nvt temp 300.0 300.0 100.0

thermo {cfg["thermo_every"]}

dump 1 all custom {cfg["dump_every"]} production.xyz id type x y z\n\ndump_modify 1 sort id\n\nrestart {cfg["restart_every"]} production.restart

run {ns_to_steps(cfg["production_ns"])}\n\nunfix 1
"""
)

print()
print("="*90)
print("DAY044 / PHASE2-A13")
print("MD INPUT FILES GENERATED")
print("="*90)
