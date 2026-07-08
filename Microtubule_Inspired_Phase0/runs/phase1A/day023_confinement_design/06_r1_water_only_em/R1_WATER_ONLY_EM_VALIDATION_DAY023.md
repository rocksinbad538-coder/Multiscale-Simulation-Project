# R1 Water-Only Energy Minimization

## Scope

The validated R1 positive-control geometry was minimized with:

- HBN frozen;
- all pyrene atoms frozen;
- both steric caps frozen;
- only the 16551 TIP4P/2005 water molecules mobile.

No molecular dynamics was performed.

## Inputs

- Coordinates:
  `runs/phase1A/day023_confinement_design/06_r1_water_only_em/r1_water_only_em_input.gro`
- Topology:
  `runs/phase1A/day023_confinement_design/06_r1_water_only_em/r1_water_only_em.top`
- Index:
  `runs/phase1A/day023_confinement_design/06_r1_water_only_em/r1_water_only_em.ndx`
- EM template:
  `runs/phase1A/day021_mobile_restraint_protocol/protocol_inputs/mdp/01_em_k10000.mdp`
- EM tolerance:
  **500.0 kJ mol^-1 nm^-1**
- Maximum steps:
  **50000**

## Energy minimization result

- Grompp return code:
  **0**
- Mdrun return code:
  **0**
- Finished mdrun:
  **YES**
- Explicit convergence message:
  **YES**
- Final maximum force:
  **489.43253
  kJ mol^-1 nm^-1**
- Final force norm:
  **16.7473
  kJ mol^-1 nm^-1**

## Energy

- Initial potential:
  **-899828.687500 kJ/mol**
- Final potential:
  **-899828.687500 kJ/mol**
- Potential change:
  **0.000000 kJ/mol**
- Initial CAP–SOL LJ:
  **14.324902 kJ/mol**
- Final CAP–SOL LJ:
  **14.324902 kJ/mol**

## Frozen-group integrity

- HBN RMS/max displacement:
  **0.000000000000/0.000000000000 nm**
- PYR RMS/max displacement:
  **0.000000000000/0.000000000000 nm**
- Cap RMS/max displacement:
  **0.000000000000/0.000000000000 nm**

## Water relaxation

- Water-O RMS displacement:
  **0.000009 nm**
- Water-O maximum displacement:
  **0.001000 nm**
- Initial/final minimum CAP–OW distance:
  **0.220189/
  0.220189 nm**
- Initial/final lumen occupancy:
  **428/
  428 waters**
- Lumen occupancy change:
  **0 waters**

## Validation gates

- `grompp_return_code_zero`: **PASS**
- `mdrun_return_code_zero`: **PASS**
- `mdrun_finished`: **PASS**
- `no_instability_signatures`: **PASS**
- `final_atom_count_is_68314`: **PASS**
- `box_is_unchanged`: **PASS**
- `HBN_is_frozen`: **PASS**
- `PYR_is_frozen`: **PASS**
- `caps_are_frozen`: **PASS**
- `energy_is_finite`: **PASS**
- `potential_did_not_increase`: **PASS**
- `force_convergence`: **PASS**
- `final_cap_water_distance`: **PASS**
- `lumen_water_retention`: **PASS**
- `water_displacement_is_finite`: **PASS**
- `water_displacement_is_local`: **PASS**
- `cap_water_energy_did_not_increase`: **PASS**

## Decision

- Decision:
  **R1_WATER_ONLY_EM_VALIDATED**
- Failed gates:
  **NONE**
- Short frozen-solute NVT preparation authorized:
  **YES**
- Molecular dynamics execution authorized:
  **NO**
- Required next step:
  `PREPARE_R1_FROZEN_SOLUTE_SHORT_NVT_SCREENING`

The next stage must prepare and statically validate a short
frozen-solute NVT screening before any dynamics are executed.
