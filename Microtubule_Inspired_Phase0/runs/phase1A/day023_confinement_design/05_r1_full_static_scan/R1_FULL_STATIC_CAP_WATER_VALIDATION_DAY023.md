# R1 Full-System Static Cap–Water Validation

## Scope

Four complete R1 topologies were constructed from the validated
68,314-atom geometry. Each topology contains:

- the accepted HBN and PYR definitions;
- 16551 TIP4P/2005 waters;
- one 163-bead lower cap;
- one 163-bead upper cap;
- zero cap charge;
- zero CAP–CAP, CAP–HBN, CAP–PYR, CAP–water-H, and CAP–water-M
  interaction;
- one explicit pure-r12 CAP–OW interaction.

No energy minimization or molecular dynamics was performed.

## Static validation

The full system was evaluated using `gmx grompp` and
`gmx mdrun -rerun`. The CAP–water energy was independently calculated
as the explicit sum of all CAP–OW C12/r^12 contributions inside the
1.0 nm static cutoff.

- `target_5kBT`: CAP–water energy=14.302773 kJ/mol; min distance=0.220189 nm; maximum water force=52.263 kJ mol^-1 nm^-1; hole barrier=551.5 kBT; **PASS**
- `target_10kBT`: CAP–water energy=28.605547 kJ/mol; min distance=0.220189 nm; maximum water force=104.527 kJ mol^-1 nm^-1; hole barrier=1103.0 kBT; **PASS**
- `target_20kBT`: CAP–water energy=57.211094 kJ/mol; min distance=0.220189 nm; maximum water force=209.054 kJ mol^-1 nm^-1; hole barrier=2206.0 kBT; **PASS**
- `target_40kBT`: CAP–water energy=114.422188 kJ/mol; min distance=0.220189 nm; maximum water force=418.108 kJ mol^-1 nm^-1; hole barrier=4411.9 kBT; **PASS**

## Selected model

- Candidate:
  **target_5kBT**
- Selection rule:
  **weakest validated repulsion satisfying a minimum 100 kBT barrier
  at the largest coverage hole**
- Pair sigma:
  **-0.170000 nm**
- Pair epsilon:
  **3.1179234818 kJ/mol**
- C12:
  **7.266286217927e-09
  kJ mol^-1 nm^12**
- Initial minimum CAP–OW distance:
  **0.220189 nm**
- Initial CAP–water energy:
  **14.302773 kJ/mol**
- Maximum initial water force from the cap:
  **52.263
  kJ mol^-1 nm^-1**
- Energy at the 0.22 nm pruning boundary:
  **0.226613 kBT**
- Barrier at the largest 0.114878 nm coverage hole:
  **551.489 kBT**
- GROMACS/analytic relative energy error:
  **1.758004e-06**

## Interaction checks

- CAP–HBN/PYR LJ energy:
  **0.000000000000 kJ/mol**
- CAP–HBN/PYR Coulomb energy:
  **0.000000000000 kJ/mol**
- CAP–water Coulomb energy:
  **0.000000000000 kJ/mol**

## Decision

- Full R1 topology validated: **YES**
- Selected cap model:
  **target_5kBT**
- Cap atoms must remain frozen: **YES**
- HBN and PYR must remain frozen during initial screening: **YES**
- Energy minimization authorized: **YES**
- Molecular dynamics authorized: **NO**
- Required next step:
  `PREPARE_AND_RUN_R1_FROZEN_SOLUTE_ENERGY_MINIMIZATION`

Selected files:

- `runs/phase1A/day023_confinement_design/05_r1_full_static_scan/selected/r1_selected_cap_model.top`
- `runs/phase1A/day023_confinement_design/05_r1_full_static_scan/selected/r1_selected_static_validation.tpr`
- `runs/phase1A/day023_confinement_design/05_r1_full_static_scan/selected/r1_selected_groups.ndx`
- `runs/phase1A/day023_confinement_design/05_r1_full_static_scan/r1_selected_cap_model.json`
