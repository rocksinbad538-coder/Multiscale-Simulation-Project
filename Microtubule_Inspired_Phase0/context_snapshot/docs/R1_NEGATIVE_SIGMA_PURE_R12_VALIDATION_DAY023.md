# R1 Negative-Sigma Pure-r12 Validation

## Correction to the previous topology audit

The previous audit treated `comb-rule = 2` as incompatible with a
pure C12 interaction. GROMACS provides a documented special case:
when sigma is negative, C6 is set to zero and C12 is calculated from
the absolute sigma value and epsilon.

Therefore, a CAP–OW pure repulsive interaction can be represented as:

- sigma = **-0.170000 nm**
- epsilon > 0
- C6 = **0**
- C12 = **4 epsilon |sigma|^12**

The previous contract is retained for provenance but is superseded by:

`runs/phase1A/day023_confinement_design/04_r1_negative_sigma_validation/r1_cap_nonbonded_model_contract_corrected.json`

## Microvalidation setup

- GROMACS executable: `/usr/local/gromacs/bin/gmx`
- Temperature reference: **300.0 K**
- kBT: **2.49433879 kJ mol^-1**
- Baseline OW sigma:
  **0.3158900000 nm**
- Baseline OW epsilon:
  **0.7749000000 kJ mol^-1**
- Negative pair sigma:
  **-0.170000 nm**
- Distances:
  **0.150, 0.170, 0.200, 0.220, 0.250, 0.300, 0.500 nm**

Each candidate was evaluated using `grompp`, `mdrun -rerun`, and the
CAP–OW short-range Lennard-Jones energy extracted from the resulting
energy file.

## Candidate results

- target_5kBT: epsilon=3.11792348 kJ mol^-1; C12=7.266286217927e-09 kJ mol^-1 nm^12; max relative error=3.063927e-06; **PASS**
- target_10kBT: epsilon=6.23584696 kJ mol^-1; C12=1.453257243585e-08 kJ mol^-1 nm^12; max relative error=3.063927e-06; **PASS**
- target_20kBT: epsilon=12.47169393 kJ mol^-1; C12=2.906514487171e-08 kJ mol^-1 nm^12; max relative error=3.063927e-06; **PASS**
- target_40kBT: epsilon=24.94338785 kJ mol^-1; C12=5.813028974342e-08 kJ mol^-1 nm^12; max relative error=3.063927e-06; **PASS**

## Decision

- Decision: **PURE_R12_NEGATIVE_SIGMA_OVERRIDE_VALIDATED**
- All candidates passed:
  **YES**
- Pure-r12 standard nonbonded override feasible:
  **YES**
- Tabulated interaction required: **NO**
- Full R1 topology built: **NO**
- Energy minimization authorized: **NO**
- MD execution authorized: **NO**
- Required next step:
  `BUILD_FULL_R1_TOPOLOGY_CANDIDATES_AND_RUN_STATIC_CAP_WATER_ENERGY_SCAN`

The next gate must evaluate the candidate interaction strengths in the
complete 68,314-atom R1 system before selecting a production value.
