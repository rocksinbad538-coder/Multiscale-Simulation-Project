# PYR4-PYR5 Hydrated Dimer TDDFT Pilot

## Calculation

- System: PYR4-PYR5 hydrated dimer, local water shell 0.50 nm
- Method: ORCA 6.1.1, PBE0/def2-SVP, TDDFT/TDA
- Roots: 20
- Status: ORCA terminated normally
- Final single-point energy: -4433.282643972761 Eh

## Key excited states

| State | Energy (eV) | Energy (cm^-1) |
|---|---:|---:|
| S1 | 2.959 | 23864.7 |
| S2 | 3.746 | 30212.2 |
| S3 | 3.757 | 30302.9 |
| S4 | 3.779 | 30479.2 |
| S5 | 3.812 | 30746.7 |

## Immediate interpretation

The hydrated PYR4-PYR5 dimer converged successfully and produced a low-energy first excited state at 2.959 eV, substantially below the monomer-like pyrene transitions around 3.75-3.81 eV. This suggests that explicit local hydration and dimer/environment effects introduce a distinct low-energy excitonic or charge-transfer-like state that must be inspected before treating the splitting as a simple two-state excitonic coupling.

## Files

- `parsed/PYR4_PYR5_water0p50_tddft_states.csv`
- `parsed/PYR4_PYR5_water0p50_absorption_block.txt`
- `PYR4_PYR5_water0p50_serial_tight.out`
- `PYR4_PYR5_water0p50_serial_tight.inp`
