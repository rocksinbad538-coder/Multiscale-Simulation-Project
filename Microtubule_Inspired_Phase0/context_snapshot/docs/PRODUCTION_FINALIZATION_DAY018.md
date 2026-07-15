# Day018 embedded TDDFT production finalization

## Production status

- ORCA calculations: 84/84
- Complete MD frames: 21/21
- Chromophores per frame: 4
- Normal ORCA terminations: 84/84
- SCF converged: 84/84
- TDDFT/TDA completed: 84/84
- Point-charge count mismatches: 0
- Non-neutral embedding files: 0

## Physical scope

The source MD trajectory has frozen h-BN and pyrene coordinates. The time-dependent variation therefore represents solvent-induced electrostatic fluctuations for fixed solute geometries.

## State tracking

The local bright state is selected from S1 and S2 using the largest 52a->53a HOMO-LUMO configuration weight rather than root number alone.

cluster  tracked_root  n_jobs
   PYR2             2      21
   PYR3             2      21
   PYR4             2      21
   PYR5             1      21

The previously observed root-ordering pattern remains unchanged across all 21 frames.

## Tracked-state statistics

cluster  n_frames  mean_energy_eV  std_energy_eV  min_energy_eV  max_energy_eV  mean_fosc  mean_character_weight  minimum_character_weight
   PYR2        21        4.088571       0.011084          4.068          4.104   0.493290               0.830708                  0.738175
   PYR3        21        4.078048       0.011902          4.052          4.096   0.501778               0.836158                  0.793849
   PYR4        21        4.090619       0.007586          4.069          4.104   0.510859               0.841515                  0.823216
   PYR5        21        3.776952       0.016430          3.750          3.797   0.648089               0.863937                  0.827679

## PYR5 site offset

- Mean offset relative to PYR2-PYR4: 308.794 meV
- Standard deviation: 18.483 meV
- Minimum: 272.667 meV
- Maximum: 340.667 meV

## Generated Hamiltonians

- Absolute CSV matrices: 21
- Absolute NPY matrices: 21
- Centered CSV matrices: 21
- Centered NPY matrices: 21

## Interpretation constraint

The tracked diagonal energies are electronically consistent, but the PYR5 offset cannot yet be assigned exclusively to solvent electrostatics because its frozen geometry differs measurably from PYR2-PYR4. Vacuum reference calculations and representative NTO analysis remain necessary.
