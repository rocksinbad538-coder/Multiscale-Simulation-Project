# Day016 site-energy trajectory

## Source

- Input: `runs/phase1A/day016_md_bath_extraction/orca_embedding_analysis/embedding_pilot_summary.csv`
- Quantity: embedded TDDFT/TDA S1 excitation energy
- Units: eV
- Chromophore order: PYR2, PYR3, PYR4, PYR5

## Site energies

 frame  PYR2  PYR3  PYR4  PYR5
     0 3.999 4.000 4.013 3.779
    10 3.998 3.996 4.020 3.791
    20 4.010 4.000 4.016 3.768

## Per-frame centered diagonal terms

 frame    PYR2    PYR3    PYR4     PYR5
     0 0.05125 0.05225 0.06525 -0.16875
    10 0.04675 0.04475 0.06875 -0.16025
    20 0.06150 0.05150 0.06750 -0.18050

## Statistics

         count      mean       std    min     25%    50%     75%    max
cluster                                                                
PYR2       3.0  4.002333  0.006658  3.998  3.9985  3.999  4.0045  4.010
PYR3       3.0  3.998667  0.002309  3.996  3.9980  4.000  4.0000  4.000
PYR4       3.0  4.016333  0.003512  4.013  4.0145  4.016  4.0180  4.020
PYR5       3.0  3.779333  0.011504  3.768  3.7735  3.779  3.7850  3.791

## Interpretation

This table is the first MD-derived embedded-TDDFT site-energy trajectory. It provides the diagonal part of the time-dependent excitonic Hamiltonian. At this stage, off-diagonal couplings are not included.
