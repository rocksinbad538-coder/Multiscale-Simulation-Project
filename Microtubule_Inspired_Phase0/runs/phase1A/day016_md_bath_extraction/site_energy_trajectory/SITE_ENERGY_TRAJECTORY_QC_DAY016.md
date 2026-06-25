# Day016 site-energy trajectory QC

- Frames: [0, 10, 20]
- Chromophores: ['PYR2', 'PYR3', 'PYR4', 'PYR5']
- All energies are embedded TDDFT/TDA S1 values in eV.
- Centered Hamiltonian diagonals use per-frame mean subtraction.

## Site-energy statistics

      count      mean       std    min     25%    50%     75%    max
PYR2    3.0  4.002333  0.006658  3.998  3.9985  3.999  4.0045  4.010
PYR3    3.0  3.998667  0.002309  3.996  3.9980  4.000  4.0000  4.000
PYR4    3.0  4.016333  0.003512  4.013  4.0145  4.016  4.0180  4.020
PYR5    3.0  3.779333  0.011504  3.768  3.7735  3.779  3.7850  3.791

## Centered row sums

-4.163336e-16
-2.220446e-16
-9.992007e-16

## Frame-to-frame S1 changes [meV]

 PYR2  PYR3  PYR4  PYR5
 -1.0  -4.0   7.0  12.0
 12.0   4.0  -4.0 -23.0

## PYR5 detuning relative to PYR2-PYR4 mean [meV]

-225.000000
-213.666667
-240.666667

## QC interpretation

The centered diagonal Hamiltonians are internally consistent: each frame sums to approximately zero after per-frame mean subtraction. PYR2-PYR4 remain tightly clustered near 4.0 eV, while PYR5 is persistently lower. This produces a stable negative PYR5 detuning that must be retained in the next excitonic Hamiltonian model.
