# Day014 Hydrated PYR4-PYR5 State Classification

## System

- PYR4-PYR5 hydrated dimer
- Local explicit water shell: 0.50 nm
- ORCA 6.1.1, PBE0/def2-SVP, TDDFT/TDA, 20 roots
- Oscillator strengths reported here use the transition electric dipole length gauge, fosc(D2).

## S1-S5 classification

| State | Energy (eV) | Wavelength (nm) | fosc(D2) | Dominant transition | Weight | Character |
|---|---:|---:|---:|---|---:|---|
| S1 | 2.959 | 419.0 | 0.000000000 | 315a->316a | 0.999988 | dark inter-chromophore CT-like |
| S2 | 3.746 | 331.0 | 0.371510312 | 314a->316a | 0.794275 | bright Frenkel-like candidate |
| S3 | 3.757 | 330.0 | 0.008126814 | 314a->318a | 0.528051 | weak/mixed transition |
| S4 | 3.779 | 328.1 | 0.001646680 | 315a->321a | 0.548333 | dark or weak transition |
| S5 | 3.812 | 325.2 | 0.424868742 | 315a->317a | 0.829674 | bright Frenkel-like candidate |

## Interpretation

The lowest state, S1 at 2.959 eV, is dark and dominated by MO315a -> MO316a. Fragment localization shows MO315a localized on PYR4 and MO316a localized on PYR5, so S1 is best interpreted as an inter-chromophore charge-transfer-like state rather than a water-localized artifact.

The bright monomer-like states are S2 and S5, with fosc(D2) ≈ 0.372 and 0.425. These are the most appropriate first-pass Frenkel-like pair for estimating the hydrated PYR4-PYR5 coupling.

Using S2 and S5:

- S2 = 3.746 eV
- S5 = 3.812 eV
- Splitting = 0.066 eV
- Preliminary hydrated Frenkel-like J_eff ≈ 33 meV

This value remains preliminary because the states may contain residual mixing. A more rigorous refinement should use transition-density or fragment-based coupling analysis.
