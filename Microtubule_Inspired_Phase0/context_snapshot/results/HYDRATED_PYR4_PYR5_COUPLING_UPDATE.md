# Hydrated PYR4-PYR5 Coupling Update

## System

- Pair: PYR4-PYR5
- Environment: explicit local water shell, 0.50 nm
- Method: ORCA 6.1.1, PBE0/def2-SVP, TDDFT/TDA
- Roots: 20

## State assignment

The hydrated dimer contains a low-energy S1 state at 2.959 eV. This state is dark and dominated by:

- MO315a -> MO316a
- excitation weight = 0.999988
- oscillator strength = 0.000000000

Fragment-localization analysis shows:

| MO | PYR4 fraction | PYR5 fraction | Water fraction |
|---|---:|---:|---:|
| 315a | 0.995615 | 0.000000 | 0.004385 |
| 316a | 0.000000 | 0.986427 | 0.013573 |

Therefore, S1 is best assigned as a dark inter-chromophore PYR4 -> PYR5 charge-transfer-like state, not as a water-localized artifact.

## Frenkel-like bright-state pair

The first bright chromophore-centered states are:

| State | Energy (eV) | fosc(D2) | Dominant transition | Character |
|---|---:|---:|---|---|
| S2 | 3.746 | 0.371510312 | 314a -> 316a | bright Frenkel-like candidate |
| S5 | 3.812 | 0.424868742 | 315a -> 317a | bright Frenkel-like candidate |

Using S2 and S5 as the Frenkel-like pair:

- Splitting = 3.812 - 3.746 = 0.066 eV
- J_eff = splitting / 2 = 0.033 eV
- J_eff = 33 meV

## Interpretation

The hydrated PYR4-PYR5 pair should not be modeled using the lowest two excited states directly, because S1 is a dark charge-transfer-like state.

A conservative Frenkel excitonic update should use the S2/S5 bright-state pair, giving a preliminary hydrated PYR4-PYR5 coupling of approximately 33 meV.

This is larger than the previous vacuum splitting-derived PYR4-PYR5 coupling estimate of approximately 11 meV, indicating that explicit local hydration and state mixing can substantially modify the effective excitonic coupling landscape.

## Recommended Hamiltonian handling

For the current Phase 1A Hamiltonian, keep both values documented:

- Vacuum first-pass PYR4-PYR5 coupling: 11 meV
- Hydrated Frenkel-like PYR4-PYR5 coupling: 33 meV

The hydrated value should be treated as a flagged refinement, not yet as a final replacement, until additional hydrated dimers or transition-density/fragment-based coupling diagnostics are completed.
