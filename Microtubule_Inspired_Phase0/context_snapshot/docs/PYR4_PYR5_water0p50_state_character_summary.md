# Day014 PYR4-PYR5 Hydrated Dimer State Character

## Source calculation

- System: PYR4-PYR5 hydrated dimer, local water shell 0.50 nm
- Method: ORCA 6.1.1, PBE0/def2-SVP, TDDFT/TDA
- Roots: 20
- Source output: Day013 hydrated dimer pilot
- Status: ORCA terminated normally

## Key result

The low-energy S1 state at 2.959 eV is dominated by a single orbital transition:

- S1: 315a -> 316a
- Weight: 0.999988
- Oscillator strength: 0.000000000

This state is dark and should not be interpreted directly as a normal bright pyrene excitonic splitting.

## Bright pyrene-like states

The first strong optical transitions are:

| State | Energy (eV) | fosc |
|---|---:|---:|
| S2 | 3.745830 | 0.371510312 |
| S5 | 3.812105 | 0.424868742 |

These are the relevant bright transitions for optical/excitonic interpretation.

## Technical interpretation

The hydrated dimer does not support a simple two-state splitting interpretation using S1/S2. The S1 state is an anomalous dark state and must be characterized further, likely by orbital visualization or fragment/localization analysis.

For coupling extraction, the immediate working hypothesis should use the bright pyrene-like manifold rather than the dark S1 state.

## Next step

Generate orbital or NTO visualizations for orbitals 315a and 316a to determine whether S1 is local, water-associated, inter-fragment charge-transfer-like, or an artifact of the hydrated cluster boundary.
