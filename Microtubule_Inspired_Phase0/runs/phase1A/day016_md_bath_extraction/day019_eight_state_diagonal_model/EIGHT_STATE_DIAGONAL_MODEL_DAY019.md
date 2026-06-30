# Day019 character-indexed eight-state diagonal model

## Basis definition

| Index | Label | Site | Family | TDDFT root |
|---:|---|---|---|---:|
| 0 | `PYR2_alternate` | PYR2 | alternate_like | S1 |
| 1 | `PYR2_bright` | PYR2 | bright_like | S2 |
| 2 | `PYR3_alternate` | PYR3 | alternate_like | S1 |
| 3 | `PYR3_bright` | PYR3 | bright_like | S2 |
| 4 | `PYR4_alternate` | PYR4 | alternate_like | S1 |
| 5 | `PYR4_bright` | PYR4 | bright_like | S2 |
| 6 | `PYR5_alternate` | PYR5 | alternate_like | S2 |
| 7 | `PYR5_bright` | PYR5 | bright_like | S1 |

The labels follow electronic character, not a globally fixed root number. PYR5 therefore uses S1 as the bright-like state, whereas PYR2-PYR4 use S2.

## Validation

- Frame/site rows consumed: 84/84
- Character-indexed state observations: 168/168
- Frames represented: 21/21
- Frame spacing: 5.000000 ps
- Bright-root mapping preserved: 84/84
- State-statistics energy SD values were validated and exported in meV.
- Hamiltonian snapshots are diagonal-only; all off-diagonal couplings are intentionally zero placeholders.

## State statistics

| State | Mean energy (eV) | SD (meV) | Range (meV) | Mean fosc |
|---|---:|---:|---:|---:|
| `PYR2_alternate` | 4.006143 | 8.061 | 27.000 | 0.011173 |
| `PYR2_bright` | 4.088571 | 10.817 | 36.000 | 0.493290 |
| `PYR3_alternate` | 4.003190 | 13.037 | 50.000 | 0.009893 |
| `PYR3_bright` | 4.078048 | 11.615 | 44.000 | 0.501778 |
| `PYR4_alternate` | 4.019238 | 11.820 | 57.000 | 0.008602 |
| `PYR4_bright` | 4.090619 | 7.403 | 35.000 | 0.510859 |
| `PYR5_alternate` | 3.847810 | 14.949 | 56.000 | 0.018624 |
| `PYR5_bright` | 3.776952 | 16.034 | 47.000 | 0.648089 |

## Aggregate fluctuation structure

- Bright-like energy SD range: 7.403-16.034 meV
- Alternate-like energy SD range: 8.061-14.949 meV
- Maximum absolute off-diagonal energy correlation: 0.730259
- PYR2 alternate/bright fluctuation correlation: 0.582864
- PYR3 alternate/bright fluctuation correlation: 0.605314
- PYR4 alternate/bright fluctuation correlation: 0.730259
- PYR5 alternate/bright fluctuation correlation: 0.717122

## Hamiltonian status

The resulting matrices provide the complete time-dependent diagonal component of the eight-state Hamiltonian. They are not yet the final excitonic Hamiltonian because intersite couplings and possible same-site interstate couplings have not been computed. No dynamical propagation should treat the current zero off-diagonal entries as physical coupling estimates.
