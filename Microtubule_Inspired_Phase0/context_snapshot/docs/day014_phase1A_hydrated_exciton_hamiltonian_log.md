# Day014 Phase 1A Log — Hydrated Excitonic Hamiltonian

## Objective

Complete the first explicit-hydration excitonic Hamiltonian for the PYR2-PYR5 chromophore segment by computing and classifying hydrated nearest-neighbor dimer TDDFT/TDA results.

## Completed work

### PYR4-PYR5 hydrated dimer state characterization

The Day013 PYR4-PYR5 explicit-water TDDFT/TDA calculation was postprocessed. The anomalous low-energy S1 state at 2.959 eV was found to be dark and dominated by MO315a -> MO316a.

Orbital localization showed:

- MO315a localized on PYR4.
- MO316a localized on PYR5.

Therefore, S1 was interpreted as an inter-chromophore charge-transfer-like state, not as a water-localized artifact. It was excluded from Frenkel coupling estimation.

The bright Frenkel-like pair S2/S5 gave:

- S2 = 3.746 eV
- S5 = 3.812 eV
- J45 ≈ 33 meV

### PYR3-PYR4 hydrated dimer

A new explicit-water PYR3-PYR4 TDDFT/TDA calculation was completed using ORCA 6.1.1, PBE0/def2-SVP, TDDFT/TDA, 20 roots, and a 0.50 nm water shell.

The calculation terminated normally.

Key results:

- S1 = 2.888 eV, dark CT-like state.
- Bright states near 3.778-3.812 eV.
- Orbital localization confirmed Frenkel-like bright-state character.
- Effective hydrated coupling: J34 ≈ 17 meV.

### PYR2-PYR3 hydrated dimer

A new explicit-water PYR2-PYR3 TDDFT/TDA calculation was completed using ORCA 6.1.1, PBE0/def2-SVP, TDDFT/TDA, 20 roots, and a 0.50 nm water shell.

The calculation terminated normally.

Key results:

- S1 = 3.139 eV, dark CT-like state dominated by MO290a -> MO291a.
- Bright states include S2 = 3.754 eV, S4 = 3.774 eV, and S5 = 3.778 eV.
- Orbital localization showed alternating chromophore-localized MOs with low water leakage.
- Effective hydrated coupling: J23 ≈ 2.37 meV.
- This coupling is weak/mixed and should be treated as an operational Frenkel-like estimate, not as a clean two-state splitting.

### Fragment localization workflow

The cube-based fragment-localization script was corrected to avoid hardcoded cube prefixes and to operate on arbitrary hydrated dimer cube directories.

Final validated PYR2-PYR3 localization table:

| MO | PYR2 fraction | PYR3 fraction | Water fraction |
|---|---:|---:|---:|
| 286a | 0.949436 | 0.000000 | 0.050564 |
| 289a | 0.000000 | 0.976026 | 0.023974 |
| 290a | 0.973880 | 0.000000 | 0.026120 |
| 291a | 0.000000 | 0.972567 | 0.027433 |
| 292a | 0.977982 | 0.000000 | 0.022018 |
| 293a | 0.000000 | 0.972007 | 0.027993 |
| 295a | 0.977719 | 0.000000 | 0.022281 |
| 296a | 0.964729 | 0.000000 | 0.035271 |

## Final hydrated nearest-neighbor couplings

| Pair | J_eff (meV) | Assignment |
|---|---:|---|
| PYR2-PYR3 | 2.37 | weak/mixed Frenkel-like effective coupling |
| PYR3-PYR4 | 17.00 | bright Frenkel-like effective coupling |
| PYR4-PYR5 | 33.00 | bright Frenkel-like effective coupling |

## Full hydrated Hamiltonian

The full four-site hydrated Hamiltonian was built for PYR2-PYR5.

Hamiltonian, eV:

|      | PYR2 | PYR3 | PYR4 | PYR5 |
|---|---:|---:|---:|---:|
| PYR2 | 3.75500 | 0.00237 | 0.00000 | 0.00000 |
| PYR3 | 0.00237 | 3.77750 | 0.01700 | 0.00000 |
| PYR4 | 0.00000 | 0.01700 | 3.77750 | 0.03300 |
| PYR5 | 0.00000 | 0.00000 | 0.03300 | 3.77750 |

Exciton eigenstates:

| State | Energy (eV) | Dominant site | Participation ratio |
|---|---:|---|---:|
| X1 | 3.740338 | PYR4 | 2.417521 |
| X2 | 3.754835 | PYR2 | 1.023333 |
| X3 | 3.777696 | PYR3 | 1.521618 |
| X4 | 3.814631 | PYR4 | 2.399793 |

## Interpretation

Explicit hydration produces a strongly asymmetric Frenkel network. PYR4-PYR5 is the dominant hydrated coupling, PYR3-PYR4 is intermediate, and PYR2-PYR3 is weak/mixed.

The resulting Hamiltonian is not a homogeneous nearest-neighbor chain. One exciton remains highly localized on PYR2, while the remaining excitonic states are distributed primarily across PYR3-PYR5.

Dark CT-like states were observed in the hydrated dimers and excluded from the Frenkel Hamiltonian construction. This distinction is essential for avoiding overestimation or misassignment of excitonic couplings.

## Key outputs

- `runs/phase1A/day014_hydrated_hamiltonian_sensitivity/hydrated_neighbor_couplings_summary.csv`
- `runs/phase1A/day014_hydrated_hamiltonian_sensitivity/full_hydrated_hamiltonian_4x4_eV.csv`
- `runs/phase1A/day014_hydrated_hamiltonian_sensitivity/full_hydrated_hamiltonian_4x4_meV.csv`
- `runs/phase1A/day014_hydrated_hamiltonian_sensitivity/full_hydrated_exciton_eigenstates_day014.csv`
- `runs/phase1A/day014_hydrated_hamiltonian_sensitivity/full_hydrated_exciton_eigenvectors_day014.csv`
- `runs/phase1A/day014_hydrated_hamiltonian_sensitivity/FULL_HYDRATED_HAMILTONIAN_DAY014.md`
- `scripts/phase1A/analyze_cube_fragment_localization_radialmask.py`
