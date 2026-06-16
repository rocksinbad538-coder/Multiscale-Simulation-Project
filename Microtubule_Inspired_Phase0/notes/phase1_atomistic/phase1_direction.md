# Phase 1 direction: atomistic chromophore-bearing model

The coarse-grained BN-like chain scaffold used in Phase 0 successfully demonstrated stable confined-water behavior and reproducible field-induced water polarization. However, topology analysis showed that the scaffold is a tubular array of independent BN-like chains rather than an atomistic h-BN network. Therefore, it is not suitable for direct DFT/TDDFT or chemical passivation.

To fulfill the master project plan, Phase 1 will move to a chemically valid atomistic model:

- atomistic h-BN nanotube or chemically defined tubular scaffold
- real chromophores instead of pseudo-dipoles
- explicit confined water
- classical MD for thermal and hydration stability
- representative fragment extraction
- DFT/TDDFT parameterization
- excitonic Hamiltonian construction
- open-system excitation dynamics

The Phase 0 model is retained as a screening and methodology-validation result, not as the final electronic-structure model.

## Phase 1A.1 chromophore decision

Pyrene is selected as the first methodological chromophore for Phase 1A.

This decision is made to validate the DFT/TDDFT → site-energy/transition-dipole → excitonic-Hamiltonian pipeline using a small, rigid, aromatic, computationally tractable chromophore.

Pyrene is not declared to be the final biological chromophore of the project. Indole/tryptophan-like chromophores remain planned as later biomimetic controls after the pyrene-based workflow is validated.

Immediate next step:

- generate standalone pyrene geometry;
- audit pyrene connectivity and bond lengths;
- place four pyrenes around the atomistic h-BN nanotube;
- audit chromophore–chromophore and chromophore–scaffold distances.
