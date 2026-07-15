# Day020 MD and Confined-Water Input Audit

## Accepted run

- Run directory: `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute`.
- Files discovered: 9.
- Trajectory files: 1.
- TPR files: 1.
- GRO files: 1.
- MDP files: 1.
- NDX files: 0.
- TOP/ITP files: 0.

## Primary files selected

- Trajectory: `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc`.
- Run input: `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.tpr`.
- Structure: `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.gro`.

## GROMACS inspection

- Executable: `/usr/local/gromacs/bin/gmx`.
- GROMACS available: True.
- Last reported frame: 200.
- First trajectory time: 0.0 ps.
- Last trajectory time: 100.0 ps.

## Structure composition

- Declared atoms in selected GRO: 68320.
- Residues in selected GRO: 16639.
- Water residues: 16634.
- Water atoms: 66536.
- Ions: 0.

## Frozen-solute evidence

- Freeze-related settings found: True.
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.mdp`: `freezegrps  = HBN PYR`
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.mdp`: `freezedim   = Y Y Y Y Y Y`

## Scientifically valid analyses for the accepted trajectory

- Confined-water density and spatial heterogeneity.
- Water orientation relative to the tube axis and local chromophore geometry.
- Water–solute contacts and hydrogen-bond occupancy.
- Snapshot-resolved electrostatic environment and disorder.
- Short-time water correlation functions, provided their sampling limitations are reported explicitly.

## Analyses not supported by the frozen-solute trajectory

- Solute RMSD or RMSF as thermal-stability metrics.
- Scaffold or chromophore conformational stability.
- Coupled water–solute structural dynamics.
- Converged long-time diffusion or residence times from the current 100 ps window.
- A microscopic spectral density derived from coupled structural fluctuations.

## Required next decision

Use the audited groups and topology to define the confined-water structural analysis, then construct a controlled restraint-release sequence for a mobile-solute trajectory. The mobile trajectory must be validated before RMSD, RMSF, structural stability, or coupled bath dynamics are reported.
