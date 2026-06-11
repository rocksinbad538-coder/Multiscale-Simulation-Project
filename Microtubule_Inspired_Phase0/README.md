# Microtubule-Inspired Nanoscale Structure — Phase 0 Modeling

## Objective

This repository documents the computational setup, model construction, simulation workflow, and analysis protocol for Phase 0 of the microtubule-inspired nanoscale structure project.

Phase 0 compares three representative baseline architectures:

1. Organic / biomimetic tubular structure
2. Hybrid organic–inorganic tubular structure
3. Predominantly inorganic tubular structure

All models will be compared under controlled conditions: common geometry, confined water, identical temperature points, comparable excitation conditions, and common analysis metrics.

## Phase 0 Purpose

The goal of Phase 0 is to identify which baseline architecture maximizes:

- physical plausibility;
- computational tractability;
- structural stability;
- confined-water persistence;
- dipolar ordering;
- dielectric relaxation response;
- measurable output contrast relative to controls.

## Initial Computational Platform

The initial setup is developed on macOS. Larger production simulations may later be migrated to a Linux workstation, HPC system, or cloud compute environment if required.

## Repository Structure

```text
Microtubule_Inspired_Phase0/
├── README.md
├── docs/
├── latex/
│   └── Phase0_Modeling_Report.tex
├── references/
├── scripts/
│   ├── geometry/
│   ├── analysis/
│   └── plotting/
├── systems/
│   ├── organic_biomimetic/
│   ├── hybrid/
│   └── inorganic/
├── controls/
├── simulations/
├── results/
├── figures/
├── logs/
├── notes/
└── environment/
    ├── macos_environment_notes.md
    ├── software_versions.md
    └── requirements_phase0.txt
```

## Planned Software Stack

### Core tools

- Python: geometry generation, analysis, plotting, and workflow automation.
- NumPy / SciPy: numerical analysis.
- pandas: tabular data organization.
- matplotlib: plotting.
- MDAnalysis or MDTraj: trajectory analysis.
- LAMMPS: atomistic/coarse-grained molecular dynamics for inorganic, hybrid, and field-driven systems.
- GROMACS: biomolecular or peptide-like molecular dynamics if required.
- VMD / OVITO: molecular and trajectory visualization.

### Later-stage tools

- Quantum ESPRESSO: DFT calculations for representative fragments or periodic units.
- Octopus: TDDFT / optical-response calculations.
- COMSOL: continuum electromagnetic, THz, microwave, and device-level modeling.

## Verified Local Software Setup

The initial Phase 0 computational environment has been verified on macOS.

### LAMMPS

The available LAMMPS executable is:

```bash
lmp_mpi
```

The installed version is:

```text
LAMMPS 10 Dec 2025
MPI-enabled build using Open MPI v5.0.7
```

Important build information:

```text
OS: Darwin 24.6.0 x86_64
Compiler: Clang C++ Apple LLVM 14.0.0 with OpenMP not enabled
C++ standard: C++17
MPI: Open MPI v5.0.7
FFT engine: mpiFFT
FFT library: KISS
Precision: double
```

This build includes several packages relevant to the project, including:

```text
MOLECULE, KSPACE, AMOEBA, ASPHERE, BROWNIAN, CG-DNA, CG-SPICA,
COLLOID, CORESHELL, DIPOLE, DPD, DRUDE, DIELECTRIC, MANYBODY,
MC, MEAM, MOFFF, OPENMP, OPT, PHONON, QEQ, REACTION, REAXFF,
RIGID, SPIN
```

This is important because Phase 0 may require:

- nanotube or nanopore models;
- confined-water simulations;
- electrostatics;
- dipolar interactions;
- external electric fields;
- inorganic or hybrid potentials;
- later exploratory spin, dielectric, or phonon modules.

The current macOS LAMMPS build is suitable for local development, testing, small/medium simulations, and workflow debugging. Large production simulations may later be migrated to a Linux workstation, HPC system, or cloud compute resource.

### GROMACS

The available GROMACS executable is:

```bash
gmx
```

The installed version is:

```text
GROMACS 2025.4
Executable: /usr/local/gromacs/bin/gmx
Data prefix: /usr/local/gromacs
Precision: mixed
Memory model: 64 bit
MPI library: thread_mpi
OpenMP support: enabled
GPU support: disabled
SIMD instructions: AVX2_256
CPU FFT library: fftw-3.3.10-sse2
```

This build is suitable for biomimetic, peptide-like, water-containing, and soft-matter systems. It will be used if the organic/biomimetic Phase 0 candidate requires a biomolecular force-field workflow.

### Python

The available Python version is:

```text
Python 3.12.12
pip 25.3
pip path: /Users/alejandro/miniforge3/lib/python3.12/site-packages/pip
Environment: base Miniforge environment
```

Python will be used for:

- geometry generation;
- workflow automation;
- trajectory analysis;
- autocorrelation functions;
- dipole orientation metrics;
- plotting;
- report-ready figure generation.

## Current Computational Role of macOS

The macOS environment will be used for:

- documentation;
- geometry generation;
- candidate model construction;
- small/medium simulation tests;
- local equilibration checks;
- analysis and plotting;
- preparation of inputs for possible Linux/HPC production runs.

Large-scale production simulations will be migrated if required by system size, number of replicas, simulation length, or computational cost.

## Installation Notes for Reproducibility

LAMMPS and GROMACS are already installed on the current macOS workstation. The following notes are included for reproducibility on future systems.

### LAMMPS

LAMMPS is selected because Phase 0 includes inorganic, hybrid, coarse-grained, and field-driven systems. It is suitable for nanotubes, confined fluids, polarizable or charged models, surface functionalization, and external-field simulations.

One possible macOS installation route is:

```bash
brew install lammps
```

Alternative builds may be required later if specific optimized packages, GPU support, or HPC-specific MPI/OpenMP settings are required.

### GROMACS

GROMACS is selected because it is highly optimized for biomolecular and soft-matter molecular dynamics, including peptide-like systems, water, ions, force-field-based simulations, trajectory analysis, and standard observables such as RMSD, RMSF, energy evolution, hydrogen bonding, and dipole-related analysis.

One possible macOS installation route is:

```bash
brew install gromacs
```

For large production simulations, optimized Linux builds with MPI, OpenMP, SIMD, and GPU acceleration may be preferred.

### Python environment

A minimal Python environment can be prepared using:

```bash
python3 -m venv phase0_env
source phase0_env/bin/activate
pip install --upgrade pip
pip install -r environment/requirements_phase0.txt
```

Alternatively, because the current system uses Miniforge, a conda environment can be created later if dependency isolation becomes necessary.

### Visualization

VMD and OVITO are useful for visual inspection of structures and trajectories. They are not required for automated analysis, but they are important for sanity checks, screenshots, and qualitative review.

## Reproducibility Principles

Each simulation must include:

- input files;
- generation scripts;
- software version information;
- run commands;
- log files;
- trajectory outputs;
- analysis scripts;
- processed data;
- figures;
- short interpretation notes.

No result should be reported without a reproducible path from input model to final figure or table.

## Phase 0 Immediate Next Steps

1. Finish documenting the local computational environment.
2. Define the initial candidate modeling matrix.
3. Select the first three representative architectures under Option B.
4. Define common geometry, confined-water setup, temperatures, controls, and first excitation conditions.
5. Build the first minimal geometry-generation scripts.

## Python Environment for Phase 0

A dedicated conda environment is used for Phase 0 analysis and workflow automation.

### Create environment

```bash
conda create -n mt_phase0 python=3.12 -y
conda activate mt_phase0
```

### Install core packages

```bash
conda install -c conda-forge numpy scipy pandas matplotlib jupyterlab ipykernel tqdm -y
conda install -c conda-forge mdanalysis mdtraj nglview -y
```

### Verify Installation

```bash
python -c "import numpy, scipy, pandas, matplotlib, MDAnalysis, mdtraj; print('Python analysis stack OK')"
```

### Register Jupyter kernel

```bash
python -m ipykernel install --user --name mt_phase0 --display-name "Python (mt_phase0)"
```

### Export reproducible environment

```bash
conda env export --no-builds > environment/mt_phase0_environment.yml
pip freeze > environment/mt_phase0_pip_freeze.txt
```

The exported environment files are stored in:

```text
environment/mt_phase0_environment.yml
environment/mt_phase0_pip_freeze.txt
```

## Generic Reference Geometry

A first generic tubular scaffold was generated to validate the common Phase 0 geometry before assigning a specific material model.

### Geometry parameters

```text
Outer diameter: 24.0 nm
Lumen diameter: 14.0 nm
Wall thickness: 5.0 nm
Segment length: 20.0 nm
Tube axis: z
```

### Generated File

```text
systems/inorganic/geometry/generic_tube_24OD_14ID_20L.xyz
```

### Inspection outputs

```text
results/geometry/generic_tube_24OD_14ID_20L_summary.csv
figures/geometry/generic_tube_24OD_14ID_20L_cross_section.png
```

### Verification
The geometry inspection confirmed:

```text
Estimated lumen diameter: 14.0 nm
Estimated outer diameter: 24.0 nm
Estimated length: 20.0 nm
Number of shell points: 17,280
```

This generic scaffold is not yet a force-field-ready molecular model. It is a reference geometric scaffold used to verify dimensions, coordinate conventions, plotting, and documentation workflow before constructing the physical Phase 0 candidates.

## Inorganic Baseline Decision

The first predominantly inorganic Phase 0 candidate will be treated as a BNNT / BN-like tubular scaffold with confined water.

A key modeling issue was identified: a single-wall BNNT cannot directly reproduce the approved common Phase 0 geometry of 24 nm outer diameter, 14 nm lumen diameter, and 5 nm wall thickness. A single-wall BNNT is effectively a rolled monolayer of h-BN, whereas the Phase 0 geometry requires a finite-thickness wall.

Therefore, the first inorganic baseline will be handled as a geometry-equivalent BN-like tubular scaffold for controlled screening. This preserves the common Phase 0 geometry and allows direct comparison against organic/biomimetic and hybrid candidates.

The detailed rationale is documented in:

```text
docs/phase0_inorganic_baseline_rationale.md
```

If the inorganic baseline shows strong measurable response, a later refinement can use explicit single-wall or multi-wall BNNT models, DFT-informed partial charges, or BN-specific force fields.

---

## Current Phase 0 Status: Inorganic BN-like Baseline

The first Phase 0 baseline system is an inorganic BN-like tubular scaffold with confined water. The initial target geometry follows the common Phase 0 dimensions:

- Outer diameter: 24 nm
- Lumen diameter: 14 nm
- Segment length: 20 nm
- Scaffold representation: BN-like alternating charged scaffold
- Confined medium: explicit TIP3P-like water

The locally tested production candidate contains:

- 17,280 scaffold atoms
- 30,000 confined water molecules
- 90,000 water atoms
- 107,280 total atoms

Scaling tests were performed up to 60,000 confined water molecules, corresponding to 197,280 total atoms. The 60,000-water system successfully passed LAMMPS read/initialization, but the 30,000-water system was selected for local relaxation and NVT testing because it is more practical on the current Mac workstation.

### Water Placement Strategy

An initial random insertion strategy was tested for confined water placement. It worked for smaller systems but became inefficient at high hydration levels because random sequential insertion becomes increasingly slow near high packing density.

A lattice-based confined-water generator was then implemented to efficiently generate larger confined-water systems. The lattice-based placement is used only to create an initial configuration. Subsequent minimization and finite-temperature MD are required to relax the artificial initial placement.

### Relaxation and Thermal Testing

The 30,000-water BN-like scaffold system has completed:

1. LAMMPS read/initialization validation.
2. Initial minimization.
3. Continued minimization.
4. Controlled low-temperature dynamic relaxation.
5. Post-relaxation FIRE minimization.
6. Corrected NVT test at 50 K.
7. NVT ramp from 50 K to 150 K.
8. NVT ramp from 150 K to 200 K.

The corrected NVT inputs report both the global default temperature and the water-group temperature. The water-group temperature, `c_twater`, is the physically relevant value because the scaffold is fixed during these initial tests.

### Stability Criteria Observed

The NVT tests completed with:

- no NaNs;
- no lost atoms;
- no dangerous neighbor builds;
- stable water-group temperature control;
- smoothly increasing confined-water MSD.

The current baseline is therefore suitable for continued short-temperature ramp testing to 250 K and 300 K, followed by trajectory analysis.

### Visualization

VMD-compatible XYZ files and LAMMPS trajectories are generated locally for qualitative inspection. Visualization is used to confirm scaffold geometry, water confinement, and overall structural integrity. Quantitative conclusions will be based on LAMMPS logs and trajectory analysis, not visual inspection alone.


## Thermal ramp and 300 K hold

The BN-like scaffold + 30,000 confined water system was advanced through a staged NVT thermal ramp and a short 300 K hold.

Completed stages:

- 50 K corrected NVT
- 50 → 150 K ramp
- 150 → 200 K ramp
- 200 → 250 K ramp
- 250 → 300 K ramp
- 300 K constant-temperature hold

Final 300 K hold metrics:

- Water-specific temperature: 299.74 K
- Water MSD: 12.349 Å²
- Dangerous builds: 0
- No NaNs detected
- No lost atoms reported
- Final fraction of water oxygen atoms inside the nominal lumen segment: 0.9983
- Final fraction radially outside the nominal lumen radius: 0.00163
- Final fraction axially outside the nominal segment: 0.000133
- Final fraction outside the outer scaffold radius: 0.0

Generated outputs include thermodynamic summaries, confinement summaries, updated figures, LAMMPS trajectories, and final data files. These results validate short-timescale numerical stability and confinement for the Phase 0 inorganic baseline. The force field remains provisional and should not yet be interpreted as a final predictive material model.

## Phase 0 status — Day 004 update

Day 004 completed the BN-like inorganic baseline controls and short field-response diagnostics.

Completed items:

- Dry BN-like scaffold 300 K tethered control.
- Hydrated BN-like scaffold-water extended contained 300 K hold.
- Water-only contained reference generated from the hydrated state.
- Scaffold-water vs water-only confinement comparison.
- Water dipole orientation and dipole autocorrelation diagnostics.
- Short axial electric-field response test at Ez = 0.01 V/Å for both scaffold-water and water-only systems.

Main result:

The BN-like scaffold preserves radial confinement of the water phase relative to the water-only contained reference. Water dipoles show near-isotropic average orientation without field, but strong axial alignment under the applied field. Under the current short-timescale provisional-force-field diagnostics, the scaffold contribution is clearer for radial confinement than for average axial dipolar ordering.

Key figures:

- `figures/day004/day004_bn_like_hydrated_and_dry_controls_dashboard.png`
- `figures/day004/day004_scaffold_water_vs_water_only_control_comparison.png`
- `figures/day004/day004_water_dipole_autocorrelation_comparison.png`
- `figures/day004/day004_fieldZ_response_scaffold_water_vs_water_only.png`

Detailed notes:

- `notes/day_003.md`
- `notes/day_004.md`

## Day 005 update: charge/polarity screening and field-response controls

Day 005 completed a matched Phase 0 charge/polarity-control matrix for the inorganic scaffold-water branch.

Starting from the validated BN-like polar scaffold-water baseline, two matched scaffold controls were generated and simulated:

- BN-neutralized scaffold-water: same BN-like geometry and B/N atom-type topology, but scaffold charges set to zero.
- carbon-like neutral scaffold-water: same tubular geometry and water configuration, but scaffold represented as a neutral carbon-like placeholder.
- water-only contained reference: same confined water loading without scaffold.

The field-free controls were compared at 20,000 steps. The BN-like polar scaffold maintained strong radial confinement, with final maximum radial water oxygen position near 73.34 Å. In contrast, BN-neutralized and carbon-like neutral controls reached approximately 95.25 Å and 97.12 Å, respectively. This indicates that the BN-like polar wall is not acting only as a geometric tube; its charge/polarity pattern is a key confinement variable.

A corresponding axial field diagnostic was completed for all four systems using \(E_z = 0.01\) V/Å. Under fieldZ, all systems showed similar axial water-dipole alignment, with final mean cos(theta_z) near 0.59 and final S_z near 0.26. However, radial confinement remained strongly architecture-dependent: BN-like polar stayed near 72.82 Å, while BN-neutralized and carbon-like neutral reached approximately 105.12 Å and 107.78 Å.

Main conclusion:

The applied axial field controls average axial water-dipole alignment, while the BN-like scaffold charge/polarity pattern controls radial confinement. The BN-like polar scaffold is therefore the strongest current inorganic Phase 0 baseline.

Key outputs:

- `results/phase0/day005_charge_polarity_screening_summary.csv`
- `figures/day005/day005_charge_polarity_controls_20k_comparison.png`
- `figures/day005/day005_fieldZ_charge_polarity_controls_comparison.png`
- `latex/Phase0_Modeling_Report.pdf`
- `latex/Phase0_Modeling_Report.tex`


## Day 006 update: BN-like polar long-stability and extended field-response validation

Day 006 extended the strongest Phase 0 inorganic baseline identified during Day 005: the BN-like polar scaffold-water system.

Two longer validation simulations were completed:

1. Field-free BN-like polar scaffold-water extension from 20k to 50k total steps.
2. BN-like polar scaffold-water fieldZ extension from 10k to 30k total fieldZ steps.

Both simulations completed cleanly with:

- 107,280 total atoms
- 17,280 scaffold atoms
- 30,000 water molecules
- timestep: 0.10 fs
- dangerous builds: 0
- no detected NaNs, lost atoms, or missing bonds

Field-free 50k result:

The BN-like polar scaffold-water system remained radially confined through 50,000 total steps. The final maximum radial water oxygen position was approximately 73.25 Å, with a final fraction inside the nominal lumen segment of approximately 0.99613. Field-free water dipole orientation remained approximately isotropic, with mean cos(theta_z) near zero and S_z near zero.

FieldZ 30k result:

The BN-like polar scaffold-water system also remained radially confined under prolonged axial field exposure. At 30,000 fieldZ steps, the final maximum radial water oxygen position was approximately 74.25 Å, with a final fraction inside the nominal lumen segment of approximately 0.99597. The axial water-dipole response strengthened substantially under prolonged field exposure, reaching mean cos(theta_z) of approximately 0.763 and S_z of approximately 0.471.

Main conclusion:

The BN-like polar scaffold-water baseline preserves radial confinement both without field and under prolonged axial field exposure. Under fieldZ, the system develops strong axial water-dipole ordering without radial destabilization over the tested timescale.

Key Day 006 outputs:

- `systems/inorganic/bn_like_scaffold_water/outputs/bn_like_scaffold_water_30000w_nvt300_hold_contained_50k.data`
- `systems/inorganic/bn_like_scaffold_water/outputs/nvt300_hold_contained_bn_like_extend_20k_to_50k_30000w.lammpstrj`
- `systems/inorganic/bn_like_scaffold_water/outputs/bn_like_scaffold_water_30000w_nvt300_fieldZ_contained_30k.data`
- `systems/inorganic/bn_like_scaffold_water/outputs/nvt300_fieldZ_contained_bn_like_extend_10k_to_30k_30000w.lammpstrj`
- `figures/day006/day006_bn_like_50k_stability_vs_day005_controls.png`
- `figures/day006/day006_bn_like_fieldZ_30k_response.png`
- `notes/day_006.md`


## Day 007 update: first lead hybrid/chromophore-bearing BN-like candidate

Day 007 transitioned the project from the validated BN-like polar scaffold-water baseline toward a first hybrid/chromophore-bearing embodiment.

A carved 12-dipole hybrid candidate was built from the validated BN-like polar 50k field-free state. The model adds 12 fixed chromophore-like dipolar pseudo-site pairs near the inner lumen wall. To avoid initial overlap artifacts, nearby water molecules were locally removed around the pseudo-sites.

Lead hybrid model:

- BN-like polar tubular scaffold
- confined water
- 12 fixed chromophore-like dipolar pseudo-site pairs
- 24 chromophore pseudo-atoms
- chromophore charge magnitude: +/- 0.10 e
- water molecules after carving: 29,987
- total atoms: 107,265

The initial uncarved model produced overlap-driven instability, but the carved model corrected this and completed stable field-free and fieldZ tests.

Completed Day 007 simulations:

1. carved hybrid field-free 5k stability test
2. carved hybrid field-free extension to 20k
3. carved hybrid fieldZ test to 10k
4. carved hybrid fieldZ extension to 20k
5. carved hybrid fieldZ extension to 30k

All successful Day 007 production runs completed with dangerous builds = 0 and no detected lost atoms, NaNs, or missing bonds.

Final field-free 20k result:

- maximum radial water oxygen position: approximately 72.80 Å
- fraction inside nominal lumen segment: approximately 0.99596
- S_z: approximately -0.00991

Final fieldZ 30k result:

- maximum radial water oxygen position: approximately 73.67 Å
- fraction inside nominal lumen segment: approximately 0.99586
- mean cos(theta_z): approximately 0.662
- S_z: approximately 0.342

Main conclusion:

The carved 12-dipole chromophore-bearing BN-like scaffold-water model preserves radial confinement under both field-free and fieldZ conditions while producing a measurable, progressively increasing axial water-dipole response under fieldZ. This model is the current lead Phase 0 hybrid embodiment candidate.

Key Day 007 outputs:

- `systems/hybrid/bn_like_chromophore_scaffold_water/outputs/bn_like_chromophore_12dipoles_carved_nvt300_contained_20k.data`
- `systems/hybrid/bn_like_chromophore_scaffold_water/outputs/bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_30k.data`
- `results/phase0/day007_lead_hybrid_candidate_summary.csv`
- `figures/day007/day007_hybrid_carved_20k_vs_bn_like_baseline.png`
- `figures/day007/day007_hybrid_carved_fieldfree_and_fieldZ_to_30k.png`
- `notes/day_007.md`

