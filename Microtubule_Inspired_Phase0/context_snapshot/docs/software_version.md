# Software Versions — Phase 0

This file records the software versions available in the initial macOS environment used for Phase 0 setup, small/medium test systems, workflow development, and post-processing.

## System

- Working platform: macOS
- Reported by LAMMPS: Darwin 24.6.0 x86_64
- Working directory: `/Users/alejandro/projects/Microtubule_Inspired_Phase0`

Additional system information should be recorded with:

```bash
sw_vers
uname -m
sysctl -n machdep.cpu.brand_string
sysctl hw.memsize
echo $SHELL
```

## LAMMPS

### Executable

```bash
lmp_mpi
```

### Version

```text
LAMMPS version: 10 Dec 2025
Git info: release / patch_10Dec2025
Executable tested with: lmp_mpi -h
```

### Build information

```text
OS: Darwin 24.6.0 x86_64
Compiler: Clang C++ Apple LLVM 14.0.0 with OpenMP not enabled
C++ standard: C++17
MPI: Open MPI v5.0.7
FFT engine: mpiFFT
FFT library: KISS
Precision: double
```

### Accelerator configuration

```text
OPENMP package API: Serial
OPENMP package precision: double
GPU acceleration: not reported in this local build
```

### Relevant installed packages for this project

The local LAMMPS build includes several packages relevant to Phase 0 and later phases:

```text
MOLECULE
KSPACE
DIPOLE
DRUDE
DIELECTRIC
MANYBODY
MEAM
REAXFF
QEQ
RIGID
SPIN
PHONON
BROWNIAN
COLLOID
CORESHELL
CG-DNA
CG-SPICA
DPD
MOFFF
REACTION
FEP
MC
OPENMP
OPT
```

### Relevance to the project

The local LAMMPS installation is suitable for:

- inorganic nanotube or nanopore models;
- hybrid organic–inorganic structures;
- confined-water simulations;
- charged or dipolar systems;
- external electric-field tests;
- long-range electrostatics using KSPACE methods;
- many-body potentials for inorganic materials when appropriate;
- ReaxFF/QEQ exploratory reactive or charge-equilibration models;
- dipolar and dielectric-response prototypes;
- spin-related exploratory models if later required.

### Local limitation

This macOS LAMMPS build is appropriate for geometry construction, testing, debugging, and small/medium simulations. Large production simulations may require migration to a Linux workstation, HPC environment, or cloud compute resource with optimized MPI/OpenMP/GPU support.

## GROMACS

### Executable

```bash
gmx
```

### Version

```text
GROMACS version: 2025.4
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

### Build information

```text
C compiler: /usr/bin/gcc AppleClang 14.0.0.14000029
C++ compiler: /usr/bin/g++ AppleClang 14.0.0.14000029
Optimization: O3
BLAS: external
LAPACK: external
```

### Relevance to the project

The local GROMACS installation is suitable for:

- biomimetic or peptide-like tubular systems;
- water and ion simulations;
- organic/biomolecular force-field simulations;
- equilibration and thermal-stability tests;
- RMSD/RMSF/energy/hydrogen-bond analyses;
- trajectory generation for later dipole, hydration, and structural analysis.

### Local limitation

This GROMACS build has OpenMP support but no GPU support. It is suitable for local model construction, equilibration tests, and moderate CPU simulations. Large production runs may be more efficient on a Linux workstation or HPC/cloud resource with GPU acceleration.

## Python

```text
Python version: 3.12.12
pip version: 25.3
pip path: /Users/alejandro/miniforge3/lib/python3.12/site-packages/pip
Environment: base Miniforge environment
```

Python will be used for:

- geometry generation;
- file conversion;
- simulation setup automation;
- trajectory analysis;
- dipole/orientation analysis;
- autocorrelation functions;
- plotting;
- report-ready figure generation.

## Summary of local readiness

The macOS environment is ready for Phase 0 setup work. LAMMPS and GROMACS are already installed and accessible. The local machine is appropriate for:

- documentation;
- model preparation;
- small/medium test simulations;
- debugging;
- analysis;
- plotting;
- initial Phase 0 workflow development.

Large production simulations should remain eligible for migration to a Linux workstation or HPC/cloud environment depending on system size, number of temperatures, number of replicas, and trajectory length.
