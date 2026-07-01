# Day020 Technical Closeout

- Date: 2026-07-01
- Project day: Day020
- Workstreams: excitonic open-system dynamics and confined-water analysis.

## 1. Combined excitonic dynamics

The coherent, pure-dephasing, and detailed-balance relaxation model was completed and numerically validated over 24 parameter conditions and 21 Hamiltonian snapshots.

- PYR4 remains the dominant direct gateway into PYR5, contributing 84.3–91.4% of the mean downward flux.
- Strong dephasing raises the transient PYR5 population at 100 ps to approximately 0.105–0.110.
- Thermal relaxation alone approaches the Gibbs sink at PYR5, whereas strong local dephasing shifts the combined stationary state toward the uniform four-state limit.
- The absolute dissipative rates remain phenomenological because no microscopic spectral density is available.

## 2. Frozen-solute MD audit

- Accepted trajectory: 201 frames over 100 ps with 0.5 ps saved-frame spacing.
- System: 68,320 atoms, including 16,634 TIP4P/2005 waters.
- HBN and all four PYR residues were frozen in all Cartesian directions.
- The trajectory supports solvent-structure and electrostatic-disorder analysis, but not solute thermal-stability metrics.

## 3. HBN architecture

- HBN atoms: 1680.
- Axial planes: 56.
- Typical plane spacing: 0.073000 nm.
- Detected continuous segments: 1.
- Detected axial gaps: 0.
- Mean wall radius: 1.199126 nm.

## 4. Confined-water density

- Integrated water count in the analyzed cylinder: 5211.768138.
- Radial depletion minimum: 1.174620 nm.
- Effective radial boundaries: 0.923345 to 1.453565 nm.
- Left mouth transition: -3.112007 to -2.198479 nm.
- Right mouth transition: 2.519318 to 3.374122 nm.

### Profile-guided populations

- Lumen core: 393.289350 waters (7.546%).
- HBN interfacial shell: 488.914171 waters (9.381%).
- Mouth transitions: 105.666141 waters (2.027%).
- Exterior solvent: 4223.898476 waters (81.045%).

- Tube-associated population: 987.869662 waters (18.955%).
- Conservation error: 2.291e-07.

## 5. Evidence classification

### Accepted

- Combined open-system mechanism and publication-quality figure set.
- One continuous HBN segment with four axially embedded PYR residues.
- Validated axial–radial water-density map and cylindrical integration.
- Predominantly exterior solvent and a hydrated lumen with a radial depletion shell at the HBN wall.

### Conditional or descriptive

- Exact interfacial occupancy because it retains moderate boundary sensitivity.
- Exact mouth occupancy and left–right asymmetry because mouth sensitivity remains high.

### Not supported by the accepted trajectory

- HBN or PYR RMSD/RMSF as thermal-stability metrics.
- Coupled water–solute conformational dynamics.
- Converged residence times or long-time diffusion.
- Microscopic bath spectral density.
