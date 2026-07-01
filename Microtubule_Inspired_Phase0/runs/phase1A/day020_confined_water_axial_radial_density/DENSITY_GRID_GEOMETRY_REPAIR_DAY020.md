# Day020 Density Grid Geometry Repair

## Correct coordinate semantics

The first DAT column and row contain lower cell edges, not cell centers. This follows from their exact agreement with the requested domain boundaries:

- Axial start: -4.009000 nm.
- Requested axial domain: [-4.009000, +4.009000] nm.
- Radial start: 0.000000 nm.
- Requested radial domain: [0, 2.699126] nm.

The physical bin centers were therefore reconstructed from consecutive edges. The radial grid now starts at exactly zero and contains no artificial negative radial boundary.

## Corrected grid

- Density cells: 160 × 54.
- Axial edges: -4.009000 to 4.009000 nm.
- Radial edges: 0.000000 to 2.699126 nm.
- Median axial spacing: 0.05011000 nm.
- Median radial spacing: 0.04998400 nm.

## Cylindrical integration

Each grid cell was assigned the exact cylindrical volume

\[
\Delta V_{ij}=\pi\left(r_{j+1}^{2}-r_j^{2}\right)\left(z_{i+1}-z_i\right).
\]

- Analysis-domain volume: 183.511018523 nm³.
- Integrated average water count in the analyzed cylinder: 5211.768137771.
- Volume-weighted mean density: 28.400300863 nm⁻³.

## Product status

- Current canonical product: `runs/phase1A/day020_confined_water_axial_radial_density/water_oxygen_axial_radial_density_canonical_v2.npz`.
- Previous product retained for traceability: `runs/phase1A/day020_confined_water_axial_radial_density/water_oxygen_axial_radial_density_canonical.npz`.

The V2 product must be used for all figures and regional cylindrical integrations.
