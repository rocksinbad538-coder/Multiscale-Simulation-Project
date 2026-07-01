# Day020 Density DAT Parser Repair

## Identified DAT structure

The GROMACS textual density output contains coordinate headers:

- Row 0, columns 1 onward: radial-bin centers.
- Column 0, rows 1 onward: axial-bin centers.
- Rows 1 onward, columns 1 onward: physical water-oxygen number density.

The negative entries detected by the original validator were axial coordinates at negative z, not negative densities.

## Canonical grid

- Raw DAT dimensions: 161 × 55.
- Density dimensions: 160 × 54.
- Axial centers: -4.009000 to 3.958890 nm.
- Radial centers: 0.000000 to 2.649140 nm.
- Median axial spacing: 0.05011000 nm.
- Median radial spacing: 0.04998400 nm.
- Minimum density: 0.000000 nm^-3.
- Maximum density: 88.541600 nm^-3.

## Cross-checks against analysis inputs

- Requested bin width: 0.05 nm.
- Requested axial half-range: 4.009 nm.
- Requested radial maximum: 2.699125700748155 nm.
- HBN axial span: 6.018000000000001 nm.
- Mean HBN wall radius: 1.1991257007481553 nm.

## Canonical product

- `runs/phase1A/day020_confined_water_axial_radial_density/water_oxygen_axial_radial_density_canonical.npz`

The NPZ file stores axial and radial centers, bin edges, and the validated 160 × 54 density matrix. It is the canonical input for subsequent plotting and regional solvent classification.

## Status

The existing density calculation is numerically valid. No GROMACS rerun and no XPM conversion are required.
