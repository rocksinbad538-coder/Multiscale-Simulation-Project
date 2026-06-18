# Phase 1A Accepted Scaffold-Mechanics Baseline

## Accepted model

- Scaffold: bonded h-BN nanotube
- Angular stiffness model: kang2000
- Geometry-referenced improper stiffness: improperGeo100
- Validation run: NVT, 100 ps, 100 K
- Validation protocol: PYR and SOL frozen; h-BN scaffold mobile
- Purpose: isolate and validate h-BN scaffold mechanics before moving to controlled hydration and electronic-structure stages

## Acceptance decision

The model `kang2000 + improperGeo100` is accepted as the current Phase 1A scaffold-mechanics baseline.

## Key validation results

The 100 ps / 100 K NVT trajectory completed successfully without crash or numerical instability.

Final PYR-HBN minimum distances:

- PYR2: 0.38235 nm
- PYR3: 0.32135 nm
- PYR4: 0.42663 nm
- PYR5: 0.70072 nm

Contact criterion:

- No PYR-HBN contacts below 0.30 nm were observed.
- The closest final contact was PYR3-HBN at 0.32135 nm.
- This is interpreted as localized adsorption-like PYR-HBN contact rather than scaffold collapse.

Radial deformation after 100 ps:

- min Δr: -0.17150 nm
- max Δr: +0.23384 nm
- mean Δr: +0.01277 nm

Energy stability:

- Potential average: -905721 kJ/mol
- Potential RMSD: 73.37 kJ/mol
- Potential total drift over 100 ps: 39.03 kJ/mol
- Relative drift: approximately 4.3e-5 with respect to the absolute potential energy

## Interpretation

The accepted scaffold does not show global nanotube collapse, uncontrolled radial deformation, or chemically unrealistic PYR-HBN overlap during the 100 ps / 100 K validation.

The remaining close PYR-HBN interaction is localized and remains above the 0.30 nm rejection threshold. Therefore, it is treated as a physically plausible adsorption-like contact at this stage, not as a mechanical failure of the h-BN scaffold.

## Files

Accepted structure and validation outputs:

- `minimized_start.gro`
- `nvt_100ps_100K_final.gro`
- `nvt_100ps_100K.log`
- `nvt_100ps_100K.mdp`
- `audit_pyrene_hbn_nvt_100ps_100K_improperGeo100_NOPBC.csv`
- `potential_100ps.xvg`

Accepted topology/parameter files are stored in:

- `parameters/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_validated/`

## Operational conclusion

This model is the accepted Phase 1A scaffold-mechanics baseline.

Next stage:

Proceed to controlled hydration / solvent-validation using the accepted h-BN scaffold model.
