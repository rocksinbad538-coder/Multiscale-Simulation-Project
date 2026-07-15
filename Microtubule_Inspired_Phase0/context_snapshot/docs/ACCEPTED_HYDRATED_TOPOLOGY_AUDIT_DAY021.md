# Day021 Accepted Hydrated Topology Audit

## Purpose

Audit the strongest accepted hydrated-topology candidate against the accepted GRO and TPR before any mobile-solute simulation is prepared.

## Accepted system

- GRO atoms: 68320.
- Accepted TPR atoms: 68320.
- HBN residues/atoms: 1/1680.
- PYR residues/atoms: 4/104.
- SOL residues/atoms: 16634/66536.
- GRO velocity records: 68320/68320.

## Candidate topology reconstruction

- Candidate topology: `parameters/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_validated/hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_hbnBonded_kang2000_improperGeo100.top`.
- Parsed topology atom total: 68320.
- Local topology files parsed: 4.
- Molecule types parsed: 3.
- Duplicate molecule types: 0.
- Missing molecule definitions: 0.
- `gmx grompp` return code: 0.
- Rebuilt TPR atoms: 68320.

## Validation

- Expected GRO composition: PASS.
- GRO versus accepted TPR atom count: PASS.
- Candidate topology total atom count: PASS.
- Topology molecule counts versus GRO: PASS.
- Topology atom totals versus GRO: PASS.
- No duplicate molecule types: PASS.
- No missing molecule types: PASS.
- Candidate grompp reconstruction: PASS.
- Rebuilt TPR atom count: PASS.

## Interpretation

The candidate topology is compositionally and structurally consistent with the accepted GRO and TPR and can be advanced to the parameter-level TPR identity audit.

The `gmx check` comparison is stored as a diagnostic because differences in coordinates, velocities, run state, or generated metadata do not by themselves prove a topology mismatch.
