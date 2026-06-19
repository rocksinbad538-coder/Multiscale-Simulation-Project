# Phase 1A Hydrated Baseline Candidate

Accepted hydrated baseline candidate:

- Scaffold: h-BN bonded scaffold
- Bonded model: kang2000 + geometry-referenced improperGeo100
- Scaffold/chromophores: frozen during water relaxation
- Water model: TIP4P/2005
- Hydration protocol:
  - strict water-only minimization
  - 2 ps / 50 K water relaxation
  - 5 ps / 75 K water relaxation
  - 10 ps / 100 K water relaxation
  - post-minimization

## Key results

Internal confined water:

- Strict internal OW after post-min: 423
- Mean-radius internal OW after post-min: 423
- Initial accepted scaffold internal OW: 435

Thus, the confined-water population is retained during staged solvent relaxation.

Spatial distribution:

- Axial internal OW profile: 38, 68, 69, 69, 69, 73, 37
- The central tube region remains uniformly hydrated.
- Lower end-bin counts are consistent with boundary/lumen cutoff effects.

OW–OW nearest-neighbor statistics:

- internal OW count: 423
- min OW–OW nearest-neighbor distance: 0.26757 nm
- mean OW–OW nearest-neighbor distance: 0.27392 nm
- max OW–OW nearest-neighbor distance: 0.28157 nm

Contact audit after post-minimization:

- min OW–HBN: 0.28139 nm
- min OW–PYR: 0.27355 nm
- min OW–solute: 0.27355 nm
- min heavy-water / heavy-solute: 0.28139 nm
- min all-atom water / solute: 0.20752 nm

Interpretation:

The remaining shortest all-atom contact is a water-hydrogen contact, not a heavy-atom clash. Heavy-atom water–solute contacts are acceptable.

## Operational conclusion

This structure is accepted as the current Phase 1A hydrated baseline candidate.

It is suitable for subsequent hydration analysis, RDF analysis, chromophore-environment inspection, and preparation of downstream electronic-structure workflows.

## Temporal confined-water stability

A frame-by-frame audit was performed over the 10 ps / 100 K water-relaxation trajectory.

Internal OW count over trajectory:

- minimum: 423
- maximum: 431
- mean: 427.10

This indicates that the confined-water population remains temporally stable during the staged relaxation protocol, with no evidence of progressive lumen drainage.

## Local pyrene hydration audit

A per-chromophore hydration audit was performed on the accepted hydrated baseline.

Per-pyrene OW counts:

- PYR2: 4 OW within 0.35 nm, 19 within 0.50 nm, 74 within 0.75 nm; min OW distance 0.27355 nm
- PYR3: 5 OW within 0.35 nm, 18 within 0.50 nm, 72 within 0.75 nm; min OW distance 0.28081 nm
- PYR4: 4 OW within 0.35 nm, 19 within 0.50 nm, 71 within 0.75 nm; min OW distance 0.28945 nm
- PYR5: 6 OW within 0.35 nm, 23 within 0.50 nm, 82 within 0.75 nm; min OW distance 0.28030 nm

Interpretation:

The four pyrene chromophores have comparable local hydration shells. PYR5 is slightly more externally hydrated, consistent with its larger radial displacement from the tube axis, but no chromophore shows anomalous dehydration or excessive water contact.
