# QM_F06 LOWER Stage-2 Boundary Decision — Day027

## Finding

The only Stage-2 blocking contact is:

- `BR4:LOWER:00:3 — HCAP:LOWER:07`
- Stage-2 distance: `2.1518932544 Å`
- Distance / vdW-sum ratio: `0.6897093764`
- Graph separation: `4`

`HCAP:LOWER:07` is an artificial hydrogen cap attached to
`A:LOWER:14:-4`. It replaces the real cut neighbor
`A:LOWER:13:-3`.

## Contact evolution

- Initial capped geometry: `2.8671187259 Å`
- Stage 1: `2.4527921522 Å`
- Stage 2: `2.1518932544 Å`

The artificial cap moved progressively toward the chemically relevant
bridge region when released during Stage 2.

## Scientific interpretation

The B–N–B–N bridge itself remained chemically connected:

- 21 intended bonds preserved;
- 0 bond-range failures;
- 0 artificial-cap bond failures;
- both inherited Stage-1 compressed contacts improved.

The remaining failure is therefore a boundary-condition artifact rather
than failure of the bridge chemistry.

Because the artificial cap lies close to the bridge, retaining it would
contaminate the local electron density and any electrostatic-potential or
charge analysis.

## Decision

`QM_F06_LOWER_STAGE2_BRIDGE_VALID_BOUNDARY_CAP07_INVALID`

## Required correction

Replace `HCAP:LOWER:07` with the real R2 atom
`A:LOWER:13:-3` and include sufficient local coordination to move the
artificial boundary away from the bridge.

## Authorization state

- Current Stage-3 energy single point: **NOT AUTHORIZED**
- Current Stage-3 ESP/charge calculation: **NOT AUTHORIZED**
- Bridge chemistry rejected: **NO**
- Local fragment-boundary redesign required: **YES**
- Existing Stage-1 and Stage-2 results retained as scientific evidence: **YES**
