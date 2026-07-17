# QM_F06 LOWER Boundary V2-A Decision — Day028

## Validated findings

The V2-A optimization terminated normally and reached full geometry
convergence.

The optimized fragment preserves:

- all 27 intended bonds;
- complete graph connectivity;
- expected valence at every atom;
- all artificial-cap bond lengths;
- absence of unintended covalent connectivity;
- absence of hard bridge–cap contacts.

## Residual blocking contact

- Pair: `HCAP:LOWER:03 — A:LOWER:13:-3`
- Distance: `2.0956573628 Å`
- Distance / vdW-sum ratio: `0.6716850522`
- Involves bridge: `NO`
- Possible unintended covalent bond: `NO`

## Interpretation

`HCAP:LOWER:03` is an artificial V1 boundary hydrogen that remained fixed
during V2-A. `A:LOWER:13:-3` is a restored real R2 atom that was mobile.

The residual contact therefore represents incomplete boundary relaxation
against a fixed artificial cap, not failure of the restored BN topology
or of the B–N–B–N bridge.

## Decision

`QM_F06_LOWER_BOUNDARY_V2A_CORE_VALID_V2B_RELAXATION_REQUIRED`

## V2-B strategy

Release:

- all hydrogen atoms;
- all artificial boundary caps;
- the three restored real heavy atoms:
  - `A:LOWER:11:-3`
  - `A:LOWER:13:-3`
  - `A:LOWER:14:-2`

Keep fixed:

- the ten original heavy atoms defining the validated bridge and scaffold
  core.

## Authorization state

- Manual cap repositioning: **NOT SELECTED**
- Additional boundary expansion: **NOT REQUIRED AT THIS GATE**
- V2-B input preparation: **AUTHORIZED**
- V2-B execution: **PENDING STATIC AUDIT AND DRY PREFLIGHT**
- Final ESP/RESP calculation: **NOT AUTHORIZED**
