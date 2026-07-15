# Day022 Technical Closeout

## Validated baseline

The current architecture is retained as:

**R0 — Open-tube validated baseline**

The following components are accepted:

- mobile HBN/PYR topology;
- progressive restraint-release protocol;
- 100 ps fully mobile NVT trajectory at 300 K;
- 144 ps temporally matched frozen-solute control;
- exact initial-coordinate identity;
- bitwise-identical initial water velocities;
- structural and numerical stability of the mobile scaffold;
- frozen-solute electronic and excitonic baseline within its declared scope.

## Matched frozen-control validation

The matched frozen-solute control completed successfully:

- duration: 144 ps;
- frames: 289;
- coordinate interval: 0.5 ps;
- HBN/PYR RMS displacement: 0.000000000000 nm;
- HBN/PYR maximum displacement: 0.000000000000 nm;
- water-atom RMS displacement: 1.405231 nm;
- instability signatures: 0;
- matched comparison authorized: YES.

The control started from the same coordinates and bitwise-identical
initial water velocities as the mobile branch.

Initial water-velocity SHA256:

`d11ac84fa2d4a2a9a594a91ac0a6e0714dc3caa9da33f0ed12617371feff722d`

## Matched hydration result

The causal comparison used the same physical branch interval:

- mobile Stage08: local 0–100 ps, corresponding to branch 44–144 ps;
- frozen control: 44–144 ps;
- frames per trajectory: 201;
- sampling interval: 0.5 ps.

Results:

- mobile mean lumen occupancy: 1.905473 waters;
- matched-frozen mean lumen occupancy: 0.457711 waters;
- absolute paired difference: 1.447761 waters;
- block-bootstrap 95% interval: [0.515714, 2.650000] waters;
- exact block sign-flip p-value: 0.00585938;
- positive block fraction: 0.900000;
- mobile/frozen zero-occupancy fractions: 0.318408/0.716418;
- mobile/frozen positive-occupancy episode counts: 20/18;
- longest mobile/frozen positive episodes: 18.5/5.0 ps.

Scientific interpretation:

**Mobility-associated transient rehydration is supported.**

The result does not demonstrate persistent confined-water
stabilization. Both systems remain substantially depleted.

## Snapshot-pair review

- matched representative pairs audited: 21/21;
- detected local hydration metrics: 8;
- strong local hydration contrast pairs: 21;
- strong lumen contrast pairs: 14;
- low-contrast internal control pairs: 0;
- candidate pairs selected for screening: 5;
- selected temporal span: 57 ps;
- selected temporal quartiles represented: 3/4.

Failed pilot gates:

- at least two low-contrast control pairs;
- at least one selected low-contrast control pair.

Decision:

**LIMITED_PAIRED_QM_PILOT_NOT_JUSTIFIED**

The selected snapshots are contrast-enriched and would not support an
unbiased mobile-versus-frozen electronic comparison.

## Electronic baseline

The existing electronic results are retained as:

**time-dependent solvent-induced site energies under frozen-solute conditions**

They describe a nonstationary dehydration trajectory and must not be
presented as a stationary confined-water equilibrium ensemble.

The existing ORCA, site-energy, excitonic-Hamiltonian, coupling, and
open-system workflows remain valid within this declared scope.

## Authorizations

- Repeat current MD simulations: NO
- Longer mobile production: NO
- Multitemperature production: NO
- Limited paired mobile-versus-frozen QM pilot: NO
- Full electronic recalculation: NO
- Publication-level causal claim: NO
- Retain frozen QM baseline: YES
- Document mobility-associated solvent kinetics: YES
- Open architectural confinement-design branch: YES

## Next design branch

The validated R0 architecture will serve as the open-tube reference and
negative control for persistent water confinement.

Candidate refinements:

1. R1 — fully capped tube as a positive confinement control;
2. R2 — partially capped tube with an axial pore;
3. R3 — polar entrance rings;
4. R4 — sparse internal polar functionalization, only if required.

The next workday will begin by defining:

- confinement success criteria;
- topology and structure-generation workflows;
- R1–R3 structural gates;
- inexpensive frozen-solute screening protocols;
- criteria for authorizing mobile validation of only the best design.

No long MD or QM calculation will begin before these gates are
formally defined.
