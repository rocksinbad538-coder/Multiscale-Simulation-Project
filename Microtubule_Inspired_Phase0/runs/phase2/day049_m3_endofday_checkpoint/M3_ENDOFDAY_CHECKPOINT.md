# M3 End-of-Day Checkpoint

Date: 2026-09-07

## Milestone state

- M1: PASS / CLOSED
- M2: PASS / CLOSED
- M3: ACTIVE — late stage
- M4: PASS / CLOSED
- M5: PASS / CLOSED
- M6: OPEN — downstream of M3

## F46X bulk collective-polarization convergence

Three independently seeded bulk TIP4P/2005 trajectories were analyzed at a
matched 300 ps production duration per replica.

Ensemble tensor:

- XX = 60.78107808389157
- YY = 47.30760460176898
- ZZ = 60.62804883516797
- XY = -4.7201911316119025
- XZ = -4.694138561668152
- YZ = 4.896823870361548
- trace = 168.71673152082852

Convergence diagnostics:

- diagonal CV = 0.13754029333160112
- diagonal max/min ratio = 1.2848056585308227
- trace between-replica CV = 0.1457103290378888
- max |off-diagonal| / mean diagonal = 0.0870718124910514

Predefined bulk isotropy criteria:

- diagonal CV <= 0.20: PASS
- diagonal max/min <= 1.5: PASS

200 ps -> 300 ps ensemble changes:

- XX: 2.5162 %
- YY: 14.6560 %
- ZZ: 3.7407 %
- trace: 4.2598 %

Decision:

F46X_BULK_COLLECTIVE_CONVERGENCE = PASS
F46X_DECISION = BULK_REFERENCE_CONVERGED

Scientific interpretation:

The strong anisotropy observed in the original single bulk trajectory was
predominantly a finite-sampling directional fluctuation. Independent replicas
sample different dominant Cartesian directions, while their equal-weight
ensemble recovers the expected isotropic bulk response.

## F46Z engineered-topology provenance

Canonical frozen-solute topology:

parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032/
hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032.top

Canonical h-BN include:

parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032/
hbn_fixed_dummy.itp

Historical Day021 reconstruction established:

- accepted TPR atoms = 68320
- rebuilt TPR atoms = 68320
- grompp return code = 0
- inputrec exactly equal = True
- topology exactly equal = True
- box exactly equal = True
- frozen HBN/PYR coordinate prefix = PASS
- gmx check return code = 0
- exact original topology identity = True
- provenance validation = PASS

F46Z_TOPOLOGY_FILE_GATE = PASS

## Remaining M3 work

The bulk reference is no longer the blocking issue.

Remaining sequence:

1. Generate independent engineered frozen-solute replicas.
2. Test temporal and replica convergence of the engineered polarization tensor.
3. Do NOT impose isotropy on the engineered system.
4. Quantify tensor anisotropy, trace, principal values/directions and collective
   correlation stability.
5. Compare converged engineered response against converged bulk reference.
6. Extract the defensible low-frequency dielectric/environmental response within
   the limitations of fixed-charge TIP4P/2005 and the frozen scaffold.
7. Close M3 and hand off effective parameters to M6.

## Overnight plan

Generate two independently velocity-seeded engineered replicas from the
validated frozen-solute system, preserving the exact canonical topology,
freeze groups, PME settings, dt, thermostat, temperature and trajectory output
cadence.

Each new replica:

- 20 ps NVT reseed equilibration
- 300 ps NVT production

The engineered convergence criterion is replica/temporal stability, not
isotropy.
