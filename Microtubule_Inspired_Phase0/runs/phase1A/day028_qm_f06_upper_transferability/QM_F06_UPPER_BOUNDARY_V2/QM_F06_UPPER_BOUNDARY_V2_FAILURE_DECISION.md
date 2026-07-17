# QM_F06 UPPER Boundary V2 Failure Decision — Day028

## Validated properties

The initial UPPER Boundary V2 preserves:

- 28 atoms;
- 27 intended covalent edges;
- complete graph connectivity;
- expected B/N/H degrees;
- valid bonded distances;
- valid restored-coordinate seam bonds;
- no hard bridge–cap contacts.

## Blocking boundary defect

The only unintended covalent-distance contact is:

`HCAP:UPPER:05 — HCAPV2:UPPER:02`

with distance:

`0.4498655562 Å`

Both artificial hydrogens terminate B–N cuts directed toward the same
omitted real R2 atom:

`A:UPPER:10:4`

Specifically:

- `A:UPPER:11:5 — A:UPPER:10:4`
- `A:UPPER:11:3 — A:UPPER:10:4`

## Chemical interpretation

Two hydrogen caps cannot independently replace the same omitted
three-coordinate nitrogen center. The cap–cap overlap is therefore a
boundary-partition artifact rather than a geometry that should be
repaired by arbitrary cap rotation.

## Selected correction

Audit and restore the shared real atom:

`A:UPPER:10:4`

Remove:

- `HCAP:UPPER:05`
- `HCAPV2:UPPER:02`

Then determine the minimum chemically complete treatment of the third
bond of `A:UPPER:10:4`.

## Authorization state

- Manual cap rotation: **REJECTED**
- Shared-atom expansion audit: **AUTHORIZED**
- Boundary V3 construction: **PENDING AUDIT**
- ORCA input preparation: **NOT AUTHORIZED**
- ORCA execution: **NOT AUTHORIZED**
