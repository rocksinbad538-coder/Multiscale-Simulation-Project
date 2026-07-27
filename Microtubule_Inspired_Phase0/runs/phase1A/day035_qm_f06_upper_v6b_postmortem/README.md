# QM_F06 UPPER V6-B Postmortem

## Status

Rejected.

The ORCA optimization converged numerically, but the optimized structure violated the structural acceptance gates.

## Root cause

The artificial boundary hydrogen HCAPV2:UPPER:03 migrated away from its intended boron atom.

The migration produced:

- loss of the canonical B-H bond
- formation of a noncanonical H-S interaction
- secondary B-P reconnectivity
- two overcoordinated atoms

Therefore the optimized structure no longer represented the intended QM fragment.

## Structural evidence

Lost bond

A:UPPER:14:2 -- HCAPV2:UPPER:03

New geometric contacts

A:UPPER:14:2 -- P:1641

HCAPV2:UPPER:03 -- S:1710

Trajectory onset

Frame 4:
HCAPV2:UPPER:03 begins approaching S:1710

Frame 58:
A:UPPER:14:2 approaches P:1641

Frame 60:
Canonical B-H bond is lost

## Scientific interpretation

The failure was not caused by numerical instability.

Instead, the chosen QM boundary was chemically incomplete.

The artificial boundary hydrogen was allowed to satisfy a neighboring sulfur atom because the adjacent canonical BN environment was absent.

## Corrective action

V7A replaces the artificial boundary with a canonical local expansion by introducing

- A:UPPER:13:1
- A:UPPER:14:0

and replacing the old cap with chemically consistent boundary passivation.

## Outcome

V6-B permanently rejected.

RESP prohibited.

Force-field generation prohibited.

Molecular dynamics prohibited.

Superseded by QM_F06 UPPER V7A.
