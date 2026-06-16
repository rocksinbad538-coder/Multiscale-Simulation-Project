# Accepted pyrene GAFF2/AM1-BCC parameterization

## Status

Accepted as preliminary MD-ready pyrene parameterization for Phase 1A.

This parameterization is suitable for initial classical MD preparation, dry-system sanity checks, hydration tests, and short NVT stability audits.

It is not labeled as RESP/DFT-final and should not be used as the final electronic-structure parameter source.

## Method

- Force field: GAFF2
- Charge route: AM1-BCC
- Net molecular charge: neutral
- Initial antechamber charge residual: +0.002 e
- Final neutralized charge residual: approximately -2.0e-06 e

## Accepted audit

- atom count: 26
- atom types: ca;ha
- frcmod exists: true
- prmtop exists: true
- inpcrd exists: true
- tleap errors: false
- frcmod missing/attention terms: false
- LEaP exit: Errors = 0; Warnings = 0; Notes = 0

## Use constraints

Accepted label:

- pyrene GAFF2/AM1-BCC MD-ready preliminary parameterization

Rejected labels:

- RESP-final
- DFT-final
- TDDFT-ready
- final excitonic parameterization
