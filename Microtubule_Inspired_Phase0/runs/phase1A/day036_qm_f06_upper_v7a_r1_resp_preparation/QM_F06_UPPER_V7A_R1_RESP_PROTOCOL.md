# QM_F06 UPPER V7-A R1 RESP Protocol — Day036

## Current formal state

Decision:

`QM_F06_UPPER_V7A_R1_RESP_PREPARATION_PASS_ESP_INPUT_PREFLIGHT_AUTHORIZED`

The UPPER V7-A R1 geometry has been optimized, independently audited,
coordinate-validated, and formally adopted.

This stage prepares, but does not execute, the ESP/RESP workflow.

## Adopted electronic-structure specification

- Method: PBE0-D4
- Basis: def2-TZVP
- Coulomb auxiliary basis: def2/J
- Exchange approximation: RIJCOSX
- SCF convergence: TightSCF
- Integration grid: DefGrid3
- Net charge: 0
- Multiplicity: 1
- Candidate ESP source: ORCA CHELPG

## Geometry and atom inventory

- Total atoms: 52
- Real transferable atoms: 37
- Artificial QM caps: 15
- Composition: B17N14H21
- Adopted geometry SHA256: `59cfd417753fbf6e5e4adf78a91761c2927824c50e24c7410010917d387574b2`

## Artificial-cap policy

Artificial caps remain included in the QM electrostatic calculation
because they are part of the finite electronic-structure model.

Their fitted charges must not be transferred directly into the full
R2 scaffold or included in transferable atom-type averages.

Final cap-charge removal or redistribution remains pending a joint
LOWER/UPPER protocol validation.

## Equivalence policy

No non-singleton RESP equivalence restraint is enforced at this stage.

Candidate groups are reported from element, node type, atom role, and
transfer status. They require independent topology and local-geometry
validation before any equality constraint may be applied.

Element identity alone is not accepted as evidence of charge
equivalence.

## Conformational scope

The current UPPER V7-A R1 geometry is an accepted reference geometry.

A transferable production charge model must later evaluate consistency
against the accepted LOWER reference and determine whether additional
conformers are required.

## Authorization state

- RESP preparation gate: PASS
- ESP input preflight: AUTHORIZED
- ESP execution: NOT AUTHORIZED
- RESP input generation: NOT AUTHORIZED
- RESP execution: NOT AUTHORIZED
- Charge adoption: NOT AUTHORIZED
- Force-field adoption: NOT AUTHORIZED
- Molecular dynamics: NOT AUTHORIZED

## Required next step

Run an independent preflight of the candidate ORCA ESP input, including:

1. executable-path validation;
2. geometry-hash validation;
3. charge and multiplicity validation;
4. CHELPG-grid validation;
5. method and basis validation;
6. output-file and working-directory isolation;
7. explicit execution authorization.
