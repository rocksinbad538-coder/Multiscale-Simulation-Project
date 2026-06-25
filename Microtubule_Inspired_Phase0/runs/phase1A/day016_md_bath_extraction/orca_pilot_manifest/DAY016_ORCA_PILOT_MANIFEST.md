# Day016 ORCA Pilot Manifest

## Purpose

This manifest defines a small pilot subset for testing the MD-to-electronic-extraction workflow before launching systematic calculations.

## Pilot frames

- frame000: beginning of accepted 100 ps trajectory
- frame010: middle of accepted 100 ps trajectory
- frame020: end of accepted 100 ps trajectory

## Pilot clusters

Monomers:

- PYR2
- PYR3
- PYR4
- PYR5

Pairs:

- PYR2_PYR3
- PYR3_PYR4
- PYR4_PYR5

## Important caveat

The current extracted clusters contain explicit nearby water molecules. These are useful for geometry/electrostatic screening but are too large for routine TDDFT across all frames and all clusters.

The production route should either:

1. Convert nearby waters into point-charge embedding, or
2. Use a smaller explicit water shell, or
3. Run only a limited number of explicit-water TDDFT pilots to calibrate a cheaper embedding approximation.

## Current scientific role

This pilot does not yet estimate sub-100 fs bath correlation times because the accepted trajectory is sampled every 0.5 ps. It is used to validate the MD-frame-to-electronic-input pipeline.
