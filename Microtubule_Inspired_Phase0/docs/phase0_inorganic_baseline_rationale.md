# Phase 0 Inorganic Baseline Rationale — BNNT / BN-like Tubular Scaffold

## Purpose

This document defines the rationale for selecting a BNNT / BN-like tubular scaffold as the first predominantly inorganic Phase 0 baseline candidate.

The objective is not yet to construct the final atomistic production model. The objective is to define a defensible starting point for comparing inorganic confinement against organic/biomimetic and hybrid tubular candidates under common Phase 0 conditions.

## Why BNNT / BN-like Scaffold?

Boron nitride nanotubes are relevant to this project because they provide:

- mechanically robust tubular confinement;
- local B–N polarity;
- dielectric contrast relative to carbon-based nanotubes;
- compatibility with confined-water simulations;
- a physically meaningful inorganic scaffold for testing water ordering, dipole relaxation, and field-induced response.

BNNTs have been studied in the context of water confinement, water transport, and electrostatic interactions between confined water and the B/N surface. Prior molecular dynamics and DFT-informed studies have used partial charges on B and N atoms to represent electrostatic effects arising from the polar B–N bond.

## Important Geometry Issue

The approved Phase 0 reference geometry is:

| Parameter | Value |
|---|---:|
| Outer diameter | 24.0 nm |
| Lumen diameter | 14.0 nm |
| Wall thickness | 5.0 nm |
| Segment length | 20.0 nm |

A single-wall BNNT does not naturally provide a 5 nm wall thickness. A single BN nanotube wall is effectively a rolled monolayer of h-BN. Therefore, matching both the 24 nm outer diameter and the 14 nm lumen diameter requires either:

1. a simplified BN-like tubular shell with finite wall thickness;
2. a multi-wall BNNT / multilayer BN shell;
3. a coarse-grained or continuum-equivalent wall model;
4. a smaller atomistic representative BNNT used only to parameterize local water-wall interactions.

## Phase 0 Recommendation

For the first Phase 0 comparison, use a geometry-equivalent BN-like tubular scaffold that preserves:

- the common outer diameter;
- the common lumen diameter;
- the 5 nm wall thickness;
- confined water;
- common temperature points;
- common analysis outputs.

This model should be treated as a controlled inorganic baseline, not as the final atomistically exact BNNT embodiment.

## Later Refinement

If the BN-like inorganic baseline shows strong confined-water ordering, dipole autocorrelation, dielectric relaxation, or field-induced response, then a more explicit model can be constructed using:

- single-wall BNNTs for local reference simulations;
- multi-wall BNNT or BN shell models;
- BN force fields such as Tersoff/extended Tersoff, ReaxFF, or literature-calibrated nonbonded B/N-water interactions;
- DFT-informed partial charges for B and N when electrostatic water-wall interactions are central.

## Initial Simulation Role

The initial inorganic model will be used to test:

- structural persistence of the tubular scaffold;
- confined-water density and residence;
- dipole orientation distribution;
- dipole autocorrelation;
- dielectric relaxation proxy;
- response to static or low-frequency electric field;
- contrast against dry and bulk-water controls.

## Technical Interpretation

The first BN-like scaffold is a Phase 0 screening model. It is intended to rank physical plausibility and output contrast, not to make final material-property claims about a specific synthesized BNNT.

