# Phase 0 Candidate Matrix — Version 1

## Purpose

Phase 0 compares three representative tubular architectures under controlled and identical conditions. The objective is to identify which class provides the strongest combination of:

- structural persistence;
- confined-water stability;
- dipolar ordering;
- dielectric relaxation response;
- field-induced measurable output;
- computational tractability;
- physical plausibility.

## Common Reference Geometry

All Phase 0 candidates should initially use the same reference geometry:

| Parameter | Value |
|---|---:|
| Outer diameter | 24.0 nm |
| Lumen diameter | 14.0 nm |
| Wall thickness | 5.0 nm |
| Segment length | 20.0 nm |
| Tube axis | z |
| Confined medium | Water |
| Temperatures | 150, 200, 250, 300 K |

## Candidate Classes

| Class | Initial candidate | Model level | Primary role |
|---|---|---|---|
| Organic / biomimetic | Peptide-like or aromatic peptide tubular scaffold with polar wall groups | Coarse-grained or simplified atomistic MD | Test biomimetic confinement, hydration-layer behavior, polar wall effects, and structural persistence |
| Hybrid organic–inorganic | BNNT or CNT scaffold with polar/peptide-like coating | Hybrid MD / simplified surface functionalization | Test whether a mechanically robust inorganic tube plus polar coating improves confined-water ordering and dielectric response |
| Predominantly inorganic | BNNT with confined water | Atomistic MD / LAMMPS | Test robust tubular confinement, confined-water ordering, dipole autocorrelation, and field-induced dielectric response |

## Initial Recommendation

The first inorganic baseline should prioritize BNNT over CNT because BNNT provides a robust tubular scaffold with local B–N polarity and a more interpretable dielectric/confined-water environment than a metallic or semiconducting CNT. CNT can remain as a later inorganic variant if conductive, plasmonic, or graphitic effects become a priority.

## Controls

| Control | Purpose |
|---|---|
| Dry / empty tube | Isolate the role of confined water |
| Water-filled unfunctionalized tube | Separate confinement from wall functionalization |
| Bulk / unconfined water | Compare confined versus unconfined water relaxation |
| Randomized wall polarity | Test whether ordered dipolar patterning matters |

## First Analysis Outputs

| Output | Purpose |
|---|---|
| Geometry stability | Determine whether the structure persists at 150–300 K |
| Water density profile | Determine whether water remains confined and spatially structured |
| Water residence | Estimate confinement persistence and exchange dynamics |
| Dipole orientation distribution | Quantify orientational ordering of confined water |
| Dipole autocorrelation | Estimate relaxation behavior |
| Dielectric relaxation proxy | Identify field-relevant relaxation timescales |
| Field-induced polarization | Estimate measurable response under electric perturbation |

## Decision Criterion

The selected Phase 0 winner should be the architecture that maximizes:

```text
physical plausibility + computational tractability + measurable output contrast relative to controls
```
