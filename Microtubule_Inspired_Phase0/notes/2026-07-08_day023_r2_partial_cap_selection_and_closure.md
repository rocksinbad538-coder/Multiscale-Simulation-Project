# Day023 — R2 Partial-Cap Architecture Selection and Closure

**Date:** 2026-07-08  
**Project:** Multiscale Simulation — Microtubule-Inspired Nanoscale Body  
**Phase:** Phase 1A — Confinement-Architecture Screening  
**Status:** Scientific and documentary closure complete

## 1. Objective of the day

The objective was to complete the frozen-solute confinement screening
for the R2 symmetric partial-cap architecture, determine whether its
lumen-water population approached a stationary hydrated state, compare
it against the closed R1 positive control, and decide which architecture
should advance to a chemically realizable end-rim design gate.

No long mobile molecular dynamics, multitemperature simulation, or
quantum recalculation was authorized during this work.

## 2. R2 static partial-cap model

The selected R2 geometry contains:

- 144 neutral steric cap beads per end;
- 288 cap beads in total;
- a symmetric central aperture;
- a potential-defined effective aperture diameter of 0.839406 nm;
- an open-area fraction of 0.142928;
- 16,565 TIP4P/2005 waters;
- 428 initially luminal waters;
- 68,332 total atoms.

The static topology and CAP–water model passed all validation gates.

The independent analytical pure-repulsive CAP–water energy and the
GROMACS energy agreed within the previously validated numerical
tolerance. CAP–solute Lennard-Jones and Coulomb interactions remained
zero by construction.

## 3. Water-only minimization

The R2 water-only minimization was performed with HBN, all four pyrenes,
and both cap assemblies frozen.

Results:

- steepest-descents convergence in 1 step;
- maximum force: 382.908780 kJ mol^-1 nm^-1;
- HBN displacement: exactly 0;
- PYR displacement: exactly 0;
- cap displacement: exactly 0;
- water-O RMS displacement: 0.000029404 nm;
- water-O maximum displacement: 0.000530000 nm;
- lumen occupancy: 428 → 428;
- minimum CAP–OW distance: 0.220189 nm;
- no instability signatures.

Decision:

`R2_WATER_ONLY_ENERGY_MINIMIZATION_VALIDATED`

## 4. Initial 20 ps frozen-solute NVT screen

The prepared and executed R2 NVT protocol used:

- 0.0005 ps timestep;
- 40,000 steps;
- 20 ps duration;
- 41 trajectory frames;
- V-rescale thermostat acting effectively on SOL;
- deterministic velocity seed 20260708;
- HBN, PYR, and CAPS frozen in all dimensions.

The trajectory completed successfully and remained thermally and
numerically stable.

Key results:

- post-5 ps temperature: 299.7600 ± 1.2273 K;
- HBN/PYR/CAPS maximum displacement: exactly 0;
- minimum CAP–OW distance: 0.169759 nm;
- initial/minimum/endpoint occupancy: 428/414/416;
- endpoint occupancy fraction: 0.971963;
- endpoint retained initial identities: 414/428;
- maximum noninitial lumen waters: 4;
- no instability signatures.

The original 20 ps stationarity gate was not passed because the
10–20 ps occupancy slope was -0.509091 waters/ps, marginally outside
the predefined ±0.500000 waters/ps threshold.

The threshold was not relaxed.

Decision:

`R2_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW`

## 5. Occupancy-transient audit

A separate reproducible audit confirmed that the only failed 20 ps gate
was the second-half occupancy-slope criterion.

The final windows showed progressive flattening:

- 12.5–20 ps slope: -0.150000 waters/ps;
- 15–20 ps slope: -0.163636 waters/ps;
- 15–20 ps net occupancy change: 0 waters.

This supported an exact checkpoint continuation as a stationarity test,
without repeating the accepted 0–20 ps trajectory.

Decision:

`R2_OCCUPANCY_TRANSIENT_CHECKPOINT_EXTENSION_JUSTIFIED`

## 6. Exact checkpoint continuation from 20 to 50 ps

The continuation was generated using:

- the exact 20 ps checkpoint;
- source checkpoint step: 40,000;
- source checkpoint time: 20.000000 ps;
- target step: 100,000;
- target time: 50.000000 ps;
- 60,000 additional steps;
- no velocity regeneration;
- no thermostat-state regeneration;
- no rerun of the 0–20 ps segment.

The extended TPR was audited against the source TPR.

The apparent residual difference was determined to be only the
nonphysical filename header printed by `gmx dump`. After removing this
header and normalizing the expected `nsteps` change:

- coordinate mismatches: 0;
- velocity mismatches: 0;
- box-state mismatches: 0;
- scalar-parameter mismatches: 0;
- physical TPR differences beyond `nsteps`: 0.

Decision:

`R2_CHECKPOINT_CONTINUATION_TO_50PS_AUTHORIZED`

The continuation was executed with 1 MPI thread and 12 OpenMP threads.
This thread count is recorded as execution metadata and does not change
the TPR, checkpoint state, or physical model.

## 7. Combined R2 0–50 ps validation

The exact checkpoint continuation completed successfully.

Execution integrity:

- final checkpoint step/time: 100,000 / 50.000000 ps;
- source frames: 41;
- continuation frames: 61;
- combined unique frames: 101;
- combined interval: 0.500000 ps;
- no instability signatures;
- HBN/PYR/CAPS displacement: exactly 0.

Thermal behavior:

- continuation temperature:
  299.8727 ± 1.3326 K;
- continuation temperature slope:
  -0.013063 K/ps;
- final-15 ps temperature:
  299.8449 ± 1.3533 K;
- final-15 ps temperature slope:
  -0.044881 K/ps.

Water-confinement behavior:

- initial occupancy: 428;
- mean occupancy: 416.6337;
- minimum occupancy: 409;
- endpoint occupancy: 411;
- endpoint occupancy fraction: 0.960280;
- endpoint retained initial identities: 409/428;
- endpoint initial-identity fraction: 0.955607;
- endpoint noninitial lumen waters: 2;
- maximum noninitial lumen waters: 4;
- minimum CAP–OW distance: 0.166949 nm.

Stationarity:

- final-20 ps mean/slope:
  412.4390 / -0.173868 waters/ps;
- final-15 ps mean/slope:
  411.6452 / +0.025806 waters/ps;
- final-10 ps mean/change/slope:
  411.7143 / 0 / +0.064935 waters/ps.

The final 10 ps interval showed no net water loss.

Decision:

`R2_FROZEN_SOLUTE_NVT_50PS_VALIDATED`

## 8. R1–R2 architecture comparison

### R1 closed positive control

- initial occupancy: 428;
- mean occupancy: 427.9109;
- minimum occupancy: 426;
- endpoint occupancy: 428;
- endpoint fraction: 1.000000;
- stationarity slope: +0.005226 waters/ps.

R1 remains the validated neutral frozen closed steric positive control.

It is not interpreted as the final chemical architecture.

### R2 partial-cap architecture

- initial occupancy: 428;
- mean occupancy: 416.6337;
- minimum occupancy: 409;
- endpoint occupancy: 411;
- endpoint fraction: 0.960280;
- endpoint relative to R1: 0.960280;
- retention penalty versus R1:
  3.9720 percentage points;
- effective aperture diameter: 0.839406 nm;
- open-area fraction: 0.142928;
- final-10 ps net occupancy change: 0;
- measurable water exchange through the aperture.

R2 sacrifices approximately 4% of the closed-control endpoint occupancy
while maintaining a highly hydrated and operationally stationary lumen
with controlled exchange.

## 9. Architecture decision

Final decision:

`R2_SELECTED_AS_PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE`

Architecture status:

- R1:
  `RETAIN_AS_CLOSED_STERIC_POSITIVE_CONTROL`
- R2:
  `PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE`
- R3:
  `DEFERRED_NOT_REQUIRED_AT_THIS_GATE`
- R4:
  `DEFERRED_NOT_REQUIRED_AT_THIS_GATE`

The R2 result validates the geometric and steric screening concept. It
does not establish chemical realizability because the current cap is
represented by ideal neutral frozen steric beads.

## 10. Scope restrictions at closure

Authorized:

- static design of a chemically realizable R2 terminal/end-rim model.

Not authorized:

- new R2 molecular dynamics;
- short mobile MD;
- long mobile MD;
- multitemperature MD;
- QM recalculation;
- excitonic reparameterization.

The accepted electronic baseline remains unchanged:

**time-dependent solvent-induced site energies under frozen-solute conditions.**

## 11. Required next step

The next session must begin with:

`BEGIN_R2_CHEMICALLY_REALIZABLE_END_RIM_DESIGN_GATE`

The immediate objective is to replace the ideal neutral frozen steric
partial cap with a chemically defensible terminal-ring or end-rim
realization while preserving:

- the validated central aperture;
- the steric exclusion envelope;
- the minimum water-contact safety distance;
- the pyrene/scaffold geometry;
- the demonstrated high-hydration regime.

No new simulation should be executed before this static
chemical-architecture gate is completed.
