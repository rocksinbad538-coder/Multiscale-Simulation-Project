# M1 — Physical Feasibility and Model Definition

## Final status

**PASS / CLOSED**

M1 is considered complete for the current project scope.

No additional structural-model-definition work is required unless a downstream calculation reveals an internal inconsistency.

---

# 1. Governing project question

The project seeks to determine whether the engineered microtubule-inspired nanoscale architecture can materially prolong a clearly defined coherence lifetime at elevated temperature relative to an appropriate control.

A negative result is scientifically acceptable.

---

# 2. Canonical engineered architecture

The accepted computational architecture consists of:

- an h-BN nanoscale scaffold;
- confined water;
- four pyrene chromophore sites;
- an engineered PYR5 low-energy sink position in the final excitonic model;
- a hydrated finite simulation environment.

The architecture is treated differently depending on the downstream method:

- frozen HBN/PYR for the canonical 300 K water-dynamics trajectory used in M3;
- explicit chromophore electronic states in M4/M5;
- continuum effective-medium representation downstream in M6.

---

# 3. Physical degrees of freedom identified

## Structural

- h-BN scaffold geometry
- confined-water region
- chromophore positions
- chromophore orientations
- inter-site distances
- cavities / confined solvent domains

## Molecular

- water translational motion
- water orientational motion
- hydrogen-bond network
- residence / confinement
- diffusion

## Electronic

- local pyrene excited states
- transition dipoles
- site energies
- inter-site excitonic couplings
- energetic disorder

## Collective / electromagnetic

- collective water polarization
- effective low-frequency dielectric response
- resonant electromagnetic response

---

# 4. Active chromophore sites

The electronic/excitonic model contains four pyrene sites.

The canonical electronic representation is:

**4_STATE_TRACKED_BRIGHT**

with sensitivity model:

**8_STATE_S1_S2**

The physical regime established downstream in M4 is:

**WEAK_COUPLING_STRONG_ENERGETIC_DETUNING**

---

# 5. Dipolar degrees of freedom

For molecular water, the canonical orientation vector is defined from the oxygen atom toward the midpoint of the two hydrogen atoms.

This molecular dipole/orientation degree of freedom is used in M3 for:

- orientational distributions;
- C1(t);
- C2(t);
- relaxation analysis.

A separate charge-weighted molecular polarization observable is used for collective response.

These two quantities must not be conflated.

---

# 6. Confined-water definition

Confined water is treated as a distinct environmental subsystem whose:

- density;
- residence;
- diffusion;
- hydrogen bonding;
- orientational relaxation;
- collective polarization

are evaluated explicitly.

These quantities characterize the environment and are not quantum coherence lifetimes.

---

# 7. Candidate coherence mechanisms considered

## Excitonic coherence

Physically defined through the pyrene electronic-state Hamiltonian and open-system dynamics.

Status downstream:
tested quantitatively in M4/M5.

## Dipolar / dielectric environmental dynamics

Physically defined through molecular orientational and collective polarization observables.

Status downstream:
tested quantitatively in M3.

## Electromagnetic / THz / microwave response

Physically defined through frequency-domain electromagnetic response.

Status downstream:
reserved for M6 after M3 parameter closure.

## Spin coherence

Conditional only.

The current high-spin h-BN calculation does NOT establish:

- a spin qubit;
- spin transport;
- long T2;
- spin protection.

A spin branch is permitted only if a physical spin-active center and corresponding Hamiltonian/environmental coupling can be independently justified.

---

# 8. Control hierarchy

The following control semantics are frozen:

## Bulk water control

Primary M3 reference.

- same TIP4P/2005 water model;
- same temperature;
- same integrator;
- same thermostat;
- same PME settings;
- same timestep;
- same output spacing;
- no HBN/PYR architecture;
- equilibrated to bulk density.

This is the primary architecture-removed reference for the current M3 molecular and collective-water comparison.

## Matched frozen reference

The historical matched frozen trajectory has:

- same HBN/PYR composition;
- same topology;
- same water content;
- same frozen-solute protocol.

It is a protocol robustness reference, not an architecture-removed control.

## Historical water-only contained model

This is a finite reflective-box feasibility-stage control and must not be relabeled as bulk water.

---

# 9. Canonical scientific boundaries

The following claims are prohibited unless supported by downstream evidence:

- water orientational relaxation = electronic T2;
- Haken-Strobl gamma_phi = measured microscopic dephasing;
- static trajectory averaging = microscopic decoherence;
- high-spin electronic state = spin coherence;
- unconverged polarization tensor = dielectric constant;
- transient slowdown = persistent coherence protection.

---

# 10. Model assumptions

The current project uses a deliberately multiscale representation.

Different physical layers are not assumed to be dynamically equivalent.

## Molecular dynamics

Classical fixed-charge force fields.

## Electronic structure

DFT / TDDFT-derived site properties.

## Exciton dynamics

Reduced open-system Hamiltonian.

## Collective dielectric response

Charge-weighted molecular polarization.

## Electromagnetics

Continuum frequency-domain representation downstream in M6.

Cross-scale parameter transfer must preserve the physical meaning and uncertainty of each observable.

---

# 11. Known model caveats

## Frozen-solute water trajectory

The canonical engineered water trajectory contains frozen HBN/PYR coordinates.

Therefore it probes water dynamics in a static engineered environment and does not include scaffold vibrational motion.

## Engineered pressure

A very large static LJ/virial contribution exists in the frozen engineered trajectory.

No pathological solvent-solvent or solvent-solute geometry was detected.

The engineered pressure is therefore not used as a thermodynamic observable.

The frozen-solute energetic contribution remains a documented topology/virial caveat.

## Fixed-charge water

TIP4P/2005 does not include explicit electronic polarization.

High-frequency electronic dielectric response must therefore not be extracted directly from the MD polarization signal.

---

# 12. M1 feasibility decision

The architecture supports quantitative testing of:

- confined-water dynamics;
- dipolar relaxation;
- collective polarization;
- pyrene excitonic coupling and transfer;
- continuum electromagnetic response.

The architecture does NOT automatically imply long-lived quantum coherence.

The project proceeds by testing each candidate mechanism independently.

---

# 13. M1 closure decision

The following M1 requirements are complete:

- physical model definition;
- canonical architecture selection;
- chromophore identification;
- active-site identification;
- dipolar degrees of freedom;
- confined-water definition;
- candidate coherence mechanisms;
- control hierarchy;
- multiscale handoff logic;
- scientific limitations;
- conditional spin rule.

Therefore:

**M1 = PASS / CLOSED**

---

# 14. Downstream status

M2:
PASS / CLOSED

M3:
ACTIVE

M4:
PASS / CLOSED

M5:
PASS / CLOSED

M6:
OPEN

Current critical path:

F46Q
→ F46R
→ M3 closure
→ M6
→ final coherence-versus-control conclusion
