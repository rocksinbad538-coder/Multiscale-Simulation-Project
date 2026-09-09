# Microtubule-Inspired Multiscale Simulation — Current Status

Last updated: 2026-09-09

---

## Governing scientific question

Determine whether the engineered microtubule-inspired nanoscale architecture materially prolongs a clearly defined coherence lifetime at elevated temperature relative to an appropriate control.

A negative result is scientifically acceptable.

---

## Current project state

### M1 — Physical Feasibility and Model Definition

**PASS / CLOSED**

The canonical engineered architecture, physical degrees of freedom, control hierarchy, candidate mechanisms, and multiscale interpretation boundaries are frozen.

---

### M2 — Molecular Dynamics, Thermal Stability, and Confined Water

**PASS / CLOSED**

Validated results include confined-water residence, confined-water diffusion, water-water hydrogen bonding, and the thermal campaign at 150, 200, 250, and 300 K.

Accepted confined-water diffusion estimate:

**D_3D ≈ 2.254 × 10^-9 m^2/s**

The engineered water-water hydrogen-bond network remains essentially bulk-like at 300 K.

No additional M2 production is required.

---

### M3 — Molecular Environmental Dynamics and Collective Polarization

**PASS / CLOSED**

M3 is scientifically and formally closed.

#### Bulk reference

The primary bulk TIP4P/2005 control is converged using three matched 300 ps replicas.

The equal-weight ensemble passes the frozen bulk isotropy gate.

#### Engineered reference

The final engineered dataset contains:

- 3 replicas;
- 500 ps per replica;
- 1001 frames per replica;
- 0.5 ps sampling;
- 16,634 waters;
- 68,320 atoms.

The final primary tensor-convergence gate passes for trace and all three ordered eigenvalues.

Minimum full-window effective sample size:

**Neff_min = 40.119**

Practical effective-sampling gate:

**PASS**

Collective-correlation convergence:

**PASS**

No additional M3 collective-polarization MD sampling is required.

#### Orientational dynamics

The apparent engineered orientational slowdown relative to bulk is predominantly transient/intermediate-time.

A persistent stationary orientational enhancement is not resolved at 300 K.

#### Engineered versus bulk collective response

The statistically controlled comparison is restricted to matched windows of at most 300 ps because the converged three-replica bulk dataset contains 300 ps per replica.

Some full-window tensor/eigenvalue differences are nominally substantial.

However, after late-window and replica-reproducibility testing:

**No persistent replica-reproducible engineered-versus-bulk collective tensor enhancement or suppression is established.**

All primary tensor observables remain classified as:

**UNRESOLVED_OR_WEAK**

A reproducible architecture-locked principal polarization direction is also not established.

#### Interpretation boundary

M3 characterizes environmental molecular and collective-polarization behavior.

M3 does not establish:

- excitonic coherence lifetime;
- spin coherence lifetime;
- quantum coherence protection;
- a microscopic excitonic spectral density;
- a general dielectric constant.

The 0.5 ps trajectory sampling restricts the collective spectral analysis to frequencies at or below 1 THz.

---

### M4 — Electronic Structure and Excitonic Hamiltonian

**PASS / CLOSED**

Canonical electronic model:

**4_STATE_TRACKED_BRIGHT**

Sensitivity model:

**8_STATE_S1_S2**

Established physical regime:

**WEAK_COUPLING_STRONG_ENERGETIC_DETUNING**

The electronic calculations provide the Hamiltonian-level foundation for the downstream coherence analysis.

---

### M5 — Excitonic Dynamics

**PASS / CLOSED**

Current evidence supports:

- local/transient coherent mixing among PYR2–PYR4;
- no supported long-range coherent transfer to PYR5;
- PYR5 as a strong low-energy sink;
- bath-assisted population access to PYR5;
- PYR4 as the dominant gateway.

The historical Haken–Strobl dephasing parameter and detailed-balance relaxation parameter are phenomenological.

They must not be interpreted as microscopic coherence lifetimes.

Static ensemble averaging is not microscopic decoherence.

---

## Current active scientific branch

### Direct coherence-lifetime comparison

**ACTIVE / NEXT**

The next central task is to determine whether the engineered environment materially prolongs a consistently defined coherence lifetime relative to an appropriate control.

Before production calculations, the branch must freeze:

1. the coherence observable;
2. the coherence-lifetime definition;
3. the microscopic system-environment coupling observable;
4. the engineered/control pairing;
5. the temperature protocol;
6. the uncertainty/statistical treatment;
7. the criterion for material and reproducible enhancement.

The preferred microscopic bridge is expected to require explicitly justified electronic site-energy or energy-gap fluctuations rather than substituting collective water-orientation or polarization observables directly for electronic decoherence.

---

## Conditional spin branch

If the excitonic coherence-lifetime comparison does not show a meaningful and reproducible enhancement, the project will evaluate the spin manifold as a distinct candidate coherence mechanism.

Existing multigeometry spin-manifold and M=6 excited-state calculations constitute electronic foundations only.

They do not establish spin-coherence protection or spin-coherence lifetime enhancement.

---

## M6 — Device-Level Electromagnetic / Functional Integration

**OPEN / DOWNSTREAM**

M6 remains downstream of the direct coherence-lifetime adjudication.

Required device-level scope includes:

- effective environmental/dielectric parameters with explicit caveats;
- frequency-domain electromagnetic model;
- resonance behavior;
- phase response;
- field localization;
- directional electromagnetic response where supported.

If neither the excitonic nor spin branch supports the original coherence-enhancement objective, M6 will include an evidence-ranked alternative-function assessment based only on reproducible results.

---

## Immediate next action

Perform a read-only audit of the existing M4/M5 electronic and dynamical outputs to determine what microscopic time-dependent Hamiltonian/environment information already exists and what is still required for the direct coherence-lifetime calculation.

No new expensive production calculation is authorized until that audit defines and freezes the coherence observable, lifetime estimator, control, temperature treatment, and acceptance criteria.
