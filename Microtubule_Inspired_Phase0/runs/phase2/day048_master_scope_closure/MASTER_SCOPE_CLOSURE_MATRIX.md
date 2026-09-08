# MASTER SCOPE CLOSURE MATRIX

## Governing objective

Determine whether the engineered microtubule-inspired nanoscale structure materially prolongs a clearly defined coherence lifetime at elevated temperature relative to an appropriate control.

A negative result is scientifically acceptable.

Every remaining task must satisfy at least one of the following:

1. directly fulfill a promised master-document deliverable;
2. provide a required validation for a promised deliverable;
3. contribute directly to the final elevated-temperature coherence/control decision.

Tasks that satisfy none of these criteria should be deferred to future work.

---

# M1 — Physical Feasibility and Model Definition

## Master-document deliverables

- structural assessment
- mechanism feasibility matrix
- simulation assumptions
- identification of active sites
- identification of dipoles
- identification of chromophores
- identification of cavities
- identification of possible spin-active centers
- final simulation protocol

## Current status

STATUS: CLOSEABLE

Existing work already establishes:
- selected engineered HBN–PYR architecture
- confined-water geometry
- pyrene chromophore sites
- dipolar degrees of freedom
- mechanism-selection logic
- conditional spin rule
- canonical simulation workflows

## Remaining closure action

PRIORITY: REQUIRED / LOW COST

- verify that all M1 deliverables are represented in the final traceability package;
- freeze one canonical model/protocol description;
- explicitly classify spin-active centers as supported, unsupported, or unresolved.

No new simulation is required for M1 unless a missing physical definition is discovered.

---

# M2 — Molecular Dynamics / Thermal Stability / Confined Water

## Master-document deliverables

- prepared MD system
- equilibrated structure
- trajectories at selected temperatures
- RMSD/RMSF
- geometric stability plots
- confined-water density
- residence analysis
- diffusion
- hydrogen-bond analysis

Recommended temperature points in the master strategy:
150 K
200 K
250 K
300 K

350 K is optional.

## Validated current state

STATUS: ADVANCED / PARTIALLY OPEN

Validated:
- canonical hydrated MD system
- major structural/geometric analysis
- confined-water density analyses
- residence / survival analysis
- diffusion / MSD analysis
- accepted 3D diffusion estimate approximately 2.254e-9 m^2/s
- multiple historical temperature-ramp products exist

## Remaining mandatory closure

PRIORITY: REQUIRED

M2-A:
Hydrogen-bond analysis.

M2-B:
Final traceability packaging of:
- density
- residence
- diffusion
- geometry
- temperature coverage

M2-C:
Verify explicitly whether the final accepted architecture has defensible results at:
150 K
200 K
250 K
300 K

Historical exploratory trajectories do not automatically satisfy this requirement unless they correspond to the accepted/canonical architecture or are clearly labeled as feasibility-stage evidence.

## Decision rule

Do not rerun accepted residence/diffusion calculations unless a downstream inconsistency is identified.

---

# M3 — Dipolar / Spectral / Dielectric Response

## Master-document deliverables

- dipole orientation distributions
- dipole autocorrelation functions
- relaxation curves
- spectral densities
- field-response tests if applicable
- preliminary effective dielectric estimates

## M3.1 Orientation

STATUS: PASS / CLOSED

## M3.2 Autocorrelation

STATUS: PASS / CLOSED

## M3.3 Relaxation

STATUS: PASS / CLOSED

Validated 300 K engineered values:

C1 integrated KWW relaxation:
5.624045366 ps

C2 integrated KWW relaxation:
2.007253750 ps

The engineered-vs-bulk apparent slowdown is transient.

Persistent stationary orientational slowdown:
NOT RESOLVED.

## M3.4 Spectral density

STATUS: OPEN WITH METHOD BOUNDARY

The canonical 0.5 ps frame spacing implies:

Nyquist frequency = 1 THz.

Permitted:
- low-frequency dipolar/environmental spectrum
- slow collective response

Not permitted:
- microscopic excitonic spectral density across unresolved ultrafast frequencies
- direct derivation of microscopic electronic gamma_phi

Required closure decision:
either
A. produce a properly normalized low-frequency dipolar spectral-density deliverable with explicit 1 THz bandwidth boundary,
or
B. classify microscopic spectral density as not extractable from the available trajectory and document the limitation.

## M3.5 Collective dielectric response

STATUS: ACTIVE

100 ps collective polarization:
NOT CONVERGED.

F46Q:
extend engineered + bulk to 200 ps.

F46R:
test:
- cumulative tensor convergence
- bulk isotropy
- 50 ps blocks
- 100 ps blocks
- collective autocorrelation
- engineered/bulk ratios

Primary convergence gates:

bulk diagonal max/min <= 1.5

bulk diagonal CV <= 0.20

If PASS:
proceed to low-frequency dielectric-response construction.

If FAIL:
do not extend MD incrementally without redesigned sampling.

---

# M4 — Electronic and Excitonic Parameterization

## Master-document deliverables

- active-site selection
- DFT/TDDFT
- site energies
- transition dipoles
- excitonic couplings
- Hamiltonian matrices
- optional optical spectra

STATUS: CLOSED / PASS

Canonical model:
4_STATE_TRACKED_BRIGHT

Sensitivity:
8_STATE_S1_S2

Physical regime:
WEAK_COUPLING_STRONG_ENERGETIC_DETUNING

Do not reopen unless a downstream inconsistency requires it.

---

# M5 — Exciton Dynamics and Thermal Sensitivity

## Master-document deliverables

- population dynamics
- coherent-like beating analysis
- transfer pathways
- dephasing sensitivity
- relaxation sensitivity
- temperature/disorder sensitivity
- excitonic feasibility classification

STATUS: CLOSED / PASS

Validated interpretation:
- local/transient PYR2–PYR4 coherent mixing
- negligible coherent access to PYR5
- bath-assisted access to PYR5
- PYR4 dominant gateway
- weak-coupling/strong-detuning regime

Maximum ensemble-mean coherent PYR5 population:
approximately 8.2e-05

Absolute bath rates remain phenomenological.

No microscopic quantum coherence lifetime is inferred from Haken–Strobl gamma_phi or phenomenological relaxation kappa.

---

# M6 — THz / Microwave / Alternative Coherence

## Master-document deliverables

- effective dielectric parameters
- COMSOL electromagnetic model
- resonance analysis
- phase-shift analysis
- field-localization maps
- optional spin feasibility analysis

STATUS: OPEN / MAJOR REMAINING WORK

## M6.1 Environmental parameter freeze

DEPENDENCY:
M3 closure.

Required:
- epsilon or susceptibility representation with explicit valid frequency window
- anisotropy if supported
- control/reference representation
- uncertainty/sensitivity bounds

## M6.2 Continuum model

REQUIRED

Build/finalize frequency-domain electromagnetic model.

Must include:
- engineered geometry
- control/reference
- physically consistent environmental parameters
- defined excitation/boundary conditions

## M6.3 Electromagnetic outputs

REQUIRED

Produce:
- resonance frequencies
- phase response
- field localization
- losses/absorption where justified
- engineered versus control comparison

## M6.4 Final coherence-mechanism gate

REQUIRED

For every surviving candidate mechanism define:

- coherence observable
- lifetime definition
- extraction method
- elevated-temperature condition
- control
- uncertainty
- enhancement magnitude
- mechanism assignment

If no meaningful enhancement exists:
state negative result explicitly.

---

# Spin Branch

STATUS: CONDITIONAL

The existing high-spin h-BN calculation is not evidence of:
- spin transport
- spin coherence
- spin qubit behavior
- long T2

Spin work is activated only if an actual spin-active center/manifold can be mapped to:

- a physical carrier
- Hamiltonian
- environmental coupling
- temperature dependence
- measurable coherence observable
- appropriate control

If this cannot be justified efficiently, classify spin mechanism as unsupported by the current model.

---

# Publication / Defense Requirements

PRIORITY: REQUIRED BEFORE FINAL REPORT

For every central result preserve:

- canonical input
- script
- output
- units
- normalization
- control definition
- uncertainty or convergence statement
- physical interpretation
- limitation
- provenance

Avoid:
- converting phenomenological parameters into measured microscopic lifetimes;
- interpreting nonstationary differences as persistent effects;
- using unconverged dielectric tensors;
- claiming spin coherence without a spin model;
- presenting exploratory Phase-0 controls as equivalent to final controls.

---

# Critical Path

CURRENT:

F46Q
→ F46R
→ M3.5 decision
→ M3.4 low-frequency spectral closure
→ M3 final closure

PARALLEL LOW-COST:

M2 hydrogen-bond closure
→ M2 traceability packaging

AFTER M3:

M6.1 environmental parameter freeze
→ M6.2 COMSOL model
→ M6.3 resonance / phase / field localization
→ M6.4 elevated-temperature coherence/control decision

FINAL:

Outcome A:
physically meaningful coherence enhancement demonstrated

or

Outcome B:
no meaningful enhancement demonstrated, with mechanism-specific explanation.

Both outcomes satisfy the scientific objective if quantitatively supported.

---

# Priority Classification

## P0 — Must complete

- F46Q completion
- F46R convergence decision
- M3.4 defensible spectral deliverable/boundary
- M3.5 dielectric-response decision
- M2 hydrogen bonds
- M2 final traceability
- verify accepted temperature coverage
- M6 environmental parameter freeze
- M6 COMSOL model
- resonance / phase / field maps
- final elevated-temperature coherence-versus-control comparison

## P1 — Publication-strengthening

- uncertainty / replicate strategy where required
- final canonical figures
- sensitivity bounds
- consolidated provenance tables
- master traceability matrix

## P2 — Future work unless needed

- additional MD replicas without a failed convergence requirement
- optional 350 K extension
- microscopic excitonic spectral density requiring new high-time-resolution simulations
- spin modeling without independently justified spin centers
- additional architecture variants
- exploratory parameter sweeps not linked to the final coherence question

