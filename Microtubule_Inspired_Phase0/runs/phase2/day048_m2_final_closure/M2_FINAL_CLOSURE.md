# M2 — Molecular Dynamics, Thermal Stability, and Confined Water

## Final status

**PASS / CLOSED**

M2 is considered scientifically complete for the current project scope.

No additional molecular-dynamics production is required for M2 unless a downstream inconsistency is discovered.

---

# 1. Scope

M2 was intended to establish:

- molecular-dynamics system preparation;
- equilibration and thermal stability;
- structural stability;
- confined-water behavior;
- residence;
- diffusion;
- hydrogen bonding;
- thermal coverage over the project temperature range.

The project thermal campaign includes:

- 150 K
- 200 K
- 250 K
- 300 K

---

# 2. Canonical 300 K hydrated engineered system

Canonical trajectory:

runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc

Associated topology:

runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.tpr

Composition:

- HBN residues: 1
- PYR residues: 4
- SOL residues: 16634
- total atoms: 68320
- trajectory frames: 201
- duration: 100 ps
- frame spacing: 0.5 ps
- temperature: 300 K

The HBN/PYR architecture is frozen in this canonical water-dynamics trajectory.

---

# 3. Confined-water transport

## Residence

Confined-water residence/survival has been quantified using explicit occupancy definitions.

The continuous survival observable is conditioned on waters that remain within the confined region from the time origin through the requested lag time.

Residence analyses are therefore interpreted as confined-region persistence rather than unrestricted bulk residence.

## Diffusion

Accepted three-dimensional confined-water diffusion estimate:

**D_3D ≈ 2.254 × 10^-9 m^2/s**

This value characterizes retained confined water under the adopted analysis definition.

The diffusion result must not be interpreted as a universal bulk-water diffusion coefficient.

---

# 4. Water hydrogen bonding

## Status

**M2.6 PASS / CLOSED**

Canonical comparison:

engineered HBN–PYR hydrated system versus equilibrated bulk TIP4P/2005 at 300 K.

Geometric hydrogen-bond definition:

- O···O distance <= 3.5 Å
- O-H···O angle >= 150 degrees

Scope:

water-water hydrogen bonds only.

Water-HBN and water-PYR contacts are not classified as hydrogen bonds without separately justified donor/acceptor chemistry.

## Engineered

Mean hydrogen bonds:

30273.069652

Hydrogen bonds per water:

1.819951284

Mean coordination per water:

3.639902567

Fraction of waters participating:

0.999782260

Mean O···O distance:

2.867050 Å

Mean O-H···O angle:

167.104209 degrees

## Bulk

Mean hydrogen bonds:

30392.532338

Hydrogen bonds per water:

1.827133121

Mean coordination per water:

3.654266242

Fraction of waters participating:

0.999807384

Mean O···O distance:

2.866537 Å

Mean O-H···O angle:

166.948319 degrees

## Engineered versus bulk

Full-trajectory H-bond-count ratio:

0.996069341

Relative change:

**-0.3931 %**

After the initial 0–25 ps transient, the engineered H-bond network remains only slightly below bulk, with an approximate persistent difference below 1%.

### Frozen interpretation

The engineered environment does not produce a major strengthening or restructuring of the water-water hydrogen-bond network at 300 K.

Hydrogen-bond enhancement is therefore not supported as a mechanism for prolonged coherence or environmental protection.

---

# 5. Thermal campaign

Canonical corrected campaign:

runs/phase2/campaign_phase5_corrected/

Temperatures:

150 K
200 K
250 K
300 K

For every temperature:

- minimization pass = True
- heating pass = True
- NVT pass = True
- production pass = True
- production return code = 0
- quantitative trajectory analysis exists
- thermodynamic analysis exists
- aligned RMSD analysis exists

All four campaigns share the same minimized starting structure.

Canonical minimized.data SHA256:

970464bed50c53ab6ff50801355df88dab596e56be59934fc03ff5e42f7f85fc

This supports a common starting architecture across the thermal campaign.

---

# 6. Thermal-control adjudication

## 150 K

Full mean temperature:

150.019805 K

Last-quarter mean:

149.946515 K

Last-quarter error:

-0.053485 K

First-to-last-quarter shift:

0.101200 K

Temperature centering:

PASS

Temperature stationarity:

PASS

Production completion:

PASS

---

## 200 K

Full mean temperature:

199.970413 K

Last-quarter mean:

199.854855 K

Last-quarter error:

-0.145145 K

First-to-last-quarter shift:

0.048261 K

Temperature centering:

PASS

Temperature stationarity:

PASS

Production completion:

PASS

---

## 250 K

Full mean temperature:

250.044539 K

Last-quarter mean:

250.189315 K

Last-quarter error:

+0.189315 K

First-to-last-quarter shift:

0.047761 K

Temperature centering:

PASS

Temperature stationarity:

PASS

Production completion:

PASS

---

## 300 K

Full mean temperature:

300.004971 K

Last-quarter mean:

300.210630 K

Last-quarter error:

+0.210630 K

First-to-last-quarter shift:

0.201184 K

Temperature centering:

PASS

Temperature stationarity:

PASS

Production completion:

PASS

---

# 7. Thermal interpretation

All four production temperatures remain centered extremely close to their requested targets.

The instantaneous temperature distributions are broad, as expected for finite molecular systems, but there is no significant systematic drift in the time-averaged thermal state.

The four-temperature campaign is therefore accepted for project-level structural and thermal comparison.

Individual instantaneous extrema must not be interpreted as evidence of thermostat failure.

---

# 8. Structural observables

The corrected Phase-5 thermal campaign contains quantitative products for:

- aligned RMSD;
- RMSF;
- shape analysis;
- radius-of-gyration-related structural observables;
- thermodynamics;
- trajectory summaries.

These outputs remain available as the canonical thermal-stability evidence.

M2 closure does not imply that every instantaneous structural fluctuation is negligible.

It means that the required project-level structural and thermal observables have been computed, adjudicated, and preserved.

---

# 9. Scientific boundaries

M2 establishes molecular and environmental behavior.

M2 does NOT directly establish:

- electronic T2;
- electronic T2*;
- microscopic pure-dephasing rates;
- excitonic coherence lifetime;
- spin coherence lifetime;
- spin transport;
- dielectric convergence;
- electromagnetic coherence.

Water residence, diffusion, hydrogen bonding, and thermal stability are environmental observables.

They may constrain downstream physical interpretation but cannot themselves be relabeled as quantum coherence observables.

---

# 10. Key negative result

At 300 K, the water-water hydrogen-bond network remains essentially bulk-like.

Together with the M3 molecular-orientational analysis, this means that the present data do not support a strong stationary water-network protection mechanism at 300 K.

This is a valid scientific result and does not invalidate the project.

---

# 11. M2 closure decision

The following M2 requirements are considered complete:

- system preparation;
- canonical hydrated structure;
- thermal simulations;
- structural analysis;
- confined-water analysis;
- residence;
- diffusion;
- hydrogen bonding;
- temperature coverage from 150 to 300 K;
- provenance and interpretation boundaries.

Therefore:

**M2 = PASS / CLOSED**

No further M2 MD production should be initiated unless required by a specific failed downstream validation gate.

---

# 12. Handoff

The active project path is now:

F46Q bulk 200 ps completion

→ F46R collective-polarization convergence test

→ M3.5 decision

→ M3.4 low-frequency spectral closure

→ M3 final closure

→ M6 dielectric / electromagnetic modeling

→ final elevated-temperature coherence-versus-control decision

