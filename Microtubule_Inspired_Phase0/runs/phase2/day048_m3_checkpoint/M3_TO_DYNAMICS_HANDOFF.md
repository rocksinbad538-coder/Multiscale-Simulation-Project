# M3 → Dynamics / Coherence Handoff

## Purpose

This document defines which M3 environmental-dynamics results can be used downstream in excitonic, dielectric, or alternative coherence models, and explicitly identifies quantities that must not be interpreted as microscopic coherence parameters.

---

## 1. Validated environmental observables

### Molecular orientational dynamics

At 300 K, the engineered hydrated HBN–PYR system gives:

C1 integrated KWW relaxation time:
5.624045366 ps

C2 integrated KWW relaxation time:
2.007253750 ps

These quantities are converged molecular orientational relaxation observables.

They characterize water rotational/environmental relaxation.

They are not electronic coherence lifetimes.

---

## 2. Engineered-versus-bulk molecular result

The apparent full-trajectory slowdown relative to bulk TIP4P/2005 is not stationary.

C1 engineered/bulk enhancement:

0–100 ps: +9.61 %
25–100 ps: +4.92 %
50–100 ps: +3.16 %
75–100 ps: +0.04 %

C2 engineered/bulk enhancement:

0–100 ps: +12.31 %
25–100 ps: +5.64 %
50–100 ps: +4.52 %
75–100 ps: +1.90 %

### Frozen interpretation

No persistent enhancement of molecular water orientational relaxation relative to bulk is resolved at 300 K.

The architecture therefore cannot presently be claimed to prolong an environmental orientational lifetime in the stationary regime.

---

## 3. What M3 can constrain downstream

The validated molecular results can constrain:

- the approximate picosecond timescale of slow environmental orientational relaxation;
- whether confinement produces persistent or transient changes in molecular water dynamics;
- the qualitative separation between molecular orientational and collective polarization timescales;
- the low-frequency environmental response regime accessible to the MD trajectory;
- whether a dielectric-response model requires additional collective sampling.

These quantities may be used as environmental timescale constraints or consistency checks.

---

## 4. What M3 cannot determine directly

The following quantities must NOT be derived directly from C1(t), C2(t), or the current collective water correlation:

- electronic T2;
- electronic T2*;
- microscopic pure-dephasing rate gamma_phi;
- excitonic bath coupling amplitude;
- microscopic excitonic spectral density J(omega);
- site-energy fluctuation variance;
- site-energy cross-correlation functions;
- spin coherence lifetime;
- spin transport lifetime.

---

## 5. Temporal-resolution boundary

The canonical MD trajectory is sampled every 0.5 ps.

Therefore:

Nyquist frequency = 1 THz.

The data can support low-frequency molecular and collective environmental response.

The trajectory does not resolve the high-frequency environmental fluctuations required for a microscopic excitonic spectral density spanning the relevant electronic transition/bath frequencies.

Any Fourier-domain analysis must respect this limit.

---

## 6. Required microscopic bridge for excitonic dephasing

A physically defensible M3 → excitonic-dephasing bridge would require an observable such as:

delta E_i(t) = E_i(t) - <E_i>

where E_i(t) is the environment-dependent excitation/site energy of chromophore i.

The relevant correlation would then be:

C_ij(t) = <delta E_i(0) delta E_j(t)>

From these fluctuations one could, under an explicitly defined open-system model, construct bath correlation functions or spectral-density information.

Water orientational correlations alone do not provide the fluctuation amplitude or chromophore-specific coupling required for this mapping.

---

## 7. Relationship to M5

M5 previously used phenomenological Haken–Strobl dephasing and detailed-balance relaxation parameters.

Those parameters remain phenomenological.

M3 does not retrospectively convert them into microscopic rates.

M3 does provide a physically measured environmental timescale against which phenomenological parameter choices may be checked for order-of-magnitude consistency, but not calibrated uniquely.

---

## 8. Relationship to collective dielectric response

The charge-weighted TIP4P/2005 total-water polarization observable has been implemented.

At 100 ps, the collective fluctuation tensor was not converged.

Bulk isotropy failed strongly:

diag CV = 0.736423

max/min diagonal ratio = 3.258218

A 100→200 ps sampling extension is currently being completed for both engineered and bulk systems.

The 200 ps F46R analysis will determine whether low-frequency collective response can be interpreted further.

---

## 9. Relationship to M6

If the 200 ps collective response converges, the resulting low-frequency anisotropic/isotropic environmental response can be used as an input or constraint for M6 dielectric/EM modeling.

If collective convergence remains inadequate, M6 must use either:

- clearly identified phenomenological dielectric parameters;
- independently validated literature/reference dielectric values;
- or a separate dedicated sampling campaign.

No dielectric constant will be inferred from unconverged polarization fluctuations.

---

## 10. Current coherence conclusion

At present there is no validated evidence that the engineered HBN–PYR environment produces a persistent enhancement of a quantum coherence lifetime at 300 K.

The existing M3 results instead establish environmental molecular relaxation timescales and constrain the low-frequency dielectric-response problem.

This distinction must remain explicit in all downstream interpretation.

---

## 11. Next decision gate

F46R will analyze the complete 0–200 ps engineered and bulk trajectories.

Primary decision:

If bulk collective polarization becomes sufficiently isotropic and stable, proceed with low-frequency dielectric-response validation.

If it remains unconverged, do not extend trajectories incrementally without a redesigned sampling strategy.

