# M3 Environmental Dynamics and Dielectric Response — Status Checkpoint

## Scope

M3 evaluates the environmental molecular and collective dynamics of the hydrated HBN–PYR engineered system and establishes which environmental observables can be used as physically defensible inputs or constraints for downstream coherence/dynamical models.

---

## M3.1 Molecular dipole orientation

**STATUS: PASS / CLOSED**

Water molecular orientation is defined using the PBC-corrected molecular dipole direction

u_mu = [r_(H1+H2)/2 - r_O] / |r_(H1+H2)/2 - r_O|.

The canonical 300 K engineered trajectory contains 16634 TIP4P/2005 waters, 201 frames, 100 ps total duration, and 0.5 ps frame spacing.

---

## M3.2 Molecular orientational correlation

**STATUS: PASS / CLOSED**

The first- and second-rank molecular orientational correlation functions C1(t) and C2(t) were calculated using all available time origins with PBC-corrected molecular reconstruction.

For the engineered HBN–PYR system at 300 K:

C1 KWW tau = 4.996633661 ps
C1 beta = 0.807471596
C1 integrated KWW relaxation time = 5.624045366 ps

C2 KWW tau = 1.444430311 ps
C2 beta = 0.640293946
C2 integrated KWW relaxation time = 2.007253750 ps

KWW behavior is strongly preferred over monoexponential decay.

These quantities describe environmental molecular orientational relaxation and are NOT quantum T2, T2*, or microscopic pure-dephasing times.

---

## M3.3 Engineered versus bulk molecular relaxation

**STATUS: PASS / CLOSED**

A physically equilibrated bulk TIP4P/2005 reference was constructed after removal of the HBN–PYR architecture and NPT equilibration to approximately 0.992 g/cm3.

Full 0–100 ps comparison initially produced:

C1 engineered/bulk ratio = 1.096100
C1 apparent change = +9.61 %

C2 engineered/bulk ratio = 1.123106
C2 apparent change = +12.31 %

However, stationarity analysis showed monotonic disappearance of this difference when progressively excluding the initial part of the trajectory.

C1 integrated relaxation enhancement:
0–100 ps: +9.61 %
25–100 ps: +4.92 %
50–100 ps: +3.16 %
75–100 ps: +0.04 %

C2 integrated relaxation enhancement:
0–100 ps: +12.31 %
25–100 ps: +5.64 %
50–100 ps: +4.52 %
75–100 ps: +1.90 %

### Frozen interpretation

The initially apparent orientational slowdown is predominantly a transient/intermediate-time effect.

A persistent stationary enhancement of molecular water orientational relaxation relative to bulk TIP4P/2005 is NOT resolved at 300 K.

---

## M3.4 Spectral-density interpretation

**STATUS: METHOD BOUNDARY DEFINED / MICROSCOPIC EXCITONIC SPECTRAL DENSITY NOT CLAIMED**

The available canonical trajectory has 0.5 ps temporal sampling, corresponding to a Nyquist frequency of 1 THz.

This sampling supports molecular orientational relaxation and low-frequency environmental response.

It does NOT resolve the ultrafast environmental fluctuations required to construct a microscopic excitonic spectral density or derive a microscopic electronic dephasing rate.

Molecular water C1(t), C2(t), or their Fourier transforms must not be interpreted as an excitonic bath spectral density.

The existing phenomenological bath models used in M5 remain phenomenological unless a separate site-energy-gap fluctuation calculation is performed.

---

## M3.5 Collective charge-weighted polarization

**STATUS: ACTIVE**

The physically relevant total water polarization observable has been implemented as a charge-weighted TIP4P/2005 molecular dipole sum:

M(t) = sum_m mu_m(t)

with molecular reconstruction under periodic boundary conditions.

The previous mean-unit-vector collective proxy is not used as a dielectric polarization observable.

At 100 ps, the bulk reference showed poor collective convergence:

bulk fluctuation tensor:
XX = 111.538485
YY = 34.232970
ZZ = 35.071467

bulk diagonal CV = 0.736423
bulk max/min diagonal ratio = 3.258218

Therefore, 100 ps is insufficient to support a converged collective dielectric-response interpretation.

The apparent engineered/bulk trace ratio of 0.800947 is NOT presently interpreted as a physical dielectric reduction.

---

## Collective-sampling extension

**STATUS: RUNNING**

Both engineered and bulk trajectories are being extended from 100 ps to 200 ps cumulative sampling.

The engineered continuation is exactly continuous with the accepted 100 ps trajectory.

No discontinuity was detected in:
- coordinates/checkpoint state
- energy
- temperature
- simulation timestep
- PME settings
- thermostat
- freeze groups
- nonbonded cutoffs

The engineered baseline contains a known very large static LJ/virial contribution and nominal pressure. This feature was already present in the accepted frozen-solute baseline and was not introduced by trajectory continuation.

No pathological short contacts were detected for:
- water–water
- water–HBN
- water–PYR
- HBN–PYR

Because the HBN/PYR solute is frozen, this static energetic feature is retained as a topology/solute-energy caveat and is not used as a thermodynamic pressure observable for the current water-dynamics analysis.

---

## F46R decision gate

After both systems reach 200 ps, the combined 0–200 ps trajectories will be evaluated using:

1. cumulative charge-weighted polarization tensors at 25 ps increments;
2. independent 50 ps blocks;
3. independent 100 ps blocks;
4. bulk isotropy diagnostics;
5. collective connected polarization autocorrelation;
6. first-zero correlation integrals;
7. engineered/bulk fluctuation ratios.

Primary bulk convergence gates:

diag_max_min_ratio <= 1.5

diag_CV <= 0.20

If these gates pass, M3.5 may proceed to the next low-frequency dielectric-response validation.

If they fail, the available trajectory will be classified as insufficient for a converged collective dielectric response, rather than extending MD indefinitely without a new sampling design.

---

## Current scientific boundary

Validated molecular water relaxation does not demonstrate prolonged quantum coherence.

At 300 K, no persistent enhancement of molecular orientational relaxation relative to bulk has been resolved.

The remaining M3 question is whether the engineered environment produces a robust and converged modification of low-frequency collective polarization fluctuations.

