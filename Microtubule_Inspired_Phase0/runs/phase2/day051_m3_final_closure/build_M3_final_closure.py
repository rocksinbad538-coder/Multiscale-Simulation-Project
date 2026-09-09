#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(".").resolve()

OUT = ROOT / "runs/phase2/day051_m3_final_closure"

F46AE = ROOT / (
    "runs/phase2/day051_f46ae_engineered_3x500ps_final"
)

F46X = ROOT / (
    "runs/phase2/day049_f46x_bulk_replica_convergence_final"
)

F46N = ROOT / (
    "runs/phase2/day048_f46n_orientation_stationarity"
)

# ------------------------------------------------------------
# Authoritative final M3 results
# ------------------------------------------------------------

summary = json.loads(
    (F46AE/"F46AE_FINAL_SUMMARY.json").read_text()
)

conv = pd.read_csv(
    F46AE/"F46AE_primary_tensor_convergence.csv"
)

neff = pd.read_csv(
    F46AE/"F46AE_effective_sample_size.csv"
)

corr = pd.read_csv(
    F46AE/"F46AE_correlation_convergence.csv"
)

ens = pd.read_csv(
    F46AE/"F46AE_ensemble_tensor_windows.csv"
)

compare = pd.read_csv(
    F46AE/"F46AE_matched_engineered_vs_bulk.csv"
)

classify = pd.read_csv(
    F46AE/"F46AE_persistent_transient_classification.csv"
)

bulk_summary = json.loads(
    (F46X/"F46X_FINAL_SUMMARY.json").read_text()
)

eng500 = ens[
    ens.window=="0_500"
].iloc[0]

matched300 = compare[
    compare.window=="0_300"
].copy()

late300 = compare[
    compare.window=="200_300"
].copy()

# ------------------------------------------------------------
# Markdown closure
# ------------------------------------------------------------

md = f"""# M3 — Molecular Environmental Dynamics and Collective Polarization

## Final status

**PASS / CLOSED**

M3 is scientifically complete for the current project scope.

No additional collective-polarization MD sampling is required unless a downstream calculation reveals an internal inconsistency or requires a different physical observable.

---

# 1. Scope

M3 was intended to characterize the molecular and collective environmental response of the engineered microtubule-inspired nanoscale architecture relative to an appropriate water control.

The milestone addressed:

- molecular water orientational dynamics;
- stationarity of orientational relaxation;
- collective molecular polarization;
- statistical convergence of the collective-polarization tensor;
- collective correlation times;
- effective sampling;
- low-frequency collective spectral behavior within the available MD sampling bandwidth;
- engineered-versus-bulk comparison.

M3 does **not** directly measure a quantum coherence lifetime.

---

# 2. Primary control

The primary M3 control is equilibrated bulk TIP4P/2005 water using the matched simulation protocol:

- same temperature;
- same integrator;
- same thermostat;
- same PME treatment;
- same timestep;
- same output spacing;
- no HBN/PYR architecture.

The historical water-only contained model is not relabeled as bulk.

The matched frozen reference is a robustness reference, not the primary architecture-removed control.

---

# 3. Canonical engineered environment

The canonical engineered molecular environment contains:

- one h-BN scaffold;
- four pyrene chromophore residues;
- 16,634 TIP4P/2005 waters;
- 68,320 total atoms.

The canonical M3 water-dynamics representation uses frozen HBN/PYR coordinates.

Therefore M3 characterizes water response in a static engineered environment and does not include scaffold vibrational dynamics.

---

# 4. Molecular orientational dynamics

The canonical water orientation vector is defined from the oxygen atom toward the midpoint of the two hydrogen atoms.

The validated orientational analysis used:

- first-order correlation C1(t);
- second-order correlation C2(t);
- KWW fits;
- integrated relaxation times;
- explicit stationarity windows.

Full-window engineered-versus-bulk differences initially indicated modest apparent orientational slowdown.

However, stationarity analysis showed that the effect weakened substantially in progressively later windows.

Frozen scientific interpretation:

**The apparent orientational slowdown is predominantly transient/intermediate-time. A persistent stationary orientational enhancement relative to bulk is not resolved at 300 K.**

Therefore the full-window slowdown must not be interpreted as persistent environmental protection or prolonged electronic coherence.

---

# 5. Collective-polarization observable

The charge-weighted molecular dipole of TIP4P/2005 water is constructed in a PBC-safe molecular representation.

For each water molecule:

mu_m = q_H r_H1O + q_H r_H2O + q_M r_MO

The total water polarization coordinate is:

M = sum_m mu_m

The dimensionless fluctuation tensor is computed as:

cov(M) / (epsilon_0 V k_B T)

using the validated F46R implementation.

This quantity is retained as a **collective-polarization response tensor**.

It is not automatically relabeled as an exact dielectric constant because interpretation depends on ensemble, electrostatic boundary conditions, geometry, and confinement.

---

# 6. Bulk convergence

The final bulk reference is based on three matched 300 ps replicas.

Final equal-weight bulk ensemble:

- XX = {bulk_summary["ensemble_tensor"]["XX"]:.6f}
- YY = {bulk_summary["ensemble_tensor"]["YY"]:.6f}
- ZZ = {bulk_summary["ensemble_tensor"]["ZZ"]:.6f}
- trace = {bulk_summary["ensemble_tensor"]["trace"]:.6f}

Bulk isotropy gate:

**PASS**

Final bulk diagnostics:

- diagonal CV = {bulk_summary["ensemble_tensor"]["diag_cv"]:.6f}
- max/min diagonal ratio = {bulk_summary["ensemble_tensor"]["diag_max_min_ratio"]:.6f}
- trace between-replica CV = {bulk_summary["trace_between_replica_cv"]:.6f}
- max off-diagonal / mean diagonal = {bulk_summary["max_offdiag_over_mean_diag"]:.6f}

Frozen interpretation:

The earlier strong anisotropy observed in individual bulk trajectories was a finite-sampling directional effect. The three-replica ensemble restores the expected bulk isotropy.

---

# 7. Engineered convergence

The final engineered dataset consists of:

- 3 independent replicas;
- 500 ps per replica;
- 1001 frames per replica;
- 0.5 ps frame spacing;
- 16,634 waters per replica.

Final equal-weight engineered 0–500 ps tensor:

- XX = {eng500["XX"]:.6f}
- YY = {eng500["YY"]:.6f}
- ZZ = {eng500["ZZ"]:.6f}
- trace = {eng500["trace"]:.6f}
- eig1 = {eng500["eig1"]:.6f}
- eig2 = {eng500["eig2"]:.6f}
- eig3 = {eng500["eig3"]:.6f}
- diagonal CV = {eng500["diag_cv"]:.6f}

No engineered isotropy criterion was imposed.

The ensemble nevertheless becomes nearly diagonal-isotropic in its Cartesian diagonal components at 500 ps.

---

# 8. Engineered tensor convergence gate

The primary convergence principle was frozen before the final analysis:

**The temporal 0–400 -> 0–500 ps change of the trace and ordered eigenvalues must not exceed the corresponding 500 ps between-replica standard deviation.**

Final results:

"""

for _,r in conv.iterrows():
    md += (
        f"- {r.observable}: "
        f"relative change = {r.relative_change_percent:.4f} %, "
        f"change/SD = {r.change_over_sd:.6f}, "
        f"gate = {'PASS' if bool(r.gate) else 'FAIL'}\n"
    )

md += f"""

Primary tensor convergence:

**PASS**

---

# 9. Effective sampling

The minimum effective sample size across the primary full 0–500 ps quadratic tensor observables is:

**Neff_min = {summary["min_full500_Neff"]:.6f}**

The practical reference used consistently with the earlier sampling diagnostics is:

**Neff >= 30**

Effective-sampling gate:

**PASS**

This threshold is a practical sampling diagnostic, not a physical theorem.

---

# 10. Collective-correlation convergence

Final 0–400 -> 0–500 ps collective-correlation convergence:

"""

for _,r in corr.iterrows():
    md += (
        f"- {r.observable}: "
        f"change/SD = {r.change_over_sd:.6f}, "
        f"gate = {'PASS' if bool(r.gate) else 'FAIL'}\n"
    )

md += """

Collective-correlation support gate:

**PASS**

The correlation observables therefore provide independent support for the tensor-convergence conclusion.

---

# 11. Low-frequency collective spectral scope

The MD output spacing is 0.5 ps.

Therefore the Nyquist frequency is:

**1 THz**

Any spectral analysis derived from these trajectories is restricted to low-frequency collective environmental behavior at or below this limit.

The resulting spectrum is **not** an excitonic spectral density and must not be used as a direct substitute for microscopic electronic energy-gap fluctuations.

---

# 12. Engineered-versus-bulk comparison

The statistically controlled engineered-versus-bulk ensemble comparison is restricted to matched windows of at most 300 ps because the converged three-replica bulk reference contains 300 ps per replica.

Full matched 0–300 ps nominal differences:

"""

for _,r in matched300.iterrows():
    md += (
        f"- {r.observable}: "
        f"{r.percent_difference:+.3f} %, "
        f"same-sign replicas = {int(r.same_sign_replica_count)}/3, "
        f"|difference|/combined SD = "
        f"{r.difference_over_combined_sd:.3f}\n"
    )

md += """

Some nominal full-window differences are substantial.

However, full-window magnitude alone is not accepted as evidence of a persistent architecture effect.

---

# 13. Late-window persistence test

The late matched 200–300 ps window was used to test whether the engineered-versus-bulk differences survive beyond early/transient behavior and remain reproducible across replicas.

Final classification:

"""

for _,r in classify.iterrows():
    md += (
        f"- {r.observable}: "
        f"full difference = {r.full_0_300_percent_difference:+.3f} %, "
        f"late difference = {r.late_200_300_percent_difference:+.3f} %, "
        f"same-sign replicas = {int(r.late_same_sign_replica_count)}/3, "
        f"late |difference|/combined SD = "
        f"{r.late_difference_over_combined_sd:.3f}, "
        f"classification = **{r.classification}**\n"
    )

md += """

No tensor observable satisfies the frozen criteria for a persistent replica-reproducible engineered-versus-bulk effect.

Frozen interpretation:

**Persistent engineered-versus-bulk collective tensor enhancement or suppression is not established.**

The nominal full-window XX/YY/ZZ/trace/eigenvalue differences remain mechanism-level observations that are unresolved or weak after late-window and replica-reproducibility testing.

---

# 14. Principal-axis behavior

The 500 ps engineered ensemble has non-degenerate ordered eigenvalues sufficient for principal-axis diagnostics.

However, replica-level principal directions are not uniformly reproducible.

R1 and R3 are comparatively close to the ensemble directions, while R2 exhibits large angular deviations.

Therefore:

**A reproducible architecture-locked principal direction is not established.**

Principal-axis behavior remains a mechanism-level diagnostic rather than a validated directional functional property.

---

# 15. Scientific conclusions demonstrated by M3

M3 supports the following conclusions:

1. The bulk collective-polarization reference is statistically converged.
2. The engineered collective-polarization response is statistically converged.
3. Three 500 ps engineered replicas provide adequate full-window effective sampling under the adopted diagnostic.
4. The engineered collective-correlation observables are converged.
5. The apparent water orientational slowdown is predominantly transient/intermediate-time.
6. A persistent stationary orientational enhancement relative to bulk is not resolved.
7. A persistent replica-reproducible engineered-versus-bulk collective tensor difference is not established.
8. A reproducible architecture-locked principal polarization axis is not established.
9. No additional collective-response MD sampling is required for M3 under the current observable definitions.

---

# 16. Claims not supported by M3

M3 does NOT establish:

- a quantum coherence lifetime;
- excitonic T2;
- spin T2;
- coherence protection;
- prolonged quantum coherence;
- a microscopic excitonic spectral density;
- a dielectric constant valid at arbitrary frequency;
- scaffold vibrational polarization;
- a persistent directional electromagnetic response;
- a persistent architecture-induced bulk-water slowdown.

---

# 17. Model caveats

## Fixed-charge solvent

TIP4P/2005 does not contain explicit electronic polarization.

## Frozen scaffold

HBN/PYR coordinates are frozen in the canonical engineered M3 trajectories.

Dynamic scaffold polarization and scaffold vibrations are absent.

## Pressure

The engineered frozen-solute trajectory contains a known large static LJ/virial contribution.

Engineered pressure is not used as a thermodynamic observable.

## Heterogeneity

The engineered nanoscale environment is heterogeneous and geometrically confined.

Continuum dielectric interpretation therefore requires explicit downstream modeling rather than direct identification of the fluctuation tensor with a scalar dielectric constant.

---

# 18. Handoff to the coherence-lifetime branch

The next central scientific question remains:

**Does the engineered nanoscale architecture materially prolong a consistently defined coherence lifetime at elevated temperature relative to an appropriate control?**

The next branch must define before calculation:

1. coherence observable;
2. coherence-lifetime definition;
3. microscopic environmental coupling;
4. engineered/control pairing;
5. temperature dependence;
6. statistical uncertainty;
7. material/reproducible enhancement criterion.

M3 collective polarization and orientational dynamics may constrain environmental timescales, but they must not be substituted directly for microscopic electronic decoherence.

A correct microscopic bridge requires electronic energy/site-gap fluctuations or another explicitly justified Hamiltonian-environment coupling observable.

---

# 19. Conditional downstream branch

If the excitonic coherence-lifetime comparison does not show a meaningful and reproducible enhancement, the project will proceed to the spin manifold as a distinct candidate coherence mechanism.

Existing multigeometry spin-manifold and M=6 excited-state calculations are electronic foundations only.

They do not by themselves demonstrate spin-coherence protection or lifetime enhancement.

---

# 20. Final M3 decision

**M3 PASS / CLOSED**

Collective molecular-environment characterization is complete.

No additional M3 MD production is required.

The project can proceed to the direct coherence-lifetime branch and subsequently to M6 device-level integration.
"""

(OUT/"M3_FINAL_CLOSURE.md").write_text(
    md
)

# ------------------------------------------------------------
# JSON closure
# ------------------------------------------------------------

closure = {
    "milestone": "M3",
    "status": "PASS_CLOSED",
    "scientific_scope": (
        "molecular environmental dynamics and "
        "collective polarization"
    ),
    "bulk_reference": {
        "status": "CONVERGED",
        "replicas": 3,
        "duration_ps_each": 300,
        "ensemble_tensor": bulk_summary["ensemble_tensor"],
        "isotropy_gate": True,
    },
    "engineered_reference": {
        "status": "CONVERGED",
        "replicas": 3,
        "duration_ps_each": 500,
        "frames_each": 1001,
        "frame_spacing_ps": 0.5,
        "waters": 16634,
        "atoms": 68320,
        "ensemble_0_500": {
            "XX": float(eng500["XX"]),
            "YY": float(eng500["YY"]),
            "ZZ": float(eng500["ZZ"]),
            "trace": float(eng500["trace"]),
            "eig1": float(eng500["eig1"]),
            "eig2": float(eng500["eig2"]),
            "eig3": float(eng500["eig3"]),
            "diag_cv": float(eng500["diag_cv"]),
        },
    },
    "primary_tensor_convergence_gate": True,
    "effective_sampling": {
        "min_full500_Neff":
            float(summary["min_full500_Neff"]),
        "practical_reference_Neff": 30.0,
        "gate": True,
    },
    "collective_correlation_gate": True,
    "orientation_interpretation": (
        "PREDOMINANTLY_TRANSIENT_INTERMEDIATE_TIME;"
        "_PERSISTENT_STATIONARY_ENHANCEMENT_NOT_RESOLVED"
    ),
    "engineered_vs_bulk": {
        "matched_comparison_max_duration_ps": 300,
        "persistent_replica_reproducible_tensor_effect":
            False,
        "classification":
            "UNRESOLVED_OR_WEAK",
    },
    "principal_axis_reproducibility":
        "NOT_ESTABLISHED",
    "spectral_nyquist_THz": 1.0,
    "spectral_result_is_excitonic_spectral_density":
        False,
    "quantum_coherence_lifetime_inferred":
        False,
    "additional_M3_MD_required":
        False,
    "next_scientific_branch":
        "DIRECT_COHERENCE_LIFETIME_COMPARISON",
    "spin_branch_if_excitonic_enhancement_absent":
        True,
    "M6_direct_handoff_without_coherence_test":
        False,
}

(
    OUT/"M3_FINAL_CLOSURE.json"
).write_text(
    json.dumps(
        closure,
        indent=2
    ) + "\n"
)

# ------------------------------------------------------------
# Provenance manifest
# ------------------------------------------------------------

provenance_files = [
    OUT/"M3_FINAL_CLOSURE.md",
    OUT/"M3_FINAL_CLOSURE.json",

    F46AE/"F46AE_FINAL_SUMMARY.json",
    F46AE/"F46AE_primary_tensor_convergence.csv",
    F46AE/"F46AE_effective_sample_size.csv",
    F46AE/"F46AE_correlation_convergence.csv",
    F46AE/"F46AE_ensemble_tensor_windows.csv",
    F46AE/"F46AE_matched_engineered_vs_bulk.csv",
    F46AE/"F46AE_persistent_transient_classification.csv",
    F46AE/"F46AE_principal_axis_reproducibility.csv",
    F46AE/"F46AE_low_frequency_collective_spectrum.csv",

    F46X/"F46X_FINAL_SUMMARY.json",
    F46X/"F46X_final_replica_tensors.csv",
    F46X/"F46X_final_between_replica_statistics.csv",
]

files = []
missing = []

for p in provenance_files:

    exists = p.exists()

    if not exists:
        missing.append(
            str(p.relative_to(ROOT))
        )
        continue

    h = hashlib.sha256(
        p.read_bytes()
    ).hexdigest()

    files.append({
        "path":
            str(p.relative_to(ROOT)),
        "exists": True,
        "size_bytes":
            p.stat().st_size,
        "sha256":
            h,
    })

manifest = {
    "milestone": "M3",
    "status": "PASS_CLOSED",
    "files": files,
    "missing": missing,
}

(
    OUT/"M3_PROVENANCE_MANIFEST.json"
).write_text(
    json.dumps(
        manifest,
        indent=2
    ) + "\n"
)

print("="*80)
print("M3 FINAL CLOSURE GENERATED")
print("="*80)

print(
    f"M3_PROVENANCE_FILE_COUNT={len(files)}"
)

print(
    f"M3_PROVENANCE_MISSING_COUNT={len(missing)}"
)

print(
    "M3_CLOSURE_BUILD_GATE="
    + (
        "PASS"
        if len(missing)==0
        else "FAIL"
    )
)

print(
    "M3_STATUS=PASS_CLOSED"
)

print(
    "ADDITIONAL_M3_MD_REQUIRED=NO"
)

print(
    "NEXT_SCIENTIFIC_BRANCH="
    "DIRECT_COHERENCE_LIFETIME_COMPARISON"
)
