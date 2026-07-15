# Matched Snapshot-Pair Review for a Limited QM Pilot

## Scope

- Representative matched pairs audited: **21**
- Snapshot time field: `time_ps`
- Strong local PYR-hydration contrast pairs:
  **21**
- Strong lumen-occupancy contrast pairs:
  **14**
- Low-contrast control pairs:
  **0**

The purpose of this review is screening only. Selection from the same
trajectory cannot support a publication-level causal claim and cannot
replace independent trajectory replicas.

## Screening definitions

A pair is considered locally informative when either:

- the maximum absolute PYR hydration-shell difference within the
  nominal 0.50 nm shell is at least 2 waters; or
- the maximum absolute PYR hydration-shell difference within the
  nominal 0.35 nm shell is at least 1 water.

A lumen contrast requires an absolute mobile-minus-frozen difference
of at least 2 waters.

A control-like pair requires:

- absolute lumen difference no greater than 1 water;
- maximum local 0.35 nm difference no greater than 1 water; and
- maximum local 0.50 nm difference no greater than 1 water.

## Selected candidate pairs

- 43.000 ps: LUMEN_REHYDRATION_CONTRAST; ΔN_lumen=5.000; max |ΔN_PYR,0.35|=5.000; max |ΔN_PYR,0.50|=8.000
- 53.500 ps: LOCAL_PYRENE_HYDRATION_CONTRAST | LUMEN_REHYDRATION_CONTRAST; ΔN_lumen=6.000; max |ΔN_PYR,0.35|=4.000; max |ΔN_PYR,0.50|=17.000
- 65.500 ps: LUMEN_REHYDRATION_CONTRAST; ΔN_lumen=4.000; max |ΔN_PYR,0.35|=5.000; max |ΔN_PYR,0.50|=8.000
- 79.000 ps: LOCAL_PYRENE_HYDRATION_CONTRAST | LUMEN_REHYDRATION_CONTRAST; ΔN_lumen=4.000; max |ΔN_PYR,0.35|=5.000; max |ΔN_PYR,0.50|=16.000
- 100.000 ps: LOCAL_PYRENE_HYDRATION_CONTRAST; ΔN_lumen=-1.000; max |ΔN_PYR,0.35|=5.000; max |ΔN_PYR,0.50|=14.000

## Pilot gates

- all_21_pairs_valid: **PASS**
- at_least_2_local_contrast_pairs: **PASS**
- at_least_2_control_like_pairs: **FAIL**
- selected_5_to_7_pairs: **PASS**
- at_least_2_selected_local_contrasts: **PASS**
- at_least_3_selected_informative_pairs: **PASS**
- at_least_3_temporal_quartiles: **PASS**
- selected_temporal_span_at_least_50ps: **PASS**
- at_least_1_selected_control: **FAIL**

## Decision

- Decision: **LIMITED_PAIRED_QM_PILOT_NOT_JUSTIFIED**
- Limited paired mobile-versus-frozen QM pilot authorized:
  **NO**
- Full electronic recalculation authorized: **NO**
- Longer mobile production authorized: **NO**
- Publication-level causal claim authorized: **NO**
- Authorized next step:
  `RETAIN_FROZEN_QM_BASELINE_AND_DOCUMENT_MOBILITY_ASSOCIATED_SOLVENT_KINETICS`

The matched-water result remains a finding of mobility-associated
transient rehydration. It is not evidence of persistent confined-water
stabilization.
