# R1–R2 Architecture Comparison and Selection

## Scope

This comparison uses the validated 50 ps frozen-solute results for the
closed R1 positive control and the partially open R2 screening
architecture.

No new molecular dynamics, minimization, topology generation, or
quantum calculation was performed.

## Architectural roles

### R1

R1 remains the neutral frozen closed steric positive control. Its role
is to demonstrate that the confinement-detection method identifies
stable lumen hydration when the end boundary is closed.

R1 is not promoted as a chemically realizable architecture.

### R2

R2 is the symmetric partial-cap screening architecture. It introduces a
central aperture while maintaining a highly hydrated lumen and allowing
measurable exchange.

R2 is also not yet chemically realizable because its cap is represented
by ideal neutral frozen steric beads.

## Hydration comparison

- R1 initial/mean/minimum/endpoint:
  **428/427.9109/
  426/428**
- R2 initial/mean/minimum/endpoint:
  **428/416.6337/
  409/411**
- R1 endpoint fraction:
  **1.000000**
- R2 endpoint fraction:
  **0.960280**
- R2 endpoint relative to R1:
  **0.960280**
- R2 retention penalty versus R1:
  **3.9720 percentage points**
- R2 endpoint initial-identity retention:
  **0.955607**
- R2 endpoint/maximum noninitial lumen waters:
  **2.0/
  4.0**

## Stationarity

- R1 stationarity slope:
  **0.005226 waters/ps**
- R2 final 20 ps mean/slope:
  **412.4390/
  -0.173868 waters, waters/ps**
- R2 final 15 ps mean/slope:
  **411.6452/
  0.025806 waters, waters/ps**
- R2 final 10 ps mean/change/slope:
  **411.7143/
  +0/
  0.064935 waters, waters, waters/ps**

## Aperture and steric safety

- R1 effective aperture diameter:
  **0.000000 nm**
- R2 effective aperture diameter:
  **0.839406 nm**
- R2 open-area fraction:
  **0.142928**
- R2 minimum CAP–OW distance:
  **0.166949 nm**

## Execution metadata

- R2 continuation OpenMP threads:
  **12**

Thread count is recorded as computational metadata. It does not modify
the TPR physical parameters or checkpoint state.

## Selection gates

- `R1_closed_positive_control_is_validated`: **PASS**
- `R2_frozen_solute_50ps_is_validated`: **PASS**
- `R2_static_partial_cap_model_is_validated`: **PASS**
- `R1_and_R2_have_matched_initial_occupancy`: **PASS**
- `R1_behaves_as_closed_positive_control`: **PASS**
- `R2_endpoint_occupancy_is_at_least_90_percent`: **PASS**
- `R2_minimum_occupancy_is_at_least_80_percent`: **PASS**
- `R2_endpoint_is_at_least_90_percent_of_R1_endpoint`: **PASS**
- `R2_final20_occupancy_is_stationary`: **PASS**
- `R2_final15_occupancy_is_stationary`: **PASS**
- `R2_final10_occupancy_is_stationary`: **PASS**
- `R2_final10_net_change_is_at_most_5_waters`: **PASS**
- `R2_demonstrates_water_exchange`: **PASS**
- `R2_aperture_is_open`: **PASS**
- `R2_open_area_fraction_is_valid`: **PASS**
- `R2_minimum_CAP_OW_distance_is_safe`: **PASS**
- `R2_HBN_remained_frozen`: **PASS**
- `R2_PYR_remained_frozen`: **PASS**
- `R2_CAPS_remained_frozen`: **PASS**

## Decision

- Decision:
  **R2_SELECTED_AS_PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE**
- R1 status:
  **RETAIN_AS_CLOSED_STERIC_POSITIVE_CONTROL**
- R2 status:
  **PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE**
- R3 status:
  **DEFERRED_NOT_REQUIRED_AT_THIS_GATE**
- R4 status:
  **DEFERRED_NOT_REQUIRED_AT_THIS_GATE**
- Chemical-realization static design authorized:
  **YES**
- New MD execution authorized:
  **NO**
- Short mobile MD authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `BEGIN_R2_CHEMICALLY_REALIZABLE_END_RIM_DESIGN_GATE`

The next stage must replace the ideal steric partial cap with a
chemically defensible end-rim or terminal-ring realization while
preserving the validated aperture and confinement envelope.
