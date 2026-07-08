# R1 NVT 20 ps Thermal-Transient Review

## Purpose

The original 20 ps screening completed without dynamical instability
and failed only the full-trajectory temperature-standard-deviation
gate.

This review determines whether the failure is confined to initial
thermalization or persists after equilibration.

## Whole-trajectory temperature

- First temperature:
  **300.056244 K**
- Second temperature:
  **179.039551 K**
- Whole-trajectory mean:
  **298.467390 K**
- Whole-trajectory standard deviation:
  **10.992822 K**
- Standard deviation after removing only the first point:
  **11.006300 K**
- Canonical temperature-fluctuation estimate from
  99303 solvent degrees of freedom:
  **1.346341 K**

## Temporal windows

- 0.00–20.00 ps: T=298.4674 ± 10.9928 K; range=179.0396–303.4521 K; slope=0.483427 K/ps; potential slope=140.453 kJ mol^-1 ps^-1
- 0.05–20.00 ps: T=298.4634 ± 11.0063 K; range=179.0396–303.4521 K; slope=0.488257 K/ps; potential slope=39.442 kJ mol^-1 ps^-1
- 0.10–20.00 ps: T=298.7627 ± 9.2476 K; range=188.0844–303.4521 K; slope=0.401694 K/ps; potential slope=-10.047 kJ mol^-1 ps^-1
- 0.25–20.00 ps: T=299.4474 ± 4.6238 K; range=249.1096–303.4521 K; slope=0.201674 K/ps; potential slope=-95.043 kJ mol^-1 ps^-1
- 0.50–20.00 ps: T=299.8987 ± 2.0430 K; range=281.8354–303.4521 K; slope=0.068951 K/ps; potential slope=-159.354 kJ mol^-1 ps^-1
- 1.00–20.00 ps: T=300.1168 ± 1.3172 K; range=296.0628–303.4521 K; slope=0.003402 K/ps; potential slope=-179.019 kJ mol^-1 ps^-1
- 2.00–20.00 ps: T=300.1337 ± 1.3244 K; range=296.0628–303.4521 K; slope=-0.002027 K/ps; potential slope=-130.592 kJ mol^-1 ps^-1
- 5.00–20.00 ps: T=300.1487 ± 1.3613 K; range=296.0628–303.4521 K; slope=-0.007659 K/ps; potential slope=-58.403 kJ mol^-1 ps^-1
- 10.00–20.00 ps: T=300.1393 ± 1.3364 K; range=296.0628–303.2159 K; slope=-0.026557 K/ps; potential slope=64.025 kJ mol^-1 ps^-1

## Confinement after 5 ps

- Minimum lumen occupancy:
  **427 waters**
- Minimum retained initially luminal waters:
  **427 waters**
- Endpoint occupancy:
  **428 waters**
- Endpoint initially luminal waters retained:
  **428 waters**
- Minimum CAP–OW distance:
  **0.170798 nm**

## Gates

- `original_mdrun_return_code_zero`: **PASS**
- `original_mdrun_finished`: **PASS**
- `original_failure_is_temperature_std_only`: **PASS**
- `post_5ps_temperature_mean_is_295_to_305K`: **PASS**
- `post_5ps_temperature_std_is_at_most_5K`: **PASS**
- `post_5ps_temperature_min_is_at_least_280K`: **PASS**
- `post_5ps_temperature_max_is_at_most_320K`: **PASS**
- `last_10ps_temperature_mean_is_295_to_305K`: **PASS**
- `last_10ps_temperature_std_is_at_most_5K`: **PASS**
- `last_10ps_temperature_slope_is_small`: **PASS**
- `post_5ps_lumen_occupancy_remains_above_90_percent`: **PASS**
- `post_5ps_initial_lumen_retention_remains_above_90_percent`: **PASS**
- `post_5ps_CAP_OW_distance_remains_above_limit`: **PASS**
- `endpoint_lumen_occupancy_is_428`: **PASS**
- `endpoint_initial_lumen_retention_is_428`: **PASS**

## Decision

- Decision:
  **R1_INITIAL_THERMALIZATION_TRANSIENT_CONFIRMED**
- Failed gates:
  **NONE**
- Preparation of a checkpoint continuation to 50 ps authorized:
  **YES**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `PREPARE_R1_30PS_CHECKPOINT_CONTINUATION_TO_50PS`

A passing result authorizes only preparation of a 30 ps continuation
from the existing 20 ps checkpoint. It does not authorize velocity
regeneration, a new independent trajectory, mobile-solute MD, or
multitemperature production.
