# Day015 Haken-Strobl Dynamics Audit

## Purpose

Audit trace conservation, Hermiticity, positivity, population conservation, coherence decay, and early PYR2 arrival times for the hydrated four-site Hamiltonian.

## Interpretation rule

The long-time uniform population reached by the pure-dephasing Haken-Strobl model should be treated as a model artifact, not as thermodynamic relaxation. Early-time transfer metrics are the physically useful diagnostic at this stage.

## Summary table

| gamma_ps | solver_success | max_trace_error | max_hermiticity_error | max_population_error | min_density_eigenvalue | final_PYR2 | final_PYR3 | final_PYR4 | final_PYR5 | final_coherence_l1 | final_purity | t_PYR2_1pct_fs | t_PYR2_5pct_fs | t_PYR2_10pct_fs | max_PYR2 | time_max_PYR2_fs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | True | 3.66374e-15 | 0 | 3.55271e-15 | -2.46734e-08 | 0.00250861 | 0.332898 | 0.647856 | 0.0167369 | 1.43774 | 1 | 70 | nan | nan | 0.0223334 | 4430 |
| 1 | True | 2.22045e-15 | 0 | 2.22045e-15 | 0 | 0.0698086 | 0.307231 | 0.313014 | 0.309947 | 0.098283 | 0.294658 | 70 | 6550 | nan | 0.0698086 | 10000 |
| 5 | True | 1.55431e-15 | 0 | 1.55431e-15 | 0 | 0.194044 | 0.267965 | 0.269221 | 0.268771 | 0.0298075 | 0.254303 | 80 | 1370 | 3320 | 0.194044 | 10000 |
| 10 | True | 8.88178e-16 | 0 | 8.88178e-16 | 0 | 0.234131 | 0.255087 | 0.255451 | 0.255331 | 0.00811165 | 0.250345 | 80 | 760 | 1810 | 0.234131 | 10000 |
| 50 | True | 1.33227e-15 | 0 | 1.33227e-15 | 0 | 0.247097 | 0.250931 | 0.250988 | 0.250984 | 0.000724225 | 0.250011 | 140 | 550 | 1190 | 0.247097 | 10000 |

## Files

- `haken_strobl_audit_summary.csv`
- `haken_strobl_audit_gamma0ps.csv`
- `haken_strobl_audit_gamma1ps.csv`
- `haken_strobl_audit_gamma5ps.csv`
- `haken_strobl_audit_gamma10ps.csv`
- `haken_strobl_audit_gamma50ps.csv`