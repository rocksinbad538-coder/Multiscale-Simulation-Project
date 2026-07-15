# Day 003 Stability Summary — BN-like Scaffold + 30,000 Confined Waters

System:

- Scaffold atoms: 17,280
- Water molecules: 30,000
- Water atoms: 90,000
- Total atoms: 107,280
- Scaffold: fixed during water equilibration
- Water model/force field: provisional Phase 0 placeholder
- Objective: numerical stability and workflow validation

Completed stages:

| Stage | Target / final water temperature | Runtime | Stability |
|---|---:|---:|---|
| Corrected NVT | ~50 K | 5,000 steps | Stable |
| Ramp | 50 → 150 K | 10,000 steps | Stable |
| Ramp | 150 → 200 K | 10,000 steps | Stable |
| Ramp | 200 → 250 K | 10,000 steps | Stable |
| Ramp | 250 → 300 K | 10,000 steps | Stable |
| Hold | 300 K | 10,000 steps | Stable |

Final 300 K hold metrics:

| Metric | Value |
|---|---:|
| Final water temperature | 299.74 K |
| Water temperature range during hold | 298.06–302.12 K |
| Final water MSD | 12.349 Å² |
| Final fraction inside nominal lumen segment | 0.9983 |
| Final fraction radially outside nominal lumen radius | 0.00163 |
| Final fraction axially outside nominal segment | 0.000133 |
| Final fraction outside outer scaffold radius | 0.0 |
| Dangerous builds | 0 |
| NaNs / lost atoms | None detected |

Interpretation:

The hydrated BN-like tubular scaffold remained numerically stable through the full 50 K to 300 K ramp and a subsequent short 300 K hold. Water remained overwhelmingly confined inside the nominal lumen. These results validate the current Phase 0 workflow and system construction for short-timescale testing. They should not yet be interpreted as final predictive material behavior because the present force field is provisional.
