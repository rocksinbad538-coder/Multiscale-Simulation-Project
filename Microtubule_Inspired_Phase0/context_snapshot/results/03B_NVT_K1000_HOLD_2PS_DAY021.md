# Day021 03b_nvt_k1000_hold_2ps

- Decision: **PASS**
- Temperature mean/last: 296.2256/298.5658 K
- Temperature range: 293.2257–299.0176 K
- Potential energy first/last: -768667.6875/-766014.0625 kJ/mol
- HBN incremental RMS/max: 0.02203817/0.11379807 nm
- PYR incremental RMS/max: 0.04088974/0.07339619 nm
- HBN cumulative RMS/max: 0.02354707/0.14006784 nm
- PYR cumulative RMS/max: 0.02620518/0.05522681 nm
- Instability signatures: none

No subsequent protocol stage was executed.

## Structural reconciliation

The maximum isolated HBN displacement was not treated as
a structural failure because the bonded-network diagnostics
showed preserved geometry:

- maximum bond-equilibrium deviation: 0.00921405 nm;
- q99 bond-equilibrium deviation: 0.00704712 nm;
- maximum angle-equilibrium deviation: 3.4410 degrees;
- no LINCS, SETTLE, fatal, or non-finite runtime values;
- only 5 of the Stage03 top-20 displacement atoms persisted
  in the Stage03b top-20 set.

Scientific decision: **PASS**.

- Rerun required: **NO**
- Stage04 authorized: **YES**
