# Day021 Working Frozen-Topology Semantic Identity

## Source-level comparison

- Raw TOP equality: False.
- Comment/whitespace-normalized equality: False.
- Unified-diff lines: 13.

## Processed TPR comparison

- `grompp`: PASS.
- `inputrec` identity: True.
- Processed topology identity: False.
- Box identity: True.

## Decision

- Working-copy semantic identity: FAIL.
- The accepted-directory topology remains the canonical provenance source regardless of working-copy equivalence.
