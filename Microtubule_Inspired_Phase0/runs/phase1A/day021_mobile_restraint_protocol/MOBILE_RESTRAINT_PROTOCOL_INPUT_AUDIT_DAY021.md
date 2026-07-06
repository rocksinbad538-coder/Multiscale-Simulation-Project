# Day021 Mobile Restraint-Protocol Input Audit

- Relevant files: 31
- TOP/ITP files: 26
- Existing position-restraint sections: 14
- Existing MDP files: 1
- Existing NDX files: 0
- Conditional/include directives: 55

## Final GRO velocities

- HBN: 0/1680 atoms with nonzero velocity; RMS speed 0.00000000 nm/ps.
- PYR: 0/104 atoms with nonzero velocity; RMS speed 0.00000000 nm/ps.
- SOL: 66536/66536 atoms with nonzero velocity; RMS speed 1.48292798 nm/ps.

## Protocol implication

The final frozen-solute GRO contains zero velocities for HBN and PYR but propagated water velocities. Direct continuation into a mobile-solute trajectory would therefore create a non-equilibrated kinetic partition. A staged protocol must explicitly regenerate or re-equilibrate velocities before full release.
