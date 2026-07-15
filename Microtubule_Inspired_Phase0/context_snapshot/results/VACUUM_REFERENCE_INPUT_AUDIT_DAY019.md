# Day019 vacuum-reference input audit

## Scope

Four vacuum-reference ORCA inputs were generated from the frame000 embedded inputs for PYR2âPYR5.

## Preserved settings

- Method line: `! wB97X-D3 def2-SVP def2/J RIJCOSX TightSCF`
- `%maxcore 4096`
- `nroots 10`
- `tda true`
- Charge/multiplicity: `0 1`
- Geometry: exactly 26 atoms, C16H10

## Deliberate modification

The only removed input line is the corresponding `%pointcharges` directive. No coordinates or electronic-structure settings were changed.

## Generated jobs

| Chromophore | Vacuum-reference job |
|---|---|
| PYR2 | `frame000_PYR2_vacuum_reference` |
| PYR3 | `frame000_PYR3_vacuum_reference` |
| PYR4 | `frame000_PYR4_vacuum_reference` |
| PYR5 | `frame000_PYR5_vacuum_reference` |

## Validation status

- Inputs generated: 4/4
- Geometry identity checks: 4/4 passed
- Composition checks: 4/4 passed
- Charge/multiplicity checks: 4/4 passed
- Residual `%pointcharges` directives: 0
- Unexpected input differences: 0
