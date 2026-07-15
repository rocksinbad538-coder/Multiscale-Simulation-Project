# PYR4-PYR5 Hydrated Dimer Orbital Cube Summary

Generated orbital cube files from the converged hydrated PYR4-PYR5 TDDFT pilot.

## Files

- `PYR4_PYR5_water0p50_serial_tight.mo315a.cube`
- `PYR4_PYR5_water0p50_serial_tight.mo316a.cube`
- `PYR4_PYR5_water0p50_serial_tight.xyz`

## State relevance

The anomalous low-energy S1 state at 2.959 eV is dominated by:

- MO 315a -> MO 316a
- weight = 0.999988

Therefore, these two cube files are the first diagnostic objects needed to determine whether S1 is localized on pyrene, water, or has charge-transfer-like character.

## Immediate interpretation

If MO315 and MO316 are both pyrene-localized, S1 may be a genuine low-energy chromophore/dimer excitation. If one orbital is water-localized or spatially separated from the other, S1 should be treated as environment-induced or charge-transfer-like and excluded from simple excitonic splitting estimates.
