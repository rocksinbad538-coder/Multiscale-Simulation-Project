# Day 004 Extended 300 K Hold Failure Diagnostic

An extended 300 K NVT hold was attempted from the previously stable 300 K hold structure.

The simulation failed at step 796 with:

`Bond atoms 32158 32160 missing on proc 0`

Interpretation:
- The prior 10,000-step 300 K hold remains valid as a short-timescale stability test.
- The extended hold exposed a limitation of the current provisional Phase 0 setup.
- The likely causes are flexible water geometry, fixed non-periodic boundaries, lack of explicit axial/radial confinement walls, and/or excessive local displacement of a water molecule.
- This is a useful diagnostic, not a project failure. Longer production-like holds require a more robust water treatment and confinement protocol.

Next corrective actions:
1. inspect the failing atom/molecule;
2. rerun with a safer timestep and rigid water constraints;
3. add explicit confinement if the lumen is intended to be closed/contained;
4. treat the result as a protocol-stability boundary for Phase 0.
