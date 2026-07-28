# QM_F06 UPPER V7-A R1 Post-QM Audit

This directory is reserved for the independent final post-QM
acceptance audit of the QM_F06 UPPER V7-A R1 optimization.

The audit must not authorize RESP unless all of the following are
confirmed:

- ORCA terminated normally.
- The geometry optimization converged.
- Atom identity and ordering are preserved.
- The system contains 52 atoms with composition B17N14H21.
- The 12 fixed atoms remain within their accepted tolerance.
- All 57 nominal edges remain present and within range.
- No nonnominal geometric reconnectivity is detected.
- No atom is overcoordinated.
- The system remains a single connected component.
- No hard contacts are present.
- The strengthened global and modified-region near-contact gates pass.
- Required V7-A canonical edges and boundary caps remain intact.
- The obsolete HCAPV2 boundary cap remains absent.
- The V6-B failure mode remains eliminated.
- The optimization trajectory is parseable and structurally acceptable.

RESP, force-field adoption, and molecular dynamics remain prohibited
until the complete post-QM gate passes.
