# Day020 HBN Axial Architecture Audit

## Purpose

This audit reconstructs the actual axial architecture of the frozen HBN scaffold directly from atomic coordinates. It does not assume that the scaffold is one continuous cylindrical wall.

## Principal axis

- HBN axis: (0.00000000, 0.00000000, 1.00000000).
- HBN axial planes: 56.
- Typical neighboring-plane spacing: 0.073000 nm.
- Segment-break threshold: 0.350000 nm.

## Detected scaffold segments

- Segment 1: -3.045214 to 3.045786 nm; 1680 atoms; mean radius 1.199126 nm.

## Detected axial gaps

- No large axial gap was detected.

## Pyrene positions

- PYR_1: z=-2.008791 nm, r=1.649167 nm, class=inside_hbn_segment.
- PYR_2: z=-0.669599 nm, r=1.649090 nm, class=inside_hbn_segment.
- PYR_3: z=0.669593 nm, r=1.649295 nm, class=inside_hbn_segment.
- PYR_4: z=2.008786 nm, r=1.768256 nm, class=inside_hbn_segment.

## Consequence for water-region classification

Water regions must be defined using the detected HBN segments and gaps. The previous continuous-cylinder partition is retained only as a diagnostic calculation and must not be used as the final physical classification unless this audit detects a single continuous segment.
