# Day019 vacuum-reference and solvent-shift analysis

## Vacuum state tracking

cluster  tracked_root  tracked_energy_eV  tracked_fosc  tracked_HOMO_LUMO_weight  alternate_root  alternate_energy_eV  alternate_HOMO_LUMO_weight  state_separation_meV
   PYR2             2              4.090      0.505315                  0.846347               1                4.001                           0                  89.0
   PYR3             2              4.081      0.514411                  0.848738               1                4.001                           0                  80.0
   PYR4             2              4.089      0.518148                  0.849718               1                4.015                           0                  74.0
   PYR5             1              3.780      0.664486                  0.877731               2                3.846                           0                  66.0

## Solvent-shift statistics

cluster  n_frames  vacuum_energy_eV  mean_embedded_energy_eV  mean_solvent_shift_meV  std_solvent_shift_meV  min_solvent_shift_meV  max_solvent_shift_meV
   PYR2        21             4.090                 4.088571               -1.428571              11.084094                  -22.0                   14.0
   PYR3        21             4.081                 4.078048               -2.952381              11.901581                  -29.0                   15.0
   PYR4        21             4.089                 4.090619                1.619048               7.586015                  -20.0                   15.0
   PYR5        21             3.780                 3.776952               -3.047619              16.430083                  -30.0                   17.0

## PYR5 offset decomposition

                                quantity     value_meV
vacuum_geometry_offset_PYR5_vs_PYR2_PYR4  3.066667e+02
mean_solvent_contribution_to_PYR5_offset  2.126984e+00
               mean_embedded_PYR5_offset  3.087937e+02
             decomposition_closure_error -1.096456e-12

The decomposition is defined as:

`embedded PYR5 offset = vacuum geometry offset + mean differential solvent shift`.

The vacuum calculation removes the TIP4P/2005 point-charge environment while preserving the frozen chromophore geometry. The vacuum offset therefore measures the fixed-geometry baseline difference, subject to residual numerical orientation or integration-grid effects. The embedded-minus-vacuum difference measures the electrostatic effect of the water point charges for each fixed chromophore geometry.

## NTO cases selected

calculation_type  frame cluster                            job  tracked_root                                reason  metric_value
        embedded      3    PYR5        frame003_PYR5_embedding             1            minimum_PYR5_solvent_shift    -30.000000
        embedded      5    PYR2        frame005_PYR2_embedding             2      minimum_tracked_character_weight      0.738175
        embedded      5    PYR5        frame005_PYR5_embedding             1            maximum_PYR5_solvent_shift     17.000000
        embedded     13    PYR5        frame013_PYR5_embedding             1          minimum_S1_S2_separation_meV     53.000000
vacuum_reference      0    PYR2 frame000_PYR2_vacuum_reference             2 vacuum_reference_for_each_chromophore      0.846347
vacuum_reference      0    PYR3 frame000_PYR3_vacuum_reference             2 vacuum_reference_for_each_chromophore      0.848738
vacuum_reference      0    PYR4 frame000_PYR4_vacuum_reference             2 vacuum_reference_for_each_chromophore      0.849718
vacuum_reference      0    PYR5 frame000_PYR5_vacuum_reference             1 vacuum_reference_for_each_chromophore      0.877731
