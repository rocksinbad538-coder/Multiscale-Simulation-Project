# Day018 pyrene structure and state-identity audit

- Jobs audited: 84
- Target orbital transition: 52a -> 53a
- Candidate roots: S1–S2

## Tracked-root counts

cluster  tracked_root  n_jobs
   PYR2             2      21
   PYR3             2      21
   PYR4             2      21
   PYR5             1      21

## Tracked-state statistics

cluster  n_jobs  mean_tracked_energy_eV  std_tracked_energy_eV  min_tracked_energy_eV  max_tracked_energy_eV  mean_tracked_fosc  mean_HOMO_LUMO_weight  min_HOMO_LUMO_weight  mean_state_separation_meV
   PYR2      21                4.088571               0.011084                  4.068                  4.104           0.493290               0.830708              0.738175                  82.428571
   PYR3      21                4.078048               0.011902                  4.052                  4.096           0.501778               0.836158              0.793849                  74.857143
   PYR4      21                4.090619               0.007586                  4.069                  4.104           0.510859               0.841515              0.823216                  71.380952
   PYR5      21                3.776952               0.016430                  3.750                  3.797           0.648089               0.863937              0.827679                  70.857143

## Structural equivalence summary

         maximum_direct_frame_RMSD_A  Kabsch_RMSD_all_atoms_vs_PYR2_A  Kabsch_RMSD_heavy_atoms_vs_PYR2_A  maximum_distance_matrix_difference_A  mean_planarity_RMSD_A
cluster                                                                                                                                                              
PYR2                             0.0                     1.406928e-15                       1.033853e-15                              0.000000               0.002598
PYR3                             0.0                     6.806661e-03                       7.351895e-03                              0.014537               0.002883
PYR4                             0.0                     8.359951e-03                       8.185885e-03                              0.014919               0.002774
PYR5                             0.0                     1.019262e-01                       9.047895e-02                              0.308059               0.009318

## Interpretation constraint

Root number alone is not used as state identity. The present tracking follows the low-lying state with the largest HOMO-to-LUMO configuration weight. Definitive diabatic assignment still requires orbital or NTO inspection for representative cases.
