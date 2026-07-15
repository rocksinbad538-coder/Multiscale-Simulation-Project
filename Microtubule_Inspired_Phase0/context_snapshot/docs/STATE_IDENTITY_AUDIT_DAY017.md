# Day017 pyrene structure and state-identity audit

- Jobs audited: 48
- Target orbital transition: 52a -> 53a
- Candidate roots: S1–S2

## Tracked-root counts

cluster  tracked_root  n_jobs
   PYR2             2      12
   PYR3             2      12
   PYR4             2      12
   PYR5             1      12

## Tracked-state statistics

cluster  n_jobs  mean_tracked_energy_eV  std_tracked_energy_eV  min_tracked_energy_eV  max_tracked_energy_eV  mean_tracked_fosc  mean_HOMO_LUMO_weight  min_HOMO_LUMO_weight  mean_state_separation_meV
   PYR2      12                4.086000               0.012828                  4.068                  4.104           0.488676               0.824557              0.738175                  79.083333
   PYR3      12                4.074417               0.012894                  4.052                  4.096           0.499139               0.834199              0.793849                  72.000000
   PYR4      12                4.090333               0.005348                  4.080                  4.099           0.510563               0.838811              0.824401                  69.333333
   PYR5      12                3.774167               0.016975                  3.750                  3.797           0.649540               0.866503              0.827679                  75.083333

## Structural equivalence summary

         maximum_direct_frame_RMSD_A  Kabsch_RMSD_all_atoms_vs_PYR2_A  Kabsch_RMSD_heavy_atoms_vs_PYR2_A  maximum_distance_matrix_difference_A  mean_planarity_RMSD_A
cluster                                                                                                                                                              
PYR2                             0.0                     1.406928e-15                       1.033853e-15                              0.000000               0.002598
PYR3                             0.0                     6.806661e-03                       7.351895e-03                              0.014537               0.002883
PYR4                             0.0                     8.359951e-03                       8.185885e-03                              0.014919               0.002774
PYR5                             0.0                     1.019262e-01                       9.047895e-02                              0.308059               0.009318

## Interpretation constraint

Root number alone is not used as state identity. The present tracking follows the low-lying state with the largest HOMO-to-LUMO configuration weight. Definitive diabatic assignment still requires orbital or NTO inspection for representative cases.
