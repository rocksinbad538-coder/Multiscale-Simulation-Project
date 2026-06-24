# Day015 Gamma-Sigma Dynamic-Disorder Map Summary

## Static reference

Reference case: gamma_phi = 1 ps^-1, sigma_E = 0 meV.

- mean max PYR2 = 0.305615
- mean t(PYR2 > 30%) = 1350.0 fs

## Main result

The map shows that transport enhancement is driven primarily by site-energy dynamic disorder amplitude rather than by additional Lindblad pure dephasing.

The strongest enhancement occurs for sigma_E around 20 meV. Enhancement is already present even when gamma_phi = 0 ps^-1, indicating that dynamic Hamiltonian fluctuations themselves provide an effective decoherence/modulation pathway.

## Best cases by maximum PYR2 population

 gamma_phi_ps  sigma_E_meV  mean_max_PYR2  mean_t30_fs  delta_max_PYR2_vs_static  t30_speedup_factor
          0.0           20       0.343320   732.500000                  0.037705            1.843003
          0.1           20       0.342794   717.500000                  0.037178            1.881533
          0.5           20       0.339653   827.500000                  0.034037            1.631420
          1.0           20       0.339067   807.500000                  0.033452            1.671827
         50.0           20       0.338322   805.833333                  0.032707            1.675284
          0.5           10       0.337771   980.833333                  0.032155            1.376381
          5.0           20       0.336740   845.833333                  0.031124            1.596059
         10.0           20       0.336226   830.833333                  0.030611            1.624875
          0.1           10       0.335909   953.333333                  0.030293            1.416084
          0.0           10       0.335438   905.833333                  0.029823            1.490340

## Best cases by fastest t(PYR2 > 30%)

 gamma_phi_ps  sigma_E_meV  mean_max_PYR2  mean_t30_fs  delta_t30_fs_vs_static  t30_speedup_factor
          0.1           20       0.342794   717.500000             -632.500000            1.881533
          0.0           20       0.343320   732.500000             -617.500000            1.843003
         50.0           20       0.338322   805.833333             -544.166667            1.675284
          1.0           20       0.339067   807.500000             -542.500000            1.671827
         25.0           20       0.335199   824.166667             -525.833333            1.638018
          0.5           20       0.339653   827.500000             -522.500000            1.631420
         10.0           20       0.336226   830.833333             -519.166667            1.624875
          5.0           20       0.336740   845.833333             -504.166667            1.596059
          0.0           10       0.335438   905.833333             -444.166667            1.490340
          1.0           10       0.333132   925.833333             -424.166667            1.458146

## Interpretation

The result supports an ENAQT-like mechanism, but specifically one dominated by site-energy fluctuations. In this regime, moderate dynamic energetic disorder helps bridge the weak PYR2-PYR3 bottleneck.

## Caveat

The disorder remains synthetic. The next physical step is to determine whether real MD-derived site-energy fluctuations fall in the 10-20 meV amplitude range with sub-100 fs to few-hundred-fs correlation times.
