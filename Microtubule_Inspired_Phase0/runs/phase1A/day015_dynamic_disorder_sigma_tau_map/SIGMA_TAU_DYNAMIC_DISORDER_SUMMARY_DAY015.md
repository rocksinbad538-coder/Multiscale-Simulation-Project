# Day015 Sigma-Tau Dynamic-Disorder Map Summary

## Static reference

- mean max PYR2 = 0.305615
- mean t(PYR2 > 30%) = 1350.0 fs

## Main result

The site-energy dynamic-disorder scan identifies a robust transport-enhancement regime.

The strongest enhancement occurs for site-energy fluctuations of approximately 10-20 meV and correlation times of approximately 0.02-0.10 ps.

## Best cases by maximum PYR2 population

 sigma_E_meV  tau_c_ps  mean_max_PYR2  mean_t30_fs  delta_max_PYR2_vs_static  t30_speedup_factor
          20      0.02       0.341821      811.250                  0.036205            1.664099
          20      0.10       0.339574      821.250                  0.033959            1.643836
          20      0.05       0.337408      876.875                  0.031792            1.539558
          20      1.00       0.336757     1161.250                  0.031141            1.162540
          20      0.20       0.336687      938.125                  0.031072            1.439041
          10      0.02       0.335198      951.250                  0.029583            1.419185
          20      0.50       0.334069      932.500                  0.028453            1.447721
          10      0.20       0.334059     1133.750                  0.028443            1.190739
          10      0.50       0.332209     1163.125                  0.026593            1.160666
          10      1.00       0.330097     1463.750                  0.024482            0.922289

## Best cases by fastest t(PYR2 > 30%)

 sigma_E_meV  tau_c_ps  mean_max_PYR2  mean_t30_fs  delta_t30_fs_vs_static  t30_speedup_factor
          20      0.02       0.341821      811.250                -538.750            1.664099
          20      0.10       0.339574      821.250                -528.750            1.643836
          10      0.10       0.328986      825.625                -524.375            1.635125
          20      0.05       0.337408      876.875                -473.125            1.539558
          20      0.50       0.334069      932.500                -417.500            1.447721
          20      0.20       0.336687      938.125                -411.875            1.439041
          10      0.02       0.335198      951.250                -398.750            1.419185
          10      0.05       0.328357     1070.625                -279.375            1.260946
           5      0.50       0.315901     1106.250                -243.750            1.220339
          10      0.20       0.334059     1133.750                -216.250            1.190739

## Interpretation

The static model already reaches PYR2, but dynamic site-energy fluctuations moderately improve PYR2 access. This is consistent with an ENAQT-like regime where environmental fluctuations help bridge the weak PYR2-PYR3 bottleneck.

The effect is not monotonic in correlation time. Fast and intermediate fluctuations are more beneficial than slow fluctuations, suggesting that slow disorder behaves more like transient energetic inhomogeneity than efficient transport-assisting noise.

## Caveat

This is still a synthetic dynamic-disorder scan. The next physical step is to replace the Ornstein-Uhlenbeck parameters with MD-derived site-energy and coupling fluctuation statistics.
