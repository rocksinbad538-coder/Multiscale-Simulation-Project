# Day018 ORCA embedding production audit

- Jobs parsed: 84
- Fully successful embedded TDDFT jobs: 84/84
- Normal ORCA terminations: 84/84
- SCF converged: 84/84
- TDDFT/TDA completed: 84/84
- Explicit error signatures: 0
- Point charges read by ORCA: 76–200
- S1 range: 3.750–4.035 eV
- Mean S1 across all sites and frames: 3.951 eV
- Standard deviation across all sites and frames: 0.102 eV

## Parsed jobs

 frame cluster  terminated_normally  scf_converged  tddft_finished  has_error_flag  n_point_charges_orca  S1_eV  S2_eV  S3_eV       f1  dipole_D  total_runtime_min
     0    PYR2                 True           True            True           False                   100  3.999  4.068  4.681 0.001619  1.231000           5.083333
     0    PYR3                 True           True            True           False                    76  4.000  4.079  4.691 0.007855  0.916531           8.550000
     0    PYR4                 True           True            True           False                    88  4.013  4.090  4.653 0.015613  0.490451           8.666667
     0    PYR5                 True           True            True           False                   120  3.779  3.845  4.347 0.652758  1.552144           8.433333
     1    PYR2                 True           True            True           False                   168  3.994  4.083  4.662 0.031321  1.448856           8.866667
     1    PYR3                 True           True            True           False                   160  4.010  4.065  4.646 0.034274  2.149699           8.900000
     1    PYR4                 True           True            True           False                   120  4.017  4.089  4.616 0.001991  0.325371           8.833333
     1    PYR5                 True           True            True           False                   180  3.756  3.839  4.303 0.660520  2.234151           8.750000
     2    PYR2                 True           True            True           False                   140  4.009  4.094  4.588 0.008972  1.380070           8.650000
     2    PYR3                 True           True            True           False                   180  4.013  4.086  4.700 0.011531  1.943994           8.733333
     2    PYR4                 True           True            True           False                   156  4.021  4.093  4.637 0.009617  1.066271           8.583333
     2    PYR5                 True           True            True           False                   172  3.779  3.847  4.431 0.637472  1.044397           9.066667
     3    PYR2                 True           True            True           False                   168  4.015  4.104  4.660 0.008823  1.118632           8.583333
     3    PYR3                 True           True            True           False                   168  4.004  4.071  4.598 0.002359  2.065674           8.950000
     3    PYR4                 True           True            True           False                   152  4.019  4.097  4.678 0.008493  1.152261           8.583333
     3    PYR5                 True           True            True           False                   176  3.750  3.852  4.394 0.615413  2.368216           8.566667
     4    PYR2                 True           True            True           False                   164  4.015  4.101  4.701 0.011137  0.576480           8.533333
     4    PYR3                 True           True            True           False                   164  4.014  4.095  4.612 0.007055  1.078524           8.500000
     4    PYR4                 True           True            True           False                   156  4.024  4.087  4.631 0.002471  0.757645           8.583333
     4    PYR5                 True           True            True           False                   164  3.761  3.815  4.318 0.648595  2.366721           8.483333
     5    PYR2                 True           True            True           False                   144  4.006  4.068  4.711 0.058606  1.996653           8.633333
     5    PYR3                 True           True            True           False                   192  3.995  4.052  4.643 0.008856  2.339854           8.816667
     5    PYR4                 True           True            True           False                   180  4.034  4.089  4.722 0.008463  1.665576           8.600000
     5    PYR5                 True           True            True           False                   188  3.797  3.871  4.351 0.647308  1.216905           8.416667
     6    PYR2                 True           True            True           False                   140  4.015  4.081  4.666 0.018233  1.975306           9.000000
     6    PYR3                 True           True            True           False                   156  4.019  4.096  4.742 0.011422  1.548970           8.566667
     6    PYR4                 True           True            True           False                   172  4.029  4.092  4.702 0.010113  1.594295           8.850000
     6    PYR5                 True           True            True           False                   164  3.776  3.845  4.267 0.659520  2.990339           8.550000
     7    PYR2                 True           True            True           False                   192  4.011  4.103  4.638 0.009094  1.795061           9.250000
     7    PYR3                 True           True            True           False                   192  3.969  4.070  4.559 0.012776  0.759346           9.133333
     7    PYR4                 True           True            True           False                   164  4.020  4.093  4.608 0.007442  1.519905           9.216667
     7    PYR5                 True           True            True           False                   196  3.790  3.859  4.381 0.659437  0.703119           8.733333
     8    PYR2                 True           True            True           False                   160  3.994  4.079  4.648 0.005736  0.970874           9.183333
     8    PYR3                 True           True            True           False                   148  4.007  4.073  4.655 0.017981  1.589176           9.300000
     8    PYR4                 True           True            True           False                   160  4.020  4.092  4.720 0.016728  0.328742           8.983333
     8    PYR5                 True           True            True           False                   156  3.793  3.871  4.358 0.655884  1.132189           8.766667
     9    PYR2                 True           True            True           False                   176  4.017  4.094  4.647 0.008134  1.396384           9.016667
     9    PYR3                 True           True            True           False                   164  4.002  4.075  4.687 0.004027  1.397275           8.700000
     9    PYR4                 True           True            True           False                   168  4.019  4.080  4.688 0.025503  1.544804           8.866667
     9    PYR5                 True           True            True           False                   168  3.750  3.840  4.309 0.623179  2.101997           8.466667
    10    PYR2                 True           True            True           False                   164  3.998  4.079  4.646 0.008701  0.789775           8.866667
    10    PYR3                 True           True            True           False                   160  3.996  4.065  4.654 0.004807  2.585030           8.900000
    10    PYR4                 True           True            True           False                   172  4.020  4.099  4.741 0.013226  0.652228           8.566667
    10    PYR5                 True           True            True           False                   152  3.791  3.859  4.398 0.663969  1.004214           8.700000
    11    PYR2                 True           True            True           False                   152  4.008  4.098  4.662 0.014138  0.490893           8.533333
    11    PYR3                 True           True            True           False                   164  3.979  4.071  4.590 0.007325  0.725963           8.566667
    11    PYR4                 True           True            True           False                   184  4.016  4.090  4.624 0.004866  1.138837           8.633333
    11    PYR5                 True           True            True           False                   188  3.787  3.854  4.365 0.659752  0.605486           8.483333
    12    PYR2                 True           True            True           False                   164  3.992  4.079  4.652 0.003766  1.517884           8.833333
    12    PYR3                 True           True            True           False                   184  4.012  4.086  4.624 0.001832  1.583738           8.983333
    12    PYR4                 True           True            True           False                   132  4.035  4.104  4.578 0.008868  0.974499           8.833333
    12    PYR5                 True           True            True           False                   160  3.797  3.870  4.400 0.644719  1.476674           9.066667
    13    PYR2                 True           True            True           False                   156  4.019  4.103  4.695 0.002870  0.978905           8.700000
    13    PYR3                 True           True            True           False                   152  4.007  4.074  4.644 0.020813  1.398414           8.683333
    13    PYR4                 True           True            True           False                   140  4.031  4.102  4.714 0.001159  0.923350           8.683333
    13    PYR5                 True           True            True           False                   172  3.783  3.836  4.387 0.661990  0.875973           8.716667
    14    PYR2                 True           True            True           False                   176  4.007  4.092  4.644 0.001844  1.247468           8.783333
    14    PYR3                 True           True            True           False                   180  4.014  4.096  4.600 0.010113  1.454747           9.200000
    14    PYR4                 True           True            True           False                   144  3.978  4.069  4.596 0.005349  1.906278           8.866667
    14    PYR5                 True           True            True           False                   180  3.785  3.843  4.320 0.634868  1.779593           9.216667
    15    PYR2                 True           True            True           False                   180  4.009  4.099  4.599 0.008699  0.592796           9.233333
    15    PYR3                 True           True            True           False                   172  4.007  4.087  4.626 0.013688  1.403129           8.733333
    15    PYR4                 True           True            True           False                   172  4.026  4.095  4.611 0.018260  2.121051           8.933333
    15    PYR5                 True           True            True           False                   188  3.797  3.861  4.390 0.617290  0.821102           8.550000
    16    PYR2                 True           True            True           False                   176  4.005  4.095  4.598 0.007981  1.794603           9.166667
    16    PYR3                 True           True            True           False                   192  3.979  4.073  4.567 0.004160  1.467141           8.983333
    16    PYR4                 True           True            True           False                   164  4.021  4.088  4.659 0.000751  1.524839           8.983333
    16    PYR5                 True           True            True           False                   196  3.787  3.848  4.284 0.640872  1.472437           8.683333
    17    PYR2                 True           True            True           False                   184  3.997  4.090  4.582 0.003951  0.640084           8.850000
    17    PYR3                 True           True            True           False                   164  4.014  4.091  4.671 0.015690  1.127828           8.616667
    17    PYR4                 True           True            True           False                   192  4.023  4.095  4.662 0.004235  1.135567           8.566667
    17    PYR5                 True           True            True           False                   192  3.756  3.815  4.317 0.658948  1.746650           8.816667
    18    PYR2                 True           True            True           False                   164  4.011  4.087  4.687 0.004741  1.719813           8.850000
    18    PYR3                 True           True            True           False                   172  4.017  4.089  4.652 0.003591  0.865701           8.566667
    18    PYR4                 True           True            True           False                   172  4.022  4.088  4.637 0.005228  1.905844           8.633333
    18    PYR5                 True           True            True           False                   200  3.753  3.835  4.308 0.647780  1.705724           8.933333
    19    PYR2                 True           True            True           False                   176  3.998  4.085  4.632 0.005388  0.624705           8.833333
    19    PYR3                 True           True            True           False                   192  4.009  4.079  4.713 0.002384  1.827814           8.716667
    19    PYR4                 True           True            True           False                   176  4.000  4.088  4.562 0.005415  1.002612           8.616667
    19    PYR5                 True           True            True           False                   192  3.781  3.851  4.358 0.649181  0.755600           8.683333
    20    PYR2                 True           True            True           False                   188  4.010  4.078  4.693 0.010874  2.030039           8.833333
    20    PYR3                 True           True            True           False                   168  4.000  4.066  4.668 0.005213  1.367918           8.716667
    20    PYR4                 True           True            True           False                   192  4.016  4.083  4.647 0.006844  1.819395           8.916667
    20    PYR5                 True           True            True           False                   188  3.768  3.848  4.366 0.670424  1.193734           8.400000
