# Day018 site-energy production QC

- Complete frames: 21
- Parsed calculations: 84
- Review threshold for observed steps: 40.0 meV
- Review threshold for S2-S1 gaps: 60.0 meV
- These thresholds are screening criteria only; they are not physical definitions of root switching.
- State identity cannot be established from excitation energy and oscillator strength alone. Transition character or NTO analysis remains necessary for definitive assignment.
- The source MD trajectory contains frozen h-BN and pyrene coordinates; the resulting variation is interpreted as solvent-induced diagonal disorder.

## Cluster statistics

cluster  n_frames  mean_S1_eV  std_S1_eV  min_S1_eV  max_S1_eV  mean_f1   std_f1  min_S2_minus_S1_meV  mean_S2_minus_S1_meV  mean_point_charges  min_point_charges  max_point_charges  maximum_observed_step_meV
   PYR2        21    4.006143   0.008260      3.992      4.019 0.011173 0.012688                 62.0             82.428571          163.428571                100                192                       27.0
   PYR3        21    4.003190   0.013359      3.969      4.019 0.009893 0.007724                 55.0             74.857143          166.666667                 76                192                       50.0
   PYR4        21    4.019238   0.012112      3.978      4.035 0.008602 0.006272                 55.0             71.380952          159.809524                 88                192                       53.0
   PYR5        21    3.776952   0.016430      3.750      3.797 0.648089 0.015305                 53.0             70.857143          175.809524                120                200                       43.0

## Frame-level diagonal disorder

 frame  time_ps  mean_PYR2_PYR4_eV  PYR5_eV  PYR5_redshift_vs_PYR2_PYR4_eV  PYR5_redshift_vs_PYR2_PYR4_meV  frame_min_eV  frame_max_eV  diagonal_spread_eV  diagonal_spread_meV
     0      0.0           4.004000    3.779                       0.225000                      225.000000         3.779         4.013               0.234                234.0
     1      5.0           4.007000    3.756                       0.251000                      251.000000         3.756         4.017               0.261                261.0
     2     10.0           4.014333    3.779                       0.235333                      235.333333         3.779         4.021               0.242                242.0
     3     15.0           4.012667    3.750                       0.262667                      262.666667         3.750         4.019               0.269                269.0
     4     20.0           4.017667    3.761                       0.256667                      256.666667         3.761         4.024               0.263                263.0
     5     25.0           4.011667    3.797                       0.214667                      214.666667         3.797         4.034               0.237                237.0
     6     30.0           4.021000    3.776                       0.245000                      245.000000         3.776         4.029               0.253                253.0
     7     35.0           4.000000    3.790                       0.210000                      210.000000         3.790         4.020               0.230                230.0
     8     40.0           4.007000    3.793                       0.214000                      214.000000         3.793         4.020               0.227                227.0
     9     45.0           4.012667    3.750                       0.262667                      262.666667         3.750         4.019               0.269                269.0
    10     50.0           4.004667    3.791                       0.213667                      213.666667         3.791         4.020               0.229                229.0
    11     55.0           4.001000    3.787                       0.214000                      214.000000         3.787         4.016               0.229                229.0
    12     60.0           4.013000    3.797                       0.216000                      216.000000         3.797         4.035               0.238                238.0
    13     65.0           4.019000    3.783                       0.236000                      236.000000         3.783         4.031               0.248                248.0
    14     70.0           3.999667    3.785                       0.214667                      214.666667         3.785         4.014               0.229                229.0
    15     75.0           4.014000    3.797                       0.217000                      217.000000         3.797         4.026               0.229                229.0
    16     80.0           4.001667    3.787                       0.214667                      214.666667         3.787         4.021               0.234                234.0
    17     85.0           4.011333    3.756                       0.255333                      255.333333         3.756         4.023               0.267                267.0
    18     90.0           4.016667    3.753                       0.263667                      263.666667         3.753         4.022               0.269                269.0
    19     95.0           4.002333    3.781                       0.221333                      221.333333         3.781         4.009               0.228                228.0
    20    100.0           4.008667    3.768                       0.240667                      240.666667         3.768         4.016               0.248                248.0

## Descriptive correlations

cluster  n  pearson_S1_vs_point_charges  pearson_S1_vs_f1  pearson_S1_vs_S2_minus_S1
   PYR2 21                     0.061617         -0.015358                  -0.195930
   PYR3 21                    -0.153146          0.143169                  -0.544802
   PYR4 21                     0.112156          0.076882                  -0.785217
   PYR5 21                    -0.103624          0.131774                  -0.454313

## State-identity review cases

 frame  time_ps cluster  S1_eV  S2_eV  S2_minus_S1_meV       f1  n_point_charges_orca  delta_time_ps  delta_S1_meV  step_review_flag  small_gap_review_flag
     1      5.0    PYR3  4.010  4.065             55.0 0.034274                   160            5.0          10.0             False                   True
     4     20.0    PYR5  3.761  3.815             54.0 0.648595                   164            5.0          11.0             False                   True
     5     25.0    PYR3  3.995  4.052             57.0 0.008856                   192            5.0         -19.0             False                   True
     5     25.0    PYR4  4.034  4.089             55.0 0.008463                   180            5.0          10.0             False                   True
     7     35.0    PYR3  3.969  4.070            101.0 0.012776                   192            5.0         -50.0              True                  False
     9     45.0    PYR5  3.750  3.840             90.0 0.623179                   168            5.0         -43.0              True                  False
    10     50.0    PYR5  3.791  3.859             68.0 0.663969                   152            5.0          41.0              True                  False
    13     65.0    PYR5  3.783  3.836             53.0 0.661990                   172            5.0         -14.0             False                   True
    14     70.0    PYR4  3.978  4.069             91.0 0.005349                   144            5.0         -53.0              True                  False
    14     70.0    PYR5  3.785  3.843             58.0 0.634868                   180            5.0           2.0             False                   True
    15     75.0    PYR4  4.026  4.095             69.0 0.018260                   172            5.0          48.0              True                  False
    17     85.0    PYR5  3.756  3.815             59.0 0.658948                   192            5.0         -31.0             False                   True

## Largest observed changes

 frame cluster  delta_time_ps  delta_S1_meV  S2_minus_S1_meV       f1  n_point_charges_orca
    14    PYR4            5.0         -53.0             91.0 0.005349                   144
     7    PYR3            5.0         -50.0            101.0 0.012776                   192
    15    PYR4            5.0          48.0             69.0 0.018260                   172
     9    PYR5            5.0         -43.0             90.0 0.623179                   168
    10    PYR5            5.0          41.0             68.0 0.663969                   152
     8    PYR3            5.0          38.0             66.0 0.017981                   148
     5    PYR5            5.0          36.0             74.0 0.647308                   188
    17    PYR3            5.0          35.0             77.0 0.015690                   164
    12    PYR3            5.0          33.0             74.0 0.001832                   184
    17    PYR5            5.0         -31.0             59.0 0.658948                   192
     3    PYR5            5.0         -29.0            102.0 0.615413                   176
    19    PYR5            5.0          28.0             70.0 0.649181                   192
