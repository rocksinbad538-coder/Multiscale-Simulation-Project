# Day017 site-energy production QC

- Complete frames: 12
- Parsed calculations: 48
- Review threshold for observed steps: 40.0 meV
- Review threshold for S2-S1 gaps: 60.0 meV
- These thresholds are screening criteria only; they are not physical definitions of root switching.
- State identity cannot be established from excitation energy and oscillator strength alone. Transition character or NTO analysis remains necessary for definitive assignment.
- The source MD trajectory contains frozen h-BN and pyrene coordinates; the resulting variation is interpreted as solvent-induced diagonal disorder.

## Cluster statistics

cluster  n_frames  mean_S1_eV  std_S1_eV  min_S1_eV  max_S1_eV  mean_f1   std_f1  min_S2_minus_S1_meV  mean_S2_minus_S1_meV  mean_point_charges  min_point_charges  max_point_charges  maximum_observed_step_meV
   PYR2        12    4.006917   0.008554      3.994      4.017 0.015104 0.015583                 62.0             79.083333          158.666667                100                192                       23.0
   PYR3        12    4.002417   0.012887      3.969      4.019 0.010680 0.008623                 55.0             72.000000          160.666667                 76                192                       50.0
   PYR4        12    4.021000   0.005705      4.013      4.034 0.010542 0.006509                 55.0             69.333333          156.666667                 88                192                       10.0
   PYR5        12    3.774167   0.016975      3.750      3.797 0.649540 0.016581                 54.0             75.083333          168.666667                120                196                       43.0

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
    20    100.0           4.008667    3.768                       0.240667                      240.666667         3.768         4.016               0.248                248.0

## Descriptive correlations

cluster  n  pearson_S1_vs_point_charges  pearson_S1_vs_f1  pearson_S1_vs_S2_minus_S1
   PYR2 12                     0.265009         -0.097311                  -0.088273
   PYR3 12                    -0.230286          0.149320                  -0.469962
   PYR4 12                     0.517049         -0.181479                  -0.705302
   PYR5 12                    -0.078392          0.513425                  -0.502672

## State-identity review cases

 frame  time_ps cluster  S1_eV  S2_eV  S2_minus_S1_meV       f1  n_point_charges_orca  delta_time_ps  delta_S1_meV  step_review_flag  small_gap_review_flag
     1      5.0    PYR3  4.010  4.065             55.0 0.034274                   160            5.0          10.0             False                   True
     4     20.0    PYR5  3.761  3.815             54.0 0.648595                   164            5.0          11.0             False                   True
     5     25.0    PYR3  3.995  4.052             57.0 0.008856                   192            5.0         -19.0             False                   True
     5     25.0    PYR4  4.034  4.089             55.0 0.008463                   180            5.0          10.0             False                   True
     7     35.0    PYR3  3.969  4.070            101.0 0.012776                   192            5.0         -50.0              True                  False
     9     45.0    PYR5  3.750  3.840             90.0 0.623179                   168            5.0         -43.0              True                  False
    10     50.0    PYR5  3.791  3.859             68.0 0.663969                   152            5.0          41.0              True                  False

## Largest observed changes

 frame cluster  delta_time_ps  delta_S1_meV  S2_minus_S1_meV       f1  n_point_charges_orca
     7    PYR3            5.0         -50.0            101.0 0.012776                   192
     9    PYR5            5.0         -43.0             90.0 0.623179                   168
    10    PYR5            5.0          41.0             68.0 0.663969                   152
     8    PYR3            5.0          38.0             66.0 0.017981                   148
     5    PYR5            5.0          36.0             74.0 0.647308                   188
     3    PYR5            5.0         -29.0            102.0 0.615413                   176
     6    PYR3            5.0          24.0             77.0 0.011422                   156
     9    PYR2            5.0          23.0             77.0 0.008134                   176
    20    PYR5           50.0         -23.0             80.0 0.670424                   188
     1    PYR5            5.0         -23.0             83.0 0.660520                   180
     2    PYR5            5.0          23.0             68.0 0.637472                   172
     6    PYR5            5.0         -21.0             69.0 0.659520                   164
