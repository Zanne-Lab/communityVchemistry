Does chemistry or community better predict mass loss?
================
Marissa Lee
10/23/2017

### Load microbial community data

### Load wood trait data

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png)

########################################## 

Wood traits as a predictor
--------------------------

*Hyp:* Stem-specific initial wood traits will predict variation in percent mass loss.

### time7

**less water content and C leads to more mass remaining after 7 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ waterperc + C, data = datasets[["time7"]])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.23820 -0.05588  0.01033  0.07087  0.19121 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.962933   0.493400   3.978  0.00022 ***
    ## waterperc   -0.003632   0.001795  -2.023  0.04831 *  
    ## C           -0.018543   0.009354  -1.982  0.05284 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09295 on 51 degrees of freedom
    ## Multiple R-squared:  0.1271, Adjusted R-squared:  0.0929 
    ## F-statistic: 3.714 on 2 and 51 DF,  p-value: 0.0312

### time 13

**larger size stems, less water content, more P and Mn, and less N leads to more mass remaining after 13 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ size + waterperc + P + Mn + N, data = datasets[["time13"]])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.277899 -0.092108  0.009127  0.079788  0.297114 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.4304321  0.1740140   8.220  1.2e-10 ***
    ## sizesmall   -0.1720630  0.0458801  -3.750 0.000484 ***
    ## waterperc   -0.0135897  0.0034603  -3.927 0.000280 ***
    ## P            0.0005436  0.0002667   2.038 0.047219 *  
    ## Mn           0.0006215  0.0002790   2.228 0.030735 *  
    ## N           -0.4242546  0.1421344  -2.985 0.004492 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1415 on 47 degrees of freedom
    ## Multiple R-squared:  0.363,  Adjusted R-squared:  0.2952 
    ## F-statistic: 5.357 on 5 and 47 DF,  p-value: 0.0005615

### time 25

**larger size stems, less water content, more P leads to more mass remaining after 25 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ size + waterperc + P, data = datasets[["time25"]])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.258325 -0.059630 -0.005782  0.068751  0.261561 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.1953820  0.1176559  10.160 1.50e-13 ***
    ## sizesmall   -0.0721459  0.0348728  -2.069   0.0440 *  
    ## waterperc   -0.0128901  0.0023561  -5.471 1.59e-06 ***
    ## P            0.0005140  0.0002134   2.408   0.0199 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.116 on 48 degrees of freedom
    ## Multiple R-squared:  0.4107, Adjusted R-squared:  0.3738 
    ## F-statistic: 11.15 on 3 and 48 DF,  p-value: 1.139e-05

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png)

### time 37

**larger size stems, less water content, more P leads to more mass remaining after 37 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ size + waterperc + P + Mn + Zn + N + 
    ##     C, data = datasets[["time37"]])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.33992 -0.08112  0.01198  0.08868  0.27249 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.3343324  0.7751083   0.431  0.66838    
    ## sizesmall   -0.1924388  0.0444506  -4.329 8.78e-05 ***
    ## waterperc   -0.0236176  0.0033363  -7.079 9.89e-09 ***
    ## P            0.0007591  0.0002639   2.877  0.00623 ** 
    ## Mn           0.0006056  0.0003035   1.995  0.05236 .  
    ## Zn          -0.0007677  0.0005331  -1.440  0.15708    
    ## N           -0.2447862  0.1443333  -1.696  0.09712 .  
    ## C            0.0255493  0.0146858   1.740  0.08906 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1324 on 43 degrees of freedom
    ## Multiple R-squared:  0.6233, Adjusted R-squared:  0.562 
    ## F-statistic: 10.16 on 7 and 43 DF,  p-value: 1.968e-07

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png)

*Hyp:* Variation in wood traits will lead to differences in decay model fit (r2), rate (k), and lagginess (alpha). Specifically, we expect samples with (a) high waterperc, (b) low density and C, (c) high P, K, Ca, Mn, Fe, Zn, and N, and (d) thicker bark (potential mech: limiting microbial colonization) to have better-fiting decay models (r2), faster decay rates (k), and less lagginess (alpha).

### r2

**greater water content and greater Zn and N leads to better-fitting decay models**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, the best model that greater water content and less C leads to better-fitting decay models

    ## 
    ## Call:
    ## lm(formula = ne.r2 ~ waterperc + barkthick + Ca + Zn + N + C, 
    ##     data = spdf.traits)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11753 -0.02941  0.00184  0.02233  0.12797 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.278e+00  6.457e-01   1.979   0.0590 .  
    ## waterperc    7.405e-03  1.509e-03   4.906 4.75e-05 ***
    ## barkthick   -3.098e-02  1.587e-02  -1.953   0.0622 .  
    ## Ca          -1.189e-05  6.783e-06  -1.753   0.0919 .  
    ## Zn           1.504e-03  5.664e-04   2.654   0.0136 *  
    ## N            2.092e-01  8.436e-02   2.480   0.0202 *  
    ## C           -1.698e-02  1.291e-02  -1.315   0.2006    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06306 on 25 degrees of freedom
    ## Multiple R-squared:  0.559,  Adjusted R-squared:  0.4532 
    ## F-statistic: 5.282 on 6 and 25 DF,  p-value: 0.001234

### k

**small size stems, greater water content, thinner bark, less Ca, more Zn, and more N lead to faster decay**

NOTE from Will: Density explains the same part of the variation in decay rates that initial water content does, only less well. (In other words, although, density gets dropped from the best model by the model selection procedure, if we remove initial water from consideration entirely, density is included in the model as the best predictor.)

So my current interpretation is that wood water rentention--related to fiber saturation point and partially captured by the density measurement--has a strong effect on long-term decomposition rates, possibly by maintaining fungal activity further into dry periods. There is also a very likely interaction between this water retention capacity with the fungal community (see results in Setting the Stage paper, Lee et al. in review).

    ## 
    ## Call:
    ## lm(formula = k ~ size + waterperc + barkthick + Ca + Zn + N + 
    ##     C, data = spdf.traits)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.092617 -0.034671  0.007271  0.029023  0.123640 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  6.021e-01  5.754e-01   1.046 0.305851    
    ## sizesmall    9.050e-02  2.405e-02   3.763 0.000958 ***
    ## waterperc    8.990e-03  1.378e-03   6.524 9.53e-07 ***
    ## barkthick   -3.862e-02  1.406e-02  -2.748 0.011198 *  
    ## Ca          -1.869e-05  5.920e-06  -3.157 0.004263 ** 
    ## Zn           1.811e-03  4.922e-04   3.680 0.001177 ** 
    ## N            2.334e-01  7.674e-02   3.041 0.005628 ** 
    ## C           -1.475e-02  1.135e-02  -1.300 0.206001    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05446 on 24 degrees of freedom
    ## Multiple R-squared:  0.7486, Adjusted R-squared:  0.6752 
    ## F-statistic: 10.21 on 7 and 24 DF,  p-value: 6.942e-06

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png)![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-2.png)

### t70

**small stem sizes, less water content, thicker bark, more Ca, less Zn, and less N lead to longer wood "70%"-lives**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, the best model indicated that large size stems, less water content, more Ca, and less Zn lead to longer wood "70%"-lives

    ## 
    ## Call:
    ## lm(formula = t70 ~ size + waterperc + barkthick + Ca + Zn + N, 
    ##     data = spdf.traits)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.55548 -0.15857 -0.00642  0.07054  0.58357 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.799e+00  4.125e-01   9.209 1.65e-09 ***
    ## sizesmall   -6.284e-01  1.314e-01  -4.783 6.54e-05 ***
    ## waterperc   -4.956e-02  7.653e-03  -6.476 8.79e-07 ***
    ## barkthick    1.752e-01  7.308e-02   2.397 0.024301 *  
    ## Ca           1.234e-04  3.203e-05   3.851 0.000725 ***
    ## Zn          -9.272e-03  2.698e-03  -3.436 0.002070 ** 
    ## N           -9.985e-01  4.258e-01  -2.345 0.027272 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3026 on 25 degrees of freedom
    ## Multiple R-squared:  0.735,  Adjusted R-squared:  0.6714 
    ## F-statistic: 11.56 on 6 and 25 DF,  p-value: 3.442e-06

### alpha--- don't interpret yet

########################################## 

Community as a predictor
------------------------

### First, filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 150 of 6128 OTUs"

*Hyp:* Stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.

### time 7

**none of the community components are significant predictors**

    ##             RMSE         R2      Avg.Bias  Max.Bias       Skill delta.RMSE
    ## Comp01 0.1096009 0.04039245 -0.0042025315 0.2086640   -7.426179   3.646601
    ## Comp02 0.1254558 0.04123537  0.0045006905 0.1767666  -40.754864  14.466035
    ## Comp03 0.1472351 0.02785978  0.0004267154 0.1727394  -93.867398  17.360158
    ## Comp04 0.1501175 0.04487863  0.0040991337 0.1600327 -101.532243   1.957664
    ## Comp05 0.1563567 0.04317337  0.0024552841 0.1777311 -118.632632   4.156223
    ##            p
    ## Comp01 0.716
    ## Comp02 0.955
    ## Comp03 0.969
    ## Comp04 0.745
    ## Comp05 0.932

    ##             RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1160401 0.03500705 0.004146464 0.2118409 -20.41985  9.7359781
    ## Comp02 0.1140335 0.03867892 0.005554219 0.2053950 -16.29114 -1.7292500
    ## Comp03 0.1141955 0.04829218 0.005694794 0.2145943 -16.62175  0.1420493
    ## Comp04 0.1137256 0.05486492 0.004972683 0.2101628 -15.66399 -0.4114765
    ## Comp05 0.1124598 0.06682365 0.006100160 0.2058893 -13.10370 -1.1129729
    ##            p
    ## Comp01 0.911
    ## Comp02 0.320
    ## Comp03 0.511
    ## Comp04 0.337
    ## Comp05 0.115

### time 13

**none of the community components are significant predictors**

    ##             RMSE         R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.2139657 0.21846772 -0.01496427 0.6260449  -74.86769 32.2375476
    ## Comp02 0.2149868 0.06885441 -0.01490427 0.5414328  -76.54081  0.4772558
    ## Comp03 0.2428530 0.06690494 -0.01653354 0.6290506 -125.27248 12.9617769
    ## Comp04 0.2693172 0.08533308 -0.01524579 0.6038045 -177.04454 10.8972377
    ## Comp05 0.2994919 0.07108880 -0.01525324 0.6220238 -242.60334 11.2041512
    ##            p
    ## Comp01 1.000
    ## Comp02 0.538
    ## Comp03 1.000
    ## Comp04 1.000
    ## Comp05 0.999

    ##             RMSE           R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.2163312 4.753407e-02 -0.02301907 0.6295907 -78.75562 33.6995201
    ## Comp02 0.2023872 5.526947e-06 -0.03447532 0.5853071 -56.45433 -6.4456583
    ## Comp03 0.2200750 8.812263e-03 -0.04106598 0.6142264 -84.99614  8.7395515
    ## Comp04 0.2244837 9.323948e-03 -0.04242462 0.6066593 -92.48239  2.0032860
    ## Comp05 0.2264918 9.961345e-03 -0.04426489 0.6070759 -95.94148  0.8945483
    ##            p
    ## Comp01 1.000
    ## Comp02 0.085
    ## Comp03 0.999
    ## Comp04 0.865
    ## Comp05 0.849

### time 25

**none of the community components are significant predictors**

    ##             RMSE         R2   Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1589171 0.07590492 0.01732259 0.3131693  -9.443048   4.615031
    ## Comp02 0.1723926 0.05607701 0.02297917 0.3430466 -28.790674   8.479611
    ## Comp03 0.1786746 0.05969270 0.02226170 0.3307871 -38.347911   3.643983
    ## Comp04 0.1813230 0.06504168 0.02194579 0.3286657 -42.479674   1.482266
    ## Comp05 0.1847188 0.06871557 0.02042740 0.3579705 -47.866326   1.872786
    ##            p
    ## Comp01 0.758
    ## Comp02 0.938
    ## Comp03 0.742
    ## Comp04 0.673
    ## Comp05 0.780

    ##             RMSE         R2      Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1568901 0.08501311 -1.251646e-03 0.2783429  -6.668957  3.2806647
    ## Comp02 0.1734482 0.03761080  1.856222e-03 0.3039221 -30.372702 10.5539613
    ## Comp03 0.1766744 0.02875939 -6.219912e-05 0.2900092 -35.267824  1.8600581
    ## Comp04 0.1753588 0.03066941  2.508066e-04 0.2874978 -33.260727 -0.7446700
    ## Comp05 0.1747009 0.03192784  1.039806e-03 0.2883656 -32.262694 -0.3751703
    ##            p
    ## Comp01 0.677
    ## Comp02 0.999
    ## Comp03 0.871
    ## Comp04 0.216
    ## Comp05 0.261

### time 37

**Comp01 is significant; if trim out OTUs that are not present in at least 20% of samples...Comp02 is significant**

    ##             RMSE        R2     Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1884529 0.1515172  0.033192872 0.2942286  5.557560 -2.8184996
    ## Comp02 0.1749664 0.2407742  0.008951682 0.2788971 18.591259 -7.1564113
    ## Comp03 0.1752001 0.2660234 -0.000556893 0.2185399 18.373655  0.1335602
    ## Comp04 0.1879489 0.2205367  0.003539710 0.2822099  6.062045  7.2766925
    ## Comp05 0.1967444 0.1897768  0.002389760 0.2656368 -2.935780  4.6797384
    ##            p
    ## Comp01 0.375
    ## Comp02 0.045
    ## Comp03 0.513
    ## Comp04 0.981
    ## Comp05 0.919

    ##             RMSE        R2     Avg.Bias  Max.Bias    Skill  delta.RMSE
    ## Comp01 0.1668533 0.2925742 -0.001261747 0.3122426 25.96599 -13.9569840
    ## Comp02 0.1615472 0.3377994 -0.007472631 0.2888952 30.59988  -3.1801339
    ## Comp03 0.1682930 0.3025985 -0.012605013 0.2854822 24.68285   4.1757949
    ## Comp04 0.1708756 0.2921546 -0.014679983 0.2863511 22.35357   1.5345464
    ## Comp05 0.1725708 0.2854071 -0.014202926 0.2902932 20.80528   0.9920909
    ##            p
    ## Comp01 0.044
    ## Comp02 0.217
    ## Comp03 0.954
    ## Comp04 0.845
    ## Comp05 0.908

Investigate the biology underlying time37-associated coefs for Comp02

By trophic mode (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By guild (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

*Hyp:* Average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### r2

**none of the community components are significant predictors**

### k

**none of the community components are significant predictors**

### t70

**none of the community components are significant predictors**

### alpha --- don't interpret yet

########################################## 

Community+traits as a predictor
-------------------------------

*Hyp:* After accounting for variation in decay due to wood traits (no models with barkthick or density), stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.

### time 7

**none of the community components are significant predictors**

    ##             RMSE          R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1115850 0.007432694 0.004235539 0.2120931  -52.60022  23.531460
    ## Comp02 0.1209657 0.004020937 0.013461663 0.1991320  -79.33627   8.406790
    ## Comp03 0.1291656 0.003573218 0.013284545 0.2371421 -104.47335   6.778618
    ## Comp04 0.1252574 0.001041763 0.011006738 0.2231303  -92.28723  -3.025653
    ## Comp05 0.1199016 0.001036001 0.011918383 0.2281983  -76.19499  -4.275843
    ##            p
    ## Comp01 0.992
    ## Comp02 0.754
    ## Comp03 0.877
    ## Comp04 0.233
    ## Comp05 0.341

### time 13

**Comp02 is significant, but if trim out OTUs that are not present in at least 20% of samples then no components are significant**

    ##             RMSE           R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1684942 7.969996e-03 -0.01227286 0.3008319  -59.80987  26.415927
    ## Comp02 0.1708079 8.029389e-05 -0.01134351 0.3428988  -64.22886   1.373151
    ## Comp03 0.1813353 1.626797e-02 -0.01398033 0.3387630  -85.09640   6.163259
    ## Comp04 0.1861811 2.670958e-03 -0.01869200 0.3402859  -95.12131   2.672318
    ## Comp05 0.2000881 9.511791e-03 -0.01941544 0.3603079 -125.35962   7.469616
    ##            p
    ## Comp01 0.999
    ## Comp02 0.594
    ## Comp03 0.976
    ## Comp04 0.846
    ## Comp05 1.000

    ##             RMSE         R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1958078 0.09587553 -0.01406066 0.3977193 -115.82077 46.9083971
    ## Comp02 0.1807162 0.02022406 -0.02532412 0.3640608  -83.83481 -7.7073210
    ## Comp03 0.1915530 0.02089810 -0.02731081 0.3899540 -106.54333  5.9965450
    ## Comp04 0.1942539 0.02389339 -0.02709579 0.3915635 -112.40896  1.4100120
    ## Comp05 0.1938709 0.02017169 -0.02696350 0.3920242 -111.57224 -0.1971538
    ##            p
    ## Comp01 1.000
    ## Comp02 0.046
    ## Comp03 0.972
    ## Comp04 0.833
    ## Comp05 0.391

Investigate the biology underlying time13-associated coefs for Comp02

By trophic mode ![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-32-1.png)![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-32-2.png)

By phylum

### time 25

**none of the community components are significant predictors**

    ##             RMSE        R2   Avg.Bias  Max.Bias     Skill delta.RMSE     p
    ## Comp01 0.1620698 0.1329636 0.01919070 0.3483804 -111.3412  45.375801 1.000
    ## Comp02 0.1795585 0.1117887 0.02724735 0.3954964 -159.4131  10.790831 0.920
    ## Comp03 0.1926807 0.1869273 0.02416387 0.4476233 -198.7147   7.308055 0.998
    ## Comp04 0.1905622 0.1539854 0.01920298 0.4451090 -192.1820  -1.099513 0.223
    ## Comp05 0.1957642 0.1553971 0.02117982 0.4710456 -208.3519   2.729837 0.955

### time 37

**none of the community components are significant predictors**

    ##             RMSE          R2     Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1519790 0.004659334  0.002664229 0.3654074  -56.17994  24.971974
    ## Comp02 0.1583613 0.009032758 -0.001947198 0.3948058  -69.57286   4.199477
    ## Comp03 0.1644513 0.007509707 -0.001848282 0.3944939  -82.86600   3.845652
    ## Comp04 0.1741376 0.004828966 -0.008279261 0.3849332 -105.04227   5.890064
    ## Comp05 0.1868349 0.013372894 -0.008065750 0.3920523 -136.03383   7.291528
    ##            p
    ## Comp01 0.994
    ## Comp02 0.844
    ## Comp03 0.884
    ## Comp04 0.992
    ## Comp05 0.996

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### r2

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was marginally significant

### k

**none of the community components are significant predictors**

### t70

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was a significant predictor

### alpha --- don't interpret yet

########################################## 

Diversity (and diversity of specific clades) as a predictor
-----------------------------------------------------------

**Note that the full community matrix was used for these analyses**

*Hyp:* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another.

### Richness

**No pattern**

### Shannon's H

**No pattern**

*Hyp:* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another.

### Saprotroph richness

**No pattern**

### Basidio richness

**No pattern**

*Hyp:* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers.

### Pathogen richness

**No pattern**

### Oomycete richness

**No pattern**

############################################## 

Diversity plus traits as a predictor
------------------------------------

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### Richness

**No pattern**

### Shannon's H

**No pattern**

### Saprotroph richness

**No pattern**

### Basidio richness

**No pattern**

### Pathogen richness

**No pattern**

### Oomycete richness

**No pattern**

############################################## 

Relationship between wood traits and community
----------------------------------------------

*Hyp:* Average initial microbial communitiy compositions will covary with initial wood traits

############################################## 

Extra pieces
------------

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
