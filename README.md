Does chemistry or community better predict mass loss?
================
Marissa Lee
10/23/2017

### Load microbial community data

### Load wood trait data

########################################## 

Wood traits as a preditor
-------------------------

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png)![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-2.png)

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

*Hyp:* Average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### First, filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 150 of 6128 OTUs"

### r2

**none of the community components are significant predictors**

    ##             RMSE           R2    Avg.Bias  Max.Bias     Skill  delta.RMSE
    ## Comp01 0.1013107 1.494708e-05 0.002024841 0.2609272 -40.32548 18.45905577
    ## Comp02 0.1104027 4.293834e-03 0.003308000 0.2568255 -66.64251  8.97443525
    ## Comp03 0.1114560 1.027904e-02 0.003371373 0.2497766 -69.83735  0.95404167
    ## Comp04 0.1119467 1.116217e-02 0.002666967 0.2514263 -71.33590  0.44020287
    ## Comp05 0.1120426 1.133987e-02 0.002626418 0.2509619 -71.62965  0.08568655
    ##            p
    ## Comp01 0.898
    ## Comp02 0.999
    ## Comp03 0.710
    ## Comp04 0.802
    ## Comp05 0.720

### k

**none of the community components are significant predictors**

    ##             RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1036580 0.01276672 0.000282797 0.2658329 -22.99604 10.9035799
    ## Comp02 0.1044388 0.01904155 0.011480363 0.2688760 -24.85590  0.7532281
    ## Comp03 0.1059335 0.01805338 0.010666194 0.2689647 -28.45534  1.4311977
    ## Comp04 0.1087908 0.02741427 0.013442320 0.2667951 -35.47840  2.6972785
    ## Comp05 0.1125616 0.01442363 0.011293413 0.2704606 -45.03280  3.4661026
    ##            p
    ## Comp01 0.884
    ## Comp02 0.568
    ## Comp03 0.636
    ## Comp04 0.875
    ## Comp05 0.964

### t70

**none of the community components are significant predictors**

    ##             RMSE         R2      Avg.Bias Max.Bias      Skill  delta.RMSE
    ## Comp01 0.5638322 0.02747152  0.0291722604 1.503069 -19.039193  9.10508360
    ## Comp02 0.5363187 0.06687807 -0.0040312802 1.258095  -7.705054 -4.87973434
    ## Comp03 0.5367893 0.08213251 -0.0004053627 1.153682  -7.894179  0.08775921
    ## Comp04 0.5583885 0.08221500 -0.0056512738 1.349964 -16.751702  4.02377323
    ## Comp05 0.5693485 0.06691379  0.0058129906 1.431625 -21.379860  1.96278895
    ##            p
    ## Comp01 0.842
    ## Comp02 0.215
    ## Comp03 0.505
    ## Comp04 0.837
    ## Comp05 0.848

### alpha --- don't interpret yet

########################################## 

Community+traits as a predictor
-------------------------------

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### r2

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was marginally significant

    ##              RMSE           R2     Avg.Bias  Max.Bias      Skill
    ## Comp01 0.07083530 0.0004108819 -0.006446032 0.1400555  -61.51372
    ## Comp02 0.09185757 0.0006934417  0.002667991 0.1556453 -171.60633
    ## Comp03 0.09573347 0.0023892283  0.004948795 0.1592020 -195.01055
    ## Comp04 0.09706111 0.0019152166  0.003816465 0.1568796 -203.24973
    ## Comp05 0.09665378 0.0024386184  0.004673572 0.1604111 -200.70985
    ##        delta.RMSE     p
    ## Comp01  27.088047 0.997
    ## Comp02  29.677682 0.999
    ## Comp03   4.219463 0.933
    ## Comp04   1.386804 0.804
    ## Comp05  -0.419656 0.304

### k

**none of the community components are significant predictors**

    ##              RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.05403597 0.01169109 0.000422351 0.1254135 -31.24173  14.560784
    ## Comp02 0.05065152 0.04533652 0.002916937 0.1244362 -15.31638  -6.263328
    ## Comp03 0.04993864 0.05990280 0.004727269 0.1237288 -12.09324  -1.407427
    ## Comp04 0.05147371 0.04410936 0.005052921 0.1270141 -19.09043   3.073904
    ## Comp05 0.05142763 0.04472957 0.005796166 0.1291158 -18.87733  -0.089509
    ##            p
    ## Comp01 0.877
    ## Comp02 0.166
    ## Comp03 0.239
    ## Comp04 0.920
    ## Comp05 0.489

### t70

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was a significant predictor

    ##             RMSE        R2   Avg.Bias  Max.Bias     Skill delta.RMSE     p
    ## Comp01 0.2564979 0.1696737 0.02216625 0.4864828  8.061448 -4.1154070 0.289
    ## Comp02 0.2615463 0.1965626 0.01788789 0.4573872  4.406754  1.9682055 0.576
    ## Comp03 0.2745169 0.1851611 0.01979289 0.4458862 -5.309638  4.9591861 0.886
    ## Comp04 0.2776540 0.1774745 0.02286057 0.4449253 -7.730289  1.1427719 0.743
    ## Comp05 0.2795150 0.1716365 0.02355805 0.4637620 -9.179284  0.6702646 0.664

### alpha --- don't interpret yet

**Note: no longer warrented based on analyses with waterperc represented as g water / g wet weight...**

Investigate the biology underlying t70-associated coefs for Comp05

By trophic mode (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By guild (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By phylum

########################################## 

Diversity (and diversity of specific clades) as a predictor
-----------------------------------------------------------

**Note that the full community matrix was used for these analyses**

*Hyp:* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another.

### Richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8943
    ## mean       1 0.000245 0.0002451  0.0305 0.8625
    ## Residuals 30 0.240983 0.0080328

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0055 0.03285 *
    ## mean       1 0.000163 0.000163  0.0198 0.88904  
    ## Residuals 30 0.246927 0.008231                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9953 0.02041 *
    ## mean       1 0.0226 0.02258  0.0925 0.76317  
    ## Residuals 30 7.3263 0.24421                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0205 0.8871
    ## mean       1 0.15961 0.159615  2.3308 0.1373
    ## Residuals 30 2.05446 0.068482

### Shannon's H

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8939
    ## mean       1 0.002429 0.0024291  0.3052 0.5848
    ## Residuals 30 0.238799 0.0079600

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1360 0.03081 *
    ## mean       1 0.006437 0.006437  0.8024 0.37750  
    ## Residuals 30 0.240653 0.008022                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2830 0.01784 *
    ## mean       1 0.3580 0.35804  1.5365 0.22475  
    ## Residuals 30 6.9908 0.23303                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.0014 0.001404  0.0196 0.8896
    ## mean       1 0.0657 0.065698  0.9174 0.3458
    ## Residuals 30 2.1484 0.071613

*Hyp:* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another.

### Saprotroph richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000006 0.0000062  0.0008 0.9780
    ## Residuals 30 0.241222 0.0080407

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.04120 0.041200  5.0515 0.03211 *
    ## mean       1 0.00241 0.002410  0.2955 0.59075  
    ## Residuals 30 0.24468 0.008156                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9893 0.02046 *
    ## mean       1 0.0152 0.01524  0.0624 0.80452  
    ## Residuals 30 7.3336 0.24445                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)  
    ## size       1 0.00140 0.001404  0.0211 0.8854  
    ## mean       1 0.22078 0.220782  3.3229 0.0783 .
    ## Residuals 30 1.99330 0.066443                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Basidio richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000072 0.0000717  0.0089 0.9254
    ## Residuals 30 0.241157 0.0080386

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1545 0.03053 *
    ## mean       1 0.007301 0.007301  0.9134 0.34686  
    ## Residuals 30 0.239789 0.007993                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.1675 0.01883 *
    ## mean       1 0.2271 0.22708  0.9565 0.33588  
    ## Residuals 30 7.1218 0.23739                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0200 0.8886
    ## mean       1 0.10539 0.105391  1.4994 0.2303
    ## Residuals 30 2.10869 0.070290

*Hyp:* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers.

### Pathogen richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8940
    ## mean       1 0.001760 0.0017599  0.2205 0.6421
    ## Residuals 30 0.239469 0.0079823

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0603 0.03197 *
    ## mean       1 0.002835 0.002835  0.3482 0.55956  
    ## Residuals 30 0.244255 0.008142                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2855 0.01782 *
    ## mean       1 0.3608 0.36079  1.5489 0.22294  
    ## Residuals 30 6.9881 0.23294                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.00140 0.001404  0.0215 0.88430  
    ## mean       1 0.25844 0.258439  3.9645 0.05564 .
    ## Residuals 30 1.95564 0.065188                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Oomycete richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0184 0.8931
    ## mean       1 0.005902 0.0059025  0.7525 0.3926
    ## Residuals 30 0.235326 0.0078442

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0038 0.03288 *
    ## mean       1 0.000077 0.000077  0.0093 0.92374  
    ## Residuals 30 0.247013 0.008234                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9831 0.02052 *
    ## mean       1 0.0076 0.00763  0.0312 0.86099  
    ## Residuals 30 7.3412 0.24471                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0192 0.8906
    ## mean       1 0.02338 0.023381  0.3202 0.5757
    ## Residuals 30 2.19070 0.073023

############################################## 

Diversity plus traits as a predictor
------------------------------------

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### Richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8943
    ## mean       1 0.000245 0.0002451  0.0305 0.8625
    ## Residuals 30 0.240983 0.0080328               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0055 0.03285 *
    ## mean       1 0.000163 0.000163  0.0198 0.88904  
    ## Residuals 30 0.246927 0.008231                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9953 0.02041 *
    ## mean       1 0.0226 0.02258  0.0925 0.76317  
    ## Residuals 30 7.3263 0.24421                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0205 0.8871
    ## mean       1 0.15961 0.159615  2.3308 0.1373
    ## Residuals 30 2.05446 0.068482

### Shannon's H

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8939
    ## mean       1 0.002429 0.0024291  0.3052 0.5848
    ## Residuals 30 0.238799 0.0079600               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1360 0.03081 *
    ## mean       1 0.006437 0.006437  0.8024 0.37750  
    ## Residuals 30 0.240653 0.008022                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2830 0.01784 *
    ## mean       1 0.3580 0.35804  1.5365 0.22475  
    ## Residuals 30 6.9908 0.23303                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.0014 0.001404  0.0196 0.8896
    ## mean       1 0.0657 0.065698  0.9174 0.3458
    ## Residuals 30 2.1484 0.071613

### Saprotroph richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000006 0.0000062  0.0008 0.9780
    ## Residuals 30 0.241222 0.0080407               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.04120 0.041200  5.0515 0.03211 *
    ## mean       1 0.00241 0.002410  0.2955 0.59075  
    ## Residuals 30 0.24468 0.008156                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9893 0.02046 *
    ## mean       1 0.0152 0.01524  0.0624 0.80452  
    ## Residuals 30 7.3336 0.24445                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)  
    ## size       1 0.00140 0.001404  0.0211 0.8854  
    ## mean       1 0.22078 0.220782  3.3229 0.0783 .
    ## Residuals 30 1.99330 0.066443                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Basidio richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000072 0.0000717  0.0089 0.9254
    ## Residuals 30 0.241157 0.0080386               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1545 0.03053 *
    ## mean       1 0.007301 0.007301  0.9134 0.34686  
    ## Residuals 30 0.239789 0.007993                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.1675 0.01883 *
    ## mean       1 0.2271 0.22708  0.9565 0.33588  
    ## Residuals 30 7.1218 0.23739                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0200 0.8886
    ## mean       1 0.10539 0.105391  1.4994 0.2303
    ## Residuals 30 2.10869 0.070290

### Pathogen richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8940
    ## mean       1 0.001760 0.0017599  0.2205 0.6421
    ## Residuals 30 0.239469 0.0079823               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0603 0.03197 *
    ## mean       1 0.002835 0.002835  0.3482 0.55956  
    ## Residuals 30 0.244255 0.008142                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2855 0.01782 *
    ## mean       1 0.3608 0.36079  1.5489 0.22294  
    ## Residuals 30 6.9881 0.23294                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.00140 0.001404  0.0215 0.88430  
    ## mean       1 0.25844 0.258439  3.9645 0.05564 .
    ## Residuals 30 1.95564 0.065188                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Oomycete richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0184 0.8931
    ## mean       1 0.005902 0.0059025  0.7525 0.3926
    ## Residuals 30 0.235326 0.0078442               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0038 0.03288 *
    ## mean       1 0.000077 0.000077  0.0093 0.92374  
    ## Residuals 30 0.247013 0.008234                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9831 0.02052 *
    ## mean       1 0.0076 0.00763  0.0312 0.86099  
    ## Residuals 30 7.3412 0.24471                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0192 0.8906
    ## mean       1 0.02338 0.023381  0.3202 0.5757
    ## Residuals 30 2.19070 0.073023

############################################## 

Extra pieces
============

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
