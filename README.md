Does chemistry or community better predict mass loss?
================
Marissa Lee
10/23/2017

``` r
#chunk options
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

#libraries
devtools::install_github("cornwell-lab-unsw/litterfitter")
library(dplyr)
library(ggplot2)
library(readr)
library(vegan)
library(knitr)
library(litterfitter)
library(magrittr)
library(tidyr)
library(gridExtra)
library(rioja)

#fxns
source("code/load_fxns.R")
source("code/curve_fitting_fxns.R")
source("code/distance_fxns.R")
source("code/otuIDs_fxns.R")
```

Load data
---------

### Microbial community data

### Wood trait data

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

    ## 
    ## Call:
    ## lm(formula = alpha ~ density + barkthick + P + K + Ca + Fe + 
    ##     Zn + N + C, data = spdf.traits)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.28961 -0.10507  0.01239  0.05847  0.37028 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  6.978e+00  2.515e+00   2.774   0.0111 *
    ## density     -1.550e+00  5.595e-01  -2.770   0.0112 *
    ## barkthick   -8.276e-02  5.373e-02  -1.540   0.1377  
    ## P           -7.629e-04  5.649e-04  -1.350   0.1906  
    ## K            6.411e-05  4.567e-05   1.404   0.1743  
    ## Ca          -4.290e-05  2.474e-05  -1.735   0.0968 .
    ## Fe          -3.027e-05  2.051e-05  -1.475   0.1543  
    ## Zn           4.064e-03  1.897e-03   2.143   0.0435 *
    ## N            3.470e-01  2.794e-01   1.242   0.2274  
    ## C           -9.831e-02  4.674e-02  -2.103   0.0471 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1989 on 22 degrees of freedom
    ## Multiple R-squared:  0.6067, Adjusted R-squared:  0.4458 
    ## F-statistic: 3.771 on 9 and 22 DF,  p-value: 0.005306

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
    ## Comp01 0.890
    ## Comp02 1.000
    ## Comp03 0.715
    ## Comp04 0.789
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
    ## Comp01 0.902
    ## Comp02 0.519
    ## Comp03 0.631
    ## Comp04 0.856
    ## Comp05 0.965

### t70

**none of the community components are significant predictors**

    ##             RMSE         R2      Avg.Bias Max.Bias      Skill  delta.RMSE
    ## Comp01 0.5638322 0.02747152  0.0291722604 1.503069 -19.039193  9.10508360
    ## Comp02 0.5363187 0.06687807 -0.0040312802 1.258095  -7.705054 -4.87973434
    ## Comp03 0.5367893 0.08213251 -0.0004053627 1.153682  -7.894179  0.08775921
    ## Comp04 0.5583885 0.08221500 -0.0056512738 1.349964 -16.751702  4.02377323
    ## Comp05 0.5693485 0.06691379  0.0058129906 1.431625 -21.379860  1.96278895
    ##            p
    ## Comp01 0.873
    ## Comp02 0.227
    ## Comp03 0.501
    ## Comp04 0.838
    ## Comp05 0.831

### alpha --- don't interpret yet

    ##             RMSE           R2   Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.3217580 5.393130e-04 0.01725191 0.7545470 -54.20726  24.180214
    ## Comp02 0.3108008 2.744262e-03 0.02498895 0.7692934 -43.88325  -3.405429
    ## Comp03 0.3317582 4.810973e-05 0.04473939 0.7738432 -63.94162   6.743021
    ## Comp04 0.3519272 1.367993e-02 0.03667556 0.8064912 -84.48099   6.079436
    ## Comp05 0.3451306 9.703285e-03 0.03466600 0.8074728 -77.42424  -1.931244
    ##            p
    ## Comp01 0.963
    ## Comp02 0.282
    ## Comp03 0.927
    ## Comp04 0.994
    ## Comp05 0.175

########################################## 

Community+traits as a predictor
-------------------------------

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### r2

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was marginally significant

``` r
rand.t.test(fit.tr.r2.cv)
```

    ##              RMSE           R2     Avg.Bias  Max.Bias      Skill
    ## Comp01 0.07083530 0.0004108819 -0.006446032 0.1400555  -61.51372
    ## Comp02 0.09185757 0.0006934417  0.002667991 0.1556453 -171.60633
    ## Comp03 0.09573347 0.0023892283  0.004948795 0.1592020 -195.01055
    ## Comp04 0.09706111 0.0019152166  0.003816465 0.1568796 -203.24973
    ## Comp05 0.09665378 0.0024386184  0.004673572 0.1604111 -200.70985
    ##        delta.RMSE     p
    ## Comp01  27.088047 0.994
    ## Comp02  29.677682 1.000
    ## Comp03   4.219463 0.911
    ## Comp04   1.386804 0.801
    ## Comp05  -0.419656 0.276

### k

**none of the community components are significant predictors**

``` r
rand.t.test(fit.tr.k.cv)
```

    ##              RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.05403597 0.01169109 0.000422351 0.1254135 -31.24173  14.560784
    ## Comp02 0.05065152 0.04533652 0.002916937 0.1244362 -15.31638  -6.263328
    ## Comp03 0.04993864 0.05990280 0.004727269 0.1237288 -12.09324  -1.407427
    ## Comp04 0.05147371 0.04410936 0.005052921 0.1270141 -19.09043   3.073904
    ## Comp05 0.05142763 0.04472957 0.005796166 0.1291158 -18.87733  -0.089509
    ##            p
    ## Comp01 0.916
    ## Comp02 0.193
    ## Comp03 0.222
    ## Comp04 0.941
    ## Comp05 0.480

### t70

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was a significant predictor

``` r
rand.t.test(fit.tr.t70.cv) 
```

    ##             RMSE        R2   Avg.Bias  Max.Bias     Skill delta.RMSE     p
    ## Comp01 0.2564979 0.1696737 0.02216625 0.4864828  8.061448 -4.1154070 0.281
    ## Comp02 0.2615463 0.1965626 0.01788789 0.4573872  4.406754  1.9682055 0.549
    ## Comp03 0.2745169 0.1851611 0.01979289 0.4458862 -5.309638  4.9591861 0.888
    ## Comp04 0.2776540 0.1774745 0.02286057 0.4449253 -7.730289  1.1427719 0.719
    ## Comp05 0.2795150 0.1716365 0.02355805 0.4637620 -9.179284  0.6702646 0.644

### alpha --- don't interpret yet

``` r
rand.t.test(fit.tr.alpha.cv)
```

    ##             RMSE           R2      Avg.Bias  Max.Bias      Skill
    ## Comp01 0.1847965 2.827379e-02 -0.0201047453 0.3426950  -25.54006
    ## Comp02 0.2460519 1.842392e-04 -0.0053048937 0.4084135 -122.56047
    ## Comp03 0.2612126 4.069172e-06  0.0020589171 0.4224961 -150.83203
    ## Comp04 0.2658597 1.990507e-04  0.0062806609 0.4345492 -159.83617
    ## Comp05 0.2667991 4.823803e-04 -0.0002445828 0.4396092 -161.67561
    ##        delta.RMSE     p
    ## Comp01 12.0446595 0.942
    ## Comp02 33.1474500 0.999
    ## Comp03  6.1616065 0.862
    ## Comp04  1.7790289 0.644
    ## Comp05  0.3533382 0.581

**Note: no longer warrented based on analyses with waterperc represented as g water / g wet weight...**

Investigate the biology underlying t70-associated coefs for Comp05

By trophic mode (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By guild (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By phylum

########################################## 

1.  Diversity (and diversity of specific clades) as a predictor

*Hyp:* Greater microbial diversity (richness, Shannon diversity, phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another.

*Hyp:* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another. **No significant relationships**

``` r
# summarize the presence of ... in each sample
sapro.df<-Calc_richOTUtype(colNam="Trophic.Mode", grepTerm="Sapro", taxAndFunguild, comm.otu)  
basidio.df<-Calc_richOTUtype(colNam="phylum", grepTerm="Basid", taxAndFunguild, comm.otu)

# create a merged df wtih spdf
saprorich.spdf<-Create_rich_spdf_DF(otutype.df=sapro.df, spdf)
basidrich.spdf<-Create_rich_spdf_DF(otutype.df=basidio.df, spdf)

# fit models
mod.sapro.r2<-lm(ne.r2~size+mean, data=saprorich.spdf)
mod.sapro.k<-lm(k~size+mean, data=saprorich.spdf)
mod.sapro.t70<-lm(t70~size+mean, data=saprorich.spdf)
mod.sapro.alpha<-lm(alpha~size+mean, data=saprorich.spdf)
mod.basid.r2<-lm(ne.r2~size+mean, data=basidrich.spdf)
mod.basid.k<-lm(k~size+mean, data=basidrich.spdf)
mod.basid.t70<-lm(t70~size+mean, data=basidrich.spdf)
mod.basid.alpha<-lm(alpha~size+mean, data=basidrich.spdf)

# evaluate models
summary(mod.sapro.r2)
```

    ## 
    ## Call:
    ## lm(formula = ne.r2 ~ size + mean, data = saprorich.spdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.20146 -0.07071 -0.00394  0.07046  0.13458 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.553e-01  6.277e-02  12.032 5.22e-13 ***
    ## sizesmall   -3.727e-03  3.263e-02  -0.114    0.910    
    ## mean        -6.339e-05  3.629e-04  -0.175    0.863    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08963 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.001612,   Adjusted R-squared:  -0.06495 
    ## F-statistic: 0.02422 on 2 and 30 DF,  p-value: 0.9761

``` r
summary(mod.sapro.k)
```

    ## 
    ## Call:
    ## lm(formula = k ~ size + mean, data = saprorich.spdf)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.144459 -0.066842 -0.002603  0.054204  0.194414 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.2450387  0.0635404   3.856 0.000566 ***
    ## sizesmall    0.0739558  0.0330253   2.239 0.032698 *  
    ## mean        -0.0000517  0.0003674  -0.141 0.889036    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09072 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.1435, Adjusted R-squared:  0.08637 
    ## F-statistic: 2.513 on 2 and 30 DF,  p-value: 0.09797

``` r
summary(mod.sapro.t70)
```

    ## 
    ## Call:
    ## lm(formula = t70 ~ size + mean, data = saprorich.spdf)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8442 -0.3241 -0.1319  0.2610  1.2428 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.5978639  0.3461056   4.617 6.86e-05 ***
    ## sizesmall   -0.4437990  0.1798894  -2.467   0.0196 *  
    ## mean         0.0006085  0.0020012   0.304   0.7632    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4942 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.1687, Adjusted R-squared:  0.1133 
    ## F-statistic: 3.044 on 2 and 30 DF,  p-value: 0.06258

``` r
summary(mod.sapro.alpha)
```

    ## 
    ## Call:
    ## lm(formula = alpha ~ size + mean, data = saprorich.spdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40959 -0.17508 -0.07703  0.23693  0.59348 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  0.612096   0.183280   3.340  0.00225 **
    ## sizesmall   -0.002206   0.095260  -0.023  0.98167   
    ## mean         0.001618   0.001060   1.527  0.13732   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2617 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.07268,    Adjusted R-squared:  0.01086 
    ## F-statistic: 1.176 on 2 and 30 DF,  p-value: 0.3224

``` r
summary(mod.basid.r2)
```

    ## 
    ## Call:
    ## lm(formula = ne.r2 ~ size + mean, data = basidrich.spdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.20146 -0.07071 -0.00394  0.07046  0.13458 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.553e-01  6.277e-02  12.032 5.22e-13 ***
    ## sizesmall   -3.727e-03  3.263e-02  -0.114    0.910    
    ## mean        -6.339e-05  3.629e-04  -0.175    0.863    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08963 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.001612,   Adjusted R-squared:  -0.06495 
    ## F-statistic: 0.02422 on 2 and 30 DF,  p-value: 0.9761

``` r
summary(mod.basid.k)
```

    ## 
    ## Call:
    ## lm(formula = k ~ size + mean, data = basidrich.spdf)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.144459 -0.066842 -0.002603  0.054204  0.194414 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.2450387  0.0635404   3.856 0.000566 ***
    ## sizesmall    0.0739558  0.0330253   2.239 0.032698 *  
    ## mean        -0.0000517  0.0003674  -0.141 0.889036    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09072 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.1435, Adjusted R-squared:  0.08637 
    ## F-statistic: 2.513 on 2 and 30 DF,  p-value: 0.09797

``` r
summary(mod.basid.t70)
```

    ## 
    ## Call:
    ## lm(formula = t70 ~ size + mean, data = basidrich.spdf)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8442 -0.3241 -0.1319  0.2610  1.2428 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.5978639  0.3461056   4.617 6.86e-05 ***
    ## sizesmall   -0.4437990  0.1798894  -2.467   0.0196 *  
    ## mean         0.0006085  0.0020012   0.304   0.7632    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4942 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.1687, Adjusted R-squared:  0.1133 
    ## F-statistic: 3.044 on 2 and 30 DF,  p-value: 0.06258

``` r
summary(mod.basid.alpha)
```

    ## 
    ## Call:
    ## lm(formula = alpha ~ size + mean, data = basidrich.spdf)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40959 -0.17508 -0.07703  0.23693  0.59348 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  0.612096   0.183280   3.340  0.00225 **
    ## sizesmall   -0.002206   0.095260  -0.023  0.98167   
    ## mean         0.001618   0.001060   1.527  0.13732   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2617 on 30 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.07268,    Adjusted R-squared:  0.01086 
    ## F-statistic: 1.176 on 2 and 30 DF,  p-value: 0.3224

``` r
# create plots
sapList<-Plot_richOTUtype(rich.spdf=saprorich.spdf, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Saprotroph")
basidList<-Plot_richOTUtype(rich.spdf = basidrich.spdf, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Basidio")

#plot
grid.arrange(sapList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             sapList[['k']] + guides(color=FALSE, shape=FALSE), 
             sapList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             basidList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             basidList[['k']] + guides(color=FALSE, shape=FALSE), 
             basidList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             ncol=3)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-25-1.png)

*Hyp:* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers. **Maybe there's some indicaiton that oomycete presence increases the likelihood of slower (k) and more laggy (alpha) decay**

``` r
# summarize the presence of ... in each sample
path.df<-Calc_richOTUtype(colNam="Trophic.Mode", grepTerm="Patho", taxAndFunguild, comm.otu)
oomy.df<-Calc_richOTUtype(colNam="kingdom", grepTerm="Protist", taxAndFunguild, comm.otu)

# create a merged df wtih spdf
pathrich.spdf<-Create_rich_spdf_DF(otutype.df=path.df, spdf)
oomyrich.spdf<-Create_rich_spdf_DF(otutype.df=oomy.df, spdf)

# # fit models
# mod.sapro.r2<-lm(ne.r2~size+mean, data=saprorich.spdf)
# mod.sapro.k<-lm(k~size+mean, data=saprorich.spdf)
# mod.sapro.t70<-lm(t70~size+mean, data=saprorich.spdf)
# mod.sapro.alpha<-lm(alpha~size+mean, data=saprorich.spdf)
# mod.basid.r2<-lm(ne.r2~size+mean, data=basidrich.spdf)
# mod.basid.k<-lm(k~size+mean, data=basidrich.spdf)
# mod.basid.t70<-lm(t70~size+mean, data=basidrich.spdf)
# mod.basid.alpha<-lm(alpha~size+mean, data=basidrich.spdf)
# 
# # evaluate models
# summary(mod.sapro.r2)
# summary(mod.sapro.k)
# summary(mod.sapro.t70)
# summary(mod.sapro.alpha)
# summary(mod.basid.r2)
# summary(mod.basid.k)
# summary(mod.basid.t70)
# summary(mod.basid.alpha)

# create plots
pathList<-Plot_richOTUtype(rich.spdf=pathrich.spdf, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Pathogen")
ooList<-Plot_richOTUtype(rich.spdf=oomyrich.spdf, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Oomycete")

# plot
grid.arrange(pathList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             pathList[['k']] + guides(color=FALSE, shape=FALSE), 
             pathList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             ooList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             ooList[['k']] + guides(color=FALSE, shape=FALSE), 
             ooList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             ncol=3)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-26-1.png)

1.  Diversity plus traits as a predictor

############################################## 

Extra pieces
============

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
