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

LOAD DATA
---------

### MICROBIAL COMMUNITY DATA

``` r
#stem sample meta data
#stemSamples<-load_stemSamples() #uncomment if the data changes
#write_csv(stemSamples, "derived_data/stemSamples.csv")
stemSamples<-read_csv("derived_data/stemSamples.csv")

#OTU table
#fung.otu<-load_matotu() #uncomment if the data changes
#comm.otu<-add_oomycetes(fung.otu) #add the oomycetes #uncomment if the data changes
#write.csv(comm.otu, "derived_data/comm_otu.csv")
comm.otu<-read.csv("derived_data/comm_otu.csv", row.names=1)

#create sequence sample meta data table
#seqSamples<-load_seqSamples(comm.otu, stemSamples) #uncomment if the data changes
#write_csv(seqSamples, "derived_data/seqSamples.csv")
seqSamples<-read_csv("derived_data/seqSamples.csv")

#taxon lookup info
#taxAndFunguild<-load_TaxAndFunguild(comm.otu) #uncomment if the data changes
#write_csv(taxAndFunguild, "derived_data/taxaAndFunguild.csv")
taxAndFunguild<-read_csv("derived_data/taxaAndFunguild.csv")

#plot_sampleEffortCurves(comm.otu)
```

### LOAD WOOD TRAIT DATA

``` r
#traits.mean<-mergeTraitData() #uncomment if the data changes
#write_csv(traits.mean, "derived_data/traits_mean.csv") 
traits.mean<-read_csv("derived_data/traits_mean.csv")
#missing data
#traits.long<-as.data.frame(gather(traits.mean, key=trait, value=value, -(1:3)))
#filter(traits.long, is.na(value))


### LOAD MASS LOSS DATA and CALCULATE % MASS REMAINING AT EACH TIMEPOINT
#initial_mass <- read_in_initial_mass() #uncomment if the data changes
#harvest_mass<-LoadHarvestFiles()
#mass.data<-bind_rows(initial_mass, harvest_mass)
#missing data
#mass.data %>% filter(is.na(totalSampleDryMass))
#plotting_df<-Calc_massRemaining(mass.data)
#matching failures
# plotting_df %>%
#   filter(is.na(pmr)) %>%
#   select(unique, species, size, time, totalSampleDryMass, notes) %>%
#   spread(key=time, value=totalSampleDryMass)
#remove NAs
#plotting_df %>% filter(!is.na(pmr)) -> plotting_df
#write_csv(plotting_df,"derived_data/plotting_df.csv")
plotting_df<-read_csv("derived_data/plotting_df.csv")

### CALCULATE DECAY TRAJECTORY FITS
#spdf <- fit_all_curves(plotting_df) #this recalculates all the curve fits, uncomment if the data changes
#indx<-select(stemSamples, code, species, size)
#spdf<-left_join(spdf, indx) #add code
#write_csv(spdf,"derived_data/mass_loss_parameters.csv")
spdf <- read_csv("derived_data/mass_loss_parameters.csv")
# ggplot(spdf,aes(x=t70,y=w.t70,col=size))+
#   geom_point()+
#   labs(x="Time to 30% mass loss (negative exponential)", 
#        y="Time to 30% mass loss (Weibull)")+
#   geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()
```

Wood traits as a preditor
-------------------------

*Hyp:* Variation in wood traits will lead to differences in decay model fit (r2), rate (k), and lagginess (alpha). Specifically, we expect samples with (a) high waterperc, (b) low density and C, (c) high P, K, Ca, Mn, Fe, Zn, and N, and (d) thicker bark (potential mech: limiting microbial colonization) to have better-fiting decay models (r2), faster decay rates (k), and less lagginess (alpha).

*Result:* - r2... greater water content and greater Zn and N leads to better-fitting decay models. Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, the best model that greater water content and less C leads to better-fitting decay models

``` r
summary(mod.select.r) # waterperc, Zn, N
```

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

``` r
#ggplot(spdf.traits, aes(x=waterperc, y=ne.r2, color=species, shape=size)) + geom_point()
```

-   k... small size stems, greater water content, thinner bark, less Ca, more Zn, and more N lead to faster decay

NOTE from Will: Density explains the same part of the variation in decay rates that initial water content does, only less well. (In other words, although, density gets dropped from the best model by the model selection procedure, if we remove initial water from consideration entirely, density is included in the model as the best predictor.)

So my current interpretation is that wood water rentention--related to fiber saturation point and partially captured by the density measurement--has a strong effect on long-term decomposition rates, possibly by maintaining fungal activity further into dry periods. There is also a very likely interaction between this water retention capacity with the fungal community (see results in Setting the Stage paper, Lee et al. in review).

``` r
summary(mod.select.k) # size, waterperc, barkthick, Ca, Zn, N
```

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

``` r
ggplot(spdf.traits, aes(x=waterperc, y=k)) + geom_point(aes(color=species)) + facet_grid(~size) +
  labs(y="k (year^-1)",x="Initial water content (% wet weight)")+geom_smooth(method="lm",se=FALSE)+theme_bw()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png)

``` r
ggplot(spdf.traits, aes(x=density, y=k)) + geom_point(aes(color=species)) + facet_grid(~size) +
  labs(y="k (year^-1)",x="Initial density (g/cm^3)")+geom_smooth(method="lm",se=FALSE)+theme_bw()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-2.png)

``` r
#ggplot(spdf.traits, aes(x=waterperc, y=density, color=species, size=k)) + geom_point() + facet_grid(~size)
```

-   t70... small stem sizes, less water content, thicker bark, more Ca, less Zn, and less N lead to longer wood "70%"-lives. Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, the best model indicated that large size stems, less water content, more Ca, and less Zn lead to longer wood "70%"-lives

``` r
summary(mod.select.t70) # size, waterperc, barkthick, Ca, Zn
```

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

``` r
#ggplot(spdf.traits, aes(x=waterperc, y=t70, color=species, size=density)) + geom_point() + facet_grid(~size)
```

-   alpha--- note: don't interpret yet

``` r
summary(mod.select.alpha) # density, Zn, C
```

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

Community as a predictor
------------------------

*Hyp:* Average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### Result:

-r2... none of the community components are significant predictors

``` r
rand.t.test(fit.r2.cv)
```

    ##             RMSE           R2    Avg.Bias  Max.Bias     Skill  delta.RMSE
    ## Comp01 0.1013107 1.494708e-05 0.002024841 0.2609272 -40.32548 18.45905577
    ## Comp02 0.1104027 4.293834e-03 0.003308000 0.2568255 -66.64251  8.97443525
    ## Comp03 0.1114560 1.027904e-02 0.003371373 0.2497766 -69.83735  0.95404167
    ## Comp04 0.1119467 1.116217e-02 0.002666967 0.2514263 -71.33590  0.44020287
    ## Comp05 0.1120426 1.133987e-02 0.002626418 0.2509619 -71.62965  0.08568655
    ##            p
    ## Comp01 0.893
    ## Comp02 0.999
    ## Comp03 0.757
    ## Comp04 0.801
    ## Comp05 0.709

-k... none of the community components are significant predictors

``` r
rand.t.test(fit.k.cv)
```

    ##              RMSE         R2    Avg.Bias  Max.Bias       Skill  delta.RMSE
    ## Comp01 0.09353450 0.05958082 0.001951813 0.2601135  -0.1449817  0.07246458
    ## Comp02 0.09959658 0.02582616 0.003520065 0.2696231 -13.5466706  6.48111904
    ## Comp03 0.10066062 0.01731354 0.003984269 0.2727955 -15.9857851  1.06835125
    ## Comp04 0.10059285 0.01921425 0.003727121 0.2728401 -15.8296569 -0.06732756
    ## Comp05 0.10067798 0.01871320 0.003714864 0.2731036 -16.0257809  0.08462470
    ##            p
    ## Comp01 0.488
    ## Comp02 0.998
    ## Comp03 0.887
    ## Comp04 0.418
    ## Comp05 0.730

-t70... none of the community components are significant predictors

``` r
rand.t.test(fit.t70.cv)
```

    ##             RMSE         R2     Avg.Bias Max.Bias      Skill delta.RMSE
    ## Comp01 0.5177469 0.07292880  0.001366029 1.280885  -0.374932 0.18729060
    ## Comp02 0.5561297 0.03410409 -0.004221057 1.345871 -15.809047 7.41343680
    ## Comp03 0.5635152 0.02574998 -0.006493032 1.393380 -18.905395 1.32801513
    ## Comp04 0.5636628 0.02712522 -0.005167941 1.382455 -18.967696 0.02619417
    ## Comp05 0.5637798 0.02665069 -0.004978126 1.383710 -19.017097 0.02076028
    ##            p
    ## Comp01 0.487
    ## Comp02 1.000
    ## Comp03 0.867
    ## Comp04 0.526
    ## Comp05 0.539

-alpha --- note: don't interpret yet

``` r
rand.t.test(fit.alpha.cv)
```

    ##             RMSE         R2     Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.2727951 0.02525264 -0.001168419 0.7257063 -10.84578 5.28332258
    ## Comp02 0.2801004 0.02756209  0.003844952 0.7085096 -16.86204 2.67793965
    ## Comp03 0.2845994 0.01564889  0.007017112 0.7064607 -20.64626 1.60620079
    ## Comp04 0.2866664 0.01260590  0.005775351 0.7015749 -22.40512 0.72629298
    ## Comp05 0.2869146 0.01259244  0.006101055 0.7009235 -22.61720 0.08659388
    ##            p
    ## Comp01 0.694
    ## Comp02 0.837
    ## Comp03 0.938
    ## Comp04 0.897
    ## Comp05 0.634

Community+traits as a predictor
-------------------------------

``` r
#run rioja::WAPLS on wood trait residuals

#make sure that the residuals and community dataframes are aligned
trait.residuals$code<-as.character(trait.residuals$code)
trait.residuals[!trait.residuals$code %in% row.names(meanOTUabund.trim),"code"] # no missing from meanOTUabund.trim
```

    ## character(0)

``` r
row.names(meanOTUabund.trim)[!row.names(meanOTUabund.trim) %in% trait.residuals$code] #missing olst from trait.residuals
```

    ## [1] "olst"

``` r
trait.residuals.trim<-trait.residuals[trait.residuals$code %in% row.names(meanOTUabund.trim),] #trim
meanOTUabund.trim2<-meanOTUabund.trim[row.names(meanOTUabund.trim) %in% trait.residuals.trim$code,] #trim
ord<-match(row.names(meanOTUabund.trim2), trait.residuals.trim$code)
trait.residuals.trim.o<-trait.residuals.trim[ord,]
sum(trait.residuals.trim.o$code != row.names(meanOTUabund.trim2)) #this needs to by 0
```

    ## [1] 0

``` r
#get rid of empty OTU cols
sum(colSums(meanOTUabund.trim2)!=0)
```

    ## [1] 3744

``` r
meanOTUabund.trim2<-meanOTUabund.trim2[,colSums(meanOTUabund.trim2)!=0]

#fit models
fit.tr.r2 <- WAPLS(meanOTUabund.trim2, trait.residuals.trim.o$r2.resid)
fit.tr.k <- WAPLS(meanOTUabund.trim2, trait.residuals.trim.o$k.resid)
fit.tr.t70 <- WAPLS(meanOTUabund.trim2, trait.residuals.trim.o$t70.resid)
fit.tr.alpha <- WAPLS(meanOTUabund.trim2, trait.residuals.trim.o$alpha.resid)

#cross-validate models using the leave-one-out method
fit.tr.r2.cv <- crossval(fit.tr.r2, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |========================                                         |  38%
      |                                                                       
      |==========================                                       |  41%
      |                                                                       
      |============================                                     |  44%
      |                                                                       
      |==============================                                   |  47%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |===================================                              |  53%
      |                                                                       
      |=====================================                            |  56%
      |                                                                       
      |=======================================                          |  59%
      |                                                                       
      |=========================================                        |  62%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |=================================================================| 100%

``` r
fit.tr.k.cv <- crossval(fit.tr.k, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |========================                                         |  38%
      |                                                                       
      |==========================                                       |  41%
      |                                                                       
      |============================                                     |  44%
      |                                                                       
      |==============================                                   |  47%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |===================================                              |  53%
      |                                                                       
      |=====================================                            |  56%
      |                                                                       
      |=======================================                          |  59%
      |                                                                       
      |=========================================                        |  62%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |=================================================================| 100%

``` r
fit.tr.t70.cv <- crossval(fit.tr.t70, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |========================                                         |  38%
      |                                                                       
      |==========================                                       |  41%
      |                                                                       
      |============================                                     |  44%
      |                                                                       
      |==============================                                   |  47%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |===================================                              |  53%
      |                                                                       
      |=====================================                            |  56%
      |                                                                       
      |=======================================                          |  59%
      |                                                                       
      |=========================================                        |  62%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |=================================================================| 100%

``` r
fit.tr.alpha.cv <- crossval(fit.tr.alpha, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |========================                                         |  38%
      |                                                                       
      |==========================                                       |  41%
      |                                                                       
      |============================                                     |  44%
      |                                                                       
      |==============================                                   |  47%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |===================================                              |  53%
      |                                                                       
      |=====================================                            |  56%
      |                                                                       
      |=======================================                          |  59%
      |                                                                       
      |=========================================                        |  62%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |=================================================================| 100%

``` r
#rand.t.test(fit.alpha.cv) #perform randomization t-test to test the significance of a cross-validated model
#screeplot(fit.r2.cv)
```

*Result:* -r2...none of the community components are significant predictors. Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was marginally significant

``` r
rand.t.test(fit.tr.r2.cv)
```

    ##              RMSE        R2     Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.07316014 0.1429569 -0.005347692 0.1618938  -72.28957 31.2591204
    ## Comp02 0.07945880 0.1883225 -0.004975824 0.1735680 -103.23289  8.6094227
    ## Comp03 0.07976155 0.2104281 -0.005255877 0.1731221 -104.78449  0.3810032
    ## Comp04 0.07959619 0.2088527 -0.004842101 0.1713720 -103.93626 -0.2073162
    ## Comp05 0.07947046 0.2110715 -0.004777034 0.1709423 -103.29250 -0.1579594
    ##            p
    ## Comp01 0.999
    ## Comp02 0.995
    ## Comp03 0.648
    ## Comp04 0.337
    ## Comp05 0.200

-k... none of the community components are significant predictors

``` r
rand.t.test(fit.tr.k.cv)
```

    ##              RMSE          R2    Avg.Bias  Max.Bias     Skill  delta.RMSE
    ## Comp01 0.05191325 0.001028608 0.003303599 0.1355113 -21.13298 10.06043001
    ## Comp02 0.05067798 0.013536731 0.002851440 0.1225496 -15.43686 -2.37949637
    ## Comp03 0.05092933 0.011624399 0.002815508 0.1227497 -16.58479  0.49598106
    ## Comp04 0.05088217 0.011202260 0.002606935 0.1222940 -16.36896 -0.09260487
    ## Comp05 0.05080847 0.011411770 0.002600391 0.1220352 -16.03213 -0.14482990
    ##            p
    ## Comp01 0.925
    ## Comp02 0.240
    ## Comp03 0.752
    ## Comp04 0.289
    ## Comp05 0.133

-t70... none of the community components are significant predictors. Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was a significant predictor

``` r
rand.t.test(fit.tr.t70.cv) 
```

    ##             RMSE         R2    Avg.Bias  Max.Bias     Skill  delta.RMSE
    ## Comp01 0.2752412 0.03752747 0.007269499 0.6434662 -5.866109  2.89125760
    ## Comp02 0.2749304 0.05477065 0.016445707 0.6775652 -5.627133 -0.11293078
    ## Comp03 0.2766449 0.05067283 0.018458392 0.6756985 -6.948652  0.62361403
    ## Comp04 0.2761444 0.05152058 0.018574458 0.6725779 -6.562016 -0.18092128
    ## Comp05 0.2759834 0.05155263 0.018441009 0.6738231 -6.437833 -0.05828514
    ##            p
    ## Comp01 0.663
    ## Comp02 0.478
    ## Comp03 0.775
    ## Comp04 0.224
    ## Comp05 0.256

-alpha --- note: don't interpret yet

``` r
rand.t.test(fit.tr.alpha.cv)
```

    ##             RMSE          R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1898628 0.003951823 -0.01449037 0.3991946 -32.51782 15.1163855
    ## Comp02 0.2027442 0.011228156 -0.01338738 0.4271724 -51.10947  6.7846166
    ## Comp03 0.2061391 0.017990710 -0.01584032 0.4289071 -56.21240  1.6744701
    ## Comp04 0.2067854 0.020289417 -0.01554990 0.4267151 -57.19334  0.3134854
    ## Comp05 0.2063872 0.020076149 -0.01537752 0.4273147 -56.58859 -0.1925460
    ##            p
    ## Comp01 0.902
    ## Comp02 0.970
    ## Comp03 0.900
    ## Comp04 0.838
    ## Comp05 0.133

*Note: no longer warrented based on analyses with waterperc represented as g water / g wet weight...* Investigate the biology underlying t70-associated coefs for Comp05

``` r
coef.comp5<-fit.tr.t70.cv$coefficients[,'Comp05']
coef.comp5.df<-data.frame(OTUId=names(coef.comp5), coefComp5=coef.comp5)

#create df of OTU taxon and guild info matched with coef value
coef.comp5.df %>%
  left_join(taxAndFunguild) -> coef.comp5.ann
```

By trophic mode (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

``` r
ggplot(coef.comp5.ann, aes(x=Trophic.Mode, y=coefComp5)) + geom_violin() + coord_flip()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-20-1.png) By guild (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

``` r
ggplot(coef.comp5.ann, aes(x=Guild, y=coefComp5)) + geom_violin() + coord_flip()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-21-1.png) By phylum

``` r
ggplot(coef.comp5.ann, aes(x=phylum, y=coefComp5)) + geom_violin() + coord_flip()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-22-1.png)

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-24-1.png)

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-25-1.png)

1.  Diversity plus traits as a predictor

############################################## 

Extra pieces
============

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2
