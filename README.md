Does chemistry or community better predict mass loss?
================
10/23/2017

### Load microbial community data

### Load wood trait data

### Load percent mass remaining (pmr) data

Plot pmr data ![](readme_files/figure-markdown_github/unnamed-chunk-5-1.png)

Calculate decay trajectory fits for each code (species+size)

Compare negative exp. vs weibull by plotting time to 70% mass remaining (t70) for each species+size ![](readme_files/figure-markdown_github/unnamed-chunk-7-1.png)

Check for missing stem-level data

    ## # A tibble: 2 x 15
    ## # Groups:   codeStem [2]
    ##   codeStem barkthick      C       Ca density        Fe       K      Mn
    ##      <chr>     <dbl>  <dbl>    <dbl>   <dbl>     <dbl>   <dbl>   <dbl>
    ## 1    acpa2        NA 50.659 3397.019      NA  3859.762 1400.18  38.376
    ## 2    lepa4        NA 48.377 2870.988      NA 14424.700  748.14 246.435
    ## # ... with 7 more variables: N <dbl>, P <dbl>, waterperc <dbl>, Zn <dbl>,
    ## #   species <chr>, size <chr>, code <chr>

    ## # A tibble: 2 x 7
    ## # Groups:   code [2]
    ##    code species  size seq_sampName  drop seq.stem codeStem
    ##   <chr>   <chr> <chr>        <chr> <chr>    <chr>    <chr>
    ## 1  acpa    acpa small        acpa2  acpa        2     <NA>
    ## 2  lepa    lepa small        lepa4  lepa        4     <NA>

    ## # A tibble: 0 x 10
    ## # Groups:   codeStem [0]
    ## # ... with 10 variables: codeStem <chr>, time0 <dbl>, time7 <dbl>,
    ## #   time13 <dbl>, time25 <dbl>, time37 <dbl>, time59 <dbl>, code <chr>,
    ## #   species <chr>, size <chr>

issue \#31 -- There are 2 unique codeStem ids that are found in the trait data (xrf sample names) and the sequence data, but not in the stemSamples data (deployment sample names). These codeStem ids are not found in the percent mass loss data. Because the main goal is to analyze decay responses, I'm going to continue to leave these codeStems out of the stemSamples dataframe. Is it possible that stem id numbers got switched? Something to follow-up on.

Figure 1. Wood species x decay params -- probably want to use weibull instead? ![](readme_files/figure-markdown_github/unnamed-chunk-9-1.png)

    ## quartz_off_screen 
    ##                 2

Figure 2. Time x percent mass remaining by wood species ![](readme_files/figure-markdown_github/unnamed-chunk-10-1.png)

    ## quartz_off_screen 
    ##                 2

########################################## 

Wood traits as a predictor
==========================

We expect initial wood traits will explain varitation in species+size decay rate (k and t70), species+size lagginess (alpha), and stem-level percent mass remaining at 7, 13, 25, and 37 months of decay. Specifically, we expect samples with (a) high water percent, (b) low density and total C, (c) high macro and micro nutrients, and (d) thicker bark (potential mech: limiting microbial colonization) to have faster decay and less lagginess.

*Hyp (species+size-level)* Species+size-level initial wood traits will predict variation decay rates and lagginess.
-------------------------------------------------------------------------------------------------------------------

    ## # A tibble: 13 x 4
    ##           term                  k                  t70               ne.r2
    ##  *       <chr>              <chr>                <chr>               <chr>
    ##  1 (Intercept)   1.328 +/- 0.54 *    -3.103 +/- 3.232  1.892 +/- 0.498 ***
    ##  2   sizesmall  0.062 +/- 0.02 ** -0.594 +/- 0.152 ***                <NA>
    ##  3   waterperc 0.005 +/- 0.002 ** -0.048 +/- 0.009 *** 0.005 +/- 0.001 ***
    ##  4   barkthick  -0.022 +/- 0.012      0.18 +/- 0.087 *                <NA>
    ##  5           C -0.023 +/- 0.009 *    0.142 +/- 0.064 *   -0.026 +/- 0.01 *
    ##  6           N 0.179 +/- 0.061 **  -1.374 +/- 0.463 **  0.182 +/- 0.055 **
    ##  7          Ca        0 +/- 0 ***          0 +/- 0 ***                <NA>
    ##  8          Zn     0.001 +/- 0 **  -0.011 +/- 0.003 **        0.001 +/- 0 
    ##  9        <NA>               <NA>                 <NA>                <NA>
    ## 10       Fstat              10.37                10.41                4.51
    ## 11       numdf                  9                    8                   5
    ## 12       dendf                 23                   24                  27
    ## 13   r.squared                0.8                 0.78                0.46

*Hyp (stem-level)* Stem-level initial wood traits will predict variation in percent mass loss at each time step.
----------------------------------------------------------------------------------------------------------------

First, we need to decide what trait data (and samples) to include in this analysis since we don't have full coverage of stem-level trait data. Density and bark thickness were only measured on small sized stems. If there is not be very much within-species variation in these traits that contribute to variation in percent mass loss than we can justify including species-level estimates of these traits in the stem-level model.

Plot the small-sized stem-level measures of density and barkthick ![](readme_files/figure-markdown_github/unnamed-chunk-12-1.png)

    ## quartz_off_screen 
    ##                 2

Compare model fits (r2) using stem and species-level data to identify how much information about percent mass remaining is lost by using species-level estimates...For density, it looks like stem-level data improves model fit a tiny bit for early percent mass remaining time points (after 7 and 13 months) but not later time points. For barkthickness, fits are about the same.

    ## $time7
    ## Cox test
    ## 
    ## Model 1: time7 ~ density_stem
    ## Model 2: time7 ~ density_code
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## fitted(M1) ~ M2  0.68487    0.56767  1.2065  0.22764  
    ## fitted(M2) ~ M1 -0.95311    0.39005 -2.4436  0.01454 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $time13
    ## Cox test
    ## 
    ## Model 1: time13 ~ density_stem
    ## Model 2: time13 ~ density_code
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2  0.24606    0.39085  0.6295   0.5290
    ## fitted(M2) ~ M1 -0.38016    0.28735 -1.3230   0.1858
    ## 
    ## $time25
    ## Cox test
    ## 
    ## Model 1: time25 ~ density_stem
    ## Model 2: time25 ~ density_code
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2  -0.3566    0.68154 -0.5232   0.6008
    ## fitted(M2) ~ M1  -0.2043    0.69623 -0.2934   0.7692
    ## 
    ## $time37
    ## Cox test
    ## 
    ## Model 1: time37 ~ density_stem
    ## Model 2: time37 ~ density_code
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2 -0.23217    0.59967 -0.3872   0.6986
    ## fitted(M2) ~ M1 -0.18912    0.60408 -0.3131   0.7542
    ## 
    ## $time59
    ## Cox test
    ## 
    ## Model 1: time59 ~ density_stem
    ## Model 2: time59 ~ density_code
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## fitted(M1) ~ M2  0.55009    0.48034  1.1452  0.25212  
    ## fitted(M2) ~ M1 -0.73549    0.31308 -2.3492  0.01881 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## $time7
    ## Cox test
    ## 
    ## Model 1: time7 ~ barkthick_stem
    ## Model 2: time7 ~ barkthick_code
    ##                  Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2 -0.137323    0.36348 -0.3778   0.7056
    ## fitted(M2) ~ M1 -0.044232    0.38730 -0.1142   0.9091
    ## 
    ## $time13
    ## Cox test
    ## 
    ## Model 1: time13 ~ barkthick_stem
    ## Model 2: time13 ~ barkthick_code
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## fitted(M1) ~ M2 -0.32781    0.16294 -2.0118  0.04424 *
    ## fitted(M2) ~ M1  0.23141    0.35594  0.6501  0.51560  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $time25
    ## Cox test
    ## 
    ## Model 1: time25 ~ barkthick_stem
    ## Model 2: time25 ~ barkthick_code
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2 -0.60720    0.38652 -1.5709   0.1162
    ## fitted(M2) ~ M1  0.29202    0.57455  0.5083   0.6113
    ## 
    ## $time37
    ## Cox test
    ## 
    ## Model 1: time37 ~ barkthick_stem
    ## Model 2: time37 ~ barkthick_code
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2 -0.72103    0.83575 -0.8627   0.3883
    ## fitted(M2) ~ M1 -0.27998    0.88766 -0.3154   0.7524
    ## 
    ## $time59
    ## Cox test
    ## 
    ## Model 1: time59 ~ barkthick_stem
    ## Model 2: time59 ~ barkthick_code
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## fitted(M1) ~ M2 -0.96335    0.73903 -1.3035   0.1924
    ## fitted(M2) ~ M1  0.11858    0.86611  0.1369   0.8911

    ##   respvars density_code_r2 density_stem_r2 barkthick_code_r2
    ## 1    time7            0.04            0.08              0.02
    ## 2   time13            0.02            0.04              0.02
    ## 3   time25            0.09            0.09              0.04
    ## 4   time37            0.08            0.07              0.10
    ## 5   time59            0.02            0.05              0.10
    ##   barkthick_stem_r2
    ## 1              0.02
    ## 2              0.00
    ## 3              0.02
    ## 4              0.09
    ## 5              0.07

Stem-level density estimates provide additional information about mass loss beyond code-level density at time7 and time 59. Stem-level bark thickness estimates are not useful (beyond code-level bark thickness) at any time point.

Compile a "stem-level" dataframe with (a) stem-level percent mass remaining values, (b) stem-level traits including waterperc and chemistry along, and (c) small species-level density and bark thickness data.

########################################## 

Community as a predictor
========================

Filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 150 of 6128 OTUs"

*Hyp (species+size-level)* Species+size-level (average) initial microbial community composition will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##                   stat      k    t70  ne.r2 trim
    ## RMSE              RMSE   0.10   0.67   0.08  yes
    ## R2                  R2   0.04   0.00   0.00  yes
    ## Avg.Bias      Avg.Bias  -0.01   0.07  -0.01  yes
    ## Max.Bias      Max.Bias   0.27   1.72   0.20  yes
    ## Skill            Skill -58.60 -23.69 -53.07  yes
    ## delta.RMSE  delta.RMSE  25.94  11.22  23.72  yes
    ## p                    p   1.00   0.98   0.99  yes
    ## RMSE1             RMSE   0.09   0.60   0.07   no
    ## R21                 R2   0.01   0.05   0.01   no
    ## Avg.Bias1     Avg.Bias   0.00   0.02   0.00   no
    ## Max.Bias1     Max.Bias   0.25   1.70   0.22   no
    ## Skill1           Skill -14.56  -1.51 -17.09   no
    ## delta.RMSE1 delta.RMSE   7.03   0.75   8.21   no
    ## p1                   p   0.77   0.55   0.88   no

*Hyp (stem-level)* Stem-level initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
---------------------------------------------------------------------------------------------------------------------------------------------------------------

Comp01 (of the non-trimmed community) is a significant predictor of percent mass remaining at 37 months.

Plot the distribution of WA-PLS scores ![](readme_files/figure-markdown_github/unnamed-chunk-19-1.png)

Who is in the top and bottom 1%?

    ##        quant kingdom        phylum                    species Trophic.Mode
    ## 1  bottom 1%   Fungi    Ascomycota         Hormonema_viticola   Saprotroph
    ## 2  bottom 1%   Fungi Basidiomycota         Hyphodontia_radula   Saprotroph
    ## 3  bottom 1%   Fungi    Ascomycota     Neosetophoma_samarorum   Saprotroph
    ## 4  bottom 1%   Fungi    Ascomycota  Neophysalospora_eucalypti unclassified
    ## 5  bottom 1%   Fungi Basidiomycota      Bensingtonia_ingoldii   Saprotroph
    ## 6  bottom 1%   Fungi Basidiomycota   Pisolithus_croceorrhizus  Symbiotroph
    ## 7  bottom 1%   Fungi    Ascomycota  Mycosphaerella_excentrica   Pathotroph
    ## 8  bottom 1%   Fungi Basidiomycota      Mycetinis_scorodonius   Saprotroph
    ## 9  bottom 1%   Fungi Basidiomycota            Odontia_fibrosa   Saprotroph
    ## 10 bottom 1%   Fungi    Ascomycota     Acremonium_cavaraeanum unclassified
    ## 11 bottom 1%   Fungi    Ascomycota    Phaeomoniella_prunicola   Saprotroph
    ## 12    top 1%   Fungi    Ascomycota Debaryomyces_vindobonensis unclassified
    ## 13    top 1%   Fungi Basidiomycota       Septobasidium_burtii   Pathotroph
    ##                   Guild
    ## 1  Undefined Saprotroph
    ## 2  Undefined Saprotroph
    ## 3  Undefined Saprotroph
    ## 4          unclassified
    ## 5  Undefined Saprotroph
    ## 6       Ectomycorrhizal
    ## 7        Plant Pathogen
    ## 8  Undefined Saprotroph
    ## 9  Undefined Saprotroph
    ## 10         unclassified
    ## 11 Undefined Saprotroph
    ## 12         unclassified
    ## 13      Animal Pathogen

Many of the bottom 1% OTUs are classified as saprotrophs. That makes sense since low WA-PLS scores indicate an association with high mass loss (i.e. less mass remaining) at time37.

But saprotrophs are also found at many points along the gradient... ![](readme_files/figure-markdown_github/unnamed-chunk-21-1.png)

    ## quartz_off_screen 
    ##                 2

Is this because there is an underlying signature of wood traits on the initial microbial community that is driving the relationship between the community and the mass remaining after 37 months? The next analysis ("Community+traits" as predictor) will test this formally. Just out of curiousity, I'd like to pull in OTU "niche" info from the boral analysis to see if there's a relationship between OTU WA-PLS scores and wood trait coeffient estimates.

Reminder of which wood traits were included in the best model to explain pmr at time37...

    ##            X            term               time7               time13
    ## 1          1     (Intercept) 1.114 +/- 0.055 ***  1.253 +/- 0.063 ***
    ## 2          2       sizesmall -0.088 +/- 0.02 *** -0.111 +/- 0.019 ***
    ## 3          3       waterperc -0.003 +/- 0.001 ** -0.006 +/- 0.001 ***
    ## 4          4 barkthick_smspp                <NA>                 <NA>
    ## 5          5               C                <NA>                 <NA>
    ## 6          6               N  -0.098 +/- 0.043 *  -0.19 +/- 0.051 ***
    ## 7          7               P                <NA>            0 +/- 0 *
    ## 8          8            <NA>                <NA>                 <NA>
    ## 9      Fstat           Fstat                6.42                 12.9
    ## 10     numdf           numdf                   5                    5
    ## 11     dendf           dendf                  49                   48
    ## 12 r.squared       r.squared                 0.4                 0.57
    ##                  time25               time37               time59
    ## 1   1.216 +/- 0.096 ***  1.359 +/- 0.126 ***    -1.323 +/- 0.718 
    ## 2   -0.099 +/- 0.034 **                 <NA>                 <NA>
    ## 3  -0.011 +/- 0.002 *** -0.017 +/- 0.003 *** -0.018 +/- 0.003 ***
    ## 4                  <NA>                 <NA>  -0.056 +/- 0.021 **
    ## 5                  <NA>                 <NA>  0.052 +/- 0.014 ***
    ## 6                  <NA>                 <NA>                 <NA>
    ## 7        0.001 +/- 0 **       0.001 +/- 0 **                 <NA>
    ## 8                  <NA>                 <NA>                 <NA>
    ## 9                  9.91                15.32                12.41
    ## 10                    5                    4                    4
    ## 11                   46                   46                   48
    ## 12                 0.52                 0.57                 0.51

More water leads to less mass remaining; more P leads to more mass remaining

Plot OTU wood trait estimates (from boral) versus signif WA-PLS score. ![](readme_files/figure-markdown_github/unnamed-chunk-23-1.png)

    ## 
    ## Call:
    ## lm(formula = coefEst ~ coefComp, data = tmp)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -8.1111 -0.8393  0.1895  1.1350  3.5687 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   0.8647     0.5203   1.662   0.0989 .
    ## coefComp     -2.4314     1.0542  -2.306   0.0226 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.783 on 136 degrees of freedom
    ## Multiple R-squared:  0.03764,    Adjusted R-squared:  0.03057 
    ## F-statistic:  5.32 on 1 and 136 DF,  p-value: 0.0226

There's a weak negative relationship between an OTU's WA-PLS score and waterperc coefficient (slope=-2.3, p=.03), suggesting that OTUs that "prefer" high-water niche space are associated with less mass remaining at time37.

########################################## 

Community+traits as a predictor
===============================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Community data doesn't improve our understanding of decay rates (k, t70) or variation in decay rate (ne.r2) beyond what is known from the trait data.

*Hyp (stem-level)* After accounting for variation in decay due to wood traits (no models with density, includes small-species level bark thickness), stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Community data doesn't improve our understanding of mass loss (pmr after 7, 13, 25, 37, and 59 months) beyond what is known from the trait data.

########################################## 

Diversity (and diversity of specific clades) as a predictor
===========================================================

**Note that the full community matrix was used for these analyses**

*Hyp-a (species+size-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (species+size-level)* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (species+size-level)* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers.

    ##        term    source respvar               coef signif
    ## 1 sizesmall Oomy.rich     t70 -0.525 +/- 0.228 *   TRUE

None of the richness variables were significant predictors of k, ne.r2, or t70

*Hyp-a (stem-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to less mass remaining esp. at early time steps because of the selection effect for fast decayers and complementarity among taxa for decay.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to more mass remaining because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (stem-level)* Greater saprotroph and basidiomycete richness will lead to less mass remaining esp. at early time steps because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to more mass remaining because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (stem-level)* Greater pathogen and oomycete richness will lead to more mass remaining because the presence of these organisms will inhibit the establishment and activity of decayers.

    ##                  term     source respvar                coef signif
    ## 1           sizesmall  Oomy.rich  time13 -0.078 +/- 0.02 ***   TRUE
    ## 2  sizesmall:sub_rich   Richness  time25  -0.001 +/- 0.001 *   TRUE
    ## 3            sub_rich   Richness  time25       0.001 +/- 0 *   TRUE
    ## 4            sub_rich  ShannonsH  time25   0.081 +/- 0.039 *   TRUE
    ## 5  sizesmall:sub_rich Sapro.rich  time25  -0.014 +/- 0.005 *   TRUE
    ## 6            sub_rich  Oomy.rich  time25   0.101 +/- 0.043 *   TRUE
    ## 7  sizesmall:sub_rich   Richness  time37  -0.002 +/- 0.001 *   TRUE
    ## 8  sizesmall:sub_rich  ShannonsH  time37   -0.16 +/- 0.074 *   TRUE
    ## 9  sizesmall:sub_rich Sapro.rich  time37  -0.022 +/- 0.009 *   TRUE
    ## 10 sizesmall:sub_rich Patho.rich  time37  -0.017 +/- 0.007 *   TRUE
    ## 11          sizesmall  Oomy.rich   time7 -0.046 +/- 0.016 **   TRUE

A couple of models were informative, particularly for time37. Saprotroph OTU richness is associated with percent mass remaining at time25 and time37. Pathogen OTU richness and Shannon's H are associated with pmr at time 37. In all cases, there is a negative interaction between stem size and diversity.

Plot the relationship between saprotroph OTU richness and pmr at time25 and time37 ![](readme_files/figure-markdown_github/unnamed-chunk-29-1.png) More saprotrophs leads to less mass remaining - but only in small stems.

Plot the relationship between pathogen OTU richness and pmr at time37 ![](readme_files/figure-markdown_github/unnamed-chunk-30-1.png) More pathotrophs leads to less mass remaining - but only in small stems.

Plot the relationship between Shannon's H and pmr at time37 ![](readme_files/figure-markdown_github/unnamed-chunk-31-1.png) Higher Shannon's H leads to less mass remaining in small stems, more mass remaining in large stems... but these are pretty weak relationships.

############################################## 

Diversity plus traits as a predictor
====================================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##   term     source respvar              coef signif
    ## 1 mean Sapro.rich       k 0.005 +/- 0.002 *   TRUE
    ## 2 mean Patho.rich   ne.r2 0.005 +/- 0.002 *   TRUE
    ## 3 mean Sapro.rich     t70 -0.04 +/- 0.017 *   TRUE

Beyond wood trait data, saprotroph richness improves models for k and t70. In addition, pathotroph richness improves the model for ne.r2.

Plot the relationship between Saprotroph richness and k ![](readme_files/figure-markdown_github/unnamed-chunk-33-1.png) Higher saprotroph richness leads to faster decay than would be expected based on wood traits alone.

Plot the relationship between Saprotroph richness and t70 ![](readme_files/figure-markdown_github/unnamed-chunk-34-1.png) Higher saprotroph richness leads to shorter times to 70% mass remaining than would be expected based on wood traits alone.

Plot the relationship between Pathotroph richness and ne.r2 ![](readme_files/figure-markdown_github/unnamed-chunk-35-1.png) Higher pathotroph richness leads to better-fitting decay models than would be expected based on wood traits alone.

*Hyp (stem-level)* After accounting for variation in decay due to wood traits, initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in percent mass loss, esp. at early time points.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##                 term     source respvar               coef signif
    ## 1 sizesmall:sub_rich   Richness  time13      0.001 +/- 0 *   TRUE
    ## 2 sizesmall:sub_rich   Richness  time25 -0.001 +/- 0.001 *   TRUE
    ## 3           sub_rich  Oomy.rich  time25  0.076 +/- 0.029 *   TRUE
    ## 4 sizesmall:sub_rich Patho.rich  time37 -0.012 +/- 0.005 *   TRUE
    ## 5          sizesmall Sapro.rich  time59 -0.221 +/- 0.106 *   TRUE
    ## 6 sizesmall:sub_rich  Oomy.rich   time7  -0.09 +/- 0.036 *   TRUE

Overall OTU richness improves model estimates for percent mass remaining at 25 months. Pathotroph richness improves model estmates for pmr at 37 months.

Plot the relationship between OTU richness and residuals of pmr at time25 ![](readme_files/figure-markdown_github/unnamed-chunk-37-1.png) In small stems, more fungal OTUs leads to less mass remaining after 25 months than would be expected based on wood traits alone. This relationship looks heavily influenced by 1 outlier w/ high OTU richness

Plot the relationship between Pathotroph OTU richness and residuals of pmr at time37 ![](readme_files/figure-markdown_github/unnamed-chunk-38-1.png) In small stems, more pathotrophs leads to less mass remaining after 37 months than would be expected based on wood traits alone.

############################################## 

Relationship between wood traits and community
==============================================

*Hyp (species+size-level)* Initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

![](readme_files/figure-markdown_github/unnamed-chunk-39-1.png)

    ## quartz_off_screen 
    ##                 2

Full community anova-like table

    ##                    term    df Fval  pval
    ## 1                  size  1.00 1.40 0.015
    ## 2             barkthick  1.00 1.60 0.007
    ## 3             waterperc  1.00 1.41 0.012
    ## 4              Residual 29.00   NA    NA
    ## 5                  <NA>    NA   NA    NA
    ## 6   Constrained inertia  1.91   NA    NA
    ## 7 Unconstrained inertia 11.16   NA    NA

Trimmed community anova-like table

    ##                    term    df Fval  pval
    ## 1                  size  1.00 1.67 0.011
    ## 2             waterperc  1.00 1.58 0.027
    ## 3             barkthick  1.00 1.54 0.030
    ## 4              Residual 29.00   NA    NA
    ## 5                  <NA>    NA   NA    NA
    ## 6   Constrained inertia  1.79   NA    NA
    ## 7 Unconstrained inertia  9.22   NA    NA

*Hyp (stem-level)* Average initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

    ## # A tibble: 6 x 14
    ##    code barkthick        C       Ca   density         Fe         K
    ##   <chr>     <dbl>    <dbl>    <dbl>     <dbl>      <dbl>     <dbl>
    ## 1  acel  2.050000 51.17750 2806.772 0.5589542 2295.93950  651.4925
    ## 2  acpa  1.680000 50.68325 5266.820 0.6999979 3522.06900 1666.8053
    ## 3  ACPA  1.680000 49.39300 3303.239 0.6999979   23.98167  744.6377
    ## 4  alli  2.503333 49.89400 7652.384 0.6188021 1942.10533 2928.1150
    ## 5  ALLI  2.503333 50.06100 4842.335 0.6188021   35.51367  775.9320
    ## 6  anba  1.380000 50.06300 6730.355 0.5047438 5951.61700 3214.2055
    ## # ... with 7 more variables: Mn <dbl>, N <dbl>, P <dbl>, waterperc <dbl>,
    ## #   Zn <dbl>, species <chr>, size <chr>

    ## barkthick_smspp       waterperc       sizesmall              Ca 
    ##        3.569130        2.608819        1.913125        2.458129 
    ##   density_smspp               C               N              Mn 
    ##        2.631552        1.709383        3.624903        2.200966 
    ##               K 
    ##        2.152031

    ##       waterperc       sizesmall               C barkthick_smspp 
    ##        2.608819        1.913125        1.709383        3.569130 
    ##               N               K              Ca   density_smspp 
    ##        3.624903        2.152031        2.458129        2.631552 
    ##              Mn 
    ##        2.200966

![](readme_files/figure-markdown_github/unnamed-chunk-42-1.png)

    ## quartz_off_screen 
    ##                 2

Full community anova-like table

    ##                     term    df Fval  pval
    ## 1        barkthick_smspp  1.00 2.33 0.001
    ## 2              waterperc  1.00 1.99 0.001
    ## 3                   size  1.00 2.13 0.001
    ## 4                     Ca  1.00 1.91 0.001
    ## 5          density_smspp  1.00 1.97 0.001
    ## 6                      C  1.00 1.96 0.001
    ## 7                      N  1.00 2.02 0.001
    ## 8                     Mn  1.00 1.62 0.002
    ## 9                      K  1.00 1.60 0.004
    ## 10              Residual 47.00   NA    NA
    ## 11                  <NA>    NA   NA    NA
    ## 12   Constrained inertia  7.55   NA    NA
    ## 13 Unconstrained inertia 17.37   NA    NA

Trimmed community anova-like table

    ##                     term    df Fval  pval
    ## 1              waterperc  1.00 2.37 0.001
    ## 2                   size  1.00 2.60 0.001
    ## 3                      C  1.00 2.70 0.001
    ## 4        barkthick_smspp  1.00 2.64 0.001
    ## 5                      N  1.00 2.02 0.005
    ## 6                      K  1.00 1.59 0.026
    ## 7                     Ca  1.00 1.93 0.003
    ## 8          density_smspp  1.00 1.68 0.011
    ## 9                     Mn  1.00 1.61 0.016
    ## 10              Residual 47.00   NA    NA
    ## 11                  <NA>    NA   NA    NA
    ## 12   Constrained inertia  7.18   NA    NA
    ## 13 Unconstrained inertia 14.86   NA    NA

############################################## 

Extra pieces
------------

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
