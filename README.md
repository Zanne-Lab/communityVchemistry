Does chemistry or community better predict mass loss?
================
10/23/2017

    ## Warning: package 'dplyr' was built under R version 3.5.1

### Load microbial community data

Figure S1. Sample-effort curves

### Load wood trait data

### Load percent mass remaining (pmr) data

Quick plot of pmr data ![](readme_files/figure-markdown_github/unnamed-chunk-6-1.png)

Calculate decay trajectory fits for each code (species+size)

Compare negative exp. vs weibull by plotting time to 70% mass remaining (t70) for each species+size ![](readme_files/figure-markdown_github/unnamed-chunk-8-1.png)

Another view... Figure S2. Comparison of negative exp and weibull models using t70 ![](readme_files/figure-markdown_github/unnamed-chunk-9-1.png)

Figure S3. Comparison of negative exp and weibull models using AIC ![](readme_files/figure-markdown_github/unnamed-chunk-10-1.png)

Check for missing stem-level data

    ## # A tibble: 2 x 15
    ## # Groups:   codeStem [2]
    ##   codeStem barkthick     C    Ca density     Fe     K    Mn     N     P
    ##   <chr>        <dbl> <dbl> <dbl>   <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 acpa2           NA  50.7 3397.      NA  3860. 1400.  38.4 0.579  78.4
    ## 2 lepa4           NA  48.4 2871.      NA 14425.  748. 246.  0.391 147. 
    ## # ... with 5 more variables: waterperc <dbl>, Zn <dbl>, species <chr>,
    ## #   size <chr>, code <chr>

    ## # A tibble: 2 x 7
    ## # Groups:   code [2]
    ##   code  species size  seq_sampName drop  seq.stem codeStem
    ##   <chr> <chr>   <chr> <chr>        <chr> <chr>    <chr>   
    ## 1 acpa  acpa    small acpa2        acpa  2        <NA>    
    ## 2 lepa  lepa    small lepa4        lepa  4        <NA>

    ## # A tibble: 0 x 10
    ## # Groups:   codeStem [0]
    ## # ... with 10 variables: codeStem <chr>, time0 <dbl>, time7 <dbl>,
    ## #   time13 <dbl>, time25 <dbl>, time37 <dbl>, time59 <dbl>, code <chr>,
    ## #   species <chr>, size <chr>

issue \#31 -- There are 2 unique codeStem ids that are found in the trait data (xrf sample names) and the sequence data, but not in the stemSamples data (deployment sample names). These codeStem ids are not found in the percent mass loss data. Because the main goal is to analyze decay responses, I'm going to continue to leave these codeStems out of the stemSamples dataframe. Is it possible that stem id numbers got switched? Something to follow-up on.

Figure 1. Wood species x decay params (weibull model)

    ## quartz_off_screen 
    ##                 2

Figure 2. Time x percent mass remaining

    ## [1]  0  7 13 25 37 59

![](readme_files/figure-markdown_github/unnamed-chunk-13-1.png)

    ## quartz_off_screen 
    ##                 2

########################################## 

Wood traits as a predictor
==========================

We expect initial wood traits will explain varitation in species+size decay rate (k and t70), species+size lagginess (alpha), and stem-level percent mass remaining at 7, 13, 25, and 37 months of decay. Specifically, we expect samples with (a) high water percent, (b) low density and total C, (c) high macro and micro nutrients, and (d) thicker bark (potential mech: limiting microbial colonization) to have faster decay and less lagginess.

*Hyp (species+size-level)* Species+size-level initial wood traits will predict variation decay rates and lagginess.
-------------------------------------------------------------------------------------------------------------------

*Hyp (stem-level)* Stem-level initial wood traits will predict variation in percent mass loss at each time step.
----------------------------------------------------------------------------------------------------------------

First, we need to decide what trait data (and samples) to include in this analysis since we don't have full coverage of stem-level trait data. Density and bark thickness were only measured on small sized stems. If there is not be very much within-species variation in these traits that contribute to variation in percent mass loss than we can justify including species-level estimates of these traits in the stem-level model.

Plot the small-sized stem-level measures of density and barkthick ![](readme_files/figure-markdown_github/unnamed-chunk-15-1.png)

    ## quartz_off_screen 
    ##                 2

Compare model fits (r2) using stem and species-level data to identify how much information about percent mass remaining is lost by using species-level estimates...For density, it looks like stem-level data improves model fit a tiny bit for early percent mass remaining time points (after 7 and 13 months) but not later time points. For barkthickness, fits are about the same.

Stem-level density estimates provide additional information about mass loss beyond code-level density at time7 and time 59. Stem-level bark thickness estimates are not useful (beyond code-level bark thickness) at any time point.

Compile a "stem-level" dataframe with (a) stem-level percent mass remaining values, (b) stem-level traits including waterperc and chemistry along, and (c) small species-level density and bark thickness data.

########################################## 

Community as a predictor
========================

Filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 150 of 6128 OTUs"

*Hyp (species+size-level)* Species+size-level (average) initial microbial community composition will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##                   stat      k    t70  ne.r2  alpha   beta  w.t70   w.r2
    ## RMSE              RMSE   0.10   0.67   0.08   0.41   1.64   0.59   0.10
    ## R2                  R2   0.05   0.00   0.00   0.01   0.01   0.00   0.02
    ## Avg.Bias      Avg.Bias  -0.01   0.07  -0.01  -0.07   0.20   0.03  -0.01
    ## Max.Bias      Max.Bias   0.27   1.73   0.20   0.79   4.02   1.39   0.23
    ## Skill            Skill -59.73 -24.66 -52.94 -78.76 -38.87 -16.61 -79.97
    ## delta.RMSE  delta.RMSE  26.38  11.65  23.67  33.70  17.84   7.98  34.15
    ## p                    p   0.99   0.97   1.00   0.98   1.00   0.92   1.00
    ## RMSE1             RMSE   0.09   0.61   0.07   0.35   1.48   0.54   0.08
    ## R21                 R2   0.01   0.04   0.01   0.00   0.01   0.06   0.00
    ## Avg.Bias1     Avg.Bias   0.00   0.02   0.00  -0.04   0.10  -0.01  -0.01
    ## Max.Bias1     Max.Bias   0.25   1.70   0.22   0.95   4.39   1.36   0.25
    ## Skill1           Skill -15.38  -2.48 -16.89 -30.59 -13.24   1.31 -27.53
    ## delta.RMSE1 delta.RMSE   7.42   1.23   8.12  14.27   6.42  -0.66  12.93
    ## p1                   p   0.80   0.59   0.86   0.93   0.73   0.45   0.94
    ##             trim
    ## RMSE         yes
    ## R2           yes
    ## Avg.Bias     yes
    ## Max.Bias     yes
    ## Skill        yes
    ## delta.RMSE   yes
    ## p            yes
    ## RMSE1         no
    ## R21           no
    ## Avg.Bias1     no
    ## Max.Bias1     no
    ## Skill1        no
    ## delta.RMSE1   no
    ## p1            no

*Hyp (stem-level)* Stem-level initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
---------------------------------------------------------------------------------------------------------------------------------------------------------------

Comp01 (of the non-trimmed community) is a significant predictor of percent mass remaining at 37 months.

Plot the distribution of WA-PLS scores ![](readme_files/figure-markdown_github/unnamed-chunk-22-1.png)

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

But saprotrophs are also found at many points along the gradient... ![](readme_files/figure-markdown_github/unnamed-chunk-24-1.png)

    ## quartz_off_screen 
    ##                 2

Is this because there is an underlying signature of wood traits on the initial microbial community that is driving the relationship between the community and the mass remaining after 37 months? The next analysis ("Community+traits" as predictor) will test this formally. Just out of curiousity, I'd like to pull in OTU "niche" info from the boral analysis to see if there's a relationship between OTU WA-PLS scores and wood trait coeffient estimates.

Reminder of which wood traits were included in the best model to explain pmr at time37...

    ##            X            term      time13     time25       time37
    ## 1          1     (Intercept)     1.3 ***    1.2 ***      1.4 ***
    ## 2          2 barkthick_smspp        <NA>       <NA>         <NA>
    ## 3          3               C        <NA>       <NA>         <NA>
    ## 4          4              Fe        <NA>  -1.6e-05  -4.1e-05 ***
    ## 5          5               K        <NA>       <NA>    -3.7e-05 
    ## 6          6               N   -0.19 ***       <NA>         <NA>
    ## 7          7               P   0.00027 * 0.00056 **   0.00057 **
    ## 8          8       sizesmall   -0.11 ***  -0.099 **         <NA>
    ## 9          9       waterperc -0.0063 *** -0.011 ***   -0.017 ***
    ## 10        10              Zn   -0.00052   -0.00081          <NA>
    ## 11        11            <NA>        <NA>       <NA>         <NA>
    ## 12     Fstat           Fstat        12.9       9.91        15.32
    ## 13     numdf           numdf           5          5            4
    ## 14     dendf           dendf          48         46           46
    ## 15 r.squared       r.squared        0.57       0.52         0.57
    ##         time59      time7
    ## 1       -1.5 *    1.1 ***
    ## 2     -0.055 *       <NA>
    ## 3    0.055 ***       <NA>
    ## 4  -2.9e-05 **   8.3e-06 
    ## 5         <NA>  1.4e-05 *
    ## 6         <NA>   -0.098 *
    ## 7         <NA>       <NA>
    ## 8         <NA> -0.088 ***
    ## 9   -0.017 ***  -0.003 **
    ## 10        <NA>       <NA>
    ## 11        <NA>       <NA>
    ## 12       11.91       6.42
    ## 13           4          5
    ## 14          49         49
    ## 15        0.49        0.4

More water leads to less mass remaining; more P leads to more mass remaining

Plot OTU wood trait estimates (from boral) versus signif WA-PLS score. ![](readme_files/figure-markdown_github/unnamed-chunk-26-1.png)

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

    ##        term     source respvar     coef signif
    ## 1 sizesmall Patho.rich   alpha  -0.56 *   TRUE
    ## 2 sizesmall  Oomy.rich     t70  -0.52 *   TRUE
    ## 3 sizesmall  Oomy.rich   w.t70 -0.58 **   TRUE

*Hyp-a (stem-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to less mass remaining esp. at early time steps because of the selection effect for fast decayers and complementarity among taxa for decay.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to more mass remaining because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (stem-level)* Greater saprotroph and basidiomycete richness will lead to less mass remaining esp. at early time steps because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to more mass remaining because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (stem-level)* Greater pathogen and oomycete richness will lead to more mass remaining because the presence of these organisms will inhibit the establishment and activity of decayers.

    ##                  term     source respvar       coef signif
    ## 1           sizesmall  Oomy.rich  time13 -0.078 ***   TRUE
    ## 2  sizesmall:sub_rich   Richness  time25  -0.0014 *   TRUE
    ## 3            sub_rich   Richness  time25  0.00087 *   TRUE
    ## 4            sub_rich  ShannonsH  time25    0.081 *   TRUE
    ## 5  sizesmall:sub_rich Sapro.rich  time25   -0.014 *   TRUE
    ## 6            sub_rich  Oomy.rich  time25      0.1 *   TRUE
    ## 7  sizesmall:sub_rich   Richness  time37  -0.0018 *   TRUE
    ## 8  sizesmall:sub_rich  ShannonsH  time37    -0.16 *   TRUE
    ## 9  sizesmall:sub_rich Sapro.rich  time37   -0.022 *   TRUE
    ## 10 sizesmall:sub_rich Patho.rich  time37   -0.017 *   TRUE
    ## 11          sizesmall  Oomy.rich   time7  -0.046 **   TRUE

A couple of models were informative, particularly for time37. Saprotroph OTU richness is associated with percent mass remaining at time25 and time37. Pathogen OTU richness and Shannon's H are associated with pmr at time 37. In all cases, there is a negative interaction between stem size and diversity.

Plot the relationship between saprotroph OTU richness and pmr at time25 and time37 ![](readme_files/figure-markdown_github/unnamed-chunk-32-1.png) More saprotrophs leads to less mass remaining - but only in small stems.

Plot the relationship between pathogen OTU richness and pmr at time37 ![](readme_files/figure-markdown_github/unnamed-chunk-33-1.png) More pathotrophs leads to less mass remaining - but only in small stems.

Plot the relationship between Shannon's H and pmr at time37 ![](readme_files/figure-markdown_github/unnamed-chunk-34-1.png) Higher Shannon's H leads to less mass remaining in small stems, more mass remaining in large stems... but these are pretty weak relationships.

############################################## 

Diversity plus traits as a predictor
====================================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##   term     source respvar     coef signif
    ## 1 mean Sapro.rich    beta -0.12 **   TRUE
    ## 2 mean Sapro.rich       k 0.0053 *   TRUE
    ## 3 mean Patho.rich   ne.r2 0.0047 *   TRUE
    ## 4 mean Sapro.rich     t70  -0.04 *   TRUE

Beyond wood trait data, saprotroph richness improves models for k and t70. In addition, pathotroph richness improves the model for ne.r2.

Plot the relationship between Saprotroph richness and k ![](readme_files/figure-markdown_github/unnamed-chunk-36-1.png) Higher saprotroph richness leads to faster decay than would be expected based on wood traits alone.

Plot the relationship between Saprotroph richness and t70 ![](readme_files/figure-markdown_github/unnamed-chunk-37-1.png) Higher saprotroph richness leads to shorter times to 70% mass remaining than would be expected based on wood traits alone.

Plot the relationship between Saprotroph richness and beta ![](readme_files/figure-markdown_github/unnamed-chunk-38-1.png) Higher saprotroph richness leads to shorter times to 70% mass remaining than would be expected based on wood traits alone.

Plot the relationship between Pathotroph richness and ne.r2 ![](readme_files/figure-markdown_github/unnamed-chunk-39-1.png) Higher pathotroph richness leads to better-fitting decay models than would be expected based on wood traits alone.

*Hyp (stem-level)* After accounting for variation in decay due to wood traits, initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in percent mass loss, esp. at early time points.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##                 term     source respvar      coef signif
    ## 1 sizesmall:sub_rich   Richness  time13 0.00075 *   TRUE
    ## 2 sizesmall:sub_rich   Richness  time25 -0.0011 *   TRUE
    ## 3           sub_rich  Oomy.rich  time25   0.076 *   TRUE
    ## 4 sizesmall:sub_rich Patho.rich  time37  -0.012 *   TRUE
    ## 5          sizesmall Sapro.rich  time59   -0.25 *   TRUE
    ## 6 sizesmall:sub_rich Sapro.rich  time59   0.013 *   TRUE
    ## 7 sizesmall:sub_rich  Oomy.rich   time7   -0.09 *   TRUE

Overall OTU richness improves model estimates for percent mass remaining at 25 months. Pathotroph richness improves model estmates for pmr at 37 months.

Plot the relationship between OTU richness and residuals of pmr at time25 ![](readme_files/figure-markdown_github/unnamed-chunk-41-1.png) In small stems, more fungal OTUs leads to less mass remaining after 25 months than would be expected based on wood traits alone. This relationship looks heavily influenced by 1 outlier w/ high OTU richness

Plot the relationship between Pathotroph OTU richness and residuals of pmr at time37 ![](readme_files/figure-markdown_github/unnamed-chunk-42-1.png) In small stems, more pathotrophs leads to less mass remaining after 37 months than would be expected based on wood traits alone.

############################################## 

Relationship between wood traits and community
==============================================

*Hyp (species+size-level)* Initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

![](readme_files/figure-markdown_github/unnamed-chunk-43-1.png)

    ## quartz_off_screen 
    ##                 2

Full community anova-like table

Trimmed community anova-like table

*Hyp (stem-level)* Average initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

    ## # A tibble: 6 x 14
    ## # Groups:   code [6]
    ##   code  barkthick     C    Ca density     Fe     K     Mn     N       P
    ##   <chr>     <dbl> <dbl> <dbl>   <dbl>  <dbl> <dbl>  <dbl> <dbl>   <dbl>
    ## 1 acel       2.05  51.2 2807.   0.559 2296.   651.  45.2  0.749   -6.79
    ## 2 acpa       1.68  50.7 5267.   0.700 3522.  1667.  39.6  0.575   69.1 
    ## 3 ACPA       1.68  49.4 3303.   0.700   24.0  745.   4.68 0.255   75.3 
    ## 4 alli       2.50  49.9 7652.   0.619 1942.  2928. 114.   0.345  256.  
    ## 5 ALLI       2.50  50.1 4842.   0.619   35.5  776.  59.6  0.180   44.7 
    ## 6 anba       1.38  50.1 6730.   0.505 5952.  3214. 169.   0.355   98.9 
    ## # ... with 4 more variables: waterperc <dbl>, Zn <dbl>, species <chr>,
    ## #   size <chr>

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

![](readme_files/figure-markdown_github/unnamed-chunk-46-1.png)

    ## quartz_off_screen 
    ##                 2

Full community anova-like table

Trimmed community anova-like table

############################################## 

Extra pieces
------------

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
