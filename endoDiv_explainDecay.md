Does endophyte diversity explain decay?
================
Marissa Lee
12/2/2018

Load libraries, functions, data

    ## calc_distances.R :
    ## calc_diversity.R :
    ## helper_fxns.R :
    ## load_decayData.R :
    ## load_microbeData.R :
    ## load_traitData.R :
    ## make_figs_decayPatterns.R :
    ## make_figs_endoComp_explainDecay.R :
    ## make_figs_endoDiv_explainDecay.R :
    ## make_figs_woodTraits_explainDecay.R :
    ## make_summaryTables.R :
    ## shape_analysis_dataframes.R :

    ## [1] "WARNING 62"
    ## [1] "Saprotroph"  "Symbiotroph" "Pathotroph" 
    ## [1] "WARNING 78"
    ## [1] "Pathotroph-Symbiotroph"            "Pathotroph-Saprotroph-Symbiotroph"
    ## [1] "WARNING 102"
    ## [1] "Pathotroph"                        "Pathotroph-Saprotroph-Symbiotroph"
    ## [1] "WARNING 173"
    ## [1] "Saprotroph"                        "Pathotroph-Saprotroph-Symbiotroph"
    ## [1] "WARNING 178"
    ## [1] "Pathotroph-Saprotroph-Symbiotroph" "Symbiotroph"                      
    ## [1] "WARNING 250"
    ## [1] "Pathotroph-Symbiotroph"            "Pathotroph-Saprotroph-Symbiotroph"

Diversity (and diversity of specific clades) as a predictor
===========================================================

**Note that the full community matrix was used for these analyses**

*Hyp-a (species+size-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (species+size-level)* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (species+size-level)* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers.

    ##        term       source respvar     coef signif
    ## 1 sizesmall Basidio.rich   alpha  -0.76 *   TRUE
    ## 2 sizesmall    Oomy.rich    beta     -1 *   TRUE
    ## 3 sizesmall    Oomy.rich       k 0.075 **   TRUE
    ## 4 sizesmall    Oomy.rich     t60 -0.86 **   TRUE
    ## 5 sizesmall    Oomy.rich   w.t60 -0.77 **   TRUE

The only term that is significant in these models is size class. Small stems had significantly fewer oomycete taxa.

*Hyp-a (stem-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to less mass remaining esp. at early time steps because of the selection effect for fast decayers and complementarity among taxa for decay.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to more mass remaining because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (stem-level)* Greater saprotroph and basidiomycete richness will lead to less mass remaining esp. at early time steps because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to more mass remaining because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (stem-level)* Greater pathogen and oomycete richness will lead to more mass remaining because the presence of these organisms will inhibit the establishment and activity of decayers.

    ##                 term     source respvar       coef signif
    ## 1          sizesmall  Oomy.rich  time13 -0.084 ***   TRUE
    ## 2 sizesmall:sub_rich   Richness  time25  -0.0016 *   TRUE
    ## 3 sizesmall:sub_rich  ShannonsH  time25    -0.14 *   TRUE
    ## 4           sub_rich  ShannonsH  time25    0.091 *   TRUE
    ## 5 sizesmall:sub_rich   Richness  time37  -0.0021 *   TRUE
    ## 6 sizesmall:sub_rich  ShannonsH  time37    -0.17 *   TRUE
    ## 7 sizesmall:sub_rich Patho.rich  time37   -0.035 *   TRUE
    ## 8          sizesmall  Oomy.rich   time7  -0.048 **   TRUE

OTU richness and Shannon's H are associated with percent mass remaining at time25 and time37. Pathogen OTU richness is associated with pmr at time 37. In all cases, there is a negative interaction between stem size and diversity.

Plot the relationship between OTU richness and Shannon's H and pmr at time25 and time37

    ## quartz_off_screen 
    ##                 2

[](output/figures/supplementary/richness_sapro_pmr.pdf)

Not especially convincing relationships. On the whole, diversity is positively correlated with mass remaining in large stems but negatively correlated in small stems.

Plot the relationship between pathogen OTU richness and pmr at time37

[](output/figures/supplementary/richness_patho_pmr.pdf)

Pathotroph richness is positively correlated with mass remaining in large stems, but negatively correlated in small stems.

############################################## 

Diversity plus traits as a predictor
====================================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##             term     source respvar     coef signif
    ## 1           mean Sapro.rich    beta -0.29 **   TRUE
    ## 2 sizesmall:mean Sapro.rich    beta   0.25 *   TRUE
    ## 3           mean Sapro.rich     t60   -0.1 *   TRUE

    ##                                  term     alpha     beta         k
    ## Richness.1                (Intercept)    -0.25      1.2    -0.016 
    ## Richness.2                       mean   0.0018  -0.0088   0.00011 
    ## Richness.3                  sizesmall   -0.078    -0.83    0.0027 
    ## Richness.4             sizesmall:mean  0.00061   0.0059  -1.7e-05 
    ## Richness.5                       <NA>      <NA>     <NA>      <NA>
    ## Richness.Fstat                  Fstat       1.5     0.71       0.1
    ## Richness.numdf                  numdf         3        3         3
    ## Richness.dendf                  dendf        29       29        29
    ## Richness.r.squared          r.squared      0.13     0.07      0.01
    ## ShannonsH.1               (Intercept)     0.44      1.5     -0.04 
    ## ShannonsH.2                      mean    -0.13    -0.43     0.012 
    ## ShannonsH.3                 sizesmall    -0.93     -1.4     0.071 
    ## ShannonsH.4            sizesmall:mean     0.27     0.41     -0.02 
    ## ShannonsH.5                      <NA>      <NA>     <NA>      <NA>
    ## ShannonsH.Fstat                 Fstat      0.82     0.19      0.14
    ## ShannonsH.numdf                 numdf         3        3         3
    ## ShannonsH.dendf                 dendf        29       29        29
    ## ShannonsH.r.squared         r.squared      0.08     0.02      0.01
    ## Sapro.rich.1              (Intercept)    -0.27    2.1 **    -0.06 
    ## Sapro.rich.2                     mean    0.038  -0.29 **   0.0082 
    ## Sapro.rich.3                sizesmall   -0.013     -1.7     0.012 
    ## Sapro.rich.4           sizesmall:mean  -0.0088    0.25 *  -0.0034 
    ## Sapro.rich.5                     <NA>      <NA>     <NA>      <NA>
    ## Sapro.rich.Fstat                Fstat      1.51     3.04      1.99
    ## Sapro.rich.numdf                numdf         3        3         3
    ## Sapro.rich.dendf                dendf        29       29        29
    ## Sapro.rich.r.squared        r.squared      0.13     0.24      0.17
    ## Basidio.rich.1            (Intercept)    0.017     0.55    -0.021 
    ## Basidio.rich.2                   mean  -0.0012   -0.039    0.0014 
    ## Basidio.rich.3              sizesmall    -0.25    -0.69     0.043 
    ## Basidio.rich.4         sizesmall:mean    0.017    0.048    -0.003 
    ## Basidio.rich.5                   <NA>      <NA>     <NA>      <NA>
    ## Basidio.rich.Fstat              Fstat      0.63     0.27      0.35
    ## Basidio.rich.numdf              numdf         3        3         3
    ## Basidio.rich.dendf              dendf        29       29        29
    ## Basidio.rich.r.squared      r.squared      0.06     0.03      0.04
    ## Patho.rich.1              (Intercept)    -0.13      1.1    -0.029 
    ## Patho.rich.2                     mean    0.017    -0.14    0.0038 
    ## Patho.rich.3                sizesmall    0.016     -1.1     0.019 
    ## Patho.rich.4           sizesmall:mean  -0.0046     0.15   -0.0027 
    ## Patho.rich.5                     <NA>      <NA>     <NA>      <NA>
    ## Patho.rich.Fstat                Fstat      0.35     0.86      0.31
    ## Patho.rich.numdf                numdf         3        3         3
    ## Patho.rich.dendf                dendf        29       29        29
    ## Patho.rich.r.squared        r.squared      0.04     0.08      0.03
    ## Oomy.rich.1               (Intercept) -2.2e-17    4e-17    -5e-18 
    ## Oomy.rich.2                 sizesmall  3.5e-17   -4e-17     5e-18 
    ## Oomy.rich.3                      <NA>      <NA>     <NA>      <NA>
    ## Oomy.rich.Fstat                 Fstat         0        0         0
    ## Oomy.rich.numdf                 numdf         1        1         1
    ## Oomy.rich.dendf                 dendf        31       31        31
    ## Oomy.rich.r.squared         r.squared         0        0         0
    ##                            ne.r2      t60      w.r2     w.t60       source
    ## Richness.1               -0.031     0.26    -0.047     0.058      Richness
    ## Richness.2              0.00029  -0.0019   0.00033  -0.00041      Richness
    ## Richness.3              -0.0038    0.022     0.027    -0.019      Richness
    ## Richness.4             -7.3e-05   -2e-04  -0.00019   0.00013      Richness
    ## Richness.5                  <NA>     <NA>      <NA>      <NA>     Richness
    ## Richness.Fstat              0.54     0.63      0.23      0.02     Richness
    ## Richness.numdf                 3        3         3         3     Richness
    ## Richness.dendf                29       29        29        29     Richness
    ## Richness.r.squared          0.05     0.06      0.02         0     Richness
    ## ShannonsH.1               -0.19     0.35    -0.089      0.34     ShannonsH
    ## ShannonsH.2               0.058     -0.1     0.026    -0.099     ShannonsH
    ## ShannonsH.3                0.12    -0.35     0.045      -0.7     ShannonsH
    ## ShannonsH.4              -0.041      0.1    -0.013       0.2     ShannonsH
    ## ShannonsH.5                 <NA>     <NA>      <NA>      <NA>    ShannonsH
    ## ShannonsH.Fstat             1.33     0.07      0.23      0.36    ShannonsH
    ## ShannonsH.numdf                3        3         3         3    ShannonsH
    ## ShannonsH.dendf               29       29        29        29    ShannonsH
    ## ShannonsH.r.squared         0.12     0.01      0.02      0.04    ShannonsH
    ## Sapro.rich.1            -0.0023    0.73 *   -0.018      0.49    Sapro.rich
    ## Sapro.rich.2             0.0016    -0.1 *   0.0024    -0.068    Sapro.rich
    ## Sapro.rich.3             0.0071    -0.34     0.013      -0.3    Sapro.rich
    ## Sapro.rich.4            -0.0026    0.061   -0.0019     0.049    Sapro.rich
    ## Sapro.rich.5                <NA>     <NA>      <NA>      <NA>   Sapro.rich
    ## Sapro.rich.Fstat            0.26     3.74      0.04      1.75   Sapro.rich
    ## Sapro.rich.numdf               3        3         3         3   Sapro.rich
    ## Sapro.rich.dendf              29       29        29        29   Sapro.rich
    ## Sapro.rich.r.squared        0.03     0.28         0      0.15   Sapro.rich
    ## Basidio.rich.1           -0.037     0.18    -0.028      0.15  Basidio.rich
    ## Basidio.rich.2           0.0032   -0.012     0.002     -0.01  Basidio.rich
    ## Basidio.rich.3            0.026    -0.11     0.041     -0.27  Basidio.rich
    ## Basidio.rich.4          -0.0029   0.0082   -0.0028     0.019  Basidio.rich
    ## Basidio.rich.5              <NA>     <NA>      <NA>      <NA> Basidio.rich
    ## Basidio.rich.Fstat          0.68      0.2      0.17      0.32 Basidio.rich
    ## Basidio.rich.numdf             3        3         3         3 Basidio.rich
    ## Basidio.rich.dendf            29       29        29        29 Basidio.rich
    ## Basidio.rich.r.squared      0.07     0.02      0.02      0.03 Basidio.rich
    ## Patho.rich.1             -0.069     0.19    -0.067      0.13    Patho.rich
    ## Patho.rich.2               0.01   -0.024    0.0087    -0.017    Patho.rich
    ## Patho.rich.3              0.058    -0.07     0.069     -0.14    Patho.rich
    ## Patho.rich.4            -0.0096    0.012    -0.009     0.018    Patho.rich
    ## Patho.rich.5                <NA>     <NA>      <NA>      <NA>   Patho.rich
    ## Patho.rich.Fstat            1.35     0.29      0.66      0.11   Patho.rich
    ## Patho.rich.numdf               3        3         3         3   Patho.rich
    ## Patho.rich.dendf              29       29        29        29   Patho.rich
    ## Patho.rich.r.squared        0.12     0.03      0.06      0.01   Patho.rich
    ## Oomy.rich.1              0.0093  4.1e-17     1e-17   3.4e-17     Oomy.rich
    ## Oomy.rich.2              -0.015   -8e-17    -5e-18    -4e-17     Oomy.rich
    ## Oomy.rich.3                 <NA>     <NA>      <NA>      <NA>    Oomy.rich
    ## Oomy.rich.Fstat             0.71        0         0         0    Oomy.rich
    ## Oomy.rich.numdf                1        1         1         1    Oomy.rich
    ## Oomy.rich.dendf               31       31        31        31    Oomy.rich
    ## Oomy.rich.r.squared         0.02        0         0         0    Oomy.rich

Beyond wood trait data, saprotroph richness improves models for beta and t60.

Plot the relationship between Saprotroph richness and beta, t60

[](output/figures/supplementary/richness_sapro_decayResids.pdf)

Higher saprotroph richness leads to smaller beta and t70 values (i.e. faster mass loss) than would be expected based on wood traits alone, especially in large stems.

*Hyp (stem-level)* After accounting for variation in decay due to wood traits, initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in percent mass loss, esp. at early time points.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##                 term       source respvar      coef signif
    ## 1 sizesmall:sub_rich     Richness  time25 -0.0012 *   TRUE
    ## 2          sizesmall Basidio.rich  time25    0.11 *   TRUE
    ## 3 sizesmall:sub_rich Basidio.rich  time25 -0.0099 *   TRUE
    ## 4 sizesmall:sub_rich   Patho.rich  time37  -0.024 *   TRUE

Beyond wood trait data, overall and basidio OTU richness improve model estimates for percent mass remaining at 25 months. Pathotroph richness improves model estmates for pmr at 37 months.

Plot the relationship between OTU richness, basidio richness and residuals of pmr at time25

[](output/figures/supplementary/richness_PMRResids.pdf)

In small stems, more OTU richness is associated with less mass remaining after 25 months than would be expected based on wood traits alone. This relationship looks heavily influenced by 1 outlier w/ high OTU (and basidio) richness

Plot the relationship between Pathotroph OTU richness and residuals of pmr at time37

[](output/figures/supplementary/richness_patho_time37Resids.pdf)

In small stems, more pathotroph OTUs is associated with less mass remaining after 37 months than would be expected based on wood traits alone.
