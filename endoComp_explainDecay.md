Does endophyte composition explain decay?
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

Figure S1. Sample-effort curves

Filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 132 of 3021 OTUs"

*Hyp (species+size-level)* Species+size-level (average) initial microbial community composition will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

None of the components are significant.

*Hyp (stem-level)* Stem-level initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
---------------------------------------------------------------------------------------------------------------------------------------------------------------

Comp01 (of the non-trimmed community) is a significant predictor of percent mass remaining at 37 months.

Plot the distribution of WA-PLS scores

[](output/figures/supplementary/wapls_score_time37.pdf)

Who is in the top and bottom 1%?

    ##        quant kingdom        phylum                           species
    ## 1  bottom 1%   Fungi    Ascomycota            Neosetophoma_samararum
    ## 2  bottom 1%   Fungi Basidiomycota          Pisolithus_croceorrhizus
    ## 3  bottom 1%   Fungi    Ascomycota Neotrimmatostroma_paraexcentricum
    ## 4  bottom 1%   Fungi    Ascomycota               Phaeoacremonium_sp3
    ## 5  bottom 1%   Fungi    Ascomycota                 Strigula_nitidula
    ## 6  bottom 1%   Fungi    Ascomycota                     Acremonium_sp
    ## 7  bottom 1%   Fungi    Ascomycota           Celerioriella_prunicola
    ## 8     top 1%   Fungi    Ascomycota             Debaryomyces_hansenii
    ## 9     top 1%   Fungi Basidiomycota                    Naganishia_sp3
    ## 10    top 1%   Fungi    Ascomycota               Marcelaria_cumingii
    ## 11    top 1%   Fungi Basidiomycota              Septobasidium_burtii
    ##             Trophic.Mode
    ## 1             Saprotroph
    ## 2            Symbiotroph
    ## 3  Pathotroph-Saprotroph
    ## 4           unclassified
    ## 5            Symbiotroph
    ## 6           unclassified
    ## 7           unclassified
    ## 8           unclassified
    ## 9           unclassified
    ## 10          unclassified
    ## 11          unclassified
    ##                                                  Guild
    ## 1                                 Undefined Saprotroph
    ## 2                                      Ectomycorrhizal
    ## 3  Animal Pathogen-Plant Pathogen-Undefined Saprotroph
    ## 4                                         unclassified
    ## 5                                           Lichenized
    ## 6                                         unclassified
    ## 7                                         unclassified
    ## 8                                         unclassified
    ## 9                                         unclassified
    ## 10                                        unclassified
    ## 11                                        unclassified

Many of the bottom 1% OTUs are classified as saprotrophs. That makes sense since low WA-PLS scores indicate an association with high mass loss (i.e. less mass remaining) at time37.

But saprotrophs are also found at many points along the gradient...

    ## quartz_off_screen 
    ##                 2

[](output/figures/supplementary/wapls_score_time37_otuCats.pdf)

Is this because there is an underlying signature of wood traits on the initial microbial community that is driving the relationship between the community and the mass remaining after 37 months? The next analysis ("Community+traits" as predictor) will test this formally. Just out of curiousity, I'd like to pull in OTU "niche" info from the boral analysis to see if there's a relationship between OTU WA-PLS scores and wood trait coeffient estimates.

Reminder of which wood traits were included in the best model to explain pmr at time37... [](output/figures/maintext/traits_explain_pmr.pdf) More water leads to less mass remaining; more P leads to more mass remaining

Plot OTU wood trait estimates (from boral) versus signif WA-PLS score.

[](output/figures/supplementary/wapls_score_time37_boral.pdf) There's a weak negative relationship between an OTU's WA-PLS score and waterperc coefficient (slope=-2.3, p=.03), suggesting that OTUs that "prefer" high-water niche space are associated with less mass remaining at time37.

########################################## 

Community+traits as a predictor
===============================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Community data doesn't improve our understanding of decay rates (k, t70) or variation in decay rate (ne.r2) beyond what is known from the trait data.

*Hyp (stem-level)* After accounting for variation in decay due to wood traits (no models with density, includes small-species level bark thickness), stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Community data doesn't improve our understanding of mass loss (pmr after 7, 13, 25, 37, and 59 months) beyond what is known from the trait data.

Relationship between wood traits and community
==============================================

*Hyp (species+size-level)* Initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

Use ordistep to find the best combination of traits to constrain variation in the endophyte community

    ## [1] "Keep 132 of 3021 OTUs"

    ## quartz_off_screen 
    ##                 2

[](output/figures/supplementary/dbRDA_code.pdf)

Anova-like tables

    ## [1] "Keep 132 of 3021 OTUs"

    ## $an.nt.code
    ##                    term    df Fval  pval
    ## 1                  size  1.00 1.40 0.012
    ## 2             barkthick  1.00 1.60 0.003
    ## 3             waterperc  1.00 1.41 0.013
    ## 4              Residual 29.00   NA    NA
    ## 5                  <NA>    NA   NA    NA
    ## 6   Constrained inertia  1.91   NA    NA
    ## 7 Unconstrained inertia 11.16   NA    NA
    ## 
    ## $an.t.code
    ##                    term    df Fval  pval
    ## 1                  size  1.00 1.67 0.009
    ## 2             waterperc  1.00 1.58 0.031
    ## 3             barkthick  1.00 1.54 0.035
    ## 4              Residual 29.00   NA    NA
    ## 5                  <NA>    NA   NA    NA
    ## 6   Constrained inertia  1.79   NA    NA
    ## 7 Unconstrained inertia  9.22   NA    NA

*Hyp (stem-level)* Average initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

    ## [1] "Keep 132 of 3021 OTUs"

    ## quartz_off_screen 
    ##                 2

Anova-like tables

    ## [1] "Keep 132 of 3021 OTUs"

    ## $an.nt.stem
    ##                     term    df Fval  pval
    ## 1        barkthick_smspp  1.00 2.33 0.001
    ## 2              waterperc  1.00 1.99 0.001
    ## 3                   size  1.00 2.13 0.001
    ## 4                     Ca  1.00 1.91 0.001
    ## 5          density_smspp  1.00 1.97 0.001
    ## 6                      C  1.00 1.96 0.001
    ## 7                      N  1.00 2.02 0.001
    ## 8                     Mn  1.00 1.62 0.001
    ## 9                      K  1.00 1.60 0.003
    ## 10              Residual 47.00   NA    NA
    ## 11                  <NA>    NA   NA    NA
    ## 12   Constrained inertia  7.55   NA    NA
    ## 13 Unconstrained inertia 17.37   NA    NA
    ## 
    ## $an.t.stem
    ##                     term    df Fval  pval
    ## 1              waterperc  1.00 2.37 0.001
    ## 2                   size  1.00 2.60 0.001
    ## 3                      C  1.00 2.70 0.001
    ## 4        barkthick_smspp  1.00 2.64 0.001
    ## 5                      N  1.00 2.02 0.003
    ## 6                      K  1.00 1.59 0.024
    ## 7                     Ca  1.00 1.93 0.003
    ## 8          density_smspp  1.00 1.68 0.011
    ## 9                     Mn  1.00 1.61 0.012
    ## 10              Residual 47.00   NA    NA
    ## 11                  <NA>    NA   NA    NA
    ## 12   Constrained inertia  7.18   NA    NA
    ## 13 Unconstrained inertia 14.86   NA    NA
