Do wood traits explain decay?
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
    ## make_figs_woodTraits_explainDecay.R :
    ## make_summaryTables.R :
    ## shape_analysis_dataframes.R :

We expect initial wood traits will explain varitation in species+size decay rate (k and t70), species+size lagginess (alpha), and stem-level percent mass remaining at 7, 13, 25, and 37 months of decay. Specifically, we expect samples with (a) high water percent, (b) low density and total C, (c) high macro and micro nutrients, and (d) thicker bark (potential mech: limiting microbial colonization) to have faster decay and less lagginess.

*Hyp (species+size-level)* Species+size-level initial wood traits will predict variation decay rates and lagginess.
-------------------------------------------------------------------------------------------------------------------

![](output/figures/maintext/traits_explain_decayparams.pdf)

Variation in w.t70 and beta (shape param) were best explained by traits (high r2); whereas traits explained much less variation in alpha (scale) and w.r2.

-   alpha: Samples that decay fast have low density. Other fast-decay attributes include: thin bark, low C, high N, and belonging to the large size class.
-   beta: Samples with a stronger S-shaped decay trajectory have low N and high density. Other attributes include: thick bark, high C, and belonging to the large size class.
-   w.r2: Samples with tightly-fitting decay models have high N, low density, and belong to the large size class.
-   w.t70: Samples with long wait times to 70% mass remaining have low N and belong to the large size class. Other attributes include: thick bark, high C, and low water content.

*Hyp (stem-level)* Stem-level initial wood traits will predict variation in percent mass loss at each time step.
----------------------------------------------------------------------------------------------------------------

First, we need to decide what trait data (and samples) to include in this analysis since we don't have full coverage of stem-level trait data. Density and bark thickness were only measured on small sized stems. If there is not be very much within-species variation in these traits that contribute to variation in percent mass loss than we can justify including species-level estimates of these traits in the stem-level model.

Plot the small-sized stem-level measures of density and barkthick

    ## quartz_off_screen 
    ##                 2

![](output/figures/supplementary/variation_densityNbarkthick.pdf)

A few species have a relatively large amount of within-species variation in density (i.e. alli) and barkthickness (i.e. leer and cota).

Compare model fits (r2) using stem and species-level data to identify how much information about percent mass remaining is lost by using species-level estimates.

    ##                       mod       pval source      meas
    ## time7.1   fitted(M1) ~ M2 0.22763856  time7   density
    ## time7.2   fitted(M2) ~ M1 0.01454204  time7   density
    ## time13.1  fitted(M1) ~ M2 0.52899319 time13   density
    ## time13.2  fitted(M2) ~ M1 0.18583285 time13   density
    ## time25.1  fitted(M1) ~ M2 0.60081755 time25   density
    ## time25.2  fitted(M2) ~ M1 0.76919241 time25   density
    ## time37.1  fitted(M1) ~ M2 0.69863046 time37   density
    ## time37.2  fitted(M2) ~ M1 0.75422808 time37   density
    ## time59.1  fitted(M1) ~ M2 0.25212110 time59   density
    ## time59.2  fitted(M2) ~ M1 0.01881453 time59   density
    ## time7.11  fitted(M1) ~ M2 0.70557781  time7 barkthick
    ## time7.21  fitted(M2) ~ M1 0.90907489  time7 barkthick
    ## time13.11 fitted(M1) ~ M2 0.04423876 time13 barkthick
    ## time13.21 fitted(M2) ~ M1 0.51560418 time13 barkthick
    ## time25.11 fitted(M1) ~ M2 0.11620024 time25 barkthick
    ## time25.21 fitted(M2) ~ M1 0.61127201 time25 barkthick
    ## time37.11 fitted(M1) ~ M2 0.38828743 time37 barkthick
    ## time37.21 fitted(M2) ~ M1 0.75244608 time37 barkthick
    ## time59.11 fitted(M1) ~ M2 0.19239411 time59 barkthick
    ## time59.21 fitted(M2) ~ M1 0.89110070 time59 barkthick

-   For density -- Stem-level data (relative to code-level data) improves estimates of pmr after 7 and 59 months.
-   For barkthickness -- Stem-level data (relative to code-level data) does not improve estimates of pmr at any time point.

Compile a "stem-level" dataframe with (a) stem-level percent mass remaining values, (b) stem-level traits including waterperc and chemistry along, and (c) small species-level density and bark thickness data. Use model selection to determine which traits to include.

![](output/figures/maintext/traits_explain_pmr.pdf)

The explanatory power of wood traits were lowest for the first decay time point (7 months). Samples had more mass remaining...

-   All time points: Low water content
-   Early-to-mid time points: Low N and belonged to the small size class
-   Mid time points: High P and low Zn
-   Late time point: High C and thin bark
