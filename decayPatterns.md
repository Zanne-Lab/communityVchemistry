How do wood species and sizes vary in decay?
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

Plot percent mass remaining (pmr) for each sample over time ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-2-1.png)

Compare the negative exp. vs weibull models by plotting time to 70% mass remaining (t70) for each species+size ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-3-1.png)

Another view... Figure S2. Comparison of negative exp and weibull models using t70

![](output/figures/supplementary/compare_t60.pdf)

Figure S3. Comparison of negative exp and weibull models using AIC

![](output/figures/supplementary/compare_aic.pdf)

Figure 1. Wood species x decay params (weibull model)

    ## quartz_off_screen 
    ##                 2

![](output/figures/maintext/decayfits.pdf)

Figure 2. Time x percent mass remaining

![](output/figures/maintext/pmr_byStem.pdf)
