Does chemistry or community better predict mass loss?
================
10/23/2017

    ## Warning in strptime(x, fmt, tz = "GMT"): unknown timezone 'default/America/
    ## New_York'

### Load microbial community data

### Load wood trait data

### Load decay data

``` r
# load mass loss data and calculate percent mass remaining (pmr)
#initial_mass <- read_in_initial_mass() #uncomment if the data changes
#harvest_mass<-LoadHarvestFiles()
#pmr<-Calc_massRemaining(initial_mass, harvest_mass)
#write.csv(pmr,"derived_data/pmr.csv")
pmr<-read.csv("derived_data/pmr.csv", row.names=1)

# average pmr for each stem and timepoint
pmr_byStem<-AvgPMR_byStem(pmr)

# calculate decay trajectory fits for each species+size
#decayfits <- fit_all_curves(pmr, stemSamples) #this recalculates all the curve fits, uncomment if the data changes
#write_csv(decayfits,"derived_data/decayfits.csv")
decayfits <- read_csv("derived_data/decayfits.csv")

# ggplot(decayfits,aes(x=t70,y=w.t70,col=size))+
#   geom_point()+
#   labs(x="Time to 30% mass loss (negative exponential)", 
#        y="Time to 30% mass loss (Weibull)")+
#   geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()

rm(pmr)
```

Check for missing stem-level data

``` r
traits.stem %>%
  filter(!codeStem %in% stemSamples$codeStem) #there are 2 stems that are in the traits df but are missing from the stem lookup table...
```

    ## # A tibble: 2 x 14
    ##   codeStem  code  size waterperc density barkthick       P       K
    ##      <chr> <chr> <chr>     <dbl>   <dbl>     <dbl>   <dbl>   <dbl>
    ## 1    acpa2  acpa small        NA      NA        NA  78.423 1400.18
    ## 2    lepa4  lepa small        NA      NA        NA 147.198  748.14
    ## # ... with 6 more variables: Ca <dbl>, Mn <dbl>, Fe <dbl>, Zn <dbl>,
    ## #   N <dbl>, C <dbl>

``` r
seqSamples %>%
  filter(!codeStem %in% stemSamples$codeStem) %>%
  separate(seq_sampName, into=c("drop","seq.stem"), 4, remove=FALSE) %>%
  filter(grepl("[1-9]", seq.stem)) #these same 2 stems are in the sequence samples but are missing from the stem lookup table...
```

    ## # A tibble: 2 x 7
    ##   seq_sampName  drop seq.stem  code codeStem species  size
    ##          <chr> <chr>    <chr> <chr>    <chr>   <chr> <chr>
    ## 1        lepa4  lepa        4  lepa     <NA>    lepa small
    ## 2        acpa2  acpa        2  acpa     <NA>    acpa small

``` r
pmr_byStem %>%
  filter(!codeStem %in% stemSamples$codeStem) # no stem samples are missing from the decay data
```

    ## # A tibble: 0 x 8
    ## # Groups:   codeStem [0]
    ## # ... with 8 variables: codeStem <chr>, code <chr>, Stem <chr>,
    ## #   time0 <dbl>, time7 <dbl>, time13 <dbl>, time25 <dbl>, time37 <dbl>

For some reason there are 2 unique codeStem ids that are found in the trait data (xrf sample names) and the sequence data, but not in the stemSamples data (deployment sample names). These codeStem ids are not found in the percent mass loss data. Because the main goal is to analyze decay responses, I'm going to continue to leave these codeStems out of the stemSamples dataframe. Is it possible that stem id numbers got switched? Something to follow-up on.

Figure 1. Wood species x decay params

``` r
#make plot of decay params by wood species
decayfits.p<-left_join(decayfits, stemSamples)

# t70
p.t70<-ggplot(decayfits.p, aes(x=reorder(Binomial, -t70), y=t70, shape=size)) + 
  geom_point() +
  #geom_errorbar(aes(ymin=t70.lower, ymax=t70.upper)) +
  coord_flip() +
  xlab("Wood species") + ylab("Time to 70% mass remaining (years)") + theme_bw() +
  scale_shape_manual(values=c('large'=19, 'small'=1)) + 
  guides(shape=F) + 
  theme(plot.margin = unit(c(0,1,0,0), "lines"))

# k
p.k<-ggplot(decayfits.p, aes(x=reorder(Binomial, -t70), y=k, shape=size)) + 
  geom_point() +
  geom_errorbar(aes(ymin=k.lower, ymax=k.upper)) +
  coord_flip() +
  xlab("Wood species") + ylab("k (year^-1)") + theme_bw() +
  scale_shape_manual(values=c('large'=19, 'small'=1)) + 
  guides(shape=F) +
  theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank(),
            plot.margin = unit(c(0,1,0,0), "lines"),
            plot.background = element_blank())

# ne.r2
p.r2<-ggplot(decayfits.p, aes(x=reorder(Binomial, -t70), y=ne.r2, shape=size)) + 
  geom_point() +
  coord_flip() +
  xlab("Wood species") + ylab("Neg. expon. R2") + theme_bw() +
  scale_shape_manual(values=c('large'=19, 'small'=1)) + 
  guides(shape=F) +
  theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank(),
            plot.margin = unit(c(0,1,0,0), "lines"),
            plot.background = element_blank())

pdf("output/decayfits.pdf", width=10, height=6)
grid.newpage()
grid.draw(cbind(ggplotGrob(p.t70), ggplotGrob(p.k), ggplotGrob(p.r2), size = "last"))
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
rm(p.t70, p.k, p.r2)
```

Figure 2. Time x percent mass remaining by wood species

``` r
#make a long version of pmr.byStem.df.w
pmr_byStem %>%
  gather(key="time",value="mean.pmr", 4:8, na.rm=TRUE) %>%
  separate(time, into=c("drop","timeNum"), sep="time") %>%
  mutate(timeNum=as.numeric(timeNum)) %>%
  mutate(species=tolower(code)) %>%
  mutate(size = case_when(code == species ~ "small", 
                          code != species ~ "large")) %>%
  select(-drop) -> pmr.byStem.df

#order the wood species by t70 to match species order in previous figure
woodSp.o<-levels(reorder(decayfits.p$species, -decayfits.p$t70))
pmr.byStem.df$species<-factor(pmr.byStem.df$species, levels=rev(woodSp.o))

p.pmrBystem<-ggplot(pmr.byStem.df, aes(x=timeNum, y=mean.pmr*100, group=codeStem, linetype=size)) + 
  geom_line() + facet_wrap(~species) +
  xlab("Time (years)") + ylab("Mass remaining (%)") + theme_bw() +
  scale_linetype_manual(values=c(2,1)) + guides(linetype=FALSE)

pdf("output/pmr_byStem.pdf", width=6, height=6)
p.pmrBystem
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
rm(decayfits.p, woodSp.o, p.pmrBystem, pmr.byStem.df)
```

########################################## 

Wood traits as a predictor
==========================

We expect initial wood traits will explain varitation in species+size decay rate (k and t70), species+size lagginess (alpha), and stem-level percent mass remaining at 7, 13, 25, and 37 months of decay. Specifically, we expect samples with (a) high water percent, (b) low density and total C, (c) high macro and micro nutrients, and (d) thicker bark (potential mech: limiting microbial colonization) to have faster decay and less lagginess.

*Hyp (species+size-level)* Species+size-level initial wood traits will predict variation decay rates and lagginess.
-------------------------------------------------------------------------------------------------------------------

*Hyp (stem-level)* Stem-level initial wood traits will predict variation in percent mass loss at each time step.
----------------------------------------------------------------------------------------------------------------

First, we need to decide what trait data (and samples) to include in this analysis since we don't have full coverage of stem-level trait data. Density and bark thickness were only measured on small sized stems. If there is not be very much within-species variation in these traits that contribute to variation in percent mass loss than we can justify including species-level estimates of these traits in the stem-level model.

Plot the small-sized stem-level measures of density and barkthick

    ## quartz_off_screen 
    ##                 2

Compare model fits (r2) using stem and species-level data to identify how much information about percent mass remaining is lost by using species-level estimates...For density, it looks like stem-level data improves model fit a tiny bit for early percent mass remaining time points (after 7 and 13 months) but not later time points. For barkthickness, fits are about the same.

    ##   respvars density_code_r2 density_stem_r2 barkthick_code_r2
    ## 1    time7            0.05            0.08              0.02
    ## 2   time13            0.03            0.07              0.00
    ## 3   time25            0.03            0.03              0.03
    ## 4   time37            0.07            0.07              0.11
    ##   barkthick_stem_r2
    ## 1              0.01
    ## 2              0.00
    ## 3              0.01
    ## 4              0.10

Compile a "stem-level" dataframe with (a) stem-level percent mass remaining values, (b) stem-level traits including waterperc and chemistry along, and (c) small species-level density and bark thickness data.

########################################## 

Community as a predictor
========================

Filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 150 of 6128 OTUs"

*Hyp (species+size-level)* Species+size-level (average) initial microbial community composition will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

*Hyp (stem-level)* Stem-level initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
---------------------------------------------------------------------------------------------------------------------------------------------------------------

Comp01 (of the non-trimmed community) is a significant predictor of percent mass remaining at 37 months. This is likely an underlying signature of wood traits on the initial microbial community that is driving the relationship between the community and the mass remaining after 37 months. Check this out by plotting OTU trait-associated coefs (from boral) versus component coef estimate.

Plot of OTU functional and phylogenetic categories versus signif WA-PLS score. Remember that higher WA-PLS scores are associated more more mass remaingin after 37 months.

    ## quartz_off_screen 
    ##                 2

Plot OTU wood trait estimates (from boral) versus signif WA-PLS score.

    ##            X            term              time7               time13
    ## 1          1     (Intercept) 2.07 +/- 0.486 ***   1.43 +/- 0.174 ***
    ## 2          2       sizesmall               <NA> -0.172 +/- 0.046 ***
    ## 3          3       waterperc -0.005 +/- 0.002 * -0.014 +/- 0.003 ***
    ## 4          4   density_smspp               <NA>                 <NA>
    ## 5          5 barkthick_smspp               <NA>                 <NA>
    ## 6          6               C -0.019 +/- 0.009 *                 <NA>
    ## 7          7               N               <NA>  -0.424 +/- 0.142 **
    ## 8          8               P               <NA>        0.001 +/- 0 *
    ## 9          9              Mn               <NA>        0.001 +/- 0 *
    ## 10        10            <NA>               <NA>                 <NA>
    ## 11     Fstat           Fstat                3.7                 5.36
    ## 12     numdf           numdf                  3                    5
    ## 13     dendf           dendf                 50                   47
    ## 14 r.squared       r.squared               0.18                 0.36
    ##                  time25               time37
    ## 1   1.195 +/- 0.118 ***    -0.437 +/- 0.904 
    ## 2    -0.072 +/- 0.035 *  -0.169 +/- 0.048 **
    ## 3  -0.013 +/- 0.002 *** -0.021 +/- 0.004 ***
    ## 4                  <NA>     0.598 +/- 0.436 
    ## 5                  <NA>  -0.075 +/- 0.026 **
    ## 6                  <NA>    0.034 +/- 0.016 *
    ## 7                  <NA>                 <NA>
    ## 8         0.001 +/- 0 *       0.001 +/- 0 **
    ## 9                  <NA>        0.001 +/- 0 *
    ## 10                 <NA>                 <NA>
    ## 11                11.15                10.25
    ## 12                    3                    8
    ## 13                   48                   42
    ## 14                 0.41                 0.66

    ## quartz_off_screen 
    ##                 2

When looking at these barkthickness OTU coefs, remember that the boral model was fit based on small species-level bark thickness data...The asco response to barkthickness seems to be driving the WA-PLS association.

########################################## 

Community+traits as a predictor
===============================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (t70, k), and lagginess (alpha).
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

*Hyp (stem-level)* After accounting for variation in decay due to wood traits (no models with density, includes small-species level bark thickness), stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

########################################## 

Diversity (and diversity of specific clades) as a predictor
===========================================================

**Note that the full community matrix was used for these analyses**

*Hyp-a (species+size-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (species+size-level)* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (species+size-level)* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers.

    ##        variable valueType                 k             ne.r2
    ## 1      Richness  estim.se          0 +/- 0           0 +/- 0 
    ## 2      Richness     numdf                 2                 2
    ## 3      Richness     dendf                30                30
    ## 4      Richness r.squared              0.14                 0
    ## 5     ShannonsH  estim.se -0.029 +/- 0.032  -0.018 +/- 0.032 
    ## 6     ShannonsH     numdf                 2                 2
    ## 7     ShannonsH     dendf                30                30
    ## 8     ShannonsH r.squared              0.17              0.01
    ## 9    Sapro.rich  estim.se  0.002 +/- 0.003       0 +/- 0.003 
    ## 10   Sapro.rich     numdf                 2                 2
    ## 11   Sapro.rich     dendf                30                30
    ## 12   Sapro.rich r.squared              0.15                 0
    ## 13 Basidio.rich  estim.se -0.001 +/- 0.002       0 +/- 0.002 
    ## 14 Basidio.rich     numdf                 2                 2
    ## 15 Basidio.rich     dendf                30                30
    ## 16 Basidio.rich r.squared              0.17                 0
    ## 17   Patho.rich  estim.se -0.001 +/- 0.002   0.001 +/- 0.002 
    ## 18   Patho.rich     numdf                 2                 2
    ## 19   Patho.rich     dendf                30                30
    ## 20   Patho.rich r.squared              0.15              0.01
    ## 21    Oomy.rich  estim.se  -0.004 +/- 0.04  -0.034 +/- 0.039 
    ## 22    Oomy.rich     numdf                 2                 2
    ## 23    Oomy.rich     dendf                30                30
    ## 24    Oomy.rich r.squared              0.14              0.03
    ##                  t70
    ## 1   0.001 +/- 0.002 
    ## 2                  2
    ## 3                 30
    ## 4               0.17
    ## 5   0.215 +/- 0.173 
    ## 6                  2
    ## 7                 30
    ## 8               0.21
    ## 9  -0.004 +/- 0.018 
    ## 10                 2
    ## 11                30
    ## 12              0.17
    ## 13  0.008 +/- 0.008 
    ## 14                 2
    ## 15                30
    ## 16              0.19
    ## 17  0.017 +/- 0.013 
    ## 18                 2
    ## 19                30
    ## 20              0.21
    ## 21 -0.038 +/- 0.216 
    ## 22                 2
    ## 23                30
    ## 24              0.17

*Hyp-a (stem-level)* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to less mass remaining esp. at early time steps because of the selection effect for fast decayers and complementarity among taxa for decay.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hyp-Alt: Greater microbial diversity will lead to more mass remaining because taxa will be allocating more of their resources to combat one another. \#\# *Hyp-b (stem-level)* Greater saprotroph and basidiomycete richness will lead to less mass remaining esp. at early time steps because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to more mass remaining because decayers will be allocating more of their resources to combat one another. \#\# *Hyp-c (stem-level)* Greater pathogen and oomycete richness will lead to more mass remaining because the presence of these organisms will inhibit the establishment and activity of decayers.

############################################## 

Diversity plus traits as a predictor
====================================

*Hyp (species+size-level)* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##        variable valueType                 k             ne.r2
    ## 1      Richness  estim.se          0 +/- 0           0 +/- 0 
    ## 2      Richness     numdf                 2                 2
    ## 3      Richness     dendf                29                29
    ## 4      Richness r.squared              0.01              0.02
    ## 5     ShannonsH  estim.se -0.007 +/- 0.018  -0.009 +/- 0.021 
    ## 6     ShannonsH     numdf                 2                 2
    ## 7     ShannonsH     dendf                29                29
    ## 8     ShannonsH r.squared              0.01              0.01
    ## 9    Sapro.rich  estim.se 0.004 +/- 0.002 *  0.002 +/- 0.002 
    ## 10   Sapro.rich     numdf                 2                 2
    ## 11   Sapro.rich     dendf                29                29
    ## 12   Sapro.rich r.squared               0.2              0.02
    ## 13 Basidio.rich  estim.se -0.001 +/- 0.001       0 +/- 0.001 
    ## 14 Basidio.rich     numdf                 2                 2
    ## 15 Basidio.rich     dendf                29                29
    ## 16 Basidio.rich r.squared              0.08              0.01
    ## 17   Patho.rich  estim.se      0 +/- 0.001   0.001 +/- 0.002 
    ## 18   Patho.rich     numdf                 2                 2
    ## 19   Patho.rich     dendf                29                29
    ## 20   Patho.rich r.squared                 0              0.01
    ## 21    Oomy.rich  estim.se -0.011 +/- 0.022  -0.029 +/- 0.025 
    ## 22    Oomy.rich     numdf                 2                 2
    ## 23    Oomy.rich     dendf                29                29
    ## 24    Oomy.rich r.squared              0.01              0.04
    ##                   t70
    ## 1        0 +/- 0.001 
    ## 2                   2
    ## 3                  29
    ## 4                   0
    ## 5   -0.001 +/- 0.101 
    ## 6                   2
    ## 7                  29
    ## 8                   0
    ## 9  -0.022 +/- 0.009 *
    ## 10                  2
    ## 11                 29
    ## 12               0.17
    ## 13   0.003 +/- 0.005 
    ## 14                  2
    ## 15                 29
    ## 16               0.02
    ## 17   0.001 +/- 0.008 
    ## 18                  2
    ## 19                 29
    ## 20                  0
    ## 21    0.07 +/- 0.122 
    ## 22                  2
    ## 23                 29
    ## 24               0.01

*Hyp (stem-level)* After accounting for variation in decay due to wood traits, initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in percent mass loss, esp. at early time points.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ##        variable valueType            time13            time25
    ## 1      Richness  estim.se          0 +/- 0           0 +/- 0 
    ## 2      Richness     numdf                 2                 2
    ## 3      Richness     dendf                50                49
    ## 4      Richness r.squared              0.01              0.02
    ## 5     ShannonsH  estim.se -0.002 +/- 0.031  -0.005 +/- 0.029 
    ## 6     ShannonsH     numdf                 2                 2
    ## 7     ShannonsH     dendf                50                49
    ## 8     ShannonsH r.squared                 0                 0
    ## 9    Sapro.rich  estim.se  0.001 +/- 0.003  -0.004 +/- 0.003 
    ## 10   Sapro.rich     numdf                 2                 2
    ## 11   Sapro.rich     dendf                50                49
    ## 12   Sapro.rich r.squared                 0              0.04
    ## 13 Basidio.rich  estim.se      0 +/- 0.001  -0.001 +/- 0.001 
    ## 14 Basidio.rich     numdf                 2                 2
    ## 15 Basidio.rich     dendf                50                49
    ## 16 Basidio.rich r.squared                 0              0.03
    ## 17   Patho.rich  estim.se      0 +/- 0.003  -0.001 +/- 0.002 
    ## 18   Patho.rich     numdf                 2                 2
    ## 19   Patho.rich     dendf                50                49
    ## 20   Patho.rich r.squared                 0                 0
    ## 21    Oomy.rich  estim.se  -0.027 +/- 0.03    0.038 +/- 0.03 
    ## 22    Oomy.rich     numdf                 2                 2
    ## 23    Oomy.rich     dendf                50                49
    ## 24    Oomy.rich r.squared              0.02              0.03
    ##               time37             time7
    ## 1           0 +/- 0           0 +/- 0 
    ## 2                  2                 2
    ## 3                 48                51
    ## 4               0.01              0.03
    ## 5   0.016 +/- 0.027   -0.016 +/- 0.02 
    ## 6                  2                 2
    ## 7                 48                51
    ## 8               0.01              0.01
    ## 9       0 +/- 0.003  -0.002 +/- 0.002 
    ## 10                 2                 2
    ## 11                48                51
    ## 12                 0              0.03
    ## 13 -0.001 +/- 0.001  -0.001 +/- 0.001 
    ## 14                 2                 2
    ## 15                48                51
    ## 16              0.03              0.01
    ## 17 -0.003 +/- 0.002  -0.001 +/- 0.002 
    ## 18                 2                 2
    ## 19                48                51
    ## 20              0.04                 0
    ## 21  0.024 +/- 0.025  -0.009 +/- 0.019 
    ## 22                 2                 2
    ## 23                48                51
    ## 24              0.02                 0

############################################## 

Relationship between wood traits and community
==============================================

*Hyp (species+size-level)* Initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

*Hyp (stem-level)* Average initial microbial communitiy compositions will covary with initial wood traits
---------------------------------------------------------------------------------------------------------

############################################## 

Extra pieces
------------

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
