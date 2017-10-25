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

#fxns
source("code/load_fxns.R")
source("code/curve_fitting_fxns.R")
source("code/distance_fxns.R")
```

### Load microbial community data

``` r
stemSamples<-load_stemSamples() #load stem sample meta data
fung.otu<-load_matotu() #load the fungal OTU table
comm.otu<-add_oomycetes(fung.otu) #add the oomycetes
seqSamples<-load_seqSamples(comm.otu, stemSamples)
#plot_sampleEffortCurves(comm.otu)
```

### Load wood trait data

``` r
traits.mean<-mergeTraitData()
traits.long<-as.data.frame(gather(traits.mean, key=trait, value=value, -(1:3)))

#missing data
filter(traits.long, is.na(value))
```

    ##   code species_lower  size     trait value
    ## 1 olst          olst small waterperc   NaN
    ## 2 eusc          eusc small         P    NA
    ## 3 eusc          eusc small         K    NA
    ## 4 eusc          eusc small        Ca    NA
    ## 5 eusc          eusc small        Mn    NA
    ## 6 eusc          eusc small        Fe    NA
    ## 7 eusc          eusc small        Zn    NA
    ## 8 eusc          eusc small         N    NA
    ## 9 eusc          eusc small         C    NA

``` r
# ggplot(traits.long, aes(x=species_lower, y=value, color=size)) + 
#   geom_point() + 
#   facet_wrap(~trait, scales="free") +
#   mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# #the Fe values for large stems are not 0, they are just small values
# subset(traits.long, trait=="Fe") %>%
#   ggplot(aes(x=species_lower, y=value, color=size)) + 
#   geom_point() +
#   facet_wrap(~size, scales="free")
```

### Load mass loss data

``` r
initial_mass <- read_in_initial_mass()
harvest_mass<-LoadHarvestFiles()
mass.data<-bind_rows(initial_mass, harvest_mass)

#look for outliers

# mass.data %>% ggplot(aes(x=time, y=totalSampleDryMass)) + geom_point(alpha=0.6)+theme_bw() 
#looks good

# mass.data %>% ggplot(aes(x=time, y=totalSampleDryMass,col=size)) + geom_point(position="jitter",alpha=0.6)+theme_bw()+scale_y_log10()
# two high values in size==small and harvest 3 are likely real, they have been checked

#mass.data[which(mass.data$totalSampleDryMass==0),]
#no longer any samples with 0 totalSampleDryMass in the dataset

#check for missing data
mass.data %>%
  filter(is.na(totalSampleDryMass)) %>%
  knitr::kable()
```

| unique  | species | size  |  time|  totalSampleDryMass|  density| fruiting | insects | drill | notes                                  |
|:--------|:--------|:------|-----:|-------------------:|--------:|:---------|:--------|:------|:---------------------------------------|
| olst1e  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1j  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1i  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2j  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2b  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2d  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4h  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1c  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst3a  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst3b  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4e  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1a  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4f  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2i  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2e  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4i  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2g  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1d  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2h  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2c  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4c  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2f  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst2a  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4b  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1b  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1k  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1g  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1h  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4g  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4a  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst3c  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst1f  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| olst4d  | olst    | small |     0|                  NA|       NA| NA       | NA      | NA    | NA                                     |
| ALLI311 | ALLI    | large |    37|                  NA|       NA|          | NA      | no    | missing from plot -- missing           |
| baae1a  | baae    | small |    37|                  NA|       NA|          | NA      | no    | missing from plot -- missing from plot |

Merge time zero with the other harvests to calculate proportion mass remaining at each time point... Matching failures, all due to missing time 0 data

``` r
#Merge time zero with the other harvests to calculate proportion mass remaining at each time point
mass.data %>%
  filter(time==0) %>%
  rename(timeZeroDensity=density) %>%
  rename(timeZeroMass=totalSampleDryMass) %>%
  select(unique,timeZeroMass,timeZeroDensity)->time_zero

mass.data %>%
  left_join(time_zero,by="unique") %>%
  mutate(pmr=totalSampleDryMass/timeZeroMass) %>%
  mutate(SpeciesCode=tolower(species)) -> plotting_df
  write_csv(plotting_df,"derived_data/plotting_df.csv")
  
  
# here are the matching failures which are currently due to the time zero adjustment for moisture
plotting_df %>%
  filter(is.na(pmr)) %>%
  select(unique, species, size, time, totalSampleDryMass, notes) %>%
  spread(key=time, value=totalSampleDryMass) %>%
  knitr::kable()
```

| unique  | species | size  | notes                                  |    0|      7|    13|    25|    37|
|:--------|:--------|:------|:---------------------------------------|----:|------:|-----:|-----:|-----:|
| ALLI311 | ALLI    | large | missing from plot -- missing           |   NA|     NA|    NA|    NA|    NA|
| baae1a  | baae    | small | missing from plot -- missing from plot |   NA|     NA|    NA|    NA|    NA|
| olst1a  | olst    | small |                                        |   NA|   3.15|    NA|    NA|    NA|
| olst1a  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1b  | olst    | small |                                        |   NA|     NA|    NA|  4.28|    NA|
| olst1b  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1c  | olst    | small | all wwe -- all wet weight excess       |   NA|     NA|    NA|    NA|  0.99|
| olst1c  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1d  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1e  | olst    | small | all crumbled once handled --           |   NA|     NA|    NA|    NA|  0.77|
| olst1e  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1f  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1g  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1h  | olst    | small |                                        |   NA|     NA|  8.64|    NA|    NA|
| olst1h  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1i  | olst    | small |                                        |   NA|  10.20|    NA|    NA|    NA|
| olst1i  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1j  | olst    | small |                                        |   NA|  10.41|    NA|    NA|    NA|
| olst1j  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst1k  | olst    | small |                                        |   NA|     NA|  6.40|    NA|    NA|
| olst1k  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2a  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2b  | olst    | small | -- bark all off                        |   NA|     NA|    NA|    NA|  4.85|
| olst2b  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2c  | olst    | small | broken, rotted                         |   NA|     NA|    NA|  4.65|    NA|
| olst2c  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2d  | olst    | small | 2 pieces                               |   NA|     NA|    NA|  3.94|    NA|
| olst2d  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2e  | olst    | small | bark flaking off                       |   NA|     NA|    NA|  5.31|    NA|
| olst2e  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2f  | olst    | small |                                        |   NA|  11.06|    NA|    NA|    NA|
| olst2f  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2g  | olst    | small | thin bark almost all gone              |   NA|     NA|    NA|  2.75|    NA|
| olst2g  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2h  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2i  | olst    | small |                                        |   NA|   5.91|    NA|    NA|    NA|
| olst2i  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst2j  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst3a  | olst    | small |                                        |   NA|     NA|    NA|  2.93|    NA|
| olst3a  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst3b  | olst    | small | --                                     |   NA|     NA|    NA|    NA|  1.53|
| olst3b  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst3c  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4a  | olst    | small | -- looks intact/slightly soft          |   NA|     NA|    NA|    NA|  2.58|
| olst4a  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4b  | olst    | small |                                        |   NA|     NA|  7.37|    NA|    NA|
| olst4b  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4c  | olst    | small |                                        |   NA|     NA|  4.60|    NA|    NA|
| olst4c  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4d  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4e  | olst    | small |                                        |   NA|     NA|  2.99|    NA|    NA|
| olst4e  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4f  | olst    | small |                                        |   NA|   5.55|    NA|    NA|    NA|
| olst4f  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4g  | olst    | small | -- intact, no rot                      |   NA|     NA|    NA|    NA|  3.32|
| olst4g  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4h  | olst    | small |                                        |   NA|     NA|  2.13|    NA|    NA|
| olst4h  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|
| olst4i  | olst    | small | NA                                     |   NA|     NA|    NA|    NA|    NA|

``` r
#remove NA rows
plotting_df %>%
  filter(!is.na(pmr)) -> plotting_df
```

### Non-linear curve fits of decay trajectories

Using `litterfitter` to apply both negative exponenial and weibull to all species/size classes

``` r
#spdf <- fit_all_curves(plotting_df) #this recalculates all the curve fits, uncomment if the data changes
#write_csv(spdf,"derived_data/mass_loss_parameters.csv")
spdf <- read_csv("derived_data/mass_loss_parameters.csv")
indx<-select(stemSamples, code, species, size)
spdf<-left_join(spdf, indx) #add code

# ggplot(spdf,aes(x=t70,y=w.t70,col=size))+
#   geom_point()+
#   labs(x="Time to 30% mass loss (negative exponential)", 
#        y="Time to 30% mass loss (Weibull)")+
#   geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()
```

### Plot microbial community distances vs decay param distances between species+size

``` r
# community composition, bray
pList<-Make_commDistvDist_Fig(distType="bray", valueCol_vec=c("k","alpha"), seqSamples, stemSamples, comm.otu, spdf)
grid.arrange(pList[['k']], pList[['alpha']])
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png)

``` r
# community composition, jaccard
pList<-Make_commDistvDist_Fig(distType="jaccard", valueCol_vec=c("k","alpha"), seqSamples, stemSamples, comm.otu, spdf)
grid.arrange(pList[['k']], pList[['alpha']])
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-2.png)

``` r
# richness
pList<-Make_commDistvDist_Fig(distType="richness", valueCol_vec=c("k","alpha"), seqSamples, stemSamples, comm.otu, spdf)
grid.arrange(pList[['k']], pList[['alpha']])
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-3.png) Not a lot of information in initial mean microbial composition distance (bray or jaccard) about differences in species+size decay trajectories. Likewise, richness is not useful.

### Plot wood trait distances vs decay param distances between species+size

``` r
pList<-Make_woodTraitDistvDist_Fig(valueCol_vec=c("k","alpha"), seqSamples, stemSamples, traits.mean, spdf)
grid.arrange(pList[['k']], pList[['alpha']])
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png) Seems like there is more information in inital wood trait distance about species+size decay trajectories than there is in the initial community composition/richness

### Plot species+size beta diversity of microbial community vs AIC of neg.exponential

``` r
#calculate pairwise community distances
comm.dist<-Calc_commDists(seqSamples, comm.otu, distType="bray") #2. calc the distances
filter(comm.dist, code1==code2) %>% #3. isolate just the distances within species+size
  select(sampID1, sampID2, code1, size, dist) %>%
  rename("code"="code1") %>% #rename the cols
  group_by(code) %>%
  summarize(mean=mean(dist), #calculate the mean community dist within code classes
            se=sd(dist)/sqrt(length(dist)),
            upper=mean+se,
            lower=mean-se) -> comm.dist.wth

#combine with decay trajectory params
spdf.sub<-select(spdf, code, neg.exp.aic, w.aic)
comm.dist.wth<-left_join(comm.dist.wth, spdf.sub)

#add back the size and species columns
comm.dist.wth$size<-"large"
comm.dist.wth[tolower(comm.dist.wth$code) == comm.dist.wth$code, "size"]<-"small"
comm.dist.wth$species<-tolower(comm.dist.wth$code)

p.negexp.aic<-ggplot(comm.dist.wth, aes(x=mean, y=neg.exp.aic, color=species, shape=size)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlab("Mean microbial community distance") + 
  ylab("Negative exponential model AIC") +
  facet_grid(~size)

p.w.aic<-ggplot(comm.dist.wth, aes(x=mean, y=w.aic, color=species, shape=size)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlab("Mean microbial community distance") + 
  ylab("Weibull model AIC") +
  facet_grid(~size)

#pdf('output/commDist_decayparam_withinCode.pdf', width=8, height=6)
grid.arrange(p.negexp.aic, p.w.aic)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png)

``` r
#dev.off()
```

Expect to see a positive relationship between R2 and within species+size class microbial community distance. In other words, expect that samples with similar initial microbial communities will have better-fitting decay models.
