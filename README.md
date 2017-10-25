Does chemistry or community better predict mass loss?
================
Marissa Lee
10/23/2017

``` r
#chunk options
knitr::opts_chunk$set(echo = TRUE)

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
fung.otu<-load_matotu()
comm.otu<-add_oomycetes(fung.otu)
```

    ## Warning: Column `seqSamp` joining factors with different levels, coercing
    ## to character vector

``` r
#plot_sampleEffortCurves(comm.otu)
```

### Load wood trait data

``` r
traits.mean<-mergeTraitData()
```

    ## Warning: Column `SampleCode` joining character vector and factor, coercing
    ## into character vector

``` r
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
ggplot(traits.long, aes(x=species_lower, y=value, color=size)) + 
  geom_point() + 
  facet_wrap(~trait, scales="free") +
  mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

    ## Warning: Removed 9 rows containing missing values (geom_point).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png)

``` r
# #the Fe values for large stems are not 0, they are just small values
# subset(traits.long, trait=="Fe") %>%
#   ggplot(aes(x=species_lower, y=value, color=size)) + 
#   geom_point() +
#   facet_wrap(~size, scales="free")
```

### Load mass loss data

``` r
initial_mass <- read_in_initial_mass()
```

    ## Parsed with column specification:
    ## cols(
    ##   `random order num` = col_integer(),
    ##   Species = col_character(),
    ##   Stem = col_integer(),
    ##   ID = col_integer(),
    ##   unique = col_character(),
    ##   `Diameter (cm)` = col_double(),
    ##   `Fresh mass (g)` = col_double(),
    ##   `Dry mass (g)` = col_double(),
    ##   notes = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   `random order num` = col_character(),
    ##   Species = col_character(),
    ##   Stem = col_integer(),
    ##   ID = col_character(),
    ##   unique = col_character(),
    ##   `Diameter.wbark (mm)` = col_double(),
    ##   `Fresh mass (g)` = col_double(),
    ##   `Diameter.nobark (mm)` = col_double(),
    ##   `Volume (g)` = col_double(),
    ##   `Dry mass wood (g)` = col_double(),
    ##   `Dry mass bark (g)` = col_double(),
    ##   `Dry mass total (g)` = col_double(),
    ##   Notes = col_character()
    ## )

    ## Joining, by = "Species"
    ## Joining, by = "Species"

``` r
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

| unique  | Species | size  |  time|  totalSampleDryMass|  density| fruiting | insects | drill | notes                                  |
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
  mutate(SpeciesCode=tolower(Species)) -> plotting_df
  write_csv(plotting_df,"derived_data/plotting_df.csv")
  
  
# here are the matching failures which are currently due to the time zero adjustment for moisture
plotting_df %>%
  filter(is.na(pmr)) %>%
  select(unique, Species, size, time, totalSampleDryMass, notes) %>%
  spread(key=time, value=totalSampleDryMass) %>%
  knitr::kable()
```

| unique  | Species | size  | notes                                  |    0|      7|    13|    25|    37|
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

### Non-linear curve fits of decay trajectories

Using `litterfitter` to apply both negative exponenial and weibull to all species/size classes

``` r
#spdf <- fit_all_curves(plotting_df) #this recalculates all the curve fits, uncomment if the data changes
#write_csv(spdf,"derived_data/mass_loss_parameters.csv")
spdf <- read_csv("derived_data/mass_loss_parameters.csv")

ggplot(spdf,aes(x=t70,y=w.t70,col=size))+
  geom_point()+
  labs(x="Time to 30% mass loss (negative exponential)", 
       y="Time to 30% mass loss (Weibull)")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png)

### Testing effects of t=0 points

replotting the last example without the t=0 points

``` r
plotting_df %>%
  filter(SpeciesCode=="ripi",size=="small",time>0, !is.na(pmr)) ->out

plot_multiple_fits(time = out$time/12,
                   mass.remaining = out$pmr,
                   bty = 'n', model = c('neg.exp', 'weibull'),
                   xlab = 'Time', ylab = 'Proportion mass remaining',iters=1000)
```

    ## Number of successful fits:  995  out of 1000 
    ## Number of successful fits:  1000  out of 1000

    ## Warning in multioptimFit(time, mass.remaining, model, iters = iters, upper
    ## = upper, : May not have found global best fit; increase iterations

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png) Checking that the fits are the same for weibull which they are

``` r
out%$%
fit_litter(time = time/12, 
        mass.remaining = pmr, model = c("weibull"), iters = 1000) ->plot_1
```

    ## Number of successful fits:  1000  out of 1000

``` r
print(plot_1)
```

    ## $optimFit
    ## $optimFit$par
    ## [1] 2.148217 1.579560
    ## 
    ## $optimFit$value
    ## [1] 0.5958932
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       19       19 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 10.29445
    ## 
    ## $fitAIC
    ## [1] -16.58889
    ## 
    ## $fitAICc
    ## [1] -16.01746
    ## 
    ## $fitBIC
    ## [1] -14.23278
    ## 
    ## $time
    ##  [1] 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 1.0833333
    ##  [8] 1.0833333 1.0833333 1.0833333 1.0833333 1.0833333 2.0833333 2.0833333
    ## [15] 2.0833333 2.0833333 2.0833333 2.0833333 3.0833333 3.0833333 3.0833333
    ## [22] 3.0833333 3.0833333 3.0833333
    ## 
    ## $mass
    ##  [1] 0.80702715 0.70635658 0.73831715 0.98531286 0.98532447 0.97733591
    ##  [7] 0.44745619 0.48033810 0.60704277 0.88422473 0.92666205 0.89016632
    ## [13] 0.15996872 0.25228526 0.40869266 0.47153137 0.60294025 0.59227957
    ## [19] 0.00000000 0.14378464 0.03073876 0.39355135 0.16262361 0.14137293
    ## 
    ## $predicted
    ##  [1] 0.8802407 0.8802407 0.8802407 0.8802407 0.8802407 0.8802407 0.7123870
    ##  [8] 0.7123870 0.7123870 0.7123870 0.7123870 0.7123870 0.3856940 0.3856940
    ## [15] 0.3856940 0.3856940 0.3856940 0.3856940 0.1703840 0.1703840 0.1703840
    ## [22] 0.1703840 0.1703840 0.1703840
    ## 
    ## $model
    ## [1] "weibull"
    ## 
    ## $nparams
    ## [1] 2
    ## 
    ## attr(,"class")
    ## [1] "litfit"

``` r
plot(plot_1)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png)

``` r
plotting_df %>%
  filter(SpeciesCode=="ripi",size=="small", !is.na(pmr)) ->out
out%$%
fit_litter(time = time/12, 
        mass.remaining = pmr, model = c("weibull"), iters = 1000) ->plot_2
```

    ## Number of successful fits:  1000  out of 1000

``` r
print(plot_2)
```

    ## $optimFit
    ## $optimFit$par
    ## [1] 2.148217 1.579560
    ## 
    ## $optimFit$value
    ## [1] 0.5958932
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       29       29 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 53.22484
    ## 
    ## $fitAIC
    ## [1] -102.4497
    ## 
    ## $fitAICc
    ## [1] -102.2391
    ## 
    ## $fitBIC
    ## [1] -98.26099
    ## 
    ## $time
    ##  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    ##  [8] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    ## [15] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    ## [22] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    ## [29] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    ## [36] 0.0000000 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333
    ## [43] 1.0833333 1.0833333 1.0833333 1.0833333 1.0833333 1.0833333 2.0833333
    ## [50] 2.0833333 2.0833333 2.0833333 2.0833333 2.0833333 3.0833333 3.0833333
    ## [57] 3.0833333 3.0833333 3.0833333 3.0833333
    ## 
    ## $mass
    ##  [1] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ##  [7] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ## [13] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ## [19] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ## [25] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ## [31] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ## [37] 0.80702715 0.70635658 0.73831715 0.98531286 0.98532447 0.97733591
    ## [43] 0.44745619 0.48033810 0.60704277 0.88422473 0.92666205 0.89016632
    ## [49] 0.15996872 0.25228526 0.40869266 0.47153137 0.60294025 0.59227957
    ## [55] 0.00000000 0.14378464 0.03073876 0.39355135 0.16262361 0.14137293
    ## 
    ## $predicted
    ##  [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ##  [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [22] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [29] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [36] 1.0000000 0.8802407 0.8802407 0.8802407 0.8802407 0.8802407 0.8802407
    ## [43] 0.7123870 0.7123870 0.7123870 0.7123870 0.7123870 0.7123870 0.3856940
    ## [50] 0.3856940 0.3856940 0.3856940 0.3856940 0.3856940 0.1703840 0.1703840
    ## [57] 0.1703840 0.1703840 0.1703840 0.1703840
    ## 
    ## $model
    ## [1] "weibull"
    ## 
    ## $nparams
    ## [1] 2
    ## 
    ## attr(,"class")
    ## [1] "litfit"

``` r
plot(plot_2)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-2.png) Conclusion: including t=0 points affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

### Plot beta diversity of microbial community vs distance in decay params

``` r
#calculate pairwise community distances
sampTab<-CreateSeqSampTab(mass.data) #1. creat sample tab for community sequences samples
```

    ## Joining, by = "code"

    ## Warning: Column `code` joining factor and character vector, coercing into
    ## character vector

``` r
comm.dist<-Calc_commDists(sampTab, comm.otu) #2. calc the distances
```

    ## Warning: Column `sampID1`/`sampID` joining factors with different levels,
    ## coercing to character vector

    ## Warning: Column `sampID2`/`sampID` joining factors with different levels,
    ## coercing to character vector

``` r
summ.comm_dist<-SummarizeCommDist_byCodePair(comm.dist) #3. summarize the distances by codePair

#calculate pairwise decay parameter distances
spdf<-mutate(spdf, minAIC = pmin(neg.exp.aic, w.aic) ) #1. find the min AIC for each species+size
spdf<-AddCodeID(sampTab) #2. add code to spdf using sampTab
```

    ## Joining, by = c("species", "size")

``` r
aic.dist<-Calc_decayParamDiffs(valueCol="minAIC", spdf, sampTab) #3. calc decay param diffs, e.g. minAIC
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
k.dist<-Calc_decayParamDiffs(valueCol="k", spdf, sampTab)
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
alpha.dist<-Calc_decayParamDiffs(valueCol="alpha", spdf, sampTab)
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
#merge community and decay param distances into the same dataframe
comm_decay.dist.aic<-MergeCommNDecaydists_byCodePair(decayparam.dist=aic.dist, summ.comm_dist)
```

    ## Joining, by = "codePair"

``` r
comm_decay.dist.k<-MergeCommNDecaydists_byCodePair(decayparam.dist=k.dist, summ.comm_dist)
```

    ## Joining, by = "codePair"

``` r
comm_decay.dist.alpha<-MergeCommNDecaydists_byCodePair(decayparam.dist=alpha.dist, summ.comm_dist)
```

    ## Joining, by = "codePair"

``` r
# delta minAIC
p.aic<-ggplot(comm_decay.dist.aic, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model AIC")

# delta k
p.k<-ggplot(comm_decay.dist.k, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model k")

# delta alpha
p.alpha<-ggplot(comm_decay.dist.alpha, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model alpha")

#pdf('output/commDist_decayparamDist_byCode.pdf', width=4, height=9)
grid.arrange(p.aic, p.k, p.alpha)
```

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png)

``` r
#dev.off()
```

Not a lot of information in mean inital community distance about species+size decay trajectories

### Plot beta diversity of wood traits vs distance in decay params

``` r
#calculate pairwise wood trait distances
traits.dist<-Calc_woodTraitDist(traits.mean)

#calculate pairwise decay parameter distances
spdf<-mutate(spdf, minAIC = pmin(neg.exp.aic, w.aic) ) #1. find the min AIC for each species+size
spdf<-AddCodeID(sampTab) #2. add code to spdf using sampTab
```

    ## Joining, by = c("species", "size", "code")

``` r
aic.dist<-Calc_decayParamDiffs(valueCol="minAIC", spdf, sampTab) #3. calc decay param diffs, e.g. minAIC
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
k.dist<-Calc_decayParamDiffs(valueCol="k", spdf, sampTab)
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
alpha.dist<-Calc_decayParamDiffs(valueCol="alpha", spdf, sampTab)
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
#merge wood trait and decay param distances into the same dataframe
wood_decay.dist.aic<-MergeWoodNDecaydists_byCodePair(decayparam.dist=aic.dist, traits.dist)
```

    ## Joining, by = "codePair"

``` r
wood_decay.dist.k<-MergeWoodNDecaydists_byCodePair(decayparam.dist=k.dist, traits.dist)
```

    ## Joining, by = "codePair"

``` r
wood_decay.dist.alpha<-MergeWoodNDecaydists_byCodePair(decayparam.dist=alpha.dist, traits.dist)
```

    ## Joining, by = "codePair"

``` r
# delta minAIC
p.aic<-ggplot(wood_decay.dist.aic, aes(x=woodTraitDist, y=decayparam_dist, color=size)) + 
  geom_point() + 
  xlab('Wood trait distance') + ylab("Delta decay model AIC")

# delta k
p.k<-ggplot(wood_decay.dist.k, aes(x=woodTraitDist, y=decayparam_dist, color=size)) + 
  geom_point() + 
  xlab('Wood trait distance') + ylab("Delta decay model k")

# delta alpha
p.alpha<-ggplot(wood_decay.dist.alpha, aes(x=woodTraitDist, y=decayparam_dist, color=size)) + 
  geom_point() + 
  xlab('Wood trait distance') + ylab("Delta decay model alpha")

#pdf('output/woodTraitDist_decayparamDist_byCode.pdf', width=4, height=9)
grid.arrange(p.aic, p.k, p.alpha)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png)

``` r
#dev.off()
```

Seems like there is more information in inital wood trait distance about species+size decay trajectories than there is in the initial mean communities distances

### Plot beta *richness* of microbial community vs distance in decay params

``` r
#calculate pairwise community distances
sampTab<-CreateSeqSampTab(mass.data) #1. creat sample tab for community sequences samples
```

    ## Joining, by = "code"

    ## Warning: Column `code` joining factor and character vector, coercing into
    ## character vector

``` r
rich.dist<-Calc_richDists(sampTab, comm.otu) #2. calc the distances
```

    ## Warning: Column `sampID1`/`sampID` joining factors with different levels,
    ## coercing to character vector

    ## Warning: Column `sampID2`/`sampID` joining factors with different levels,
    ## coercing to character vector

``` r
summ.rich_dist<-SummarizeCommDist_byCodePair(rich.dist) #3. summarize the distances by codePair
 
#calculate pairwise decay parameter distances
spdf<-mutate(spdf, minAIC = pmin(neg.exp.aic, w.aic) ) #1. find the min AIC for each species+size
spdf<-AddCodeID(sampTab) #2. add code to spdf using sampTab
```

    ## Joining, by = c("species", "size", "code")

``` r
aic.dist<-Calc_decayParamDiffs(valueCol="minAIC", spdf, sampTab) #3. calc decay param diffs, e.g. minAIC
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
k.dist<-Calc_decayParamDiffs(valueCol="k", spdf, sampTab)
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
alpha.dist<-Calc_decayParamDiffs(valueCol="alpha", spdf, sampTab)
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
#merge community and decay param distances into the same dataframe
rich_decay.dist.aic<-MergeCommNDecaydists_byCodePair(decayparam.dist=aic.dist, summ.rich_dist)
```

    ## Joining, by = "codePair"

``` r
rich_decay.dist.k<-MergeCommNDecaydists_byCodePair(decayparam.dist=k.dist, summ.rich_dist)
```

    ## Joining, by = "codePair"

``` r
rich_decay.dist.alpha<-MergeCommNDecaydists_byCodePair(decayparam.dist=alpha.dist, summ.rich_dist)
```

    ## Joining, by = "codePair"

``` r
# delta minAIC
p.aic<-ggplot(rich_decay.dist.aic, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial richness distance') + ylab("Delta decay model AIC")

# delta k
p.k<-ggplot(rich_decay.dist.k, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial richness distance') + ylab("Delta decay model k")

# delta alpha
p.alpha<-ggplot(rich_decay.dist.alpha, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial richness distance') + ylab("Delta decay model alpha")

# #pdf('output/richDist_decayparamDist_byCode.pdf', width=4, height=9)
grid.arrange(p.aic, p.k, p.alpha)
```

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png)

``` r
# #dev.off()
```

Really not a lot of information in mean inital richness distance about species+size decay trajectories

### Look at 3-way relationship between species+size distances: wood trait, microbial community, and k

``` r
comm_wood_decay.dist.k<-left_join(comm_decay.dist.k, wood_decay.dist.k)
```

    ## Joining, by = c("codePair", "size", "decayparam_dist")

``` r
#what is up with the missing data in woodTraitDist?
select(comm_wood_decay.dist.k, size, mean_comm_dist, woodTraitDist, decayparam_dist) %>%
  filter(complete.cases(size, comm_wood_decay.dist.k, mean_comm_dist, woodTraitDist, decayparam_dist)) -> df.k

p1<-ggplot(df.k, aes(x=mean_comm_dist, y=decayparam_dist, color=size)) + 
  geom_point() + xlab("Mean microbial community distance") + ylab("Decay rate (k) distance")

p2<-ggplot(df.k, aes(x=woodTraitDist, y=decayparam_dist, color=size)) + 
  geom_point() +  xlab("Wood trait distance") + ylab("Decay rate (k) distance")

p3<-ggplot(df.k, aes(x=mean_comm_dist, y=woodTraitDist, color=size)) + 
  geom_point() +  xlab("Mean microbial community distance") + ylab("Wood trait distance")

#pdf('output/triDist_k_byCode.pdf', width=6, height=6)
grid.arrange(p1,p2,p3, ncol=2)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-1.png)

``` r
#dev.off()

cor(df.k[,-1])
```

    ##                 mean_comm_dist woodTraitDist decayparam_dist
    ## mean_comm_dist      1.00000000    0.07362931      0.09170396
    ## woodTraitDist       0.07362931    1.00000000      0.26622256
    ## decayparam_dist     0.09170396    0.26622256      1.00000000

Another view that emphasizes the relationship between wood trait distance and decay rate distances between species+size classes. Mean microbial community distance is not correlated with either variable.

### Plot species+size beta diversity of microbial community vs AIC of neg.exponential

``` r
#calculate pairwise community distances
sampTab<-CreateSeqSampTab(mass.data) #1. creat sample tab for community sequences samples
```

    ## Joining, by = "code"

    ## Warning: Column `code` joining factor and character vector, coercing into
    ## character vector

``` r
comm.dist<-Calc_commDists(sampTab, comm.otu) #2. calc the distances
```

    ## Warning: Column `sampID1`/`sampID` joining factors with different levels,
    ## coercing to character vector

    ## Warning: Column `sampID2`/`sampID` joining factors with different levels,
    ## coercing to character vector

``` r
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
```

    ## Joining, by = "code"

``` r
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

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_errorbarh).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-13-1.png)

``` r
#dev.off()
```

Expected to see a positive relationship between AIC and within species+size class microbial community distance. In other words, expected that samples with similar initial microbial communities would have better-fitting decay models.

The large size samples show no pattern The small size samples show (maybe) the opposite pattern
