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
library(litterfitter)
library(magrittr)
library(tidyr)
source("code/load_fxns.R")
source("code/curve_fitting_fxns.R")
source("code/distance_fxns.R")
library(gridExtra)
```

### Load microbial community data

``` r
fung.otu<-load_matotu()
comm.otu<-add_oomycetes(fung.otu)
```

    ## Joining, by = "seqSamp"

    ## Warning: Column `seqSamp` joining factors with different levels, coercing
    ## to character vector

``` r
#plot_sampleEffortCurves(comm.otu)
```

### Load wood trait data

``` r
traits.mean<-mergeTraitData()
```

    ## Joining, by = "SampleCode"

    ## Warning: Column `SampleCode` joining character vector and factor, coercing
    ## into character vector

    ## Joining, by = "code"

    ## Joining, by = "code"
    ## Joining, by = "code"

``` r
traits.long<-as.data.frame(gather(traits.mean, key=trait, value=value, -(1:3)))

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
#initial mass
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
# initial_mass %>%
#   mutate(SpeciesCode=tolower(Species))%>%
#   ggplot(aes(y=totalSampleDryMass,x=SpeciesCode,fill=size))+
#   geom_violin()+ 
#   theme(axis.text.x=element_text(angle=90,hjust=1)) + 
#   scale_y_log10()

#mass at harvest
harvest_mass<-LoadHarvestFiles()
```

    ## Warning in read.samp3(): NAs introduced by coercion

``` r
#create a complete sample mass df for all time points
mass.data<-bind_rows(initial_mass, harvest_mass)
```

Identify outliers in mass loss data

``` r
mass.data %>% 
  ggplot(aes(x=time, y=totalSampleDryMass)) + geom_point(alpha=0.6)+theme_bw() 
```

    ## Warning: Removed 33 rows containing missing values (geom_point).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png)

``` r
mass.data[which.max(mass.data$totalSampleDryMass),]
```

    ## # A tibble: 1 x 10
    ##    unique Species  size  time totalSampleDryMass density
    ##     <chr>   <chr> <chr> <dbl>              <dbl>   <dbl>
    ## 1 ALLI111    ALLI large    37            1242.64    0.74
    ## # ... with 4 more variables: fruiting <chr>, insects <chr>, drill <chr>,
    ## #   notes <chr>

``` r
outlier.uniques<-as.character(mass.data[which.max(mass.data$totalSampleDryMass),"unique"])
outlier.uniques
```

    ## [1] "ALLI111"

...another view...might want to check out those two high mass value outliers from harvest 3

``` r
mass.data %>%
  ggplot(aes(x=time, y=totalSampleDryMass,col=size)) +
  geom_point(position="jitter",alpha=0.6)+theme_bw()+scale_y_log10()
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 33 rows containing missing values (geom_point).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-6-1.png)

``` r
#Might want to check out those two high mass value outliers from harvest 3
filter(mass.data, size=="small", time==25) %>%
  arrange(-totalSampleDryMass)->tmp
tmp[1:2,]
```

    ## # A tibble: 2 x 10
    ##   unique Species  size  time totalSampleDryMass density fruiting insects
    ##    <chr>   <chr> <chr> <dbl>              <dbl>   <dbl>    <chr>   <chr>
    ## 1 alli3l    alli small    25              18.37    0.65   hyphae        
    ## 2 hase1k    hase small    25              16.04    0.60                 
    ## # ... with 2 more variables: drill <chr>, notes <chr>

``` r
outlier.uniques<-c(outlier.uniques, tmp$unique[1:2])
```

Are these real 0's?

``` r
mass.data[which(mass.data$totalSampleDryMass==0),]
```

    ## # A tibble: 5 x 10
    ##    unique Species  size  time totalSampleDryMass density fruiting insects
    ##     <chr>   <chr> <chr> <dbl>              <dbl>   <dbl>    <chr>   <chr>
    ## 1  ripi1j    ripi small    37                  0     NaN                4
    ## 2 ALLI311    ALLI large    37                  0     NaN             <NA>
    ## 3  hase2b    hase small    37                  0     NaN                4
    ## 4  baae1a    baae small    37                  0     NaN             <NA>
    ## 5  eute1e    eute small    37                  0     NaN                4
    ## # ... with 2 more variables: drill <chr>, notes <chr>

Remove outliers from mass.data and merge time zero with the other harvests to calculate proportion mass remaining at each time point

``` r
mass.data<-mass.data[!mass.data$unique %in% outlier.uniques,]

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
```

### Non-linear curve fits of decay trajectories

Using `litterfitter` to apply both negative exponenial and weibull to all species/size classes

``` r
#spdf <- fit_all_curves(plotting_df) #this recalculates all the curve fits, uncomment if the data changes
spdf <- read_csv("derived_data/mass_loss_parameters.csv")
#write_csv(spdf,"derived_data/mass_loss_parameters.csv")

ggplot(spdf,aes(x=t70,y=w.t70,col=size))+
  geom_point()+
  labs(x="Time to 30% mass loss (negative exponential)", 
       y="weibull time to 30% mass loss")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png)

``` r
# #Here is an example where the neg.exp and the weibull curves are almost identical
# plotting_df %>%
#   filter(SpeciesCode=="eute",size=="small")  ->one_example
# plot_multiple_fits(time = one_example$time/12,
#                    mass.remaining = one_example$pmr,
#                    bty = 'n', model = c('neg.exp', 'weibull'),
#                    xlab = 'Time', ylab = 'Proportion mass remaining',iters=1000)
# #and one where they are pretty different:
# plotting_df %>%
#   filter(SpeciesCode=="ripi",size=="small") -> another_example
# plot_multiple_fits(time = another_example$time/12,
#                    mass.remaining = another_example$pmr,
#                    bty = 'n', model = c('neg.exp', 'weibull'),
#                    xlab = 'Time', ylab = 'Proportion mass remaining',iters=1000)
```

### Testing effects of t=0 points

replotting the last example without the t=0 points

``` r
#why filtering for ripi small?

plotting_df %>%
  filter(SpeciesCode=="ripi",size=="small",time>0) ->out

plot_multiple_fits(time = out$time/12,
                   mass.remaining = out$pmr,
                   bty = 'n', model = c('neg.exp', 'weibull'),
                   xlab = 'Time', ylab = 'Proportion mass remaining',iters=1000)
```

    ## Number of successful fits:  988  out of 1000 
    ## Number of successful fits:  1000  out of 1000

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png) Checking that the fits are the same for weibull which they are

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
    ## [1] 2.167791 1.576424
    ## 
    ## $optimFit$value
    ## [1] 0.5690418
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       22       22 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 10.84774
    ## 
    ## $fitAIC
    ## [1] -17.69548
    ## 
    ## $fitAICc
    ## [1] -17.12405
    ## 
    ## $fitBIC
    ## [1] -15.33937
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
    ## [13] 0.22928849 0.25228526 0.40869266 0.47153137 0.60294025 0.59227957
    ## [19] 0.00000000 0.14378464 0.03073876 0.39355135 0.16262361 0.14137293
    ## 
    ## $predicted
    ##  [1] 0.8813822 0.8813822 0.8813822 0.8813822 0.8813822 0.8813822 0.7153111
    ##  [8] 0.7153111 0.7153111 0.7153111 0.7153111 0.7153111 0.3909110 0.3909110
    ## [15] 0.3909110 0.3909110 0.3909110 0.3909110 0.1750647 0.1750647 0.1750647
    ## [22] 0.1750647 0.1750647 0.1750647
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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png)

``` r
#why filtering for ripi small?

plotting_df %>%
  filter(SpeciesCode=="ripi",size=="small") ->out
out%$%
fit_litter(time = time/12, 
        mass.remaining = pmr, model = c("weibull"), iters = 1000) ->plot_2
```

    ## Number of successful fits:  1000  out of 1000

    ## Warning in multioptimFit(time, mass.remaining, model, iters = iters, upper
    ## = upper, : May not have found global best fit; increase iterations

``` r
print(plot_2)
```

    ## $optimFit
    ## $optimFit$par
    ## [1] 2.167791 1.576424
    ## 
    ## $optimFit$value
    ## [1] 0.5690418
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       32       32 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 54.60807
    ## 
    ## $fitAIC
    ## [1] -105.2161
    ## 
    ## $fitAICc
    ## [1] -105.0056
    ## 
    ## $fitBIC
    ## [1] -101.0274
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
    ## [49] 0.22928849 0.25228526 0.40869266 0.47153137 0.60294025 0.59227957
    ## [55] 0.00000000 0.14378464 0.03073876 0.39355135 0.16262361 0.14137293
    ## 
    ## $predicted
    ##  [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ##  [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [22] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [29] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [36] 1.0000000 0.8813822 0.8813822 0.8813822 0.8813822 0.8813822 0.8813822
    ## [43] 0.7153111 0.7153111 0.7153111 0.7153111 0.7153111 0.7153111 0.3909110
    ## [50] 0.3909110 0.3909110 0.3909110 0.3909110 0.3909110 0.1750647 0.1750647
    ## [57] 0.1750647 0.1750647 0.1750647 0.1750647
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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-2.png) Conclusion: including t=0 points affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-1.png)

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-13-1.png)

``` r
#dev.off()
```

Seems like there is more information in inital wood trait distance about species+size decay trajectories than there is in the initial mean communities distances

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-14-1.png)

``` r
#dev.off()

cor(df.k[,-1])
```

    ##                 mean_comm_dist woodTraitDist decayparam_dist
    ## mean_comm_dist      1.00000000    0.07362931      0.09170396
    ## woodTraitDist       0.07362931    1.00000000      0.26622256
    ## decayparam_dist     0.09170396    0.26622256      1.00000000

Another view that emphasizes the relationship between wood trait distance and decay rate distances between species+size classes. Mean microbial community distance is not correlated with either variable.

### Plot with species+size beta diversity of microbial community vs AIC of neg.exponential

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

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-15-1.png)

``` r
#dev.off()
```

Expected to see a positive relationship between AIC and within species+size class microbial community distance. In other words, expected that samples with similar initial microbial communities would have better-fitting decay models.

The large size samples show no pattern The small size samples show (maybe) the opposite pattern
