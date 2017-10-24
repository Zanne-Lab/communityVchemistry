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
mass.data<-rbind(initial_mass, harvest_mass)
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

    ## # A tibble: 1 x 9
    ##    unique Species  size  time totalSampleDryMass   density
    ##     <chr>   <chr> <chr> <dbl>              <dbl>     <dbl>
    ## 1 ALLI111    ALLI large    37           1242.641 0.7351301
    ## # ... with 3 more variables: fruiting <chr>, insects <chr>, drill <chr>

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

    ## # A tibble: 2 x 9
    ##   unique Species  size  time totalSampleDryMass   density fruiting insects
    ##    <chr>   <chr> <chr> <dbl>              <dbl>     <dbl>    <chr>   <chr>
    ## 1 alli3l    alli small    25              18.37 0.6488873   hyphae        
    ## 2 hase1k    hase small    25              16.04 0.6005241                 
    ## # ... with 1 more variables: drill <chr>

``` r
outlier.uniques<-c(outlier.uniques, tmp$unique[1:2])
```

Are these real 0's?

``` r
mass.data[which(mass.data$totalSampleDryMass==0),]
```

    ## # A tibble: 10 x 9
    ##     unique Species  size  time totalSampleDryMass density fruiting insects
    ##      <chr>   <chr> <chr> <dbl>              <dbl>   <dbl>    <chr>   <chr>
    ##  1  ripi1j    ripi small    37                  0     NaN                4
    ##  2  eute2b    eute small    37                  0     NaN             <NA>
    ##  3  acel2f    acel small    37                  0     NaN                3
    ##  4  olst1c    olst small    37                  0     NaN             <NA>
    ##  5  eusc3j    eusc small    37                  0     NaN                3
    ##  6  olst1e    olst small    37                  0     NaN             <NA>
    ##  7 ALLI311    ALLI large    37                  0     NaN             <NA>
    ##  8  hase2b    hase small    37                  0     NaN                4
    ##  9  baae1a    baae small    37                  0     NaN             <NA>
    ## 10  eute1e    eute small    37                  0     NaN                4
    ## # ... with 1 more variables: drill <chr>

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


ggplot(spdf,aes(x=t70,y=w.t70,col=size))+geom_point()+labs(x="Time to 30% mass loss (negative exponential)", y="weibull time to 30% mass loss")+geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()
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

    ## Number of successful fits:  986  out of 1000 
    ## Number of successful fits:  1000  out of 1000

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png) Checking that the fits are the same for weibull which they are

``` r
out%$%
fit_litter(time = time/12, 
        mass.remaining = pmr, model = c("weibull"), iters = 1000) ->plot_1
```

    ## Number of successful fits:  999  out of 1000

``` r
print(plot_1)
```

    ## $optimFit
    ## $optimFit$par
    ## [1] 2.167896 1.575979
    ## 
    ## $optimFit$value
    ## [1] 0.5691124
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       38       38 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 10.84625
    ## 
    ## $fitAIC
    ## [1] -17.6925
    ## 
    ## $fitAICc
    ## [1] -17.12107
    ## 
    ## $fitBIC
    ## [1] -15.33639
    ## 
    ## $time
    ##  [1] 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 1.0833333
    ##  [8] 1.0833333 1.0833333 1.0833333 1.0833333 1.0833333 2.0833333 2.0833333
    ## [15] 2.0833333 2.0833333 2.0833333 2.0833333 3.0833333 3.0833333 3.0833333
    ## [22] 3.0833333 3.0833333 3.0833333
    ## 
    ## $mass
    ##  [1] 0.8074037 0.7059619 0.7383556 0.9853129 0.9853245 0.9773359 0.4475656
    ##  [8] 0.4802992 0.6065266 0.8842247 0.9266620 0.8901663 0.2295541 0.2520352
    ## [15] 0.4086021 0.4715314 0.6029402 0.5922796 0.0000000 0.1444158 0.0307866
    ## [22] 0.3935514 0.1626236 0.1413729
    ## 
    ## $predicted
    ##  [1] 0.8813256 0.8813256 0.8813256 0.8813256 0.8813256 0.8813256 0.7152553
    ##  [8] 0.7152553 0.7152553 0.7152553 0.7152553 0.7152553 0.3909324 0.3909324
    ## [15] 0.3909324 0.3909324 0.3909324 0.3909324 0.1751357 0.1751357 0.1751357
    ## [22] 0.1751357 0.1751357 0.1751357
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

    ## Number of successful fits:  999  out of 1000

``` r
print(plot_2)
```

    ## $optimFit
    ## $optimFit$par
    ## [1] 2.167896 1.575979
    ## 
    ## $optimFit$value
    ## [1] 0.5691124
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       18       18 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 54.60435
    ## 
    ## $fitAIC
    ## [1] -105.2087
    ## 
    ## $fitAICc
    ## [1] -104.9982
    ## 
    ## $fitBIC
    ## [1] -101.02
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
    ##  [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ##  [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [22] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [29] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [36] 1.0000000 0.8074037 0.7059619 0.7383556 0.9853129 0.9853245 0.9773359
    ## [43] 0.4475656 0.4802992 0.6065266 0.8842247 0.9266620 0.8901663 0.2295541
    ## [50] 0.2520352 0.4086021 0.4715314 0.6029402 0.5922796 0.0000000 0.1444158
    ## [57] 0.0307866 0.3935514 0.1626236 0.1413729
    ## 
    ## $predicted
    ##  [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ##  [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [22] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [29] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [36] 1.0000000 0.8813256 0.8813256 0.8813256 0.8813256 0.8813256 0.8813256
    ## [43] 0.7152553 0.7152553 0.7152553 0.7152553 0.7152553 0.7152553 0.3909324
    ## [50] 0.3909324 0.3909324 0.3909324 0.3909324 0.3909324 0.1751357 0.1751357
    ## [57] 0.1751357 0.1751357 0.1751357 0.1751357
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
join.dist.aic<-MergeCommNDecaydists_byCodePair(decayparam.dist=aic.dist, summ.comm_dist)
```

    ## Joining, by = "codePair"

``` r
join.dist.k<-MergeCommNDecaydists_byCodePair(decayparam.dist=k.dist, summ.comm_dist)
```

    ## Joining, by = "codePair"

``` r
join.dist.alpha<-MergeCommNDecaydists_byCodePair(decayparam.dist=alpha.dist, summ.comm_dist)
```

    ## Joining, by = "codePair"

``` r
# delta minAIC
p.aic<-ggplot(join.dist.aic, aes(x=mean_comm_dist, y=abs(decayparam_dist), color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model AIC")

# delta k
p.k<-ggplot(join.dist.k, aes(x=mean_comm_dist, y=abs(decayparam_dist), color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model k")

# delta alpha
p.alpha<-ggplot(join.dist.alpha, aes(x=mean_comm_dist, y=abs(decayparam_dist), color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model alpha")

grid.arrange(p.aic, p.k, p.alpha)
```

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-1.png)

``` r
#combine delta decay param columns
j1<-rename(join.dist.aic, "aic_dist"="decayparam_dist") %>%
  select(codePair, size, mean_comm_dist, aic_dist)
j2<-rename(join.dist.k, "k_dist"="decayparam_dist") %>%
  select(codePair, size, mean_comm_dist, k_dist)
j3<-rename(join.dist.alpha, "alpha_dist"="decayparam_dist") %>%
  select(codePair, size, mean_comm_dist, alpha_dist)

left_join(j1, j2) %>%
  left_join(j3) -> join.dist
```

    ## Joining, by = c("codePair", "size", "mean_comm_dist")

    ## Joining, by = c("codePair", "size", "mean_comm_dist")

``` r
join.dist<-mutate(join.dist, aic_dist=abs(aic_dist)) %>%
  mutate(k_dist=abs(k_dist)) %>%
  mutate(alpha_dist=abs(alpha_dist))
```

### Plot beta diversity of wood traits vs distance in decay params

``` r
#calculate pairwise wood trait distances

### LEFT OFF HERE
calc_woodTraitDist <- function(traits.mean){
  
  # calculate wood functional trait distance in multivariate space 
  
  # identify rows with no missing values
  x <- complete.cases(traits.mean[,-(1:2)]) 
  traits.mean1<-traits.mean[x,-(1:2)]
  View(traits.mean1)
  # did you lose any species doing that?
  length(unique(traits.mean$species_lower)); length(unique(traits.mean[x,]$species_lower)) #yes, lost 1 species
  #unique(meta$species)[!unique(meta$species) %in% unique(traits.mean[x,]$species)] #this species is missing
  # Olax stricta is missing because it doesn't have a waterperc value and it is only represented in small stem samples 
  
  #log-transform and scale, do PCA and take the first 3 axis, measures euc
  
  traits.scaled <- apply(log(traits.mean1+10), 2, scale)
  pc <- princomp(traits.scaled)  # stats package
  # pc=rda(traits.scaled) # vegan package
  pc.scores <- pc$scores[, 1:3]  # the first 3 axes
  
  # make a unique identifier for each row
  uniqNames<-paste(traits.mean[x,]$species, traits.mean[x,]$size, sep="_")
  row.names(pc.scores)<-uniqNames
  pc.dist.mat <- dist(pc.scores, method = "euclidean", diag=TRUE, upper=TRUE)
  mat.traitDist<-as.matrix(pc.dist.mat)
  
  #make it long
  traitDist.l <- extract_uniquePairDists(dist.mat=mat.traitDist) #make it long
  colnames(traitDist.l)<-c("samp1_spSize","samp2_spSize","woodTraitDist")
  #separate species and size identifiers
  x<-separate(traitDist.l, col=samp1_spSize, into=c("samp1_sp","samp1_size"), sep="_")
  xx<-separate(x, col=samp2_spSize, into=c("samp2_sp","samp2_size"), sep="_")
  xx$samp1_sp<-gsub(" ", "_", xx$samp1_sp)
  xx$samp2_sp<-gsub(" ", "_", xx$samp2_sp)#fix species names
  traitDist.l<-xx
  
  #add speciesXsize vs itself
  length(unique(traitDist.l$samp1_sp)); length(unique(traitDist.l$samp2_sp)) #missing 2 species (1 is ok)
  length(unique(c(traitDist.l$samp1_sp, traitDist.l$samp2_sp))) #missing only 1 species, ok
  allSp<-unique(c(traitDist.l$samp1_sp, traitDist.l$samp2_sp))
  sp<-allSp
  size<-unique(traitDist.l$samp1_size)
  rep(sp, times=2)
  rep(size, each=length(sp))
  sameSp<-data.frame(samp1_sp=rep(sp, times=2), 
                     samp1_size=rep(size, each=length(sp)),
                     samp2_sp=rep(sp, times=2), 
                     samp2_size=rep(size, each=length(sp)),
                     woodTraitDist=rep(0,length(rep(sp, times=2)))
  )
  traitDist.l<-rbind(traitDist.l,sameSp)
  #View(traitDist.l)
  
  traitDist <- traitDist.l
  return(traitDist)
}

#within species size class
```

### Plot with species+size beta diversity of microbial community vs AIC of neg.exponential
