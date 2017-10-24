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
trait.means<-mergeTraitData()
```

    ## Joining, by = "SampleCode"

    ## Warning: Column `SampleCode` joining character vector and factor, coercing
    ## into character vector

    ## Joining, by = "code"

    ## Joining, by = "code"
    ## Joining, by = "code"

``` r
traits.long<-as.data.frame(gather(trait.means, key=trait, value=value, -(1:3)))

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

    ## # A tibble: 9 x 9
    ##    unique Species  size  time totalSampleDryMass density fruiting insects
    ##     <chr>   <chr> <chr> <dbl>              <dbl>   <dbl>    <chr>   <chr>
    ## 1  eute2b    eute small    37                  0     NaN             <NA>
    ## 2  acel2f    acel small    37                  0     NaN                3
    ## 3  olst1c    olst small    37                  0     NaN             <NA>
    ## 4  eusc3j    eusc small    37                  0     NaN                3
    ## 5  olst1e    olst small    37                  0     NaN             <NA>
    ## 6 ALLI311    ALLI large    37                  0     NaN             <NA>
    ## 7  hase2b    hase small    37                  0     NaN                4
    ## 8  baae1a    baae small    37                  0     NaN             <NA>
    ## 9  eute1e    eute small    37                  0     NaN                4
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

    ## Number of successful fits:  991  out of 1000 
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
    ## [1] 2.225099 1.501404
    ## 
    ## $optimFit$value
    ## [1] 0.5348948
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       15       15 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 10.61798
    ## 
    ## $fitAIC
    ## [1] -17.23595
    ## 
    ## $fitAICc
    ## [1] -16.63595
    ## 
    ## $fitBIC
    ## [1] -14.96497
    ## 
    ## $time
    ##  [1] 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 0.5833333 1.0833333
    ##  [8] 1.0833333 1.0833333 1.0833333 1.0833333 1.0833333 2.0833333 2.0833333
    ## [15] 2.0833333 2.0833333 2.0833333 2.0833333 3.0833333 3.0833333 3.0833333
    ## [22] 3.0833333 3.0833333
    ## 
    ## $mass
    ##  [1] 0.8074037 0.7059619 0.7383556 0.9853129 0.9853245 0.9773359 0.4475656
    ##  [8] 0.4802992 0.6065266 0.8842247 0.9266620 0.8901663 0.2295541 0.2520352
    ## [15] 0.4086021 0.4715314 0.6029402 0.5922796 0.1444158 0.0307866 0.3935514
    ## [22] 0.1626236 0.1413729
    ## 
    ## $predicted
    ##  [1] 0.8746090 0.8746090 0.8746090 0.8746090 0.8746090 0.8746090 0.7122151
    ##  [8] 0.7122151 0.7122151 0.7122151 0.7122151 0.7122151 0.4041833 0.4041833
    ## [15] 0.4041833 0.4041833 0.4041833 0.4041833 0.1955488 0.1955488 0.1955488
    ## [22] 0.1955488 0.1955488
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

    ## Number of successful fits:  998  out of 1000

``` r
print(plot_2)
```

    ## $optimFit
    ## $optimFit$par
    ## [1] 2.225099 1.501404
    ## 
    ## $optimFit$value
    ## [1] 0.5348948
    ## 
    ## $optimFit$counts
    ## function gradient 
    ##       36       36 
    ## 
    ## $optimFit$convergence
    ## [1] 0
    ## 
    ## $optimFit$message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## 
    ## $logLik
    ## [1] 55.02769
    ## 
    ## $fitAIC
    ## [1] -106.0554
    ## 
    ## $fitAICc
    ## [1] -105.8411
    ## 
    ## $fitBIC
    ## [1] -101.9003
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
    ## [57] 3.0833333 3.0833333 3.0833333
    ## 
    ## $mass
    ##  [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ##  [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [22] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [29] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [36] 1.0000000 0.8074037 0.7059619 0.7383556 0.9853129 0.9853245 0.9773359
    ## [43] 0.4475656 0.4802992 0.6065266 0.8842247 0.9266620 0.8901663 0.2295541
    ## [50] 0.2520352 0.4086021 0.4715314 0.6029402 0.5922796 0.1444158 0.0307866
    ## [57] 0.3935514 0.1626236 0.1413729
    ## 
    ## $predicted
    ##  [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ##  [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [15] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [22] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [29] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [36] 1.0000000 0.8746090 0.8746090 0.8746090 0.8746090 0.8746090 0.8746090
    ## [43] 0.7122151 0.7122151 0.7122151 0.7122151 0.7122151 0.7122151 0.4041832
    ## [50] 0.4041832 0.4041832 0.4041832 0.4041832 0.4041832 0.1955488 0.1955488
    ## [57] 0.1955488 0.1955488 0.1955488
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

### Plot beta diversity of microbial community vs beta diversity of decay params

``` r
#####
#calculate pairwise community distances

#1. create sample table for the community data
sampTab<-unique(mass.data[,c("Species","size")])
sampTab<-rename(sampTab, "code"="Species")
seq_sampName<-row.names(comm.otu)
indx<-data.frame(seq_sampName=seq_sampName, code=substr(seq_sampName, 1, 4))
sampTab<-left_join(indx, sampTab)
```

    ## Joining, by = "code"

    ## Warning: Column `code` joining factor and character vector, coercing into
    ## character vector

``` r
#2. calc dist, make it long format
dist <- vegdist(as.matrix(comm.otu), method = "bray", upper=TRUE, diag=TRUE)  #calculate jaccard index
dist.mat<-as.matrix(dist)
dist.df<-extract_uniquePairDists(dist.mat) #make the distance matrix long and add metadata
dist.df<-rename(dist.df, "sampID1"="sp1", "sampID2"="sp2")

#3. annotate distances with sample info
indx<-rename(sampTab, "sampID"="seq_sampName")
left_join(dist.df, indx, by=c("sampID1" = "sampID")) %>%
  rename("code1"="code",
         "size1"="size") -> dist.df1 # use the index to add info for 1st sample
```

    ## Warning: Column `sampID1`/`sampID` joining factors with different levels,
    ## coercing to character vector

``` r
left_join(dist.df1, indx, by=c("sampID2" = "sampID")) %>%
  rename("code2"="code",
         "size2"="size") -> dist.df2 # use the index to add info for the 2nd sample
```

    ## Warning: Column `sampID2`/`sampID` joining factors with different levels,
    ## coercing to character vector

``` r
#4. remove distances between size classes
comm.dist<-subset(dist.df2, size1 == size2)

#####
#calculate pairwise min AIC distances


#1. find the min AIC for each species+size
spdf<-mutate(spdf, minAIC = pmin(neg.exp.aic, w.aic) )

#2. identify spdf by code
sampTab<-mutate(sampTab, species=tolower(code))
indx<-unique(sampTab[,c("code","species","size")])
spdf<-left_join(spdf, indx) #put code into the dataframe
```

    ## Joining, by = c("species", "size")

``` r
#fix an NA code
spdf[is.na(spdf$code),"code"]<-"eusc"

#3. calc the pairwise differences in minAIC
v <- spdf$minAIC
z <- outer(v,v,'-') # sp1 - sp2 = dist
colnames(z)<-spdf$code
row.names(z)<-spdf$code
dist.df<-extract_uniquePairDists(z) #make the distance matrix long and add metadata
dist.df<-rename(dist.df, "code1"="sp1", "code2"="sp2")

#4. add back species + size identifiers
indx<-unique(sampTab[,c("code","species","size")])
left_join(dist.df, indx, by=c("code1" = "code")) %>%
  rename("species1"="species",
         "size1"="size") -> dist.df1 # use the index to add info for 1st sample
```

    ## Warning: Column `code1`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
left_join(dist.df1, indx, by=c("code2" = "code")) %>%
  rename("species2"="species",
         "size2"="size") -> dist.df2 # use the index to add info for the 2nd sample
```

    ## Warning: Column `code2`/`code` joining factor and character vector,
    ## coercing into character vector

``` r
#4. remove distances between size classes
aic.dist<-subset(dist.df2, size1 == size2)



#####
#merge community and AIC distances into the same dataframe
#community
comm.dist$codePair<-paste(comm.dist$code1, comm.dist$code2, sep="_")
comm.df<-comm.dist[,c("codePair","dist")]
comm.df<-rename(comm.df, "comm_dist"="dist")
summ.comm_dist<-group_by(comm.df, codePair) %>%
  summarize(mean=mean(comm_dist),
            se=sd(comm_dist)/sqrt(length(comm_dist)),
            lower=mean-se,
            upper=mean+se)
summ.comm_dist<-rename(summ.comm_dist, "mean_comm_dist"="mean",
                       "lower_comm_dist"="lower",
                       "upper_comm_dist"="upper")
unique(summ.comm_dist$codePair)
```

    ##   [1] "acel_acel" "acel_acpa" "acel_alli" "acel_baae" "acel_basp"
    ##   [6] "acel_cali" "acel_cota" "acel_eute" "acel_excu" "acel_hase"
    ##  [11] "acel_isan" "acel_jasc" "acel_leer" "acel_lepa" "acel_mede"
    ##  [16] "acel_olst" "acel_peli" "acel_penu" "acel_pepu" "acel_ripi"
    ##  [21] "acpa_acel" "acpa_acpa" "ACPA_ACPA" "acpa_alli" "ACPA_ALLI"
    ##  [26] "acpa_anba" "acpa_baae" "acpa_basp" "acpa_cali" "acpa_cota"
    ##  [31] "acpa_eute" "acpa_excu" "ACPA_EXCU" "acpa_hase" "acpa_isan"
    ##  [36] "acpa_jasc" "acpa_leer" "acpa_lepa" "acpa_mede" "ACPA_MEDE"
    ##  [41] "acpa_olst" "acpa_peli" "acpa_penu" "acpa_pepu" "acpa_ripi"
    ##  [46] "alli_acel" "alli_acpa" "alli_alli" "ALLI_ALLI" "alli_anba"
    ##  [51] "alli_baae" "alli_basp" "alli_cali" "alli_cota" "alli_eute"
    ##  [56] "alli_excu" "ALLI_EXCU" "alli_hase" "alli_isan" "alli_jasc"
    ##  [61] "alli_leer" "alli_lepa" "alli_mede" "ALLI_MEDE" "alli_olst"
    ##  [66] "alli_peli" "alli_penu" "alli_pepu" "alli_ripi" "anba_acel"
    ##  [71] "anba_acpa" "ANBA_ACPA" "anba_alli" "ANBA_ALLI" "anba_anba"
    ##  [76] "ANBA_ANBA" "anba_baae" "ANBA_BAAE" "anba_basp" "anba_cali"
    ##  [81] "anba_cota" "anba_eute" "anba_excu" "ANBA_EXCU" "anba_hase"
    ##  [86] "anba_isan" "anba_jasc" "ANBA_JASC" "anba_leer" "anba_lepa"
    ##  [91] "anba_mede" "ANBA_MEDE" "anba_olst" "anba_peli" "ANBA_PELI"
    ##  [96] "anba_penu" "anba_pepu" "ANBA_PEPU" "anba_ripi" "baae_acel"
    ## [101] "baae_acpa" "BAAE_ACPA" "baae_alli" "BAAE_ALLI" "baae_baae"
    ## [106] "BAAE_BAAE" "baae_basp" "baae_cali" "baae_cota" "baae_eute"
    ## [111] "baae_excu" "BAAE_EXCU" "baae_hase" "baae_isan" "baae_jasc"
    ## [116] "baae_leer" "baae_lepa" "baae_mede" "BAAE_MEDE" "baae_olst"
    ## [121] "baae_peli" "baae_penu" "baae_pepu" "BAAE_PEPU" "baae_ripi"
    ## [126] "basp_acel" "basp_acpa" "basp_alli" "basp_baae" "basp_basp"
    ## [131] "basp_cali" "basp_cota" "basp_eute" "basp_excu" "basp_hase"
    ## [136] "basp_isan" "basp_jasc" "basp_leer" "basp_lepa" "basp_mede"
    ## [141] "basp_olst" "basp_peli" "basp_penu" "basp_pepu" "basp_ripi"
    ## [146] "cali_acel" "cali_acpa" "cali_alli" "cali_baae" "cali_basp"
    ## [151] "cali_cali" "cali_cota" "cali_eute" "cali_excu" "cali_hase"
    ## [156] "cali_isan" "cali_jasc" "cali_leer" "cali_lepa" "cali_mede"
    ## [161] "cali_olst" "cali_peli" "cali_penu" "cali_pepu" "cali_ripi"
    ## [166] "cota_acel" "cota_acpa" "cota_alli" "cota_baae" "cota_basp"
    ## [171] "cota_cali" "cota_cota" "cota_eute" "cota_excu" "cota_hase"
    ## [176] "cota_isan" "cota_jasc" "cota_leer" "cota_lepa" "cota_mede"
    ## [181] "cota_olst" "cota_peli" "cota_penu" "cota_pepu" "cota_ripi"
    ## [186] "EUSC_ACPA" "EUSC_ALLI" "EUSC_ANBA" "EUSC_BAAE" "EUSC_EUSC"
    ## [191] "EUSC_EXCU" "EUSC_JASC" "EUSC_MEDE" "EUSC_PELI" "EUSC_PEPU"
    ## [196] "eute_acel" "eute_acpa" "EUTE_ACPA" "eute_alli" "EUTE_ALLI"
    ## [201] "EUTE_ANBA" "eute_baae" "EUTE_BAAE" "eute_basp" "eute_cali"
    ## [206] "eute_cota" "EUTE_EUSC" "eute_eute" "EUTE_EUTE" "eute_excu"
    ## [211] "EUTE_EXCU" "eute_hase" "eute_isan" "eute_jasc" "EUTE_JASC"
    ## [216] "eute_leer" "eute_lepa" "EUTE_LEPA" "eute_mede" "EUTE_MEDE"
    ## [221] "eute_olst" "eute_peli" "EUTE_PELI" "eute_penu" "eute_pepu"
    ## [226] "EUTE_PEPU" "eute_ripi" "excu_acel" "excu_acpa" "excu_alli"
    ## [231] "excu_baae" "excu_basp" "excu_cali" "excu_cota" "excu_eute"
    ## [236] "excu_excu" "EXCU_EXCU" "excu_hase" "excu_isan" "excu_jasc"
    ## [241] "excu_leer" "excu_lepa" "excu_mede" "EXCU_MEDE" "excu_olst"
    ## [246] "excu_peli" "excu_penu" "excu_pepu" "excu_ripi" "hase_acel"
    ## [251] "hase_acpa" "hase_alli" "hase_baae" "hase_basp" "hase_cali"
    ## [256] "hase_cota" "hase_eute" "hase_excu" "hase_hase" "hase_isan"
    ## [261] "hase_jasc" "hase_leer" "hase_lepa" "hase_mede" "hase_olst"
    ## [266] "hase_peli" "hase_penu" "hase_pepu" "hase_ripi" "isan_acel"
    ## [271] "isan_acpa" "isan_alli" "isan_baae" "isan_basp" "isan_cali"
    ## [276] "isan_cota" "isan_eute" "isan_excu" "isan_hase" "isan_isan"
    ## [281] "isan_jasc" "isan_leer" "isan_lepa" "isan_mede" "isan_olst"
    ## [286] "isan_peli" "isan_penu" "isan_pepu" "isan_ripi" "jasc_acel"
    ## [291] "jasc_acpa" "JASC_ACPA" "jasc_alli" "JASC_ALLI" "jasc_baae"
    ## [296] "JASC_BAAE" "jasc_basp" "jasc_cali" "jasc_cota" "jasc_eute"
    ## [301] "jasc_excu" "JASC_EXCU" "jasc_hase" "jasc_isan" "jasc_jasc"
    ## [306] "JASC_JASC" "jasc_leer" "jasc_lepa" "jasc_mede" "JASC_MEDE"
    ## [311] "jasc_olst" "jasc_peli" "jasc_penu" "jasc_pepu" "JASC_PEPU"
    ## [316] "jasc_ripi" "leer_acel" "leer_acpa" "leer_alli" "leer_baae"
    ## [321] "leer_basp" "leer_cali" "leer_cota" "leer_eute" "leer_excu"
    ## [326] "leer_hase" "leer_isan" "leer_jasc" "leer_leer" "leer_lepa"
    ## [331] "leer_mede" "leer_olst" "leer_peli" "leer_penu" "leer_pepu"
    ## [336] "leer_ripi" "lepa_acel" "lepa_acpa" "LEPA_ACPA" "lepa_alli"
    ## [341] "LEPA_ALLI" "LEPA_ANBA" "lepa_baae" "LEPA_BAAE" "lepa_basp"
    ## [346] "lepa_cali" "lepa_cota" "LEPA_EUSC" "lepa_eute" "lepa_excu"
    ## [351] "LEPA_EXCU" "lepa_hase" "lepa_isan" "lepa_jasc" "LEPA_JASC"
    ## [356] "lepa_leer" "lepa_lepa" "LEPA_LEPA" "lepa_mede" "LEPA_MEDE"
    ## [361] "lepa_olst" "lepa_peli" "LEPA_PELI" "lepa_penu" "lepa_pepu"
    ## [366] "LEPA_PEPU" "lepa_ripi" "mede_acel" "mede_acpa" "mede_alli"
    ## [371] "mede_baae" "mede_basp" "mede_cali" "mede_cota" "mede_eute"
    ## [376] "mede_excu" "mede_hase" "mede_isan" "mede_jasc" "mede_leer"
    ## [381] "mede_lepa" "mede_mede" "MEDE_MEDE" "mede_olst" "mede_peli"
    ## [386] "mede_penu" "mede_pepu" "mede_ripi" "olst_acel" "olst_acpa"
    ## [391] "olst_alli" "olst_baae" "olst_basp" "olst_cali" "olst_cota"
    ## [396] "olst_eute" "olst_excu" "olst_hase" "olst_isan" "olst_jasc"
    ## [401] "olst_leer" "olst_lepa" "olst_mede" "olst_olst" "olst_peli"
    ## [406] "olst_penu" "olst_pepu" "olst_ripi" "peli_acel" "peli_acpa"
    ## [411] "PELI_ACPA" "peli_alli" "PELI_ALLI" "peli_baae" "PELI_BAAE"
    ## [416] "peli_basp" "peli_cali" "peli_cota" "peli_eute" "peli_excu"
    ## [421] "PELI_EXCU" "peli_hase" "peli_isan" "peli_jasc" "PELI_JASC"
    ## [426] "peli_leer" "peli_lepa" "peli_mede" "PELI_MEDE" "peli_olst"
    ## [431] "peli_peli" "PELI_PELI" "peli_penu" "peli_pepu" "PELI_PEPU"
    ## [436] "peli_ripi" "penu_acel" "penu_acpa" "penu_alli" "penu_baae"
    ## [441] "penu_basp" "penu_cali" "penu_cota" "penu_eute" "penu_excu"
    ## [446] "penu_hase" "penu_isan" "penu_jasc" "penu_leer" "penu_lepa"
    ## [451] "penu_mede" "penu_olst" "penu_peli" "penu_penu" "penu_pepu"
    ## [456] "penu_ripi" "pepu_acel" "pepu_acpa" "PEPU_ACPA" "pepu_alli"
    ## [461] "PEPU_ALLI" "pepu_baae" "pepu_basp" "pepu_cali" "pepu_cota"
    ## [466] "pepu_eute" "pepu_excu" "PEPU_EXCU" "pepu_hase" "pepu_isan"
    ## [471] "pepu_jasc" "pepu_leer" "pepu_lepa" "pepu_mede" "PEPU_MEDE"
    ## [476] "pepu_olst" "pepu_peli" "pepu_penu" "pepu_pepu" "PEPU_PEPU"
    ## [481] "pepu_ripi" "ripi_acel" "ripi_acpa" "ripi_alli" "ripi_baae"
    ## [486] "ripi_basp" "ripi_cali" "ripi_cota" "ripi_eute" "ripi_excu"
    ## [491] "ripi_hase" "ripi_isan" "ripi_jasc" "ripi_leer" "ripi_lepa"
    ## [496] "ripi_mede" "ripi_olst" "ripi_peli" "ripi_penu" "ripi_pepu"
    ## [501] "ripi_ripi"

``` r
#aic
aic.dist$codePair<-paste(aic.dist$code1, aic.dist$code2, sep="_")
aic.dist$codePair_rev<-paste(aic.dist$code2, aic.dist$code1, sep="_")
aic.df<-aic.dist[,c("codePair","codePair_rev","dist")]
aic.df<-rename(aic.df, "aic_dist"="dist")
aic.df_forward<-aic.df[,c("codePair","aic_dist")]
aic.df_reverse<-aic.df[,c("codePair_rev","aic_dist")]
uniqAICcodePairs<-unique(c(aic.df$codePair, aic.df$codePair_rev))

#why are all these code pairs missing from aic.df?
missingCodePairs<-summ.comm_dist$codePair[!unique(summ.comm_dist$codePair) %in% uniqAICcodePairs]
# df<-data.frame(codePair=missingCodePairs)
#  df %>%
#   separate(codePair, into=c("code1","code2"), remove=FALSE) %>%
#   filter(code1 != code2) #1. aic does not include code pairs within the same species+size

summ.comm_dist %>%
  left_join(aic.df_forward) %>% #join with the forward versions of aic.df$codePair
  rename("aic_dist_forward"="aic_dist") %>%
  left_join(aic.df_reverse, by=c("codePair"="codePair_rev")) %>%
  rename("aic_dist_reverse"="aic_dist") -> summ.comm_dist
```

    ## Joining, by = "codePair"

``` r
filter(summ.comm_dist, is.na(aic_dist_forward) & !is.na(aic_dist_reverse)) %>%
  mutate(aic_dist=aic_dist_reverse) -> reverse.rows
filter(summ.comm_dist, !is.na(aic_dist_forward) & is.na(aic_dist_reverse)) %>%
  mutate(aic_dist=aic_dist_forward) -> forward.rows
summ.comm_dist<-bind_rows(reverse.rows, forward.rows)
summ.comm_dist<-summ.comm_dist[,c("codePair","mean_comm_dist","lower_comm_dist","upper_comm_dist","aic_dist")]

#add size back
summ.comm_dist$size<-"large"
summ.comm_dist[tolower(summ.comm_dist$codePair) == summ.comm_dist$codePair,"size"]<-"small"

ggplot(summ.comm_dist, aes(x=mean_comm_dist, y=abs(aic_dist), color=size)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lower_comm_dist, xmax=upper_comm_dist)) +
  xlab('Mean microbial community distance') + ylab("Delta decay model AIC")
```

    ## Warning: Removed 5 rows containing missing values (geom_errorbarh).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-12-1.png)
