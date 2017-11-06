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
source("code/otuIDs_fxns.R")
```

Load microbial community data

``` r
stemSamples<-load_stemSamples() #load stem sample meta data
write.csv(stemSamples, "derived_data/stemSamples.csv")

fung.otu<-load_matotu() #load the fungal OTU table
comm.otu<-add_oomycetes(fung.otu) #add the oomycetes
write.csv(comm.otu, "derived_data/comm_otu.csv")

seqSamples<-load_seqSamples(comm.otu, stemSamples) #create sequence sample meta data table
write.csv(seqSamples, "derived_data/seqSamples.csv")

#plot_sampleEffortCurves(comm.otu)

#load taxon lookup info
taxAndFunguild<-load_TaxAndFunguild(comm.otu)
```

Load wood trait data

``` r
traits.mean<-mergeTraitData()
write.csv(traits.mean, "derived_data/traits_mean.csv")

#missing data
traits.long<-as.data.frame(gather(traits.mean, key=trait, value=value, -(1:3)))
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

Load mass loss data

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

| unique  | species | size  |  time|  totalSampleDryMass| density | fruiting | insects | drill |  bark\_density|  xylem\_density|  total\_density| notes                                  |
|:--------|:--------|:------|-----:|-------------------:|:--------|:---------|:--------|:------|--------------:|---------------:|---------------:|:---------------------------------------|
| ALLI311 | ALLI    | large |    37|                  NA| NA      |          | NA      | no    |             NA|              NA|              NA| missing from plot -- missing           |
| baae1a  | baae    | small |    37|                  NA| NA      |          | NA      | no    |             NA|              NA|              NA| missing from plot -- missing from plot |

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

| unique  | species | size  | notes                                  |   37|
|:--------|:--------|:------|:---------------------------------------|----:|
| ALLI311 | ALLI    | large | missing from plot -- missing           |   NA|
| baae1a  | baae    | small | missing from plot -- missing from plot |   NA|

``` r
#remove NA rows
plotting_df %>%
  filter(!is.na(pmr)) -> plotting_df
```

Non-linear curve fits of decay trajectories Using `litterfitter` to apply both negative exponenial and weibull to all species/size classes

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

### Correlation BETWEEN species+size dissimilarities in (a) initial microbial community composition, (b) initial wood traits, and (b) decay trajectory params

``` r
# community composition (bray) vs delta k
pList<-Make_commDistvDist_Fig(distType="bray", valueCol_vec=c("ne.r2", "k","alpha"), 
                              seqSamples, stemSamples, comm.otu, spdf)
p.commVdecay.k<-pList[['k']] + guides(color=FALSE)

# wood trait vs delta k
pList<-Make_woodTraitDistvDist_Fig(valueCol_vec=c("ne.r2","k","alpha"), seqSamples, stemSamples, traits.mean, spdf)
p.traitVdecay.k<-pList[['k']] + guides(color=FALSE)

# wood trait vs community composition (bray)
traits.dist<-Calc_woodTraitDist(traits.mean)
comm.dist<-Calc_commDists(seqSamples, comm.otu, distType="bray")
summ.comm_dist<-SummarizeCommDist_byCodePair(comm.dist)
trait_comm.dist<-MergeCommNWoodTraitdists_byCodePair(traits.dist, summ.comm_dist)
p.commVtrait<-ggplot(trait_comm.dist, aes(x=mean_comm_dist, y=woodTraitDist, color=size)) + geom_point() + 
  guides(color=FALSE) + xlab("Microb comm distance (bray)")

# decay param == k
grid.arrange(p.commVdecay.k, p.traitVdecay.k, p.commVtrait, ncol=2)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png) *Check out code/initialDist\_vs\_decayDist\_btwCode.Rmd for full set of plots.*

H1. Greater average distances between species+size initial microbial communitiy compositions will lead to greater differences in decay model fit (r2), rate (k -- shown here), and lagginess (alpha).
**There's not a lot of explanatory power because all microbial compositions are really different between species+size**

H2. Greater average distances between species+size initial wood traits will lead to greater differences in decay model fit (r2), rate (k -- shown here), and lagginess (alpha).
**Looks like there's some evidence to support this hypothesis**

H3. Greater average distances between species+size initial wood traits will correspond with greater distances between species+size initial microbial community compositions.
**All microbial compositions are really different between species+size. Variation in these composition differences do not appear to correspond to wood trait distances**

### Plot richness of key players vs decay param distances BETWEEN species+size

``` r
# summarize the presence of ... in each sample
sapro.df<-Calc_richOTUtype(colNam="Trophic.Mode", grepTerm="Sapro", taxAndFunguild, comm.otu)  
sapList<-Plot_richOTUtype(otutype.df=sapro.df, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Saprotroph", spdf)

basidio.df<-Calc_richOTUtype(colNam="phylum", grepTerm="Basid", taxAndFunguild, comm.otu)
basidList<-Plot_richOTUtype(otutype.df=basidio.df, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Basidio", spdf)

path.df<-Calc_richOTUtype(colNam="Trophic.Mode", grepTerm="Patho", taxAndFunguild, comm.otu)
pathList<-Plot_richOTUtype(otutype.df=path.df, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Pathogen", spdf)

oomy.df<-Calc_richOTUtype(colNam="kingdom", grepTerm="Protist", taxAndFunguild, comm.otu)
ooList<-Plot_richOTUtype(otutype.df=oomy.df, 
                        valueCol_vec=c("ne.r2", "k","alpha"), 
                        otutypeNam="Oomycete", spdf)

#plot
grid.arrange(sapList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             sapList[['k']] + guides(color=FALSE, shape=FALSE), 
             sapList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             basidList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             basidList[['k']] + guides(color=FALSE, shape=FALSE), 
             basidList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             pathList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             pathList[['k']] + guides(color=FALSE, shape=FALSE), 
             pathList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             ooList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             ooList[['k']] + guides(color=FALSE, shape=FALSE), 
             ooList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             ncol=3)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png)

H1a. Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
H1b. Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another. **Maybe there's some indication that more saprotroph richness leads to better-fitting decay models (ne.r2)**

H2. Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers. **Maybe there's some indicaiton that more oomycete richness/presence leads to slower (k) and more laggy (alpha) decay**

### Plot frequency of positive/negative taxa 'interactions' vs decay param distances BETWEEN species+size

``` r
#identify the negative and positive residually-correlated OTU pairs in a given community
#this takes forever.... just read in the derived-data file unless input data changes
#IDcorrelatedOTUsinEachSampl(taxaAndFunguild, comm.otu)
pairsPresent.df<-read.csv(file="derived_data/pairsPresent.csv", row.names=1)


#merge the pairsPresent dataframes with the decay params data and plot
posList<-Plot_signifCor(sign='pos', 
                        valueCol_vec=c("ne.r2", "k", "alpha"), 
                        pairsPresent.df, spdf, seqSamples)
negList<-Plot_signifCor(sign='neg', 
                        valueCol_vec=c("ne.r2", "k", "alpha"), 
                        pairsPresent.df, spdf, seqSamples)

grid.arrange(posList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             posList[['k']] + guides(color=FALSE, shape=FALSE), 
             posList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             negList[['ne.r2']] + guides(color=FALSE, shape=FALSE), 
             negList[['k']] + guides(color=FALSE, shape=FALSE), 
             negList[['alpha']] + guides(color=FALSE, shape=FALSE),
             
             ncol=3)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png)

H1. A higher frequency of postively-correlated OTUs (after controlling for wood trait environmental influences) in a sample community will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha). This could be because these OTUs facilitate each other thereby making it more difficult for later-colonizing OTUs that specialize in decay to establish. \[Alternatively, a high frequency of positively-correlated taxa may indicate that these samples are dominated by co-occurring taxa that are in combat with one another...\] **No pattern apparent**

H2. A higher frequency of negatively-correlated OTUs (after controlling for wood trait environmental influences) in a sample community will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha). This could be because the occurance of negatively-correlated OTUs in the same sample indicates that there are negative biotic interactions among members of the initial community, thereby making it easier for later-colonizing OTUs that specialize in decay to establish. **No pattern apparent**

### Plot microbial beta diversity WITHIN species+size vs R2 of neg.exponential

``` r
#calculate pairwise community distances WITHIN code
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
comm.dist.wth<-left_join(comm.dist.wth, spdf)

p.ne.r2<-ggplot(comm.dist.wth, aes(x=mean, y=ne.r2, color=species, shape=size)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlab("Microbial community distance within Code") + 
  ylab("Negative exponential R2") 

p.w.r2<-ggplot(comm.dist.wth, aes(x=mean, y=w.r2, color=species, shape=size)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlab("Microbial community distance within Code") + 
  ylab("Weibull model R2") 

p.k<-ggplot(comm.dist.wth, aes(x=mean, y=k, color=species, shape=size)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlab("Microbial community distance within Code") + 
  ylab("k") 
p.k
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png)

``` r
#pdf('output/commDist_decayparam_withinCode.pdf', width=8, height=6)
grid.arrange(p.ne.r2 + guides(color=FALSE, shape=FALSE), 
             p.w.r2 + guides(color=FALSE, shape=FALSE))
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-2.png)

``` r
#dev.off()
```

H1. Expect that samples with similar initial microbial communities will have better-fitting decay models. **No pattern apparent**
