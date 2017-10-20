Does chemistry or community better predict mass loss?
================
Marissa Lee
10/18/2017

``` r
devtools::install_github("cornwell-lab-unsw/litterfitter")
source("code/read_initial.R")
source("code/functions.R")
library(dplyr)
library(ggplot2)
library(readr)
library(vegan)
library(litterfitter)
```

Load harvest data

``` r
knitr::opts_chunk$set(echo = TRUE)

#initial mass
initial_mass <- read_in_initial_mass()

initial_mass %>%
  mutate(SpeciesCode=tolower(Species))%>%
  ggplot(aes(y=totalSampleDryMass,x=SpeciesCode,fill=size))+
  geom_violin()+ 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + 
  scale_y_log10()
```

    ## Warning: Removed 33 rows containing non-finite values (stat_ydensity).

![](readme_files/figure-markdown_github-ascii_identifiers/setup-1.png)

``` r
#mass at harvests
harvest_data<-LoadHarvestFiles()
```

    ## Warning in read.samp3(): NAs introduced by coercion

``` r
harvest_data1<-CalcTotalDryMass(harvest_data)
harvest_data2<-CalcDensity(harvest_data1)
harvest_mass<-ReorgDataFrame(harvest_data2)
```

A first pass at plotting to look for outliers

``` r
#create a complete sample mass df for all time points
mass.data<-rbind(initial_mass, harvest_mass)
mass.data %>%
  ggplot(aes(x=time, y=totalSampleDryMass)) + geom_point(alpha=0.6)+theme_bw() 
```

    ## Warning: Removed 34 rows containing missing values (geom_point).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-2-1.png)

``` r
#funky outlier
max(mass.data$totalSampleDryMass, na.rm=TRUE)
```

    ## [1] 1242.641

``` r
mass.data[which.max(mass.data$totalSampleDryMass),]
```

    ## # A tibble: 1 x 9
    ##    unique Species  size  time totalSampleDryMass   density
    ##     <chr>   <chr> <chr> <dbl>              <dbl>     <dbl>
    ## 1 ALLI111    ALLI large    37           1242.641 0.7351301
    ## # ... with 3 more variables: fruiting <chr>, insects <chr>, drill <chr>

and another view. Might want to check out those two high mass value outliers from harvest 3.

``` r
mass.data %>%
  ggplot(aes(x=time, y=totalSampleDryMass,col=size)) + geom_point(position="jitter",alpha=0.6)+theme_bw()+scale_y_log10()
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 34 rows containing missing values (geom_point).

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png) Looking for outliers in the other direction, but I guess these are real 0's?

``` r
#funky outlier
min(mass.data$totalSampleDryMass, na.rm=TRUE)
```

    ## [1] 0

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

``` r
#colnames(mass.data)
#unique(mass.data$unique)
```

Merge time zero with the other harvests to calculate proportion mass remaining at each time point

``` r
mass.data %>%
  filter(time==0) %>%
  rename(timeZeroDensity=density) %>%
  rename(timeZeroMass=totalSampleDryMass) %>%
  select(unique,timeZeroMass,timeZeroDensity)->time_zero

mass.data %>%
  left_join(time_zero,by="unique") %>%
  mutate(pmr=totalSampleDryMass/timeZeroMass) %>%
  filter(pmr<2) %>% #remove outlier noted above
  mutate(SpeciesCode=tolower(Species)) -> plotting_df

plotting_df%>%
  ggplot(aes(x=time, y=pmr,col=SpeciesCode,shape=size)) + 
  geom_point()+geom_smooth(method="lm",se=FALSE) ->gg1
gg1
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png)

Interactive plotting for fun

``` r
#library(plotly)
#ggplotly(gg1)
```

An example using litterFitter to remember how it works.

``` r
plotting_df %>%
  filter(SpeciesCode=="eute",size=="small") ->out

plot_multiple_fits(time = out$time/12,
                   mass.remaining = out$pmr,
                   bty = 'n', model = c('neg.exp', 'weibull'),
                   xlab = 'Time', ylab = 'Proportion mass remaining',iters=1000)
```

    ## Number of successful fits:  994  out of 1000 
    ## Number of successful fits:  1000  out of 1000

    ## Warning in multioptimFit(time, mass.remaining, model, iters = iters, upper
    ## = upper, : May not have found global best fit; increase iterations

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png)

and looking at branches seperately

``` r
plotting_df %>%
  mutate(branch=as.factor(substr(plotting_df$unique, 5, 5))) %>%
   filter(SpeciesCode=="eute",size=="large") %>%
    ggplot(aes(x=time, y=pmr,col=branch,shape=size))+geom_point()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png)

``` r
  geom_smooth(se=FALSE) ->gg1
```

applying to all species/size classes

``` r
plotting_df %>%
  mutate(sp_size=paste(SpeciesCode,size,sep="_")) ->out

  fits <- lapply(split(out,factor(out$sp_size)),function(x){
  fit_litter(time = x$time/12, 
        mass.remaining = x$pmr, model = c("neg.exp"), iters = 1000)
  })
  
  spdf<-data.frame(k=unlist(lapply(fits,function(x)x$optimFit$par)))
  spdf$sp_size<-row.names(spdf)
  
  spdf %>%
    separate(sp_size, c("species", "size"))%>%
    ggplot(aes(x=size,y=k,col=species))+geom_point()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png)

``` r
spdf$sp_size<-row.names(spdf)
  
  
  w.fits <- lapply(split(out,factor(out$sp_size)),function(x){
  fit_litter(time = x$time/12, 
        mass.remaining = x$pmr, model = c("weibull"), iters = 1000)
  })
  
spdf<-data.frame(k=unlist(lapply(fits,
                                 function(x)x$optimFit$par)),
                 t70=unlist(lapply(fits,function(x)             time_to_prop_mass_remaining(x,threshold.mass=0.70))),
                 neg.exp.aic=unlist(lapply(fits,
                                 function(x)x$fitAICc)),
                 w.t70=unlist(lapply(w.fits,function(x)             time_to_prop_mass_remaining(x,threshold.mass=0.70))),
                 w.aic=unlist(lapply(w.fits,
                                 function(x)x$fitAICc)),
                 alpha=unlist(lapply(w.fits,
                                 function(x)x$optimFit$par[[2]]))#check parameter order
                 
                 )

spdf$sp_size<-rownames(spdf)


spdf %>%
  separate(sp_size,sep="_",into = c("species","size"))->spdf

ggplot(spdf,aes(x=t70,y=w.t70,col=size))+geom_point()+labs(x="Time to 30% mass loss (negative exponential)", y="weibull time to 30% mass loss")+geom_abline(slope=1,intercept=0,linetype="dashed")+theme_bw()
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-10-1.png)

``` r
write_csv(spdf,"derived_data/mass_loss_parameters.csv")

spdf[which(-spdf$neg.exp.aic+spdf$w.aic < -2),]
```

``` r
plotting_df %>%
  filter(SpeciesCode=="ripi",size=="small") ->out

plot_multiple_fits(time = out$time/12,
                   mass.remaining = out$pmr,
                   bty = 'n', model = c('neg.exp', 'weibull'),
                   xlab = 'Time', ylab = 'Proportion mass remaining',iters=1000)
```

    ## Number of successful fits:  991  out of 1000 
    ## Number of successful fits:  999  out of 1000

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png)
