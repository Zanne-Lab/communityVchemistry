Does chemistry or community better predict mass loss?
================
Marissa Lee
10/23/2017

### Load microbial community data

### Load wood trait data

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png)

########################################## 

Wood traits as a predictor
--------------------------

*Hyp:* Stem-specific initial wood traits will predict variation in percent mass loss.

``` r
#function to pull out samples with complete trait info
CreateTraitPMRpair<-function(timePoint, traits.codeStem, pmr.byStem.df.w){
  
  #make a dataframe using the current time point's pmr and remove NAs
  pmr.byStem.df.w %>%
    select_("codeStem",timePoint) %>%
    rename_("curr.time"=timePoint) %>%
    filter(!is.na(curr.time)) -> pmr.noNAs
  
  #subset the trait matrix using these unique codeStems
  traits.codeStem %>%
    filter(codeStem %in% pmr.noNAs$codeStem) -> curr.traits
  
  #most complete trait set
  curr.traits %>%
    filter(!is.na(waterperc) & !is.na(P) & !is.na(K) & !is.na(Ca) & !is.na(Mn) & !is.na(Fe) & !is.na(Zn) & !is.na(N) & !is.na(C)) -> curr.traits
  #curr.traits<-curr.traits[complete.cases(curr.traits),] #so few if you do this..

  #get rid of pmr rows for which there is missing trait data
  pmr.noNAs %>%
    filter(codeStem %in% curr.traits$codeStem) -> curr.pmr
    
  #merge the dataframes
  curr.df<-left_join(curr.pmr, curr.traits) 
  
  #add code and species and size
  curr.df<-separate(curr.df, col=codeStem, into=c("code","Stem"), sep=4, remove=FALSE)
  curr.df$species<-tolower(curr.df$code)
  curr.df$size<-"large"
  curr.df[curr.df$code == tolower(curr.df$code),"size"]<-"small"
  
  return(curr.df)
    
  }

timeList<-list('time7','time13','time25','time37')

#very few complete cases of wood traits on stem-level
#density and barkthick only on small stems
datasets<-lapply(timeList, CreateTraitPMRpair, traits.codeStem, pmr.byStem.df.w)
names(datasets)<-c('time7','time13','time25','time37')

#fit full models (no barkthick or density because these weren't measured on small stems)
mod.full.t7<-lm(curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + C, data=datasets[['time7']])
mod.full.t13<-lm(curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + C, data=datasets[['time13']])
mod.full.t25<-lm(curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + C, data=datasets[['time25']])
mod.full.t37<-lm(curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + C, data=datasets[['time37']])

#do stepwise model selection
mod.select.t7<-step(mod.full.t7, direction="backward")
```

    ## Start:  AIC=-243.89
    ## curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + 
    ##     C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - P          1  0.000018 0.39267 -245.88
    ## - Ca         1  0.000389 0.39304 -245.83
    ## - Fe         1  0.001448 0.39410 -245.69
    ## - size       1  0.003001 0.39565 -245.47
    ## - Mn         1  0.004194 0.39685 -245.31
    ## - Zn         1  0.005655 0.39831 -245.11
    ## - K          1  0.005962 0.39861 -245.07
    ## - N          1  0.009073 0.40173 -244.65
    ## - C          1  0.010164 0.40282 -244.51
    ## <none>                   0.39265 -243.89
    ## - waterperc  1  0.056020 0.44867 -238.68
    ## 
    ## Step:  AIC=-245.88
    ## curr.time ~ size + waterperc + K + Ca + Mn + Fe + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Ca         1  0.000371 0.39304 -247.83
    ## - Fe         1  0.001433 0.39410 -247.69
    ## - size       1  0.003061 0.39573 -247.46
    ## - Mn         1  0.004419 0.39709 -247.28
    ## - Zn         1  0.005994 0.39866 -247.06
    ## - K          1  0.006341 0.39901 -247.02
    ## - N          1  0.009457 0.40213 -246.60
    ## - C          1  0.010279 0.40295 -246.49
    ## <none>                   0.39267 -245.88
    ## - waterperc  1  0.057242 0.44991 -240.53
    ## 
    ## Step:  AIC=-247.83
    ## curr.time ~ size + waterperc + K + Mn + Fe + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Fe         1  0.001299 0.39434 -249.65
    ## - size       1  0.003444 0.39649 -249.36
    ## - Mn         1  0.005941 0.39898 -249.02
    ## - Zn         1  0.006283 0.39932 -248.98
    ## - K          1  0.007530 0.40057 -248.81
    ## - N          1  0.010920 0.40396 -248.35
    ## - C          1  0.013585 0.40663 -248.00
    ## <none>                   0.39304 -247.83
    ## - waterperc  1  0.061470 0.45451 -241.99
    ## 
    ## Step:  AIC=-249.65
    ## curr.time ~ size + waterperc + K + Mn + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Mn         1  0.007355 0.40169 -250.66
    ## - Zn         1  0.007746 0.40209 -250.60
    ## - K          1  0.010032 0.40437 -250.30
    ## - C          1  0.012733 0.40707 -249.94
    ## <none>                   0.39434 -249.65
    ## - size       1  0.015546 0.40989 -249.57
    ## - N          1  0.021488 0.41583 -248.79
    ## - waterperc  1  0.065194 0.45953 -243.39
    ## 
    ## Step:  AIC=-250.66
    ## curr.time ~ size + waterperc + K + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - K          1  0.008072 0.40977 -251.58
    ## - Zn         1  0.008147 0.40984 -251.57
    ## - size       1  0.010681 0.41238 -251.24
    ## - C          1  0.013147 0.41484 -250.92
    ## - N          1  0.014710 0.41640 -250.71
    ## <none>                   0.40169 -250.66
    ## - waterperc  1  0.060186 0.46188 -245.12
    ## 
    ## Step:  AIC=-251.58
    ## curr.time ~ size + waterperc + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - size       1  0.005825 0.41559 -252.82
    ## - Zn         1  0.011080 0.42085 -252.14
    ## - N          1  0.011625 0.42139 -252.07
    ## <none>                   0.40977 -251.58
    ## - C          1  0.038036 0.44780 -248.79
    ## - waterperc  1  0.053430 0.46320 -246.96
    ## 
    ## Step:  AIC=-252.82
    ## curr.time ~ waterperc + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Zn         1  0.013591 0.42918 -253.08
    ## <none>                   0.41559 -252.82
    ## - N          1  0.017199 0.43279 -252.63
    ## - C          1  0.035175 0.45077 -250.43
    ## - waterperc  1  0.048414 0.46401 -248.87
    ## 
    ## Step:  AIC=-253.08
    ## curr.time ~ waterperc + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - N          1  0.011423 0.44061 -253.66
    ## <none>                   0.42918 -253.08
    ## - C          1  0.031256 0.46044 -251.29
    ## - waterperc  1  0.045071 0.47425 -249.69
    ## 
    ## Step:  AIC=-253.66
    ## curr.time ~ waterperc + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## <none>                   0.44061 -253.66
    ## - C          1  0.033950 0.47456 -251.66
    ## - waterperc  1  0.035363 0.47597 -251.50

``` r
mod.select.t13<-step(mod.full.t13, direction="backward")
```

    ## Start:  AIC=-193.9
    ## curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + 
    ##     C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - K          1  0.000049 0.90187 -195.90
    ## - Ca         1  0.000404 0.90222 -195.88
    ## - C          1  0.001761 0.90358 -195.80
    ## - Fe         1  0.011103 0.91292 -195.25
    ## - Zn         1  0.016941 0.91876 -194.92
    ## <none>                   0.90182 -193.90
    ## - P          1  0.060108 0.96193 -192.48
    ## - size       1  0.065550 0.96737 -192.18
    ## - Mn         1  0.077968 0.97979 -191.51
    ## - N          1  0.100012 1.00183 -190.33
    ## - waterperc  1  0.280175 1.18199 -181.56
    ## 
    ## Step:  AIC=-195.9
    ## curr.time ~ size + waterperc + P + Ca + Mn + Fe + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Ca         1  0.000367 0.90223 -197.88
    ## - C          1  0.002414 0.90428 -197.76
    ## - Fe         1  0.011820 0.91369 -197.21
    ## - Zn         1  0.017711 0.91958 -196.87
    ## <none>                   0.90187 -195.90
    ## - P          1  0.068543 0.97041 -194.02
    ## - size       1  0.075385 0.97725 -193.65
    ## - Mn         1  0.079536 0.98140 -193.42
    ## - N          1  0.104874 1.00674 -192.07
    ## - waterperc  1  0.300277 1.20214 -182.67
    ## 
    ## Step:  AIC=-197.88
    ## curr.time ~ size + waterperc + P + Mn + Fe + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - C          1  0.002073 0.90431 -199.76
    ## - Fe         1  0.012035 0.91427 -199.18
    ## - Zn         1  0.017350 0.91958 -198.87
    ## <none>                   0.90223 -197.88
    ## - P          1  0.071888 0.97412 -195.81
    ## - size       1  0.075079 0.97731 -195.64
    ## - Mn         1  0.084494 0.98673 -195.13
    ## - N          1  0.106146 1.00838 -193.98
    ## - waterperc  1  0.308175 1.21041 -184.31
    ## 
    ## Step:  AIC=-199.76
    ## curr.time ~ size + waterperc + P + Mn + Fe + Zn + N
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Fe         1  0.011648 0.91595 -201.08
    ## - Zn         1  0.016865 0.92117 -200.78
    ## <none>                   0.90431 -199.76
    ## - size       1  0.074599 0.97891 -197.56
    ## - P          1  0.081694 0.98600 -197.17
    ## - Mn         1  0.084434 0.98874 -197.03
    ## - N          1  0.110341 1.01465 -195.66
    ## - waterperc  1  0.306395 1.21070 -186.29
    ## 
    ## Step:  AIC=-201.08
    ## curr.time ~ size + waterperc + P + Mn + Zn + N
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Zn         1   0.02559 0.94155 -201.62
    ## <none>                   0.91595 -201.08
    ## - P          1   0.09503 1.01099 -197.85
    ## - Mn         1   0.09993 1.01589 -197.59
    ## - N          1   0.20110 1.11705 -192.56
    ## - size       1   0.26624 1.18220 -189.56
    ## - waterperc  1   0.31820 1.23415 -187.28
    ## 
    ## Step:  AIC=-201.62
    ## curr.time ~ size + waterperc + P + Mn + N
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## <none>                   0.94155 -201.62
    ## - P          1  0.083191 1.02474 -199.13
    ## - Mn         1  0.099404 1.04095 -198.30
    ## - N          1  0.178484 1.12003 -194.42
    ## - size       1  0.281754 1.22330 -189.74
    ## - waterperc  1  0.308982 1.25053 -188.58

``` r
mod.select.t25<-step(mod.full.t25, direction="backward")
```

    ## Start:  AIC=-211.1
    ## curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + 
    ##     C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - C          1   0.00160 0.58940 -212.96
    ## - Ca         1   0.00712 0.59492 -212.47
    ## - K          1   0.00850 0.59630 -212.35
    ## - N          1   0.01762 0.60542 -211.56
    ## - Mn         1   0.01981 0.60761 -211.37
    ## <none>                   0.58780 -211.10
    ## - Fe         1   0.02320 0.61100 -211.08
    ## - Zn         1   0.02482 0.61262 -210.95
    ## - P          1   0.07163 0.65943 -207.12
    ## - size       1   0.07740 0.66520 -206.66
    ## - waterperc  1   0.34091 0.92871 -189.31
    ## 
    ## Step:  AIC=-212.95
    ## curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Ca         1   0.00564 0.59503 -214.46
    ## - K          1   0.01324 0.60264 -213.80
    ## - N          1   0.01884 0.60823 -213.32
    ## - Mn         1   0.01887 0.60826 -213.32
    ## <none>                   0.58940 -212.96
    ## - Zn         1   0.02360 0.61300 -212.91
    ## - Fe         1   0.02407 0.61346 -212.87
    ## - P          1   0.07074 0.66013 -209.06
    ## - size       1   0.07951 0.66891 -208.37
    ## - waterperc  1   0.34120 0.93059 -191.21
    ## 
    ## Step:  AIC=-214.46
    ## curr.time ~ size + waterperc + P + K + Mn + Fe + Zn + N
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - K          1   0.00916 0.60419 -215.67
    ## - Mn         1   0.01425 0.60929 -215.23
    ## - N          1   0.01453 0.60957 -215.21
    ## - Fe         1   0.02113 0.61616 -214.65
    ## - Zn         1   0.02215 0.61719 -214.56
    ## <none>                   0.59503 -214.46
    ## - P          1   0.06513 0.66017 -211.06
    ## - size       1   0.07398 0.66901 -210.37
    ## - waterperc  1   0.33914 0.93417 -193.01
    ## 
    ## Step:  AIC=-215.67
    ## curr.time ~ size + waterperc + P + Mn + Fe + Zn + N
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - N          1   0.01061 0.61479 -216.76
    ## - Mn         1   0.01081 0.61500 -216.74
    ## - Fe         1   0.01570 0.61989 -216.33
    ## <none>                   0.60419 -215.67
    ## - Zn         1   0.02643 0.63062 -215.44
    ## - size       1   0.06550 0.66969 -212.31
    ## - P          1   0.10153 0.70572 -209.59
    ## - waterperc  1   0.33197 0.93616 -194.90
    ## 
    ## Step:  AIC=-216.76
    ## curr.time ~ size + waterperc + P + Mn + Fe + Zn
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Mn         1   0.00386 0.61866 -218.44
    ## - Fe         1   0.00724 0.62203 -218.15
    ## - Zn         1   0.01883 0.63363 -217.19
    ## <none>                   0.61479 -216.76
    ## - size       1   0.05514 0.66994 -214.29
    ## - P          1   0.09108 0.70588 -211.58
    ## - waterperc  1   0.34787 0.96267 -195.44
    ## 
    ## Step:  AIC=-218.44
    ## curr.time ~ size + waterperc + P + Fe + Zn
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Fe         1   0.00815 0.62681 -219.75
    ## - Zn         1   0.02147 0.64013 -218.66
    ## <none>                   0.61866 -218.44
    ## - size       1   0.05171 0.67037 -216.26
    ## - P          1   0.08872 0.70738 -213.47
    ## - waterperc  1   0.37301 0.99167 -195.90
    ## 
    ## Step:  AIC=-219.75
    ## curr.time ~ size + waterperc + P + Zn
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Zn         1   0.01947 0.64628 -220.16
    ## <none>                   0.62681 -219.75
    ## - size       1   0.05480 0.68161 -217.40
    ## - P          1   0.08547 0.71228 -215.11
    ## - waterperc  1   0.40065 1.02746 -196.06
    ## 
    ## Step:  AIC=-220.16
    ## curr.time ~ size + waterperc + P
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## <none>                   0.64628 -220.16
    ## - size       1   0.05763 0.70391 -217.72
    ## - P          1   0.07810 0.72438 -216.23
    ## - waterperc  1   0.40299 1.04927 -196.96

``` r
mod.select.t37<-step(mod.full.t37, direction="backward")
```

    ## Start:  AIC=-195.11
    ## curr.time ~ size + waterperc + P + K + Ca + Mn + Fe + Zn + N + 
    ##     C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Fe         1   0.00536 0.72772 -196.73
    ## - N          1   0.01048 0.73285 -196.38
    ## - K          1   0.01449 0.73686 -196.10
    ## - Ca         1   0.02338 0.74575 -195.49
    ## - Zn         1   0.02646 0.74883 -195.27
    ## <none>                   0.72237 -195.11
    ## - Mn         1   0.03937 0.76173 -194.40
    ## - C          1   0.04210 0.76447 -194.22
    ## - size       1   0.06412 0.78649 -192.77
    ## - P          1   0.09567 0.81804 -190.77
    ## - waterperc  1   0.65231 1.37468 -164.29
    ## 
    ## Step:  AIC=-196.73
    ## curr.time ~ size + waterperc + P + K + Ca + Mn + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - K          1   0.01129 0.73901 -197.95
    ## - Ca         1   0.02137 0.74909 -197.26
    ## <none>                   0.72772 -196.73
    ## - N          1   0.02946 0.75718 -196.71
    ## - Zn         1   0.03451 0.76224 -196.37
    ## - C          1   0.04598 0.77370 -195.61
    ## - Mn         1   0.04605 0.77377 -195.60
    ## - P          1   0.10355 0.83127 -191.95
    ## - size       1   0.21432 0.94205 -185.57
    ## - waterperc  1   0.67829 1.40601 -165.15
    ## 
    ## Step:  AIC=-197.95
    ## curr.time ~ size + waterperc + P + Ca + Mn + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## - Ca         1   0.01523 0.75424 -198.91
    ## - Zn         1   0.02917 0.76818 -197.97
    ## <none>                   0.73901 -197.95
    ## - N          1   0.04293 0.78194 -197.07
    ## - Mn         1   0.05794 0.79695 -196.10
    ## - C          1   0.06759 0.80660 -195.48
    ## - P          1   0.09346 0.83248 -193.87
    ## - size       1   0.32624 1.06525 -181.30
    ## - waterperc  1   0.85388 1.59289 -160.78
    ## 
    ## Step:  AIC=-198.91
    ## curr.time ~ size + waterperc + P + Mn + Zn + N + C
    ## 
    ##             Df Sum of Sq     RSS     AIC
    ## <none>                   0.75424 -198.91
    ## - Zn         1   0.03638 0.79062 -198.50
    ## - N          1   0.05045 0.80470 -197.60
    ## - C          1   0.05309 0.80733 -197.44
    ## - Mn         1   0.06984 0.82408 -196.39
    ## - P          1   0.14515 0.89940 -191.93
    ## - size       1   0.32876 1.08300 -182.46
    ## - waterperc  1   0.87901 1.63326 -161.50

``` r
#save the residuals
traitResid.t7<- data.frame(codeStem = datasets[['time7']]$codeStem, trait.resid = mod.select.t7$residuals)
traitResid.t13<- data.frame(codeStem = datasets[['time13']]$codeStem, trait.resid = mod.select.t13$residuals)
traitResid.t25<- data.frame(codeStem = datasets[['time25']]$codeStem, trait.resid = mod.select.t25$residuals)
traitResid.t37<- data.frame(codeStem = datasets[['time37']]$codeStem, trait.resid = mod.select.t37$residuals)
```

### time7

**less water content and C leads to more mass remaining after 7 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ waterperc + C, data = datasets[["time7"]])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.23820 -0.05588  0.01033  0.07087  0.19121 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.962933   0.493400   3.978  0.00022 ***
    ## waterperc   -0.003632   0.001795  -2.023  0.04831 *  
    ## C           -0.018543   0.009354  -1.982  0.05284 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09295 on 51 degrees of freedom
    ## Multiple R-squared:  0.1271, Adjusted R-squared:  0.0929 
    ## F-statistic: 3.714 on 2 and 51 DF,  p-value: 0.0312

### time 13

**larger size stems, less water content, more P and Mn, and less N leads to more mass remaining after 13 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ size + waterperc + P + Mn + N, data = datasets[["time13"]])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.277899 -0.092108  0.009127  0.079788  0.297114 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.4304321  0.1740140   8.220  1.2e-10 ***
    ## sizesmall   -0.1720630  0.0458801  -3.750 0.000484 ***
    ## waterperc   -0.0135897  0.0034603  -3.927 0.000280 ***
    ## P            0.0005436  0.0002667   2.038 0.047219 *  
    ## Mn           0.0006215  0.0002790   2.228 0.030735 *  
    ## N           -0.4242546  0.1421344  -2.985 0.004492 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1415 on 47 degrees of freedom
    ## Multiple R-squared:  0.363,  Adjusted R-squared:  0.2952 
    ## F-statistic: 5.357 on 5 and 47 DF,  p-value: 0.0005615

### time 25

**larger size stems, less water content, more P leads to more mass remaining after 25 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ size + waterperc + P, data = datasets[["time25"]])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.258325 -0.059630 -0.005782  0.068751  0.261561 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.1953820  0.1176559  10.160 1.50e-13 ***
    ## sizesmall   -0.0721459  0.0348728  -2.069   0.0440 *  
    ## waterperc   -0.0128901  0.0023561  -5.471 1.59e-06 ***
    ## P            0.0005140  0.0002134   2.408   0.0199 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.116 on 48 degrees of freedom
    ## Multiple R-squared:  0.4107, Adjusted R-squared:  0.3738 
    ## F-statistic: 11.15 on 3 and 48 DF,  p-value: 1.139e-05

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png)

### time 37

**larger size stems, less water content, more P leads to more mass remaining after 37 months**

    ## 
    ## Call:
    ## lm(formula = curr.time ~ size + waterperc + P + Mn + Zn + N + 
    ##     C, data = datasets[["time37"]])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.33992 -0.08112  0.01198  0.08868  0.27249 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.3343324  0.7751083   0.431  0.66838    
    ## sizesmall   -0.1924388  0.0444506  -4.329 8.78e-05 ***
    ## waterperc   -0.0236176  0.0033363  -7.079 9.89e-09 ***
    ## P            0.0007591  0.0002639   2.877  0.00623 ** 
    ## Mn           0.0006056  0.0003035   1.995  0.05236 .  
    ## Zn          -0.0007677  0.0005331  -1.440  0.15708    
    ## N           -0.2447862  0.1443333  -1.696  0.09712 .  
    ## C            0.0255493  0.0146858   1.740  0.08906 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1324 on 43 degrees of freedom
    ## Multiple R-squared:  0.6233, Adjusted R-squared:  0.562 
    ## F-statistic: 10.16 on 7 and 43 DF,  p-value: 1.968e-07

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-8-1.png)

*Hyp:* Variation in wood traits will lead to differences in decay model fit (r2), rate (k), and lagginess (alpha). Specifically, we expect samples with (a) high waterperc, (b) low density and C, (c) high P, K, Ca, Mn, Fe, Zn, and N, and (d) thicker bark (potential mech: limiting microbial colonization) to have better-fiting decay models (r2), faster decay rates (k), and less lagginess (alpha).

### r2

**greater water content and greater Zn and N leads to better-fitting decay models**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, the best model that greater water content and less C leads to better-fitting decay models

    ## 
    ## Call:
    ## lm(formula = ne.r2 ~ waterperc + barkthick + Ca + Zn + N + C, 
    ##     data = spdf.traits)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11753 -0.02941  0.00184  0.02233  0.12797 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.278e+00  6.457e-01   1.979   0.0590 .  
    ## waterperc    7.405e-03  1.509e-03   4.906 4.75e-05 ***
    ## barkthick   -3.098e-02  1.587e-02  -1.953   0.0622 .  
    ## Ca          -1.189e-05  6.783e-06  -1.753   0.0919 .  
    ## Zn           1.504e-03  5.664e-04   2.654   0.0136 *  
    ## N            2.092e-01  8.436e-02   2.480   0.0202 *  
    ## C           -1.698e-02  1.291e-02  -1.315   0.2006    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06306 on 25 degrees of freedom
    ## Multiple R-squared:  0.559,  Adjusted R-squared:  0.4532 
    ## F-statistic: 5.282 on 6 and 25 DF,  p-value: 0.001234

### k

**small size stems, greater water content, thinner bark, less Ca, more Zn, and more N lead to faster decay**

NOTE from Will: Density explains the same part of the variation in decay rates that initial water content does, only less well. (In other words, although, density gets dropped from the best model by the model selection procedure, if we remove initial water from consideration entirely, density is included in the model as the best predictor.)

So my current interpretation is that wood water rentention--related to fiber saturation point and partially captured by the density measurement--has a strong effect on long-term decomposition rates, possibly by maintaining fungal activity further into dry periods. There is also a very likely interaction between this water retention capacity with the fungal community (see results in Setting the Stage paper, Lee et al. in review).

    ## 
    ## Call:
    ## lm(formula = k ~ size + waterperc + barkthick + Ca + Zn + N + 
    ##     C, data = spdf.traits)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.092617 -0.034671  0.007271  0.029023  0.123640 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  6.021e-01  5.754e-01   1.046 0.305851    
    ## sizesmall    9.050e-02  2.405e-02   3.763 0.000958 ***
    ## waterperc    8.990e-03  1.378e-03   6.524 9.53e-07 ***
    ## barkthick   -3.862e-02  1.406e-02  -2.748 0.011198 *  
    ## Ca          -1.869e-05  5.920e-06  -3.157 0.004263 ** 
    ## Zn           1.811e-03  4.922e-04   3.680 0.001177 ** 
    ## N            2.334e-01  7.674e-02   3.041 0.005628 ** 
    ## C           -1.475e-02  1.135e-02  -1.300 0.206001    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05446 on 24 degrees of freedom
    ## Multiple R-squared:  0.7486, Adjusted R-squared:  0.6752 
    ## F-statistic: 10.21 on 7 and 24 DF,  p-value: 6.942e-06

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-1.png)![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-11-2.png)

### t70

**small stem sizes, less water content, thicker bark, more Ca, less Zn, and less N lead to longer wood "70%"-lives**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, the best model indicated that large size stems, less water content, more Ca, and less Zn lead to longer wood "70%"-lives

    ## 
    ## Call:
    ## lm(formula = t70 ~ size + waterperc + barkthick + Ca + Zn + N, 
    ##     data = spdf.traits)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.55548 -0.15857 -0.00642  0.07054  0.58357 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.799e+00  4.125e-01   9.209 1.65e-09 ***
    ## sizesmall   -6.284e-01  1.314e-01  -4.783 6.54e-05 ***
    ## waterperc   -4.956e-02  7.653e-03  -6.476 8.79e-07 ***
    ## barkthick    1.752e-01  7.308e-02   2.397 0.024301 *  
    ## Ca           1.234e-04  3.203e-05   3.851 0.000725 ***
    ## Zn          -9.272e-03  2.698e-03  -3.436 0.002070 ** 
    ## N           -9.985e-01  4.258e-01  -2.345 0.027272 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3026 on 25 degrees of freedom
    ## Multiple R-squared:  0.735,  Adjusted R-squared:  0.6714 
    ## F-statistic: 11.56 on 6 and 25 DF,  p-value: 3.442e-06

### alpha--- don't interpret yet

########################################## 

Community as a predictor
------------------------

### First, filter community matrix to include only taxa that are present in a least 20% of all the samples. This step removes taxa that may not contribute much to our understanding of the relationship between species’ multivariate abundance and environment.

    ## [1] "Keep 150 of 6128 OTUs"

*Hyp:* Stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.

``` r
#fxns
reformatMatrix<-function(commmat){
  newmat<-matrix(as.numeric(as.matrix(commmat)), ncol=dim(commmat)[2], nrow=dim(commmat)[1])
  return(newmat)
}
CreateCommPMRpair<-function(timePoint, comm.mat, pmr.byStem.df.w){
    
  #make a dataframe using the current time point's pmr and remove NAs
  pmr.byStem.df.w %>%
    select_("codeStem",timePoint) %>%
    rename_("curr.time"=timePoint) %>%
    filter(!is.na(curr.time)) -> pmr.noNAs
  
  #subset the community matrix using these unique codeStems
  curr.comm<-comm.mat[row.names(comm.mat) %in% pmr.noNAs$codeStem,]
  curr.comm<-curr.comm[,colSums(curr.comm) != 0] #get rid of any OTU columns with all 0s
    
  #get rid of pmr rows for which there is missing community data
  pmr.noNAs %>%
    filter(codeStem %in% row.names(curr.comm)) -> curr.pmr
    
    #make sure the row orders match
    ord<-match(row.names(curr.comm), curr.pmr$codeStem)
    curr.pmr<-curr.pmr[ord,]
    
    #reformat the community matrix
    curr.comm.reform<-reformatMatrix(curr.comm)
    row.names(curr.comm.reform)<-row.names(curr.comm)
    colnames(curr.comm.reform)<-colnames(curr.comm)
    
    modelDat.list<-list(pmr=curr.pmr, comm=curr.comm.reform)
    
    return(modelDat.list)
    
  }

timeList<-list('time7','time13','time25','time37')
datasets<-lapply(timeList, CreateCommPMRpair, comm.mat = comm.otu.trimmed, pmr.byStem.df.w)
names(datasets)<-c('time7','time13','time25','time37')
datasets.nt<-lapply(timeList, CreateCommPMRpair, comm.mat = comm.otu, pmr.byStem.df.w)
names(datasets.nt)<-c('time7','time13','time25','time37')

#fit models
fit.t7 <- WAPLS(datasets[['time7']][['comm']], datasets[['time7']][['pmr']]$curr.time)
fit.t13 <- WAPLS(datasets[['time13']][['comm']], datasets[['time13']][['pmr']]$curr.time)
fit.t25 <- WAPLS(datasets[['time25']][['comm']], datasets[['time25']][['pmr']]$curr.time)
fit.t37 <- WAPLS(datasets[['time37']][['comm']], datasets[['time37']][['pmr']]$curr.time)
fit.t7.notrim <- WAPLS(datasets.nt[['time7']][['comm']], datasets.nt[['time7']][['pmr']]$curr.time)
fit.t13.notrim <- WAPLS(datasets.nt[['time13']][['comm']], datasets.nt[['time13']][['pmr']]$curr.time)
fit.t25.notrim <- WAPLS(datasets.nt[['time25']][['comm']], datasets.nt[['time25']][['pmr']]$curr.time)
fit.t37.notrim <- WAPLS(datasets.nt[['time37']][['comm']], datasets.nt[['time37']][['pmr']]$curr.time)

#cross-validate models using the leave-one-out method
fit.t7.cv <- crossval(fit.t7, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |==========                                                       |  15%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |============                                                     |  18%
      |                                                                       
      |=============                                                    |  20%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  23%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |=================                                                |  26%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |===========================                                      |  42%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  45%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  55%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  58%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |================================================                 |  74%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |==================================================               |  77%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |====================================================             |  80%
      |                                                                       
      |=====================================================            |  82%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=======================================================          |  85%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t13.cv <- crossval(fit.t13, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |========                                                         |  13%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |=============                                                    |  21%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  24%
      |                                                                       
      |=================                                                |  25%
      |                                                                       
      |==================                                               |  27%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |====================                                             |  30%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  70%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  73%
      |                                                                       
      |================================================                 |  75%
      |                                                                       
      |==================================================               |  76%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |====================================================             |  79%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |=========================================================        |  87%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t25.cv <- crossval(fit.t25, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |========                                                         |  13%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |=============                                                    |  21%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  24%
      |                                                                       
      |=================                                                |  25%
      |                                                                       
      |==================                                               |  27%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |====================                                             |  30%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  70%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  73%
      |                                                                       
      |================================================                 |  75%
      |                                                                       
      |==================================================               |  76%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |====================================================             |  79%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |=========================================================        |  87%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t37.cv <- crossval(fit.t37, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   7%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |=========                                                        |  13%
      |                                                                       
      |==========                                                       |  15%
      |                                                                       
      |===========                                                      |  16%
      |                                                                       
      |============                                                     |  18%
      |                                                                       
      |=============                                                    |  20%
      |                                                                       
      |==============                                                   |  21%
      |                                                                       
      |===============                                                  |  23%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |=================                                                |  26%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |===================                                              |  30%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |=====================                                            |  33%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |=======================                                          |  36%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  39%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  61%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |==========================================                       |  64%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |============================================                     |  67%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |==============================================                   |  70%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |================================================                 |  74%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |==================================================               |  77%
      |                                                                       
      |===================================================              |  79%
      |                                                                       
      |====================================================             |  80%
      |                                                                       
      |=====================================================            |  82%
      |                                                                       
      |======================================================           |  84%
      |                                                                       
      |=======================================================          |  85%
      |                                                                       
      |========================================================         |  87%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  93%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t7.nt.cv <- crossval(fit.t7.notrim, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |==========                                                       |  15%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |============                                                     |  18%
      |                                                                       
      |=============                                                    |  20%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  23%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |=================                                                |  26%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |===========================                                      |  42%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  45%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  55%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  58%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |================================================                 |  74%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |==================================================               |  77%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |====================================================             |  80%
      |                                                                       
      |=====================================================            |  82%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=======================================================          |  85%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t13.nt.cv <- crossval(fit.t13.notrim, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |========                                                         |  13%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |=============                                                    |  21%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  24%
      |                                                                       
      |=================                                                |  25%
      |                                                                       
      |==================                                               |  27%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |====================                                             |  30%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  70%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  73%
      |                                                                       
      |================================================                 |  75%
      |                                                                       
      |==================================================               |  76%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |====================================================             |  79%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |=========================================================        |  87%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t25.nt.cv <- crossval(fit.t25.notrim, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |========                                                         |  13%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |=============                                                    |  21%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  24%
      |                                                                       
      |=================                                                |  25%
      |                                                                       
      |==================                                               |  27%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |====================                                             |  30%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  70%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  73%
      |                                                                       
      |================================================                 |  75%
      |                                                                       
      |==================================================               |  76%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |====================================================             |  79%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |=========================================================        |  87%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
fit.t37.nt.cv <- crossval(fit.t37.notrim, cv.method="loo")
```

    ## Cross-validating:
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   7%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |=========                                                        |  13%
      |                                                                       
      |==========                                                       |  15%
      |                                                                       
      |===========                                                      |  16%
      |                                                                       
      |============                                                     |  18%
      |                                                                       
      |=============                                                    |  20%
      |                                                                       
      |==============                                                   |  21%
      |                                                                       
      |===============                                                  |  23%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |=================                                                |  26%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |===================                                              |  30%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |=====================                                            |  33%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |=======================                                          |  36%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |==========================                                       |  39%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  61%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |==========================================                       |  64%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |============================================                     |  67%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |==============================================                   |  70%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |================================================                 |  74%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |==================================================               |  77%
      |                                                                       
      |===================================================              |  79%
      |                                                                       
      |====================================================             |  80%
      |                                                                       
      |=====================================================            |  82%
      |                                                                       
      |======================================================           |  84%
      |                                                                       
      |=======================================================          |  85%
      |                                                                       
      |========================================================         |  87%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |=============================================================    |  93%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |=================================================================| 100%

``` r
#perform randomization t-test to test the significance of a cross-validated model....
```

### time 7

**none of the community components are significant predictors**

    ##             RMSE         R2      Avg.Bias  Max.Bias       Skill delta.RMSE
    ## Comp01 0.1096009 0.04039245 -0.0042025315 0.2086640   -7.426179   3.646601
    ## Comp02 0.1254558 0.04123537  0.0045006905 0.1767666  -40.754864  14.466035
    ## Comp03 0.1472351 0.02785978  0.0004267154 0.1727394  -93.867398  17.360158
    ## Comp04 0.1501175 0.04487863  0.0040991337 0.1600327 -101.532243   1.957664
    ## Comp05 0.1563567 0.04317337  0.0024552841 0.1777311 -118.632632   4.156223
    ##            p
    ## Comp01 0.715
    ## Comp02 0.957
    ## Comp03 0.977
    ## Comp04 0.729
    ## Comp05 0.930

    ##             RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1160401 0.03500705 0.004146464 0.2118409 -20.41985  9.7359781
    ## Comp02 0.1140335 0.03867892 0.005554219 0.2053950 -16.29114 -1.7292500
    ## Comp03 0.1141955 0.04829218 0.005694794 0.2145943 -16.62175  0.1420493
    ## Comp04 0.1137256 0.05486492 0.004972683 0.2101628 -15.66399 -0.4114765
    ## Comp05 0.1124598 0.06682365 0.006100160 0.2058893 -13.10370 -1.1129729
    ##            p
    ## Comp01 0.919
    ## Comp02 0.316
    ## Comp03 0.518
    ## Comp04 0.361
    ## Comp05 0.112

### time 13

**none of the community components are significant predictors**

    ##             RMSE         R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.2139657 0.21846772 -0.01496427 0.6260449  -74.86769 32.2375476
    ## Comp02 0.2149868 0.06885441 -0.01490427 0.5414328  -76.54081  0.4772558
    ## Comp03 0.2428530 0.06690494 -0.01653354 0.6290506 -125.27248 12.9617769
    ## Comp04 0.2693172 0.08533308 -0.01524579 0.6038045 -177.04454 10.8972377
    ## Comp05 0.2994919 0.07108880 -0.01525324 0.6220238 -242.60334 11.2041512
    ##            p
    ## Comp01 1.000
    ## Comp02 0.565
    ## Comp03 0.999
    ## Comp04 0.999
    ## Comp05 0.997

    ##             RMSE           R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.2163312 4.753407e-02 -0.02301907 0.6295907 -78.75562 33.6995201
    ## Comp02 0.2023872 5.526947e-06 -0.03447532 0.5853071 -56.45433 -6.4456583
    ## Comp03 0.2200750 8.812263e-03 -0.04106598 0.6142264 -84.99614  8.7395515
    ## Comp04 0.2244837 9.323948e-03 -0.04242462 0.6066593 -92.48239  2.0032860
    ## Comp05 0.2264918 9.961345e-03 -0.04426489 0.6070759 -95.94148  0.8945483
    ##            p
    ## Comp01 0.998
    ## Comp02 0.083
    ## Comp03 0.999
    ## Comp04 0.880
    ## Comp05 0.855

### time 25

**none of the community components are significant predictors**

    ##             RMSE         R2   Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1589171 0.07590492 0.01732259 0.3131693  -9.443048   4.615031
    ## Comp02 0.1723926 0.05607701 0.02297917 0.3430466 -28.790674   8.479611
    ## Comp03 0.1786746 0.05969270 0.02226170 0.3307871 -38.347911   3.643983
    ## Comp04 0.1813230 0.06504168 0.02194579 0.3286657 -42.479674   1.482266
    ## Comp05 0.1847188 0.06871557 0.02042740 0.3579705 -47.866326   1.872786
    ##            p
    ## Comp01 0.755
    ## Comp02 0.924
    ## Comp03 0.796
    ## Comp04 0.680
    ## Comp05 0.768

    ##             RMSE         R2      Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1568901 0.08501311 -1.251646e-03 0.2783429  -6.668957  3.2806647
    ## Comp02 0.1734482 0.03761080  1.856222e-03 0.3039221 -30.372702 10.5539613
    ## Comp03 0.1766744 0.02875939 -6.219912e-05 0.2900092 -35.267824  1.8600581
    ## Comp04 0.1753588 0.03066941  2.508066e-04 0.2874978 -33.260727 -0.7446700
    ## Comp05 0.1747009 0.03192784  1.039806e-03 0.2883656 -32.262694 -0.3751703
    ##            p
    ## Comp01 0.706
    ## Comp02 0.997
    ## Comp03 0.888
    ## Comp04 0.213
    ## Comp05 0.265

### time 37

**Comp01 is significant** **if trim out OTUs that are not present in at least 20% of samples...Comp02 is marginally significant**

    ##             RMSE        R2     Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1884529 0.1515172  0.033192872 0.2942286  5.557560 -2.8184996
    ## Comp02 0.1749664 0.2407742  0.008951682 0.2788971 18.591259 -7.1564113
    ## Comp03 0.1752001 0.2660234 -0.000556893 0.2185399 18.373655  0.1335602
    ## Comp04 0.1879489 0.2205367  0.003539710 0.2822099  6.062045  7.2766925
    ## Comp05 0.1967444 0.1897768  0.002389760 0.2656368 -2.935780  4.6797384
    ##            p
    ## Comp01 0.393
    ## Comp02 0.043
    ## Comp03 0.497
    ## Comp04 0.986
    ## Comp05 0.936

    ##             RMSE        R2     Avg.Bias  Max.Bias    Skill  delta.RMSE
    ## Comp01 0.1668533 0.2925742 -0.001261747 0.3122426 25.96599 -13.9569840
    ## Comp02 0.1615472 0.3377994 -0.007472631 0.2888952 30.59988  -3.1801339
    ## Comp03 0.1682930 0.3025985 -0.012605013 0.2854822 24.68285   4.1757949
    ## Comp04 0.1708756 0.2921546 -0.014679983 0.2863511 22.35357   1.5345464
    ## Comp05 0.1725708 0.2854071 -0.014202926 0.2902932 20.80528   0.9920909
    ##            p
    ## Comp01 0.037
    ## Comp02 0.186
    ## Comp03 0.946
    ## Comp04 0.859
    ## Comp05 0.908

Investigate the biology underlying time37-associated coefs for Comp02

By trophic mode (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By guild (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

*Hyp:* Average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### r2

**none of the community components are significant predictors**

    ##             RMSE           R2    Avg.Bias  Max.Bias     Skill  delta.RMSE
    ## Comp01 0.1013107 1.494708e-05 0.002024841 0.2609272 -40.32548 18.45905577
    ## Comp02 0.1104027 4.293834e-03 0.003308000 0.2568255 -66.64251  8.97443525
    ## Comp03 0.1114560 1.027904e-02 0.003371373 0.2497766 -69.83735  0.95404167
    ## Comp04 0.1119467 1.116217e-02 0.002666967 0.2514263 -71.33590  0.44020287
    ## Comp05 0.1120426 1.133987e-02 0.002626418 0.2509619 -71.62965  0.08568655
    ##            p
    ## Comp01 0.894
    ## Comp02 1.000
    ## Comp03 0.702
    ## Comp04 0.808
    ## Comp05 0.720

### k

**none of the community components are significant predictors**

    ##             RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.1036580 0.01276672 0.000282797 0.2658329 -22.99604 10.9035799
    ## Comp02 0.1044388 0.01904155 0.011480363 0.2688760 -24.85590  0.7532281
    ## Comp03 0.1059335 0.01805338 0.010666194 0.2689647 -28.45534  1.4311977
    ## Comp04 0.1087908 0.02741427 0.013442320 0.2667951 -35.47840  2.6972785
    ## Comp05 0.1125616 0.01442363 0.011293413 0.2704606 -45.03280  3.4661026
    ##            p
    ## Comp01 0.895
    ## Comp02 0.522
    ## Comp03 0.606
    ## Comp04 0.877
    ## Comp05 0.953

### t70

**none of the community components are significant predictors**

    ##             RMSE         R2      Avg.Bias Max.Bias      Skill  delta.RMSE
    ## Comp01 0.5638322 0.02747152  0.0291722604 1.503069 -19.039193  9.10508360
    ## Comp02 0.5363187 0.06687807 -0.0040312802 1.258095  -7.705054 -4.87973434
    ## Comp03 0.5367893 0.08213251 -0.0004053627 1.153682  -7.894179  0.08775921
    ## Comp04 0.5583885 0.08221500 -0.0056512738 1.349964 -16.751702  4.02377323
    ## Comp05 0.5693485 0.06691379  0.0058129906 1.431625 -21.379860  1.96278895
    ##            p
    ## Comp01 0.853
    ## Comp02 0.227
    ## Comp03 0.512
    ## Comp04 0.840
    ## Comp05 0.833

### alpha --- don't interpret yet

########################################## 

Community+traits as a predictor
-------------------------------

*Hyp:* After accounting for variation in decay due to wood traits (no models with barkthick or density), stem-specific initial microbial communitiy compositions will predict variation in percent mass loss, particularly in the early stages of decay.

### time 7

**none of the community components are significant predictors**

    ##             RMSE          R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1115850 0.007432694 0.004235539 0.2120931  -52.60022  23.531460
    ## Comp02 0.1209657 0.004020937 0.013461663 0.1991320  -79.33627   8.406790
    ## Comp03 0.1291656 0.003573218 0.013284545 0.2371421 -104.47335   6.778618
    ## Comp04 0.1252574 0.001041763 0.011006738 0.2231303  -92.28723  -3.025653
    ## Comp05 0.1199016 0.001036001 0.011918383 0.2281983  -76.19499  -4.275843
    ##            p
    ## Comp01 0.993
    ## Comp02 0.769
    ## Comp03 0.878
    ## Comp04 0.247
    ## Comp05 0.328

### time 13

**Comp02 is significant** **if trim out OTUs that are not present in at least 20% of samples...no components are significant**

    ##             RMSE           R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1684942 7.969996e-03 -0.01227286 0.3008319  -59.80987  26.415927
    ## Comp02 0.1708079 8.029389e-05 -0.01134351 0.3428988  -64.22886   1.373151
    ## Comp03 0.1813353 1.626797e-02 -0.01398033 0.3387630  -85.09640   6.163259
    ## Comp04 0.1861811 2.670958e-03 -0.01869200 0.3402859  -95.12131   2.672318
    ## Comp05 0.2000881 9.511791e-03 -0.01941544 0.3603079 -125.35962   7.469616
    ##            p
    ## Comp01 0.998
    ## Comp02 0.592
    ## Comp03 0.986
    ## Comp04 0.820
    ## Comp05 0.999

    ##             RMSE         R2    Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1958078 0.09587553 -0.01406066 0.3977193 -115.82077 46.9083971
    ## Comp02 0.1807162 0.02022406 -0.02532412 0.3640608  -83.83481 -7.7073210
    ## Comp03 0.1915530 0.02089810 -0.02731081 0.3899540 -106.54333  5.9965450
    ## Comp04 0.1942539 0.02389339 -0.02709579 0.3915635 -112.40896  1.4100120
    ## Comp05 0.1938709 0.02017169 -0.02696350 0.3920242 -111.57224 -0.1971538
    ##            p
    ## Comp01 1.000
    ## Comp02 0.052
    ## Comp03 0.968
    ## Comp04 0.860
    ## Comp05 0.417

Investigate the biology underlying time13-associated coefs for Comp02

By trophic mode ![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-32-1.png)![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-32-2.png)

By phylum

### time 25

**none of the community components are significant predictors**

    ##             RMSE        R2   Avg.Bias  Max.Bias     Skill delta.RMSE     p
    ## Comp01 0.1620698 0.1329636 0.01919070 0.3483804 -111.3412  45.375801 1.000
    ## Comp02 0.1795585 0.1117887 0.02724735 0.3954964 -159.4131  10.790831 0.939
    ## Comp03 0.1926807 0.1869273 0.02416387 0.4476233 -198.7147   7.308055 0.999
    ## Comp04 0.1905622 0.1539854 0.01920298 0.4451090 -192.1820  -1.099513 0.207
    ## Comp05 0.1957642 0.1553971 0.02117982 0.4710456 -208.3519   2.729837 0.952

### time 37

**none of the community components are significant predictors**

    ##             RMSE          R2     Avg.Bias  Max.Bias      Skill delta.RMSE
    ## Comp01 0.1519790 0.004659334  0.002664229 0.3654074  -56.17994  24.971974
    ## Comp02 0.1583613 0.009032758 -0.001947198 0.3948058  -69.57286   4.199477
    ## Comp03 0.1644513 0.007509707 -0.001848282 0.3944939  -82.86600   3.845652
    ## Comp04 0.1741376 0.004828966 -0.008279261 0.3849332 -105.04227   5.890064
    ## Comp05 0.1868349 0.013372894 -0.008065750 0.3920523 -136.03383   7.291528
    ##            p
    ## Comp01 0.994
    ## Comp02 0.821
    ## Comp03 0.914
    ## Comp04 0.988
    ## Comp05 0.998

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial communitiy compositions will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### r2

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was marginally significant

    ##              RMSE           R2     Avg.Bias  Max.Bias      Skill
    ## Comp01 0.07083530 0.0004108819 -0.006446032 0.1400555  -61.51372
    ## Comp02 0.09185757 0.0006934417  0.002667991 0.1556453 -171.60633
    ## Comp03 0.09573347 0.0023892283  0.004948795 0.1592020 -195.01055
    ## Comp04 0.09706111 0.0019152166  0.003816465 0.1568796 -203.24973
    ## Comp05 0.09665378 0.0024386184  0.004673572 0.1604111 -200.70985
    ##        delta.RMSE     p
    ## Comp01  27.088047 0.993
    ## Comp02  29.677682 1.000
    ## Comp03   4.219463 0.925
    ## Comp04   1.386804 0.791
    ## Comp05  -0.419656 0.303

### k

**none of the community components are significant predictors**

    ##              RMSE         R2    Avg.Bias  Max.Bias     Skill delta.RMSE
    ## Comp01 0.05403597 0.01169109 0.000422351 0.1254135 -31.24173  14.560784
    ## Comp02 0.05065152 0.04533652 0.002916937 0.1244362 -15.31638  -6.263328
    ## Comp03 0.04993864 0.05990280 0.004727269 0.1237288 -12.09324  -1.407427
    ## Comp04 0.05147371 0.04410936 0.005052921 0.1270141 -19.09043   3.073904
    ## Comp05 0.05142763 0.04472957 0.005796166 0.1291158 -18.87733  -0.089509
    ##            p
    ## Comp01 0.879
    ## Comp02 0.180
    ## Comp03 0.194
    ## Comp04 0.935
    ## Comp05 0.490

### t70

**none of the community components are significant predictors**

Note: This result changed when I changed waterperc to g water/g wet weight. When waterperc was in terms of g/g dry weight, Comp05 was a significant predictor

    ##             RMSE        R2   Avg.Bias  Max.Bias     Skill delta.RMSE     p
    ## Comp01 0.2564979 0.1696737 0.02216625 0.4864828  8.061448 -4.1154070 0.291
    ## Comp02 0.2615463 0.1965626 0.01788789 0.4573872  4.406754  1.9682055 0.560
    ## Comp03 0.2745169 0.1851611 0.01979289 0.4458862 -5.309638  4.9591861 0.893
    ## Comp04 0.2776540 0.1774745 0.02286057 0.4449253 -7.730289  1.1427719 0.735
    ## Comp05 0.2795150 0.1716365 0.02355805 0.4637620 -9.179284  0.6702646 0.637

### alpha --- don't interpret yet

**Note: no longer warrented based on analyses with waterperc represented as g water / g wet weight...**

Investigate the biology underlying t70-associated coefs for Comp05

By trophic mode (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By guild (note that apparently empty cateogies have at least 1 data point, the violin plot just doesn't show it)

By phylum

########################################## 

Diversity (and diversity of specific clades) as a predictor
-----------------------------------------------------------

**Note that the full community matrix was used for these analyses**

*Hyp:* Greater microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because of the selection effect for fast decayers and complementarity among taxa for decay.
Hyp-Alt: Greater microbial diversity will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because taxa will be allocating more of their resources to combat one another.

### Richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8943
    ## mean       1 0.000245 0.0002451  0.0305 0.8625
    ## Residuals 30 0.240983 0.0080328

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0055 0.03285 *
    ## mean       1 0.000163 0.000163  0.0198 0.88904  
    ## Residuals 30 0.246927 0.008231                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9953 0.02041 *
    ## mean       1 0.0226 0.02258  0.0925 0.76317  
    ## Residuals 30 7.3263 0.24421                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0205 0.8871
    ## mean       1 0.15961 0.159615  2.3308 0.1373
    ## Residuals 30 2.05446 0.068482

### Shannon's H

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8939
    ## mean       1 0.002429 0.0024291  0.3052 0.5848
    ## Residuals 30 0.238799 0.0079600

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1360 0.03081 *
    ## mean       1 0.006437 0.006437  0.8024 0.37750  
    ## Residuals 30 0.240653 0.008022                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2830 0.01784 *
    ## mean       1 0.3580 0.35804  1.5365 0.22475  
    ## Residuals 30 6.9908 0.23303                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.0014 0.001404  0.0196 0.8896
    ## mean       1 0.0657 0.065698  0.9174 0.3458
    ## Residuals 30 2.1484 0.071613

*Hyp:* Greater saprotroph and basidiomycete richness will lead to better-fitting decay models (ne.r2), faster decay (k), and less lagginess (alpha) because the community does not need to wait for the arrival of key decayers to act on the wood substrate.
Hyp-Alt: Greater saprotroph and basidiomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because decayers will be allocating more of their resources to combat one another.

### Saprotroph richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000006 0.0000062  0.0008 0.9780
    ## Residuals 30 0.241222 0.0080407

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.04120 0.041200  5.0515 0.03211 *
    ## mean       1 0.00241 0.002410  0.2955 0.59075  
    ## Residuals 30 0.24468 0.008156                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9893 0.02046 *
    ## mean       1 0.0152 0.01524  0.0624 0.80452  
    ## Residuals 30 7.3336 0.24445                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)  
    ## size       1 0.00140 0.001404  0.0211 0.8854  
    ## mean       1 0.22078 0.220782  3.3229 0.0783 .
    ## Residuals 30 1.99330 0.066443                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Basidio richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000072 0.0000717  0.0089 0.9254
    ## Residuals 30 0.241157 0.0080386

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1545 0.03053 *
    ## mean       1 0.007301 0.007301  0.9134 0.34686  
    ## Residuals 30 0.239789 0.007993                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.1675 0.01883 *
    ## mean       1 0.2271 0.22708  0.9565 0.33588  
    ## Residuals 30 7.1218 0.23739                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0200 0.8886
    ## mean       1 0.10539 0.105391  1.4994 0.2303
    ## Residuals 30 2.10869 0.070290

*Hyp:* Greater pathogen and oomycete richness will lead to worse-fitting decay models (ne.r2), slower decay (k), and more lagginess (alpha) because the presence of these organisms will inhibit the establishment and activity of decayers.

### Pathogen richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8940
    ## mean       1 0.001760 0.0017599  0.2205 0.6421
    ## Residuals 30 0.239469 0.0079823

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0603 0.03197 *
    ## mean       1 0.002835 0.002835  0.3482 0.55956  
    ## Residuals 30 0.244255 0.008142                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2855 0.01782 *
    ## mean       1 0.3608 0.36079  1.5489 0.22294  
    ## Residuals 30 6.9881 0.23294                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.00140 0.001404  0.0215 0.88430  
    ## mean       1 0.25844 0.258439  3.9645 0.05564 .
    ## Residuals 30 1.95564 0.065188                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Oomycete richness

**No pattern**

    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0184 0.8931
    ## mean       1 0.005902 0.0059025  0.7525 0.3926
    ## Residuals 30 0.235326 0.0078442

    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0038 0.03288 *
    ## mean       1 0.000077 0.000077  0.0093 0.92374  
    ## Residuals 30 0.247013 0.008234                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9831 0.02052 *
    ## mean       1 0.0076 0.00763  0.0312 0.86099  
    ## Residuals 30 7.3412 0.24471                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0192 0.8906
    ## mean       1 0.02338 0.023381  0.3202 0.5757
    ## Residuals 30 2.19070 0.073023

############################################## 

Diversity plus traits as a predictor
------------------------------------

*Hyp:* After accounting for variation in decay due to wood traits, average initial microbial diversity (richness, Shannon diversity, ... add phylogenetic diversity) will predict variation in decay model fit (r2), rate (k), and lagginess (alpha).

### Richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8943
    ## mean       1 0.000245 0.0002451  0.0305 0.8625
    ## Residuals 30 0.240983 0.0080328               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0055 0.03285 *
    ## mean       1 0.000163 0.000163  0.0198 0.88904  
    ## Residuals 30 0.246927 0.008231                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9953 0.02041 *
    ## mean       1 0.0226 0.02258  0.0925 0.76317  
    ## Residuals 30 7.3263 0.24421                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0205 0.8871
    ## mean       1 0.15961 0.159615  2.3308 0.1373
    ## Residuals 30 2.05446 0.068482

### Shannon's H

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8939
    ## mean       1 0.002429 0.0024291  0.3052 0.5848
    ## Residuals 30 0.238799 0.0079600               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1360 0.03081 *
    ## mean       1 0.006437 0.006437  0.8024 0.37750  
    ## Residuals 30 0.240653 0.008022                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2830 0.01784 *
    ## mean       1 0.3580 0.35804  1.5365 0.22475  
    ## Residuals 30 6.9908 0.23303                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.0014 0.001404  0.0196 0.8896
    ## mean       1 0.0657 0.065698  0.9174 0.3458
    ## Residuals 30 2.1484 0.071613

### Saprotroph richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000006 0.0000062  0.0008 0.9780
    ## Residuals 30 0.241222 0.0080407               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.04120 0.041200  5.0515 0.03211 *
    ## mean       1 0.00241 0.002410  0.2955 0.59075  
    ## Residuals 30 0.24468 0.008156                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9893 0.02046 *
    ## mean       1 0.0152 0.01524  0.0624 0.80452  
    ## Residuals 30 7.3336 0.24445                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)  
    ## size       1 0.00140 0.001404  0.0211 0.8854  
    ## mean       1 0.22078 0.220782  3.3229 0.0783 .
    ## Residuals 30 1.99330 0.066443                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Basidio richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0179 0.8944
    ## mean       1 0.000072 0.0000717  0.0089 0.9254
    ## Residuals 30 0.241157 0.0080386               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.1545 0.03053 *
    ## mean       1 0.007301 0.007301  0.9134 0.34686  
    ## Residuals 30 0.239789 0.007993                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.1675 0.01883 *
    ## mean       1 0.2271 0.22708  0.9565 0.33588  
    ## Residuals 30 7.1218 0.23739                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0200 0.8886
    ## mean       1 0.10539 0.105391  1.4994 0.2303
    ## Residuals 30 2.10869 0.070290

### Pathogen richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0181 0.8940
    ## mean       1 0.001760 0.0017599  0.2205 0.6421
    ## Residuals 30 0.239469 0.0079823               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0603 0.03197 *
    ## mean       1 0.002835 0.002835  0.3482 0.55956  
    ## Residuals 30 0.244255 0.008142                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  6.2855 0.01782 *
    ## mean       1 0.3608 0.36079  1.5489 0.22294  
    ## Residuals 30 6.9881 0.23294                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.00140 0.001404  0.0215 0.88430  
    ## mean       1 0.25844 0.258439  3.9645 0.05564 .
    ## Residuals 30 1.95564 0.065188                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Oomycete richness

**No pattern**

    ## $r2
    ## Analysis of Variance Table
    ## 
    ## Response: ne.r2
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## size       1 0.000144 0.0001441  0.0184 0.8931
    ## mean       1 0.005902 0.0059025  0.7525 0.3926
    ## Residuals 30 0.235326 0.0078442               
    ## 
    ## $k
    ## Analysis of Variance Table
    ## 
    ## Response: k
    ##           Df   Sum Sq  Mean Sq F value  Pr(>F)  
    ## size       1 0.041200 0.041200  5.0038 0.03288 *
    ## mean       1 0.000077 0.000077  0.0093 0.92374  
    ## Residuals 30 0.247013 0.008234                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $t70
    ## Analysis of Variance Table
    ## 
    ## Response: t70
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## size       1 1.4641 1.46411  5.9831 0.02052 *
    ## mean       1 0.0076 0.00763  0.0312 0.86099  
    ## Residuals 30 7.3412 0.24471                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $alpha
    ## Analysis of Variance Table
    ## 
    ## Response: alpha
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## size       1 0.00140 0.001404  0.0192 0.8906
    ## mean       1 0.02338 0.023381  0.3202 0.5757
    ## Residuals 30 2.19070 0.073023

############################################## 

Relationship between wood traits and community
----------------------------------------------

*Hyp:* Average initial microbial communitiy compositions will covary with initial wood traits

``` r
#average OTU abundances by code
meanOTUabund.trim2<-AverageOTUabund_byCode(comm.otu.trimmed, seqSamples)

#get rid of missing data from traits.mean
traits.mean %>%
  filter(!is.na(waterperc)) %>% #get NaN row
  filter(code %in% row.names(meanOTUabund.trim2)) -> traits.mean.trim #get rid of rows for which there is missing community data... this is eusc

#make sure dimensions of the community matrix match
meanOTUabund.trim3<-meanOTUabund.trim2[row.names(meanOTUabund.trim2) %in% traits.mean.trim$code,]

#make a envVars matrix for model fitting
envVars<-as.matrix(traits.mean.trim[,-(1:3)])

#fit models
fit.cca <- cca(meanOTUabund.trim3 ~ envVars)
fit.cca
```

    ## Call: cca(formula = meanOTUabund.trim3 ~ envVars)
    ## 
    ##               Inertia Proportion Rank
    ## Total          6.3027     1.0000     
    ## Constrained    2.4529     0.3892   11
    ## Unconstrained  3.8499     0.6108   20
    ## Inertia is mean squared contingency coefficient 
    ## 
    ## Eigenvalues for constrained axes:
    ##   CCA1   CCA2   CCA3   CCA4   CCA5   CCA6   CCA7   CCA8   CCA9  CCA10 
    ## 0.4651 0.3617 0.3376 0.3269 0.2671 0.1797 0.1409 0.1104 0.1072 0.0829 
    ##  CCA11 
    ## 0.0733 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
    ## 0.5121 0.4635 0.3576 0.2992 0.2897 0.2709 0.2157 0.2003 
    ## (Showed only 8 of all 20 unconstrained eigenvalues)

``` r
plot(fit.cca)
```

![](readme_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-62-1.png)

``` r
anova(fit.cca)
```

    ## Permutation test for cca under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = meanOTUabund.trim3 ~ envVars)
    ##          Df ChiSquare      F Pr(>F)
    ## Model    11    2.4529 1.1584  0.139
    ## Residual 20    3.8499

############################################## 

Extra pieces
------------

1.  *code/testing\_time\_zero.Rmd* -- Including t=0 points to fit decay model affects the liklihood and the model selection criteria, but the curve fits are identical with this formulation. Excluding the t=0 fits has an effect of prefering simpler models, which is the same effect as increasing the penalty for model complexity.

2.  *code/initialDist\_vs\_decayDist\_btwCode.Rmd* -- No apparent relationship between species+size dissimilarities in initial microbial community composition (bray and jaccard) and decay trajectory params

3.  *code/boralOTUpairs\_vs\_decay.Rmd* -- No apparent relationship between frequency of boral-ID'd positively/negatively correlated OTU pairs and decay params

4.  *code/withinInitialDist\_vs\_decayR2.Rmd* -- No apparent relationship between initial microbial diversity WITHIN species+size and decay model R2

5.  *code/unexpectedTaxa.Rmd* -- Mycorrhizal fungi and animal-associated fungi that somehow made it into our OTU table
