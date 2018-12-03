How do wood species and sizes vary in decay?
================
Marissa Lee
12/2/2018

``` r
#chunk options

knitr::opts_chunk$set(echo = FALSE, message=FALSE, warning=FALSE)
```

Load libraries, functions, data

Plot percent mass remaining (pmr) for each sample over time ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-2-1.png)

Calculate decay trajectory fits for each code (species+size)

For the weibull model, alpha = shape, beta = scale.

Compare the negative exp. vs weibull models by plotting time to 70% mass remaining (t70) for each species+size ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-4-1.png)

Another view... Figure S2. Comparison of negative exp and weibull models using t70 ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-5-1.png)

Figure S3. Comparison of negative exp and weibull models using AIC ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-6-1.png)

Figure 1. Wood species x decay params (weibull model) ![](decayPatterns_files/figure-markdown_github/unnamed-chunk-7-1.png)

Figure 2. Time x percent mass remaining

    ## [1]  0  7 13 25 37 59

![](decayPatterns_files/figure-markdown_github/unnamed-chunk-8-1.png)
