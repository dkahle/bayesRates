<!-- README.md is generated from README.Rmd. Please edit that file -->
bayesRates
==========

**bayesRates** is an R package that allows users to (1) perform two-sample tests with binomial and Poisson data from a Bayesian perspective and (2) determine sample sizes for designing such procedures.

``` r
library(bayesRates)
#> Loading required package: ggplot2
```

Binomial tests
==============

Suppose that we have two coins, a quarter and a half dollar, and we're interested in determining whether or not they flip heads with the same probability. After flipping each 50 times, the quarter flips heads 26 times and the half dollar flips heads 34 times. Should we conclude the likelihoods are the same?

One way to address the problem is with a standard test of proportions using `prop.test()` (`binom.test()` is similar, but uses the exact different procedure):

``` r
prop.test(x = c(26, 34), n = c(50, 50), correct = FALSE)
#> 
#>  2-sample test for equality of proportions without continuity
#>  correction
#> 
#> data:  c(26, 34) out of c(50, 50)
#> X-squared = 2.6667, df = 1, p-value = 0.1025
#> alternative hypothesis: two.sided
#> 95 percent confidence interval:
#>  -0.34945868  0.02945868
#> sample estimates:
#> prop 1 prop 2 
#>   0.52   0.68
```

This function uses the standard two-sample proportions test using pooling:

``` r
x <- 26; nx <- 50; px <- x/nx
y <- 34; ny <- 50; py <- y/ny
pp <- (x+y)/(nx+ny)

t <- (px-py)/sqrt(pp*(1-pp)/nx + pp*(1-pp)/ny)
2*pnorm(t) # two-sided p-value using the clt
#> [1] 0.1024704
```

Installation
------------

-   From Github: `devtools::install_github("dkahle/bayesRates")`

<!-- * From CRAN: `install.packages("bayesRates")` -->
