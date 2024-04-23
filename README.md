
<!-- README.md is generated from README.Rmd. Please edit that file -->

# netidmtpreg

<!-- badges: start -->

[![R-CMD-check](https://github.com/qmarcou/netidmtpreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/qmarcou/netidmtpreg/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of netidmtpreg is to …

## Installation

You can install the development version of netidmtpreg like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

If working from source through a git clone, the package should be
installed locally to get tests running `future::multisession` planning
working. This can be accomplished using `devtools`

``` r
devtools::dev_mode(on = TRUE)
devtools::install_local(force = TRUE) # force package update
testthat::test_package("netidmtpreg")
# or testthat::test_check("netidmtpreg") to run R CMD check
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(netidmtpreg)
#> Loading required package: doParallel
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
#> Loading required package: survival
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
