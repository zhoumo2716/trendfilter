
<!-- README.md is generated from README.Rmd. Please edit that file -->

# trendfilter

<!-- badges: start -->

[![R-CMD-check](https://github.com/glmgen/trendfilter/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glmgen/trendfilter/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/trendfilter)](https://CRAN.R-project.org/package=trendfilter)
<!-- badges: end -->

The goal of trendfilter is to **solve** nonparametric regression.

## Installation

You can install the development version of trendfilter from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("glmgen/trendfilter")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(trendfilter)
library(ggplot2)
x <- 1:100 / 101 * 2 * pi
y <- sin(x) + .2 * rnorm(100)
out <- trendfilter(y, x, nlambda = 15)
plot(out) +
  geom_point(data = data.frame(x = x, y = y), aes(x, y), color = "black")
```

<img src="man/figures/README-example-1.svg" width="100%" />
