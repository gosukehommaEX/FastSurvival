
# FastSurvival

<!-- badges: start -->
<!-- badges: end -->

The goal of FastSurvival is to ...

## Installation

You can install the development version of FastSurvival from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/FastSurvival")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(survival)
library(FastSurvival)
library(microbenchmark)

microbenchmark(
  surv.diff.result  = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
  FastLRtest.result = FastLRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2)
)

# Unit: microseconds
#               expr     min      lq      mean   median      uq    max neval
#   surv.diff.result 586.201 620.401 719.57500 647.2005 705.001 1935.6   100
#  FastLRtest.result  18.001  20.951  28.76396  28.1510  31.101  195.0   100
```

