
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
#               expr     min      lq      mean   median       uq      max neval
#   surv.diff.result 609.001 752.001 811.44694 782.1005 858.6510 1287.801   100
#  FastLRtest.result  60.400  79.851  97.57995  93.0515 107.5515  214.700   100
```

