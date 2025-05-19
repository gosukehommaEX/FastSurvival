
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

survdiff.result  = survdiff(Surv(futime, fustat) ~ rx, data = ovarian)
FastLRtest.result = FastLRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2)

print(list(survdiff.chisq = survdiff.result[['chisq']], FastLRtest.chisq = FastLRtest.result ^ 2))
# $survdiff.chisq
# [1] 1.06274
# 
# $FastLRtest.chisq
# [1] 1.06274

microbenchmark(
  survdiff.result  = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
  FastLRtest.result = FastLRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2)
)

# Unit: microseconds
#               expr     min      lq      mean  median       uq      max neval
#    survdiff.result 598.301 668.101 827.18999 742.851 835.8005 1601.901   100
#  FastLRtest.result  56.002  72.750  96.40398  89.151 102.3005  259.701   100
```

