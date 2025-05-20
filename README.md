

# FastSurvival


<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/FastSurvival)](http://cran.r-project.org/package=FastSurvival)
<!-- badges: end -->

The goal of FastSurvival is to (1) implement log-rank test, (2) estimate hazard ratio, and (3) generate simulation datasets, much faster.


## Installation

You can install the development version of FastSurvival from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/FastSurvival")
```


## Usage

``` r
library(survival)
library(FastSurvival)
library(microbenchmark)

## Log-rank test
# Compare the log-rank test result by FastSurvival package with that by survdiff() function in survival package
survdiff.result  = survdiff(Surv(futime, fustat) ~ rx, data = ovarian)
FastLRtest.result = FastLRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2)
print(list(survdiff.chisq = survdiff.result[['chisq']], FastLRtest.chisq = FastLRtest.result ^ 2))
# $survdiff.chisq
# [1] 1.06274
# 
# $FastLRtest.chisq
# [1] 1.06274

# Computing time comparison
microbenchmark(
  survdiff.result  = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
  FastLRtest.result = FastLRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2)
)
# Unit: microseconds
#               expr     min      lq      mean  median       uq      max neval
#    survdiff.result 598.301 668.101 827.18999 742.851 835.8005 1601.901   100
#  FastLRtest.result  56.002  72.750  96.40398  89.151 102.3005  259.701   100

## Hazard ratio estimation
# Compare the hazard ratio estimation by FastSurvival package with that by coxph() function in survival package
coxph.result          = data.frame(method = 'CoxPH', HR = exp(-coxph(Surv(futime, fustat) ~ rx, data = ovarian)[['coefficients']][[1]]))
FastHRest.PY.result   = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'PY')
FastHRest.MST.result  = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'MST')
FastHRest.Pike.result = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'Pike')
FastHRest.Peto.result = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'Peto')
FastHRest.LR.result   = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'LR')
print(rbind(coxph.result, FastHRest.PY.result, FastHRest.MST.result, FastHRest.Pike.result, FastHRest.Peto.result, FastHRest.LR.result))
#   method       HR
# 1  CoxPH 1.815535
# 2     PY 1.832895
# 3    MST 1.363636
# 4   Pike 1.810070
# 5   Peto 1.825070
# 6     LR 1.813369

# Computing time comparison
microbenchmark(
  coxph.result          = data.frame(method = 'CoxPH', HR = exp(-coxph(Surv(futime, fustat) ~ rx, data = ovarian)[['coefficients']][[1]])),
  FastHRest.PY.result   = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'PY'),
  FastHRest.MST.result  = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'MST'),
  FastHRest.Pike.result = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'Pike'),
  FastHRest.Peto.result = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'Peto'),
  FastHRest.LR.result   = FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'LR')
)
# Unit: microseconds
#                   expr    min      lq     mean  median      uq    max neval
#           coxph.result 1054.5 1261.40 1602.396 1407.90 1640.55 3398.7   100
#    FastHRest.PY.result  110.9  139.85  185.484  154.00  174.40  686.1   100
#   FastHRest.MST.result  145.2  172.10  215.489  193.80  221.20  869.1   100
#  FastHRest.Pike.result  153.3  183.90  292.961  207.70  237.20 6543.7   100
#  FastHRest.Peto.result  154.9  183.15  225.609  208.15  237.40  552.6   100
#    FastHRest.LR.result  150.4  182.80  246.916  205.35  233.00  886.2   100
```


