# Speed comparison

## Purpose

FastSurvival is designed for repeated evaluation inside large simulation
loops. This vignette shows how to benchmark each estimation and testing
function against an established reference and reports representative
results. The benchmark code is shown but not executed when the vignette
is built, because timing many microbenchmark replicates would exceed the
build-time limits. To reproduce the numbers, run the code blocks
interactively. The same code is collected in the
`tools/benchmark_speed.R` script shipped in the package sources.

The reported figures are median times from microbenchmark replicates on
a single desktop machine. Absolute timings depend on hardware, sample
size, and event rate, so the ratios matter more than the raw values.

``` r

library(FastSurvival)
library(survival)
library(microbenchmark)
```

## Setup

The key to the speed gain is that the analysis functions accept
pre-sorted vectors. Inside a simulation loop the data are sorted once
and reused, so the sort cost is paid a single time rather than on every
call. We build a single two-group dataset of 500 subjects with
[`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
and prepare the sorted vectors, the binary arm indicator, and the
restriction horizon used by the time-restricted methods.

``` r

dataset <- simdata_fast(
  nsim     = 1,
  n        = 500,
  a.time   = c(0, 12.5),
  a.rate   = 40,
  e.median = list(5.811, 4.3),
  seed     = 1
)

# Sort once and reuse, the intended pattern for the pre-sorted fast path.
ord <- order(dataset$tte)
t_s <- dataset$tte[ord]
e_s <- dataset$event[ord]
g_s <- dataset$group[ord]

# Control is group 1, treatment is group 2.
arm <- as.integer(dataset$group == 2)

# Restriction horizon within both arms' follow-up.
tau <- floor(min(tapply(t_s, g_s, max)))

# Factor arm for the nphRCT reference used in the rmw_fast benchmark.
df_rmw <- data.frame(
  tte   = dataset$tte,
  event = dataset$event,
  arm   = factor(ifelse(dataset$group == 1, "control", "treatment"),
                 levels = c("control", "treatment"))
)
```

## survfit_fast vs survfit + summary

``` r

microbenchmark(
  fast = survfit_fast(t_s, e_s, t_eval = tau, presorted = TRUE),
  ref  = summary(survfit(Surv(tte, event) ~ 1, data = dataset), times = tau),
  times = 1000
)
```

## survdiff_fast vs survdiff

``` r

microbenchmark(
  fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 1, presorted = TRUE),
  ref  = survdiff(Surv(tte, event) ~ group, data = dataset),
  times = 1000
)
```

## coxph_fast vs coxph

``` r

microbenchmark(
  fast = coxph_fast(t_s, e_s, g_s, control = 1, side = 1, presorted = TRUE),
  ref  = coxph(Surv(tte, event) ~ I(group == 2), data = dataset),
  times = 1000
)
```

## rmst_fast vs survRM2::rmst2

``` r

microbenchmark(
  fast = rmst_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                   presorted = TRUE),
  ref  = survRM2::rmst2(time = dataset$tte, status = dataset$event,
                        arm = arm, tau = tau),
  times = 1000
)
```

## survdiff_fast(weight = “fh”) vs nph::logrank.test

``` r

microbenchmark(
  fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 1,
                       weight = "fh", rho = 0, gamma = 1, presorted = TRUE),
  ref  = nph::logrank.test(dataset$tte, dataset$event, dataset$group,
                           rho = 0, gamma = 1),
  times = 1000
)
```

## wmst_fast vs survWMST::wmst

The window mean survival time is benchmarked against `wmst()` from the
survWMST package. survWMST is distributed on GitHub (pauknemj/survWMST),
not CRAN, so this benchmark is shown as a static block rather than a
live chunk, and the vignette carries no undeclared dependency. Install
survWMST with `remotes::install_github("pauknemj/survWMST")` and run the
block to reproduce it.

``` r

microbenchmark(
  fast = wmst_fast(t_s, e_s, g_s, control = 1, tau1 = 0, tau2 = tau,
                   side = 1, presorted = TRUE),
  ref  = survWMST::wmst(time = dataset$tte, status = dataset$event,
                        arm = arm, tau0 = 0, tau1 = tau),
  times = 1000
)
```

## milestone_fast vs survfit + summary

``` r

microbenchmark(
  fast = milestone_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                        presorted = TRUE),
  ref  = summary(survfit(Surv(tte, event) ~ group, data = dataset),
                 times = tau),
  times = 1000
)
```

## medsurv_fast vs nph::nphparams

``` r

microbenchmark(
  fast = medsurv_fast(t_s, e_s, g_s, control = 1, side = 1,
                      method = "nph", presorted = TRUE),
  ref  = nph::nphparams(time = dataset$tte, event = dataset$event,
                        group = as.integer(dataset$group == 2),
                        param_type = "Q", param_par = 0.5),
  times = 1000
)
```

## maxcombo_fast vs nph::logrank.maxtest

``` r

microbenchmark(
  fast = maxcombo_fast(t_s, e_s, g_s, control = 1, side = 1,
                       rho = c(0, 0, 1), gamma = c(0, 1, 0), presorted = TRUE),
  ref  = nph::logrank.maxtest(dataset$tte, dataset$event,
                              as.integer(dataset$group == 2)),
  times = 1000
)
```

## rmw_fast vs nphRCT::wlrt

[`rmw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md)
combines a standard and a modestly-weighted log-rank statistic, so the
reference computes both weighted log-rank components with nphRCT.

``` r

microbenchmark(
  fast = rmw_fast(t_s, e_s, g_s, control = 1, side = 1, s_star = 0.5,
                  presorted = TRUE),
  ref  = {
    nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw,
                 method = "mw", s_star = 1)
    nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw,
                 method = "mw", s_star = 0.5)
  },
  times = 1000
)
```

## wkm_fast vs nphsim::wkm.Stat

The weighted Kaplan-Meier (Pepe-Fleming) test is benchmarked against
`wkm.Stat()` from the nphsim package. nphsim is distributed on GitHub
(keaven/nphsim), not CRAN, so this benchmark is shown as a static block.
Install nphsim with `remotes::install_github("keaven/nphsim")` and run
the block to reproduce it.

``` r

microbenchmark(
  fast = wkm_fast(t_s, e_s, g_s, control = 1, side = 1, weight = "PF",
                  presorted = TRUE),
  ref  = nphsim::wkm.Stat(survival = dataset$tte, cnsr = 1 - dataset$event,
                          trt = ifelse(dataset$group == 1,
                                       "control", "experimental")),
  times = 1000
)
```

## ahsw_fast vs survAH::ah2

``` r

microbenchmark(
  fast = ahsw_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                   presorted = TRUE),
  ref  = survAH::ah2(time = dataset$tte, status = dataset$event,
                     arm = arm, tau = tau),
  times = 1000
)
```

## ahr_fast vs AHR::ahrKM

The Kalbfleisch-Prentice average hazard ratio is benchmarked against
`ahrKM()` from the AHR package, the reference implementation used by
Dormuth et al. (2024). Because AHR has been archived on CRAN, this
benchmark is shown as a static block rather than a live chunk. Install
AHR with `remotes::install_github("cran/AHR")` and run the block to
reproduce it.

``` r

microbenchmark(
  fast = ahr_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                  presorted = TRUE),
  ref  = AHR::ahrKM(tau, Surv(tte, event) ~ group, dataset),
  times = 1000
)
```

## Representative results

The table below summarizes representative median timings on the n = 500
two-group dataset generated above, with `presorted = TRUE` and one-sided
tests where applicable. The exact values will differ on your machine,
but the order of magnitude of the speedup is stable.

| Function | Replaces | Approximate speed gain |
|----|----|----|
| [`survfit_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md) | [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) + [`summary()`](https://rdrr.io/r/base/summary.html) at one time point | ~40x |
| [`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md) | [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) | ~40x |
| [`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md) | [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) (point estimate + Wald CI) | ~30x |
| [`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md) | [`survRM2::rmst2()`](https://rdrr.io/pkg/survRM2/man/rmst2.html) | ~40x |
| `survdiff_fast(weight = "fh")` | [`nph::logrank.test()`](https://rdrr.io/pkg/nph/man/logrank.test.html) | ~350x |
| [`wmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md) | [`survWMST::wmst()`](https://rdrr.io/pkg/survWMST/man/wmst.html) | ~560x |
| [`milestone_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md) | [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) + [`summary()`](https://rdrr.io/r/base/summary.html) at a milestone | ~20x |
| [`medsurv_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md) | [`nph::nphparams()`](https://rdrr.io/pkg/nph/man/nphparams.html) | ~30x |
| [`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md) | [`nph::logrank.maxtest()`](https://rdrr.io/pkg/nph/man/logrank.maxtest.html) | ~320x |
| [`rmw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md) | [`nphRCT::wlrt()`](https://rdrr.io/pkg/nphRCT/man/wlrt.html) (two components) | ~80x |
| [`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md) | [`survAH::ah2()`](https://rdrr.io/pkg/survAH/man/ah2.html) | ~410x |
| [`ahr_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md) | [`AHR::ahrKM()`](https://rdrr.io/pkg/AHR/man/ahrKM.html) | ~130x |

## Why it is faster

Each function avoids the overhead that the standard implementations
incur on every call. The standard functions parse a formula, build an S3
model object, and construct intermediate vectors before producing the
result, which is appropriate for interactive use but wasteful when the
same operation is repeated thousands of times. The FastSurvival
functions take plain vectors, do the core computation in a single C++
pass over the data, and return a lightweight numeric vector. When the
input is already sorted the sort cost is avoided entirely. In a
simulation loop these savings accumulate across every iteration.
