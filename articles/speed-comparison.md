# Speed comparison

## Purpose

FastSurvival is designed for repeated evaluation inside large simulation
loops. This vignette shows how to benchmark each function against its
standard counterpart and reports representative results. The benchmark
code is shown but not executed when the vignette is built, because
timing many microbenchmark replicates would exceed the build-time
limits. To reproduce the numbers, run the code blocks interactively.

The reported figures are median times from microbenchmark replicates on
a single desktop machine. Absolute timings depend on hardware, sample
size, and event rate, so the ratios matter more than the raw values.

``` r

library(FastSurvival)
library(survival)
library(microbenchmark)
# Comparison packages used in the benchmarks below
library(survRM2)
library(nph)
```

## Setup

The key to the speed gain is that the analysis functions accept
pre-sorted vectors. Inside a simulation loop the data are sorted once
and reused, so the sort cost is paid a single time rather than on every
call. We build a single phase-3-sized dataset for the benchmarks.

``` r

set.seed(1)
n  <- 600
tt <- rexp(n, rate = 0.05)
ev <- rbinom(n, 1, 0.8)
gp <- rep(1:2, each = n / 2)

ord <- order(tt)
t_s <- tt[ord]
e_s <- ev[ord]
g_s <- gp[ord]
```

## survfit_fast vs survfit + summary

``` r

microbenchmark(
  fast = survfit_fast(t_s, e_s, t_eval = 20, presorted = TRUE),
  base = summary(survfit(Surv(tt, ev) ~ 1), times = 20),
  times = 1000
)
```

## survdiff_fast vs survdiff

``` r

microbenchmark(
  fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 2, presorted = TRUE),
  base = survdiff(Surv(tt, ev) ~ gp),
  times = 1000
)
```

## coxph_fast vs coxph

``` r

microbenchmark(
  fast = coxph_fast(t_s, e_s, g_s, control = 1, presorted = TRUE),
  base = coxph(Surv(tt, ev) ~ I(gp == 2)),
  times = 1000
)
```

## rmst_fast vs survRM2::rmst2

``` r

arm <- as.integer(gp == 2)
microbenchmark(
  fast = rmst_fast(t_s, e_s, g_s, control = 1, tau = 20, presorted = TRUE),
  base = survRM2::rmst2(time = tt, status = ev, arm = arm, tau = 20),
  times = 1000
)
```

## survdiff_fast(weight = “fh”) vs nph::logrank.test

``` r

microbenchmark(
  fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 2,
                       weight = "fh", rho = 0, gamma = 1, presorted = TRUE),
  base = nph::logrank.test(tt, ev, gp, rho = 0, gamma = 1),
  times = 1000
)
```

## Representative results

The table below summarizes representative median timings on a typical
phase-3 dataset (n = 600, event rate 80%) with `presorted = TRUE`. The
exact values will differ on your machine, but the order of magnitude of
the speedup is stable.

| Function | Replaces | Approximate speed gain |
|----|----|----|
| [`survfit_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md) | [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) + [`summary()`](https://rdrr.io/r/base/summary.html) at one time point | ~50x |
| [`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md) | [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) | ~30x |
| [`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md) | [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) (point estimate + Wald CI) | ~30x |
| [`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md) | [`survRM2::rmst2()`](https://rdrr.io/pkg/survRM2/man/rmst2.html) | ~50x |
| `survdiff_fast(weight = "fh")` | [`nph::logrank.test()`](https://rdrr.io/pkg/nph/man/logrank.test.html) | a few hundred x |

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

## References

Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
estimator. *Manuscript under review.*

Collett, D. (2014). *Modelling Survival Data in Medical Research* (3rd
ed.). Chapman and Hall/CRC.
