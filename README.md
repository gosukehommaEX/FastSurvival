# FastSurvival

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/FastSurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/FastSurvival/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

FastSurvival provides fast alternatives to the standard survival analysis
functions in the
[survival](https://cran.r-project.org/package=survival) package.
Every function is designed for repeated evaluation inside large simulation
loops — adaptive sample-size re-estimation, probability-of-success
calculations, regional consistency evaluation in multi-regional trials —
where the standard iterative or object-building overhead of `survfit()`,
`survdiff()`, and `coxph()` becomes a computational bottleneck.

## Functions

| Function | Replaces | Speed gain |
|----------|----------|------------|
| `survfit_fast()` | `survfit()` + `summary()` at a single time point | ~63x |
| `survdiff_fast()` | `survdiff()` | ~8x |
| `coxph_fast()` | `coxph()` (point estimate + Wald CI) | ~17x |
| `simdata_fast()` | Custom simulation scripts | — |

Speed gains are median times from 1,000 microbenchmark replicates on a
typical phase-3 trial dataset (n = 600, event rate 80%). Results vary by
hardware and sample size.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("gosukehommaEX/FastSurvival")
```

## Quick start

```r
library(FastSurvival)
library(survival)

# ----------------------------------------------------------------
# survfit_fast(): Kaplan-Meier at a single time point
# ----------------------------------------------------------------
ord <- order(ovarian$futime)
t_s <- ovarian$futime[ord]
e_s <- ovarian$fustat[ord]

survfit_fast(t_s, e_s, t_eval = 500, conf.type = "log")
#>       surv    std.err      lower      upper
#> 0.59607843 0.09992615 0.42914950 0.82793874

# ----------------------------------------------------------------
# survdiff_fast(): log-rank test
# ----------------------------------------------------------------
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
              control = 1, side = 2)
#> [1] 1.06274

# ----------------------------------------------------------------
# coxph_fast(): hazard ratio via PiHE (closed-form)
# ----------------------------------------------------------------
coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
#>       coef  exp(coef)   se(coef)  lower .95  upper .95
#> -0.5963800  0.5508019  0.5868442  0.1743704  1.7398754

# ----------------------------------------------------------------
# simdata_fast(): clinical trial simulation
# ----------------------------------------------------------------

# Two-group trial, simple exponential, no dropout
df <- simdata_fast(
  nsim     = 1000,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.median = list(18, 24),
  seed     = 1
)
head(df)

# Two-group trial, piecewise exponential (delayed treatment effect)
df2 <- simdata_fast(
  nsim     = 1000,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.hazard = list(c(0.08, 0.08), c(0.08, 0.04)),
  e.time   = c(0, 6, Inf),
  seed     = 2
)
```

## Design principles

**survfit_fast** evaluates the Kaplan-Meier estimator at a single specified
time point using binary search (`findInterval()`), processing only the
observations up to that point. This makes it orders of magnitude faster than
`survfit()` when the same sorted data are evaluated repeatedly at a fixed
landmark time inside a simulation loop.

**survdiff_fast** computes the log-rank statistic using the same
rank-based at-risk counting as `survdiff()` but avoids the overhead of
S3 object construction and formula parsing. It returns either a one-sided
Z-score (`side = 1`) or a two-sided chi-square statistic (`side = 2`).

**coxph_fast** implements the Pike-Halley Estimator (PiHE) proposed by
Homma (2025). PiHE anchors at the Pike closed-form estimate and applies a
single analytic Halley correction to the Cox partial likelihood score. The
residual error relative to the Cox MLE satisfies
|PiHE - Cox| = O_p(n^{-3/2}), compared with O_p(n^{-1/2}) for Peto and
Pike. On the `pharmacoSmoking` dataset (tie rate 77.5%), PiHE reproduces
the Breslow-based Cox MLE to within 1.10e-08. The Wald confidence interval
uses the observed information at the Pike anchor as the variance estimate.

**simdata_fast** generates individual patient data for one- or two-group
time-to-event trials. Accrual times follow a piecewise uniform distribution.
Survival and dropout times follow either a simple or piecewise exponential
distribution, selected automatically based on whether a scalar or vector
hazard is supplied. All random number generation uses
[dqrng](https://cran.r-project.org/package=dqrng) for speed.

## Using the functions together

A typical simulation workflow combines all four functions:

```r
library(FastSurvival)

results <- vector("list", 1000L)

df <- simdata_fast(
  nsim     = 1000,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.hazard = list(0.08, 0.05),
  d.hazard = list(0.01, 0.01),
  seed     = 42
)

for (s in seq_len(1000L)) {
  d <- df[df$sim == s, ]
  results[[s]] <- coxph_fast(d$tte, d$event, d$group, control = 1)
}

hr_mat <- do.call(rbind, results)
colMeans(hr_mat)
```

## Reference

Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
estimator. *Manuscript under review.*

## License

MIT © 2025 Gosuke Homma
