# Validation of FastSurvival

## Purpose

FastSurvival is built for speed, but a fast estimator is only useful if
it returns the same answer as the established implementation. This
vignette checks the numerical agreement between each FastSurvival
function and a reference from a well-known package: the survival package
for the Kaplan-Meier estimate, the log-rank test, the Cox hazard ratio,
and milestone survival; survRM2 for the restricted mean survival time;
and survAH for the average hazard with survival weight. The weighted
log-rank and max-combo tests are checked against the nph package. The
same comparisons form the basis of the automated test suite shipped with
the package.

``` r

library(FastSurvival)
```

## Kaplan-Meier survival

[`survfit_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)
evaluates the Kaplan-Meier estimate at a single time point. We compare
against `summary(survfit(...))` at the same time on the `ovarian`
dataset.

``` r

library(survival)

ord <- order(ovarian$futime)
t_s <- ovarian$futime[ord]
e_s <- ovarian$fustat[ord]

fast <- survfit_fast(t_s, e_s, t_eval = 500, conf.type = "log")

fit  <- survfit(Surv(futime, fustat) ~ 1, data = ovarian)
ref  <- summary(fit, times = 500)

data.frame(
  quantity = c("survival", "std.err"),
  fast     = c(unclass(fast)["surv"], unclass(fast)["std.err"]),
  survival = c(ref$surv, ref$std.err),
  row.names = NULL
)
#>   quantity       fast   survival
#> 1 survival 0.59607843 0.59607843
#> 2  std.err 0.09992615 0.09992615
```

## Log-rank test

[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
returns a two-sided chi-square statistic when `side = 2`, matching the
statistic from
[`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html).

``` r

fast_lr <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                         control = 1, side = 2)

ref_lr  <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)

c(fast = as.numeric(fast_lr), survival = ref_lr$chisq)
#>     fast survival 
#>  1.06274  1.06274
```

## Weighted log-rank test

[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
also computes Fleming-Harrington G(rho, gamma) weighted log-rank tests
through the `weight = "fh"` argument. The
[nph](https://cran.r-project.org/package=nph) package provides
[`logrank.test()`](https://rdrr.io/pkg/nph/man/logrank.test.html) with
the same family, so we compare the chi-square statistic across several
weight choices on the `ovarian` data.

``` r

fh_grid <- data.frame(rho = c(0, 1, 0, 1), gamma = c(1, 0, 0, 1))

do.call(rbind, lapply(seq_len(nrow(fh_grid)), function(i) {
  r <- fh_grid$rho[i]
  g <- fh_grid$gamma[i]
  fast <- as.numeric(survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                                   control = 1, side = 2,
                                   weight = "fh", rho = r, gamma = g))
  nph_chisq <- nph::logrank.test(ovarian$futime, ovarian$fustat, ovarian$rx,
                                 rho = r, gamma = g)$test$Chisq
  data.frame(rho = r, gamma = g, fast = fast, nph = nph_chisq)
}))
#>   rho gamma         fast          nph
#> 1   0     1 0.0001020735 0.0001020735
#> 2   1     0 1.6848546117 1.6848546117
#> 3   0     0 1.0627398613 1.0627398613
#> 4   1     1 0.0033228087 0.0033228087
```

The two implementations agree to numerical precision across the weight
family, including the ordinary log-rank test recovered at
`rho = 0, gamma = 0`.

## Cox hazard ratio

[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
returns the Pike-Halley Estimator, a closed-form approximation to the
Cox partial likelihood maximizer. On data with a moderate hazard ratio
the two agree closely. We report the log hazard ratio from both.

``` r

fast_hr <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)

ref_cox <- coxph(Surv(futime, fustat) ~ I(rx == 2), data = ovarian)

c(fast = unclass(fast_hr)["coef"], cox = unname(coef(ref_cox)))
#> fast.coef       cox 
#>  -0.59638  -0.59638
```

## Restricted mean survival time

[`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
integrates the Kaplan-Meier survival curve up to a horizon. We compare
the two-group RMST difference against
[`survRM2::rmst2()`](https://rdrr.io/pkg/survRM2/man/rmst2.html).

``` r

library(survRM2)

arm  <- as.integer(ovarian$rx == 2)
tau  <- 500

fast_rmst <- rmst_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                       control = 1, tau = tau)

ref_rmst  <- rmst2(time = ovarian$futime, status = ovarian$fustat,
                   arm = arm, tau = tau)

c(fast = unclass(fast_rmst)["diff"],
  survRM2 = ref_rmst$unadjusted.result[1, 1])
#> fast.diff   survRM2 
#>  98.88034  98.88034
```

## Milestone survival

[`milestone_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md)
compares Kaplan-Meier survival between two groups at a milestone
timepoint. The per-group survival probabilities match those from
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).

``` r

tstar <- 500

fast_ms <- milestone_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                          tau = tstar, method = "wald")

fit_g <- survfit(Surv(futime, fustat) ~ rx, data = ovarian)
ref_g <- summary(fit_g, times = tstar)

ms_surv <- fast_ms$surv

data.frame(
  group    = c("control", "treatment"),
  fast     = c(ms_surv["control"], ms_surv["treatment"]),
  survival = ref_g$surv,
  row.names = NULL
)
#>       group      fast  survival
#> 1   control 0.5384615 0.5384615
#> 2 treatment 0.6581197 0.6581197
```

## Average hazard with survival weight

[`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md)
computes the average hazard with survival weight of Uno and Horiguchi.
We compare the per-group average hazard against
[`survAH::ah2()`](https://rdrr.io/pkg/survAH/man/ah2.html).

``` r

library(survAH)

arm <- as.integer(ovarian$rx == 2)
tau <- 500

fast_ah <- ahsw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                     control = 1, tau = tau)

ref_ah  <- ah2(time = ovarian$futime, status = ovarian$fustat,
               arm = arm, tau = tau)

data.frame(
  quantity = c("AH (control)", "AH (treatment)"),
  fast     = unclass(fast_ah)[c("ah.ctrl", "ah.trt")],
  survAH   = c(ref_ah$ah["AH (arm0)", "Est."],
               ref_ah$ah["AH (arm1)", "Est."]),
  row.names = NULL
)
#>         quantity        fast      survAH
#> 1   AH (control) 0.001235076 0.001235076
#> 2 AH (treatment) 0.000723445 0.000723445
```

The average hazard is the ratio of the cumulative event probability to
the restricted mean survival time, so for the `ovarian` data, where
follow-up is measured in days, the values are on the order of 1e-03 per
day in both groups.

## Max-combo test

[`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md)
computes the max-combo test, the most extreme of a set of
Fleming-Harrington weighted log-rank statistics. The `nph` package
provides
[`logrank.maxtest()`](https://rdrr.io/pkg/nph/man/logrank.maxtest.html),
whose default weight set is FH(0, 0), FH(0, 1), and FH(1, 0). We request
the same three weights from
[`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md)
and compare the component Z-scores. The individual Z-scores are kept in
the `z` attribute, signed so that a negative value indicates benefit for
the treatment group under the package convention. Because `nph` orients
the contrast in the opposite direction, the signs are mirrored, so we
compare absolute values.

``` r

rho   <- c(0, 0, 1)
gamma <- c(0, 1, 0)

fast_mc <- maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                         control = 1, side = 2, rho = rho, gamma = gamma)

nph_mc <- nph::logrank.maxtest(ovarian$futime, ovarian$fustat, ovarian$rx)

data.frame(
  weight = c("FH(0,0)", "FH(0,1)", "FH(1,0)"),
  fast   = abs(attr(fast_mc, "z")),
  nph    = abs(nph_mc$tests$z),
  row.names = NULL
)
#>    weight       fast        nph
#> 1 FH(0,0) 1.03089275 1.03089275
#> 2 FH(0,1) 0.01010314 0.01010314
#> 3 FH(1,0) 1.29801950 1.29801950
```

The component statistics agree in absolute value. The final max-combo
statistic is the largest of these, which for the two-sided test is
`max(abs(Z))` and is therefore always non-negative. The p-values from
the two packages are close but not identical, since they use different
numerical methods for the joint distribution (a multivariate normal
integral here, a Bonferroni-style combination in `nph`).

## Summary

Across all functions the FastSurvival results reproduce the reference
values from the established packages. The point estimates and test
statistics agree to numerical precision for the Kaplan-Meier, log-rank,
RMST, milestone, and average-hazard quantities, and the closed-form Cox
hazard ratio agrees with the partial-likelihood maximizer to the order
expected for the Pike-Halley approximation. This agreement is verified
continuously by the package test suite.

## References

Tang, Y. (2021). Some new confidence intervals for Kaplan-Meier based
estimators from one and two sample survival data. *Statistics in
Medicine*, 40(23), 4961-4976.

Uno, H., & Horiguchi, M. (2023). Ratio and difference of average hazard
with survival weight: new measures to quantify survival benefit of new
therapy. *Statistics in Medicine*, 42(7), 936-952.

Karrison, T. G. (2016). Versatile tests for comparing survival curves
based on weighted log-rank statistics. *The Stata Journal*, 16(3),
678-690.
