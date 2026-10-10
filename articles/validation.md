# Validation of FastSurvival

## Purpose

FastSurvival is built for speed, but a fast estimator is only useful if
it returns the same answer as the established implementation. This
vignette checks the numerical agreement between each FastSurvival
function and a reference, on a real clinical-trial dataset. We use the
`gbsg` data from the survival package, the German Breast Cancer Study
Group cohort of 686 patients, with recurrence-free survival time
`rfstime` (in days), the event indicator `status`, and the treatment
indicator `hormon` (0 = no hormonal therapy, 1 = hormonal therapy).
Throughout we take the no-hormone arm as the control (`control = 0`) and
report one-sided tests of treatment benefit (`side = 1`).

The references are: the survival package for the Kaplan-Meier estimate,
the log-rank test, the Cox hazard ratio, milestone survival, and the
median survival comparison; survRM2 for the restricted mean survival
time; survAH for the average hazard with survival weight; nph for the
weighted log-rank test and the max-combo test; and nphRCT for the robust
modestly-weighted test. The window mean survival time and the weighted
Kaplan-Meier statistic, which have no direct package counterpart on
CRAN, are checked against a Kaplan-Meier integral computed from
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) and against
the restricted mean survival time identity, respectively. The same
comparisons form the basis of the automated test suite shipped with the
package.

``` r

library(FastSurvival)
```

## Reference data

``` r

library(survival)

# German Breast Cancer Study Group cohort
str(gbsg[, c("rfstime", "status", "hormon")])
#> 'data.frame':    686 obs. of  3 variables:
#>  $ rfstime: int  1838 403 1603 177 1855 842 293 42 564 1093 ...
#>  $ status : int  0 1 0 0 0 1 1 0 1 1 ...
#>  $ hormon : int  0 0 0 0 1 0 0 1 1 0 ...
table(hormon = gbsg$hormon)
#> hormon
#>   0   1 
#> 440 246
```

## Kaplan-Meier survival

[`survfit_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)
evaluates the Kaplan-Meier estimate at a single time point. We compare
against `summary(survfit(...))` at the same time (1000 days).

``` r

# survfit_fast assumes time-sorted input; the C++ core groups tied times in the
# survival::survfit risk-set convention, so the data only need ordering by time.
ord  <- order(gbsg$rfstime)
fast <- survfit_fast(gbsg$rfstime[ord], gbsg$status[ord],
                     t_eval = 1000, conf.type = "log-log")

fit  <- survfit(Surv(rfstime, status) ~ 1, data = gbsg)
ref  <- summary(fit, times = 1000)

data.frame(
  quantity = c("survival", "std.err"),
  fast     = c(unclass(fast)["surv"], unclass(fast)["std.err"]),
  survival = c(ref$surv, ref$std.err),
  row.names = NULL
)
#>   quantity       fast   survival
#> 1 survival 0.65775040 0.65775040
#> 2  std.err 0.01905146 0.01905146
```

## Log-rank test

[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
with `side = 1` returns the signed Z-score of the log-rank test (Mantel,
1966). Its square is the chi-square statistic from
[`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html).

``` r

fast_lr <- survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                         control = 0, side = 1)

ref_lr  <- survdiff(Surv(rfstime, status) ~ hormon, data = gbsg)

c(fast = as.numeric(fast_lr)^2, survival = ref_lr$chisq)
#>     fast survival 
#> 8.564781 8.564781
```

## Weighted log-rank test

[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
also computes Fleming-Harrington G(rho, gamma) weighted log-rank tests
through the `weight = "fh"` argument. The
[nph](https://cran.r-project.org/package=nph) package provides
`logrank.test()` with the same family, so we compare the chi-square
statistic (the squared one-sided Z) across several weight choices.

``` r

fh_grid <- data.frame(rho = c(0, 1, 0, 1), gamma = c(1, 0, 0, 1))

do.call(rbind, lapply(seq_len(nrow(fh_grid)), function(i) {
  r <- fh_grid$rho[i]
  g <- fh_grid$gamma[i]
  fast <- as.numeric(survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                                   control = 0, side = 1,
                                   weight = "fh", rho = r, gamma = g))
  nph_chisq <- nph::logrank.test(gbsg$rfstime, gbsg$status, gbsg$hormon,
                                 rho = r, gamma = g)$test$Chisq
  data.frame(rho = r, gamma = g, fast = fast^2, nph = nph_chisq)
}))
#>   rho gamma     fast      nph
#> 1   0     1 5.110660 5.110660
#> 2   1     0 8.713791 8.713791
#> 3   0     0 8.564781 8.564781
#> 4   1     1 5.881310 5.881310
```

The two implementations agree to numerical precision across the weight
family, including the ordinary log-rank test recovered at
`rho = 0, gamma = 0`.

## Cox hazard ratio

[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
returns the Pike-Halley Estimator (Homma, 2025), a closed-form
approximation to the maximizer of the Cox partial likelihood with the
Breslow method for ties. The reference is therefore
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) with
`ties = "breslow"` (its default is the Efron method). We report the log
hazard ratio from both; the sign convention is the same and `side` does
not affect the point estimate.

``` r

fast_hr <- coxph_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                      control = 0, side = 1)

ref_cox <- coxph(Surv(rfstime, status) ~ hormon, data = gbsg,
                 ties = "breslow")

c(fast = unclass(fast_hr)["coef"], cox = unname(coef(ref_cox)))
#>  fast.coef        cox 
#> -0.3638987 -0.3638988
cox_diff <- abs(unname(unclass(fast_hr)["coef"]) - unname(coef(ref_cox)))
```

Because the Pike-Halley Estimator is a closed-form approximation rather
than the exact partial-likelihood maximizer, the two log hazard ratios
are not expected to be identical. Their absolute difference here is
2.1e-09, which reflects the approximation and is not a sign of error.

With `strata`,
[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
approximates the stratified Cox model, in which each stratum has its own
baseline hazard. Stratifying by tumor grade, the estimate is compared
with [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) using
[`strata()`](https://rdrr.io/pkg/survival/man/strata.html) and the
Breslow method for ties, which is the partial likelihood that the
Pike-Halley Estimator targets.

``` r

fast_hr_s <- coxph_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                        control = 0, strata = gbsg$grade)

ref_cox_s <- coxph(Surv(rfstime, status) ~ hormon + strata(grade),
                   data = gbsg, ties = "breslow")

c(fast = unclass(fast_hr_s)["coef"], cox = unname(coef(ref_cox_s)))
#>  fast.coef        cox 
#> -0.3396679 -0.3396679
```

## Restricted mean survival time

[`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
integrates the Kaplan-Meier survival curve up to a horizon (Royston and
Parmar, 2013). We compare the two-group RMST difference against
[`survRM2::rmst2()`](https://rdrr.io/pkg/survRM2/man/rmst2.html) at 1000
days.

``` r

library(survRM2)

tau <- 1000

fast_rmst <- rmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                       control = 0, tau = tau, side = 1)

ref_rmst  <- rmst2(time = gbsg$rfstime, status = gbsg$status,
                   arm = gbsg$hormon, tau = tau)

c(fast = unclass(fast_rmst)["diff"],
  survRM2 = ref_rmst$unadjusted.result[1, 1])
#> fast.diff   survRM2 
#>  51.14447  51.14447
```

## Window mean survival time

[`wmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md)
integrates the Kaplan-Meier curve between a lower and an upper window
limit, generalizing the restricted mean survival time (Paukner and
Chappell, 2021). CRAN has no dedicated WMST package, so we validate the
per-group windowed area directly against a Kaplan-Meier step-function
integral computed from
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) over the
same window. A full comparison against the `survWMST` package, which is
distributed on GitHub, is provided in `tools/compare_wmst_survwmst.R` of
the package’s GitHub repository.

``` r

tau1 <- 200
tau2 <- 1000

fast_wm <- wmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, side = 1, tau1 = tau1, tau2 = tau2)

# Integral of the Kaplan-Meier step function over [lo, hi] for one group.
km_window_area <- function(gi, lo, hi) {
  sf <- survfit(Surv(rfstime, status) ~ 1, data = gbsg[gbsg$hormon == gi, ])
  tk <- c(0, sf$time)
  sk <- c(1, sf$surv)
  brk  <- sort(unique(c(lo, hi, sf$time[sf$time > lo & sf$time < hi])))
  area <- 0
  for (b in seq_len(length(brk) - 1L)) {
    u  <- brk[b]
    su <- sk[max(which(tk <= u))]
    area <- area + su * (brk[b + 1L] - u)
  }
  area
}

data.frame(
  group    = c("control", "treatment"),
  fast     = unclass(fast_wm)[c("wmst.control", "wmst.treatment")],
  survfit  = c(km_window_area(0, tau1, tau2), km_window_area(1, tau1, tau2)),
  row.names = NULL
)
#>       group     fast  survfit
#> 1   control 628.0064 628.0064
#> 2 treatment 678.3393 678.3393
```

## Weighted Kaplan-Meier test

[`wkm_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md)
computes the weighted Kaplan-Meier (Pepe-Fleming) test, the weighted
integral of the difference between the two Kaplan-Meier curves. With the
constant weight the weighted difference reduces exactly to the
difference in restricted mean survival time over the observed range,
which gives a clean internal check against
[`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
evaluated at the largest observed time.

``` r

tmax <- max(gbsg$rfstime)

fast_wk <- wkm_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                    control = 0, side = 1, weight = "constant")

# tau = tmax extends past the shorter group's follow-up on purpose, to cover
# the whole observed range as the constant-weight test does; rmst_fast()
# warns about this, so the warning is suppressed here.
ref_rm  <- suppressWarnings(
  rmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
            control = 0, tau = tmax, side = 1))

c(wkm.constant = unclass(fast_wk)["wdiff"], rmst = unclass(ref_rm)["diff"])
#> wkm.constant.wdiff          rmst.diff 
#>           256.7207           256.7207
```

With the default Pepe-Fleming weight,
[`wkm_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md)
reproduces the weighted Kaplan-Meier statistic of the `nphsim` package.
Because `nphsim` is distributed only on GitHub, that comparison is kept
in `tools/compare_wkm_nphsim.R` (which also bundles a survival-only
reproduction of the statistic) and in the test suite, rather than in
this vignette.

## Milestone survival

[`milestone_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md)
compares Kaplan-Meier survival between two groups at a milestone
timepoint, here with the complementary log-log intervals and the MOVER
interval of the difference of Tang (2021). The per-group survival
probabilities match those from
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).

``` r

tstar <- 1000

fast_ms <- milestone_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                          control = 0, tau = tstar, method = "loglog",
                          side = 1)

fit_g <- survfit(Surv(rfstime, status) ~ hormon, data = gbsg)
ref_g <- summary(fit_g, times = tstar)

data.frame(
  group    = c("control", "treatment"),
  fast     = fast_ms$surv[c("control", "treatment")],
  survival = ref_g$surv,
  row.names = NULL
)
#>       group      fast  survival
#> 1   control 0.6208762 0.6208762
#> 2 treatment 0.7230082 0.7230082
```

## Median survival time

[`medsurv_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md)
estimates the Kaplan-Meier median survival time and, for two groups,
their difference. The point estimate is the Kaplan-Meier median (the
first time at which the product-limit curve reaches 0.5), the same
convention as
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html), so the
per-group medians match those reported by
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) exactly.

``` r

fast_med <- medsurv_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                         control = 0, side = 1, method = "nph")

fit_med <- survfit(Surv(rfstime, status) ~ hormon, data = gbsg)
med_ref <- summary(fit_med)$table[, "median"]

data.frame(
  group    = c("control", "treatment"),
  fast     = unclass(fast_med)[c("median.control", "median.treatment")],
  survival = c(med_ref[1], med_ref[2]),
  row.names = NULL
)
#>       group fast survival
#> 1   control 1528     1528
#> 2 treatment 2018     2018
```

We use [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) as
the reference here, rather than
[`nph::nphparams()`](https://rdrr.io/pkg/nph/man/nphparams.html),
because the two define the median on different survival curves:
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) and
[`medsurv_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md)
use the Kaplan-Meier (product-limit) median, whereas
[`nph::nphparams()`](https://rdrr.io/pkg/nph/man/nphparams.html) reads
the median off the Nelson-Aalen curve, S(t) = exp(-H(t)), which lies
above the Kaplan-Meier curve and so can cross 0.5 at a later time on
tied data. The `method = "nph"` standard error reproduces the
[`nph::nphparams()`](https://rdrr.io/pkg/nph/man/nphparams.html)
standard error to numerical precision when the two medians coincide;
that comparison is run on tie-free scenarios in
`tools/compare_medsurv_nphparams.R` and in the package test suite.

## Average hazard with survival weight

[`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md)
computes the average hazard with survival weight of Uno and Horiguchi.
We compare the per-group average hazard against
[`survAH::ah2()`](https://rdrr.io/pkg/survAH/man/ah2.html).

``` r

library(survAH)

tau <- 1000

fast_ah <- ahsw_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, tau = tau, side = 1)

ref_ah  <- ah2(time = gbsg$rfstime, status = gbsg$status,
               arm = gbsg$hormon, tau = tau)

data.frame(
  quantity = c("AH (control)", "AH (treatment)"),
  fast     = unclass(fast_ah)[c("ah.ctrl", "ah.trt")],
  survAH   = c(ref_ah$ah["AH (arm0)", "Est."],
               ref_ah$ah["AH (arm1)", "Est."]),
  row.names = NULL
)
#>         quantity         fast       survAH
#> 1   AH (control) 0.0004585857 0.0004585857
#> 2 AH (treatment) 0.0003155277 0.0003155277
```

The average hazard is the ratio of the cumulative event probability to
the restricted mean survival time, so for the `gbsg` data, where
follow-up is measured in days, the values are on the order of 1e-04 per
day in both groups.

## Average hazard ratio

[`ahr_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md)
computes the Kalbfleisch-Prentice average hazard ratio between two
groups over a restricted interval. The reference implementation is
`ahrKM()` from the AHR package, but that package has been archived on
CRAN, so here we check the point estimates against a direct
survival-based computation: each group’s Kaplan-Meier curve from
[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html), integrated
to form the group shares of the total hazard.

``` r

tau <- 1000

fast_ahr <- ahr_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, tau = tau, side = 1)

g  <- gbsg$hormon
ev <- c(gbsg$rfstime[g == 0 & gbsg$status == 1],
        gbsg$rfstime[g == 1 & gbsg$status == 1])
grid <- sort(unique(c(0, ev[ev <= tau], tau)))

km_on_grid <- function(gi) {
  sf <- survfit(Surv(rfstime, status) ~ 1, data = gbsg[g == gi, ])
  approxfun(sf$time, sf$surv, method = "constant",
            yleft = 1, rule = 2, f = 0)(grid)
}

S0  <- km_on_grid(0)
S1  <- km_on_grid(1)
m   <- length(grid)
dS0 <- S0 - c(1, S0[-m])
GL  <- S0[m] * S1[m]
ref_theta_ctrl <- -sum(S1 * dS0) / (1 - GL)

data.frame(
  quantity  = c("theta (control)", "theta (treatment)", "AHR"),
  fast      = c(fast_ahr$theta[[1]], fast_ahr$theta[[2]], fast_ahr$ahr),
  reference = c(ref_theta_ctrl, 1 - ref_theta_ctrl,
                (1 - ref_theta_ctrl) / ref_theta_ctrl),
  row.names = NULL
)
#>            quantity      fast reference
#> 1   theta (control) 0.5974108 0.5974108
#> 2 theta (treatment) 0.4025892 0.4025892
#> 3               AHR 0.6738899 0.6738899
```

The point estimates match the survival-based reference to numerical
precision. The average hazard ratio was additionally cross-checked
against the archived AHR package, the reference implementation used by
Dormuth et al. (2024): the two group shares, the average hazard ratio,
the variances, and both test statistics agree to within about 1e-14.
That external comparison is reproducible with the
`tools/compare_ahr_ahrKM.R` script of the package’s GitHub repository.

## Max-combo test

[`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md)
computes the max-combo test, the most extreme of a set of
Fleming-Harrington weighted log-rank statistics (Karrison, 2016). The
`nph` package provides `logrank.maxtest()`, whose default weight set is
FH(0, 0), FH(0, 1), and FH(1, 0). We request the same three weights from
[`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md)
and compare the component Z-scores. The individual Z-scores are kept in
the `z` attribute, signed so that a negative value indicates benefit for
the treatment group under the package convention. Because `nph` orients
the contrast in the opposite direction, the signs are mirrored, so we
compare absolute values.

``` r

rho   <- c(0, 0, 1)
gamma <- c(0, 1, 0)

fast_mc <- maxcombo_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                         control = 0, side = 1, rho = rho, gamma = gamma)

nph_mc <- nph::logrank.maxtest(gbsg$rfstime, gbsg$status, gbsg$hormon)

data.frame(
  weight = c("FH(0,0)", "FH(0,1)", "FH(1,0)"),
  fast   = abs(attr(fast_mc, "z")),
  nph    = abs(nph_mc$tests$z),
  row.names = NULL
)
#>    weight     fast      nph
#> 1 FH(0,0) 2.926565 2.926565
#> 2 FH(0,1) 2.260677 2.260677
#> 3 FH(1,0) 2.951913 2.951913
```

The component statistics agree in absolute value. The final max-combo
statistic is the largest of these, which for the one-sided test is the
most extreme component Z. The p-values are not compared here, because
[`nph::logrank.maxtest()`](https://rdrr.io/pkg/nph/man/logrank.maxtest.html)
reports two-sided p-values and the call above is one-sided.

## Robust modestly-weighted log-rank test

[`rmw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md)
computes the robust modestly-weighted (rMW) test of Magirr and Öhrn, the
maximum of the standard log-rank statistic and a single
modestly-weighted log-rank statistic. Each component is a weighted
log-rank test, so we check the two component Z-scores against
[nphRCT](https://cran.r-project.org/package=nphRCT), the implementation
used by the method’s authors. The modestly-weighted component uses the
survival-threshold parameterization, so `s_star = 1` recovers the
standard log-rank test and `s_star = 0.5` gives the modestly-weighted
component. As with the max-combo test the package orients a negative Z
toward treatment benefit while `nphRCT` uses the opposite sign, so we
compare absolute values.

``` r

gbsg_df <- data.frame(
  time  = gbsg$rfstime,
  event = gbsg$status,
  arm   = factor(ifelse(gbsg$hormon == 0, "control", "experimental"),
                 levels = c("control", "experimental"))
)

fit_rmw <- rmw_fast(gbsg_df$time, gbsg_df$event, gbsg_df$arm,
                    control = "control", side = 1, s_star = 0.5)

z_lr_nph <- nphRCT::wlrt(Surv(time, event) ~ arm, data = gbsg_df,
                         method = "mw", s_star = 1)$z
z_mw_nph <- nphRCT::wlrt(Surv(time, event) ~ arm, data = gbsg_df,
                         method = "mw", s_star = 0.5)$z

data.frame(
  component = c("log-rank (s_star = 1)", "modestly-weighted (s_star = 0.5)"),
  fast      = abs(attr(fit_rmw, "z")),
  nphRCT    = abs(c(z_lr_nph, z_mw_nph)),
  row.names = NULL
)
#>                          component     fast   nphRCT
#> 1            log-rank (s_star = 1) 2.926565 2.926565
#> 2 modestly-weighted (s_star = 0.5) 2.773749 2.773749
```

The component statistics agree in absolute value. The null correlation
of the two components, reported in the `corr` attribute, is checked
against an independent implementation in the package test suite. The
combined statistic and one-sided p-value then follow from the bivariate
normal distribution of the two components, evaluated with
[`mvtnorm::pmvnorm`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html).

## Analysis cutoffs

The simulation layer determines the calendar time of each analysis with
[`cutoff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md),
which combines a target number of events, a planned calendar time, a
maximum calendar time, a minimum time after the previous look, and a
minimum follow-up after a given number of enrolled subjects. The same
rule is implemented one simulated trial at a time by
`get_analysis_date()` of the
[simtrial](https://cran.r-project.org/package=simtrial) package. We
compare the two on 50 simulated trials for a look at “150 events, but
not before month 18 and not before 6 months after the 250th enrolled
subject, and in any case by month 30”.

``` r

sim_c <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
                      a.rate = 300 / 12, e.median = list(12, 16),
                      d.hazard = 0.01, seed = 2026)
cut_c <- cutoff_fast(sim_c, event.looks = 150, time.looks = 18,
                     max.time = 30, min.enrolled = 250, min.followup = 6)
ref_c <- vapply(seq_len(50), function(s) {
  d <- sim_c[sim_c$sim == s, ]
  x <- data.frame(enroll_time = d$accrual_time,
                  cte = d$accrual_time + d$tte,
                  fail = d$event, stratum = "All")
  simtrial::get_analysis_date(x, planned_calendar_time = 18,
                              target_event_overall = 150,
                              max_extension_for_target_event = 30,
                              min_n_overall = 250, min_followup = 6)
}, numeric(1))
c(max_abs_difference = max(abs(cut_c[, 1] - ref_c)))
#> max_abs_difference 
#>                  0
```

The cutoffs agree exactly. The two functions differ only when an event
target is never met in the simulated data, which
[`cutoff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md)
reports as an unreached look rather than as the time of the last
observed event.

## Treatment switching

[`switch_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md)
changes only the outcomes after each subject’s switch. For an
exponential survival time with hazard `lambda`, a switch at time `s`
after accrual that multiplies the remaining time by `f` gives the
survival function `exp(-lambda t)` for `t <= s` and
`exp(-lambda s - lambda (t - s) / f)` for `t > s`. We let every control
subject still alive at calendar month 12 switch with `f = 2`, so that
`s = 12 - a` for a subject accrued at time `a`, and compare the
empirical survival of the control group with the average of this
survival function over the accrual times.

``` r

lam   <- log(2) / 12
sim_s <- simdata_fast(nsim = 200, n = c(250, 250), a.time = c(0, 6),
                      a.prop = 1, e.hazard = list(lam, lam / 1.5),
                      seed = 2027)
sw    <- switch_fast(sim_s, group = 1, when = "cutoff",
                     cutoff = rep(12, 200), aft.factor = 2)
ctl   <- sw$group == 1
s_sw  <- 12 - sw$accrual_time[ctl]
t_grid <- c(6, 12, 18, 24, 36)
data.frame(
  t         = t_grid,
  empirical = round(sapply(t_grid, function(t) mean(sw$surv_time[ctl] > t)), 4),
  analytic  = round(sapply(t_grid, function(t) {
    mean(ifelse(t <= s_sw, exp(-lam * t),
                exp(-lam * s_sw - lam * (t - s_sw) / 2)))
  }), 4)
)
#>    t empirical analytic
#> 1  6    0.7052   0.7071
#> 2 12    0.5456   0.5460
#> 3 18    0.4594   0.4591
#> 4 24    0.3860   0.3861
#> 5 36    0.2728   0.2730
```

The empirical survival matches the analytic survival function to Monte
Carlo error (the standard error is at most about 0.002 with 50,000
control subjects). The vignette on treatment switching also checks the
switching model with the rank-preserving structural failure time
estimator of the rpsftm package.

## Summary

Across all functions the FastSurvival results reproduce the reference
values on the `gbsg` data. The point estimates and test statistics agree
to numerical precision for the Kaplan-Meier, log-rank, weighted
log-rank, RMST, milestone, median, and average-hazard quantities; the
window mean survival time matches a Kaplan-Meier integral and the
weighted Kaplan-Meier statistic matches the RMST identity; the average
hazard ratio matches a survival-based reference; and the closed-form Cox
hazard ratio, with or without strata, agrees with the partial-likelihood
maximizer to the order expected for the Pike-Halley approximation. In
the simulation layer, the analysis cutoffs agree with simtrial and the
switched survival times follow their analytic distribution. This
agreement is verified continuously by the package test suite.

## References

Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from
incomplete observations. *Journal of the American Statistical
Association*, 53(282), 457-481.

Mantel, N. (1966). Evaluation of survival data and two new rank order
statistics arising in its consideration. *Cancer Chemotherapy Reports*,
50(3), 163-170.

Fleming, T. R., & Harrington, D. P. (1991). *Counting Processes and
Survival Analysis*. New York: John Wiley & Sons.

Cox, D. R. (1972). Regression models and life-tables. *Journal of the
Royal Statistical Society. Series B (Methodological)*, 34(2), 187-202.

Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
estimator. *Manuscript under review.*

Royston, P., & Parmar, M. K. B. (2013). Restricted mean survival time:
an alternative to the hazard ratio for the design and analysis of
randomized trials with a time-to-event outcome. *BMC Medical Research
Methodology*, 13, 152.

Paukner, M., & Chappell, R. (2021). Window mean survival time.
*Statistics in Medicine*, 40(25), 5521-5533.

Pepe, M. S., & Fleming, T. R. (1989). Weighted Kaplan-Meier statistics:
a class of distance tests for censored survival data. *Biometrics*,
45(2), 497-507.

Pepe, M. S., & Fleming, T. R. (1991). Weighted Kaplan-Meier statistics:
large sample and optimality considerations. *Journal of the Royal
Statistical Society. Series B (Methodological)*, 53(2), 341-352.

Tang, Y. (2021). Some new confidence intervals for Kaplan-Meier based
estimators from one and two sample survival data. *Statistics in
Medicine*, 40(23), 4961-4976.

Uno, H., Claggett, B., Tian, L., et al. (2014). Moving beyond the hazard
ratio in quantifying the between-group difference in survival analysis.
*Journal of Clinical Oncology*, 32(22), 2380-2385.

Uno, H., & Horiguchi, M. (2023). Ratio and difference of average hazard
with survival weight: new measures to quantify survival benefit of new
therapy. *Statistics in Medicine*, 42(7), 936-952.

Kalbfleisch, J. D., & Prentice, R. L. (1981). Estimation of the average
hazard ratio. *Biometrika*, 68(1), 105-112.

Dormuth, I., Pauly, M., Rauch, G., & Herrmann, C. (2024). Sample size
calculation under nonproportional hazards using average hazard ratios.
*Biometrical Journal*, 66(6), e202300271.

Karrison, T. G. (2016). Versatile tests for comparing survival curves
based on weighted log-rank statistics. *The Stata Journal*, 16(3),
678-690.

Magirr, D., & Burman, C.-F. (2019). Modestly weighted logrank tests.
*Statistics in Medicine*, 38(20), 3782-3790.

Magirr, D., & Öhrn, F. (2026). Robust modestly weighted log-rank tests.
*Pharmaceutical Statistics*, 25(1), e70066.
