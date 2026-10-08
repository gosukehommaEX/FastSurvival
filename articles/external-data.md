# Using your own data generator

## Overview

[`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
generates piecewise-exponential survival and dropout times, subgroups,
multi-arm trials, and correlated endpoints from an illness-death model.
Other data-generating models, such as Weibull or cure-model survival,
copula-based dependence, or resampling from historical data, are easy to
write in a few lines of R. The analysis side of FastSurvival does not
depend on how the data were generated:
[`cutoff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md),
[`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md),
[`switch_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md),
and
[`simsummary_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md)
read a data frame with a fixed set of columns. This vignette generates
trials from a Weibull model with a cure fraction in the treatment group,
a model that
[`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
does not provide, and runs the same workflow as for data from
[`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md).

The columns that the functions read are

| Function | Required columns |
|----|----|
| [`cutoff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md) | `sim`, `accrual_time`, and the time and event columns named by `tte.col` and `event.col` (`tte` and `event` by default) |
| [`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md), [`pairwise_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/pairwise_fast.md) | `sim`, `group`, `accrual_time`, `tte`, `event` |
| [`switch_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md) | `sim`, `group`, `accrual_time`, and the latent times `surv_time` and `dropout_time` |

where `tte` is the observed time from accrual, `event` is 1 for an event
and 0 for censoring, and `accrual_time` is the calendar time of
enrollment. The latent columns are needed only by
[`switch_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md),
which changes the survival time after the switch and recomputes the
observed columns.

``` r

library(FastSurvival)
```

## A Weibull model with a cure fraction

The control group has Weibull survival with shape 1.3 and scale 16
months. In the treatment group 15% of the subjects are cured and the
others have Weibull survival with the same shape and scale 18 months, so
the survival curves separate slowly and the treatment curve levels off
at 0.15. Subjects are enrolled uniformly over 12 months and drop out at
a constant rate of 0.01 per month. The generator below writes the
columns listed above for all simulated trials at once.

``` r

gen_cure_weibull <- function(nsim, n_per, accrual, shape, scale, cure,
                             drop_rate, seed) {
  set.seed(seed)
  n_tot <- 2L * n_per
  N     <- nsim * n_tot
  sim   <- rep(seq_len(nsim), each = n_tot)
  group <- rep(rep(1:2, each = n_per), times = nsim)
  accrual_time <- stats::runif(N, 0, accrual)
  surv_time    <- stats::rweibull(N, shape = shape, scale = scale[group])
  surv_time[stats::runif(N) < cure[group]] <- Inf   # cured subjects
  dropout_time <- stats::rexp(N, rate = drop_rate)
  tte <- pmin(surv_time, dropout_time)
  data.frame(sim, group, accrual_time, surv_time, dropout_time, tte,
             event = as.integer(surv_time <= dropout_time),
             calendar_time = accrual_time + tte)
}

shape <- 1.3
scale <- c(16, 18)
cure  <- c(0, 0.15)

dat <- gen_cure_weibull(nsim = 2000, n_per = 200, accrual = 12,
                        shape = shape, scale = scale, cure = cure,
                        drop_rate = 0.01, seed = 2026)
head(dat)
#>   sim group accrual_time surv_time dropout_time        tte event calendar_time
#> 1   1     1    8.3840817  3.939239  300.0990358  3.9392386     1     12.323320
#> 2   1     1    6.6783661 19.459257    0.1160362  0.1160362     0      6.794402
#> 3   1     1    1.6816795 22.747133   42.4435221 22.7471334     1     24.428813
#> 4   1     1    3.4286797 14.702976   19.9134914 14.7029759     1     18.131656
#> 5   1     1    6.6644281 17.846507  166.8662675 17.8465068     1     24.510935
#> 6   1     1    0.3015738  5.973659   24.0272766  5.9736589     1      6.275233
```

A quick check of the generator compares the proportion of latent
survival times beyond a few time points with the survival function of
the model, `cure + (1 - cure) * exp(-(t / scale)^shape)`.

``` r

surv_model <- function(t, g) cure[g] + (1 - cure[g]) * exp(-(t / scale[g])^shape)
t_chk <- c(6, 12, 24, 36)
data.frame(
  t               = t_chk,
  control_sim     = sapply(t_chk, function(t) mean(dat$surv_time[dat$group == 1] > t)),
  control_model   = surv_model(t_chk, 1),
  treatment_sim   = sapply(t_chk, function(t) mean(dat$surv_time[dat$group == 2] > t)),
  treatment_model = surv_model(t_chk, 2)
)
#>    t control_sim control_model treatment_sim treatment_model
#> 1  6   0.7549500    0.75623042     0.8188400       0.8188069
#> 2 12   0.5025250    0.50258723     0.6206900       0.6210314
#> 3 24   0.1838250    0.18377917     0.3485375       0.3486846
#> 4 36   0.0565925    0.05671565     0.2223875       0.2224537
```

## Analysis times

The interim analysis takes place at 150 events or at month 24, whichever
comes first, and the final analysis at 280 events, but not earlier than
9 months after the interim and not later than month 48.
[`cutoff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md)
returns the calendar time of both looks for every simulated trial.

``` r

cut <- cutoff_fast(dat, event.looks = c(150, 280), max.time = c(24, 48),
                   min.gap = c(NA, 9))
head(cut)
#>      look1    look2
#> 1 17.72435 30.83964
#> 2 17.01267 34.59461
#> 3 18.24693 32.53216
#> 4 18.15349 40.04166
#> 5 17.81107 38.90341
#> 6 18.86204 37.84150
colMeans(cut)
#>    look1    look2 
#> 17.06437 34.10453
```

## Log-rank and max-combo tests

Because the treatment effect builds up late and ends in a plateau, the
max-combo test of Fleming-Harrington weights is compared with the
log-rank test. The nominal one-sided levels 0.002 and 0.024 are close to
those of an O’Brien-Fleming-type spending function at these looks and
are used for illustration. With `mc.alpha` set to these levels, the
max-combo p-value is integrated only when the Bonferroni bounds do not
decide the comparison with the level, which leaves the decisions
unchanged.

``` r

alpha_look <- c(0.002, 0.024)
set.seed(1)
res <- analysis_fast(dat, control = 1, cutoff.looks = cut,
                     stat = c("logrank", "maxcombo"), side = 1,
                     mc.alpha = alpha_look)

oc_lr <- simsummary_fast(res, p.col = "logrank.p",  alpha = alpha_look)
oc_mc <- simsummary_fast(res, p.col = "maxcombo.p", alpha = alpha_look)
data.frame(
  test             = c("Log-rank", "Max-combo"),
  reject_interim   = c(oc_lr[oc_lr$look == "1", "prob.stop.efficacy"],
                       oc_mc[oc_mc$look == "1", "prob.stop.efficacy"]),
  power            = c(oc_lr[oc_lr$look == "overall", "cum.reject"],
                       oc_mc[oc_mc$look == "overall", "cum.reject"])
)
#>        test reject_interim  power
#> 1  Log-rank         0.2445 0.9600
#> 2 Max-combo         0.2325 0.9635
mean(res$maxcombo.p.exact, na.rm = TRUE)
#> [1] 0.09275
```

The last value is the share of the max-combo p-values that required the
multivariate normal integral.

## Crossover after the interim analysis

[`switch_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md)
works on the same data. Here every control subject who is still on study
at the interim analysis switches to the experimental treatment, which
multiplies the remaining survival time by 1.5. The final analysis is
then triggered by the events of the modified data. The outcomes observed
by the interim analysis are unchanged, so the interim cutoffs and
statistics are the same as before.

``` r

sw <- switch_fast(dat, group = 1, when = "cutoff",
                  cutoff = cut[, 1, drop = FALSE], aft.factor = 1.5)
cut_sw <- cutoff_fast(sw, event.looks = c(150, 280), max.time = c(24, 48),
                      min.gap = c(NA, 9))
all.equal(cut_sw[, 1], cut[, 1])
#> [1] TRUE

res_sw <- analysis_fast(sw, control = 1, cutoff.looks = cut_sw,
                        stat = "logrank", side = 1)
all.equal(res_sw$logrank.z[res_sw$look == 1], res$logrank.z[res$look == 1])
#> [1] TRUE

oc_sw <- simsummary_fast(res_sw, p.col = "logrank.p", alpha = alpha_look)
c(log_rank_power_without_switching = oc_lr[oc_lr$look == "overall", "cum.reject"],
  log_rank_power_with_crossover    = oc_sw[oc_sw$look == "overall", "cum.reject"])
#> log_rank_power_without_switching    log_rank_power_with_crossover 
#>                           0.9600                           0.6145
```

## Remarks

Any generator can be used in the same way, provided that it writes one
row per subject with the columns above and that the simulation
identifiers group the rows of each trial. The data do not need to be
sorted, but generating them grouped by `sim`, as
[`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
does, avoids a sorting step in the analysis functions. For an
illness-death model,
[`switch_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md)
reads the columns `e1_surv_time`, `e2_surv_time`, `dropout_time`, and
`intermediate` instead, as described in its help page.
