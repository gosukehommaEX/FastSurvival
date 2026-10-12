# Compare Milestone Survival Probabilities Between Two Groups

Compares the Kaplan-Meier survival probabilities of two groups at a
prespecified milestone timepoint. The point estimate of interest is the
difference in milestone survival, treatment minus control. Three
inference methods are provided. The `"wald"` method uses the unpooled
Greenwood variance directly. The `"loglog"` and `"mover"` methods build
the confidence interval for the difference with the method of variance
estimates recovery (MOVER), recovering the variance from the one-sample
complementary log-log and log transformed confidence intervals
respectively. When the estimate of a group is 0 or 1 at `tau`, its
one-sample interval degenerates to the estimate, so the interval of the
difference is still given; the `"loglog"` test statistic is then `NA`.
See Tang (2021) for the MOVER difference interval and Tang (2022) for
the use of milestone survival in trial design.

## Usage

``` r
milestone_fast(
  time,
  event,
  group,
  control,
  side = 2,
  conf.level = 0.95,
  tau,
  method = c("wald", "loglog", "mover"),
  presorted = FALSE,
  data = NULL
)
```

## Arguments

- time:

  A numeric vector of follow-up times. Alternatively, a formula
  `Surv(time, event) ~ group`, with the variables in `data`; see
  [`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).

- event:

  An integer vector of event indicators, 1 for an event and 0 for a
  censored observation.

- group:

  A vector with exactly two distinct values identifying the group.

- control:

  The value of `group` that denotes the control group. The other value
  is the treatment group and the difference is reported as treatment
  minus control. With the formula interface it defaults to the first
  level of the group variable.

- side:

  1 for a one-sided test in the direction of treatment benefit
  (treatment milestone survival larger than control) or 2 for a
  two-sided test (default 2). The confidence interval is always reported
  as a two-sided interval at `conf.level`.

- conf.level:

  The confidence level for the reported intervals.

- tau:

  The milestone timepoint at which the survival probabilities are
  compared. A single positive number. A warning is given when it exceeds
  the largest observed time of either group, because the Kaplan-Meier
  curve is not estimated beyond that time.

- method:

  The inference method for the difference in milestone survival, one of
  `"wald"`, `"loglog"`, or `"mover"`.

- presorted:

  Logical. If `TRUE` the input is assumed to be sorted by `time` in
  ascending order and the internal sort is skipped. This is intended for
  repeated calls inside simulation loops. The order is checked, and an
  error is given when it does not hold.

- data:

  A data frame in which the variables of a formula given as `time` are
  evaluated. Used only with the formula interface; see
  [`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).

## Value

An object of class `"milestone_fast"`, a list with the per-group
milestone survival estimates and standard errors, the difference
estimate with its confidence interval, the test statistic, and the
p-value. The object also has the class `"fastsurv"`, whose
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods are
described in
[`as.data.frame.fastsurv`](https://gosukehommaEX.github.io/FastSurvival/reference/as.data.frame.fastsurv.md).

## Details

The `"wald"` statistic is the difference divided by its Greenwood
standard error. The `"loglog"` statistic is the difference of the
complementary log-log transforms of the two estimates divided by its
standard error, so it is negative when the treatment group has the
higher milestone survival, and its one-sided p-value is the lower tail.
The `"mover"` statistic is the upper `(1 - conf.level) / 2` quantile of
the standard normal distribution (1.96 for `conf.level = 0.95`) times
the difference divided by the half-width of the MOVER interval on the
side of zero. It is obtained by inverting the interval at `conf.level`,
so the statistic and its p-value depend on `conf.level`, and the
two-sided p-value is below `1 - conf.level` exactly when the interval
excludes zero. For `"wald"` and `"mover"` a positive statistic favors
treatment and the one-sided p-value is the upper tail.

## References

Tang, Y. (2021). Some new confidence intervals for Kaplan-Meier based
estimators from one and two sample survival data. *Statistics in
Medicine*, *40*(23), 4961-4976.

Tang, Y. (2022). Complex survival trial design by the product
integration method. *Statistics in Medicine*, *41*(4), 798-814.

## Examples

``` r
set.seed(1)
time <- c(rexp(50, 0.1), rexp(50, 0.07))
event <- rep(1, 100)
group <- rep(c(0, 1), each = 50)
milestone_fast(time, event, group, control = 0, tau = 10, method = "loglog")
#> Milestone survival (two-group)
#> 
#>   tau = 10,  control = 0
#>   method = loglog,  alternative = two.sided
#> 
#>           survival std.err lower 95% upper 95%
#> control       0.42  0.0698    0.2829    0.5510
#> treatment     0.56  0.0702    0.4124    0.6842
#> 
#>                                      Est. lower 95% upper 95%      z Pr(>|z|)
#> difference (treatment - control)  0.14000  -0.05738   0.32500 -1.395    0.163

# Formula interface
d <- data.frame(time, event, group)
milestone_fast(Surv(time, event) ~ group, data = d, tau = 10)
#> Milestone survival (two-group)
#> 
#>   tau = 10,  control = 0
#>   method = wald,  alternative = two.sided
#> 
#>           survival std.err lower 95% upper 95%
#> control       0.42  0.0698    0.2832    0.5568
#> treatment     0.56  0.0702    0.4224    0.6976
#> 
#>                                      Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  0.14000  -0.05403   0.33403 1.414    0.157
```
