# Compare Milestone Survival Probabilities Between Two Groups

Compares the Kaplan-Meier survival probabilities of two groups at a
prespecified milestone timepoint. The point estimate of interest is the
difference in milestone survival, treatment minus control. Three
inference methods are provided. The `"wald"` method uses the unpooled
Greenwood variance directly. The `"loglog"` and `"mover"` methods build
the confidence interval for the difference with the method of variance
estimates recovery (MOVER), recovering the variance from the one-sample
complementary log-log and log transformed confidence intervals
respectively. See Tang (2021) for the MOVER difference interval and Tang
(2022) for the use of milestone survival in trial design.

## Usage

``` r
milestone_fast(
  time,
  status,
  group,
  tau,
  method = c("wald", "loglog", "mover"),
  side = c("two.sided", "upper", "lower"),
  conf.level = 0.95,
  presorted = FALSE
)
```

## Arguments

- time:

  A numeric vector of follow-up times.

- status:

  An integer vector of event indicators, 1 for an event and 0 for a
  censored observation.

- group:

  A vector with exactly two distinct values identifying the group. The
  smaller value (or the first factor level) is treated as the control
  group and the other as the treatment group.

- tau:

  The milestone timepoint at which the survival probabilities are
  compared. A single positive number.

- method:

  The inference method for the difference in milestone survival, one of
  `"wald"`, `"loglog"`, or `"mover"`.

- side:

  The alternative hypothesis for the difference, one of `"two.sided"`,
  `"upper"` (treatment survival larger), or `"lower"` (treatment
  survival smaller). The confidence interval is always reported as a
  two-sided interval at `conf.level`.

- conf.level:

  The confidence level for the reported intervals.

- presorted:

  Logical. If `TRUE` the input is assumed to be sorted by `time` in
  ascending order and the internal sort is skipped. This is intended for
  repeated calls inside simulation loops.

## Value

An object of class `"milestone_fast"`, a list with the per-group
milestone survival estimates and standard errors, the difference
estimate with its confidence interval, the test statistic, and the
p-value.

## References

Tang Y (2021). Some new confidence intervals for Kaplan-Meier based
estimators from one and two sample survival data. Statistics in
Medicine, 40(23), 4961-4976.

Tang Y (2022). Complex survival trial design by the product integration
method. Statistics in Medicine, 41(4), 798-814.

## Examples

``` r
set.seed(1)
time <- c(rexp(50, 0.1), rexp(50, 0.07))
status <- rep(1, 100)
group <- rep(c(0, 1), each = 50)
milestone_fast(time, status, group, tau = 10, method = "loglog")
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
```
