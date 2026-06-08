# Fast non-parametric median survival time and between-group difference

Estimates the Kaplan-Meier median survival time for a single group, or
the difference in median survival time between a treatment group and a
control group, together with standard errors, confidence intervals and a
Wald test for the difference. The estimator is designed for repeated
evaluation inside simulation loops, with the single scan over the sorted
data performed in C++.

## Usage

``` r
medsurv_fast(
  time,
  event,
  group = NULL,
  control = NULL,
  side = 2,
  conf.level = 0.95,
  conf.type = "log",
  method = c("km", "nph"),
  bw = NULL,
  presorted = FALSE
)
```

## Arguments

- time:

  Numeric vector of event or censoring times.

- event:

  Integer vector, 1 for an event and 0 for censoring.

- group:

  Optional grouping vector with exactly two distinct levels for a
  two-group comparison. If omitted, a single-group median is returned.

- control:

  The level of `group` that denotes the control group. Required when
  `group` is supplied.

- side:

  Either 2 for a two-sided test or 1 for a one-sided test of treatment
  superiority (difference greater than 0).

- conf.level:

  Confidence level for the intervals.

- conf.type:

  Confidence interval type for each group median, either `"plain"` or
  `"log"`. The interval for the difference is always on the plain scale.

- method:

  Variance method, either `"km"` (Greenwood increment with a kernel
  hazard) or `"nph"` (local constant hazard with a `sum 1 / (Y - k)^2`
  increment, matching
  [`nph::nphparams`](https://rdrr.io/pkg/nph/man/nphparams.html)).

- bw:

  Optional kernel bandwidth for the hazard at the median, used only when
  `method = "km"`. Either a single value applied to both groups or one
  value per group. If omitted, a Silverman type default is used.

- presorted:

  Logical; set to `TRUE` when `time`, `event` and `group` are already
  sorted in ascending order of `time`, to skip the internal sort.
  Defaults to `FALSE`.

## Value

A named numeric vector of class `"medsurv_fast"`. For a single group the
elements are the median, its standard error and confidence limits. For
two groups the elements are the control and treatment medians, their
difference, the standard errors, the confidence limits and the Wald
statistics for the difference.

## Details

The median in each group is the first event time at which the
Kaplan-Meier estimate drops to 0.5 or below, matching the convention
used by `survfit`. The point estimate is the same for both variance
methods.

Two variance methods are available through the `method` argument. With
`method = "km"` the variance of the estimated median follows the
counting process delta method, var(median) = greenwood_sum / hazard^2,
where the instantaneous hazard at the median is obtained by a
Ramlau-Hansen kernel smoother with an Epanechnikov kernel and bandwidth
`bw`. With `method = "nph"` the variance reproduces the computation used
by [`nph::nphparams`](https://rdrr.io/pkg/nph/man/nphparams.html) with
`param_type = "Q"` and `haz_method = "local"`: a local constant hazard
at the median and a `sum 1 / (Y - k)^2` variance increment. The two
methods give the same median but generally different standard errors,
since the variance of a quantile depends on the local hazard estimate.
When the Kaplan-Meier and Nelson-Aalen medians coincide, which is the
usual case, `method = "nph"` reproduces the
[`nph::nphparams`](https://rdrr.io/pkg/nph/man/nphparams.html) standard
error and p-value to numerical precision.

The difference is computed as treatment minus control, so a positive
difference indicates a longer median survival time under treatment.

## References

Fleming, T. R., & Harrington, D. P. (1991). Counting Processes and
Survival Analysis. New York: John Wiley & Sons.

## Examples

``` r
set.seed(1)
n <- 200
g <- rep(0:1, each = n / 2)
tt <- rexp(n, rate = ifelse(g == 0, 0.1, 0.07))
cc <- rexp(n, rate = 0.02)
time <- pmin(tt, cc)
event <- as.integer(tt <= cc)
medsurv_fast(time, event, group = g, control = 0)
#> Median survival time (two-group)
#> 
#>   control = 0
#>   method = km,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    8.1851  0.9708    6.4873   10.3273
#> treatment 10.7251  1.4586    8.2156   14.0012
#> 
#>                                     Est. lower 95% upper 95%    z Pr(>|z|)
#> difference (treatment - control)  2.5400   -0.8942    5.9742 1.45    0.147
medsurv_fast(time, event, group = g, control = 0, method = "nph")
#> Median survival time (two-group)
#> 
#>   control = 0
#>   method = nph,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    8.1851  0.9813    6.4710   10.3533
#> treatment 10.7251  1.6181    7.9796   14.4152
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  2.540    -1.169     6.249 1.342     0.18
```
