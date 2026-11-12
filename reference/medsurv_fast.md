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
  bw = NULL,
  presorted = FALSE,
  data = NULL
)
```

## Arguments

- time:

  Numeric vector of event or censoring times. Alternatively, a formula
  `Surv(time, event) ~ group` or `Surv(time, event) ~ 1`, with the
  variables in `data`; see
  [`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).

- event:

  Integer vector, 1 for an event and 0 for censoring.

- group:

  Optional grouping vector with exactly two distinct levels for a
  two-group comparison. If omitted, a single-group median is returned.

- control:

  The level of `group` that denotes the control group. Required when
  `group` is supplied. With the formula interface it defaults to the
  first level of the group variable.

- side:

  Either 2 for a two-sided test or 1 for a one-sided test of treatment
  superiority (difference greater than 0).

- conf.level:

  Confidence level for the intervals.

- conf.type:

  Confidence interval type for each group median, either `"plain"` or
  `"log"`. The interval for the difference is always on the plain scale.

- bw:

  Optional kernel bandwidth for the hazard at the median. Either a
  single value applied to both groups or one value per group. If
  omitted, the default described in Details is computed for each group.

- presorted:

  Logical; set to `TRUE` when `time`, `event` and `group` are already
  sorted in ascending order of `time`, to skip the internal sort.
  Defaults to `FALSE`. The order is checked, and an error is given when
  it does not hold.

- data:

  A data frame in which the variables of a formula given as `time` are
  evaluated. Used only with the formula interface; see
  [`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).
  With the native pipe, `d |> medsurv_fast(formula, ...)` uses the data
  frame `d` in the same way.

## Value

A named numeric vector of class `"medsurv_fast"`. For a single group the
elements are the median, its standard error and confidence limits. For
two groups the elements are the control and treatment medians, their
difference, the standard errors, the confidence limits and the Wald
statistics for the difference. The object also has the class
`"fastsurv"`, whose
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods are
described in
[`as.data.frame.fastsurv`](https://gosukehommaEX.github.io/FastSurvival/reference/as.data.frame.fastsurv.md).

## Details

The median in each group is the first event time at which the
Kaplan-Meier estimate drops to 0.5 or below. When the estimate equals
0.5 exactly over an interval, the midpoint of that interval is returned
(the next event time, or the last observed time when no later event
exists, closes the interval). This is the convention used by `survfit`,
including its tolerance of `sqrt(.Machine$double.eps)` when comparing
the estimate with 0.5. Collett (2023, Section 2.4) defines the median in
the same way. When the estimate stays at 0.5 up to an infinite observed
time (for example a subject without a finite event or dropout time in
simulated data), the midpoint is not defined, and the median of that
group and the statistics that depend on it are `NA`.

The variance of the estimated median follows from the delta method for a
percentile (Collett, 2023, Section 2.5): var(median) = var(S(median)) /
f^2, where f is the density of the survival time at the median. The
variance of S(median) is Greenwood's, S(median)^2 times greenwood_sum,
and the density is estimated as S(median) times the hazard at the
median, so that var(median) = greenwood_sum / hazard^2. The hazard at
the median is the kernel smoother of Ramlau-Hansen (1983) applied to the
Nelson-Aalen increments, with the Epanechnikov kernel 0.75 (1 - u^2) on
\[-1, 1\] and bandwidth `bw`. Collett (2023) estimates the density
instead by a difference quotient of the Kaplan-Meier estimate.

The default bandwidth of a group is (40 sqrt(pi))^(1/5) s m^(-1/5),
about 2.34 s m^(-1/5), where m is the number of events in the group and
s the standard deviation of their times. This is the bandwidth that
minimizes the asymptotic mean integrated squared error of an
Epanechnikov kernel estimate of a normal density, applied to the event
times as a simple default. When a group has fewer than two events, or
their times do not vary, the range of the event times divided by m^(1/3)
is used, or 1 when that range is 0.

The difference is computed as treatment minus control, so a positive
difference indicates a longer median survival time under treatment.

## References

Collett, D. (2023). *Modelling Survival Data in Medical Research* (4th
ed.). Chapman and Hall/CRC.

Ramlau-Hansen, H. (1983). Smoothing counting process intensities by
means of kernel functions. *The Annals of Statistics*, *11*(2), 453-466.

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
#>   conf.type = log,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    8.1851  0.9241    6.5603   10.2124
#> treatment 10.7251  1.7252    7.8249   14.7002
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  2.540    -1.296     6.376 1.298    0.194
medsurv_fast(time, event, group = g, control = 0, bw = 10)
#> Median survival time (two-group)
#> 
#>   control = 0
#>   conf.type = log,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    8.1851  1.0445    6.3738   10.5112
#> treatment 10.7251  1.7302    7.8177   14.7137
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  2.540    -1.421     6.501 1.257    0.209

# Formula interface
d <- data.frame(time, event, g)
medsurv_fast(Surv(time, event) ~ g, data = d)
#> Median survival time (two-group)
#> 
#>   control = 0
#>   conf.type = log,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    8.1851  0.9241    6.5603   10.2124
#> treatment 10.7251  1.7252    7.8249   14.7002
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  2.540    -1.296     6.376 1.298    0.194
d |> medsurv_fast(Surv(time, event) ~ g)
#> Median survival time (two-group)
#> 
#>   control = 0
#>   conf.type = log,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    8.1851  0.9241    6.5603   10.2124
#> treatment 10.7251  1.7252    7.8249   14.7002
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  2.540    -1.296     6.376 1.298    0.194
```
