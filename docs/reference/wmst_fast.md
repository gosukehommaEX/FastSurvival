# Fast window mean survival time and between-group difference

Estimates the window mean survival time (WMST) for a single group, or
the difference in WMST between a treatment group and a control group,
together with standard errors, confidence intervals and a Wald test for
the difference. The WMST is the area under the Kaplan-Meier curve
between the window limits `tau1` and `tau2`. With `tau1 = 0` it is the
restricted mean survival time at `tau2`, and the result then matches
[`survRM2::rmst2`](https://rdrr.io/pkg/survRM2/man/rmst2.html). The
single scan over the sorted data is performed in C++ for use inside
simulation loops.

## Usage

``` r
wmst_fast(
  time,
  event,
  group = NULL,
  control = NULL,
  tau1 = 0,
  tau2 = NULL,
  side = 2,
  conf.level = 0.95,
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
  two-group comparison. If omitted, a single-group WMST is returned.

- control:

  The level of `group` that denotes the control group. Required when
  `group` is supplied.

- tau1:

  Lower limit of the window. Defaults to 0.

- tau2:

  Upper limit of the window. If omitted, the largest time common to both
  groups is used (the minimum over groups of the maximum observed time),
  or the maximum observed time for a single group.

- side:

  Either 2 for a two-sided test or 1 for a one-sided test of treatment
  superiority (difference greater than 0).

- conf.level:

  Confidence level for the intervals.

- presorted:

  Logical; set to `TRUE` when `time`, `event` and `group` are already
  sorted in ascending order of `time`, to skip the internal sort.
  Defaults to `FALSE`.

## Value

A named numeric vector of class `"wmst_fast"`. For a single group the
elements are the WMST, its standard error and confidence limits. For two
groups the elements are the control and treatment WMST values, their
difference, the standard errors, the confidence limits and the Wald
statistics for the difference.

## Details

The variance of the WMST estimate uses the counting process (Greenwood)
form in which each event time at or before `tau2` contributes the
squared remaining window area. Event times at or before `tau1`
contribute the full window area, since the hazard before `tau1` scales
the survival curve throughout the window.

The difference is computed as treatment minus control, so a positive
difference indicates a longer mean survival time within the window under
treatment.

## References

Paukner, M., & Chappell, R. (2021). Window mean survival time.
Statistics in Medicine, 40(25), 5521-5533.

## Examples

``` r
set.seed(1)
n <- 200
g <- rep(0:1, each = n / 2)
tt <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.07))
cc <- rexp(n, rate = 0.02)
time <- pmin(tt, cc)
event <- as.integer(tt <= cc)
wmst_fast(time, event, group = g, control = 0, tau1 = 2, tau2 = 12)
#> Window mean survival time (two-group)
#> 
#>   window = [2, 12],  control = 0
#>   alternative = two.sided
#> 
#>             wmst std.err lower 95% upper 95%
#> control   5.8371   0.402    5.0491    6.6251
#> treatment 6.6057   0.402    5.8177    7.3936
#> 
#>                                     Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  0.7686   -0.3458    1.8829 1.352    0.176
```
