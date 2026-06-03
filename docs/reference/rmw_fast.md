# Robust Modestly-Weighted Log-Rank Test for Two-Group Survival Data

Computes the robust modestly-weighted (rMW) log-rank test of Magirr and
Ohrn, which combines the standard log-rank test with a single
modestly-weighted log-rank test. The test statistic is the maximum of
the two standardized components, evaluated against their joint null
distribution. Because the standard log-rank statistic is included as one
of the two components and the components are strongly correlated under
the null, the multiplicity adjustment is small, so the test loses little
power relative to the log-rank test in worst-case scenarios while
gaining substantial power under delayed effects. The two component
numerators, their variances, and their null covariance are computed in a
single pass by the C++ core `rmw_core`.

## Usage

``` r
rmw_fast(time, event, group, control, side, s_star = 0.5, presorted = FALSE)
```

## Arguments

- time:

  A numeric vector of follow-up times for all subjects.

- event:

  An integer or numeric vector of event indicators (1 = event, 0 =
  censored), aligned with `time`.

- group:

  A vector of group labels aligned with `time`.

- control:

  A scalar value indicating which level of `group` represents the
  control group.

- side:

  An integer, either 1 or 2. If `side = 1`, the statistic is the minimum
  of the two standardized components and a one-sided p-value is
  returned. If `side = 2`, the statistic is the maximum of their
  absolute values and a two-sided p-value is returned.

- s_star:

  A single numeric value in `(0, 1]`, the survival-probability threshold
  of the modestly-weighted component. The weight is capped at
  `1 / s_star`. Defaults to 0.5. A value of 1 caps the weight at one, so
  the modestly-weighted component equals the log-rank component and the
  test reduces to the ordinary log-rank test.

- presorted:

  A logical value. If `TRUE`, `time`, `event`, and `group` are assumed
  to be already sorted by ascending `time`, and the internal
  [`order()`](https://rdrr.io/r/base/order.html) call is skipped. If
  `FALSE` (default), sorting is handled internally.

## Value

An object of class `"rmw_fast"`, a length-two numeric vector
`c(statistic, p.value)` with attributes `z` (the named component
Z-scores `c(logrank, mwlrt)`), `corr` (the 2 by 2 null correlation
matrix of the two components), `s_star`, `O1` (the observed number of
events in the treatment group), `side`, and `n` (the total sample size).
Returns `NA` statistic and p-value (still with class `"rmw_fast"`) when
either component variance is zero (e.g., all events in one group).

## Details

The first component is the ordinary log-rank statistic, with weight one
at every event time. The second component is a modestly-weighted
log-rank statistic with weight `min(1 / S(t-), 1 / s_star)`, where
`S(t-)` is the left-continuous pooled Kaplan-Meier estimate just prior
to each event time and `s_star` is a survival-probability threshold.
This is the survival-threshold parameterization of the modestly-weighted
test of Magirr and Burman, in which the weight is capped at
`1 / s_star`. It differs from the timepoint parameterization exposed by
[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
with `weight = "mwlrt"`, where the cap is derived from a timepoint
`t_star`. The choice `s_star = 0.5` caps the weight at 2 and is the
value used in the original rMW article.

Writing the two standardized components as `Z_lr` and `Z_mw`, under the
null hypothesis of equal survival the pair is asymptotically bivariate
normal with zero means, unit variances, and correlation
`rho = C / sqrt(V_lr V_mw)`, where `C` is the covariance of the two
numerators returned by the C++ core. When `side = 1` the statistic is
`min(Z_lr, Z_mw)`, so a protective treatment effect (fewer events than
expected) yields a small, negative value, and the one-sided p-value is
`P(min(Z_lr, Z_mw) <= observed)` under the joint null. When `side = 2`
the statistic is `max(abs(Z_lr), abs(Z_mw))` and the two-sided p-value
is `P(max(abs(Z_lr), abs(Z_mw)) >= observed)`. The joint normal
probability is evaluated with
[`mvtnorm::pmvnorm`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html),
using the exact `TVPACK` algorithm for the one-sided half-space and the
deterministic `Miwa` algorithm for the two-sided rectangle.

When `presorted = TRUE`, the input vectors are assumed to be sorted by
ascending `time` and the internal
[`order()`](https://rdrr.io/r/base/order.html) call is skipped, which
avoids one O(n log n) pass in simulation loops where the data are
already generated in time order.

## References

Magirr, D., & Ohrn, F. (2026). Robust modestly weighted log-rank tests.
*Pharmaceutical Statistics*, *25*(1), e70066.

Magirr, D., & Burman, C.-F. (2019). Modestly weighted logrank tests.
*Statistics in Medicine*, *38*(20), 3782-3790.

## See also

[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
for the standard and weighted log-rank tests.
[`print.rmw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.rmw_fast.md)
for the print method.

## Examples

``` r
library(survival)

# One-sided robust modestly-weighted test with s_star = 0.5
fit <- rmw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                control = 1, side = 1, s_star = 0.5)
fit
#> Robust modestly-weighted log-rank test (two-group)
#> 
#>   N = 26,  s_star = 0.5
#> 
#>                         Z
#> log-rank          -1.0309
#> modestly-weighted -0.7583
#> 
#>   Null correlation = 0.9821
#>  min Z = -1.031,  one-sided p-value = 0.169

# The log-rank component matches survdiff_fast with weight = "logrank"
z_lr <- as.numeric(survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                                 control = 1, side = 1))
cat("rmw_fast log-rank component:", attr(fit, "z")[["logrank"]], "\n")
#> rmw_fast log-rank component: -1.030893 
cat("survdiff_fast log-rank Z   :", z_lr, "\n")
#> survdiff_fast log-rank Z   : -1.030893 

# Two-sided test
rmw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
         control = 1, side = 2, s_star = 0.5)
#> Robust modestly-weighted log-rank test (two-group)
#> 
#>   N = 26,  s_star = 0.5
#> 
#>                         Z
#> log-rank          -1.0309
#> modestly-weighted -0.7583
#> 
#>   Null correlation = 0.9821
#>  Max |Z| = 1.031,  two-sided p-value = 0.338

# presorted = TRUE: sort once outside, reuse inside a loop
ord <- order(ovarian$futime)
rmw_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
         control = 1, side = 1, s_star = 0.5, presorted = TRUE)
#> Robust modestly-weighted log-rank test (two-group)
#> 
#>   N = 26,  s_star = 0.5
#> 
#>                         Z
#> log-rank          -1.0309
#> modestly-weighted -0.7583
#> 
#>   Null correlation = 0.9821
#>  min Z = -1.031,  one-sided p-value = 0.169
```
