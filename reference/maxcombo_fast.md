# Fast Max-Combo Weighted Log-Rank Test for Two-Group Survival Data

Computes the max-combo test, the maximum over a set of
Fleming-Harrington weighted log-rank statistics, for comparing survival
between two groups under non-proportional hazards. The C++ backend
evaluates every weighted numerator and the full between-scheme
covariance matrix in a single scan over the pooled sorted data, and the
p-value is obtained from the multivariate normal distribution implied by
the correlation of the component statistics. The test is robust to the
shape of the hazard difference because the most extreme of several
complementary weights is taken, with the multiplicity accounted for
through the joint distribution.

## Usage

``` r
maxcombo_fast(
  time,
  event,
  group,
  control,
  side = 1,
  rho = c(0, 0, 1, 1),
  gamma = c(0, 1, 0, 1),
  presorted = FALSE,
  abseps = 1e-05,
  maxpts = 25000
)
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

  An integer, either 1 or 2. If `side = 1` (default), the one-sided
  max-combo test for treatment benefit is computed. If `side = 2`, the
  two-sided test based on the maximum absolute component is computed.

- rho:

  A numeric vector of Fleming-Harrington first parameters, one per
  component weight. Defaults to `c(0, 0, 1, 1)`.

- gamma:

  A numeric vector of Fleming-Harrington second parameters, one per
  component weight, aligned with `rho`. Defaults to `c(0, 1, 0, 1)`. The
  default pairs are the standard four-weight max-combo: G(0,0) for
  proportional hazards, G(0,1) for late differences, G(1,0) for early
  differences, and G(1,1) for middle differences.

- presorted:

  A logical value. If `TRUE`, `time`, `event`, and `group` are assumed
  to be sorted in ascending order of `time`, and the internal
  [`order()`](https://rdrr.io/r/base/order.html) call is skipped. If
  `FALSE` (default), sorting is handled internally.

- abseps:

  A single positive numeric value, the absolute error tolerance passed
  to the multivariate normal integration. Defaults to 1e-5. Larger
  values speed up the four-or-more-weight case at the cost of p-value
  precision.

- maxpts:

  A single positive integer, the maximum number of function evaluations
  for the quasi-Monte-Carlo integration used when four or more weights
  are supplied. Defaults to 25000.

## Value

An object of class `"maxcombo_fast"`, a named numeric vector of length
two with elements `statistic` (the max-combo statistic; min_k Z_k when
`side = 1`, so a negative value favours treatment, and max_k abs(Z_k)
when `side = 2`) and `p.value`. The component Z-scores are stored in the
attribute `z`, their correlation matrix in `corr`, the
Fleming-Harrington parameters in `rho` and `gamma`, the requested
`side`, and the total sample size in `n`. Returns `NA` values (still
with class `"maxcombo_fast"`) when any component variance is zero or not
finite.

## Details

Each component is a Fleming-Harrington G(rho, gamma) weighted log-rank
statistic Z_k = U_k / sqrt(V_kk), where U_k = sum w_k (O_1 - E_1) and
the weight is w_k = S(t-)^rho_k (1 - S(t-))^gamma_k evaluated from the
left-continuous pooled Kaplan-Meier estimate S(t-). The between-scheme
covariance is V_ab = sum w_a w_b v, with v the hypergeometric variance
increment shared by all schemes, so the diagonal of V reproduces the
single-scheme weighted variances and the off-diagonal entries give the
correlation matrix R of the component Z-scores. The sign convention
matches
[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md):
a component Z is negative when the treatment group is favoured.

The max-combo statistic and its p-value depend on `side`. When
`side = 1`, the statistic is the most negative component, min_k Z_k, so
that a negative value favours the treatment group in the same way as
[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
with `side = 1`. The one-sided p-value is 1 - P(G_1 \>= m, ..., G_K \>=
m) for G distributed as multivariate normal with mean zero and
correlation R, where m = min_k Z_k. When `side = 2`, the statistic is
max_k abs(Z_k) and the p-value is 1 - P(-m \<= G_1 \<= m, ..., -m \<=
G_K \<= m).

The joint normal probability is evaluated by dimension. With a single
weight the univariate normal is used. With two or three weights the
deterministic
[`TVPACK`](https://rdrr.io/pkg/mvtnorm/man/algorithms.html) algorithm is
used. With four or more weights the quasi-Monte-Carlo
[`GenzBretz`](https://rdrr.io/pkg/mvtnorm/man/algorithms.html) algorithm
is used, whose precision is governed by `abseps` and `maxpts`. In a
simulation study the Monte Carlo error of the estimated rejection rate
is driven by the number of simulated trials rather than by the precision
of each individual p-value, so `abseps` can be loosened to speed up the
four-weight case with negligible effect on the operating
characteristics.

When `presorted = TRUE`, the inputs are assumed to be sorted in
ascending order of `time` and the internal
[`order()`](https://rdrr.io/r/base/order.html) call is skipped, which is
useful inside simulation loops where the data are generated in sorted
order.

## References

Lin, R. S., Lin, J., Roychoudhury, S., et al. (2020). Alternative
analysis methods for time to event endpoints under nonproportional
hazards: a comparative analysis. *Statistics in Biopharmaceutical
Research*, *12*(2), 187-198.

Karrison, T. G. (2016). Versatile tests for comparing survival curves
based on weighted log-rank statistics. *The Stata Journal*, *16*(3),
678-690.

## See also

[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
for the single-scheme weighted log-rank test.

## Examples

``` r
library(survival)

# Standard four-weight max-combo, one-sided
fit <- maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
fit["statistic"]
#> statistic 
#> -1.298019 
fit["p.value"]
#>   p.value 
#> 0.1511426 

# Two-sided test
maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 1, side = 2)
#> Max-combo weighted log-rank test (two-group)
#> 
#>   N = 26
#> 
#>               Z
#> FH(0,0) -1.0309
#> FH(0,1)  0.0101
#> FH(1,0) -1.2980
#> FH(1,1) -0.0576
#> 
#>  Max-combo statistic = 1.298 (two-sided),  p-value = 0.3023

# Custom weight set: proportional plus late-difference only
maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 1,
              rho = c(0, 0), gamma = c(0, 1))
#> Max-combo weighted log-rank test (two-group)
#> 
#>   N = 26
#> 
#>               Z
#> FH(0,0) -1.0309
#> FH(0,1)  0.0101
#> 
#>  Max-combo statistic = -1.031 (one-sided),  p-value = 0.2045

# presorted = TRUE: sort once outside, reuse inside a loop
ord <- order(ovarian$futime)
maxcombo_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
              control = 1, presorted = TRUE)
#> Max-combo weighted log-rank test (two-group)
#> 
#>   N = 26
#> 
#>               Z
#> FH(0,0) -1.0309
#> FH(0,1)  0.0101
#> FH(1,0) -1.2980
#> FH(1,1) -0.0576
#> 
#>  Max-combo statistic = -1.298 (one-sided),  p-value = 0.1512

# \donttest{
# Cross-check against simtrial::maxcombo (one-sided)
if (requireNamespace("simtrial", quietly = TRUE)) {
  df <- data.frame(
    stratum   = "All",
    treatment = ifelse(ovarian$rx == 2, "experimental", "control"),
    tte       = ovarian$futime,
    event     = ovarian$fustat
  )
  simtrial::maxcombo(df, rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1),
                     return_corr = TRUE)$p_value
}
#> [1] 0.1512125
# }
```
