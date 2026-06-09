# Fast Log-Rank Test for Two-Group Survival Data

Computes the log-rank test statistic for comparing survival curves
between two groups. Returns either a one-sided Z-score or a two-sided
chi-square statistic. The C++ backend uses a two-pointer merge scan over
pooled sorted vectors, eliminating the
[`rank()`](https://rdrr.io/r/base/rank.html) call that dominates the
pure-R implementation. When a `strata` argument is supplied, the
stratified log-rank test is computed instead, matching
[`survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html) with a
[`strata()`](https://rdrr.io/pkg/survival/man/strata.html) term. A
non-default `weight` argument selects a weighted log-rank test
(Fleming-Harrington, modestly-weighted, Gehan-Breslow, or Tarone-Ware)
for non-proportional hazards.

## Usage

``` r
survdiff_fast(
  time,
  event,
  group,
  control,
  side = 2,
  presorted = FALSE,
  strata = NULL,
  weight = c("logrank", "fh", "mwlrt", "gehan", "tarone-ware"),
  rho = 0,
  gamma = 0,
  t_star = NULL
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

  An integer, either 1 or 2 (default 2). If `side = 1`, returns the
  standardized log-rank statistic (Z-score), defined as
  `(O_1 - E_1) / sqrt(V_1)` for the treatment group, so the Z-score is
  negative when the treatment group has fewer events than expected (a
  protective treatment effect). If `side = 2`, returns the chi-square
  statistic (Z^2).

- presorted:

  A logical value. If `TRUE`, `time`, `event`, and `group` (and `strata`
  when supplied) are assumed to be already sorted in the required order,
  and the internal [`order()`](https://rdrr.io/r/base/order.html) call
  is skipped. If `FALSE` (default), sorting is handled internally. See
  Details for the required order in the stratified case.

- strata:

  An optional vector of stratum labels aligned with `time`. If `NULL`
  (default), the ordinary log-rank test is computed and the behavior is
  identical to earlier versions of this function. If supplied, the
  stratified log-rank test is computed, matching
  [`survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html) with a
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) term. Any
  type that supports equality comparison is accepted. May be combined
  with a non-default `weight` to obtain a stratified weighted log-rank
  test.

- weight:

  A character string naming the weight scheme. `"logrank"` (default) is
  the ordinary unweighted log-rank test and reproduces the behavior of
  earlier versions of this function exactly. `"fh"` is the
  Fleming-Harrington G(rho, gamma) test with weight
  `S(t-)^rho (1 - S(t-))^gamma`. `"mwlrt"` is the modestly-weighted
  log-rank test of Magirr and Burman with weight
  `1 / max(S(t-), S(t_star))`. `"gehan"` is the Gehan-Breslow test with
  weight equal to the at-risk count, and `"tarone-ware"` uses the square
  root of the at-risk count. Here `S(t-)` is the left-continuous pooled
  Kaplan-Meier estimate just prior to each event time.

- rho:

  A numeric Fleming-Harrington first parameter, used only when
  `weight = "fh"`. Defaults to 0.

- gamma:

  A numeric Fleming-Harrington second parameter, used only when
  `weight = "fh"`. Defaults to 0. The pair `rho = 0, gamma = 0`
  reproduces the ordinary log-rank test, and `rho = 0, gamma = 1` is the
  Fleming-Harrington G(0, 1) test for delayed effects.

- t_star:

  A single non-negative numeric value, the timepoint of the
  modestly-weighted log-rank test. Required only when
  `weight = "mwlrt"`. The weight is capped at `1 / S(t_star)`, where
  `S(t_star)` is the smallest pooled Kaplan-Meier value at or after
  `t_star`. A value of 0 yields the ordinary log-rank test.

## Value

An object of class `"survdiff_fast"`, which is a length-one numeric
value with attributes `O0`, `E0`, `O1`, `E1`, `V1`, `side`, and `n`,
plus `strata` (the number of strata) when `strata` is supplied, or
`weight` (the scheme name) when a non-default `weight` is used. The
numeric value is the Z-score when `side = 1`, or the chi-square
statistic when `side = 2`. For a weighted test the value is U / sqrt(V)
(or its square), `V1` holds the weighted variance V, `O0` and `O1` hold
the raw observed event counts, and `E0` and `E1` are `NA` because a
single unweighted expected count is not defined for a weighted test.
Returns `NA_real_` (still with class `"survdiff_fast"`) when the
variance is zero (e.g., all events in one group).

## Details

The log-rank statistic is computed as:

Z = (O_1 - E_1) / sqrt(V_1)

where O_1 is the observed number of events in the treatment group, E_1
is the expected number under the null hypothesis of equal survival, and
V_1 is the hypergeometric variance. Tied event times are handled
correctly: all subjects sharing the same event time form a tied block,
and the block is processed atomically in the two-pointer merge.

When `strata` is `NULL` (default), the ordinary two-group log-rank test
is computed by the C++ core `logrank_core`, which walks the pooled
sorted data with a single two-pointer scan, maintaining running at-risk
counts per group. No rank vector is constructed, so the dominant O(n log
n) cost of [`rank()`](https://rdrr.io/r/base/rank.html) in the pure-R
version is removed.

When `strata` is supplied, the stratified log-rank test is computed by
the C++ core `stratified_logrank_core`. The contributions O_1, E_1, and
V_1 are accumulated within each stratum and then summed across strata,
so the overall statistic is Z = (sum O_1 - sum E_1) / sqrt(sum V_1). A
stratum that contains only one group contributes zero to all three
totals, the same convention used by
[`survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html). This is
the standard stratified log-rank test, equivalent to a log-rank test
that conditions on the stratum at each event time.

When a non-default `weight` is requested, the weighted log-rank test is
computed by the C++ core `weighted_logrank_core`. The statistic is Z = U
/ sqrt(V) with U = sum w_j (O_1j - E_1j) and V = sum w_j^2 n_0j n_1j O_j
(n_j - O_j) / (n_j^2 (n_j - 1)), where the weight w_j is one of the
schemes named by `weight`. The Fleming-Harrington and modestly-weighted
schemes use the left-continuous pooled Kaplan-Meier estimate S(t-),
initialized at 1 and updated after each event time, so the first event
always has S(t-) = 1. The modestly-weighted scheme determines its weight
cap in a first pass over the event times before accumulating U and V in
a second pass; the other schemes accumulate in a single pass. When both
`weight` and `strata` are supplied, the stratified weighted log-rank
test is computed by the C++ core `stratified_weighted_logrank_core`:
each stratum is an independent weighted log-rank test whose weights come
from that stratum's own pooled Kaplan-Meier estimate, and the
per-stratum U and V are summed before standardizing once as Z = sum U /
sqrt(sum V).

When `presorted = TRUE`, the input vectors are assumed to be sorted and
the internal [`order()`](https://rdrr.io/r/base/order.html) call is
skipped. In the unstratified case the assumed order is ascending `time`.
In the stratified case the assumed order is by stratum first and by
ascending `time` within each stratum, so that rows of the same stratum
are contiguous. When `presorted = FALSE` (default), sorting is handled
internally. In simulation loops where the data are generated in the
required order, setting `presorted = TRUE` avoids one O(n log n) pass.

The returned object has class `"survdiff_fast"` and is a length-one
numeric value (Z-score or chi-square) with the underlying counts O_0,
E_0, O_1, E_1, V_1, the requested `side`, and the total sample size
stored as attributes. When `strata` is supplied, the number of strata is
also stored as the attribute `strata`. A
[`print()`](https://rdrr.io/r/base/print.html) method formats the result
similarly to `print(survival::survdiff(...))`, displaying observed and
expected event counts for both the control and treatment groups.

## References

Gehan, E. A. (1965). A generalized Wilcoxon test for comparing
arbitrarily single-censored samples. *Biometrika*, *52*, 203-223.

Mantel, N. (1966). Evaluation of survival data and two new rank order
statistics arising in its consideration. *Cancer Chemotherapy Reports*,
*50*(3), 163-170.

Tarone, R. E., & Ware, J. (1977). On distribution-free tests for
equality of survival distributions. *Biometrika*, *64*, 156-160.

Fleming, T. R., & Harrington, D. P. (1991). *Counting Processes and
Survival Analysis*. New York: John Wiley & Sons.

Magirr, D., & Burman, C.-F. (2019). Modestly weighted logrank tests.
*Statistics in Medicine*, *38*(20), 3782-3790.

## See also

[`survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html) for the
standard implementation.
[`print.survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.survdiff_fast.md)
for the print method.

## Examples

``` r
library(survival)

# Two-sided test: compare with survdiff
fit <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 2)
fit
#> Log-rank test (two-group)
#> 
#>   N = 26,  control = 2
#> 
#>           Observed Expected (O-E)^2/E (O-E)^2/V
#> control          5   6.7665    0.4612    1.0627
#> treatment        7   5.2335    0.5962    1.0627
#> 
#>  Chi-square = 1.063 on 1 df,  p-value = 0.3026  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

chisq_ref <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)$chisq
cat("survdiff_fast chi-square:", as.numeric(fit), "\n")
#> survdiff_fast chi-square: 1.06274 
cat("survdiff      chi-square:", chisq_ref,       "\n")
#> survdiff      chi-square: 1.06274 

# One-sided test (Z-score)
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1)
#> Log-rank test (two-group)
#> 
#>   N = 26,  control = 2
#> 
#>           Observed Expected (O-E)^2/E (O-E)^2/V
#> control          5   6.7665    0.4612    1.0627
#> treatment        7   5.2335    0.5962    1.0627
#> 
#>  Z = 1.031,  one-sided p-value = 0.8487  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# presorted = TRUE: sort once outside, reuse inside a loop
ord <- order(ovarian$futime)
survdiff_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
              control = 2, side = 2, presorted = TRUE)
#> Log-rank test (two-group)
#> 
#>   N = 26,  control = 2
#> 
#>           Observed Expected (O-E)^2/E (O-E)^2/V
#> control          5   6.7665    0.4612    1.0627
#> treatment        7   5.2335    0.5962    1.0627
#> 
#>  Chi-square = 1.063 on 1 df,  p-value = 0.3026  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Stratified log-rank test: compare with survdiff + strata()
fit_str <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2,
                         side = 2, strata = ovarian$resid.ds)
chisq_str <- survdiff(Surv(futime, fustat) ~ rx + strata(resid.ds),
                      data = ovarian)$chisq
cat("stratified survdiff_fast:", as.numeric(fit_str), "\n")
#> stratified survdiff_fast: 1.279643 
cat("stratified survdiff     :", chisq_str,           "\n")
#> stratified survdiff     : 1.279643 

# Weighted log-rank tests for non-proportional hazards
# Fleming-Harrington G(0, 1), emphasizing late differences
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1,
              weight = "fh", rho = 0, gamma = 1)
#> Fleming-Harrington weighted log-rank test (two-group)
#> 
#>   N = 26,  control = 2
#> 
#>           Observed
#> control          5
#> treatment        7
#> 
#>  Z = -0.0101,  one-sided p-value = 0.496  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Modestly-weighted log-rank test with t_star = 365
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1,
              weight = "mwlrt", t_star = 365)
#> Modestly-weighted log-rank test (two-group)
#> 
#>   N = 26,  control = 2
#> 
#>           Observed
#> control          5
#> treatment        7
#> 
#>  Z = 0.7583,  one-sided p-value = 0.7759  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Stratified weighted log-rank test: Fleming-Harrington G(0,1) within strata
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1,
              weight = "fh", rho = 0, gamma = 1, strata = ovarian$resid.ds)
#> Stratified Fleming-Harrington weighted log-rank test (two-group, 2 strata)
#> 
#>   N = 26,  control = 2
#> 
#>           Observed
#> control          5
#> treatment        7
#> 
#>  Z = 0.5667,  one-sided p-value = 0.7145  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# \donttest{
library(microbenchmark)
microbenchmark(
  survdiff_fast = survdiff_fast(ovarian$futime, ovarian$fustat,
                                ovarian$rx, 2, side = 2),
  survdiff      = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
  times = 1000
)
#> Unit: microseconds
#>           expr   min      lq      mean  median      uq     max neval cld
#>  survdiff_fast  23.5   53.25   84.0894   73.35   97.50   729.6  1000  a 
#>       survdiff 694.2 1262.15 1725.7884 1589.85 1814.05 16551.0  1000   b
# }
```
