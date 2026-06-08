# Fast Restricted Mean Survival Time (Single Group or Two-Group Comparison)

Computes the restricted mean survival time (RMST) up to a horizon `tau`
from the Kaplan-Meier estimator. With a single group (the default), it
returns the RMST with its Greenwood-type standard error and a Wald
confidence interval. When a `group` is supplied, it additionally returns
the two-group contrasts: the RMST difference (treatment minus control)
and the RMST ratio (treatment over control), each with a standard error,
confidence interval, and two-sided test. The C++ backend integrates the
survival step function in a single scan and is reused once per group, so
the function is suitable for simulation loops with `presorted = TRUE`.

## Usage

``` r
rmst_fast(
  time,
  event,
  group = NULL,
  control = NULL,
  side = 2,
  conf.level = 0.95,
  tau,
  presorted = FALSE
)
```

## Arguments

- time:

  A numeric vector of event or censoring times.

- event:

  An integer or numeric vector of event indicators (1 = event, 0 =
  censored), aligned with `time`.

- group:

  An optional vector of group labels aligned with `time`. If `NULL`
  (default), a single-group RMST is computed. If supplied, the two-group
  contrasts are returned and `control` must be given.

- control:

  A scalar value indicating which level of `group` represents the
  control group. Required when `group` is supplied.

- side:

  1 for a one-sided test in the direction of treatment benefit (longer
  restricted mean survival in the treatment group) or 2 for a two-sided
  test (default 2). Applies to the two-group contrasts; the confidence
  intervals are always two-sided at `conf.level`.

- conf.level:

  A single numeric value in (0, 1) specifying the confidence level.
  Defaults to 0.95.

- tau:

  A single positive numeric value specifying the restriction horizon.

- presorted:

  A logical value. If `TRUE`, the inputs are assumed to be sorted in
  ascending order of `time`. If `FALSE` (default), sorting is handled
  internally.

## Value

An object of class `"rmst_fast"`, a named numeric vector. In
single-group mode it has length 4 with elements `rmst`, `std.err`,
`lower`, and `upper`. In two-group mode it contains the per-group RMST
(`rmst.ctrl`, `rmst.trt`), the difference contrast (`diff`, `se.diff`,
`diff.lower`, `diff.upper`, `z.diff`, `p.diff`), and the ratio contrast
(`ratio`, `ratio.lower`, `ratio.upper`, `z.ratio`, `p.ratio`). The
restriction horizon, confidence level, and test side are stored as
attributes `tau`, `conf.level`, and `side`; in two-group mode the
`control` label is also stored. Returns `NA_real_` values (still with
class `"rmst_fast"`) when `n` is zero in single-group mode.

## Details

The RMST is the area under the Kaplan-Meier curve from 0 to `tau`:

RMST(tau) = integral over \[0, tau\] of S(u) du.

The variance follows the Greenwood-type estimator

Var\[RMST\] = sum\_{t_i \<= tau} A_i^2 d_i / (n_i (n_i - d_i)),

where A_i = integral over \[t_i, tau\] of S(u) du is the area to the
right of event time t_i. This matches the restricted mean reported by
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) and by the
`survRM2` package.

When `group` is supplied, the treatment group is the level of `group`
that is not equal to `control`. Writing the treatment and control RMST
as r1 and r0 with variances v1 and v0, the difference contrast is diff =
r1 - r0 with Var\[diff\] = v1 + v0, using the independence of the two
groups. The ratio contrast is formed on the log scale by the delta
method, Var\[log(r1 / r0)\] = v1 / r1^2 + v0 / r0^2, with the confidence
interval exponentiated back to the ratio scale. These match the
unadjusted contrasts reported by `survRM2`.

When `presorted = TRUE`, the input vectors are assumed to be sorted in
ascending order of `time`; splitting into groups preserves the ascending
order within each group, so no re-sorting is performed. When
`presorted = FALSE` (default), sorting is handled internally. In
simulation loops where the data are generated in sorted order,
`presorted = TRUE` avoids one O(n log n) pass.

## References

Royston, P., & Parmar, M. K. B. (2013). Restricted mean survival time:
an alternative to the hazard ratio for the design and analysis of
randomized trials with a time-to-event outcome. *BMC Medical Research
Methodology*, *13*, 152.

Uno, H., Claggett, B., Tian, L., et al. (2014). Moving beyond the hazard
ratio in quantifying the between-group difference in survival analysis.
*Journal of Clinical Oncology*, *32*(22), 2380-2385.

## See also

[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)
for the Kaplan-Meier estimate at a single time point.

## Examples

``` r
set.seed(42)
t_raw <- rexp(100, rate = 1 / 10)
e_raw <- rbinom(100, 1, 0.7)

# Single-group RMST
rmst_fast(t_raw, e_raw, tau = 10)
#> Restricted mean survival time (single-group)
#> 
#>   tau = 10
#> 
#>            RMST std.err lower 95% upper 95%
#> tau = 10 7.0341  0.3452    6.3575    7.7107

# Single-group, pre-sorted (sort once, reuse in a loop)
ord <- order(t_raw)
rmst_fast(t_raw[ord], e_raw[ord], tau = 10, presorted = TRUE)
#> Restricted mean survival time (single-group)
#> 
#>   tau = 10
#> 
#>            RMST std.err lower 95% upper 95%
#> tau = 10 7.0341  0.3452    6.3575    7.7107

# Two-group comparison (difference and ratio)
set.seed(7)
n <- 200
time  <- c(rexp(n, 0.10), rexp(n, 0.07))
event <- rbinom(2 * n, 1, 0.8)
group <- rep(0:1, each = n)
rmst_fast(time, event, group = group, control = 0, tau = 10)
#> Restricted mean survival time (two-group)
#> 
#>   tau = 10,  control = 0
#>   alternative = two.sided
#> 
#>             RMST
#> control   6.7572
#> treatment 7.8384
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)   
#> difference (treatment - control) 1.0813    0.4094    1.7531 3.154  0.00161 **
#> ratio (treatment / control)      1.1600    1.0572    1.2729 3.133  0.00173 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# \donttest{
# Validation against survRM2
if (requireNamespace("survRM2", quietly = TRUE)) {
  survRM2::rmst2(time, event, group, tau = 10)$unadjusted.result
}
#>                           Est. lower .95 upper .95           p
#> RMST (arm=1)-(arm=0) 1.0812605 0.4094144 1.7531066 0.001608610
#> RMST (arm=1)/(arm=0) 1.1600165 1.0571618 1.2728783 0.001727894
#> RMTL (arm=1)/(arm=0) 0.6665677 0.5119289 0.8679185 0.002596875
# }
```
