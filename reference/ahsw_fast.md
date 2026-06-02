# Fast Average Hazard with Survival Weight (Two-Group Comparison)

Computes the average hazard with survival weight (AHSW) of Uno and
Horiguchi for two groups and the between-group contrasts. The average
hazard on the window from 0 to `tau` is the ratio of the cumulative
event probability at `tau` to the restricted mean survival time at
`tau`, both based on the Kaplan-Meier estimate. The function returns the
per-group average hazard, the ratio of average hazards (RAH, treatment
over control) on the log scale, and the difference of average hazards
(DAH, treatment minus control) on the identity scale, each with a
confidence interval and a two-sided test. The C++ backend walks the
pooled sorted data once per group, so the function is suitable for
simulation loops with `presorted = TRUE`.

## Usage

``` r
ahsw_fast(time, event, group, control, tau, conf.int = 0.95, presorted = FALSE)
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

- tau:

  A single positive numeric value, the truncation time point for the
  average hazard. Both groups must have positive Kaplan-Meier survival
  at `tau`.

- conf.int:

  A single numeric value in (0, 1) specifying the confidence level.
  Defaults to 0.95.

- presorted:

  A logical value. If `TRUE`, `time`, `event`, and `group` are assumed
  to be sorted in ascending order of `time`, and the internal
  [`order()`](https://rdrr.io/r/base/order.html) call is skipped. If
  `FALSE` (default), sorting is handled internally.

## Value

An object of class `"ahsw_fast"`, a named numeric vector containing the
per-group average hazards (`ah.ctrl`, `ah.trt`), the ratio contrast
(`rah`, `rah.lower`, `rah.upper`, `p.rah`), and the difference contrast
(`dah`, `dah.lower`, `dah.upper`, `p.dah`). The truncation time and
confidence level are stored as attributes `tau` and `conf.int`, and the
`control` label is also stored. Returns `NA` values (still with class
`"ahsw_fast"`) when either group has zero survival at `tau` or a
non-finite variance.

## Details

The average hazard with survival weight is

AH(tau) = (1 - S(tau)) / integral over \[0, tau\] of S(u) du,

the ratio of the cumulative event probability at `tau` to the restricted
mean survival time at `tau`. It can be read as a general censoring-free
incidence rate on the window from 0 to `tau` and stays interpretable
under non-proportional hazards. This is a different quantity from the
average hazard ratio of Kalbfleisch, which averages the time-varying
ratio of hazards rather than forming a single average hazard per group
and then contrasting.

Writing the treatment and control average hazards as a1 and a0, the
ratio contrast is RAH = a1 / a0, formed on the log scale with variance
v_Q1 / n1 + v_Q0 / n0, and the difference contrast is DAH = a1 - a0,
with variance v_U1 / n1 + v_U0 / n0, using the independence of the two
groups. The per-group variance terms v_Q (log scale) and v_U (identity
scale) follow the asymptotic variance of Uno and Horiguchi, computed
from the Nelson-Aalen increments, the running restricted mean survival
time and the at-risk fraction. The confidence interval for RAH is
exponentiated from the log scale, and the two-sided p-values are based
on the normal approximation.

When `presorted = TRUE`, the inputs are assumed to be sorted in
ascending order of `time`, so the internal
[`order()`](https://rdrr.io/r/base/order.html) call is skipped.
Splitting into groups preserves the ascending order within each group.

## References

Uno, H., & Horiguchi, M. (2023). Ratio and difference of average hazard
with survival weight: new measures to quantify survival benefit of new
therapy. *Statistics in Medicine*, *42*(7), 936-952.

## See also

[`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
for the restricted mean survival time, and
[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)
for the Kaplan-Meier estimate at a time point.

## Examples

``` r
library(survival)

# Average hazard contrasts on the ovarian data
ahsw_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1, tau = 600)
#> Average hazard with survival weight (two-group)
#> 
#>   tau = 600,  control = 1
#> 
#>               AH
#> control   0.0011
#> treatment 0.0008
#> 
#>                                        Est.  lower 95%  upper 95% p.value
#> ratio (treatment / control)       0.7548735  0.2622750  2.1726579   0.602
#> difference (treatment - control) -0.0002646 -0.0012914  0.0007622   0.613

# presorted = TRUE: sort once outside, reuse inside a loop
ord <- order(ovarian$futime)
ahsw_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
          control = 1, tau = 600, presorted = TRUE)
#> Average hazard with survival weight (two-group)
#> 
#>   tau = 600,  control = 1
#> 
#>               AH
#> control   0.0011
#> treatment 0.0008
#> 
#>                                        Est.  lower 95%  upper 95% p.value
#> ratio (treatment / control)       0.7548735  0.2622750  2.1726579   0.602
#> difference (treatment - control) -0.0002646 -0.0012914  0.0007622   0.613

# \donttest{
# Validation against survAH
if (requireNamespace("survAH", quietly = TRUE)) {
  arm <- as.numeric(ovarian$rx == 2)
  survAH::ah2(time = ovarian$futime, status = ovarian$fustat,
              arm = arm, tau = 600)
}
#> 
#> The time window: [eta, tau] = [0, 600] was specified. Warning: The normal approximation may be questionable with the specified tau. A smaller value of tau would be recommended for this data.
#> 
#> Number of observations: 
#>      Total N Event by tau Censor by tau At risk at tau
#> arm0      13            6             2              5
#> arm1      13            5             2              6
#> 
#> 
#> Average Hazard (AH) by arm: 
#>            Est. Lower 0.95 Upper 0.95
#> AH (arm0) 0.001          0      0.002
#> AH (arm1) 0.001          0      0.002
#> 
#> 
#> Between-group contrast: 
#>                               Est. Lower 0.95 Upper 0.95 P-value
#> Ratio of AH (arm1/arm0)      0.755      0.262      2.173   0.602
#> Difference of AH (arm1-arm0) 0.000     -0.001      0.001   0.613
# }
```
