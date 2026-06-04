# Fast Kalbfleisch-Prentice average hazard ratio for two groups

Estimates the average hazard ratio of Kalbfleisch and Prentice (1981)
between two groups over the time interval from 0 to `tau`, based on the
Kaplan-Meier estimator of each group's survival function. This is the
Kaplan-Meier (unweighted) special case of the estimator implemented in
the archived `AHR` package, restricted to two groups and recoded in C++
for use inside simulation loops.

## Usage

``` r
ahr_fast(
  time,
  status,
  group,
  tau = NULL,
  null.ahr = 1,
  conf.level = 0.95,
  presorted = FALSE
)
```

## Arguments

- time:

  vector of right-censored event times

- status:

  0/1 (or logical) event indicators, 1 for an event

- group:

  vector with exactly two distinct values identifying the groups. The
  smaller value (or first factor level) is the reference group and the
  average hazard ratio is reported for the other group relative to it.

- tau:

  upper limit of the interval over which the average hazard ratio is
  computed. If `NULL` (default) the largest time observed in both groups
  is used.

- null.ahr:

  value of the average hazard ratio under the null hypothesis used for
  the Z statistic and p-value (default 1)

- conf.level:

  confidence level for the confidence interval (default 0.95)

- presorted:

  if `TRUE`, assume `time` is already sorted in ascending order so that
  each group's observations are also ascending; this skips the internal
  sort (default `FALSE`)

## Value

An object of class `"ahr_fast"`, a list with elements `ahr` (the average
hazard ratio, comparison vs reference), `log.ahr`, `se.loghr`, `lower`,
`upper`, `conf.level`, `z` and `p.value` (the primary test on the theta
/ group-share scale, as in Dormuth et al. 2024 eq. 5), `z.loghr` and
`p.value.loghr` (the equivalent test on the `log(ahr)` scale),
`se.theta` (standard error of the tested comparison-group share),
`null.share`, `null.ahr`, `theta` (the two group shares), `var.theta1`,
`var.theta2`, `tau`, `n` (the two group sizes) and `groups`.

## Details

The estimator works with the group shares of the total hazard. Writing
`S1` and `S2` for the two survival functions, the reference-group share
is `theta1 = -integral(S2 dS1) / (1 - S1(tau) S2(tau))` over `[0, tau]`,
the comparison-group share is `theta2 = 1 - theta1`, and the average
hazard ratio is `ahr = theta2 / theta1`. Under proportional hazards with
hazard ratio `psi`, `ahr` estimates `psi`. A value above 1 indicates
higher hazard (worse survival) in the comparison group. The variance of
`theta1` is the direct Greenwood-based estimator. The primary test is on
the theta (group-share) scale, as in Kalbfleisch and Prentice (1981) and
Dormuth et al. (2024, eq. 5): the comparison-group share is compared
with its null value (0.5 when `null.ahr = 1`). An equivalent test and a
confidence interval on the `log(ahr)` scale are also reported.

This is distinct from
[`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md),
which estimates the Uno-Horiguchi average hazard with survival weight.

## References

J. D. Kalbfleisch and R. L. Prentice. Estimation of the average hazard
ratio. Biometrika, 68(1):105-112, 1981.

I. Dormuth, M. Pauly, G. Rauch, and C. Herrmann. Sample size calculation
under nonproportional hazards using average hazard ratios. Biometrical
Journal, 66(6):e202300271, 2024.

## See also

[`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md),
[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)

## Examples

``` r
set.seed(1)
n <- 200
time1 <- rexp(n, 0.1)
time2 <- rexp(n, 0.18)
cens <- rexp(2 * n, 0.05)
obs <- pmin(c(time1, time2), cens)
status <- as.integer(c(time1, time2) <= cens)
group <- rep(c(0, 1), each = n)
ahr_fast(obs, status, group, tau = 8)
#> 
#> Kalbfleisch-Prentice average hazard ratio
#> Comparison group '1' vs reference '0' over [0, 8]
#> 
#>   AHR          : 2.248
#>   95% CI      : 1.699, 2.973
#>   log(AHR)     : 0.81  (se 0.1427)
#>   Z (H0: AHR = 1) : 6.317  [theta scale]
#>   p-value      : 2.661e-10  [theta scale]
#>   (log scale: Z = 5.676, p = 1.379e-08)
#> 
#> Group shares of total hazard (theta):
#>              0        1
#> theta   0.3079   0.6921
#> n     200.0000 200.0000
#> 
```
