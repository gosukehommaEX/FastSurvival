# Fast weighted Kaplan-Meier (Pepe-Fleming) two-sample test

Computes the weighted Kaplan-Meier (WKM) statistic of Pepe and Fleming
for comparing two survival curves. The statistic is the weighted
integral of the difference between the Kaplan-Meier estimates of the
treatment and control groups over the observed range, standardized to a
Wald z statistic. Unlike weighted log-rank tests, this test targets the
integrated difference in survival and is sensitive to differences even
when the hazard functions cross. The single scan over the sorted data is
performed in C++ for use inside simulation loops.

## Usage

``` r
wkm_fast(
  time,
  event,
  group,
  control = NULL,
  side = 2,
  conf.level = 0.95,
  weight = c("PF", "sqrtPF", "constant"),
  presorted = FALSE
)
```

## Arguments

- time:

  Numeric vector of event or censoring times.

- event:

  Integer vector, 1 for an event and 0 for censoring.

- group:

  Grouping vector with exactly two distinct levels.

- control:

  The level of `group` that denotes the control group.

- side:

  Either 2 for a two-sided test or 1 for a one-sided test of treatment
  superiority (weighted difference greater than 0).

- conf.level:

  Confidence level for the interval of the weighted difference.

- weight:

  Weight function, one of `"PF"` (Pepe-Fleming combined censoring
  weight), `"sqrtPF"` (its square root) or `"constant"` (weight 1).

- presorted:

  Logical; set to `TRUE` when `time`, `event` and `group` are already
  sorted in ascending order of `time`, to skip the internal sort.
  Defaults to `FALSE`.

## Value

A named numeric vector of class `"wkm_fast"` with the weighted
integrated difference, its standard error and confidence limits, and the
Wald statistics for the test.

## Details

The default weight is the Pepe-Fleming combined censoring weight w(t) =
n G1(t) G2(t) / (n1 G1(t) + n2 G2(t)), where G1 and G2 are the
Kaplan-Meier estimates of the censoring survival functions of the two
groups. This weight gives less weight to regions with heavy censoring
and stabilizes the variance in the tail. The choice `weight = "sqrtPF"`
uses its square root, and `weight = "constant"` uses a weight of 1, in
which case the numerator reduces to the difference in restricted mean
survival time over the observed range. With `weight = "PF"` the result
reproduces `nphsim::wkm.Stat` for data without tied times.

The weighted difference is computed as treatment minus control, so a
positive value and a positive z indicate longer survival under
treatment.

## References

Pepe, M. S., & Fleming, T. R. (1989). Weighted Kaplan-Meier statistics:
a class of distance tests for censored survival data. Biometrics, 45(2),
497-507.

Pepe, M. S., & Fleming, T. R. (1991). Weighted Kaplan-Meier statistics:
large sample and optimality considerations. Journal of the Royal
Statistical Society. Series B (Methodological), 53(2), 341-352.

## Examples

``` r
set.seed(1)
n <- 300
g <- rep(0:1, each = n / 2)
tt <- c(rexp(n / 2, log(2) / 12), rexp(n / 2, log(2) / 16))
cc <- rexp(n, rate = 0.02)
time <- pmin(tt, cc)
event <- as.integer(tt <= cc)
wkm_fast(time, event, group = g, control = 0)
#> Weighted Kaplan-Meier test (Pepe-Fleming, two-group)
#> 
#>   control = 0
#>   weight = PF,  alternative = two.sided
#> 
#>                                             Est. lower 95% upper 95%     z
#> weighted difference (treatment - control) 2.8030    0.4869    5.1192 2.372
#>                                           Pr(>|z|)  
#> weighted difference (treatment - control)   0.0177 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
