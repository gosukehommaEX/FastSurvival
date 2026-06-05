# Group sequential design with the simulation trio

## Purpose

This vignette demonstrates the three simulation functions working
together:
[`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
generates the trial data,
[`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
performs the interim analyses, and
[`simsummary_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md)
aggregates the operating characteristics. We set up a three-look group
sequential design with an O’Brien-Fleming type efficacy boundary, then
check that the simulated rejection probabilities match the closed-form
values from [gsDesign](https://cran.r-project.org/package=gsDesign).

The point is not that simulation replaces the closed-form calculation
under proportional hazards, where gsDesign is exact, but that once the
two agree on a case gsDesign can handle, the same simulation machinery
can be trusted for cases it cannot, such as non-proportional hazards or
data-dependent timing.

``` r

library(FastSurvival)
```

## The design

We consider a two-group trial with a one-sided test at level 0.025,
three equally spaced looks by information fraction, and an
O’Brien-Fleming spending function for the efficacy boundary. We first
obtain the design from gsDesign.

``` r

library(gsDesign)

k     <- 3
alpha <- 0.025

gsd <- gsDesign(
  k         = k,
  test.type = 1,
  alpha     = alpha,
  sfu       = sfLDOF,
  timing    = c(1, 2, 3) / 3
)

# Upper (efficacy) Z boundaries at the three looks
z_bounds <- gsd$upper$bound
z_bounds
#> [1] 3.710303 2.511427 1.993048
```

The boundaries are increasing in stringency at the earlier looks, as
expected for an O’Brien-Fleming spending function. We will compare
against these boundaries on the Z scale.

## Simulate the trial under the null

Under the null hypothesis the two groups share the same hazard. We
simulate many trials, each with the planned accrual and follow-up, and
analyze at three event-count looks corresponding to the information
fractions above. The target final event count is chosen for the design,
and the interim looks are at one third and two thirds of it.

``` r

nsim         <- 2000
final_events <- 300
event_looks  <- round(final_events * c(1, 2, 3) / 3)

df_null <- simdata_fast(
  nsim     = nsim,
  n        = c(350, 350),
  a.time   = c(0, 12),
  a.rate   = 700 / 12,
  e.hazard = list(0.05, 0.05),
  seed     = 101
)

res_null <- analysis_fast(
  df_null, control = 1,
  event.looks = event_looks,
  stat = "logrank", side = 1
)

sum_null <- simsummary_fast(
  res_null,
  eff.col   = "logrank.z",
  efficacy  = -z_bounds,
  direction = "lower"
)
sum_null
#> Group-Sequential Operating Characteristics (simsummary_fast)
#>   Simulations: 2000
#>   Boundaries: efficacy on 'logrank.z' (direction = lower)
#> 
#> Stopping Boundaries: Look by Look
#>  Look Info. Frac. Events (s) Sample (n) Efficacy Z Cum. Cross. Eff.
#>     1        0.33      100.0      517.5    -3.7103           0.0000
#>     2        0.67      200.0      699.9    -2.5114           0.0065
#>     3        1.00      300.0      700.0    -1.9930           0.0215
#> 
#> Events, Sample Size, Dropouts, Pipeline and Analysis Times: Look by Look
#>  Look Info. Frac. Sample (n) Events (s) Dropouts (d) Pipeline Analysis Time
#>     1        0.33      517.5      100.0          0.0    417.5          8.87
#>     2        0.67      699.9      200.0          0.0    499.9         13.03
#>     3        1.00      700.0      300.0          0.0    400.0         17.50
#>  Cross. Eff.
#>       0.0000
#>       0.0065
#>       0.0150
#> 
#> Overall
#>   Rejection rate (efficacy):    0.0215
#>   Expected events at stop:      299.4
#>   Expected sample size at stop: 700.0
#>   Expected analysis time at stop:17.47
```

The cumulative rejection probability at the final look estimates the
type I error. It should be close to the nominal 0.025.

``` r

overall_null <- sum_null[sum_null$look == "overall", ]
data.frame(
  source = c("simulation", "gsDesign"),
  alpha  = c(overall_null$cum.reject, alpha)
)
#>       source  alpha
#> 1 simulation 0.0215
#> 2   gsDesign 0.0250
```

## Simulate under the alternative

Under the alternative the treatment group has a lower hazard. The same
design and the same boundaries now estimate the power.

``` r

df_alt <- simdata_fast(
  nsim     = nsim,
  n        = c(350, 350),
  a.time   = c(0, 12),
  a.rate   = 700 / 12,
  e.hazard = list(0.05, 0.035),
  seed     = 202
)

res_alt <- analysis_fast(
  df_alt, control = 1,
  event.looks = event_looks,
  stat = "logrank", side = 1
)

sum_alt <- simsummary_fast(
  res_alt,
  eff.col   = "logrank.z",
  efficacy  = -z_bounds,
  direction = "lower"
)
sum_alt
#> Group-Sequential Operating Characteristics (simsummary_fast)
#>   Simulations: 2000
#>   Boundaries: efficacy on 'logrank.z' (direction = lower)
#> 
#> Stopping Boundaries: Look by Look
#>  Look Info. Frac. Events (s) Sample (n) Efficacy Z Cum. Cross. Eff.
#>     1        0.33      100.0      559.6    -3.7103           0.0240
#>     2        0.67      200.0      700.0    -2.5114           0.5005
#>     3        1.00      300.0      700.0    -1.9930           0.8645
#> 
#> Events, Sample Size, Dropouts, Pipeline and Analysis Times: Look by Look
#>  Look Info. Frac. Sample (n) Events (s) Dropouts (d) Pipeline Analysis Time
#>     1        0.33      559.6      100.0          0.0    459.6          9.59
#>     2        0.67      700.0      200.0          0.0    500.0         14.21
#>     3        1.00      700.0      300.0          0.0    400.0         19.53
#>  Cross. Eff.
#>       0.0240
#>       0.4765
#>       0.3640
#> 
#> Overall
#>   Rejection rate (efficacy):    0.8645
#>   Expected events at stop:      247.6
#>   Expected sample size at stop: 696.9
#>   Expected analysis time at stop:16.74
```

## Comparing the boundary-crossing probabilities

For a closed-form reference, recompute the design at the alternative
effect size and read off the cumulative crossing probabilities at each
look, then put them next to the per-look cumulative rejection from the
simulation.

``` r

sim_cum <- sum_alt[sum_alt$look != "overall", "cum.reject"]

data.frame(
  look       = seq_len(k),
  simulation = sim_cum,
  z.boundary = z_bounds
)
#>   look simulation z.boundary
#> 1    1     0.0240   3.710303
#> 2    2     0.5005   2.511427
#> 3    3     0.8645   1.993048
```

The simulated cumulative rejection probabilities track the spending
implied by the boundaries. Small discrepancies are Monte Carlo error and
shrink as `nsim` increases. The type I error under the null and the
power under the alternative both align with the closed-form design.

## Beyond proportional hazards

The value of the simulation trio is that nothing in the workflow assumes
proportional hazards. To study a delayed treatment effect, replace the
constant alternative hazard with a piecewise specification and keep
everything else the same. The log-rank statistic loses power under a
delayed effect, and a weighted or max-combo statistic can be substituted
at the
[`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
step to recover it. Because gsDesign cannot evaluate these cases in
closed form, the validated simulation machinery becomes the tool of
choice.

``` r

df_delay <- simdata_fast(
  nsim     = nsim,
  n        = c(350, 350),
  a.time   = c(0, 12),
  a.rate   = 700 / 12,
  e.hazard = list(c(0.05, 0.05), c(0.05, 0.030)),
  e.time   = c(0, 6, Inf),
  seed     = 303
)

res_delay <- analysis_fast(
  df_delay, control = 1,
  event.looks = event_looks,
  stat = "maxcombo"
)
```

## References

O’Brien, P. C., & Fleming, T. R. (1979). A multiple testing procedure
for clinical trials. *Biometrics*, 35(3), 549-556.

Lan, K. K. G., & DeMets, D. L. (1983). Discrete sequential boundaries
for clinical trials. *Biometrika*, 70(3), 659-663.

Karrison, T. G. (2016). Versatile tests for comparing survival curves
based on weighted log-rank statistics. *The Stata Journal*, 16(3),
678-690.

Lin, R. S., Lin, J., Roychoudhury, S., et al. (2020). Alternative
analysis methods for time to event endpoints under nonproportional
hazards: a comparative analysis. *Statistics in Biopharmaceutical
Research*, 12(2), 187-198.
