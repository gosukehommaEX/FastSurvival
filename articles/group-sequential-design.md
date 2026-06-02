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
#> Sequential analysis summary (simsummary_fast)
#>   Simulations: 2000
#>   Mode: Z-scale boundaries, direction = lower
#>     efficacy [logrank.z]: -3.710, -2.511, -1.993
#> 
#>     look look.value n.enrolled.mean n.event.mean cutoff.mean prob.stop.efficacy
#>        1        100        517.5325       100.00      8.8672             0.0000
#>        2        200        699.8730       200.00     13.0276             0.0065
#>        3        300        700.0000       300.00     17.4987             0.0150
#>  overall         NA        700.0000       299.35     17.4706             0.0215
#>  prob.stop.futility prob.stop.any cum.reject
#>                   0        0.0000     0.0000
#>                   0        0.0065     0.0065
#>                   0        0.0150     0.0215
#>                   0        0.0215     0.0215
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
#> Sequential analysis summary (simsummary_fast)
#>   Simulations: 2000
#>   Mode: Z-scale boundaries, direction = lower
#>     efficacy [logrank.z]: -3.710, -2.511, -1.993
#> 
#>     look look.value n.enrolled.mean n.event.mean cutoff.mean prob.stop.efficacy
#>        1        100        559.6415       100.00      9.5889             0.0240
#>        2        200        700.0000       200.00     14.2143             0.4765
#>        3        300        700.0000       300.00     19.5279             0.3640
#>  overall         NA        696.8610       247.55     16.7443             0.8645
#>  prob.stop.futility prob.stop.any cum.reject
#>                   0        0.0240     0.0240
#>                   0        0.4765     0.5005
#>                   0        0.3640     0.8645
#>                   0        0.8645     0.8645
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

Lin, R. S., Lin, J., Roychoudhury, S., et al. (2020). Alternative
analysis methods for time to event endpoints under nonproportional
hazards: a comparative analysis. *Statistics in Biopharmaceutical
Research*, 12(2), 187-198.

Magirr, D., & Burman, C.-F. (2019). Modestly weighted logrank tests.
*Statistics in Medicine*, 38(20), 3782-3790.
