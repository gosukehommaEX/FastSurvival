# Fast Simulation of Two-Group Time-to-Event Trial Data

Simulates time-to-event trial data for one or two groups across many
simulated trials, with piecewise accrual, piecewise-exponential survival
and dropout, and optional subgroups defined by a prevalence
specification. The entire generation pipeline (accrual, survival,
dropout, derived columns, and two-group interleaving) runs in a single
C++ kernel that materializes the output data frame once, avoiding
intermediate R-level vector operations and copies. The random-number
stream is consumed in the same order as a per-group reference
implementation, so results are reproducible from `seed`.

## Usage

``` r
simdata_fast(
  nsim = 1000,
  n,
  alloc = c(1, 1),
  a.time,
  a.rate = NULL,
  a.prop = NULL,
  e.hazard = NULL,
  e.median = NULL,
  e.time = NULL,
  d.hazard = NULL,
  d.median = NULL,
  d.time = NULL,
  seed = NULL,
  prevalence = NULL,
  fixed.alloc = FALSE,
  h01.hazard = NULL,
  h01.median = NULL,
  h01.time = NULL,
  h02.hazard = NULL,
  h02.median = NULL,
  h02.time = NULL,
  h12.hazard = NULL,
  h12.median = NULL,
  h12.time = NULL,
  switch.prop = NULL,
  h12.switch.hazard = NULL,
  h12.switch.median = NULL,
  h12.switch.time = NULL,
  switch.clock = "reset"
)
```

## Arguments

- nsim:

  Number of simulated trials.

- n:

  Either a single total sample size (split by `alloc`) or a length-two
  vector of per-group sample sizes.

- alloc:

  A length-two allocation ratio, used when `n` is scalar.

- a.time:

  A numeric vector of accrual-interval breakpoints.

- a.rate:

  Absolute accrual rates (subjects per unit time), interpreted in one of
  two ways. With length `length(a.time) - 1` the accrual period is fully
  specified and the rates must accrue exactly `sum(n)` subjects (an
  inconsistent total is an error). With length `length(a.time)` the
  final rate applies to an open last interval whose end time is computed
  so the total is `sum(n)`. Supply exactly one of `a.rate` and `a.prop`.

- a.prop:

  Accrual proportions, one per accrual interval (length
  `length(a.time) - 1`), giving the fraction of subjects enrolled in
  each interval. Values are normalized to sum to one and distribute the
  fixed total `sum(n)`. Unlike `a.rate` this carries no rate, so the
  accrual period must be fully specified by `a.time`. Supply exactly one
  of `a.rate` and `a.prop`.

- e.hazard:

  Survival hazard(s). A scalar or vector for one group, or a two-element
  list for two groups; per-cell lists are used with subgroups.

- e.median:

  Survival median(s); an alternative to `e.hazard`.

- e.time:

  Survival breakpoints for piecewise hazards (last element `Inf`).

- d.hazard:

  Dropout hazard(s), same structure as `e.hazard`.

- d.median:

  Dropout median(s); an alternative to `d.hazard`.

- d.time:

  Dropout breakpoints for piecewise hazards.

- seed:

  Optional integer seed for the `dqrng` generator.

- prevalence:

  Optional subgroup prevalence specification (numeric vector, list of
  vectors, array, or a named `control`/`treatment` list for
  group-specific prevalence).

- fixed.alloc:

  Logical; when `TRUE` subgroup sizes are deterministic rather than
  drawn.

- h01.hazard:

  Transition hazard(s) for the non-terminal (intermediate) event (state
  0 to state 1) in the illness-death model. A scalar or vector for one
  group, or a two-element list for two groups. Supplying any of `h01.*`
  or `h02.*` activates the illness-death model with two correlated
  endpoints, and is mutually exclusive with `e.hazard` / `e.median`.

- h01.median:

  Median(s) for the intermediate event; an alternative to `h01.hazard`.

- h01.time:

  Breakpoints for a piecewise `h01.hazard` (last element `Inf`).

- h02.hazard:

  Transition hazard(s) for the terminal event without a prior
  intermediate event (state 0 to state 2). Same group and piecewise
  conventions as `h01.hazard`.

- h02.median:

  Median(s) for the direct terminal event; an alternative to
  `h02.hazard`.

- h02.time:

  Breakpoints for a piecewise `h02.hazard`.

- h12.hazard:

  Transition hazard(s) for the terminal event after an intermediate
  event (state 1 to state 2) for subjects who do not switch. Defaults to
  `h02.hazard`, which gives the Fleischer maximal-independence model
  (Fleischer Theorem 1 when there is no switching).

- h12.median:

  Median(s) for the post-event terminal event; an alternative to
  `h12.hazard`.

- h12.time:

  Breakpoints for a piecewise `h12.hazard`, measured from the
  intermediate-event time (clock-reset).

- switch.prop:

  Probability that a subject with an intermediate event switches
  treatment, a scalar or a two-element list (control, treatment).
  Defaults to zero (no switching); the treatment group is typically left
  at zero.

- h12.switch.hazard:

  Transition hazard(s) from state 1 to state 2 for subjects who switch.
  Required when any `switch.prop` is positive.

- h12.switch.median:

  Median(s) for the post-switch terminal event; an alternative to
  `h12.switch.hazard`.

- h12.switch.time:

  Breakpoints for a piecewise `h12.switch.hazard`, measured from the
  switch (intermediate-event) time (clock-reset).

- switch.clock:

  Time origin for the post-event hazards. Currently only `"reset"`
  (measured from the intermediate event) is implemented.

## Value

A `data.frame` with `nsim * sum(n)` rows. The columns are `sim`,
`group`, any subgroup columns, `accrual_time`, `surv_time`,
`dropout_time`, `tte`, `event`, and `calendar_time`. In the
illness-death model the columns are instead `sim`, `group`,
`accrual_time`, `e1_surv_time`, `e2_surv_time`, `dropout_time`,
`e1_tte`, `e1_event`, `e2_tte`, `e2_event`, `e1_calendar_time`,
`e2_calendar_time`, `intermediate`, `switched`, and `switch_time`, where
`e1` is the first (state-0 exit) endpoint and `e2` is the terminal
endpoint. In oncology `e1` is progression-free survival, `e2` is overall
survival, and `intermediate` flags progression.

## Details

For each subject the observed time-to-event is
`tte = pmin(surv_time, dropout_time)` and `event` is 1 when the survival
time occurs first. The calendar time of the observed event is
`accrual_time + tte`.

The total enrolled is fixed at `sum(n)`. With `a.rate` the rates are
absolute (subjects per unit time): when the accrual period is fully
specified the rates must accrue exactly `sum(n)`, and when one extra
rate is given the end of the final interval is solved so the total is
met. With `a.prop` the values are relative proportions that distribute
`sum(n)` across the fully specified intervals. Each accrual interval
receives a deterministic number of subjects (the rate or proportion
times the group total, rounded to keep the per-group total exact),
placed uniformly within the interval.

Survival and dropout are exponential when a single hazard (or median) is
supplied and piecewise-exponential when a vector is supplied together
with the corresponding `e.time` or `d.time` breakpoints, whose last
element must be `Inf`. Group-specific parameters are supplied as a
two-element list (control first, treatment second).

When `prevalence` is supplied the trial has subgroups. A numeric vector
defines a single factor; a list of numeric vectors defines several
independent factors; a multi-dimensional array defines the joint
distribution of correlated factors. Per-cell hazards may be supplied as
a list with one element per cell. With `fixed.alloc = TRUE` the subgroup
sizes are deterministic; otherwise subgroup membership is drawn from the
prevalence distribution.

## See also

[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)

## Examples

``` r
# One-group simulation, simple exponential, no dropout
df1 <- simdata_fast(
  nsim     = 100,
  n        = 50,
  a.time   = c(0, 12),
  a.rate   = 50 / 12,
  e.median = 18,
  seed     = 1
)
head(df1)
#>   sim group accrual_time surv_time dropout_time       tte event calendar_time
#> 1   1     1    0.3819659 71.464787          Inf 71.464787     1     71.846753
#> 2   1     1    4.3863151 50.157102          Inf 50.157102     1     54.543417
#> 3   1     1    1.8796625  1.363115          Inf  1.363115     1      3.242777
#> 4   1     1    4.1827192  4.388424          Inf  4.388424     1      8.571144
#> 5   1     1    2.6645346 27.617923          Inf 27.617923     1     30.282457
#> 6   1     1   10.4718355  9.350214          Inf  9.350214     1     19.822049

# Accrual rate with the final interval computed from the total: 20 per unit
# time for the first 12 units, then 30 per unit time until 500 are enrolled
df1b <- simdata_fast(
  nsim     = 100,
  n        = 500,
  a.time   = c(0, 12),
  a.rate   = c(20, 30),
  e.median = 18,
  seed     = 1
)
head(df1b)
#>   sim group accrual_time surv_time dropout_time      tte event calendar_time
#> 1   1     1    0.3819659  35.69947          Inf 35.69947     1      36.08143
#> 2   1     1    4.3863151  19.97617          Inf 19.97617     1      24.36248
#> 3   1     1    1.8796625  30.96883          Inf 30.96883     1      32.84849
#> 4   1     1    4.1827192  45.91519          Inf 45.91519     1      50.09791
#> 5   1     1    2.6645346  15.75625          Inf 15.75625     1      18.42079
#> 6   1     1   10.4718355  29.59420          Inf 29.59420     1      40.06603

# Accrual by proportion: 30% enrolled in [0, 6], 70% in [6, 12]
df1c <- simdata_fast(
  nsim     = 100,
  n        = 50,
  a.time   = c(0, 6, 12),
  a.prop   = c(0.3, 0.7),
  e.median = 18,
  seed     = 1
)
head(df1c)
#>   sim group accrual_time surv_time dropout_time       tte event calendar_time
#> 1   1     1    0.1909830 71.464787          Inf 71.464787     1     71.655770
#> 2   1     1    2.1931575 50.157102          Inf 50.157102     1     52.350260
#> 3   1     1    0.9398312  1.363115          Inf  1.363115     1      2.302946
#> 4   1     1    2.0913596  4.388424          Inf  4.388424     1      6.479784
#> 5   1     1    1.3322673 27.617923          Inf 27.617923     1     28.950190
#> 6   1     1    5.2359177  9.350214          Inf  9.350214     1     14.586131

# Two-group simulation, simple exponential, with dropout
df2 <- simdata_fast(
  nsim     = 100,
  n        = c(60, 60),
  a.time   = c(0, 6, 12),
  a.rate   = c(8, 12),
  e.median = list(18, 24),
  d.hazard = list(0.01, 0.01),
  seed     = 2
)
head(df2)
#>   sim group accrual_time surv_time dropout_time       tte event calendar_time
#> 1   1     1     3.356372 49.398104    338.23794 49.398104     1     52.754476
#> 2   1     1     4.642096 15.082166     68.86648 15.082166     1     19.724261
#> 3   1     1     1.722016 10.150280     34.10893 10.150280     1     11.872297
#> 4   1     1     2.639623 19.258733     52.33605 19.258733     1     21.898356
#> 5   1     1     1.618623  4.988974    111.28931  4.988974     1      6.607597
#> 6   1     1     4.107536 60.757715     46.11484 46.114837     0     50.222374

# One factor with three levels: single subgroup column
df3 <- simdata_fast(
  nsim       = 100,
  n          = c(150, 150),
  a.time     = c(0, 12),
  a.rate     = 300 / 12,
  e.hazard   = list(list(0.10, 0.08, 0.06), 0.05),
  prevalence = c(0.5, 0.3, 0.2),
  seed       = 3
)
head(df3)
#>   sim group subgroup accrual_time surv_time dropout_time       tte event
#> 1   1     1        1     1.777177  5.115282          Inf  5.115282     1
#> 2   1     1        2     6.189047 14.841849          Inf 14.841849     1
#> 3   1     1        1     2.772930 16.503950          Inf 16.503950     1
#> 4   1     1        3     5.795480  7.764437          Inf  7.764437     1
#> 5   1     1        3     6.756990 42.481282          Inf 42.481282     1
#> 6   1     1        2    11.895699  8.678879          Inf  8.678879     1
#>   calendar_time
#> 1       6.89246
#> 2      21.03090
#> 3      19.27688
#> 4      13.55992
#> 5      49.23827
#> 6      20.57458

# Two independent factors (2 x 2): columns subgroup1 and subgroup2.
# Four cells in column-major order: (1,1), (2,1), (1,2), (2,2).
df4 <- simdata_fast(
  nsim       = 100,
  n          = 200,
  a.time     = c(0, 12),
  a.rate     = 200 / 12,
  e.hazard   = list(0.10, 0.08, 0.07, 0.05),
  prevalence = list(c(0.5, 0.5), c(0.6, 0.4)),
  seed       = 4
)
head(df4)
#>   sim group subgroup1 subgroup2 accrual_time  surv_time dropout_time        tte
#> 1   1     1         1         1    9.1022179  5.9978017          Inf  5.9978017
#> 2   1     1         1         1    2.8634934  0.7745617          Inf  0.7745617
#> 3   1     1         1         2    7.8409159  4.3302474          Inf  4.3302474
#> 4   1     1         1         1   11.1231621 18.1107495          Inf 18.1107495
#> 5   1     1         2         1    0.1179968 21.9852009          Inf 21.9852009
#> 6   1     1         1         2    8.8916664 15.3635872          Inf 15.3635872
#>   event calendar_time
#> 1     1     15.100020
#> 2     1      3.638055
#> 3     1     12.171163
#> 4     1     29.233912
#> 5     1     22.103198
#> 6     1     24.255254

# Two correlated factors via a joint-distribution array (2 x 2)
df5 <- simdata_fast(
  nsim       = 100,
  n          = 200,
  a.time     = c(0, 12),
  a.rate     = 200 / 12,
  e.hazard   = 0.08,
  prevalence = array(c(0.40, 0.10, 0.15, 0.35), dim = c(2, 2)),
  seed       = 5
)
head(df5)
#>   sim group subgroup1 subgroup2 accrual_time surv_time dropout_time       tte
#> 1   1     1         1         2     2.799254 27.395432          Inf 27.395432
#> 2   1     1         2         1    11.012518  2.739611          Inf  2.739611
#> 3   1     1         1         1     9.539372 16.660831          Inf 16.660831
#> 4   1     1         1         2     8.875943 25.432814          Inf 25.432814
#> 5   1     1         1         1     6.182424  5.134396          Inf  5.134396
#> 6   1     1         1         1     7.800613 10.776088          Inf 10.776088
#>   event calendar_time
#> 1     1      30.19469
#> 2     1      13.75213
#> 3     1      26.20020
#> 4     1      34.30876
#> 5     1      11.31682
#> 6     1      18.57670

# Two correlated endpoints, no switching (reduces to Fleischer Theorem 1
# because h12 defaults to h02). In oncology e1 is PFS and e2 is OS; here the
# control has faster intermediate events and faster direct terminal events.
dfid <- simdata_fast(
  nsim       = 100,
  n          = c(150, 150),
  a.time     = c(0, 12),
  a.rate     = 300 / 12,
  h01.median = list(8, 12),
  h02.median = list(24, 30),
  seed       = 6
)
head(dfid)
#>   sim group accrual_time e1_surv_time e2_surv_time dropout_time    e1_tte
#> 1   1     1     3.210288     1.262139    86.027025          Inf  1.262139
#> 2   1     1     5.488111    13.830983    91.509522          Inf 13.830983
#> 3   1     1    10.250797     2.811514     6.206717          Inf  2.811514
#> 4   1     1     9.997081    10.769714    19.784649          Inf 10.769714
#> 5   1     1     6.210510    15.100728    56.583233          Inf 15.100728
#> 6   1     1     7.528979     7.054905    10.651538          Inf  7.054905
#>   e1_event    e2_tte e2_event e1_calendar_time e2_calendar_time intermediate
#> 1        1 86.027025        1         4.472427         89.23731            1
#> 2        1 91.509522        1        19.319094         96.99763            1
#> 3        1  6.206717        1        13.062311         16.45751            1
#> 4        1 19.784649        1        20.766795         29.78173            1
#> 5        1 56.583233        1        21.311238         62.79374            1
#> 6        1 10.651538        1        14.583884         18.18052            1
#>   switched switch_time
#> 1        0          NA
#> 2        0          NA
#> 3        0          NA
#> 4        0          NA
#> 5        0          NA
#> 6        0          NA

# With treatment switching: 40 percent of control subjects who reach the
# intermediate event switch and then follow a more favorable post-event hazard
dfsw <- simdata_fast(
  nsim              = 100,
  n                 = c(150, 150),
  a.time            = c(0, 12),
  a.rate            = 300 / 12,
  h01.median        = list(8, 12),
  h02.median        = list(24, 30),
  switch.prop       = list(0.4, 0),
  h12.switch.median = list(36, 36),
  seed              = 7
)
head(dfsw)
#>   sim group accrual_time e1_surv_time e2_surv_time dropout_time     e1_tte
#> 1   1     1     7.559415   14.4396084   36.8236195          Inf 14.4396084
#> 2   1     1     1.437074    0.6125605    0.6125605          Inf  0.6125605
#> 3   1     1     4.659691   61.4042624   61.4042624          Inf 61.4042624
#> 4   1     1     8.729861    0.6111453   83.4156438          Inf  0.6111453
#> 5   1     1     1.966126   10.6181068  101.0631930          Inf 10.6181068
#> 6   1     1    11.790063    3.0603097   87.5269704          Inf  3.0603097
#>   e1_event      e2_tte e2_event e1_calendar_time e2_calendar_time intermediate
#> 1        1  36.8236195        1        21.999024        44.383035            1
#> 2        1   0.6125605        1         2.049635         2.049635            0
#> 3        1  61.4042624        1        66.063953        66.063953            0
#> 4        1  83.4156438        1         9.341006        92.145505            1
#> 5        1 101.0631930        1        12.584233       103.029319            1
#> 6        1  87.5269704        1        14.850373        99.317034            1
#>   switched switch_time
#> 1        0          NA
#> 2        0          NA
#> 3        0          NA
#> 4        1   0.6111453
#> 5        0          NA
#> 6        0          NA
```
