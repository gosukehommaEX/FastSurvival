# Fast Sequential Analysis of Simulated Trial Data

Performs interim or sequential analyses of simulated two-group
time-to-event data at one or more analysis times ("looks"). Each look is
defined either by a target cumulative number of events
(information-based timing) or by a calendar time (calendar-based
timing). At every look the data are administratively censored at the
corresponding calendar cutoff, and the requested statistics are computed
for each simulated trial by reusing
[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md),
[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md),
[`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md),
[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md),
[`maxcombo_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md),
and
[`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md).
Optionally the same statistics are also reported within each subgroup.
The censoring and time sorting are handled by a C++ backend, and the
analysis cores are called with `presorted = TRUE`, so each look avoids a
redundant sort.

## Usage

``` r
analysis_fast(
  data,
  control,
  event.looks = NULL,
  time.looks = NULL,
  stat = "logrank",
  tau = NULL,
  t.eval = NULL,
  conf.int = 0.95,
  side = 2,
  by.subgroup = FALSE,
  weight = c("logrank", "fh", "mwlrt", "gehan", "tarone-ware"),
  rho = 0,
  gamma = 0,
  t_star = NULL,
  strata = NULL,
  mc.rho = c(0, 0, 1, 1),
  mc.gamma = c(0, 1, 0, 1),
  abseps = 1e-05,
  maxpts = 25000
)
```

## Arguments

- data:

  A data frame from
  [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
  for a two-group trial, containing at least `sim`, `group`,
  `accrual_time`, `tte`, and `event`.

- control:

  A scalar value indicating which level of `group` represents the
  control group.

- event.looks:

  A numeric vector of target cumulative event counts, one per look.
  Mutually exclusive with `time.looks`.

- time.looks:

  A numeric vector of calendar times, one per look. Mutually exclusive
  with `event.looks`.

- stat:

  A character vector naming the statistics to compute. Any subset of
  `"logrank"`, `"coxph"`, `"rmst"`, `"km"`, `"maxcombo"`, and `"ahsw"`.
  Defaults to `"logrank"`.

- tau:

  A single positive numeric value, the restriction horizon for `"rmst"`
  and the truncation time for `"ahsw"`. Required only when `"rmst"` or
  `"ahsw"` is requested.

- t.eval:

  A single positive numeric value, the landmark time for `"km"`.
  Required only when `"km"` is requested.

- conf.int:

  A single numeric value in (0, 1), the confidence level for `"coxph"`,
  `"rmst"`, and `"ahsw"`. Defaults to 0.95.

- side:

  An integer, either 1 or 2. When `side = 2` (default), two-sided
  p-values `2 pnorm(-|z|)` are reported for log-rank, Cox, and RMST, and
  the two-sided max-combo test is used. When `side = 1`, the one-sided
  p-value in the direction of treatment benefit is reported for each of
  those statistics, and the one-sided max-combo test is used. The test
  statistics themselves are always reported with their natural sign, so
  the choice of `side` affects only the p-value columns. For log-rank
  and Cox the benefit direction is a negative Z (the one-sided p-value
  is the lower tail `pnorm(z)`); for RMST it is a positive Z (the upper
  tail `pnorm(-z)`). The AHSW p-values are always two-sided. For
  group-sequential boundary comparisons (for example with gsDesign or
  rpact), align the sign of the reported Z with the boundary convention
  before comparing.

- by.subgroup:

  A logical value. When `TRUE`, the analysis is also reported within
  each subgroup level, and the output gains a `population` column (long
  form). When `FALSE` (default), only the whole-population analysis is
  returned and no `population` column is added. Requires at least one
  subgroup column in `data`.

- weight:

  A character string naming the weight scheme for the `"logrank"`
  statistic. `"logrank"` (default) is the ordinary unweighted test.
  `"fh"` is the Fleming-Harrington G(rho, gamma) test, `"mwlrt"` is the
  modestly-weighted log-rank test (requires `t_star`), `"gehan"` is
  Gehan-Breslow, and `"tarone-ware"` is Tarone-Ware. See
  [`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md).

- rho:

  A numeric Fleming-Harrington first parameter, used only when
  `weight = "fh"`. Defaults to 0.

- gamma:

  A numeric Fleming-Harrington second parameter, used only when
  `weight = "fh"`. Defaults to 0.

- t_star:

  A single non-negative numeric value, the timepoint of the
  modestly-weighted log-rank test. Required only when
  `weight = "mwlrt"`.

- strata:

  An optional character vector naming one or more subgroup columns of
  `data` to use as the stratification variable for the `"logrank"`
  statistic. `NULL` (default) gives the unstratified test. When several
  columns are named their interaction defines the strata. Stratification
  applies only to the `"logrank"` statistic; the other statistics ignore
  it.

- mc.rho:

  A numeric vector of Fleming-Harrington first parameters for the
  `"maxcombo"` statistic, one per component weight. Defaults to
  `c(0, 0, 1, 1)`.

- mc.gamma:

  A numeric vector of Fleming-Harrington second parameters for the
  `"maxcombo"` statistic, aligned with `mc.rho`. Defaults to
  `c(0, 1, 0, 1)`.

- abseps:

  A single positive numeric value, the absolute error tolerance passed
  to the multivariate normal integration of the `"maxcombo"` p-value.
  Defaults to 1e-5.

- maxpts:

  A single positive integer, the maximum number of function evaluations
  for the quasi-Monte-Carlo integration used by the `"maxcombo"` p-value
  when four or more weights are supplied. Defaults to 25000.

## Value

A data frame. When `by.subgroup = FALSE`, it has `nsim * length(looks)`
rows. When `by.subgroup = TRUE`, it has
`nsim * length(looks) * (1 + total subgroup levels)` rows and an extra
`population` column placed after `look.value`. The common columns are
`sim`, `look` (1-based look index), `look.type` (`"event"` or `"time"`),
`look.value` (the requested event count or calendar time), optionally
`population`, `cutoff` (the calendar time used, `NA` when an event
target was not reached), `reached`, `n.enrolled`, and `n.event`,
followed by the columns of the requested statistics. A statistic that
cannot be computed for a row (no events, or an empty group) is `NA`. The
statistic columns are `logrank.z`, `logrank.chisq`, and `logrank.p` for
`"logrank"`; `cox.coef`, `cox.hr`, `cox.se`, `cox.z`, `cox.p`,
`cox.lower`, and `cox.upper` for `"coxph"`; `rmst.ctrl`, `rmst.trt`,
`rmst.diff`, `rmst.diff.lower`, `rmst.diff.upper`, `rmst.z`, and
`rmst.p` for `"rmst"`; `km.surv.ctrl` and `km.surv.trt` for `"km"`;
`maxcombo.stat` and `maxcombo.p` for `"maxcombo"`; and `ahsw.ah.ctrl`,
`ahsw.ah.trt`, `ahsw.rah`, `ahsw.rah.lower`, `ahsw.rah.upper`,
`ahsw.p.rah`, `ahsw.dah`, `ahsw.dah.lower`, `ahsw.dah.upper`, and
`ahsw.p.dah` for `"ahsw"`. The Z columns `logrank.z`, `cox.z`, and
`rmst.z` carry the natural sign of each test, and the p-value columns
follow `side` except for the AHSW p-values, which are two-sided.

## Details

The input `data` is the data frame returned by
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
for a two-group trial. The columns `sim`, `group`, `accrual_time`,
`tte`, and `event` are required.

For a look at calendar time `cutoff`, each subject with accrual time `a`
contributes only if enrolled by then (`a <= cutoff`). The observed time
at the look is `min(tte, cutoff - a)` and the observed event indicator
is the original `event` when the natural event or dropout occurred on or
before `cutoff` (`a + tte <= cutoff`), and zero otherwise
(administrative censoring at the look).

When `event.looks` is supplied, the calendar cutoff for a target of `d`
events is the calendar time of the `d`-th event in that simulated trial,
counted over the whole trial population. If a simulation contains fewer
than `d` events, the target is never reached: the full data are used,
`reached` is `FALSE`, and `cutoff` is `NA`. When `time.looks` is
supplied, the cutoff is the specified calendar time and `reached` is
always `TRUE`. In both cases the cutoff is determined once on the whole
population and then used for the overall analysis and for every subgroup
analysis at that look.

Exactly one of `event.looks` and `time.looks` must be supplied.

The statistics are selected with `stat`, which may name one or more of
`"logrank"`, `"coxph"`, `"rmst"`, `"km"`, `"maxcombo"`, and `"ahsw"`.

The `"logrank"` statistic is configurable. By default it is the ordinary
unweighted, unstratified two-group log-rank test and reproduces the
behavior of earlier versions of this function exactly. A non-default
`weight` selects a weighted log-rank test (Fleming-Harrington,
modestly-weighted, Gehan-Breslow, or Tarone-Ware) for non-proportional
hazards, and a non-`NULL` `strata` selects the stratified test, summing
the per-stratum contributions. The two options combine to give the
stratified weighted log-rank test. Whatever configuration is chosen, the
result is written to the same `logrank.z`, `logrank.chisq`, and
`logrank.p` columns; the unweighted unstratified case is the `K = 1`
single-stratum degenerate form of the general statistic. The columns
`rho`, `gamma`, and `t_star` parametrize `weight`, and `strata` names
one or more subgroup columns of `data` used as the stratification
variable. The stratification is determined on the whole cut data,
independently of the `population` marginalization, so a stratified
overall analysis is the canonical primary test; in a single-subgroup
population the stratum is constant and the stratified test degenerates
to the ordinary one within that subset.

The `"maxcombo"` statistic is the max-combo test of `mc.rho` and
`mc.gamma` Fleming-Harrington weights. Its `maxcombo.stat` is the most
extreme component (`min` of the component Z-scores when `side = 1`, so
that a negative value favors treatment, and the maximum absolute
component when `side = 2`), and `maxcombo.p` is the joint
multivariate-normal p-value, which already follows `side`.

The `"ahsw"` statistic is the average hazard with survival weight of Uno
and Horiguchi on the window from 0 to `tau`. It reports the per-group
average hazards, the ratio (RAH) and difference (DAH) contrasts with
their confidence intervals, and two-sided p-values for both contrasts.
The AHSW p-values are always two-sided and do not depend on `side`,
matching
[`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md).

When `by.subgroup = TRUE`, the output is given in long form with a
`population` column. Each `(sim, look)` produces one row for the whole
trial (`population = "overall"`) plus one row per subgroup level of each
subgroup factor in `data`. Subgroup factors are the columns named
`subgroup` or `subgroup1`, `subgroup2`, and so on. Each factor is
marginalized separately, so a factor with levels 1 and 2 yields the
populations `subgroup_1` and `subgroup_2` (or `subgroup1_1`,
`subgroup1_2`, and so on for numbered factors). The look cutoff is
always determined on the whole population, so subgroup rows at a given
look share the same `cutoff`, `reached`, and `look.value`; their
`n.enrolled` and `n.event` are the counts within that subgroup. When
`by.subgroup = FALSE` (default) the output has no `population` column
and one row per `(sim, look)`, matching the whole-population analysis.

## See also

[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md),
[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md),
[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md),
[`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md),
[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md),
[`maxcombo_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md),
[`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md).

## Examples

``` r
df <- simdata_fast(
  nsim       = 50,
  n          = c(150, 150),
  a.time     = c(0, 12),
  a.rate     = 300 / 12,
  e.hazard   = list(list(0.10, 0.07), 0.05),
  prevalence = c(0.5, 0.5),
  seed       = 1
)

# Whole-population analysis at two event-based looks
res1 <- analysis_fast(df, control = 1, event.looks = c(80, 140))
head(res1)
#>   sim look look.value   cutoff reached n.enrolled n.event logrank.z
#> 1   1    1         80 11.76464    TRUE        291      80 -1.171526
#> 2   1    2        140 16.08387    TRUE        300     140 -1.555453
#> 3   2    1         80 11.40147    TRUE        285      80 -2.780165
#> 4   2    2        140 16.62229    TRUE        300     140 -4.862092
#> 5   3    1         80 11.48259    TRUE        286      80 -3.631606
#> 6   3    2        140 16.06186    TRUE        300     140 -5.043861
#>   logrank.chisq    logrank.p
#> 1      1.372472 2.413876e-01
#> 2      2.419433 1.198383e-01
#> 3      7.729316 5.433132e-03
#> 4     23.639938 1.161516e-06
#> 5     13.188562 2.816629e-04
#> 6     25.440534 4.562306e-07

# Stratified log-rank on the subgroup factor
res2 <- analysis_fast(df, control = 1, time.looks = 24,
                      stat = "logrank", strata = "subgroup")
head(res2)
#>   sim look look.value cutoff reached n.enrolled n.event logrank.z logrank.chisq
#> 1   1    1         24     24    TRUE        300     208 -3.059687      9.361684
#> 2   2    1         24     24    TRUE        300     192 -5.125810     26.273925
#> 3   3    1         24     24    TRUE        300     194 -5.532218     30.605431
#> 4   4    1         24     24    TRUE        300     206 -4.736340     22.432915
#> 5   5    1         24     24    TRUE        300     194 -3.741408     13.998135
#> 6   6    1         24     24    TRUE        300     203 -1.584396      2.510311
#>      logrank.p
#> 1 2.215684e-03
#> 2 2.962617e-07
#> 3 3.162074e-08
#> 4 2.176126e-06
#> 5 1.829920e-04
#> 6 1.131036e-01

# Fleming-Harrington G(0, 1) weighted log-rank for delayed effects
res3 <- analysis_fast(df, control = 1, time.looks = 24,
                      stat = "logrank", weight = "fh", rho = 0, gamma = 1)
head(res3)
#>   sim look look.value cutoff reached n.enrolled n.event logrank.z logrank.chisq
#> 1   1    1         24     24    TRUE        300     208 -3.156345      9.962516
#> 2   2    1         24     24    TRUE        300     192 -4.348303     18.907740
#> 3   3    1         24     24    TRUE        300     194 -4.895489     23.965814
#> 4   4    1         24     24    TRUE        300     206 -4.627069     21.409771
#> 5   5    1         24     24    TRUE        300     194 -3.584443     12.848231
#> 6   6    1         24     24    TRUE        300     203 -1.224151      1.498546
#>      logrank.p
#> 1 1.597596e-03
#> 2 1.371949e-05
#> 3 9.806151e-07
#> 4 3.708763e-06
#> 5 3.377984e-04
#> 6 2.208952e-01

# Max-combo and AHSW together
res4 <- analysis_fast(df, control = 1, time.looks = 24,
                      stat = c("maxcombo", "ahsw"), tau = 18)
head(res4)
#>   sim look look.value cutoff reached n.enrolled n.event maxcombo.stat
#> 1   1    1         24     24    TRUE        300     208      3.192382
#> 2   2    1         24     24    TRUE        300     192      5.332615
#> 3   3    1         24     24    TRUE        300     194      5.574249
#> 4   4    1         24     24    TRUE        300     206      4.927248
#> 5   5    1         24     24    TRUE        300     194      3.840689
#> 6   6    1         24     24    TRUE        300     203      1.467686
#>     maxcombo.p ahsw.ah.ctrl ahsw.ah.trt  ahsw.rah ahsw.rah.lower ahsw.rah.upper
#> 1 3.151767e-03   0.07881581  0.05237260 0.6644936      0.5034745      0.8770091
#> 2 9.680876e-08   0.08345052  0.04012456 0.4808185      0.3581851      0.6454384
#> 3 2.486058e-08   0.09089574  0.03990436 0.4390125      0.3278702      0.5878301
#> 4 1.424037e-06   0.09257119  0.04918586 0.5313301      0.4019580      0.7023413
#> 5 2.042602e-04   0.08147862  0.04592903 0.5636943      0.4213291      0.7541639
#> 6 2.403671e-01   0.07356423  0.06111041 0.8307082      0.6176444      1.1172711
#>     ahsw.p.rah    ahsw.dah ahsw.dah.lower ahsw.dah.upper   ahsw.p.dah
#> 1 3.890394e-03 -0.02644321    -0.04446913   -0.008417282 4.038005e-03
#> 2 1.091478e-06 -0.04332597    -0.06151709   -0.025134854 3.040474e-06
#> 3 3.250633e-08 -0.05099138    -0.06965300   -0.032329748 8.534762e-08
#> 4 8.920363e-06 -0.04338533    -0.06288035   -0.023890304 1.289882e-05
#> 5 1.135479e-04 -0.03554959    -0.05370932   -0.017389860 1.246293e-04
#> 6 2.199664e-01 -0.01245382    -0.03244248    0.007534848 2.220316e-01
```
