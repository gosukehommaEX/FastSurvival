# Pairwise comparisons of each experimental arm against a shared control

Runs the two-group analysis of
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
once for each experimental arm against a common control arm, on
multi-arm simulated data such as the output of
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
with a length-greater-than-two `n`. The contrasts are analyzed at a
shared set of looks, and the results are stacked into one long data
frame with an `arm` column, optionally with a Bonferroni-adjusted
p-value across the contrasts.

## Usage

``` r
pairwise_fast(
  data,
  control,
  event.looks = NULL,
  time.looks = NULL,
  primary = NULL,
  arms = NULL,
  stat = "logrank",
  adjust = c("none", "bonferroni"),
  p.col = NULL,
  ...
)
```

## Arguments

- data:

  A data frame of simulated trial data with columns `sim`, `group`,
  `accrual_time`, `tte`, and `event`, as produced by
  [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md).

- control:

  The group label of the control arm.

- event.looks:

  A numeric vector of target cumulative event counts for the primary
  contrast, one per look. Mutually exclusive with `time.looks`; requires
  `primary`.

- time.looks:

  A numeric vector of calendar times, one per look. Mutually exclusive
  with `event.looks`.

- primary:

  The group label of the experimental arm whose control-versus-arm
  comparison defines the shared calendar cutoff. Required with
  `event.looks` and ignored with `time.looks`.

- arms:

  A vector of experimental arm labels to compare against the control.
  Defaults to every group other than `control`.

- stat:

  The test statistic passed to
  [`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
  (for example `"logrank"`, `"coxph"`, or `"rmst"`).

- adjust:

  Multiplicity adjustment across contrasts. Either `"none"` (default) or
  `"bonferroni"`, which adds a `p.adj` column.

- p.col:

  The name of the p-value column to adjust. By default it is the column
  matching the chosen `stat` (for example `"logrank.p"`).

- ...:

  Further arguments passed to
  [`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md),
  such as `side`, `conf.level`, `tau`, `rho`, and `gamma`. The argument
  `by.subgroup` is not supported.

## Value

A data frame with one row per contrast, look, and simulation. The first
column is `arm` (the experimental arm compared against `control`),
followed by the columns returned by
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md):
`sim`, `look`, `look.value`, `cutoff`, `reached`, the enrollment and
event counts, and the statistic columns for the chosen `stat`. When
`adjust = "bonferroni"` a `p.adj` column is appended.

## Details

Two timing regimes are supported, exactly one of which must be
requested.

With `time.looks`, every contrast is analyzed at the same fixed calendar
time or times, so all contrasts share the same data cutoff by
construction.

With `event.looks`, the analysis is event-driven and `primary` names the
experimental arm whose control-versus-arm comparison defines the data
cutoff. The primary contrast is analyzed at the requested cumulative
event counts, its per-simulation calendar cutoffs are recorded, and then
every contrast, including the primary one, is analyzed at those same
cutoffs. This reproduces the standard design in which the primary
event-driven analysis fixes a single data cutoff at which all
comparisons are performed. For a simulation in which the primary event
target is not reached at a look, that look is marked with
`reached = FALSE` and `NA` statistics for every contrast.

The Bonferroni option multiplies each p-value by the number of contrasts
and caps it at one, controlling the family-wise error rate across the
control-versus-arm comparisons at each look. Multiplicity across looks
is a separate matter handled by group-sequential boundaries in
[`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md),
not by this adjustment.

This is a single-endpoint helper: it reads the `tte` and `event` columns
and does not support subgroups. Comparisons for a second endpoint are
obtained by calling `pairwise_fast` again on that endpoint's columns.

## See also

[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md),
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)

## Examples

``` r
# Three-arm trial: control (group 1) and two experimental arms.
dfk <- simdata_fast(
  nsim     = 100,
  n        = c(120, 120, 120),
  a.time   = c(0, 12),
  a.rate   = 360 / 12,
  e.median = list(12, 16, 20),
  seed     = 8
)

# Fixed calendar look at month 30, Bonferroni across the two contrasts.
pw <- pairwise_fast(dfk, control = 1, time.looks = 30,
                    stat = "logrank", side = 1, adjust = "bonferroni")
head(pw)
#>   arm sim look look.value cutoff reached n.enrolled n.event n.dropout
#> 1   2   1    1         30     30    TRUE        240     173         0
#> 2   2   2    1         30     30    TRUE        240     160         0
#> 3   2   3    1         30     30    TRUE        240     176         0
#> 4   2   4    1         30     30    TRUE        240     176         0
#> 5   2   5    1         30     30    TRUE        240     168         0
#> 6   2   6    1         30     30    TRUE        240     160         0
#>   n.pipeline  logrank.z logrank.chisq  logrank.p      p.adj
#> 1         67 -2.0654860    4.26623257 0.01943852 0.03887704
#> 2         80 -2.1030267    4.42272144 0.01773172 0.03546343
#> 3         64 -0.6031565    0.36379774 0.27320230 0.54640460
#> 4         64 -2.1304840    4.53896204 0.01656584 0.03313168
#> 5         72  0.1039376    0.01080303 0.54139057 1.00000000
#> 6         80 -2.0645649    4.26242820 0.01948209 0.03896419

# Event-driven: arm 3 is the primary contrast and its 200th control-plus-arm-3
# event fixes the cutoff at which both contrasts are analyzed.
pw2 <- pairwise_fast(dfk, control = 1, event.looks = 200, primary = 3,
                     stat = "logrank", side = 1, adjust = "bonferroni")
head(pw2)
#>   arm sim look look.value   cutoff reached n.enrolled n.event n.dropout
#> 1   2   1    1        200 45.50440    TRUE        240     210        30
#> 2   2   2    1        200 47.69751    TRUE        240     208        32
#> 3   2   3    1        200 45.98303    TRUE        240     210        30
#> 4   2   4    1        200 44.93080    TRUE        240     210        30
#> 5   2   5    1        200 43.31435    TRUE        240     209        31
#> 6   2   6    1        200 55.37777    TRUE        240     220        20
#>   n.pipeline  logrank.z logrank.chisq  logrank.p      p.adj
#> 1          0 -1.6385354     2.6847984 0.05065502 0.10131005
#> 2          0 -2.5714433     6.6123205 0.00506378 0.01012756
#> 3          0 -1.1258351     1.2675046 0.13011767 0.26023534
#> 4          0 -1.7611765     3.1017428 0.03910426 0.07820853
#> 5          0 -0.3822314     0.1461008 0.35114487 0.70228974
#> 6          0 -1.3978726     1.9540478 0.08107566 0.16215132
```
