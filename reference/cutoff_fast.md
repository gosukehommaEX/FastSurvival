# Fast Per-Simulation Analysis Cutoffs from Combined Trigger Rules

Computes the calendar time of each analysis look in every simulated
trial from a combination of trigger rules: a target number of events, a
planned calendar time, a maximum calendar time, a minimum time after the
previous look, and a minimum follow-up after a given number of subjects
are enrolled. The cutoffs are computed for all simulated trials in a
single C++ pass and are returned as a matrix that can be passed to the
`cutoff.looks` argument of
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
and
[`pairwise_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/pairwise_fast.md),
or to the `cutoff` argument of
[`switch_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md).
Separating the determination of the analysis time from the analysis
itself lets the same cutoffs drive several endpoints (for example
overall survival analyzed at the progression-free survival event-driven
looks), several contrasts, or a modification of the data after an
interim analysis.

## Usage

``` r
cutoff_fast(
  data,
  event.looks = NULL,
  time.looks = NULL,
  max.time = NULL,
  min.gap = NULL,
  min.enrolled = NULL,
  min.followup = NULL,
  event.subset = NULL,
  tte.col = "tte",
  event.col = "event"
)
```

## Arguments

- data:

  A data frame of simulated trials, such as the output of
  [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md),
  with columns `sim`, `accrual_time`, and the columns named by `tte.col`
  and `event.col`.

- event.looks:

  Target cumulative event counts. Either a vector with one positive
  whole number per look (`NA` to omit the event condition at a look), or
  a matrix of positive whole numbers with one row per simulated trial
  (in the order of the sorted distinct values of `data$sim`) and one
  column per look.

- time.looks:

  Planned calendar times, one per look (`NA` to omit).

- max.time:

  Maximum calendar times, one per look (`NA` to omit). The cutoff never
  exceeds this value.

- min.gap:

  Minimum time after the previous look, one per look (`NA` to omit). At
  the first look it is measured from time zero.

- min.enrolled:

  Number of enrolled subjects after which the minimum follow-up
  `min.followup` must elapse, one per look (`NA` to omit).

- min.followup:

  Minimum follow-up after the `min.enrolled`-th subject is enrolled, one
  per look. Defaults to zero when `min.enrolled` is supplied, and
  requires `min.enrolled`.

- event.subset:

  An optional logical vector with one element per row of `data`. When
  supplied, only the events of the rows with `TRUE` are counted for
  `event.looks`.

- tte.col:

  A single character string naming the observed-time column used to
  count events. Defaults to `"tte"`.

- event.col:

  A single character string naming the event-indicator column (0 or 1)
  used to count events. Defaults to `"event"`.

## Value

A numeric matrix with one row per simulated trial and one column per
look, holding the calendar cutoffs (`NA` for a look that is not
reached). The row names are the simulation identifiers, and the
attribute `"look.value"` holds, for each look, the event target when a
common target is supplied, otherwise the planned calendar time,
otherwise `NA`.
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
reports this attribute in its `look.value` column.

## Details

For look `l` of a simulated trial the candidate calendar times are

- `t_cal`, the planned calendar time `time.looks[l]`;

- `t_event`, the calendar time (accrual time plus observed time) of the
  `d`-th counted event, where `d` is `event.looks[l]` (or the entry of
  the `event.looks` matrix for that simulation and look);

- `t_gap`, the cutoff of the previous look plus `min.gap[l]` (zero plus
  `min.gap[1]` at the first look);

- `t_enr`, the accrual time of the `min.enrolled[l]`-th enrolled subject
  plus `min.followup[l]`;

and the cutoff is `min(max(t_cal, t_event, t_gap, t_enr), max.time[l])`,
where a condition that is not supplied (or is `NA` at that look) is left
out. This is the rule of `get_analysis_date()` in the simtrial package
(with `max.time` playing the role of its
`max_extension_for_target_event`), applied to every simulated trial at
once. For example, `event.looks = 300` with `max.time = 48` is "300
events or month 48, whichever comes first", and `event.looks = 300` with
`time.looks = 36` is "300 events but not before month 36".

An event or enrollment target that is not met in the simulated data has
no finite calendar time. Unless a `max.time` caps the look, its cutoff
is then `NA`, which
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
reports as `reached = FALSE` with the statistics of the full data, as it
does for an unreached `event.looks` target. In this respect the function
differs from simtrial, which uses the last observed event time for an
unmet target.

Events are counted on the columns named by `tte.col` and `event.col`, so
a look can be triggered by an endpoint other than the one analyzed (for
example the `e1_tte` and `e1_event` columns of an illness-death
simulation). By default every subject's events are counted;
`event.subset` restricts the count to a subset of rows, such as the
control group (`data$group == 1`), one subgroup, or the two arms of the
primary contrast in a multi-arm trial. When the subjects tie at the
target event time, all events at that time are included in the analysis,
so the analyzed event count can exceed the target, as in
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md).

A matrix `event.looks` with one row per simulated trial gives
simulation-specific event targets, for example targets re-estimated at
an interim analysis.

Each look is determined by its own conditions; only `min.gap` refers to
the previous look. A look with `max.time` as its only condition is
placed at that time. When the previous look is not reached, `t_gap` is
infinite, so a look with `min.gap` is not reached either unless its
`max.time` caps it. The looks are not reordered: when the conditions of
a later look give an earlier cutoff than the previous look (for example
`event.looks = c(150, 220)` with `max.time = c(NA, 30)`), a warning is
given. Supplying `min.gap` (for example 0) at the later looks keeps them
in order unless their `max.time` is earlier.

## See also

[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md),
[`pairwise_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/pairwise_fast.md),
[`switch_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/switch_fast.md),
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md).

## Examples

``` r
df <- simdata_fast(nsim = 20, n = c(150, 150), a.time = c(0, 12),
                   a.rate = 300 / 12, e.median = list(12, 18),
                   d.hazard = 0.01, seed = 1)

# 150 events or month 30, whichever comes first, then a final look at
# 220 events but at least 6 months after the interim
cut <- cutoff_fast(df, event.looks = c(150, 220), max.time = c(30, NA),
                   min.gap = c(NA, 6))
head(cut)
#>      look1    look2
#> 1 21.68353 56.12726
#> 2 19.56962 37.92612
#> 3 21.45020 45.57526
#> 4 22.33437 44.48730
#> 5 25.01527 61.62116
#> 6 22.71343 47.55300

res <- analysis_fast(df, control = 1, cutoff.looks = cut)
head(res)
#>   sim look look.value   cutoff reached n.enrolled n.event n.dropout n.pipeline
#> 1   1    1        150 21.68353    TRUE        300     150        34        116
#> 2   1    2        220 56.12726    TRUE        300     220        57         23
#> 3   2    1        150 19.56962    TRUE        300     150        25        125
#> 4   2    2        220 37.92612    TRUE        300     220        38         42
#> 5   3    1        150 21.45020    TRUE        300     150        35        115
#> 6   3    2        220 45.57526    TRUE        300     220        45         35
#>    logrank.z logrank.chisq    logrank.p
#> 1 -3.5830408     12.838182 0.0003396175
#> 2 -2.7397778      7.506383 0.0061480726
#> 3 -0.5443804      0.296350 0.5861797207
#> 4 -2.1257515      4.518820 0.0335239525
#> 5 -3.6000536     12.960386 0.0003181515
#> 6 -2.9619756      8.773300 0.0030567199

# Looks triggered by the events of the control group only
cut_c <- cutoff_fast(df, event.looks = 100, event.subset = df$group == 1)
head(cut_c)
#>      look1
#> 1 28.57429
#> 2 29.06672
#> 3 25.51201
#> 4 26.65380
#> 5 32.54351
#> 6 30.44908
```
