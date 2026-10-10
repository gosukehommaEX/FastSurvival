# Fast Treatment Switching in Simulated Trial Data

Applies a treatment-switching rule to simulated trial data, such as the
output of
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md),
by changing only the outcomes that occur after each subject's switch.
Subjects of the selected groups switch either at their intermediate
event (for example disease progression), at an opening time after an
interim analysis (crossover at a milestone), or at the later of the two,
and the remaining time to the terminal event after the switch is either
multiplied by an acceleration factor or redrawn from a new (piecewise)
exponential hazard. All simulated trials are processed at once with
vectorized operations. Because the outcomes observed before the switch
are kept exactly, an analysis at or before the opening time gives the
same result before and after switching, so a crossover that depends on
an interim result can be simulated by analyzing the interim, deciding
per simulated trial, and switching afterward.

## Usage

``` r
switch_fast(
  data,
  group,
  prob = 1,
  when = c("intermediate", "cutoff", "later"),
  cutoff = NULL,
  delay = 0,
  sims = NULL,
  aft.factor = NULL,
  hazard = NULL,
  median = NULL,
  time = NULL,
  seed = NULL,
  stream = NULL
)
```

## Arguments

- data:

  A data frame from
  [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
  (single-endpoint, multi-arm, or illness-death), containing at least
  `sim`, `group`, `accrual_time`, and the latent columns described in
  Details.

- group:

  A vector of group labels whose subjects may switch, for example the
  control group `1`.

- prob:

  A single probability in \[0, 1\] that an eligible subject switches.
  Defaults to 1.

- when:

  A character string, one of `"intermediate"`, `"cutoff"`, and
  `"later"`, giving the switch time (see Details).

- cutoff:

  Per-simulation calendar times at which switching opens, used with
  `when = "cutoff"` or `"later"`. Either a numeric vector with one
  element per simulated trial or a one-column matrix from
  [`cutoff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md).
  As in the `cutoff.looks` argument of
  [`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md),
  the names of a vector (or the row names of a matrix) are matched to
  the values of `data$sim`; without names the elements are taken in the
  order of the sorted distinct values of `data$sim`. `NA` disables
  switching in that simulated trial.

- delay:

  A single non-negative number added to `cutoff`, for example the time
  needed to implement the crossover after an interim analysis. Defaults
  to 0.

- sims:

  An optional logical vector with one element per simulated trial (in
  the order of the sorted distinct values of `data$sim`). Switching is
  applied only in the simulated trials with `TRUE`. Requires
  `when = "cutoff"` or `"later"`.

- aft.factor:

  A single positive acceleration factor applied to the remaining time
  after the switch. Exactly one of `aft.factor`, `hazard`, and `median`
  must be supplied.

- hazard:

  Post-switch hazard, a single value or a vector for a piecewise
  exponential distribution with breakpoints `time`, measured from the
  switch.

- median:

  A single post-switch median; an alternative to a single `hazard`.

- time:

  Breakpoints for a piecewise `hazard`, starting at 0 and ending with
  `Inf`, measured from the switch.

- seed:

  Optional integer seed for the `dqrng` generator. If `NULL` (default),
  the draws come from the current state of that generator, which
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not control;
  supply `seed` for reproducible results (see
  [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)).

- stream:

  Optional non-negative whole number selecting a `dqrng` stream;
  requires `seed`. See
  [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md).

## Value

The data frame `data` with the latent and observed columns of the
switched subjects updated, and the columns `switched` (0 or 1) and
`switch_time` (time from accrual to the switch, `NA` for subjects who
did not switch) added or updated.

## Details

The function works on the latent columns of
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md).
For a single-endpoint simulation these are `surv_time` and
`dropout_time`; for an illness-death simulation they are `e1_surv_time`
(time to the first event, such as progression-free survival),
`e2_surv_time` (time to the terminal event, such as overall survival),
`dropout_time`, and `intermediate`. The observed columns (`tte`,
`event`, and `calendar_time`, or their `e1_` and `e2_` counterparts) are
recomputed from the modified latent times, and the columns `switched`
and `switch_time` (time from accrual to the switch) record the switches.
The columns `sim`, `group`, `accrual_time`, and the subgroup columns are
unchanged.

The switch time `s`, measured from accrual, is

- `when = "intermediate"`: the intermediate event time `e1_surv_time` of
  a subject whose intermediate event occurred (`intermediate == 1`);
  illness-death data only;

- `when = "cutoff"`: the opening time `cutoff + delay` of the subject's
  simulated trial minus the accrual time (zero for a subject accrued
  after the opening);

- `when = "later"`: the later of the intermediate event time and the
  opening time, for a subject whose intermediate event occurred;
  illness-death data only.

A subject of the selected `group` switches with probability `prob` when
the switch occurs before the terminal event and before dropout
(`s < e2_surv_time`, or `s < surv_time` for single-endpoint data, and
`s < dropout_time`) and the subject has not already switched
(`switched == 1`, for example at progression through the `switch.prop`
argument of
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)).
With `when = "cutoff"` or `"later"`, the switch never occurs before the
opening time, so every outcome observed by calendar time
`cutoff + delay` is left unchanged; subjects accrued after the opening
are also eligible, as in a protocol that opens crossover from that time
on.

The effect of switching on the terminal time `T > s` is

- `aft.factor = f`: `T_new = s + f * (T - s)`, the causal accelerated
  failure time model of the rank-preserving structural failure time
  (RPSFT) method, under which `f > 1` prolongs the remaining survival;

- `hazard` (or `median`): `T_new = s + R`, where the remaining time `R`
  is drawn from the (piecewise) exponential distribution with hazard
  `hazard` and breakpoints `time` measured from the switch.

In illness-death data with `when = "cutoff"`, a switch can occur before
the intermediate event. The acceleration factor is then applied in the
same way to the first-event time when it lies after the switch, which
keeps the first event no later than the terminal event; a redrawn hazard
is not available in this case because the latent time to the
intermediate event after the terminal event is not stored.

`sims` restricts the switching to some simulated trials, for example
those in which the interim analysis crossed an efficacy boundary. It is
available with `when = "cutoff"` or `"later"` only, so that the decision
cannot change outcomes observed before it is made.

The random numbers (one uniform per row when `prob < 1`, and one
exponential per row when `hazard` or `median` is used) are drawn with
`dqrng`, in row order, so the result is reproducible from `seed` and
`stream`, or from the state left by a seeded
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
call.

## References

Robins, J. M., & Tsiatis, A. A. (1991). Correcting for non-compliance in
randomized trials using rank preserving structural failure time models.
*Communications in Statistics - Theory and Methods*, *20*(8), 2609-2631.

## See also

[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md),
[`cutoff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/cutoff_fast.md),
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md).

## Examples

``` r
# Crossover at an interim analysis: control subjects still on study at the
# 150th event switch to the experimental treatment, which multiplies their
# remaining survival time by 1.5.
df <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
                   a.rate = 300 / 12, e.median = list(12, 18),
                   d.hazard = 0.01, seed = 1)
ia <- cutoff_fast(df, event.looks = 150)
dfs <- switch_fast(df, group = 1, when = "cutoff", cutoff = ia,
                   aft.factor = 1.5)
mean(dfs$switched[dfs$group == 1])
#> [1] 0.3404

# The interim analysis is unchanged; the final analysis is diluted.
r0 <- analysis_fast(df,  control = 1, cutoff.looks = ia)
r1 <- analysis_fast(dfs, control = 1, cutoff.looks = ia)
all.equal(r0$logrank.z, r1$logrank.z)
#> [1] TRUE
fa0 <- analysis_fast(df,  control = 1, event.looks = 220)
fa1 <- analysis_fast(dfs, control = 1, event.looks = 220)
c(mean(fa0$logrank.z), mean(fa1$logrank.z))
#> [1] -3.078266 -2.112271

# Switching at progression in an illness-death simulation: 50 percent of
# the control subjects who progress switch, with a post-progression median
# of 15 instead of 10.
id <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
                   a.rate = 300 / 12, h01.median = list(6, 9),
                   h02.median = list(20, 24), h12.median = list(10, 10),
                   seed = 2)
ids <- switch_fast(id, group = 1, prob = 0.5, when = "intermediate",
                   median = 15)
table(ids$group, ids$switched)
#>    
#>        0    1
#>   1 4587 2913
#>   2 7500    0
```
