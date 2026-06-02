# Core stratified log-rank computation on stratum-blocked sorted vectors

Internal C++ function that computes the stratified log-rank totals O1,
E1, and V1 from a dataset sorted by stratum first and by time within
stratum. Each stratum is processed as an independent contiguous block
using the same single-pass two-group scan as `logrank_core`, and the
per-stratum O1, E1, and V1 are accumulated into the returned totals.
Tied event times are processed atomically within each stratum. A stratum
that contains only one group contributes zero to all three totals,
matching the behaviour of
[`survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html) with a
[`strata()`](https://rdrr.io/pkg/survival/man/strata.html) term. Not
intended to be called directly by users; use
[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
with the `strata` argument instead.

## Usage

``` r
stratified_logrank_core(time_sorted, event_sorted, j_sorted, strata_sorted)
```

## Arguments

- time_sorted:

  A numeric vector of follow-up times sorted by stratum first and by
  ascending time within each stratum.

- event_sorted:

  An integer vector of event indicators (1 = event, 0 = censored),
  aligned with `time_sorted`.

- j_sorted:

  An integer vector of group indicators (1 = treatment, 0 = control),
  aligned with `time_sorted`.

- strata_sorted:

  An integer vector of stratum identifiers, aligned with `time_sorted`.
  Rows of the same stratum must be contiguous, which is guaranteed when
  the data are sorted by stratum first.

## Value

A numeric vector of length 3: `c(O1, E1, V1)`, the stratified observed
events, expected events, and variance for the treatment group.
