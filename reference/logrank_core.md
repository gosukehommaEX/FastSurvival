# Core log-rank computation on pooled sorted vectors (C++ backend)

Internal C++ function that computes O1, E1, and V1 for the log-rank test
directly from a pooled sorted dataset and a group indicator. Eliminates
the R-level group-splitting copies (`time[is1]`, `time[!is1]`,
`event[is1]`, `event[!is1]`) of the previous design by walking the
pooled vector once and updating per-group at-risk counters in C++. Tied
event times are processed atomically. Not intended to be called directly
by users; use
[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
instead.

## Usage

``` r
logrank_core(time_sorted, event_sorted, j_sorted)
```

## Arguments

- time_sorted:

  A numeric vector of pooled follow-up times sorted in ascending order.

- event_sorted:

  An integer vector of event indicators (1 = event, 0 = censored),
  aligned with `time_sorted`.

- j_sorted:

  An integer vector of group indicators (1 = treatment, 0 = control),
  aligned with `time_sorted`.

## Value

A numeric vector of length 3: `c(O1, E1, V1)`.
