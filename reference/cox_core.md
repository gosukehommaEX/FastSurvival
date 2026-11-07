# Core two-group Cox model fit (C++ backend)

Internal C++ function that fits the Cox proportional hazards model with
a single treatment indicator by Newton-Raphson iteration on the risk-set
counts at the distinct event times, as described in
[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md).
Not intended to be called directly by users; use
[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
instead.

## Usage

``` r
cox_core(
  time_sorted,
  event_sorted,
  j_sorted,
  strata_sorted,
  efron,
  iter_max = 20L
)
```

## Arguments

- time_sorted:

  A numeric vector of pooled follow-up times sorted in ascending order
  (by stratum, then time, when `strata_sorted` is not empty).

- event_sorted:

  An integer vector of event indicators (1 = event, 0 = censored),
  aligned with `time_sorted`.

- j_sorted:

  An integer vector of group indicators (1 = treatment, 0 = control),
  aligned with `time_sorted`.

- strata_sorted:

  An integer vector of stratum codes aligned with `time_sorted`, with
  the rows of each stratum contiguous, or an integer vector of length 0
  for the unstratified model.

- efron:

  `TRUE` for the Efron handling of tied event times, `FALSE` for the
  Breslow handling.

- iter_max:

  The maximum number of iterations (20 in
  [`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)).

## Value

A numeric vector of length 6: the estimate of the log hazard ratio, its
standard error, the log partial likelihood at 0 and at the estimate, the
number of iterations, and a status code: 0 converged; 1 not converged
within `iter_max` iterations, in which case the estimate is the last
iterate; 2 the estimate is infinite or cannot be computed, in which case
the other elements except the number of iterations are `NA`.
