# Core PiHE hazard ratio computation (C++ backend)

Internal C++ function that computes the Pike-Halley Estimator (PiHE) for
the hazard ratio. Accepts pooled sorted vectors plus an integer group
indicator, performs group splitting and the two-pointer merge scan
entirely in C++, and returns the quantities needed for the Halley
correction and Wald interval. Not intended to be called directly by
users; use
[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
instead.

Implementation notes: Pass 1 accumulates pooled scalars (O_T, O_C, E_T,
E_C) and stores per-distinct-event-time summaries in a single struct
array (better cache locality than four parallel vectors) with a single
`reserve(n)` call (no event-count prepass over the input). Pass 2 walks
the saved summaries once at the Pike anchor theta_0 to compute U_0, I_0,
and J_0.

## Usage

``` r
pihe_core(time_sorted, event_sorted, j_sorted)
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

A numeric vector of length 4: `c(theta_0, U_0, I_0, J_0)`. Returns a
length-4 vector of `NA_real_` when the estimate cannot be computed.
