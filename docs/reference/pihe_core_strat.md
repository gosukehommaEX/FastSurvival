# Core stratified PiHE hazard ratio computation (C++ backend)

Internal C++ function that computes the quantities of the stratified
Pike-Halley Estimator. The at-risk sets are restarted in every stratum
and the observed and expected event totals, the score, the information,
and the curvature term are summed over strata, which gives the Pike
anchor and its Halley correction for the stratified Breslow partial
likelihood. Not intended to be called directly by users; use
[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
with `strata` instead.

## Usage

``` r
pihe_core_strat(time_sorted, event_sorted, j_sorted, strata_sorted)
```

## Arguments

- time_sorted:

  A numeric vector of pooled follow-up times, sorted by stratum and in
  ascending order of time within stratum.

- event_sorted:

  An integer vector of event indicators (1 = event, 0 = censored),
  aligned with `time_sorted`.

- j_sorted:

  An integer vector of group indicators (1 = treatment, 0 = control),
  aligned with `time_sorted`.

- strata_sorted:

  An integer vector of stratum codes aligned with `time_sorted`; rows of
  the same stratum must be contiguous.

## Value

A numeric vector of length 4: `c(theta_0, U_0, I_0, J_0)`. Returns a
length-4 vector of `NA_real_` when the estimate cannot be computed.
