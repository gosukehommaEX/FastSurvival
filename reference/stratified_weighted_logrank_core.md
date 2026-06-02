# Core stratified weighted log-rank computation (C++ backend)

Internal C++ function that computes the stratified weighted log-rank
totals from a dataset sorted by stratum first and by time within
stratum. Each stratum is processed as an independent weighted log-rank
test using the same logic as `weighted_logrank_core`, with weights
derived from the within-stratum pooled Kaplan-Meier estimate, and the
per-stratum numerator U and variance V are accumulated into the returned
totals. The overall statistic is Z = sum U / sqrt(sum V). This matches
the stratified weighted log-rank test of Magirr, which sums the
per-stratum U and V and standardizes once. Not intended to be called
directly by users; use
[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
with both `weight` and `strata` instead.

## Usage

``` r
stratified_weighted_logrank_core(
  time_sorted,
  event_sorted,
  j_sorted,
  strata_sorted,
  scheme,
  rho,
  gamma,
  t_star
)
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
  Rows of the same stratum must be contiguous.

- scheme:

  An integer weight-scheme code: 0 = Fleming-Harrington, 1 =
  modestly-weighted, 2 = Gehan-Breslow, 3 = Tarone-Ware.

- rho:

  A numeric Fleming-Harrington first parameter (scheme 0).

- gamma:

  A numeric Fleming-Harrington second parameter (scheme 0).

- t_star:

  A numeric timepoint for the modestly-weighted scheme (scheme 1).

## Value

A numeric vector of length 3: `c(O1, U, V)`, where O1 is the pooled
observed number of events in the treatment group, U is the summed
weighted numerator, and V is the summed weighted variance.

## Details

The scheme codes match `weighted_logrank_core`: 0 = Fleming-Harrington
G(rho, gamma); 1 = modestly-weighted with cap `1 / S(t_star)` computed
within each stratum; 2 = Gehan-Breslow; 3 = Tarone-Ware. The
left-continuous pooled Kaplan-Meier estimate is restarted at 1 at the
beginning of every stratum, so the weights of each stratum depend only
on that stratum's pooled data. For the modestly-weighted scheme the
weight cap is determined in a first pass within each stratum before
accumulation.
