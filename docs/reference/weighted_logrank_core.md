# Core weighted log-rank computation on pooled sorted vectors (C++ backend)

Internal C++ function that computes the weighted log-rank numerator U
and variance V directly from a pooled sorted dataset and a group
indicator. Four weight schemes are supported, selected by `scheme`:
Fleming- Harrington G(rho, gamma), the modestly-weighted log-rank test,
Gehan- Breslow, and Tarone-Ware. The weights are functions of the
left-continuous pooled Kaplan-Meier estimate (for Fleming-Harrington and
the modestly- weighted test) or of the at-risk count (for Gehan-Breslow
and Tarone-Ware). Tied event times are processed atomically. Not
intended to be called directly by users; use
[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
with the `weight` argument instead.

## Usage

``` r
weighted_logrank_core(
  time_sorted,
  event_sorted,
  j_sorted,
  scheme,
  rho,
  gamma,
  t_star
)
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

- scheme:

  An integer weight-scheme code: 0 = Fleming-Harrington, 1 =
  modestly-weighted, 2 = Gehan-Breslow, 3 = Tarone-Ware.

- rho:

  A numeric Fleming-Harrington first parameter (used when `scheme = 0`).

- gamma:

  A numeric Fleming-Harrington second parameter (used when
  `scheme = 0`).

- t_star:

  A numeric timepoint for the modestly-weighted scheme (used when
  `scheme = 1`).

## Value

A numeric vector of length 3: `c(O1, U, V)`, where O1 is the observed
number of events in the treatment group, U is the weighted log-rank
numerator sum w (O1 - E1), and V is the weighted variance.

## Details

The scheme codes are: 0 = Fleming-Harrington G(rho, gamma) with weight
`S_minus^rho * (1 - S_minus)^gamma`; 1 = modestly-weighted log-rank with
weight `min(1 / S_minus, max_weight)`, where `max_weight` is the
reciprocal of the smallest right-continuous pooled Kaplan-Meier value at
or after `t_star` (and is 1 when `t_star = 0`); 2 = Gehan-Breslow with
weight `n_j`; 3 = Tarone-Ware with weight `sqrt(n_j)`. Here `S_minus` is
the left-continuous pooled Kaplan-Meier estimate just prior to each
event time, initialized at 1. The modestly-weighted scheme requires the
largest weight to be known before accumulation, so its `max_weight` is
obtained in a first pass over the event times; the other schemes
accumulate in a single pass.
