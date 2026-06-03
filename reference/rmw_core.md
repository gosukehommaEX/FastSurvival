# Core robust modestly-weighted log-rank computation (C++ backend)

Internal C++ function that computes, in a single pass over a pooled
sorted dataset, the two component statistics of the robust
modestly-weighted (rMW) log-rank test of Magirr and Ohrn together with
their null covariance. The first component is the standard log-rank
statistic (weight one at every event time); the second is a
modestly-weighted log-rank statistic with weight
`min(1 / S(t-), 1 / s_star)`, where `S(t-)` is the left-continuous
pooled Kaplan-Meier estimate just prior to each event time and `s_star`
is a survival-probability threshold. The rMW test rejects for an extreme
value of the maximum of the two standardized components, so the
covariance of the two numerators under the null is required to recover
the joint distribution. Tied event times are processed atomically. Not
intended to be called directly by users; use
[`rmw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md)
instead.

## Usage

``` r
rmw_core(time_sorted, event_sorted, j_sorted, s_star)
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

- s_star:

  A numeric survival-probability threshold for the modestly-weighted
  component. The weight is capped at `1 / s_star`, so a value of `0.5`
  caps the weight at 2.

## Value

A numeric vector of length 6: `c(O1, U_lr, V_lr, U_mw, V_mw, C)`, where
O1 is the observed number of events in the treatment group, U_lr and
V_lr are the standard log-rank numerator and variance, U_mw and V_mw are
the modestly-weighted numerator and variance, and C is the null
covariance of the two numerators.

## Details

The standard log-rank numerator is `U_lr = sum (d1 - e1)` with variance
`V_lr = sum var_d`, where `d1` is the number of treatment-group events
at an event time, `e1 = d n1 / n` is the expected count, and
`var_d = d n1 n0 (n - d) / (n^2 (n - 1))` is the hypergeometric variance
of `d1`. The modestly-weighted numerator is `U_mw = sum w (d1 - e1)`
with variance `V_mw = sum w^2 var_d`. Because the standard log-rank
weight is one, the null covariance of the two numerators reduces to
`C = sum w var_d`. The modestly-weighted weight is capped at
`1 / s_star`, a constant, so unlike the timepoint parameterization used
by `weighted_logrank_core` no first pass is needed and all quantities
accumulate in one scan. When `s_star` is not positive the cap is treated
as infinite, giving the uncapped `1 / S(t-)` weight. Setting
`s_star = 1` caps the weight at one, so the modestly-weighted component
equals the standard log-rank component and the two are perfectly
correlated.
