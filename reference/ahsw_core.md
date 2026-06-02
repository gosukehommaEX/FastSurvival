# Core average hazard with survival weight computation for one group (C++ backend)

Internal C++ function that computes the average hazard with survival
weight (AHSW) of Uno and Horiguchi for a single group, together with the
building blocks needed for its asymptotic variance. The average hazard
on the window from 0 to `tau` is the ratio of the cumulative event
probability at `tau` to the restricted mean survival time at `tau`, both
based on the Kaplan-Meier estimate. The function walks a pooled sorted
dataset once, accumulating the Nelson-Aalen increment, the Kaplan-Meier
survival, the running restricted mean survival time and the at-risk
fraction at each distinct event time at or before `tau`, and forms the
two variance sums used for the ratio (log scale) and the difference
(identity scale). Not intended to be called directly by users; use
[`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md)
instead.

## Usage

``` r
ahsw_core(time_sorted, event_sorted, tau)
```

## Arguments

- time_sorted:

  A numeric vector of follow-up times sorted in ascending order.

- event_sorted:

  An integer vector of event indicators (1 = event, 0 = censored),
  aligned with `time_sorted`.

- tau:

  A numeric truncation time point.

## Value

A numeric vector of length 6:
`c(F_tau, RMST_tau, AH, v_Q, v_U, surv_tau)`, where F_tau is the
cumulative event probability at `tau`, RMST_tau is the restricted mean
survival time at `tau`, AH is the average hazard F_tau / RMST_tau, v_Q
is the variance sum on the log scale, v_U is the variance sum on the
identity scale, and surv_tau is the Kaplan-Meier survival at `tau`. The
variance sums must be divided by the group size by the caller.

## Details

Writing F for the cumulative event probability at `tau`, R for the
restricted mean survival time at `tau`, R(t) for the running restricted
mean survival time, dH(t) for the Nelson-Aalen increment and G(t) for
the at-risk fraction, the average hazard is F / R and the two variance
sums are the sum over event times of dH(t) (kernel(t))^2 / G(t). For the
average hazard on the log scale the kernel is 1 / F - R(t) / R, matching
the `v_Q` term of the reference implementation, and for the average
hazard on the identity scale the kernel is 1 / R - F R(t) / R^2,
matching the `v_U` term. The survival function is the right-continuous
Kaplan-Meier product, and the running restricted mean survival time
integrates the left-hand survival over each grid interval, exactly as in
the reference.
