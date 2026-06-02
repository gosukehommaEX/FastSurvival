# Core Kaplan-Meier computation (C++ backend)

Internal C++ function that computes the Kaplan-Meier survival estimate
and the Greenwood variance sum at a specified time point. Binary search
locates the cutoff index; a single left-to-right scan accumulates the KM
product and Greenwood sum over event positions only. Accepts the event
vector as a numeric (double) vector to avoid an
[`as.integer()`](https://rdrr.io/r/base/integer.html) copy in R. Uses
`std::log1p(-1/n_risk)` for the KM log-product, which is more accurate
than `std::log(1 - 1/n_risk)` near n_risk = 1 and avoids one subtraction
per event. Not intended to be called directly by users; use
[`survfit_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)
instead.

## Usage

``` r
km_core(t_sorted, e_sorted, t_eval)
```

## Arguments

- t_sorted:

  A numeric vector of follow-up times sorted in ascending order.

- e_sorted:

  A numeric vector of event indicators (1 = event, 0 = censored),
  aligned with `t_sorted`. Accepted as numeric to avoid an
  integer-conversion copy on the R side.

- t_eval:

  A single numeric value specifying the evaluation time point.

## Value

A numeric vector of length 2: `c(surv, gw_sum)`, where `surv` is the KM
estimate and `gw_sum` is the Greenwood variance sum used to compute
SE\[S(t)\] = surv \* sqrt(gw_sum). Returns `c(1.0, 0.0)` when no events
occur up to `t_eval`.
