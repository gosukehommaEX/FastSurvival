# Core Kaplan-Meier computation (C++ backend)

Internal C++ function that computes the Kaplan-Meier survival estimate
and the Greenwood variance sum at a specified time point. Binary search
locates the cutoff index; a single left-to-right scan groups tied times
by distinct time, applying one Kaplan-Meier factor and one Greenwood
increment per distinct event time. Grouping the ties follows the
risk-set convention of
[`survival::survfit`](https://rdrr.io/pkg/survival/man/survfit.html) and
makes the result independent of the order of events and censorings
within a tied time. Accepts the event vector as a numeric (double)
vector to avoid an [`as.integer()`](https://rdrr.io/r/base/integer.html)
copy in R. Uses `std::log1p(-d/n_risk)` for the KM log-product, which is
more accurate than `std::log(1 - d/n_risk)` near a full drop and avoids
one subtraction per distinct event time. Not intended to be called
directly by users; use
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
