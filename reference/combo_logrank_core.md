# Core max-combo weighted log-rank computation on pooled sorted vectors (C++ backend)

Internal C++ function that computes a vector of Fleming-Harrington
weighted log-rank numerators and their full covariance matrix in a
single pass over a pooled sorted dataset. Each weight is a
Fleming-Harrington G(rho, gamma) weight evaluated from the
left-continuous pooled Kaplan-Meier estimate, with the pairs supplied
through `rho_vec` and `gamma_vec`. The variance increment at each event
time does not depend on the weight, so it is computed once and reused
for every numerator and every pair of weights. The diagonal of the
returned covariance matrix holds the per-scheme variances, and the
off-diagonal entries hold the between-scheme covariances used to form
the correlation matrix for the max-combo p-value. Tied event times are
processed atomically. Not intended to be called directly by users; use
the max-combo R wrapper instead.

## Usage

``` r
combo_logrank_core(time_sorted, event_sorted, j_sorted, rho_vec, gamma_vec)
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

- rho_vec:

  A numeric vector of Fleming-Harrington first parameters, one per
  weight scheme.

- gamma_vec:

  A numeric vector of Fleming-Harrington second parameters, one per
  weight scheme, aligned with `rho_vec`.

## Value

A list with three elements: `O1`, the observed number of events in the
treatment group; `U`, a numeric vector of weighted log-rank numerators
(one per weight scheme); and `V`, the numeric covariance matrix of the
numerators.

## Details

For a weight indexed by k the numerator is `U_k = sum w_k (O1 - E1)` and
the variance is `V_kk = sum w_k^2 v`, where `v` is the hypergeometric
variance increment `d n1 (n - n1) (n - d) / (n^2 (n - 1))`. The
between-scheme covariance is `V_ab = sum w_a w_b v`. Here `S_minus` is
the left-continuous pooled Kaplan-Meier estimate just prior to each
event time, initialized at 1, and the weight is
`S_minus^rho (1 - S_minus)^gamma`. The convention is the same as the
single-scheme weighted core: the numerator is negative when the
treatment group is favored.
