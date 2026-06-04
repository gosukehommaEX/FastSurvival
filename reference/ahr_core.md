# Two-group average hazard ratio core (internal)

Internal Rcpp core for
[`ahr_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md).
Computes the Kalbfleisch-Prentice average hazard ratio over the interval
from 0 to `tau` for two groups, with the direct Greenwood-based variance
of the first group share `theta1`. The two groups' observations must
each be sorted in ascending time order and `grid` must be sorted and
unique.

## Usage

``` r
ahr_core(time1, status1, time2, status2, grid)
```

## Arguments

- time1:

  sorted right-censored times for the reference group

- status1:

  0/1 event indicators for the reference group

- time2:

  sorted right-censored times for the comparison group

- status2:

  0/1 event indicators for the comparison group

- grid:

  common evaluation times (sorted, unique, including 0 and tau)

## Value

A named numeric vector with elements theta1, theta2, ahr, var.theta1,
var.theta2, GL and valid.
