# Kaplan-Meier Curves as a Data Frame

Returns the Kaplan-Meier curves of a `kmcurve_fast` object as a data
frame, one row per group and distinct observed time, for drawing the
curves with other graphics software or for tabulation.

## Usage

``` r
# S3 method for class 'kmcurve_fast'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  conf.level = 0.95,
  conf.type = "log",
  ...
)
```

## Arguments

- x:

  A `kmcurve_fast` object from
  [`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md).

- row.names:

  `NULL` or a character vector of row names.

- optional:

  Ignored; the column names are fixed.

- conf.level:

  Confidence level of the pointwise limits. Defaults to 0.95.

- conf.type:

  Type of the pointwise limits: `"plain"`, `"log"` (default), or
  `"log-log"`, as in
  [`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md).

- ...:

  Further arguments, currently ignored.

## Value

A data frame with the columns `group` (the group label, as printed by
[`print.kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.kmcurve_fast.md)),
`time`, `n.risk`, `n.event`, `n.censor`, `surv`, `std.err`, `lower`, and
`upper`. In the two-group case the rows of the control group come first.

## Details

For each group the rows are the distinct observed times (event and
censoring times) in ascending order, preceded by a row at time 0 with
survival 1, so that a step function drawn through the rows starts at 1
and ends at the largest observed time. At each time, `n.risk` is the
number of subjects with an observed time at or after it, `n.event` and
`n.censor` are the numbers of events and censored times at it, and
`surv` is the Kaplan-Meier estimate just after it. `std.err` is the
Greenwood standard error on the survival scale, the `std.err` reported
by [`summary()`](https://rdrr.io/r/base/summary.html) of a
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
object, and `lower` and `upper` are the pointwise confidence limits of
the type `conf.type`, as drawn by
[`plot.kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/plot.kmcurve_fast.md).
At a time at which every subject at risk fails, the estimate is 0 and
`std.err` and the limits are 0, as in
[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md),
where [`summary()`](https://rdrr.io/r/base/summary.html) reports a `NaN`
standard error and `NA` limits. Apart from the row at time 0, the times,
counts, and survival estimates are those reported by
[`summary()`](https://rdrr.io/r/base/summary.html) of a
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
object with `censored = TRUE`.

## See also

[`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md),
[`plot.kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/plot.kmcurve_fast.md)

## Examples

``` r
set.seed(1)
n <- 50
t0 <- rexp(n, log(2) / 12)
t1 <- rexp(n, log(2) / 18)
cens <- runif(2 * n, 0, 30)
time  <- pmin(c(t0, t1), cens)
event <- as.integer(c(t0, t1) <= cens)
group <- rep(c(0, 1), each = n)

km <- as.data.frame(kmcurve_fast(time, event, group, control = 0))
head(km)
#>   group      time n.risk n.event n.censor      surv    std.err     lower
#> 1     0 0.0000000     50       0        0 1.0000000 0.00000000 1.0000000
#> 2     0 0.6452054     50       1        0 0.9800000 0.01979899 0.9419529
#> 3     0 1.0290310     49       1        0 0.9600000 0.02771281 0.9071919
#> 4     0 1.3938266     48       0        1 0.9600000 0.02771281 0.9071919
#> 5     0 1.8363654     47       1        0 0.9395745 0.03382296 0.8755673
#> 6     0 2.4201832     46       1        0 0.9191489 0.03876760 0.8462217
#>      upper
#> 1 1.000000
#> 2 1.000000
#> 3 1.000000
#> 4 1.000000
#> 5 1.000000
#> 6 0.998361
```
