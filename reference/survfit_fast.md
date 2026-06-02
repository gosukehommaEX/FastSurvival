# Fast Kaplan-Meier Survival Probability at a Specified Time Point

Computes the Kaplan-Meier survival probability at a specified time
point, together with a standard error and confidence interval based on
Greenwood's variance formula. The C++ backend performs binary search for
the evaluation cutoff and accumulates the Kaplan-Meier product and
Greenwood sum in a single scan over event positions only, without
constructing intermediate vectors.

## Usage

``` r
survfit_fast(
  t_sorted,
  e_sorted,
  t_eval,
  conf.int = 0.95,
  conf.type = "log",
  presorted = TRUE
)
```

## Arguments

- t_sorted:

  A numeric vector of event or censoring times. Must be sorted in
  ascending order when `presorted = TRUE`.

- e_sorted:

  An integer or numeric vector of event indicators (1 = event, 0 =
  censored), aligned with `t_sorted`.

- t_eval:

  A single numeric value specifying the time point at which the survival
  probability is evaluated.

- conf.int:

  A single numeric value in (0, 1) specifying the confidence level.
  Defaults to 0.95.

- conf.type:

  A character string specifying the confidence interval type. Must be
  one of `"plain"`, `"log"`, or `"log-log"`. Defaults to `"log"`.

- presorted:

  A logical value. If `TRUE` (default), `t_sorted` and `e_sorted` are
  assumed to be sorted in ascending order of time. If `FALSE`, the
  vectors are sorted internally before computation.

## Value

An object of class `"survfit_fast"`, which is a named numeric vector of
length 4 with elements `surv`, `std.err`, `lower`, and `upper`,
representing the Kaplan-Meier survival estimate, the Greenwood standard
error SE\[S(t)\], and the lower and upper confidence limits at `t_eval`.
The evaluation time, confidence level, and confidence interval type are
stored as attributes `t_eval`, `conf.int`, and `conf.type`. Returns a
vector of `NA_real_` values (still with class `"survfit_fast"`) when `n`
is zero.

## Details

The Kaplan-Meier estimate at time `t_eval` is defined as the
product-limit estimator evaluated at the largest observed event time
less than or equal to `t_eval`. If `t_eval` is smaller than the first
observed event time, `S(t) = 1` and the standard error is zero.

The standard error is estimated by Greenwood's formula:

SE\[S(t)\] = S(t) \* sqrt(sum\_{t_i \<= t, d_i \> 0} d_i / (n_i \*
(n_i - d_i)))

where d_i is the number of events and n_i is the number at risk at time
t_i. The output field `std.err` follows the convention of
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html), which
reports SE\[S(t)\] / S(t) (i.e., the standard error on the log scale)
when `conf.type != "plain"`, and SE\[S(t)\] when `conf.type = "plain"`.
This function always returns SE\[S(t)\] (the standard error on the
survival scale).

When `S(t_eval) = 0` (all subjects have experienced the event by
`t_eval`), the standard error is zero and the confidence interval
collapses to `[0, 0]`, consistent with
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html).

When `presorted = TRUE` (default), `t_sorted` and `e_sorted` are assumed
to be sorted in ascending order of time. When `presorted = FALSE`, the
vectors are sorted internally before computation.

Three confidence interval types are supported via `conf.type`:

- `"plain"`: Linear interval on the survival scale, S(t) +/- z \* SE.
  The bounds are clipped to \[0, 1\].

- `"log"`: Interval on the log scale (default in
  [`survfit`](https://rdrr.io/pkg/survival/man/survfit.html)), S(t) \*
  exp(+/- z \* SE / S(t)).

- `"log-log"`: Interval on the complementary log-log scale, S(t)^exp(+/-
  z \* SE / (S(t) \* log(S(t)))).

The returned object has class `"survfit_fast"` and is a named numeric
vector of length 4 with the evaluation time `t_eval`, the confidence
level `conf.int`, and the confidence interval type `conf.type` stored as
attributes. A [`print()`](https://rdrr.io/r/base/print.html) method
formats the result similarly to
`print(summary(survival::survfit(...)))`.

## References

Kaplan, E. L., & Meier, P. (1958). Nonparametric Estimation from
Incomplete Observations. Journal of the American Statistical
Association, 53(282), 457–481.

## See also

[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) for the
standard Kaplan-Meier estimator.
[`print.survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.survfit_fast.md)
for the print method.

## Examples

``` r
set.seed(42)
t_raw <- rexp(100, rate = 1 / 10)
e_raw <- rbinom(100, 1, 0.7)

# presorted = TRUE (default): sort once outside, reuse inside a loop
ord <- order(t_raw)
t_s <- t_raw[ord]
e_s <- e_raw[ord]
survfit_fast(t_s, e_s, t_eval = 10, conf.type = "plain")
#> Kaplan-Meier survival estimate (single time point)
#> 
#>        survival std.err lower 95% upper 95%
#> t = 10   0.4744  0.0543    0.3679    0.5809
#> 
#>  Confidence interval type: plain
survfit_fast(t_s, e_s, t_eval = 10, conf.type = "log")
#> Kaplan-Meier survival estimate (single time point)
#> 
#>        survival std.err lower 95% upper 95%
#> t = 10   0.4744  0.0543     0.379    0.5938
#> 
#>  Confidence interval type: log
survfit_fast(t_s, e_s, t_eval = 10, conf.type = "log-log")
#> Kaplan-Meier survival estimate (single time point)
#> 
#>        survival std.err lower 95% upper 95%
#> t = 10   0.4744  0.0543    0.3651    0.5759
#> 
#>  Confidence interval type: log-log

# presorted = FALSE: sort internally, convenient for one-off calls
survfit_fast(t_raw, e_raw, t_eval = 10, presorted = FALSE)
#> Kaplan-Meier survival estimate (single time point)
#> 
#>        survival std.err lower 95% upper 95%
#> t = 10   0.4744  0.0543     0.379    0.5938
#> 
#>  Confidence interval type: log

# Validation against survival::survfit
library(survival)
fit <- survfit(Surv(t_raw, e_raw) ~ 1, conf.type = "plain")
summary(fit, times = 10)
#> Call: survfit(formula = Surv(t_raw, e_raw) ~ 1, conf.type = "plain")
#> 
#>  time n.risk n.event survival std.err lower 95% CI upper 95% CI
#>    10     37      46    0.474  0.0543        0.368        0.581
```
