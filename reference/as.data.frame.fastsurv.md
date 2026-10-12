# Data Frame, Estimates, and Confidence Limits of a Result

The results of the estimation and testing functions have, besides their
own class, the class `"fastsurv"`. Its
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
returns every estimate and test of a result as a row of a data frame
with the same columns for all the functions, so that results can be
stacked with [`rbind()`](https://rdrr.io/r/base/cbind.html).
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the estimates and
[`confint()`](https://rdrr.io/r/stats/confint.html) their confidence
limits.

## Usage

``` r
# S3 method for class 'fastsurv'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'fastsurv'
coef(object, ...)

# S3 method for class 'fastsurv'
confint(object, parm, level = NULL, ...)
```

## Arguments

- x, object:

  A result of an estimation or testing function, of class `"fastsurv"`.

- row.names:

  `NULL` or a character vector of row names.

- optional:

  Ignored; the column names are fixed.

- ...:

  Further arguments, currently ignored.

- parm:

  The estimates for which limits are returned, as terms or as positions
  among the estimates. All estimates by default.

- level:

  `NULL` (default) or the confidence level at which the result was
  computed.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with the columns described in Details.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector, and [`confint()`](https://rdrr.io/r/stats/confint.html) a matrix
with one row per estimate and two columns labeled by the percentages of
the limits, as
[`stats::confint()`](https://rdrr.io/r/stats/confint.html).

## Details

The data frame has the columns `term`, `estimate`, `std.error`,
`conf.low`, `conf.high`, `statistic`, and `p.value`, with one row for
each quantity in the table below. The terms are the names of the
elements of the result. The estimate, its standard error, and its
confidence limits are on the same scale. The statistic and the p-value
are those of the test that the function reports for that quantity,
following its `side`, as printed by the
[`print()`](https://rdrr.io/r/base/print.html) method. Entries that the
function does not compute are `NA`.

|  |  |
|----|----|
| [`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md) | `surv` |
| [`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md) | `logrank`, or the name of the weight (test only) |
| [`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md) | `coef` (the log hazard ratio) |
| [`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md) | `rmst`; or `rmst.ctrl`, `rmst.trt`, `diff`, `log.ratio` |
| [`wmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md) | `wmst`; or `wmst.control`, `wmst.treatment`, `diff` |
| [`medsurv_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md) | `median`; or `median.control`, `median.treatment`, `diff` |
| [`milestone_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md) | `surv.control`, `surv.treatment`, `diff` |
| [`maxcombo_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md) | `maxcombo` and the component statistics, such as `FH(0,1)` (test only) |
| [`rmw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md) | `rmw`, `logrank`, `mwlrt` (test only) |
| [`wkm_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md) | `wdiff` |
| [`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md) | `ah.ctrl`, `ah.trt`, `log.rah`, `dah` |
| [`ahr_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md) | `log.ahr`, `theta` |

Ratios are reported on the log scale, on which their confidence
intervals are symmetric and their tests are computed: `coef` of
[`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
is the log hazard ratio (as from
[`coef()`](https://rdrr.io/r/stats/coef.html) of
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)), and
`log.ratio`, `log.rah`, and `log.ahr` are the logarithms of the ratio of
restricted mean survival times, the ratio of average hazards, and the
average hazard ratio. Take [`exp()`](https://rdrr.io/r/base/Log.html) of
the estimate and of the confidence limits for the ratio.
[`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
and
[`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md)
do not store the standard error of a log ratio or of the difference of
average hazards; it is recovered from the confidence limits, which are
symmetric on these scales, as `(conf.high - conf.low) / (2 q)`, where
`q` is the normal quantile of the confidence level, and the statistic of
`log.rah` and `dah` is the estimate divided by this standard error.

The rows of a test only carry the statistic: the value of
[`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
(the Z-score when `side = 1`, the chi-square when `side = 2`), the
max-combo and rMW statistics with their components, which have no
p-value of their own. For
[`milestone_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md),
the statistic of `diff` is that of `method`, which for `"loglog"` is on
the complementary log-log scale; a standard error is given only for
`method = "wald"`, whose interval is symmetric. For
[`ahr_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md),
`theta` is the treatment group's share of the total hazard with its
standard error and the primary test of
[`ahr_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md),
and `log.ahr` carries the equivalent test on the log scale.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the `estimate`
column named by `term`, and
[`confint()`](https://rdrr.io/r/stats/confint.html) the confidence
limits as a matrix, for the rows that are estimates. Both give an error
for a test only. The limits are those computed by the function at its
`conf.level`; they cannot be recomputed at another level, so `level`
must be `NULL` or equal to that `conf.level`.

## See also

[`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).

## Examples

``` r
library(survival)

fit <- coxph_fast(Surv(futime, fustat) ~ rx, data = ovarian)
as.data.frame(fit)
#>   term estimate std.error  conf.low conf.high statistic   p.value
#> 1 coef -0.59638 0.5868442 -1.746574 0.5538135 -1.016249 0.3095107
coef(fit)
#>     coef 
#> -0.59638 
exp(confint(fit))
#>          2.5 %   97.5 %
#> coef 0.1743704 1.739875

as.data.frame(rmst_fast(Surv(futime, fustat) ~ rx, data = ovarian,
                        tau = 600))
#>        term    estimate  std.error     conf.low   conf.high statistic
#> 1 rmst.ctrl 427.5384615         NA           NA          NA        NA
#> 2  rmst.trt 534.9059829         NA           NA          NA        NA
#> 3      diff 107.3675214 62.5917172 -15.30999001 230.0450327  1.715363
#> 4 log.ratio   0.2240467  0.1418789  -0.05403079   0.5021243  1.579141
#>      p.value
#> 1         NA
#> 2         NA
#> 3 0.08627868
#> 4 0.11430378

# Results of a loop (with the vector interface) stacked into a data frame
res <- lapply(1:3, function(s) {
  set.seed(s)
  i <- sample(nrow(ovarian), replace = TRUE)
  fit_s <- coxph_fast(ovarian$futime[i], ovarian$fustat[i], ovarian$rx[i],
                      control = 1)
  cbind(sim = s, as.data.frame(fit_s))
})
do.call(rbind, res)
#>   sim term   estimate std.error  conf.low conf.high  statistic    p.value
#> 1   1 coef -0.9047674 0.5449332 -1.972817 0.1632820 -1.6603272 0.09684865
#> 2   2 coef -0.1123616 0.6784779 -1.442154 1.2174305 -0.1656084 0.86846515
#> 3   3 coef -0.5031459 0.5086887 -1.500157 0.4938656 -0.9891037 0.32261240
```
