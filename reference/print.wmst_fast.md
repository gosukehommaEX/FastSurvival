# Print Method for wmst_fast Objects

Formats and prints a `wmst_fast` object in the same layout as the other
two-group summaries in the package. The header shows the window, the
control label and the inference settings. The body shows the per-group
window mean survival time with its confidence interval, followed, for a
two-group object, by the difference contrast (treatment minus control)
with a confidence interval, the test statistic, and the p-value.

## Usage

``` r
# S3 method for class 'wmst_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"wmst_fast"` returned by
  [`wmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`wmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md)

## Examples

``` r
set.seed(1)
n <- 200
g <- rep(0:1, each = n / 2)
tt <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.07))
cc <- rexp(n, 0.02)
time <- pmin(tt, cc)
event <- as.integer(tt <= cc)
print(wmst_fast(time, event, group = g, control = 0, tau1 = 2, tau2 = 12))
#> Window mean survival time (two-group)
#> 
#>   window = [2, 12],  control = 0
#>   alternative = two.sided
#> 
#>             wmst std.err lower 95% upper 95%
#> control   5.8371   0.402    5.0491    6.6251
#> treatment 6.6057   0.402    5.8177    7.3936
#> 
#>                                     Est. lower 95% upper 95%     z Pr(>|z|)
#> difference (treatment - control)  0.7686   -0.3458    1.8829 1.352    0.176
```
