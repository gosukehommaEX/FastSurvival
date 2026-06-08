# Print Method for medsurv_fast Objects

Formats and prints a `medsurv_fast` object in the same layout as the
other two-group estimation summaries in the package. The header shows
the control label and the inference settings. The body shows the
per-group median survival with its confidence interval, followed, for a
two-group object, by the difference contrast (treatment minus control)
with a confidence interval, the test statistic, and the p-value.

## Usage

``` r
# S3 method for class 'medsurv_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"medsurv_fast"` returned by
  [`medsurv_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`medsurv_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md)

## Examples

``` r
set.seed(1)
time <- c(rexp(50, 0.1), rexp(50, 0.07))
status <- rep(1, 100)
group <- rep(c(0, 1), each = 50)
print(medsurv_fast(time, status, group, control = 0))
#> Median survival time (two-group)
#> 
#>   control = 0
#>   method = km,  alternative = two.sided
#> 
#>            median std.err lower 95% upper 95%
#> control    7.5518  1.3257    5.3533   10.6532
#> treatment 11.6195  1.9604    8.3480   16.1732
#> 
#>                                     Est. lower 95% upper 95%     z Pr(>|z|)  
#> difference (treatment - control)  4.0677   -0.5706    8.7061 1.719   0.0856 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
