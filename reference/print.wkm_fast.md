# Print Method for wkm_fast Objects

Formats and prints a `wkm_fast` object in the same layout as the other
two-group summaries in the package. The header shows the control label
and the inference settings. The body shows the weighted integrated
survival difference (treatment minus control) with a confidence
interval, the test statistic, and the p-value.

## Usage

``` r
# S3 method for class 'wkm_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"wkm_fast"` returned by
  [`wkm_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`wkm_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md)

## Examples

``` r
set.seed(1)
n <- 200
g <- rep(0:1, each = n / 2)
tt <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.07))
cc <- rexp(n, 0.02)
time <- pmin(tt, cc)
event <- as.integer(tt <= cc)
print(wkm_fast(time, event, group = g, control = 0))
#> Weighted Kaplan-Meier test (Pepe-Fleming, two-group)
#> 
#>   control = 0
#>   weight = PF,  alternative = two.sided
#> 
#>                                             Est. lower 95% upper 95%     z
#> weighted difference (treatment - control) 2.5230    0.4581    4.5878 2.395
#>                                           Pr(>|z|)  
#> weighted difference (treatment - control)   0.0166 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
