# Print Method for rmw_fast Objects

Formats and prints an `rmw_fast` object, showing the standardized
log-rank and modestly-weighted component Z-scores, their null
correlation, the survival-probability threshold `s_star`, the combined
test statistic, and the corresponding p-value.

## Usage

``` r
# S3 method for class 'rmw_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"rmw_fast"` returned by
  [`rmw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`rmw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md)

## Examples

``` r
library(survival)
fit <- rmw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                control = 1, side = 1, s_star = 0.5)
print(fit)
#> Robust modestly-weighted log-rank test (two-group)
#> 
#>   N = 26,  s_star = 0.5
#> 
#>                         Z
#> log-rank          -1.0309
#> modestly-weighted -0.7583
#> 
#>   Null correlation = 0.9821
#>  min Z = -1.031,  one-sided p-value = 0.169
```
