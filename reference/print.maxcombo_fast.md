# Print Method for maxcombo_fast Objects

Formats and prints a `maxcombo_fast` object. It shows the total sample
size, a table of the Fleming-Harrington component Z-scores, and the
max-combo statistic with its p-value. For a one-sided test the statistic
is the most negative component, so a negative value favours the
treatment group, matching the sign convention of
[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md).
For a two-sided test the statistic is the largest component in absolute
value.

## Usage

``` r
# S3 method for class 'maxcombo_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"maxcombo_fast"` returned by
  [`maxcombo_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`maxcombo_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md)

## Examples

``` r
library(survival)
fit <- maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                     control = 1, side = 1)
print(fit)
#> Max-combo weighted log-rank test (two-group)
#> 
#>   N = 26
#> 
#>               Z
#> FH(0,0) -1.0309
#> FH(0,1)  0.0101
#> FH(1,0) -1.2980
#> FH(1,1) -0.0576
#> 
#>  Max-combo statistic = -1.298 (one-sided),  p-value = 0.1512
```
