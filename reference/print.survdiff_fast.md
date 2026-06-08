# Print Method for survdiff_fast Objects

Formats and prints a `survdiff_fast` object similarly to
`print(survival::survdiff(...))`, showing the observed and expected
event counts for the control and treatment groups, the per-group
contributions `(O-E)^2 / E` and `(O-E)^2 / V`, the test statistic, and
the corresponding p-value. For a weighted log-rank test the header names
the weight scheme (and notes stratification), and only the observed
event counts are shown alongside the weighted statistic, since a single
unweighted expected count is not defined for a weighted test.

## Usage

``` r
# S3 method for class 'survdiff_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"survdiff_fast"` returned by
  [`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)

## Examples

``` r
library(survival)
fit <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                     control = 1, side = 2)
print(fit)
#> Log-rank test (two-group)
#> 
#>   N = 26,  control = 1
#> 
#>           Observed Expected (O-E)^2/E (O-E)^2/V
#> control          7   5.2335    0.5962    1.0627
#> treatment        5   6.7665    0.4612    1.0627
#> 
#>  Chi-square = 1.063 on 1 df,  p-value = 0.3026  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
