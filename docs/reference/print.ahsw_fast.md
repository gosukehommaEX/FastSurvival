# Print Method for ahsw_fast Objects

Formats and prints an `ahsw_fast` object. The header shows the
truncation time and the control label. The body shows the per-group
average hazard with survival weight, followed by the between-group
contrasts: the ratio of average hazards (treatment over control) and the
difference of average hazards (treatment minus control), each with a
confidence interval and a two-sided p-value.

## Usage

``` r
# S3 method for class 'ahsw_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"ahsw_fast"` returned by
  [`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md)

## Examples

``` r
library(survival)
fit <- ahsw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
                 control = 1, tau = 600)
print(fit)
#> Average hazard with survival weight (two-group)
#> 
#>   tau = 600,  control = 1
#> 
#>               AH
#> control   0.0011
#> treatment 0.0008
#> 
#>                                        Est.  lower 95%  upper 95% p.value
#> ratio (treatment / control)       0.7548735  0.2622750  2.1726579   0.602
#> difference (treatment - control) -0.0002646 -0.0012914  0.0007622   0.613
```
