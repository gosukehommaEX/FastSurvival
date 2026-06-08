# Print Method for coxph_fast Objects

Formats and prints a `coxph_fast` object similarly to
`summary(survival::coxph(...))`, showing the point estimate of the log
hazard ratio, the hazard ratio, the standard error on the log scale, the
Wald z-statistic, the corresponding two-sided p-value, and the Wald
confidence interval for the hazard ratio.

## Usage

``` r
# S3 method for class 'coxph_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"coxph_fast"` returned by
  [`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments passed to
  [`format`](https://rdrr.io/r/base/format.html) (currently unused).

## Value

Invisibly returns `x`.

## See also

[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)

## Examples

``` r
library(survival)
fit <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
print(fit)
#> Pike-Halley estimator for the hazard ratio (two-group)
#> 
#>   control = 1
#>   alternative = two.sided
#> 
#> Coefficients:
#>          coef exp(coef) se(coef)      z Pr(>|z|)
#> group -0.5964    0.5508   0.5868 -1.016     0.31
#> 
#> Hazard ratio and 95% Wald confidence interval:
#>       exp(coef) exp(-coef) lower .95 upper .95
#> group    0.5508     1.8155    0.1744    1.7399
```
