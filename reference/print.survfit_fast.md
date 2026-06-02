# Print Method for survfit_fast Objects

Formats and prints a `survfit_fast` object similarly to
`print(summary(survival::survfit(...), times = t_eval))`, showing the
Kaplan-Meier survival estimate, the Greenwood standard error on the
survival scale, and the confidence interval at the requested evaluation
time.

## Usage

``` r
# S3 method for class 'survfit_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"survfit_fast"` returned by
  [`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)

## Examples

``` r
set.seed(42)
t_raw <- rexp(100, rate = 1 / 10)
e_raw <- rbinom(100, 1, 0.7)
ord   <- order(t_raw)
fit   <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = 10)
print(fit)
#> Kaplan-Meier survival estimate (single time point)
#> 
#>        survival std.err lower 95% upper 95%
#> t = 10   0.4744  0.0543     0.379    0.5938
#> 
#>  Confidence interval type: log
```
