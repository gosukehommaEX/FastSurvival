# Print Method for rmst_fast Objects

Formats and prints an `rmst_fast` object. In single-group mode it shows
the restricted mean survival time, its Greenwood standard error, and the
Wald confidence interval at the requested horizon. In two-group mode it
shows the per-group restricted mean survival times together with the
difference (treatment minus control) and ratio (treatment over control)
contrasts, each with a Wald z-statistic and two-sided p-value.

## Usage

``` r
# S3 method for class 'rmst_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"rmst_fast"` returned by
  [`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)

## Examples

``` r
set.seed(42)
t_raw <- rexp(100, rate = 1 / 10)
e_raw <- rbinom(100, 1, 0.7)

# Single-group
print(rmst_fast(t_raw, e_raw, tau = 10))
#> Restricted mean survival time (single-group)
#> 
#>   tau = 10
#> 
#>            RMST std.err lower 95% upper 95%
#> tau = 10 7.0341  0.3452    6.3575    7.7107

# Two-group comparison
set.seed(7)
n <- 200
time  <- c(rexp(n, 0.10), rexp(n, 0.07))
event <- rbinom(2 * n, 1, 0.8)
group <- rep(0:1, each = n)
print(rmst_fast(time, event, group = group, control = 0, tau = 10))
#> Restricted mean survival time (two-group)
#> 
#>   tau = 10,  control = 0
#> 
#>             RMST
#> control   6.7572
#> treatment 7.8384
#> 
#>                                    Est. lower 95% upper 95%     z Pr(>|z|)   
#> difference (treatment - control) 1.0813    0.4094    1.7531 3.154  0.00161 **
#> ratio (treatment / control)      1.1600    1.0572    1.2729 3.133  0.00173 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
