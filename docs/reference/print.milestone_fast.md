# Print Method for milestone_fast Objects

Formats and prints a `milestone_fast` object in the same layout as the
other two-group estimation summaries in the package. The header shows
the milestone timepoint, the control label, and the inference settings.
The body shows the per-group milestone survival with its confidence
interval, followed by the difference contrast (treatment minus control)
with a confidence interval, the test statistic, and the p-value. The
p-value follows the alternative recorded in the object.

## Usage

``` r
# S3 method for class 'milestone_fast'
print(x, digits = max(1L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"milestone_fast"` returned by
  [`milestone_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md).

- digits:

  Number of significant digits to display. Defaults to the global option
  `getOption("digits")`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`milestone_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md)

## Examples

``` r
set.seed(1)
time <- c(rexp(50, 0.1), rexp(50, 0.07))
status <- rep(1, 100)
group <- rep(c(0, 1), each = 50)
print(milestone_fast(time, status, group, control = 0, tau = 10, method = "loglog"))
#> Milestone survival (two-group)
#> 
#>   tau = 10,  control = 0
#>   method = loglog,  alternative = two.sided
#> 
#>           survival std.err lower 95% upper 95%
#> control       0.42  0.0698    0.2829    0.5510
#> treatment     0.56  0.0702    0.4124    0.6842
#> 
#>                                      Est. lower 95% upper 95%      z Pr(>|z|)
#> difference (treatment - control)  0.14000  -0.05738   0.32500 -1.395    0.163
```
