# Formula Interface of the Estimation and Testing Functions

The estimation and testing functions and
[`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md)
take vectors of times, event indicators, and group labels, which suits
repeated calls inside a simulation loop. For the analysis of a single
data set they also accept a model formula with a `Surv(time, event)`
response and a data frame, as the functions of the survival package do:


    coxph_fast(Surv(time, status) ~ arm, data = d)
    survdiff_fast(Surv(time, status) ~ arm + strata(region), data = d)
    rmst_fast(Surv(time, status) ~ 1, data = d, tau = 12)

The formula is given as the first argument (by position, or by the name
of that argument, such as `time`; there is no argument named `formula`),
and the data frame as `data` or as the second argument. The other
arguments (`side`, `tau`, `conf.level`, and so on) are given by name as
with the vector interface, and the result is the same as that of the
vector interface with the vectors evaluated from the formula.

## Details

The response must be `Surv(time, event)` with right-censored data, where
`event` is 1 (or `TRUE`) for an event and 0 (or `FALSE`) for a censored
time; an expression such as `Surv(time, status == 2)` may be used. A
status coded 1 (censored) and 2 (event), as in data sets of the survival
package, gives an error that suggests this expression. The call to
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) is read from the
formula and not evaluated, so the survival package need not be attached.
Other forms of [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html),
such as counting-process data, are not supported.

The right-hand side is `1` for a single group, or a single group
variable with exactly two distinct values for a two-group comparison.
Which of the two is allowed follows the vector interface of each
function:
[`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md)
takes only `~ 1`;
[`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md),
[`wmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md),
[`medsurv_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md),
and
[`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md)
take `~ 1` or a group variable; the other functions require a group
variable.
[`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
and
[`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md)
also accept [`strata()`](https://rdrr.io/pkg/survival/man/strata.html)
terms, such as `arm + strata(region)` or `arm + strata(region, sex)`;
several stratification variables define a stratum for each combination
of their values. Interactions, further covariates, and other formula
operators are not supported, because the functions compare two groups
without adjustment.

With a group variable, `control` defaults to its first level: the first
level of a factor that occurs in the data, and otherwise the smallest
value. This is the reference level of a model formula with the group
variable as a factor. `control` can be given to choose the other group
as the control group. A numeric group variable is treated as two groups,
whereas
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
treats it as a continuous covariate, so the coefficient of
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
equals the log hazard ratio of the two groups only when the two values
differ by 1.

The variables are evaluated in `data`, and variables not found there in
the environment of the formula. Missing values in the variables of the
formula are an error, as in the vector interface; remove the incomplete
rows first, for example with
[`na.omit()`](https://rdrr.io/r/stats/na.fail.html).

The formula interface adds the cost of evaluating the formula to every
call, so the vector interface, with `presorted = TRUE` where available,
remains the choice inside simulation loops.

## See also

[`as.data.frame.fastsurv`](https://gosukehommaEX.github.io/FastSurvival/reference/as.data.frame.fastsurv.md)
for the results as data frames.

## Examples

``` r
library(survival)

# The same hazard ratio from the formula and from the vectors
coxph_fast(Surv(futime, fustat) ~ rx, data = ovarian)
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
coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
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

# Stratified log-rank test
survdiff_fast(Surv(futime, fustat) ~ rx + strata(resid.ds), data = ovarian)
#> Stratified log-rank test (two-group, 2 strata)
#> 
#>   N = 26,  control = 1
#> 
#>           Observed Expected (O-E)^2/E (O-E)^2/V
#> control          7   5.0952    0.7121    1.2796
#> treatment        5   6.9048    0.5255    1.2796
#> 
#>  Chi-square = 1.28 on 1 df,  p-value = 0.258  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Single-group and two-group restricted mean survival time
rmst_fast(Surv(futime, fustat) ~ 1, data = ovarian, tau = 600)
#> Restricted mean survival time (single-group)
#> 
#>   tau = 600
#> 
#>               RMST std.err lower 95% upper 95%
#> tau = 600 480.6607 33.0567  415.8707  545.4507
rmst_fast(Surv(futime, fustat) ~ rx, data = ovarian, tau = 600, control = 2)
#> Restricted mean survival time (two-group)
#> 
#>   tau = 600,  control = 2
#>   alternative = two.sided
#> 
#>               RMST
#> control   534.9060
#> treatment 427.5385
#> 
#>                                       Est. lower 95% upper 95%      z Pr(>|z|)
#> difference (treatment - control) -107.3675 -230.0450   15.3100 -1.715   0.0863
#> ratio (treatment / control)         0.7993    0.6052    1.0555 -1.579   0.1143
#>                                   
#> difference (treatment - control) .
#> ratio (treatment / control)       
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
