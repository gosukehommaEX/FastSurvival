# Fast Cox Proportional Hazards Estimation for Two Groups

Fits the Cox proportional hazards model with a single treatment
indicator to a two-group parallel trial and returns the maximum partial
likelihood estimate of the hazard ratio, the standard error of its
logarithm, and a Wald-type confidence interval, using output names
consistent with `summary(survival::coxph(...))`. The C++ backend forms
the risk sets in a single pass over the pooled sorted data and runs the
Newton-Raphson iteration on the numbers at risk and of events at the
distinct event times, so that an iteration takes time in proportion to
the number of distinct event times rather than the number of subjects.

## Usage

``` r
coxph_fast(
  time,
  event,
  group,
  control,
  side = 2,
  conf.level = 0.95,
  presorted = FALSE,
  strata = NULL,
  data = NULL,
  ties = c("efron", "breslow")
)
```

## Arguments

- time:

  A numeric vector of follow-up times for all subjects (pooled over both
  groups). Alternatively, a formula `Surv(time, event) ~ group`,
  optionally with
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) terms, with
  the variables in `data`; see
  [`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).

- event:

  An integer or numeric vector of event indicators (1 = event, 0 =
  censored), aligned with `time`.

- group:

  A vector of group labels aligned with `time`. Any type that supports
  equality comparison is accepted.

- control:

  A scalar value indicating which level of `group` represents the
  control group. Subjects with `group != control` are treated as the
  treatment group. With the formula interface it defaults to the first
  level of the group variable.

- side:

  1 for a one-sided test in the direction of treatment benefit (hazard
  ratio below 1, i.e. a negative coefficient) or 2 for a two-sided test
  (default 2). The returned vector has no p-value: `side` is stored as
  an attribute, and the [`print()`](https://rdrr.io/r/base/print.html)
  method reports the p-value that follows it. The confidence interval is
  always two-sided at `conf.level`.

- conf.level:

  A single numeric value in (0, 1) specifying the confidence level for
  the Wald interval. Defaults to 0.95.

- presorted:

  A logical value. If `TRUE`, `time`, `event`, and `group` are assumed
  to be already sorted in ascending order of `time` (with `strata`,
  sorted by stratum and by time within stratum, so that the rows of each
  stratum are contiguous), and the internal
  [`order()`](https://rdrr.io/r/base/order.html) call is skipped. If
  `FALSE` (default), sorting is handled internally. The order is
  checked, and an error is given when it does not hold.

- strata:

  An optional vector of stratum labels aligned with `time`, without
  missing values. When supplied, the stratified estimator described in
  Details is computed. Several stratification factors can be combined
  with [`interaction()`](https://rdrr.io/r/base/interaction.html).
  `NULL` (default) gives the unstratified estimator.

- data:

  A data frame in which the variables of a formula given as `time` are
  evaluated. Used only with the formula interface; see
  [`formula_interface`](https://gosukehommaEX.github.io/FastSurvival/reference/formula_interface.md).
  With the native pipe, `d |> coxph_fast(formula, ...)` uses the data
  frame `d` in the same way.

- ties:

  A character string naming the handling of tied event times, `"efron"`
  (default, as in
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)) or
  `"breslow"`. Without tied event times the two give the same result.

## Value

An object of class `"coxph_fast"`, which is a named numeric vector of
length 5 with elements matching the column names of
`summary(coxph(...))$coefficients` and `summary(coxph(...))$conf.int`:

- `coef`:

  Estimate of the log hazard ratio.

- `exp(coef)`:

  Estimate of the hazard ratio.

- `se(coef)`:

  Standard error of `coef`, from the observed information at the
  estimate.

- `lower .95`:

  Lower bound of the Wald confidence interval for the hazard ratio. The
  label reflects `conf.level` (e.g., `"lower .90"` when
  `conf.level = 0.90`).

- `upper .95`:

  Upper bound of the Wald confidence interval.

The attributes `ties`, `iter` (the number of Newton-Raphson iterations),
and `converged` record the fit. With `strata`, the number of strata is
stored in the attribute `strata`. The elements are `NA_real_` (still
with class `"coxph_fast"`) when the estimate cannot be computed (no
events, or an infinite estimate as described in Details), with
`converged = NA`, and when the iteration does not converge in 20
iterations, with `converged = FALSE` and a warning. The object also has
the class `"fastsurv"`, whose
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods are
described in
[`as.data.frame.fastsurv`](https://gosukehommaEX.github.io/FastSurvival/reference/as.data.frame.fastsurv.md).

## Details

Let t_k (k = 1, ..., K) denote the distinct observed event times in the
pooled sample. At each t_k, let n_T_k and n_C_k be the numbers at risk
in the treatment and control groups just before t_k (subjects censored
at t_k are at risk at t_k), and let O_T_k, O_C_k, and O_k = O_T_k +
O_C_k be the numbers of events in the treatment group, in the control
group, and in total at t_k. With the log hazard ratio beta and r =
exp(beta), the log partial likelihood of Cox (1972) with the Breslow
(1974) handling of tied event times is

l(beta) = sum_k \[O_T_k beta - O_k log(n_C_k + n_T_k r)\].

The Efron (1977) handling treats the O_k tied events at t_k as occurring
one after another and replaces the k-th term by O_T_k beta minus the
sum, over m = 0, ..., O_k - 1, of log(n_C_k + n_T_k r - (m / O_k)
(O_C_k + O_T_k r)). Without tied event times the two coincide. With
`strata`, the risk sets are formed within each stratum and l(beta) is
the sum over strata, so that the hazard ratio is common to all strata
and each stratum has its own baseline hazard.

The estimate maximizes l(beta) by Newton-Raphson iteration from beta =
0, which Therneau and Grambsch (2000, Appendix E) work through on a
small data set. A step after which l(beta) decreases is halved. From the
first halved step on, the iteration also keeps the interval between the
largest value at which the score l'(beta) has been found positive and
the smallest value at which it has been found negative, and a trial
value outside this interval is replaced by its midpoint; this keeps the
iteration convergent when the hazard ratio is extreme and one group is
small. The iteration stops when the relative change in l(beta) is at
most 1e-9, and at most 20 iterations are taken. The starting value, the
convergence criterion, and the maximum number of iterations are those of
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) with
the default
[`coxph.control()`](https://rdrr.io/pkg/survival/man/coxph.control.html),
so that the two estimates agree to rounding error when no step is
halved; when a step is halved, they can differ within the convergence
tolerance.

The standard error of the log hazard ratio is 1 / sqrt(I), where I =
-l”(beta) is the observed information at the estimate, and the Wald
confidence interval for the hazard ratio is
`exp(coef -/+ z * se(coef))`, where `z` is the standard normal quantile
for `conf.level`.

Let O_T be the total number of events in the treatment group, L the
number of events at times when no control subject is at risk, and H the
number of events at times when a treatment-group subject is at risk. The
score decreases in beta, from O_T - L as beta tends to -Inf to O_T - H
as beta tends to Inf, so the estimate is finite if and only if L \< O_T
\< H. Otherwise (for example, when all events are in one group) the
estimate is infinite and the result is `NA`, whereas
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
returns a large finite value, usually with a warning.

Times are tied only when they are exactly equal. By default
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) also
treats as tied the times that differ only by rounding error
(`timefix = TRUE` in
[`coxph.control()`](https://rdrr.io/pkg/survival/man/coxph.control.html));
the two give the same result when no distinct times are that close, and
`coxph.control(timefix = FALSE)` gives the comparable fit otherwise.

The returned object has class `"coxph_fast"` and is a named numeric
vector of length 5. A [`print()`](https://rdrr.io/r/base/print.html)
method formats the result similarly to `summary(coxph(...))`.

## References

Breslow, N. (1974). Covariance analysis of censored survival data.
*Biometrics*, *30*(1), 89-99.

Cox, D. R. (1972). Regression models and life-tables. *Journal of the
Royal Statistical Society. Series B (Methodological)*, *34*(2), 187-202.

Efron, B. (1977). The efficiency of Cox's likelihood function for
censored data. *Journal of the American Statistical Association*,
*72*(359), 557-565.

Therneau, T. M., & Grambsch, P. M. (2000). *Modeling Survival Data:
Extending the Cox Model*. New York: Springer.

## See also

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) for the Cox model
with any covariates.
[`print.coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.coxph_fast.md)
for the print method.

## Examples

``` r
library(survival)

# Compare coxph_fast with coxph on the ovarian dataset.
# coxph() treats rx as numeric with rx = 1 as the reference (control),
# so set control = 1 for a consistent comparison.
fit_fast <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
fit_fast
#> Cox proportional hazards model for two groups
#> 
#>   control = 1
#>   ties = efron
#>   alternative = two.sided
#> 
#> Coefficients:
#>          coef exp(coef) se(coef)      z Pr(>|z|)
#> group -0.5964    0.5508   0.5870 -1.016     0.31
#> 
#> Hazard ratio and 95% Wald confidence interval:
#>       exp(coef) exp(-coef) lower .95 upper .95
#> group    0.5508     1.8155    0.1743    1.7404

fit_cox <- coxph(Surv(futime, fustat) ~ rx, data = ovarian)
c(coxph_fast = unname(fit_fast["coef"]), coxph = unname(coef(fit_cox)))
#> coxph_fast      coxph 
#>   -0.59638   -0.59638 
c(coxph_fast = unname(fit_fast["se(coef)"]),
  coxph = sqrt(fit_cox$var[1, 1]))
#> coxph_fast      coxph 
#>  0.5869895  0.5869895 

# The Breslow handling of tied event times
coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1,
           ties = "breslow")
#> Cox proportional hazards model for two groups
#> 
#>   control = 1
#>   ties = breslow
#>   alternative = two.sided
#> 
#> Coefficients:
#>          coef exp(coef) se(coef)      z Pr(>|z|)
#> group -0.5964    0.5508   0.5870 -1.016     0.31
#> 
#> Hazard ratio and 95% Wald confidence interval:
#>       exp(coef) exp(-coef) lower .95 upper .95
#> group    0.5508     1.8155    0.1743    1.7404

# Stratified by residual disease, compared with coxph(... + strata())
coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1,
           strata = ovarian$resid.ds)
#> Stratified Cox proportional hazards model for two groups (2 strata)
#> 
#>   control = 1
#>   ties = efron
#>   alternative = two.sided
#> 
#> Coefficients:
#>          coef exp(coef) se(coef)      z Pr(>|z|)
#> group -0.6626    0.5155   0.5955 -1.113    0.266
#> 
#> Hazard ratio and 95% Wald confidence interval:
#>       exp(coef) exp(-coef) lower .95 upper .95
#> group    0.5155     1.9399    0.1605    1.6562
coef(coxph(Surv(futime, fustat) ~ rx + strata(resid.ds), data = ovarian))
#>        rx 
#> -0.662632 

# presorted = TRUE: sort once outside, reuse inside a loop
ord <- order(ovarian$futime)
coxph_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
           control = 1, presorted = TRUE)
#> Cox proportional hazards model for two groups
#> 
#>   control = 1
#>   ties = efron
#>   alternative = two.sided
#> 
#> Coefficients:
#>          coef exp(coef) se(coef)      z Pr(>|z|)
#> group -0.5964    0.5508   0.5870 -1.016     0.31
#> 
#> Hazard ratio and 95% Wald confidence interval:
#>       exp(coef) exp(-coef) lower .95 upper .95
#> group    0.5508     1.8155    0.1743    1.7404

# Formula interface, and the result as a data frame (log hazard ratio)
fit_f <- coxph_fast(Surv(futime, fustat) ~ rx, data = ovarian)
as.data.frame(fit_f)
#>   term estimate std.error  conf.low conf.high statistic   p.value
#> 1 coef -0.59638 0.5869895 -1.746858 0.5540983 -1.015998 0.3096304
exp(confint(fit_f))
#>          2.5 %   97.5 %
#> coef 0.1743207 1.740371

# The same with the native pipe
ovarian |> coxph_fast(Surv(futime, fustat) ~ rx) |> as.data.frame()
#>   term estimate std.error  conf.low conf.high statistic   p.value
#> 1 coef -0.59638 0.5869895 -1.746858 0.5540983 -1.015998 0.3096304

# \donttest{
# Speed comparison against coxph()
if (requireNamespace("microbenchmark", quietly = TRUE)) {
  microbenchmark::microbenchmark(
    coxph_fast = coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 1),
    coxph      = coxph(Surv(futime, fustat) ~ rx, data = ovarian),
    times = 1000
  )
}
#> Unit: microseconds
#>        expr      min       lq       mean   median       uq      max neval cld
#>  coxph_fast   74.970   86.070   98.00606  102.466  106.564  187.450  1000  a 
#>       coxph 1515.918 1551.214 1612.96558 1566.958 1588.970 7905.719  1000   b
# }
```
