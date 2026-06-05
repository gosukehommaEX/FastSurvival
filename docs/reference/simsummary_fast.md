# Summarize Operating Characteristics from Sequential Analysis Output

Aggregates the per-simulation, per-look output of
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
into operating characteristics: the rejection rate, the
futility-stopping rate, the distribution of the stopping look, and the
expected analysis timing (events and calendar time at stopping). A
group-sequential design is summarized by applying the supplied per-look
boundaries in sequence, and a fixed design is the single-look (`K = 1`)
degenerate case of the same logic. This function consumes boundaries
computed elsewhere (for example by gsDesign or rpact) and does not
compute or spend alpha itself; it estimates the stopping probabilities
by Monte Carlo over the simulated trials rather than by the analytic
numerical integration used by those packages, so the two agree only up
to Monte Carlo error and converge as the number of simulations grows.

## Usage

``` r
simsummary_fast(
  data,
  eff.col = NULL,
  efficacy = NULL,
  fut.col = eff.col,
  futility = NULL,
  direction = c("lower", "upper"),
  p.col = NULL,
  alpha = NULL
)
```

## Arguments

- data:

  A data frame from
  [`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md),
  containing at least `sim` and `look`, the statistic columns named by
  the boundary arguments, and, when available, `n.event` and `cutoff`
  for the analysis-timing summaries.

- eff.col:

  A single character naming the statistic column compared with
  `efficacy`. Required with the Z mode.

- efficacy:

  A numeric vector of efficacy boundaries, one per look, on the scale of
  `eff.col`. Entries may be `NA` to omit the efficacy test at that look.
  Required with the Z mode.

- fut.col:

  A single character naming the statistic column compared with
  `futility`. Defaults to `eff.col`. Used only in the Z mode and only
  when `futility` is supplied; it may differ from `eff.col` so that
  futility is judged on a different scale (for example a log hazard
  ratio).

- futility:

  A numeric vector of futility boundaries, one per look, on the scale of
  `fut.col`, or `NULL` for no futility stopping. Entries may be `NA` to
  omit a futility boundary at that look. Used only in the Z mode.

- direction:

  A single string, either `"lower"` (default) or `"upper"`, giving the
  direction in which the boundaries are crossed. See Details.

- p.col:

  A single character naming the p-value column compared with `alpha`.
  Selects the p mode; use either the Z mode (`eff.col` with `efficacy`)
  or the p mode (`p.col` with `alpha`), not both.

- alpha:

  A numeric vector of per-look nominal significance levels, one per
  look. Entries may be `NA` to omit the test at a look. Required with
  `p.col`.

## Value

An object of class `"simsummary_fast"`: a data frame with one row per
population and look plus an `overall` summary row appended after each
population's looks. The columns are `population`, `look` (the look
index, or `"overall"` on the summary row), optionally `look.value`,
`n.enrolled.mean` and `n.event.mean` (the mean enrolled and event counts
at that look, or at the stopping look on the summary row),
`n.dropout.mean` and `n.pipeline.mean` (the mean dropout count and
pipeline count `n.enrolled - n.event - n.dropout`, when those columns
are present in `data`), `cutoff.mean` (the mean calendar time,
likewise), `prob.stop.efficacy`, `prob.stop.futility`, `prob.stop.any`,
and `cum.reject`. On the `overall` row `prob.stop.efficacy` is the total
rejection rate, `prob.stop.futility` the total futility rate,
`prob.stop.any` their sum, and `cum.reject` again the total rejection
rate; its timing columns are the expected counts and calendar time at
the stopping look. The number of simulations is stored in the attribute
`nsim` and the boundary settings in the attribute `boundary`.

## Details

The boundaries are supplied per look, one value for each distinct value
of the `look` column in `data`. Two boundary modes are available and
exactly one must be used.

In the Z mode the efficacy boundary is compared with the column named by
`eff.col` and the futility boundary with the column named by `fut.col`.
The two columns may differ, which lets the efficacy and futility rules
live on different scales. For example a beta-spending efficacy boundary
can be applied to the standardized statistic `logrank.z` while a
futility boundary expressed as a log hazard ratio is applied to
`cox.coef`. The crossing direction is set by `direction`: with `"lower"`
(the natural sign of `logrank.z`, `cox.z`, and `cox.coef`, where
treatment benefit is negative) the efficacy boundary is crossed when the
efficacy statistic is at or below `efficacy` and the futility boundary
is crossed when the futility statistic is at or above `futility`; with
`"upper"` the inequalities are reversed. Either boundary vector may
contain `NA` at some looks to omit that rule there, so efficacy-only and
futility-only looks are expressed by placing `NA` in the other vector. A
look with `NA` on both rules can never stop the trial.

In the p mode the p-value in the column named by `p.col` is compared
with the per-look nominal level `alpha`, rejecting when the p-value is
at or below the level. Futility is not used in the p mode, and `alpha`
may contain `NA` to omit the efficacy test at a look.

For each simulated trial the looks are examined in order. The trial
stops at the first look whose statistic crosses a boundary. Crossing the
efficacy boundary is a rejection of the null hypothesis; crossing the
futility boundary (Z mode only) is a stop without rejection. When both
are crossed at the same look the efficacy stop takes precedence. A look
whose relevant statistic is `NA`, or whose boundary is `NA`, triggers no
crossing of that rule and the trial continues. A trial that reaches the
final look without crossing the efficacy boundary does not reject.

Each look's `prob.stop.efficacy` is the marginal probability of stopping
for efficacy for the first time at that look, that is the probability of
not stopping at any earlier look and crossing the efficacy boundary at
this one. This is the stage-wise rejection contribution of a
group-sequential design; the cumulative sum `cum.reject` is the
cumulative power up to and including that look, and the total over all
looks is `rejection.rate`. These are the same quantities that gsDesign
and rpact report, estimated here by simulation.

The rejection rate is the type I error under a null data-generating
truth and the power under an alternative truth, but because the function
does not know the truth used to generate `data` it is reported neutrally
as the rejection rate and its interpretation is left to the user.

When `data` carries a `population` column (the long form produced by
[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
with `by.subgroup = TRUE`), the same boundaries are applied within each
population and the output has one block of rows per population.

## See also

[`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md),
[`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md),
[`print.simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/print.simsummary_fast.md).

## Examples

``` r
df <- simdata_fast(
  nsim     = 200,
  n        = c(150, 150),
  a.time   = c(0, 12),
  a.rate   = 300 / 12,
  e.hazard = list(0.05, 0.035),
  seed     = 1
)

res <- analysis_fast(df, control = 1, event.looks = c(60, 105, 150),
                     stat = c("logrank", "coxph"), side = 1)

# Efficacy on the standardized log-rank Z, futility on the log hazard ratio,
# with a futility-only first look and efficacy-only later looks
simsummary_fast(res,
                eff.col = "logrank.z",
                efficacy = c(NA, -2.96, -1.97),
                fut.col = "cox.coef",
                futility = c(log(1.2), NA, NA),
                direction = "lower")
#> Group-Sequential Operating Characteristics (simsummary_fast)
#>   Simulations: 200
#>   Boundaries: efficacy on 'logrank.z' (direction = lower), futility on 'cox.coef'
#> 
#> Stopping Boundaries: Look by Look
#>  Look Info. Frac. Events (s) Sample (n) Efficacy Z Futility Z Cum. Cross. Eff.
#>     1        0.40       60.0      286.1         NA     0.1823           0.0000
#>     2        0.70      105.0      300.0    -2.9600         NA           0.1300
#>     3        1.00      150.0      300.0    -1.9700         NA           0.5700
#> 
#> Events, Sample Size, Dropouts, Pipeline and Analysis Times: Look by Look
#>  Look Info. Frac. Sample (n) Events (s) Dropouts (d) Pipeline Analysis Time
#>     1        0.40      286.1       60.0          0.0    226.1         11.54
#>     2        0.70      300.0      105.0          0.0    195.0         16.50
#>     3        1.00      300.0      150.0          0.0    150.0         22.85
#>  Cross. Eff. Cross. Fut.
#>       0.0000      0.0100
#>       0.1300      0.0000
#>       0.4400      0.0000
#> 
#> Overall
#>   Rejection rate (efficacy):    0.5700
#>   Futility-stop rate:           0.0100
#>   Expected events at stop:      143.2
#>   Expected sample size at stop: 299.9
#>   Expected analysis time at stop:21.90

# p-value boundaries instead
simsummary_fast(res, p.col = "logrank.p",
                alpha = c(0.0006, 0.0151, 0.0245))
#> Group-Sequential Operating Characteristics (simsummary_fast)
#>   Simulations: 200
#>   Boundaries: nominal p-value on 'logrank.p'
#> 
#> Stopping Boundaries: Look by Look
#>  Look Info. Frac. Events (s) Sample (n) Nominal p Cum. Cross. Eff.
#>     1        0.40       60.0      286.1    0.0006           0.0550
#>     2        0.70      105.0      300.0    0.0151           0.3800
#>     3        1.00      150.0      300.0    0.0245           0.5900
#> 
#> Events, Sample Size, Dropouts, Pipeline and Analysis Times: Look by Look
#>  Look Info. Frac. Sample (n) Events (s) Dropouts (d) Pipeline Analysis Time
#>     1        0.40      286.1       60.0          0.0    226.1         11.54
#>     2        0.70      300.0      105.0          0.0    195.0         16.50
#>     3        1.00      300.0      150.0          0.0    150.0         22.85
#>  Cross. Eff.
#>       0.0550
#>       0.3250
#>       0.2100
#> 
#> Overall
#>   Rejection rate (efficacy):    0.5900
#>   Expected events at stop:      130.4
#>   Expected sample size at stop: 299.4
#>   Expected analysis time at stop:20.16
```
