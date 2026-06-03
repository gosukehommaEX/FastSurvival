# Changelog

## FastSurvival 0.2.0

- New estimation and testing functions:
  - [`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md):
    restricted mean survival time for a single group or a two-group
    comparison (difference and ratio contrasts), integrating the
    Kaplan-Meier survival step function in a single C++ scan.
  - [`milestone_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md):
    two-group comparison of Kaplan-Meier survival at a milestone
    timepoint, with Wald, log-log, and MOVER inference methods.
  - [`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md):
    max-combo test over a set of Fleming-Harrington weighted log-rank
    statistics, with the joint p-value obtained from the implied
    multivariate normal distribution.
  - [`rmw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md):
    robust modestly-weighted log-rank test of Magirr and Ohrn, the
    maximum of the standard log-rank and a modestly-weighted log-rank
    statistic, with the joint p-value obtained from the implied
    bivariate normal distribution.
  - [`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md):
    average hazard with survival weight of Uno and Horiguchi, reporting
    the ratio (RAH) and difference (DAH) contrasts.
- [`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
  gains weighted log-rank tests (Fleming-Harrington, modestly-weighted,
  Gehan-Breslow, Tarone-Ware) and stratified and stratified-weighted
  variants, all sharing the single-scan C++ backend.
- New simulation layer:
  - [`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
    extended with optional subgroups defined by a prevalence
    specification, and rewritten so the entire generation pipeline runs
    in a single C++ kernel that materializes the output data frame once.
  - [`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md):
    interim or sequential analysis of simulated data at one or more
    looks, defined by target event counts or calendar times, computed by
    a fused C++ kernel that reuses the analysis cores of the standalone
    functions. Supports subgroup analyses.
  - [`simsummary_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md):
    operating-characteristic summary (rejection and futility rates,
    stopping-look distribution, expected timing) from
    [`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
    output and supplied group-sequential boundaries.
- Each estimation and testing function has a corresponding
  [`print()`](https://rdrr.io/r/base/print.html) method, and the print
  methods share a unified display format.

## FastSurvival 0.1.0

CRAN release: 2026-05-27

- Initial release.
- Core computations implemented in `C++` via `Rcpp` for use inside large
  simulation loops.
- [`survfit_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md):
  single-time-point Kaplan-Meier estimator with Greenwood standard error
  and plain / log / log-log confidence intervals. The C++ backend
  locates the evaluation cutoff via binary search and accumulates the
  Kaplan-Meier product and Greenwood variance sum in a single scan over
  event positions. Returns an object of class `"survfit_fast"` with a
  [`print()`](https://rdrr.io/r/base/print.html) method.
- [`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md):
  log-rank test returning a one-sided Z-score or a two-sided chi-square
  statistic. The C++ backend uses a two-pointer merge scan over pooled
  sorted vectors, eliminating the rank construction,
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html), and reverse
  cumulative sum operations of the standard implementation. Returns an
  object of class `"survdiff_fast"` with a
  [`print()`](https://rdrr.io/r/base/print.html) method.
- [`coxph_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md):
  closed-form hazard ratio estimator via the Pike-Halley Estimator
  method with Wald confidence interval. The C++ backend performs group
  splitting, at-risk counting, and per-distinct-event-time accumulation
  in a single pass. Returns an object of class `"coxph_fast"` with a
  [`print()`](https://rdrr.io/r/base/print.html) method.
- [`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md):
  clinical trial data simulator supporting one- and two-group designs,
  piecewise uniform accrual, and simple and piecewise exponential
  survival and dropout times. C++ backends handle piecewise sampling and
  two-group interleaving, and random number generation uses `dqrng`.
