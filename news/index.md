# Changelog

## FastSurvival 0.2.0

- New estimation and testing functions:
  - [`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md):
    restricted mean survival time for a single group or a two-group
    comparison (difference and ratio contrasts), integrating the
    Kaplan-Meier survival step function in a single C++ scan.
  - [`wmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md):
    window mean survival time over an interval, generalizing
    [`rmst_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
    (which is the special case with a lower window limit of zero), for a
    single group or a two-group difference, computed in the same single
    C++ scan with a Greenwood-type variance in which each event time
    contributes its squared remaining window area.
  - [`milestone_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md):
    two-group comparison of Kaplan-Meier survival at a milestone
    timepoint, with Wald, log-log, and MOVER inference methods.
  - [`medsurv_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md):
    median survival time for a single group or a two-group difference,
    with a native kernel-hazard variance method and an `nph`-compatible
    local-constant-hazard method that reproduces the median comparison
    of the `nph` package to numerical precision; the point estimate is
    the same under both methods.
  - [`maxcombo_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md):
    max-combo test over a set of Fleming-Harrington weighted log-rank
    statistics, with the joint p-value obtained from the implied
    multivariate normal distribution.
  - [`rmw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md):
    robust modestly-weighted log-rank test of Magirr and Öhrn, the
    maximum of the standard log-rank and a modestly-weighted log-rank
    statistic, with the joint p-value obtained from the implied
    bivariate normal distribution.
  - [`wkm_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md):
    weighted Kaplan-Meier (Pepe-Fleming) test, the weighted integrated
    difference between two Kaplan-Meier curves, with Pepe-Fleming,
    square-root, and constant weights, reproducing the weighted
    Kaplan-Meier statistic of the `nphsim` package.
  - [`ahsw_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md):
    average hazard with survival weight of Uno and Horiguchi, reporting
    the ratio (RAH) and difference (DAH) contrasts.
  - [`ahr_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md):
    Kalbfleisch-Prentice average hazard ratio between two groups over a
    restricted interval, the estimator used by Dormuth et al. (2024) for
    sample-size calculation under non-proportional hazards, with a test
    on the group-share scale and an equivalent test and confidence
    interval on the log scale.
- [`survdiff_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md)
  gains weighted log-rank tests (Fleming-Harrington, modestly-weighted,
  Gehan-Breslow, Tarone-Ware) and stratified and stratified-weighted
  variants, all sharing the single-scan C++ backend.
- New simulation layer:
  - [`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md)
    extended with optional subgroups defined by a prevalence
    specification and a flexible accrual specification: `a.rate` gives
    absolute accrual rates (with the end of an open final interval
    solved from the total when a trailing rate is supplied) and `a.prop`
    gives accrual proportions, with deterministic per-interval accrual
    counts. The entire generation pipeline runs in a single C++ kernel
    that materializes the output data frame once.
  - [`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md):
    interim or sequential analysis of simulated data at one or more
    looks, defined by target event counts or calendar times, computed by
    a fused C++ kernel that reuses the analysis cores of the standalone
    functions. Supports subgroup analyses.
  - [`simsummary_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md):
    operating-characteristic summary (rejection and futility rates,
    stopping-look distribution, expected timing) from
    [`analysis_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md)
    output and supplied group-sequential boundaries, with a
    [`print()`](https://rdrr.io/r/base/print.html) method that lays the
    results out as a group-sequential design report.
- New visualization layer:
  - [`gen_scenario_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/gen_scenario_fast.md):
    assembles one or more two-group scenarios into a `scenario_fast`
    object for design-stage exploration, with a
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
    draws the analytic survival curves and the piecewise hazard ratio of
    each scenario and a [`print()`](https://rdrr.io/r/base/print.html)
    method that summarizes the medians, the start and end hazard ratios,
    and whether the curves cross.
  - [`kmcurve_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md):
    builds the Kaplan-Meier curves of a single trial realization (for
    example one replicate of
    [`simdata_fast()`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md))
    into a `kmcurve_fast` object, with a
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
    draws the curves with optional restricted-mean shading and a
    smoothed time-varying hazard-ratio panel, and a
    [`print()`](https://rdrr.io/r/base/print.html) method that
    summarizes the events and medians.
- Each estimation and testing function has a corresponding
  [`print()`](https://rdrr.io/r/base/print.html) method, and the print
  methods share a unified display format.
- New vignettes accompany the analysis, simulation, and visualization
  layers: validation against established packages, a speed comparison, a
  group-sequential design reproduction, a log-rank versus RMST
  comparison under nonproportional hazards, the Freidlin-Korn
  strong-null investigation, a correlated PFS and OS group-sequential
  design under the Fleischer model, and a multiregional
  regional-consistency evaluation.

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
