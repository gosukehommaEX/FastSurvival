# FastSurvival 0.1.0

* Initial release.
* `survfit_fast()`: single-time-point Kaplan-Meier estimator with Greenwood
  standard error and plain / log / log-log confidence intervals.
* `survdiff_fast()`: log-rank test returning a one-sided Z-score or a
  two-sided chi-square statistic.
* `coxph_fast()`: closed-form hazard ratio estimator via the Pike-Halley
  method (PiHE) with Wald confidence interval.
* `simdata_fast()`: clinical trial data simulator supporting one- and
  two-group designs, piecewise uniform accrual, simple and piecewise
  exponential survival and dropout times, and fast random number generation
  via `dqrng`.
