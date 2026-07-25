# FastSurvival: Fast Survival Analysis and Simulation for Clinical Trials

FastSurvival provides fast alternatives to the standard survival
analysis functions in the survival package, together with a simulation
layer for designing and evaluating time-to-event trials. Every function
is designed for repeated evaluation inside large simulation loops, such
as adaptive sample-size re-estimation, probability-of-success
calculations, and regional consistency evaluation in multi-regional
clinical trials. Core computations are implemented in C++ via Rcpp.

## Details

The estimation and testing functions return S3-class objects that are
internally named numeric vectors, so a result can be used directly in
arithmetic, subsetting, and aggregation after stripping the class with
[`unclass`](https://rdrr.io/r/base/class.html). Each class has a
[`print()`](https://rdrr.io/r/base/print.html) method that formats the
result similarly to the corresponding survival package output.

Estimation and testing:

- [`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md):

  Kaplan-Meier survival probability, standard error, and confidence
  interval at a single specified time point. The C++ backend locates the
  evaluation cutoff by binary search and accumulates the Kaplan-Meier
  product and the Greenwood sum in a single scan over event positions
  only.

- [`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md):

  Log-rank test for two-group survival data, including
  Fleming-Harrington, modestly-weighted, Gehan-Breslow, and Tarone-Ware
  weights, and stratified and stratified-weighted variants. The C++
  backend uses a two-pointer merge scan over pooled sorted vectors.

- [`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md):

  Closed-form hazard ratio estimator via the Pike-Halley Estimator
  method, with Wald confidence interval. The estimator anchors at the
  Pike estimate and applies a single analytic Halley correction to the
  Cox partial likelihood score.

- [`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md):

  Restricted mean survival time for a single group or a two-group
  comparison, with difference and ratio contrasts.

- [`wmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wmst_fast.md):

  Window mean survival time over an interval, generalizing
  [`rmst_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmst_fast.md)
  to a nonzero lower window limit.

- [`milestone_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/milestone_fast.md):

  Two-group comparison of Kaplan-Meier survival at a milestone
  timepoint, with Wald, log-log, and MOVER inference methods.

- [`medsurv_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/medsurv_fast.md):

  Median survival time for a single group or a two-group difference,
  with a kernel-hazard variance method and a local-constant-hazard
  method.

- [`maxcombo_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/maxcombo_fast.md):

  Max-combo test over a set of Fleming-Harrington weighted log-rank
  statistics, with the joint p-value obtained from the implied
  multivariate normal distribution.

- [`rmw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/rmw_fast.md):

  Robust modestly-weighted log-rank test, the maximum of the standard
  log-rank and a modestly-weighted log-rank statistic.

- [`wkm_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/wkm_fast.md):

  Weighted Kaplan-Meier (Pepe-Fleming) test, the weighted integrated
  difference between two Kaplan-Meier curves.

- [`ahsw_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahsw_fast.md):

  Average hazard with survival weight, reporting the ratio and
  difference contrasts.

- [`ahr_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/ahr_fast.md):

  Kalbfleisch-Prentice average hazard ratio over a restricted interval,
  with a test on the group-share scale and an equivalent test and
  confidence interval on the log scale.

Simulation and sequential analysis:

- [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md):

  Individual patient data simulator for one-group, two-group, and K-arm
  time-to-event trials. Supports piecewise uniform accrual, piecewise
  exponential survival and dropout times, subgroups, and correlated
  two-endpoint illness-death generation. Random number generation uses
  dqrng.

- [`analysis_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/analysis_fast.md):

  Interim and final analyses of simulated trials at a set of information
  times or calendar cutoffs, returning one row of test statistics per
  replicate and look.

- [`pairwise_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/pairwise_fast.md):

  Pairwise comparisons of each experimental arm against a shared control
  in a multi-arm trial.

- [`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md):

  Aggregation of replicate-level results into operating characteristics,
  given efficacy and futility boundaries computed elsewhere.

Visualization:

- [`gen_scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/gen_scenario_fast.md):

  Design-stage assembly of survival, dropout, and accrual scenarios,
  with [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods.

- [`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md):

  Analysis-stage Kaplan-Meier curves with a risk table, with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods.

## References

Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
estimator. Manuscript under review.

Collett, D. (2014). *Modelling Survival Data in Medical Research* (3rd
ed.). Chapman and Hall/CRC.

## See also

Useful links:

- <https://github.com/gosukehommaEX/FastSurvival>

- <https://gosukehommaEX.github.io/FastSurvival/>

- Report bugs at <https://github.com/gosukehommaEX/FastSurvival/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
