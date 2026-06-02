# FastSurvival: Fast Survival Analysis Functions for Simulation Studies

FastSurvival provides fast alternatives to the standard survival
analysis functions in the survival package. Every function is designed
for repeated evaluation inside large simulation loops, such as adaptive
sample-size re-estimation, probability-of-success calculations, and
regional consistency evaluation in multi-regional clinical trials. Core
computations are implemented in C++ via Rcpp for maximum performance
inside simulation loops.

## Details

The package exports four functions. The three analysis functions return
S3-class objects (`"survfit_fast"`, `"survdiff_fast"`, `"coxph_fast"`)
with [`print()`](https://rdrr.io/r/base/print.html) methods that format
the results similarly to the corresponding survival package output. Each
object is internally a numeric vector, so it can be used directly in
arithmetic, subsetting, and aggregation after stripping the class with
[`unclass`](https://rdrr.io/r/base/class.html).

- [`survfit_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survfit_fast.md):

  Kaplan-Meier survival probability, standard error, and confidence
  interval at a single specified time point. The C++ backend locates the
  evaluation cutoff via binary search and accumulates the Kaplan-Meier
  product and Greenwood sum in a single scan over event positions only,
  without constructing intermediate vectors, making it approximately 50
  times faster than
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) plus
  [`summary()`](https://rdrr.io/r/base/summary.html) for repeated
  single-time-point evaluations.

- [`survdiff_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/survdiff_fast.md):

  Log-rank test statistic (one-sided Z-score or two-sided chi-square)
  for two-group survival data. The C++ backend uses a two-pointer merge
  scan over pooled sorted vectors, eliminating the
  [`rank()`](https://rdrr.io/r/base/rank.html),
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html), and
  `rev(cumsum(rev(...)))` overhead of the standard implementation,
  making it approximately 40 times faster than
  [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html).

- [`coxph_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/coxph_fast.md):

  Closed-form hazard ratio estimator via the Pike-Halley Estimator
  method, with Wald confidence interval. The Pike-Halley Estimator
  anchors at the Pike estimate and applies a single analytic Halley
  correction to the Cox partial likelihood score, achieving residual
  error of order O_p(n^{-3/2}) relative to the Cox maximum likelihood
  estimate. The C++ backend performs group splitting, at-risk counting,
  and per-distinct-event-time accumulation in a single pass, making it
  approximately 30 times faster than
  [`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

- [`simdata_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simdata_fast.md):

  Clinical trial data simulator for one- and two-group time-to-event
  trials. Supports piecewise uniform accrual, simple and piecewise
  exponential survival and dropout times. C++ backends handle piecewise
  sampling and two-group interleaving, and random number generation uses
  dqrng, making it approximately 4 times faster than an equivalent
  pure-R implementation.

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
