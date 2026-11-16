## Update

This is an update from version 1.2.0 to 2.0.0. The major version changes
because several functions return different results for the same input, and
two arguments were removed. The changes that affect existing code are:

* `coxph_fast()` now returns the maximum partial likelihood estimate of the
  Cox model, which agrees with `survival::coxph()`, instead of the
  closed-form Pike-Halley approximation of earlier versions. Tied event times
  are handled by the Efron method by default (new argument `ties`), and the
  result is NA when the estimate is infinite or the iteration does not
  converge. The same holds for the Cox statistic of `analysis_fast()` and
  `pairwise_fast()`.
* The variance of the weighted Kaplan-Meier test of `wkm_fast()` and of the
  `"wkm"` statistic of `analysis_fast()` is now the pooled estimator of Pepe
  and Fleming (1989), so the standard errors and p-values change.
* The default kernel bandwidth of `medsurv_fast()` and of the `"medsurv"`
  statistic of `analysis_fast()` is now the normal-reference bandwidth of
  the Epanechnikov kernel, so the standard errors of the median change.
* The argument `method` of `medsurv_fast()` and the argument
  `medsurv.method` of `analysis_fast()` were removed; the arguments after
  them move up by one position.
* In the two-group mode, `rmst_fast()` also returns the standard errors of
  the restricted mean survival times of the two groups, so the later
  elements of its result have moved by two positions.

The release also adds support for the native pipe in the formula interface
and the simulation functions, the choice of the endpoint columns in
`analysis_fast()` and `pairwise_fast()`, further treatment-switching rules
in `simdata_fast()`, and a vignette on simulation studies. See NEWS.md for
the full list of changes.

## Notes for the reviewer

If the incoming check reports possibly misspelled words in the DESCRIPTION,
"Kalbfleisch" and "Pepe" are author surnames, used to name the
Kalbfleisch-Prentice average hazard ratio and the Pepe-Fleming weighted
Kaplan-Meier test. The spelling is correct.

As in the previous releases, no example uses \dontrun{}. Examples that exceed
the 5-second limit are wrapped in \donttest{}, and examples that use Suggests
packages other than the recommended package 'survival' are guarded with
requireNamespace(). The package was checked with --run-donttest.

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.1
* win-builder: R-release (R 4.6.2)
* win-builder: R-devel (2026-11-14 r90890)
* GitHub Actions (R-CMD-check workflow):
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - windows-latest (R release)
  - macos-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no reverse dependencies on CRAN (checked on 2026-11-15 with
`tools::package_dependencies()`, `reverse = TRUE`, `which = "all"`).
