## Update

This is an update from version 1.1.0 to 1.2.0. It corrects errors and
documentation found while preparing an article on the package, and adds no
new features:

* `print.simsummary_fast()` failed for a selection of the columns of a
  `simsummary_fast()` result;
* the weighted Kaplan-Meier test returned NaN or an infinite value, instead of
  NA, when an observed time was infinite (for example a cure fraction without
  dropout in simulated data);
* the documentation of `maxcombo_fast()` and `analysis_fast()` described the
  precision of the max-combo p-value and the `mc.alpha` shortcut incorrectly.
  The correlation matrix of the default weights is singular, so the
  quasi-Monte-Carlo integration does not reach the requested tolerance. This
  is now stated, and new tests compare the p-values with an independent
  numerical integral.

See NEWS.md for the full list of changes.

## Notes for the reviewer

The checks below reported no NOTE. If the incoming check reports possibly
misspelled words in the DESCRIPTION, "Kalbfleisch" and "Pepe" are author
surnames, used to name the
Kalbfleisch-Prentice average hazard ratio and the Pepe-Fleming weighted
Kaplan-Meier test. The spelling is correct.

As in the previous releases, no example uses \dontrun{}. Examples that exceed
the 5-second limit are wrapped in \donttest{}, and examples that use Suggests
packages other than the recommended package 'survival' are guarded with
requireNamespace(). The package was checked with --run-donttest.

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.0
* win-builder: R-release (R 4.6.1)
* win-builder: R-devel (2026-10-05 r90641)
* GitHub Actions (R-CMD-check workflow):
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - windows-latest (R release)
  - macos-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies.
