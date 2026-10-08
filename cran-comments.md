## Update

This is an update from version 1.0.0 to 1.1.0. It adds features for the
simulation of clinical trials:

* `cutoff_fast()` (new) computes the calendar time of each analysis in every
  simulated trial from combined event, calendar-time, and enrollment rules;
* `switch_fast()` (new) applies treatment switching to simulated data;
* `analysis_fast()` and `pairwise_fast()` gain a `cutoff.looks` argument,
  `analysis_fast()` gains an `mc.alpha` argument that avoids most of the
  multivariate normal integrals of the max-combo p-values when only the
  decisions at given levels are needed, and `simdata_fast()` gains a `stream`
  argument for reproducible simulation in batches.

It also fixes an error of `analysis_fast()` for the two-sided max-combo test
with two or three weights. A new vignette uses the 'rpsftm' package, which is
added to Suggests and used only when it is installed; the 'rpact' package is no
longer used and is removed from Suggests. See NEWS.md for the full
list of changes.

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

* Local: Windows 11 x64 (build 26200), R 4.6.0 [TO BE CONFIRMED]
* win-builder: R-release (R 4.6.1) [TO BE CONFIRMED]
* win-builder: R-devel [TO BE CONFIRMED]
* GitHub Actions (R-CMD-check workflow) [TO BE CONFIRMED]:
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - windows-latest (R release)
  - macos-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies.
