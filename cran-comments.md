## Update

This is an update from version 0.2.0 to 1.0.0. The release fixes several
bugs found in a review of the whole package (among them a wrong one-sided
p-value tail for the log-log milestone test in `analysis_fast()`, incorrect
dropout and pipeline counts in the event-driven mode of `pairwise_fast()`,
and a sample-size split in `simdata_fast()` that could lose a subject),
adds input validation to the analysis functions, and adds a stratified
version of the closed-form hazard ratio estimator in `coxph_fast()`. It also
resolves the two issues reported on GitHub (#1 and #2). See NEWS.md for the
full list of changes.

## Notes for the reviewer

The checks below reported no NOTE. If the incoming check reports possibly
misspelled words in the DESCRIPTION, "Kalbfleisch" and "Pepe", both are author
surnames, used to name the Kalbfleisch-Prentice average hazard ratio and the
Pepe-Fleming weighted Kaplan-Meier test. The spelling is correct.

As in the previous release, no example uses \dontrun{}. Examples that exceed
the 5-second limit are wrapped in \donttest{}, and examples that use Suggests
packages other than the recommended package 'survival' are guarded with
requireNamespace(). The package was checked with --run-donttest.

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.0
* win-builder: R-release (R 4.6.1)
* GitHub Actions (R-CMD-check workflow):
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - windows-latest (R release)
  - macos-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies.
