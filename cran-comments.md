## Update

This is an update from version 1.1.0 to 1.2.0. It is submitted soon after
version 1.1.0 (published on 2026-10-08) because it corrects results that the
published version returns without a warning, found in an independent review
of the package while preparing an article on it. It adds no new features. The
main corrections are:

* the median survival time of `medsurv_fast()` and `analysis_fast()` was
  infinite, with a p-value of 0, when a Kaplan-Meier curve stayed at 0.5 up
  to an infinite observed time;
* `simdata_fast()` truncated a group size such as `90 * 0.7` to the integer
  below and, with fixed subgroup sizes, wrote past the end of a buffer in the
  C++ code;
* the modestly-weighted log-rank test reduced to the ordinary log-rank test
  when the pooled Kaplan-Meier estimate reached 0 before `t_star`;
* `print.simsummary_fast()` printed boundaries next to the wrong looks for
  some selections of rows, and failed for others;
* the weighted Kaplan-Meier test returned NaN or an infinite value, instead of
  NA, when an observed time was infinite;
* several functions accepted a factor event indicator or negative times and
  returned wrong results; these inputs are now errors.

The documentation was also corrected, including the precision of the
max-combo p-value. See NEWS.md for the full list of changes.

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
