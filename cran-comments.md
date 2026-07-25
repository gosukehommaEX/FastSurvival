## Update

This is an update from version 0.1.0 to 0.2.0. It adds several survival
analysis functions (`rmst_fast()`, `wmst_fast()`, `milestone_fast()`,
`medsurv_fast()`, `maxcombo_fast()`, `rmw_fast()`, `wkm_fast()`, `ahsw_fast()`,
`ahr_fast()`), extends `survdiff_fast()` with weighted and stratified log-rank
tests, adds a simulation and sequential-analysis layer (`simdata_fast()`
subgroups and correlated two-endpoint illness-death simulation,
`analysis_fast()`, `pairwise_fast()`, `simsummary_fast()`), and adds a
visualization layer (`gen_scenario_fast()`, `kmcurve_fast()`). See NEWS.md for
the full list of changes.

## Notes for the reviewer

The incoming check reports one NOTE on possibly misspelled words in the
DESCRIPTION, "Kalbfleisch" and "Pepe". Both are author surnames, used to name
the Kalbfleisch-Prentice average hazard ratio and the Pepe-Fleming weighted
Kaplan-Meier test. The spelling is correct.

Following the feedback on the 0.1.0 submission, software names in the title
and description are wrapped in single quotes ('C++', 'Rcpp', 'survival'), and
no example uses \dontrun{}. Examples that exceed the 5-second limit or rely on
Suggests packages are wrapped in \donttest{} and guarded with
requireNamespace(), so they do not fail when those packages are absent. The
package was checked with --run-donttest to confirm this.

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.0
* win-builder: R Under development (unstable) (2026-07-24 r90297 ucrt)
* GitHub Actions (R-CMD-check workflow):
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - windows-latest (R release)
  - macos-latest (R release)

## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE is the DESCRIPTION spelling note described above.

## Downstream dependencies

There are no downstream dependencies.
