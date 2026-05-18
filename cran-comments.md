## Resubmission

This is a resubmission. In response to the previous review by
Benjamin Altmann, the following changes have been made:

* DESCRIPTION: wrapped the software name `C++` in single quotes in the
  Description field.
* `coxph_fast.Rd`: replaced `\dontrun{}` with `\donttest{}` in the example
  that uses `microbenchmark` with `times = 1000`, which cannot be executed
  in under 5 seconds.
* `survdiff_fast.Rd`: applied the same `\dontrun{}` to `\donttest{}`
  replacement for consistency.

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.0
* GitHub Actions (R-CMD-check workflow):
  - ubuntu-latest (R release)
  - ubuntu-latest (R devel)
  - windows-latest (R release)
  - macos-latest (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies (initial submission).
