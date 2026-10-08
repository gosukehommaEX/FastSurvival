## Update

This is an update from version 1.0.0 to 1.1.0. It adds features for the
simulation of clinical trials:

* `cutoff_fast()` (new) computes the calendar time of each analysis in every
  simulated trial from combined event, calendar-time, and enrollment rules;
* `switch_fast()` (new) applies treatment switching to simulated data;
* `analysis_fast()` and `pairwise_fast()` gain a `cutoff.looks` argument, and
  `simdata_fast()` gains a `stream` argument for reproducible simulation in
  batches.

It also fixes an error of `analysis_fast()` for the two-sided max-combo test
with two or three weights. A new vignette uses the 'rpsftm' package, which is
added to Suggests and used only when it is installed. See NEWS.md for the full
list of changes.

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

[TO BE FILLED IN FROM THE CHECK RESULTS]

## R CMD check results

[TO BE FILLED IN FROM THE CHECK RESULTS]

## Downstream dependencies

[TO BE CONFIRMED]
