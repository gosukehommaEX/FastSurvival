# Scripts for the software article

These scripts produce the comparisons and benchmarks reported in the article
on FastSurvival for The R Journal. They are not part of the package build
(`tools/` is listed in `.Rbuildignore`). Install FastSurvival 1.1.0 from CRAN
(`install.packages("FastSurvival", type = "source")`), open
`FastSurvival.Rproj`, and run each script from the package root, for example
`source("tools/paper/data/bench_gsd.R")`. Each script sources
`data/machine_info.R`, which stops if another version of FastSurvival is
installed and writes the computing environment to
`tools/paper/output/machine_info.csv`. Each script writes its results and the
session information to `tools/paper/output/`, and the article reads only these
output files.

| Script | Content | Packages used |
|--------|---------|---------------|
| `data/bench_gsd.R` | Two-arm group-sequential design with two event-driven looks: power, analysis times, and elapsed time per simulated trial | gsDesign, simtrial, TrialSimulator |
| `data/bench_crossover.R` | Crossover after a positive PFS analysis in an illness-death model: PFS and OS power, analysis times, and elapsed time | TrialSimulator |
| `data/bench_functions.R` | Agreement of the analysis functions with reference implementations on CRAN (gbsg data), of `cutoff_fast()` with `simtrial::get_analysis_date()`, and of `analysis_fast()` with `survdiff()` and `rmst2()` on the censored data; median time per call against the references | survival, survRM2, nph, nphRCT, survAH, simtrial, microbenchmark |
| `data/bench_scaling.R` | Elapsed time and data size against the number of simulated trials, and identical results of sequential and parallel batched runs with dqrng streams; max-combo timing with and without `mc.alpha` | gsDesign, parallel |

The operating characteristics of the packages should agree within Monte Carlo
error. The output tables give each estimate with its Monte Carlo standard error
(for a power estimate `p` from `nsim` trials, `sqrt(p * (1 - p) / nsim)`) and
the difference from FastSurvival in units of the combined standard error.
simtrial and TrialSimulator run fewer simulated trials (2,000) than
FastSurvival (10,000) because they are slower; the elapsed time is compared per
simulated trial, with every package running sequentially.
