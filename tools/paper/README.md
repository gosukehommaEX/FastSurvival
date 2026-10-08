# Scripts for the software article

These scripts reproduce the comparisons and benchmarks planned for the article
on FastSurvival. They are not part of the package build (`tools/` is listed in
`.Rbuildignore`). Run them from the package root, after installing the package
with `devtools::install()`. Each script writes its results and the session
information to `tools/paper/output/`.

| Script | Content | Packages used |
|--------|---------|---------------|
| `bench_gsd.R` | Two-arm group-sequential design with two event-driven looks: power, analysis times, and elapsed time per simulated trial | gsDesign, simtrial, TrialSimulator |
| `bench_crossover.R` | Crossover after a positive PFS analysis in an illness-death model: PFS and OS power, analysis times, and elapsed time | TrialSimulator |
| `bench_scaling.R` | Elapsed time and data size against the number of simulated trials, and identical results of sequential and parallel batched runs with dqrng streams; max-combo timing with and without `mc.alpha` | gsDesign, parallel |

The operating characteristics of the packages should agree within Monte Carlo
error (the standard error of a power estimate `p` from `nsim` trials is
`sqrt(p * (1 - p) / nsim)`). TrialSimulator and simtrial run fewer simulated
trials because they are slower; the elapsed time is compared per simulated
trial.
