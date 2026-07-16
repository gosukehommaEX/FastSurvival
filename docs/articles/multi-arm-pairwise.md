# Multi-arm trials: K-arm generation and pairwise comparisons against a shared control

## Overview

A multi-arm trial randomizes subjects to a shared control and two or
more experimental arms, and each experimental arm is compared against
the same control. This is the standard screening design for several
candidate doses or regimens, where the goal is to identify which
experimental arms beat control while controlling the family-wise error
rate across the comparisons.

`simdata_fast` generates such a trial directly. When `n` is a vector of
length greater than two together with a per-arm survival list, each arm
is generated with the validated single-group kernel over a common
accrual window, and the arms are stacked into one data frame with a
`group` column labeled 1 to `length(n)`. The analysis path keeps its
two-group contract, so a multi-arm design is analyzed as a set of
pairwise contrasts. Each experimental arm is compared against the
control by subsetting the generated data to those two arms and calling
`analysis_fast` once per contrast. Because every contrast reuses the
same control arm, the pairwise statistics are positively correlated, and
that correlation is what a Dunnett-type adjustment exploits to be less
conservative than a Bonferroni correction.

This vignette generates a three-arm trial (one control and two
experimental arms), runs the two pairwise log-rank tests, and reports
the operating characteristics of the screening design under a Bonferroni
correction, including the disjunctive power of finding at least one
winning arm and the probability of selecting the truly better
experimental arm. It then confirms empirically that the family-wise
error rate is controlled under the global null.

## Design parameters

``` r

nsim <- 2000
seed <- 20260711

# Three arms: group 1 = control, groups 2 and 3 = experimental.
# Per-arm sample size (balanced 1:1:1 randomization).
n_arm <- c(150, 150, 150)

# Median survival by arm, in months. The two experimental arms have an
# increasing benefit over control.
mst <- c(control = 12, dose_low = 16, dose_high = 20)

# Piecewise-uniform accrual: constant rate over 12 months enrolling all 450.
a_time <- c(0, 12)
a_rate <- sum(n_arm) / 12

# Single calendar analysis time, in months.
analysis_month <- 30

# One-sided family-wise significance level and the number of experimental arms.
alpha_total <- 0.025
n_contrast  <- length(n_arm) - 1L

# Implied exponential hazard ratios of each experimental arm versus control.
data.frame(
  arm = names(mst)[-1],
  median = mst[-1],
  HR_vs_control = round(mst["control"] / mst[-1], 3),
  row.names = NULL
)
#>         arm median HR_vs_control
#> 1  dose_low     16          0.75
#> 2 dose_high     20          0.60
```

The experimental arm with the longer median has the smaller hazard ratio
and is the stronger candidate.

## Data generation

A single `simdata_fast` call generates all three arms. The per-arm
medians are supplied as an `e.median` list with one element per arm, and
the arms share the accrual window and the master seed.

``` r

df <- simdata_fast(
  nsim     = nsim,
  n        = n_arm,
  a.time   = a_time,
  a.rate   = a_rate,
  e.median = list(mst["control"], mst["dose_low"], mst["dose_high"]),
  seed     = seed
)

# One row per subject per simulation; the group column labels the arm.
table(df$group[df$sim == 1])
#> 
#>   1   2   3 
#> 150 150 150
```

## Pairwise analysis against the shared control

Each experimental arm is compared against the control by restricting the
data to that arm and the control and calling `analysis_fast` with
`control = 1`. On a two-arm subset the analysis is the ordinary
two-group log-rank test, one-sided in the direction of treatment
benefit, with the natural sign convention in which a negative Z favors
the experimental arm.

``` r

control_arm <- 1L
exp_arms    <- setdiff(sort(unique(df$group)), control_arm)

# Matrix of one-sided log-rank statistics and p-values, one column per contrast.
Z <- matrix(NA_real_, nrow = nsim, ncol = n_contrast)
P <- matrix(NA_real_, nrow = nsim, ncol = n_contrast)
colnames(Z) <- colnames(P) <- paste0("arm", exp_arms, "_vs_ctrl")

for (k in seq_along(exp_arms)) {
  sub <- df[df$group %in% c(control_arm, exp_arms[k]), ]
  res <- analysis_fast(
    sub,
    control    = control_arm,
    time.looks = analysis_month,
    stat       = "logrank",
    side       = 1
  )
  Z[, k] <- res$logrank.z
  P[, k] <- res$logrank.p
}

head(round(P, 4))
#>      arm2_vs_ctrl arm3_vs_ctrl
#> [1,]       0.0016        0e+00
#> [2,]       0.0089        0e+00
#> [3,]       0.0111        0e+00
#> [4,]       0.0010        2e-04
#> [5,]       0.0000        0e+00
#> [6,]       0.1123        7e-04
```

## The `pairwise_fast()` convenience wrapper

The loop above is packaged by `pairwise_fast`, which runs the
control-versus-arm analysis for every experimental arm and stacks the
results into one long data frame with an `arm` column. With `time.looks`
all contrasts share the fixed calendar cutoff, and
`adjust = "bonferroni"` appends a `p.adj` column that multiplies each
p-value by the number of contrasts.

``` r

pw <- pairwise_fast(
  df,
  control    = control_arm,
  time.looks = analysis_month,
  stat       = "logrank",
  side       = 1,
  adjust     = "bonferroni"
)
head(pw[, c("arm", "sim", "cutoff", "logrank.z", "logrank.p", "p.adj")])
#>   arm sim cutoff logrank.z    logrank.p        p.adj
#> 1   2   1     30 -2.949901 1.589381e-03 3.178762e-03
#> 2   2   2     30 -2.369475 8.906683e-03 1.781337e-02
#> 3   2   3     30 -2.286564 1.111065e-02 2.222130e-02
#> 4   2   4     30 -3.101323 9.632898e-04 1.926580e-03
#> 5   2   5     30 -4.009756 3.039077e-05 6.078154e-05
#> 6   2   6     30 -1.214419 1.122938e-01 2.245877e-01
```

The marginal per-contrast rejection rates from the wrapper reproduce the
manual loop above.

``` r

round(tapply(pw$logrank.p <= alpha_total, pw$arm, mean), 3)
#>     2     3 
#> 0.539 0.941
```

An event-driven screening design instead lets a designated primary
contrast fix the data cutoff. Here arm 3 is the primary comparison, its
target of 200 control-plus-arm-3 events defines the per-simulation
calendar cutoff, and both contrasts are analyzed at that same cutoff.

``` r

pw_ev <- pairwise_fast(
  df,
  control     = control_arm,
  event.looks = 200,
  primary     = 3,
  stat        = "logrank",
  side        = 1,
  adjust      = "bonferroni"
)
head(pw_ev[, c("arm", "sim", "cutoff", "reached", "logrank.z", "p.adj")])
#>   arm sim   cutoff reached logrank.z        p.adj
#> 1   2   1 30.41747    TRUE -2.948070 0.0031976438
#> 2   2   2 28.86407    TRUE -2.304046 0.0212200415
#> 3   2   3 30.85221    TRUE -2.390590 0.0168213177
#> 4   2   4 30.28744    TRUE -3.230967 0.0012337209
#> 5   2   5 30.39027    TRUE -3.860549 0.0001131327
#> 6   2   6 28.89168    TRUE -1.090690 0.2754092176
```

The two contrasts at a given simulation carry the same `cutoff`, which
is the calendar time at which the primary contrast reached its event
target.

## Operating characteristics under the alternative

The screening design rejects the control-versus-arm null for an
experimental arm when its one-sided p-value falls below the
multiplicity-adjusted threshold. A Bonferroni correction spends the
family-wise level equally across the two contrasts, giving a
per-contrast threshold of `alpha_total / n_contrast`. The marginal
(unadjusted) power of each contrast is shown alongside the adjusted
power for reference.

``` r

thr_marginal   <- alpha_total
thr_bonferroni <- alpha_total / n_contrast

reject_marginal   <- P <= thr_marginal
reject_bonferroni <- P <= thr_bonferroni

power_table <- data.frame(
  contrast          = colnames(P),
  power_marginal    = round(colMeans(reject_marginal), 3),
  power_bonferroni  = round(colMeans(reject_bonferroni), 3),
  row.names = NULL
)
power_table
#>       contrast power_marginal power_bonferroni
#> 1 arm2_vs_ctrl          0.539            0.431
#> 2 arm3_vs_ctrl          0.941            0.899
```

The disjunctive power is the probability of declaring at least one
experimental arm superior to control under the Bonferroni rule. The
probability of correctly selecting the stronger experimental arm is the
probability that the arm with the larger true effect is among the
rejected arms.

``` r

any_reject   <- rowSums(reject_bonferroni) > 0
# The stronger experimental arm is the last column (largest median).
stronger_col <- n_contrast
select_stronger <- reject_bonferroni[, stronger_col]

data.frame(
  quantity = c("Disjunctive power (any arm beats control)",
               "Power to declare the stronger arm",
               "Expected number of arms declared"),
  value = c(round(mean(any_reject), 3),
            round(mean(select_stronger), 3),
            round(mean(rowSums(reject_bonferroni)), 3)),
  row.names = NULL
)
#>                                    quantity value
#> 1 Disjunctive power (any arm beats control) 0.911
#> 2         Power to declare the stronger arm 0.899
#> 3          Expected number of arms declared 1.330
```

## Correlation between the contrasts

The two pairwise log-rank statistics share the control arm, so they are
positively correlated. This correlation is the reason a Dunnett
adjustment is less conservative than Bonferroni, since Bonferroni
ignores it and treats the comparisons as if they were independent.

``` r

round(stats::cor(Z, use = "complete.obs"), 3)
#>              arm2_vs_ctrl arm3_vs_ctrl
#> arm2_vs_ctrl        1.000        0.508
#> arm3_vs_ctrl        0.508        1.000
```

For balanced allocation with a shared control the theoretical
correlation of the two standardized statistics is close to one half, and
the empirical value here reproduces that magnitude to simulation error.

## Family-wise error rate under the global null

Regenerating the trial with all three arms sharing the control median
gives the global null, under which no experimental arm truly beats
control. The Bonferroni rule should reject at least one arm no more than
the nominal family-wise level.

``` r

df0 <- simdata_fast(
  nsim     = nsim,
  n        = n_arm,
  a.time   = a_time,
  a.rate   = a_rate,
  e.median = list(mst["control"], mst["control"], mst["control"]),
  seed     = seed + 1
)

P0 <- matrix(NA_real_, nrow = nsim, ncol = n_contrast)
for (k in seq_along(exp_arms)) {
  sub0 <- df0[df0$group %in% c(control_arm, exp_arms[k]), ]
  res0 <- analysis_fast(sub0, control = control_arm,
                        time.looks = analysis_month, stat = "logrank", side = 1)
  P0[, k] <- res0$logrank.p
}

data.frame(
  quantity = c("Per-contrast type I error (nominal)",
               "Family-wise error, Bonferroni"),
  value = c(round(mean(P0 <= alpha_total), 3),
            round(mean(rowSums(P0 <= thr_bonferroni) > 0), 3)),
  row.names = NULL
)
#>                              quantity value
#> 1 Per-contrast type I error (nominal) 0.024
#> 2       Family-wise error, Bonferroni 0.022
```

The per-contrast error is near the nominal one-sided level, and the
Bonferroni family-wise error is at or below the target, confirming that
the screening rule controls the error rate. A Dunnett procedure would
recover some of the conservatism by using the shared-control
correlation, at the cost of the closed-form simplicity of the Bonferroni
threshold.

## Remarks

The multi-arm generation reuses the validated single-group kernel once
per arm, so every survival, dropout, and accrual draw comes from the
same primitives that back the two-group path, and the whole study runs
in a few seconds. The analysis keeps its two-group contract, so
multiplicity and arm selection are handled as a thin post-processing
layer over the pairwise statistics rather than inside the analysis
kernel. This keeps the design flexible. Any of the statistics offered by
`analysis_fast`, not only the log-rank test, can drive the pairwise
comparisons, and any multiple-comparison procedure can be applied to the
resulting p-values.

Multi-arm mode shares the master seed across arms and generates each arm
over a common accrual window, which is adequate for
operating-characteristic screening. It does not reproduce a single
randomization list with a common accrual stream, and it does not combine
with subgroups or the illness-death model, both of which remain
two-group.

### References

Dunnett, C. W. (1955). A multiple comparison procedure for comparing
several treatments with a control. *Journal of the American Statistical
Association*, 50(272), 1096-1121.
