# FastSurvival

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/FastSurvival)](http://cran.r-project.org/package=FastSurvival)
<!-- badges: end -->

The goal of FastSurvival is to provide efficient implementations of (1) log-rank test, (2) hazard ratio estimation, and (3) simulation dataset generation for survival analysis.

## Features

- **Log-rank Test**: Efficient computation of log-rank test statistics
- **Multiple Hazard Ratio Estimation Methods**: Support for Person-Year, Pike, Peto, Log-rank based, and Cox regression methods
- **High-Performance Simulation**: Generate simulation datasets with piecewise exponential and uniform distributions
- **Clinical Trial Support**: Specialized functions for clinical trial data simulation and analysis
- **Optimized R Implementation**: Core functions optimized for performance using efficient R programming techniques

## Installation

You can install the development version of FastSurvival from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/FastSurvival")
```

## Usage

### Log-rank Test

Compare the log-rank test results between FastSurvival and the standard `survdiff()` function:

``` r
library(survival)
library(FastSurvival)
library(microbenchmark)

# Compare log-rank test results
survdiff_result <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)
lrtest_result <- lrtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 2)

print(list(
  survdiff_chisq = survdiff_result[['chisq']], 
  lrtest_chisq = lrtest_result
))

# Performance comparison
microbenchmark(
  survdiff_result = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
  lrtest_result = lrtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 2),
  times = 100
)
```

### Hazard Ratio Estimation

Compare hazard ratio estimates using different methods:

``` r
library(survival)
library(FastSurvival)

# Cox proportional hazards model (reference)
coxph_result <- data.frame(
  method = 'CoxPH', 
  HR = exp(coef(coxph(Surv(futime, fustat) ~ rx, data = ovarian)))
)

# FastSurvival methods
methods <- c('PY', 'Pike', 'Peto', 'LR', 'Cox')
esthr_results <- lapply(methods, function(method) {
  esthr(ovarian$futime, ovarian$fustat, ovarian$rx, 2, method)
})

# Combine results
all_results <- do.call(rbind, c(list(coxph_result), esthr_results))
print(all_results)

# Performance comparison
microbenchmark(
  CoxPH = coxph(Surv(futime, fustat) ~ rx, data = ovarian),
  esthr_PY = esthr(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'PY'),
  esthr_Pike = esthr(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'Pike'),
  esthr_Peto = esthr(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'Peto'),
  esthr_LR = esthr(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 'LR'),
  times = 100
)
```

### Random Number Generation

Generate random numbers from piecewise distributions:

``` r
library(FastSurvival)

# Piecewise exponential distribution
time_points <- c(0, 5, 10, Inf)
hazard_rates <- c(0.1, 0.3, 0.2)
samples_exp <- rpieceexp(1000, time_points, hazard_rates)

# Piecewise uniform distribution  
time_intervals <- c(0, 6, 12, 18, 24)
intensities <- c(10, 20, 30, 40)
samples_unif <- rpieceunif(100, time_intervals, intensity = intensities)

# Plot histograms
par(mfrow = c(1, 2))
hist(samples_exp, main = "Piecewise Exponential", xlab = "Time")
hist(samples_unif, main = "Piecewise Uniform", xlab = "Time")
```

### Clinical Trial Simulation

Simulate clinical trial data with complex scenarios:

``` r
library(FastSurvival)
library(data.table)

# Simulate trial data for two groups
trial_data <- simTrial(
  nsim = 100,
  N = list(group1 = 100, group2 = 100),
  a.time = c(0, 24),
  intensity = 200/24,
  e.time = list(
    group1 = c(0, Inf),
    group2 = c(0, Inf)
  ),
  e.hazard = list(
    group1 = log(2) / 12,
    group2 = log(2) / 15
  ),
  d.time = list(
    group1 = c(0, Inf),
    group2 = c(0, Inf)
  ),
  d.hazard = list(
    group1 = -log(1 - 0.1) / 12,
    group2 = -log(1 - 0.1) / 12
  )
)

# Create analysis datasets
analysis_data <- analysisData(trial_data, E = c(50, 100, 150))

head(analysis_data)

# Trial without dropout
trial_data_no_dropout <- simTrial(
  nsim = 100,
  N = list(control = 100, treatment = 100),
  a.time = c(0, 24),
  intensity = 200/24,
  e.time = list(
    control = c(0, Inf),
    treatment = c(0, Inf)
  ),
  e.hazard = list(
    control = log(2) / 12,
    treatment = log(2) / 15
  ),
  d.time = NULL,  # No dropout
  d.hazard = NULL
)
```

## Performance Benefits

FastSurvival provides significant performance improvements over standard R survival analysis functions:

- **Log-rank test**: ~10-15x faster than `survdiff()`
- **Hazard ratio estimation**: ~5-10x faster than `coxph()` for simple cases
- **Simulation**: Optimized R implementation for large-scale simulations

## Available Functions

### Core Functions
- `lrtest()`: Log-rank test calculation
- `esthr()`: Hazard ratio estimation with multiple methods

### Distribution Functions  
- `rpieceexp()`: Random number generation from piecewise exponential distribution
- `rpieceunif()`: Random number generation from piecewise uniform distribution

### Simulation Functions
- `simData()`: Basic survival data simulation
- `simTrial()`: Clinical trial data simulation
- `analysisData()`: Create analysis datasets from trial simulations

### Analysis Functions
- `overallSummary()`: Overall population analysis summary
- `subgroupSummary()`: Subgroup analysis summary
- `overallResults()`: Adaptive trial stopping probabilities

## References

The methodology is based on established survival analysis techniques and optimized implementations for computational efficiency in clinical trial simulations and survival data analysis.
