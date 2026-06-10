# Print Method for Scenario Sets

Prints a compact one-row-per-scenario summary of a `scenario_fast`
object: the scenario label, the null flag, the total sample size, the
control and experimental medians, the hazard ratio at the start and end
of the time window, and whether the two survival curves cross. The
medians and hazard ratios are read from the analytic survival curves, so
they are defined for both exponential and piecewise-exponential
scenarios.

## Usage

``` r
# S3 method for class 'scenario_fast'
print(x, tmax = NULL, ...)
```

## Arguments

- x:

  A `scenario_fast` object from
  [`gen_scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/gen_scenario_fast.md).

- tmax:

  Upper limit of the time window used to read medians, the end-of-
  window hazard ratio, and the crossing flag. Defaults to a value
  derived from the scenarios.

- ...:

  Further arguments, currently ignored.

## Value

The input `x`, invisibly.

## See also

[`gen_scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/gen_scenario_fast.md),
[`plot.scenario_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/plot.scenario_fast.md)
