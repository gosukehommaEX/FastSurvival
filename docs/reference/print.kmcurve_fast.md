# Print Method for Kaplan-Meier Curves

Prints a compact per-group summary of a `kmcurve_fast` object: the group
label and role, the sample size, the number of events, and the median
survival read from the Kaplan-Meier curve. In the two-group case the
constant Cox proportional-hazards estimate of the hazard ratio is
reported below the table.

## Usage

``` r
# S3 method for class 'kmcurve_fast'
print(x, ...)
```

## Arguments

- x:

  A `kmcurve_fast` object from
  [`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md).

- ...:

  Further arguments, currently ignored.

## Value

The input `x`, invisibly.

## See also

[`kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/kmcurve_fast.md),
[`plot.kmcurve_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/plot.kmcurve_fast.md)
