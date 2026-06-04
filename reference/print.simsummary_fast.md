# Print Method for Sequential Analysis Summaries

Formats an object returned by
[`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md)
for display. The boundary settings are shown first, then a
per-population table whose look rows carry the stage-wise stopping
probabilities and analysis timing and whose `overall` row carries the
totals.

## Usage

``` r
# S3 method for class 'simsummary_fast'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `"simsummary_fast"` from
  [`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md).

- digits:

  A single positive integer, the number of significant digits used when
  rounding the printed values. Defaults to 4.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Details

The probability columns are printed as proportions rounded to `digits`
significant digits, and the timing columns (mean enrolled, mean events,
mean dropouts, mean pipeline, mean calendar time) are rounded for
readability. When the object holds more than one population, each
population is printed as its own block with a heading. The underlying
object is an ordinary data frame, so the unrounded values remain
available by subsetting it directly.

## See also

[`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md).
