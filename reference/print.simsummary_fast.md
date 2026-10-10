# Print Method for Sequential Analysis Summaries

Formats an object returned by
[`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md)
in the style of a group-sequential design report. After a header with
the simulation count and the boundary settings, two look-by-look tables
are shown: a stopping-boundary table (information fraction, events,
sample size, the efficacy and futility boundaries, and the cumulative
efficacy crossing probability) and an analysis-timing table (sample
size, events, dropouts, pipeline, analysis time, and the per-look
efficacy and futility crossing probabilities). An overall block reports
the rejection rate and the expected counts and timing at the stopping
look.

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

  A single positive integer, the number of decimal places used for the
  printed probabilities. Defaults to 4.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Details

Column labels follow the convention of group-sequential design software:
`Events (s)`, `Sample (n)`, `Dropouts (d)`, `Pipeline` (the enrolled
count minus events minus dropouts), `Analysis Time` (mean calendar
time), and `Info. Frac.` (the information fraction, computed as the mean
events at a look divided by the mean events at the final look, or from
`look.value` when no event count is available). Probabilities are
printed to `digits` decimal places and counts and times to fewer.
Because the summary is a Monte Carlo estimate under a single
data-generating truth, it does not carry the separate null and
alternative columns or the alpha and beta spending of an analytic design
report. When the summary has several blocks (arms of
[`pairwise_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/pairwise_fast.md)
output or subgroup populations), each block is printed under its own
heading. The object is a data frame with an additional class, so the
unrounded values remain available by subsetting it. A selection of rows
that keeps, in every block, all the looks in their original order (for
example one population or one arm) is printed as a report. Any other
selection of rows, such as some of the looks or rows in a different
order, is printed as an ordinary data frame, because the boundaries and
the information fractions refer to all the looks of the design. A
selection of columns keeps the class but not the boundary settings, and
is printed as an ordinary data frame, as is a selection without any look
row.

## See also

[`simsummary_fast`](https://gosukehommaEX.github.io/FastSurvival/reference/simsummary_fast.md).
