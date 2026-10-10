# Print an ahr_fast object

The p-value column is labeled `Pr(>|z|)` for a two-sided test and
`Pr(<z)` for a one-sided test, whose p-value is the lower tail in the
direction of treatment benefit.

## Usage

``` r
# S3 method for class 'ahr_fast'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  an object of class `"ahr_fast"`

- digits:

  number of significant digits to print

- ...:

  further arguments (currently ignored)

## Value

`x`, invisibly
