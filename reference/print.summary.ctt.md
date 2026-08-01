# Print the Full Report for a Combined CTT Analysis

Displays the complete classical test theory (CTT) report for an object
of class `"summary.ctt"` produced by
[`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md):
the function call, the full item-level statistics table (including
flags, if present), the test-level reliability summary, and the
total-score frequency distribution, following the sectioned reporting
style used by `print.summary.est_irt()` in irtQ.

## Usage

``` r
# S3 method for class 'summary.ctt'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"summary.ctt"`, as returned by
  [`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md).

- digits:

  Number of decimal places used when rounding numeric values for
  display. Default is `3`. Note that the statistics bundled into `x` are
  already rounded to 3 decimal places upstream, so a `digits` value
  greater than 3 here cannot recover precision that was already
  discarded; `digits` is only useful for displaying the report at 3 or
  fewer decimal places.

- ...:

  Additional arguments passed to or from other methods (currently not
  used).

## Value

`x`, invisibly.

## See also

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md),
[`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
