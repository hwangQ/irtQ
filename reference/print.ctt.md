# Print a Combined CTT Analysis

Prints a condensed, headline-style summary of an object of class `"ctt"`
returned by [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md),
following the same brief/full print-summary split used by
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
(compare `print.est_irt()` vs.
[`summary.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md)/`print.summary.est_irt()`).

## Usage

``` r
# S3 method for class 'ctt'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"ctt"`, as returned by
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md).

- digits:

  Number of decimal places used when rounding numeric values for
  display. Default is `3`. Note that the item- and test-level statistics
  bundled into `x` are already rounded to 3 decimal places before
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) returns
  them, so a `digits` value greater than 3 here cannot recover precision
  that was already discarded upstream.

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
