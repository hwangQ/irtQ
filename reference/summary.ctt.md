# Summarize a Combined CTT Analysis

Prepares the full, detailed report for an object of class `"ctt"`
returned by [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md).
Following the convention used by
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)/[`summary.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md),
this method returns an object of class `"summary.ctt"` whose own
[`print.summary.ctt()`](https://hwangQ.github.io/irtQ/reference/print.summary.ctt.md)
method displays the complete item-level table, test-level reliability
summary, and total-score frequency distribution.

## Usage

``` r
# S3 method for class 'ctt'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"ctt"`, as returned by
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md).

- ...:

  Additional arguments passed to or from other methods (currently not
  used).

## Value

An object of class `"summary.ctt"`: a list with the same `item`, `crit`,
`alpha`, `freq`, and `call` elements as `object` (see
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md)'s **Value**),
to be displayed by
[`print.summary.ctt()`](https://hwangQ.github.io/irtQ/reference/print.summary.ctt.md).

## See also

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md),
[`print.ctt()`](https://hwangQ.github.io/irtQ/reference/print.ctt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
