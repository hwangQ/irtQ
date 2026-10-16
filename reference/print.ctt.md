# Print a Combined CTT Analysis

Prints a condensed, headline-style summary of an object of class `"ctt"`
returned by [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md),
following the same brief/full print-summary split used by
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
(compare
[`print.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md)
vs.
[`summary.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md)/[`print.summary.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md)).

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

  Number of decimal places used for display. Default is `3`. The
  statistics in `x` are already rounded to three decimal places by
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md), so a value
  above 3 adds no precision.

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

## Examples

``` r
# simulate the responses of 300 examinees to 15 dichotomous 3PLM items
set.seed(1)
x <- shape_df(
  par.drm = list(a = rep(1.5, 15), b = seq(-1.5, 1.5, length.out = 15),
                 g = rep(0.2, 15)),
  cats = 2, model = "3PLM"
)
dat <- simdat(x = x, theta = rnorm(300), D = 1)

# run the CTT analysis
out <- ctt(data = dat)

# print the condensed report
print(out)
#> 
#> Call:
#> ctt(data = dat)
#> 
#> Classical Test Theory (CTT) Analysis
#>  Number of items: 15
#>  Number of examinees: 300
#>  Cronbach's alpha: 0.675
#>  SEM: 1.648
#>  Mean difficulty: 0.605
#>  Mean discrimination (raw / corrected): 0.424 / 0.284
#>  Flagged items: 0 of 15
#> 
#> Use summary() for the full item-level and frequency-distribution report.
```
