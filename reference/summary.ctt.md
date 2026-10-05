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

# create the full report object and print it
summary(out)
#> 
#> Call:
#> ctt(data = dat)
#> 
#> Item-Level Statistics
#>   item  cats  difficulty  discrimination_raw  discrimination_corrected
#>     V1     2       0.887               0.309                     0.205
#>     V2     2       0.837               0.369                     0.250
#>     V3     2       0.813               0.479                     0.365
#>     V4     2       0.763               0.489                     0.364
#>     V5     2       0.790               0.329                     0.195
#>     V6     2       0.697               0.403                     0.257
#>     V7     2       0.637               0.514                     0.375
#>     V8     2       0.630               0.498                     0.356
#>     V9     2       0.560               0.376                     0.215
#>    V10     2       0.500               0.497                     0.349
#>    V11     2       0.467               0.470                     0.319
#>    V12     2       0.457               0.459                     0.307
#>    V13     2       0.407               0.382                     0.223
#>    V14     2       0.313               0.395                     0.247
#>    V15     2       0.323               0.386                     0.236
#>   alpha_removed  flag
#>           0.668      
#>           0.664      
#>           0.651      
#>           0.649      
#>           0.670      
#>           0.663      
#>           0.646      
#>           0.649      
#>           0.669      
#>           0.650      
#>           0.654      
#>           0.656      
#>           0.668      
#>           0.664      
#>           0.666      
#> 
#> Flagging thresholds: difficulty in [0.1, 0.95], discrimination >= 0.2
#> 0 of 15 item(s) flagged.
#> 
#> Test-Level Reliability Summary
#>   n_examinee  n_item  alpha  alpha_std    sem  mean_difficulty
#>          300      15  0.675      0.675  1.648            0.605
#>   mean_discrimination_raw  mean_discrimination_corrected
#>                     0.424                          0.284
#> 
#> Total-Score Frequency Distribution
#>   score  freq    pct  cum_pct
#>       3     6   2.00     2.00
#>       4    16   5.33     7.33
#>       5    16   5.33    12.67
#>       6    24   8.00    20.67
#>       7    32  10.67    31.33
#>       8    27   9.00    40.33
#>       9    37  12.33    52.67
#>      10    42  14.00    66.67
#>      11    38  12.67    79.33
#>      12    23   7.67    87.00
#>      13    22   7.33    94.33
#>      14    10   3.33    97.67
#>      15     7   2.33   100.00
#> 
```
