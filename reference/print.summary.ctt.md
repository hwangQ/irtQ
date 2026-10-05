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

# print the full report with values rounded to 2 decimal places
print(summary(out), digits = 2)
#> 
#> Call:
#> ctt(data = dat)
#> 
#> Item-Level Statistics
#>   item  cats  difficulty  discrimination_raw  discrimination_corrected
#>     V1     2        0.89                0.31                      0.20
#>     V2     2        0.84                0.37                      0.25
#>     V3     2        0.81                0.48                      0.36
#>     V4     2        0.76                0.49                      0.36
#>     V5     2        0.79                0.33                      0.20
#>     V6     2        0.70                0.40                      0.26
#>     V7     2        0.64                0.51                      0.38
#>     V8     2        0.63                0.50                      0.36
#>     V9     2        0.56                0.38                      0.22
#>    V10     2        0.50                0.50                      0.35
#>    V11     2        0.47                0.47                      0.32
#>    V12     2        0.46                0.46                      0.31
#>    V13     2        0.41                0.38                      0.22
#>    V14     2        0.31                0.40                      0.25
#>    V15     2        0.32                0.39                      0.24
#>   alpha_removed  flag
#>            0.67      
#>            0.66      
#>            0.65      
#>            0.65      
#>            0.67      
#>            0.66      
#>            0.65      
#>            0.65      
#>            0.67      
#>            0.65      
#>            0.65      
#>            0.66      
#>            0.67      
#>            0.66      
#>            0.67      
#> 
#> Flagging thresholds: difficulty in [0.1, 0.95], discrimination >= 0.2
#> 0 of 15 item(s) flagged.
#> 
#> Test-Level Reliability Summary
#>   n_examinee  n_item  alpha  alpha_std   sem  mean_difficulty
#>          300      15   0.68       0.68  1.65              0.6
#>   mean_discrimination_raw  mean_discrimination_corrected
#>                      0.42                           0.28
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
