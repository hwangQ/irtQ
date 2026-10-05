# Simulated Single-Item Format CAT Data

A simulated dataset containing an item pool, sparse response data, and
examinee ability estimates, designed for single-item computerized
adaptive testing (CAT).

## Usage

``` r
simCAT_DC
```

## Format

A list of length three:

- item.prm:

  A data frame in item metadata format containing 100 dichotomous items.

  - Items 1-90: Generated and calibrated under the IRT 2PL model.

  - Items 91-100: Generated under the IRT 3PL model but calibrated using
    the 2PL model.

- res.dat:

  A matrix of item responses from 10,000 examinees (rows) to 100 items
  (columns). `NA` marks an item that was not administered to the
  examinee. The columns have no names; they follow the row order of
  `item.prm`, whose `id` values are `V1` to `V100`.

- score:

  A numeric vector of ability estimates for the 10,000 examinees.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# structure of the data
str(simCAT_DC, max.level = 1)
#> List of 3
#>  $ item.prm:'data.frame':    100 obs. of  6 variables:
#>  $ res.dat : num [1:10000, 1:100] NA NA NA NA NA NA NA NA NA NA ...
#>  $ score   : num [1:10000] 0.846 -0.914 -0.583 -0.314 -2.036 ...

# \donttest{
# item fit of the first five items, using the ability estimates as scores
x <- simCAT_DC$item.prm[1:5, ]
data <- simCAT_DC$res.dat[, 1:5]
irtfit(
  x = x, score = simCAT_DC$score, data = data, group.method = "equal.freq",
  n.width = 10, loc.theta = "average", range.score = c(-4, 4), D = 1
)
#> 
#> Call:
#> irtfit.default(x = x, score = simCAT_DC$score, data = data, group.method = "equal.freq", 
#>     n.width = 10, loc.theta = "average", range.score = c(-4, 
#>         4), D = 1)
#> 
#> Significance level for chi-square fit statistic: 0.05 
#> 
#> Item fit statistics: 
#>    id      X2      G2  df.X2  df.G2  crit.val.X2  crit.val.G2   p.X2   p.G2
#> 1  V1  18.570  20.901      8     10        15.51        18.31  0.017  0.022
#> 2  V2  28.242  29.106      8     10        15.51        18.31  0.000  0.001
#> 3  V3  16.746  17.621      8     10        15.51        18.31  0.033  0.062
#> 4  V4  14.266  14.650      8     10        15.51        18.31  0.075  0.145
#> 5  V5  50.441  55.780      8     10        15.51        18.31  0.000  0.000
#>    outfit  infit     N  overSR.prop
#> 1   0.892  0.945  1629          0.2
#> 2   0.955  0.961  2739          0.2
#> 3   0.939  0.963  5999          0.1
#> 4   0.946  0.964  5999          0.1
#> 5   0.889  0.922  5071          0.4
#> 
#> Caution is needed in interpreting infit and outfit statistics for non-Rasch models. 
# }
```
