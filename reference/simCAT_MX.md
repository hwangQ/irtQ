# Simulated Mixed-Item Format CAT Data

A simulated dataset for computerized adaptive testing (CAT), containing
an item pool, sparse response data, and examinee ability estimates. The
item pool includes both dichotomous and polytomous items.

## Usage

``` r
simCAT_MX
```

## Format

A list of length three:

- item.prm:

  A data frame in item metadata format consisting of 200 dichotomous
  items and 30 polytomous items.

  - Dichotomous items: Calibrated using the IRT 2PL model.

  - Polytomous items: Calibrated using the Generalized Partial Credit
    Model (GPCM), with four score categories (0, 1, 2, 3).

- res.dat:

  A matrix of item responses from 30,000 examinees (rows) to 230 items
  (columns). `NA` marks an item that was not administered to the
  examinee. The columns are named `Item.dc.1` to `Item.dc.200` and
  `Item.py.1` to `Item.py.30`, not by the `id` column of `item.prm`
  (`V1` to `V230`); the columns follow the row order of `item.prm`.

- score:

  A numeric vector of ability estimates for the 30,000 examinees.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# structure of the data
str(simCAT_MX, max.level = 1)
#> List of 3
#>  $ item.prm:'data.frame':    230 obs. of  7 variables:
#>  $ res.dat : num [1:30000, 1:230] NA NA NA NA NA NA NA 1 NA NA ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ score   : num [1:30000] -0.303 -0.672 -0.735 1.769 -0.91 ...

# \donttest{
# item fit of three dichotomous items and two polytomous items
loc <- c(1:3, 201:202)
x <- simCAT_MX$item.prm[loc, ]
data <- simCAT_MX$res.dat[, loc]
irtfit(
  x = x, score = simCAT_MX$score, data = data, group.method = "equal.freq",
  n.width = 10, loc.theta = "average", range.score = c(-4, 4), D = 1
)
#> 
#> Call:
#> irtfit.default(x = x, score = simCAT_MX$score, data = data, group.method = "equal.freq", 
#>     n.width = 10, loc.theta = "average", range.score = c(-4, 
#>         4), D = 1)
#> 
#> Significance level for chi-square fit statistic: 0.05 
#> 
#> Item fit statistics: 
#>      id        X2        G2  df.X2  df.G2  crit.val.X2  crit.val.G2  p.X2  p.G2
#> 1    V1    86.094    88.522      8     10        15.51        18.31     0     0
#> 2    V2    85.780    86.229      8     10        15.51        18.31     0     0
#> 3    V3   201.658   179.886      8     10        15.51        18.31     0     0
#> 4  V201   910.867  1154.259     26     30        38.89        43.77     0     0
#> 5  V202  1480.850  1868.760     26     30        38.89        43.77     0     0
#>    outfit  infit      N  overSR.prop
#> 1   0.940  0.945   2558         0.60
#> 2   1.018  1.016   2018         0.60
#> 3   1.124  1.090  11041         0.50
#> 4   0.776  0.807   7384         0.75
#> 5   0.725  0.760  13254         0.80
#> 
#> Caution is needed in interpreting infit and outfit statistics for non-Rasch models. 
# }
```
