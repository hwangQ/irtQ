# Traditional IRT Item Fit Statistics

This function computes traditional IRT item fit statistics, including
the \\\chi^{2}\\ fit statistic (e.g., Bock, 1960; Yen, 1981), the
log-likelihood ratio \\\chi^{2}\\ fit statistic (\\G^{2}\\; McKinley &
Mills, 1985), and the infit and outfit statistics (Ames et al., 2015).
It also returns contingency tables used to compute the \\\chi^{2}\\ and
\\G^{2}\\ statistics.

## Usage

``` r
irtfit(x, ...)

# Default S3 method
irtfit(
  x,
  score,
  data,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  D = 1,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  pcm.loc = NULL,
  ...
)

# S3 method for class 'est_item'
irtfit(
  x,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  pcm.loc = NULL,
  ...
)

# S3 method for class 'est_irt'
irtfit(
  x,
  score,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  pcm.loc = NULL,
  ...
)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.); or an object of class `est_irt`
  obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md), or
  `est_item` from
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

  See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
  or [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- ...:

  Further arguments passed to or from other methods.

- score:

  A numeric vector containing examinees' ability estimates (theta
  values).

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- group.method:

  A character string specifying the method used to group examinees along
  the ability scale when computing the \\\chi^{2}\\ and \\G^{2}\\ fit
  statistics. Available options are:

  - `"equal.width"`: Divides the ability scale into intervals of equal
    width.

  - `"equal.freq"`: Divides the examinees into groups with
    (approximately) equal numbers of examinees.

  Note that `"equal.freq"` does not guarantee exactly equal group sizes
  due to ties in ability estimates. Default is `"equal.width"`. The
  number of groups and the range of the ability scale are controlled by
  the `n.width` and `range.score` arguments, respectively.

- n.width:

  An integer specifying the number of intervals (groups) into which the
  ability scale is divided for computing the fit statistics. Default is
  10.

- loc.theta:

  A character string indicating the point on the ability scale at which
  the expected category probabilities are calculated for each group.
  Available options are:

  - `"average"`: Uses the average ability estimate of examinees within
    each group.

  - `"middle"`: Uses the midpoint of each group's ability interval.

  Default is `"average"`.

- range.score:

  A numeric vector of length two specifying the lower and upper bounds
  of the ability scale. Ability estimates below the lower bound or above
  the upper bound are truncated to the respective bound. If `NULL`, the
  observed minimum and maximum of the `score` vector are used. Note that
  this range restriction is independent of the grouping method specified
  in `group.method`. Default is `NULL`.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for the
  hypothesis tests of the \\\chi^{2}\\ and \\G^{2}\\ item fit
  statistics. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.
  See **Details** below.

- overSR:

  A numeric threshold used to identify ability groups (intervals) whose
  standardized residuals exceed the specified value. This is used to
  compute the proportion of misfitting groups per item. Default is 2.

- min.collapse:

  An integer specifying the minimum expected frequency required for a
  cell in the contingency table. Neighboring groups will be merged if
  any expected cell frequency falls below this threshold when computing
  the \\\chi^{2}\\ and \\G^{2}\\ statistics. Default is 1.

- pcm.loc:

  A vector of integers indicating the locations (indices) of partial
  credit model (PCM) items for which slope parameters are fixed.

## Value

This function returns an object of class `irtfit`, which includes the
following components:

- fit_stat:

  A data frame containing the results of three IRT item fit statistics -
  \\\chi^{2}\\, \\G^{2}\\, infit, and outfit - for all evaluated items.
  Each row corresponds to one item, and the columns include: the item
  ID; \\\chi^{2}\\ statistic; \\G^{2}\\ statistic; degrees of freedom
  for \\\chi^{2}\\ and \\G^{2}\\; critical values and p-values for both
  statistics; outfit and infit values; the number of examinees used to
  compute these statistics; and the proportion of ability groups (prior
  to cell collapsing) that have standardized residuals greater than the
  threshold specified in the `overSR` argument.

- contingency.fitstat:

  A list of contingency tables used to compute the \\\chi^{2}\\ and
  \\G^{2}\\ fit statistics for all items. Note that the cell-collapsing
  strategy is applied to these tables to ensure sufficient expected
  frequencies.

- contingency.plot:

  A list of contingency tables used to generate raw and standardized
  residual plots (Hambleton et al., 1991) via the
  [`plot.irtfit()`](https://hwangQ.github.io/irtQ/reference/plot.irtfit.md).
  Note that these tables are based on the original, uncollapsed
  groupings.

- individual.info:

  A list of data frames containing individual residuals and
  corresponding variance values. This information is used to compute
  infit and outfit statistics.

- item_df:

  A data frame containing the item metadata provided in the argument
  `x`.

- ancillary:

  A list of ancillary information used during the item fit analysis.

## Details

To compute the \\\chi^2\\ and \\G^2\\ item fit statistics, the
`group.method` argument determines how the ability scale is divided into
groups:

- `"equal.width"`: Examinees are grouped based on intervals of equal
  width a long the ability scale.

- `"equal.freq"`: Examinees are grouped such that each group contains
  (approximately) the same number of individuals.

Note that `"equal.freq"` does not guarantee *exactly* equal frequencies
across all groups, since grouping is based on quantiles.

When dividing the ability scale into intervals to compute the \\\chi^2\\
and \\G^2\\ fit statistics, the intervals should be:

- **Wide enough** to ensure that each group contains a sufficient number
  of examinees (to avoid unstable estimates),

- **Narrow enough** to ensure that examinees within each group are
  relatively homogeneous in ability (Hambleton et al., 1991).

If you want to divide the ability scale into a number of groups other
than the default of 10, specify the desired number using the `n.width`
argument. For reference:

- Yen (1981) used 10 fixed-width groups,

- Bock (1960) allowed for flexibility in the number of groups.

Regarding degrees of freedom (*df*):

- The \\\chi^2\\ statistic is approximately chi-square distributed with
  degrees of freedom equal to the number of ability groups minus the
  number of item parameters (Ames et al., 2015).

- The \\G^2\\ statistic is approximately chi-square distributed with
  degrees of freedom equal to the number of ability groups (Ames et al.,
  2015; Muraki & Bock, 2003).

  Note that if `"DRM"` is specified for an item in the item metadata
  set, the item is treated as a `"3PLM"` when computing the degrees of
  freedom for the \\\chi^{2}\\ fit statistic.

Note that infit and outfit statistics should be interpreted with caution
when applied to non-Rasch models. The returned object - particularly the
contingency tables - can be passed to
[`plot.irtfit()`](https://hwangQ.github.io/irtQ/reference/plot.irtfit.md)
to generate raw and standardized residual plots (Hambleton et al.,
1991).

## Methods (by class)

- `irtfit(default)`: Default method for computing traditional IRT item
  fit statistics using a data frame `x` that contains item metadata.

- `irtfit(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

- `irtfit(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

## References

Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional Module on
Item-Fit Statistics for Item Response Theory Models. *Educational
Measurement: Issues and Practice, 34*(3), 39-48.

Bock, R.D. (1960), *Methods and applications of optimal scaling*. Chapel
Hill, NC: L.L. Thurstone Psychometric Laboratory.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).*Fundamentals
of item response theory*. Newbury Park, CA: Sage.

McKinley, R., & Mills, C. (1985). A comparison of several
goodness-of-fit statistics. *Applied Psychological Measurement, 9*,
49-57.

Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test
scoring for rating scale data (Computer Software). Chicago, IL:
Scientific Software International. URL http://www.ssicentral.com

Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric
procedure for assessing goodness-of-fit in item response theory.
*Applied Measurement in Education, 21*(1), 22-40.

Yen, W. M. (1981). Using simulation results to choose a latent trait
model. *Applied Psychological Measurement, 5*, 245-262.

## See also

[`plot.irtfit()`](https://hwangQ.github.io/irtQ/reference/plot.irtfit.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
## Example 1
## Use the simulated CAT data
# Identify items with more than 10,000 responses
over10000 <- which(colSums(simCAT_MX$res.dat, na.rm = TRUE) > 10000)

# Select items with more than 10,000 responses
x <- simCAT_MX$item.prm[over10000, ]

# Extract response data for the selected items
data <- simCAT_MX$res.dat[, over10000]

# Extract examinees' ability estimates
score <- simCAT_MX$score

# Compute item fit statistics
fit1 <- irtfit(
  x = x, score = score, data = data, group.method = "equal.width",
  n.width = 10, loc.theta = "average", range.score = NULL, D = 1, alpha = 0.05,
  missing = NA, overSR = 2
)

# View the fit statistics
fit1$fit_stat
#>      id       X2       G2 df.X2 df.G2 crit.val.X2 crit.val.G2 p.X2 p.G2 outfit
#> 1   V37  157.239  137.111     7     9       14.07       16.92    0    0  1.129
#> 2   V40  134.853  121.555     8    10       15.51       18.31    0    0  1.087
#> 3   V59  182.470  166.750     8    10       15.51       18.31    0    0  1.087
#> 4   V63  172.870  153.429     7     9       14.07       16.92    0    0  1.091
#> 5   V66  216.148  191.918     8    10       15.51       18.31    0    0  1.196
#> 6   V74  158.851  145.358     8    10       15.51       18.31    0    0  1.105
#> 7   V75  251.708  226.279     8    10       15.51       18.31    0    0  1.161
#> 8   V84  220.177  190.837     8    10       15.51       18.31    0    0  1.099
#> 9   V89  216.934  190.497     8    10       15.51       18.31    0    0  1.151
#> 10  V97  162.835  146.236     8    10       15.51       18.31    0    0  1.086
#> 11 V107  212.665  187.084     8    10       15.51       18.31    0    0  1.112
#> 12 V150  241.878  210.551     8    10       15.51       18.31    0    0  1.112
#> 13 V157  166.947  148.510     8    10       15.51       18.31    0    0  1.139
#> 14 V166  161.125  135.352     8    10       15.51       18.31    0    0  1.087
#> 15 V172  257.947  232.775     8    10       15.51       18.31    0    0  1.109
#> 16 V176  217.621  191.815     8    10       15.51       18.31    0    0  1.112
#> 17 V185  159.682  135.802     7     9       14.07       16.92    0    0  1.184
#> 18 V186  190.901  170.650     8    10       15.51       18.31    0    0  1.118
#> 19 V187  217.861  192.502     8    10       15.51       18.31    0    0  1.090
#> 20 V201  707.671  998.945    26    30       38.89       43.77    0    0  0.776
#> 21 V202 1236.908 1746.382    20    24       31.41       36.42    0    0  0.725
#> 22 V211  935.227 1280.437    23    27       35.17       40.11    0    0  0.769
#> 23 V214 2157.262 2669.312    26    30       38.89       43.77    0    0  0.702
#> 24 V217 1322.558 1758.534    17    21       27.59       32.67    0    0  0.683
#> 25 V220 2091.594 3119.452    26    30       38.89       43.77    0    0  0.693
#> 26 V222 1006.092 1360.961    20    24       31.41       36.42    0    0  0.657
#> 27 V223 1090.650 1372.033    23    27       35.17       40.11    0    0  0.699
#> 28 V226  988.500 1425.971    23    27       35.17       40.11    0    0  0.765
#>    infit     N overSR.prop
#> 1  1.078 17939       0.500
#> 2  1.058 25834       0.800
#> 3  1.064 21249       0.800
#> 4  1.062 19778       0.700
#> 5  1.092 26369       0.900
#> 6  1.068 25948       0.800
#> 7  1.093 24606       0.700
#> 8  1.066 22926       0.800
#> 9  1.081 25991       0.700
#> 10 1.058 22275       0.800
#> 11 1.072 22510       0.800
#> 12 1.075 21040       0.700
#> 13 1.074 25736       0.800
#> 14 1.057 18980       0.600
#> 15 1.077 23044       0.700
#> 16 1.072 23353       0.800
#> 17 1.095 15965       0.500
#> 18 1.072 24494       0.800
#> 19 1.065 20929       0.800
#> 20 0.807  7384       0.725
#> 21 0.760 13254       0.750
#> 22 0.793  8375       0.850
#> 23 0.728  8348       0.825
#> 24 0.722 11879       0.600
#> 25 0.722  8159       0.825
#> 26 0.725 16423       0.750
#> 27 0.754 17179       0.800
#> 28 0.790  9487       0.650

# View the contingency tables used to compute fit statistics
fit1$contingency.fitstat
#> $V37
#>   total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1  obs.prop.0 obs.prop.1
#> 1   814        642        172  660.044874   153.9551 0.788697789  0.2113022
#> 2  1537       1001        536 1038.234261   498.7657 0.651268705  0.3487313
#> 3  2630       1388       1242 1400.800535  1229.1995 0.527756654  0.4722433
#> 4  4386       1667       2719 1663.315681  2722.6843 0.380072959  0.6199270
#> 5  3834       1081       2753  943.666706  2890.3333 0.281950965  0.7180490
#> 6  3107        616       2491  470.658865  2636.3411 0.198261989  0.8017380
#> 7  1147        145       1002   97.867617  1049.1324 0.126416739  0.8735833
#> 8   381         47        334   18.705807   362.2942 0.123359580  0.8766404
#> 9   103          1        102    2.614117   100.3859 0.009708738  0.9902913
#>   exp.prob.0 exp.prob.1     raw.rsd.0     raw.rsd.1
#> 1 0.81086594  0.1891341 -0.0221681499  0.0221681499
#> 2 0.67549399  0.3245060 -0.0242252840  0.0242252840
#> 3 0.53262378  0.4673762 -0.0048671235  0.0048671235
#> 4 0.37923294  0.6207671  0.0008400179 -0.0008400179
#> 5 0.24613112  0.7538689  0.0358198473 -0.0358198473
#> 6 0.15148338  0.8485166  0.0467786079 -0.0467786079
#> 7 0.08532486  0.9146751  0.0410918769 -0.0410918769
#> 8 0.04909661  0.9509034  0.0742629740 -0.0742629740
#> 9 0.02537977  0.9746202 -0.0156710371  0.0156710371
#> 
#> $V40
#>    total obs.freq.0 obs.freq.1 exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    269        225         44  249.56106   19.43894 0.83643123  0.1635688
#> 2   1097        906        191  952.18873  144.81127 0.82588879  0.1741112
#> 3   2264       1693        571 1785.98317  478.01683 0.74779152  0.2522085
#> 4   4469       2929       1540 3033.53507 1435.46493 0.65540389  0.3445961
#> 5   5123       2632       2491 2706.50075 2416.49925 0.51376147  0.4862385
#> 6   4778       1881       2897 1823.57709 2954.42291 0.39367936  0.6063206
#> 7   3957       1098       2859  999.18951 2957.81049 0.27748294  0.7225171
#> 8   2654        484       2170  414.03275 2239.96725 0.18236624  0.8176338
#> 9    768         98        670   66.98114  701.01886 0.12760417  0.8723958
#> 10   455         26        429   23.89007  431.10993 0.05714286  0.9428571
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.92773626 0.07226374 -0.091305037  0.091305037
#> 2  0.86799337 0.13200663 -0.042104586  0.042104586
#> 3  0.78886182 0.21113818 -0.041070305  0.041070305
#> 4  0.67879505 0.32120495 -0.023391154  0.023391154
#> 5  0.52830388 0.47169612 -0.014542407  0.014542407
#> 6  0.38166117 0.61833883  0.012018190 -0.012018190
#> 7  0.25251188 0.74748812  0.024971062 -0.024971062
#> 8  0.15600330 0.84399670  0.026362944 -0.026362944
#> 9  0.08721503 0.91278497  0.040389136 -0.040389136
#> 10 0.05250564 0.94749436  0.004637215 -0.004637215
#> 
#> $V59
#>    total obs.freq.0 obs.freq.1 exp.freq.0  exp.freq.1 obs.prop.0 obs.prop.1
#> 1    110         88         22  100.54210    9.457899  0.8000000  0.2000000
#> 2    494        389        105  424.47426   69.525740  0.7874494  0.2125506
#> 3   1645       1215        430 1299.36185  345.638149  0.7386018  0.2613982
#> 4   3933       2605       1328 2718.12089 1214.879109  0.6623443  0.3376557
#> 5   4239       2416       1823 2431.36905 1807.630949  0.5699457  0.4300543
#> 6   4744       2183       2561 2136.01650 2607.983504  0.4601602  0.5398398
#> 7   3568       1325       2243 1171.24358 2396.756418  0.3713565  0.6286435
#> 8   1969        566       1403  450.92230 1518.077696  0.2874556  0.7125444
#> 9    438        100        338   66.23718  371.762820  0.2283105  0.7716895
#> 10   109         20         89    9.88696   99.113040  0.1834862  0.8165138
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.91401910  0.0859809 -0.114019098  0.114019098
#> 2  0.85925964  0.1407404 -0.071810243  0.071810243
#> 3  0.78988562  0.2101144 -0.051283800  0.051283800
#> 4  0.69110625  0.3088937 -0.028761986  0.028761986
#> 5  0.57357137  0.4264286 -0.003625631  0.003625631
#> 6  0.45025643  0.5497436  0.009903774 -0.009903774
#> 7  0.32826334  0.6717367  0.043093167 -0.043093167
#> 8  0.22901082  0.7709892  0.058444741 -0.058444741
#> 9  0.15122644  0.8487736  0.077084063 -0.077084063
#> 10 0.09070605  0.9092939  0.092780185 -0.092780185
#> 
#> $V63
#>   total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1   748        623        125  647.359334  100.64067  0.8328877  0.1671123
#> 2  1331       1008        323 1035.496631  295.50337  0.7573253  0.2426747
#> 3  2540       1663        877 1721.385605  818.61440  0.6547244  0.3452756
#> 4  4481       2439       2042 2438.934149 2042.06585  0.5442981  0.4557019
#> 5  4133       1759       2374 1687.815345 2445.18466  0.4255988  0.5744012
#> 6  3855       1228       2627 1120.874746 2734.12525  0.3185473  0.6814527
#> 7  1956        513       1443  377.938553 1578.06145  0.2622699  0.7377301
#> 8   633        146        487   78.386234  554.61377  0.2306477  0.7693523
#> 9   101         16         85    7.415968   93.58403  0.1584158  0.8415842
#>   exp.prob.0 exp.prob.1     raw.rsd.0     raw.rsd.1
#> 1 0.86545366  0.1345463 -3.256595e-02  3.256595e-02
#> 2 0.77798395  0.2220161 -2.065863e-02  2.065863e-02
#> 3 0.67771087  0.3222891 -2.298646e-02  2.298646e-02
#> 4 0.54428345  0.4557165  1.469553e-05 -1.469553e-05
#> 5 0.40837536  0.5916246  1.722348e-02 -1.722348e-02
#> 6 0.29075869  0.7092413  2.778865e-02 -2.778865e-02
#> 7 0.19322012  0.8067799  6.904982e-02 -6.904982e-02
#> 8 0.12383291  0.8761671  1.068148e-01 -1.068148e-01
#> 9 0.07342543  0.9265746  8.499042e-02 -8.499042e-02
#> 
#> $V66
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    714        660         54  684.298856   29.70114 0.92436975 0.07563025
#> 2   1602       1387        215 1441.152929  160.84707 0.86579276 0.13420724
#> 3   2560       1979        581 2060.105629  499.89437 0.77304688 0.22695313
#> 4   4577       2877       1700 2993.539060 1583.46094 0.62857767 0.37142233
#> 5   5127       2269       2858 2270.443275 2856.55673 0.44255900 0.55744100
#> 6   4753       1372       3381 1241.490009 3511.50999 0.28865979 0.71134021
#> 7   3828        666       3162  515.668506 3312.33149 0.17398119 0.82601881
#> 8   2353        234       2119  152.320754 2200.67925 0.09944751 0.90055249
#> 9    585         34        551   15.905649  569.09435 0.05811966 0.94188034
#> 10   270          9        261    3.552967  266.44703 0.03333333 0.96666667
#>    exp.prob.0 exp.prob.1     raw.rsd.0     raw.rsd.1
#> 1  0.95840176 0.04159824 -0.0340320110  0.0340320110
#> 2  0.89959609 0.10040391 -0.0338033267  0.0338033267
#> 3  0.80472876 0.19527124 -0.0316818862  0.0316818862
#> 4  0.65403956 0.34596044 -0.0254618876  0.0254618876
#> 5  0.44284051 0.55715949 -0.0002815047  0.0002815047
#> 6  0.26120135 0.73879865  0.0274584453 -0.0274584453
#> 7  0.13470964 0.86529036  0.0392715502 -0.0392715502
#> 8  0.06473470 0.93526530  0.0347128116 -0.0347128116
#> 9  0.02718914 0.97281086  0.0309305142 -0.0309305142
#> 10 0.01315914 0.98684086  0.0201741951 -0.0201741951
#> 
#> $V74
#>    total obs.freq.0 obs.freq.1 exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    262        236         26  250.67682   11.32318 0.90076336 0.09923664
#> 2    838        733        105  763.54729   74.45271 0.87470167 0.12529833
#> 3   1969       1575        394 1667.13159  301.86841 0.79989843 0.20010157
#> 4   4220       3029       1191 3154.21550 1065.78450 0.71777251 0.28222749
#> 5   5095       2952       2143 3043.63129 2051.36871 0.57939156 0.42060844
#> 6   4869       2178       2691 2122.31567 2746.68433 0.44731978 0.55268022
#> 7   4334       1318       3016 1232.88261 3101.11739 0.30410706 0.69589294
#> 8   3071        655       2416  525.77970 2545.22030 0.21328557 0.78671443
#> 9    887        117        770   82.12432  804.87568 0.13190530 0.86809470
#> 10   403         21        382   21.71724  381.28276 0.05210918 0.94789082
#>    exp.prob.0 exp.prob.1   raw.rsd.0   raw.rsd.1
#> 1  0.95678174 0.04321826 -0.05601838  0.05601838
#> 2  0.91115429 0.08884571 -0.03645262  0.03645262
#> 3  0.84668948 0.15331052 -0.04679106  0.04679106
#> 4  0.74744443 0.25255557 -0.02967192  0.02967192
#> 5  0.59737611 0.40262389 -0.01798455  0.01798455
#> 6  0.43588328 0.56411672  0.01143650 -0.01143650
#> 7  0.28446761 0.71553239  0.01963945 -0.01963945
#> 8  0.17120798 0.82879202  0.04207760 -0.04207760
#> 9  0.09258660 0.90741340  0.03931870 -0.03931870
#> 10 0.05388892 0.94611108 -0.00177974  0.00177974
#> 
#> $V75
#>    total obs.freq.0 obs.freq.1 exp.freq.0  exp.freq.1 obs.prop.0  obs.prop.1
#> 1    198        197          1  194.90459    3.095412  0.9949495 0.005050505
#> 2    556        518         38  536.13800   19.861999  0.9316547 0.068345324
#> 3   1254       1108        146 1164.94404   89.055955  0.8835726 0.116427432
#> 4   3346       2744        602 2902.30143  443.698567  0.8200837 0.179916318
#> 5   4800       3439       1361 3604.96323 1195.036774  0.7164583 0.283541667
#> 6   4828       2830       1998 2860.76808 1967.231922  0.5861640 0.413835957
#> 7   4440       1822       2618 1807.42060 2632.579401  0.4103604 0.589639640
#> 8   3457        994       2463  849.94075 2607.059249  0.2875325 0.712467457
#> 9   1109        210        899  143.92794  965.072058  0.1893598 0.810640216
#> 10   618         80        538   44.78631  573.213694  0.1294498 0.870550162
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.98436660 0.01563340  0.010582891 -0.010582891
#> 2  0.96427698 0.03572302 -0.032622304  0.032622304
#> 3  0.92898249 0.07101751 -0.045409924  0.045409924
#> 4  0.86739433 0.13260567 -0.047310650  0.047310650
#> 5  0.75103401 0.24896599 -0.034575672  0.034575672
#> 6  0.59253688 0.40746312 -0.006372841  0.006372841
#> 7  0.40707671 0.59292329  0.003283649 -0.003283649
#> 8  0.24586079 0.75413921  0.041671753 -0.041671753
#> 9  0.12978173 0.87021827  0.059578051 -0.059578051
#> 10 0.07246975 0.92753025  0.056980087 -0.056980087
#> 
#> $V84
#>    total obs.freq.0 obs.freq.1 exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    148        126         22  137.42177   10.57823  0.8513514  0.1486486
#> 2    667        534        133  583.41002   83.58998  0.8005997  0.1994003
#> 3   1945       1468        477 1564.90287  380.09713  0.7547558  0.2452442
#> 4   4143       2815       1328 2923.28176 1219.71824  0.6794593  0.3205407
#> 5   4394       2524       1870 2550.49492 1843.50508  0.5744197  0.4255803
#> 6   5181       2275       2906 2302.80495 2878.19505  0.4391044  0.5608956
#> 7   3748       1281       2467 1160.49271 2587.50729  0.3417823  0.6582177
#> 8   1929        518       1411  410.15413 1518.84587  0.2685329  0.7314671
#> 9    601        126        475   85.18522  515.81478  0.2096506  0.7903494
#> 10   170         40        130   14.22200  155.77800  0.2352941  0.7647059
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.92852544 0.07147456 -0.077174091  0.077174091
#> 2  0.87467769 0.12532231 -0.074077987  0.074077987
#> 3  0.80457731 0.19542269 -0.049821526  0.049821526
#> 4  0.70559540 0.29440460 -0.026136074  0.026136074
#> 5  0.58044946 0.41955054 -0.006029794  0.006029794
#> 6  0.44447113 0.55552887 -0.005366714  0.005366714
#> 7  0.30962986 0.69037014  0.032152425 -0.032152425
#> 8  0.21262526 0.78737474  0.055907657 -0.055907657
#> 9  0.14173913 0.85826087  0.067911454 -0.067911454
#> 10 0.08365885 0.91634115  0.151635272 -0.151635272
#> 
#> $V89
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    376        358         18  362.432617   13.56738 0.95212766 0.04787234
#> 2   1091        954        137 1002.634148   88.36585 0.87442713 0.12557287
#> 3   2171       1737        434 1844.427583  326.57242 0.80009212 0.19990788
#> 4   4349       3068       1281 3206.237171 1142.76283 0.70544953 0.29455047
#> 5   5102       2850       2252 2873.633453 2228.36655 0.55860447 0.44139553
#> 6   4852       1918       2934 1858.613969 2993.38603 0.39530091 0.60469909
#> 7   4193       1047       3146  955.855905 3237.14410 0.24970188 0.75029812
#> 8   2796        457       2339  348.234296 2447.76570 0.16344778 0.83655222
#> 9    767         96        671   45.798485  721.20152 0.12516297 0.87483703
#> 10   294         19        275    9.465732  284.53427 0.06462585 0.93537415
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.96391654 0.03608346 -0.011788876  0.011788876
#> 2  0.91900472 0.08099528 -0.044577588  0.044577588
#> 3  0.84957512 0.15042488 -0.049482996  0.049482996
#> 4  0.73723550 0.26276450 -0.031785967  0.031785967
#> 5  0.56323666 0.43676334 -0.004632194  0.004632194
#> 6  0.38306141 0.61693859  0.012239495 -0.012239495
#> 7  0.22796468 0.77203532  0.021737204 -0.021737204
#> 8  0.12454732 0.87545268  0.038900466 -0.038900466
#> 9  0.05971119 0.94028881  0.065451780 -0.065451780
#> 10 0.03219637 0.96780363  0.032429482 -0.032429482
#> 
#> $V97
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    572        485         87  511.987531   60.01247  0.8479021  0.1520979
#> 2   1239        986        253 1022.768395  216.23161  0.7958031  0.2041969
#> 3   2483       1762        721 1838.311851  644.68815  0.7096255  0.2903745
#> 4   4495       2757       1738 2782.912757 1712.08724  0.6133482  0.3866518
#> 5   4298       2083       2215 2083.695592 2214.30441  0.4846440  0.5153560
#> 6   4535       1707       2828 1622.134751 2912.86525  0.3764057  0.6235943
#> 7   3008        891       2117  735.214020 2272.78598  0.2962101  0.7037899
#> 8   1357        301       1056  218.740833 1138.25917  0.2218128  0.7781872
#> 9    249         53        196   25.020167  223.97983  0.2128514  0.7871486
#> 10    39          6         33    2.254606   36.74539  0.1538462  0.8461538
#>    exp.prob.0 exp.prob.1     raw.rsd.0     raw.rsd.1
#> 1  0.89508310  0.1049169 -0.0471809980  0.0471809980
#> 2  0.82547893  0.1745211 -0.0296758636  0.0296758636
#> 3  0.74035918  0.2596408 -0.0307337296  0.0307337296
#> 4  0.61911296  0.3808870 -0.0057647958  0.0057647958
#> 5  0.48480586  0.5151941 -0.0001618408  0.0001618408
#> 6  0.35769234  0.6423077  0.0187133956 -0.0187133956
#> 7  0.24441955  0.7555804  0.0517905519 -0.0517905519
#> 8  0.16119442  0.8388056  0.0606183988 -0.0606183988
#> 9  0.10048260  0.8995174  0.1123688056 -0.1123688056
#> 10 0.05781041  0.9421896  0.0960357451 -0.0960357451
#> 
#> $V107
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    147        119         28  135.452087   11.54791  0.8095238  0.1904762
#> 2    620        508        112  534.542775   85.45723  0.8193548  0.1806452
#> 3   1808       1353        455 1415.191511  392.80849  0.7483407  0.2516593
#> 4   4094       2605       1489 2722.474468 1371.52553  0.6362970  0.3637030
#> 5   4277       2210       2067 2247.527136 2029.47286  0.5167173  0.4832827
#> 6   4820       1913       2907 1857.437044 2962.56296  0.3968880  0.6031120
#> 7   3783       1141       2642  973.927861 2809.07214  0.3016125  0.6983875
#> 8   2267        495       1772  371.722870 1895.27713  0.2183502  0.7816498
#> 9    559         98        461   54.746938  504.25306  0.1753131  0.8246869
#> 10   135         19        116    7.180443  127.81956  0.1407407  0.8592593
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.92144277 0.07855723 -0.111918960  0.111918960
#> 2  0.86216577 0.13783423 -0.042810927  0.042810927
#> 3  0.78273867 0.21726133 -0.034397959  0.034397959
#> 4  0.66499132 0.33500868 -0.028694301  0.028694301
#> 5  0.52549150 0.47450850 -0.008774173  0.008774173
#> 6  0.38536038 0.61463962  0.011527584 -0.011527584
#> 7  0.25744855 0.74255145  0.044163928 -0.044163928
#> 8  0.16397127 0.83602873  0.054378972 -0.054378972
#> 9  0.09793728 0.90206272  0.077375782 -0.077375782
#> 10 0.05318847 0.94681153  0.087552276 -0.087552276
#> 
#> $V150
#>    total obs.freq.0 obs.freq.1 exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    280        225         55  252.09604   27.90396 0.80357143  0.1964286
#> 2    877        686        191  727.72118  149.27882 0.78221209  0.2177879
#> 3   2167       1524        643 1602.51519  564.48481 0.70327642  0.2967236
#> 4   4349       2601       1748 2654.70473 1694.29527 0.59806852  0.4019315
#> 5   4286       2049       2237 2002.84093 2283.15907 0.47806813  0.5219319
#> 6   4520       1606       2914 1509.46067 3010.53933 0.35530973  0.6446903
#> 7   2977        816       2161  653.46383 2323.53617 0.27410144  0.7258986
#> 8   1313        287       1026  181.86735 1131.13265 0.21858340  0.7814166
#> 9    232         48        184   19.20296  212.79704 0.20689655  0.7931034
#> 10    39          1         38    1.77352   37.22648 0.02564103  0.9743590
#>    exp.prob.0 exp.prob.1   raw.rsd.0   raw.rsd.1
#> 1  0.90034302 0.09965698 -0.09677159  0.09677159
#> 2  0.82978469 0.17021531 -0.04757261  0.04757261
#> 3  0.73950863 0.26049137 -0.03623221  0.03623221
#> 4  0.61041728 0.38958272 -0.01234875  0.01234875
#> 5  0.46729840 0.53270160  0.01076973 -0.01076973
#> 6  0.33395148 0.66604852  0.02135826 -0.02135826
#> 7  0.21950414 0.78049586  0.05459730 -0.05459730
#> 8  0.13851283 0.86148717  0.08007057 -0.08007057
#> 9  0.08277137 0.91722863  0.12412518 -0.12412518
#> 10 0.04547487 0.95452513 -0.01983384  0.01983384
#> 
#> $V157
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    573        510         63  539.603559   33.39644 0.89005236  0.1099476
#> 2   1507       1265        242 1319.518800  187.48120 0.83941606  0.1605839
#> 3   2532       1906        626 1974.430854  557.56915 0.75276461  0.2472354
#> 4   4573       2851       1722 2923.395659 1649.60434 0.62344194  0.3765581
#> 5   5126       2342       2784 2309.368158 2816.63184 0.45688646  0.5431135
#> 6   4731       1482       3249 1349.075807 3381.92419 0.31325301  0.6867470
#> 7   3701        708       2993  599.373605 3101.62640 0.19129965  0.8087004
#> 8   2189        265       1924  187.280219 2001.71978 0.12105984  0.8789402
#> 9    537         45        492   21.449704  515.55030 0.08379888  0.9162011
#> 10   267          7        260    5.562012  261.43799 0.02621723  0.9737828
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.94171651 0.05828349 -0.051664152  0.051664152
#> 2  0.87559310 0.12440690 -0.036177041  0.036177041
#> 3  0.77979102 0.22020898 -0.027026404  0.027026404
#> 4  0.63927305 0.36072695 -0.015831108  0.015831108
#> 5  0.45052051 0.54947949  0.006365947 -0.006365947
#> 6  0.28515659 0.71484341  0.028096426 -0.028096426
#> 7  0.16194910 0.83805090  0.029350553 -0.029350553
#> 8  0.08555515 0.91444485  0.035504697 -0.035504697
#> 9  0.03994358 0.96005642  0.043855300 -0.043855300
#> 10 0.02083151 0.97916849  0.005385721 -0.005385721
#> 
#> $V166
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    779        658        121  669.116761  109.88324  0.8446727  0.1553273
#> 2    863        646        217  675.691289  187.30871  0.7485516  0.2514484
#> 3   2396       1603        793 1644.067798  751.93220  0.6690317  0.3309683
#> 4   3797       2068       1729 2085.249368 1711.75063  0.5446405  0.4553595
#> 5   3700       1597       2103 1579.414109 2120.58589  0.4316216  0.5683784
#> 6   3979       1333       2646 1221.380479 2757.61952  0.3350088  0.6649912
#> 7   2277        591       1686  478.965498 1798.03450  0.2595520  0.7404480
#> 8    950        194        756  129.789957  820.21004  0.2042105  0.7957895
#> 9    213         49        164   18.553175  194.44683  0.2300469  0.7699531
#> 10    26          5         21    1.347068   24.65293  0.1923077  0.8076923
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.85894321  0.1410568 -0.014270553  0.014270553
#> 2  0.78295630  0.2170437 -0.034404738  0.034404738
#> 3  0.68617187  0.3138281 -0.017140150  0.017140150
#> 4  0.54918340  0.4508166 -0.004542894  0.004542894
#> 5  0.42686868  0.5731313  0.004752943 -0.004752943
#> 6  0.30695664  0.6930434  0.028052154 -0.028052154
#> 7  0.21034936  0.7896506  0.049202680 -0.049202680
#> 8  0.13662101  0.8633790  0.067589519 -0.067589519
#> 9  0.08710411  0.9128959  0.142942842 -0.142942842
#> 10 0.05181031  0.9481897  0.140497386 -0.140497386
#> 
#> $V172
#>    total obs.freq.0 obs.freq.1 exp.freq.0  exp.freq.1 obs.prop.0 obs.prop.1
#> 1     66         58          8   62.11728    3.882718  0.8787879  0.1212121
#> 2    417        331         86  373.90608   43.093921  0.7937650  0.2062350
#> 3   1332       1029        303 1113.35331  218.646690  0.7725225  0.2274775
#> 4   3639       2510       1129 2722.70719  916.292814  0.6897499  0.3102501
#> 5   4962       2957       2005 3079.63866 1882.361340  0.5959291  0.4040709
#> 6   4827       2378       2449 2324.55470 2502.445303  0.4926455  0.5073545
#> 7   4181       1484       2697 1433.87328 2747.126721  0.3549390  0.6450610
#> 8   2720        752       1968  622.87207 2097.127926  0.2764706  0.7235294
#> 9    683        155        528   93.16829  589.831706  0.2269400  0.7730600
#> 10   217         22        195   18.82221  198.177788  0.1013825  0.8986175
#>    exp.prob.0 exp.prob.1   raw.rsd.0   raw.rsd.1
#> 1   0.9411709 0.05882906 -0.06238307  0.06238307
#> 2   0.8966573 0.10334274 -0.10289227  0.10289227
#> 3   0.8358508 0.16414917 -0.06332831  0.06332831
#> 4   0.7482020 0.25179797 -0.05845210  0.05845210
#> 5   0.6206446 0.37935537 -0.02471557  0.02471557
#> 6   0.4815734 0.51842662  0.01107216 -0.01107216
#> 7   0.3429498 0.65705016  0.01198917 -0.01198917
#> 8   0.2289971 0.77100291  0.04747350 -0.04747350
#> 9   0.1364104 0.86358961  0.09052958 -0.09052958
#> 10  0.0867383 0.91326170  0.01464419 -0.01464419
#> 
#> $V176
#>    total obs.freq.0 obs.freq.1 exp.freq.0  exp.freq.1 obs.prop.0 obs.prop.1
#> 1     34         30          4   32.64428    1.355724  0.8823529  0.1176471
#> 2    253        208         45  234.86038   18.139617  0.8221344  0.1778656
#> 3    953        773        180  839.06830  113.931704  0.8111228  0.1888772
#> 4   3167       2387        780 2552.06750  614.932500  0.7537101  0.2462899
#> 5   4775       3145       1630 3285.34265 1489.657353  0.6586387  0.3413613
#> 6   4843       2649       2194 2652.26790 2190.732104  0.5469750  0.4530250
#> 7   4425       1761       2664 1747.81331 2677.186687  0.3979661  0.6020339
#> 8   3381        998       2383  884.31622 2496.683785  0.2951789  0.7048211
#> 9   1007        215        792  158.16623  848.833767  0.2135055  0.7864945
#> 10   515         70        445   51.38197  463.618027  0.1359223  0.8640777
#>    exp.prob.0 exp.prob.1     raw.rsd.0     raw.rsd.1
#> 1  0.96012578 0.03987422 -0.0777728360  0.0777728360
#> 2  0.92830191 0.07169809 -0.1061675203  0.1061675203
#> 3  0.88044942 0.11955058 -0.0693266487  0.0693266487
#> 4  0.80583123 0.19416877 -0.0521210925  0.0521210925
#> 5  0.68802987 0.31197013 -0.0293911303  0.0293911303
#> 6  0.54764978 0.45235022 -0.0006747668  0.0006747668
#> 7  0.39498606 0.60501394  0.0029800422 -0.0029800422
#> 8  0.26155463 0.73844537  0.0336243078 -0.0336243078
#> 9  0.15706677 0.84293323  0.0564386963 -0.0564386963
#> 10 0.09977082 0.90022918  0.0361515080 -0.0361515080
#> 
#> $V185
#>   total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1   840        615        225  639.165283   200.8347 0.73214286  0.2678571
#> 2  1540        878        662  912.443545   627.5565 0.57012987  0.4298701
#> 3  2613       1114       1499 1126.441351  1486.5586 0.42632989  0.5736701
#> 4  4208       1288       2920 1173.914686  3034.0853 0.30608365  0.6939163
#> 5  3287        640       2647  541.119390  2745.8806 0.19470642  0.8052936
#> 6  2362        316       2046  220.376072  2141.6239 0.13378493  0.8662151
#> 7   777         81        696   37.603521   739.3965 0.10424710  0.8957529
#> 8   253         17        236    6.565049   246.4350 0.06719368  0.9328063
#> 9    85          1         84    1.061403    83.9386 0.01176471  0.9882353
#>   exp.prob.0 exp.prob.1     raw.rsd.0     raw.rsd.1
#> 1 0.76091105  0.2390889 -0.0287681938  0.0287681938
#> 2 0.59249581  0.4075042 -0.0223659385  0.0223659385
#> 3 0.43109122  0.5689088 -0.0047613284  0.0047613284
#> 4 0.27897212  0.7210279  0.0271115291 -0.0271115291
#> 5 0.16462409  0.8353759  0.0300823273 -0.0300823273
#> 6 0.09330062  0.9066994  0.0404843050 -0.0404843050
#> 7 0.04839578  0.9516042  0.0558513240 -0.0558513240
#> 8 0.02594881  0.9740512  0.0412448655 -0.0412448655
#> 9 0.01248709  0.9875129 -0.0007223891  0.0007223891
#> 
#> $V186
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    426        372         54  395.970241   30.02976  0.8732394  0.1267606
#> 2   1230       1021        209 1064.755950  165.24405  0.8300813  0.1699187
#> 3   2495       1849        646 1943.114708  551.88529  0.7410822  0.2589178
#> 4   4490       2839       1651 2925.008568 1564.99143  0.6322940  0.3677060
#> 5   4430       2247       2183 2198.053318 2231.94668  0.5072235  0.4927765
#> 6   5172       1838       3334 1776.638854 3395.36115  0.3553751  0.6446249
#> 7   3699        946       2753  782.351535 2916.64847  0.2557448  0.7442552
#> 8   1849        328       1521  240.488195 1608.51181  0.1773932  0.8226068
#> 9    558         76        482   43.538047  514.46195  0.1362007  0.8637993
#> 10   145         16        129    5.996887  139.00311  0.1103448  0.8896552
#>    exp.prob.0 exp.prob.1   raw.rsd.0   raw.rsd.1
#> 1  0.92950761 0.07049239 -0.05626817  0.05626817
#> 2  0.86565524 0.13434476 -0.03557394  0.03557394
#> 3  0.77880349 0.22119651 -0.03772133  0.03772133
#> 4  0.65144957 0.34855043 -0.01915558  0.01915558
#> 5  0.49617456 0.50382544  0.01104891 -0.01104891
#> 6  0.34351099 0.65648901  0.01186410 -0.01186410
#> 7  0.21150352 0.78849648  0.04424127 -0.04424127
#> 8  0.13006392 0.86993608  0.04732926 -0.04732926
#> 9  0.07802517 0.92197483  0.05817554 -0.05817554
#> 10 0.04135784 0.95864216  0.06898699 -0.06898699
#> 
#> $V187
#>    total obs.freq.0 obs.freq.1  exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1    137        103         34  123.162022   13.83798  0.7518248  0.2481752
#> 2    568        445        123  475.142273   92.85773  0.7834507  0.2165493
#> 3   1766       1241        525 1343.242718  422.75728  0.7027180  0.2972820
#> 4   4066       2585       1481 2664.223102 1401.77690  0.6357600  0.3642400
#> 5   4257       2257       2000 2274.322306 1982.67769  0.5301856  0.4698144
#> 6   4671       1964       2707 1927.387137 2743.61286  0.4204667  0.5795333
#> 7   3358       1117       2241  996.377436 2361.62256  0.3326385  0.6673615
#> 8   1689        476       1213  346.536448 1342.46355  0.2818236  0.7181764
#> 9    350         89        261   47.241324  302.75868  0.2542857  0.7457143
#> 10    67         12         55    5.494661   61.50534  0.1791045  0.8208955
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.89899286  0.1010071 -0.147168044  0.147168044
#> 2  0.83651809  0.1634819 -0.053067382  0.053067382
#> 3  0.76061309  0.2393869 -0.057895084  0.057895084
#> 4  0.65524425  0.3447558 -0.019484285  0.019484285
#> 5  0.53425471  0.4657453 -0.004069135  0.004069135
#> 6  0.41262837  0.5873716  0.007838335 -0.007838335
#> 7  0.29671752  0.7032825  0.035920954 -0.035920954
#> 8  0.20517256  0.7948274  0.076651007 -0.076651007
#> 9  0.13497521  0.8650248  0.119310503 -0.119310503
#> 10 0.08200987  0.9179901  0.097094610 -0.097094610
#> 
#> $V201
#>    total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1     38         33          1          1          3   21.74482   11.54804
#> 2    271        170         92          3          6  133.47491   87.44788
#> 3    442        252        135         30         25  185.48923  145.04012
#> 4    934        442        354         57         81  323.56372  300.62245
#> 5   1782        559        540        202        481  474.30099  535.05881
#> 6   1282        178        421        150        533  250.91952  341.75891
#> 7   1088        102        203        152        631  155.24163  250.08344
#> 8    760          0        114         79        567   72.59432  141.37026
#> 9    584          0          0         54        530   34.85516   83.03033
#> 10   203          0          0          0        203    6.28567   19.40046
#>    exp.freq.2 exp.freq.3 obs.prop.0 obs.prop.1 obs.prop.2 obs.prop.3 exp.prob.0
#> 1    1.911186   2.795956  0.8684211 0.02631579 0.02631579 0.07894737 0.57223205
#> 2   17.854230  32.222982  0.6273063 0.33948339 0.01107011 0.02214022 0.49252733
#> 3   35.342701  76.127947  0.5701357 0.30542986 0.06787330 0.05656109 0.41965890
#> 4   87.041260 222.772564  0.4732334 0.37901499 0.06102784 0.08672377 0.34642797
#> 5  188.101033 584.539168  0.3136925 0.30303030 0.11335578 0.26992144 0.26616217
#> 6  145.059911 544.261660  0.1388456 0.32839314 0.11700468 0.41575663 0.19572506
#> 7  125.546273 557.128657  0.0937500 0.18658088 0.13970588 0.57996324 0.14268532
#> 8   85.793742 460.241679  0.0000000 0.15000000 0.10394737 0.74605263 0.09551884
#> 9   61.638014 404.476505  0.0000000 0.00000000 0.09246575 0.90753425 0.05968348
#> 10  18.660136 158.653732  0.0000000 0.00000000 0.00000000 1.00000000 0.03096389
#>    exp.prob.1 exp.prob.2 exp.prob.3   raw.rsd.0    raw.rsd.1    raw.rsd.2
#> 1  0.30389579 0.05029437  0.0735778  0.29618900 -0.277579996 -0.023978577
#> 2  0.32268591 0.06588277  0.1189040  0.13477894  0.016797486 -0.054812657
#> 3  0.32814506 0.07996086  0.1722352  0.15047684 -0.022715198 -0.012087560
#> 4  0.32186558 0.09319193  0.2385145  0.12680544  0.057149408 -0.032164090
#> 5  0.30025747 0.10555614  0.3280242  0.04753031  0.002772832  0.007799645
#> 6  0.26658261 0.11315126  0.4245411 -0.05687950  0.061810524  0.003853424
#> 7  0.22985610 0.11539180  0.5120668 -0.04893532 -0.043275218  0.024314087
#> 8  0.18601350 0.11288650  0.6055812 -0.09551884 -0.036013496 -0.008939134
#> 9  0.14217522 0.10554454  0.6925968 -0.05968348 -0.142175216 -0.013078790
#> 10 0.09556878 0.09192185  0.7815455 -0.03096389 -0.095568783 -0.091921852
#>       raw.rsd.3
#> 1   0.005369571
#> 2  -0.096763772
#> 3  -0.115674087
#> 4  -0.151790754
#> 5  -0.058102788
#> 6  -0.008784446
#> 7   0.067896455
#> 8   0.140471474
#> 9   0.214937491
#> 10  0.218454524
#> 
#> $V202
#>   total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1   775        775          0          0          0 559.447133  185.03778
#> 2   770        550        220          0          0 468.824055  237.01445
#> 3  1457        840        524         93          0 719.826935  519.40759
#> 4  2171       1016        910        236          9 786.348925  823.81848
#> 5  3525        647       1463        899        516 811.352068 1253.13473
#> 6  2698        169        641        875       1013 331.323585  762.47827
#> 7  1486         33        171        387        895  86.799790  290.99760
#> 8   372          1         18         48        305   8.596281   42.25366
#>   exp.freq.2 exp.freq.3  obs.prop.0 obs.prop.1 obs.prop.2  obs.prop.3
#> 1   26.24556   4.269527 1.000000000  0.0000000 0.00000000 0.000000000
#> 2   51.38471  12.776778 0.714285714  0.2857143 0.00000000 0.000000000
#> 3  160.72471  57.040762 0.576527111  0.3596431 0.06382979 0.000000000
#> 4  370.11946 190.713138 0.467987103  0.4191617 0.10870567 0.004145555
#> 5  830.00430 630.508904 0.183546099  0.4150355 0.25503546 0.146382979
#> 6  752.48309 851.715059 0.062638992  0.2375834 0.32431431 0.375463306
#> 7  418.36399 689.838615 0.022207268  0.1150740 0.26043069 0.602288022
#> 8   90.68049 230.469575 0.002688172  0.0483871 0.12903226 0.819892473
#>   exp.prob.0 exp.prob.1 exp.prob.2  exp.prob.3   raw.rsd.0   raw.rsd.1
#> 1 0.72186727  0.2387584 0.03386524 0.005509067  0.27813273 -0.23875842
#> 2 0.60886241  0.3078110 0.06673339 0.016593218  0.10542330 -0.02209669
#> 3 0.49404731  0.3564911 0.11031209 0.039149459  0.08247980  0.00315196
#> 4 0.36220586  0.3794650 0.17048340 0.087845757  0.10578124  0.03969669
#> 5 0.23017080  0.3554992 0.23546221 0.178867774 -0.04662470  0.05953625
#> 6 0.12280340  0.2826087 0.27890404 0.315683862 -0.06016441 -0.04502530
#> 7 0.05841170  0.1958261 0.28153701 0.464225178 -0.03620443 -0.08075209
#> 8 0.02310828  0.1135851 0.24376475 0.619541869 -0.02042011 -0.06519800
#>     raw.rsd.2    raw.rsd.3
#> 1 -0.03386524 -0.005509067
#> 2 -0.06673339 -0.016593218
#> 3 -0.04648230 -0.039149459
#> 4 -0.06177773 -0.083700202
#> 5  0.01957325 -0.032484795
#> 6  0.04541027  0.059779444
#> 7 -0.02110632  0.138062843
#> 8 -0.11473249  0.200350604
#> 
#> $V211
#>   total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1   135        101         22         12          0   65.69728   43.98025
#> 2   369        234         99         28          8  146.94513  124.79548
#> 3   948        380        370        148         50  290.36831  314.35598
#> 4  1392        415        499        296        182  318.05316  429.25679
#> 5  2142        218        652        739        533  344.68686  580.71800
#> 6  1299        111        213        416        559  136.64613  290.23285
#> 7   914          0        146        232        536   60.21927  160.32232
#> 8   837          0          0        211        626   32.38527  108.94207
#> 9   339          0          0          0        339    6.25990   28.42711
#>   exp.freq.2 exp.freq.3 obs.prop.0 obs.prop.1 obs.prop.2 obs.prop.3 exp.prob.0
#> 1   17.43195   7.890518 0.74814815  0.1629630 0.08888889 0.00000000 0.48664648
#> 2   62.15991  35.099478 0.63414634  0.2682927 0.07588076 0.02168022 0.39822529
#> 3  199.60071 143.674999 0.40084388  0.3902954 0.15611814 0.05274262 0.30629569
#> 4  339.78365 304.906407 0.29813218  0.3584770 0.21264368 0.13074713 0.22848646
#> 5  573.81736 642.777779 0.10177404  0.3043884 0.34500467 0.24883287 0.16091824
#> 6  361.54607 510.574950 0.08545035  0.1639723 0.32024634 0.43033102 0.10519332
#> 7  250.33427 443.124141 0.00000000  0.1597374 0.25382932 0.58643326 0.06588542
#> 8  214.93717 480.735492 0.00000000  0.0000000 0.25209080 0.74790920 0.03869208
#> 9   75.71219 228.600798 0.00000000  0.0000000 0.00000000 1.00000000 0.01846578
#>   exp.prob.1 exp.prob.2 exp.prob.3   raw.rsd.0   raw.rsd.1    raw.rsd.2
#> 1 0.32577965  0.1291256 0.05844828  0.26150167 -0.16281669 -0.040236697
#> 2 0.33819914  0.1684550 0.09512054  0.23592106 -0.06990646 -0.092574281
#> 3 0.33159913  0.2105493 0.15155591  0.09454819  0.05869622 -0.054431128
#> 4 0.30837413  0.2440974 0.21904196  0.06964572  0.05010288 -0.031453771
#> 5 0.27111018  0.2678886 0.30008300 -0.05914419  0.03327825  0.077116077
#> 6 0.22342791  0.2783265 0.39305231 -0.01974298 -0.05945562  0.041919881
#> 7 0.17540735  0.2738887 0.48481854 -0.06588542 -0.01566993 -0.020059378
#> 8 0.13015779  0.2567947 0.57435543 -0.03869208 -0.13015779 -0.004703902
#> 9 0.08385578  0.2233398 0.67433864 -0.01846578 -0.08385578 -0.223339804
#>     raw.rsd.3
#> 1 -0.05844828
#> 2 -0.07344032
#> 3 -0.09881329
#> 4 -0.08829483
#> 5 -0.05125013
#> 6  0.03727871
#> 7  0.10161473
#> 8  0.17355377
#> 9  0.32566136
#> 
#> $V214
#>    total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1    351        317         34          0          0 253.794941  42.834200
#> 2    218        172         24         22          0 134.918384  29.668287
#> 3    370        265         41         46         18 191.133607  52.208898
#> 4   1303       1016        187         75         25 534.637002 179.995633
#> 5   2248        713        373        799        363 667.221101 282.639875
#> 6   2130         85        175        854       1016 424.998919 226.369418
#> 7   1269         15         24        253        977 162.265520 107.710347
#> 8    265          0          9         83        173  21.564078  17.504119
#> 9    150          0          0         13        137   7.204947   7.270299
#> 10    44          0          0          0         44   1.110013   1.439040
#>    exp.freq.2 exp.freq.3 obs.prop.0 obs.prop.1 obs.prop.2 obs.prop.3 exp.prob.0
#> 1   38.467628   15.90323 0.90313390 0.09686610 0.00000000 0.00000000 0.72306251
#> 2   34.714497   18.69883 0.78899083 0.11009174 0.10091743 0.00000000 0.61889167
#> 3   75.883798   50.77370 0.71621622 0.11081081 0.12432432 0.04864865 0.51657732
#> 4  322.449701  265.91766 0.77973906 0.14351497 0.05755948 0.01918649 0.41031236
#> 5  637.080262  661.05876 0.31717082 0.16592527 0.35542705 0.16147687 0.29680654
#> 6  641.571021  837.06064 0.03990610 0.08215962 0.40093897 0.47699531 0.19953001
#> 7  380.439550  618.58458 0.01182033 0.01891253 0.19936958 0.76989756 0.12786881
#> 8   75.604336  150.32747 0.00000000 0.03396226 0.31320755 0.65283019 0.08137388
#> 9   39.036491   96.48826 0.00000000 0.00000000 0.08666667 0.91333333 0.04803298
#> 10   9.926924   31.52402 0.00000000 0.00000000 0.00000000 1.00000000 0.02522758
#>    exp.prob.1 exp.prob.2 exp.prob.3   raw.rsd.0   raw.rsd.1   raw.rsd.2
#> 1  0.12203476  0.1095944 0.04530835  0.18007139 -0.02516866 -0.10959438
#> 2  0.13609306  0.1592408 0.08577446  0.17009915 -0.02600132 -0.05832338
#> 3  0.14110513  0.2050913 0.13722621  0.19963890 -0.03029432 -0.08076702
#> 4  0.13813940  0.2474672 0.20408109  0.36942671  0.00537557 -0.18990768
#> 5  0.12572948  0.2833987 0.29406529  0.02036428  0.04019578  0.07202835
#> 6  0.10627672  0.3012071 0.39298622 -0.15962391 -0.02411710  0.09973192
#> 7  0.08487813  0.2997948 0.48745830 -0.11604848 -0.06596560 -0.10042518
#> 8  0.06605328  0.2852994 0.56727346 -0.08137388 -0.03209102  0.02790817
#> 9  0.04846866  0.2602433 0.64325508 -0.04803298 -0.04846866 -0.17357661
#> 10 0.03270546  0.2256119 0.71645506 -0.02522758 -0.03270546 -0.22561190
#>      raw.rsd.3
#> 1  -0.04530835
#> 2  -0.08577446
#> 3  -0.08857756
#> 4  -0.18489460
#> 5  -0.13258842
#> 6   0.08400909
#> 7   0.28243926
#> 8   0.08555673
#> 9   0.27007825
#> 10  0.28354494
#> 
#> $V217
#>   total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3  exp.freq.0 exp.freq.1
#> 1   909        872         32          5          0  652.823460 148.604849
#> 2  2243       1552        485        192         14 1348.046871 424.894714
#> 3  3118       1497        671        600        350 1335.829332 616.252804
#> 4  3119        525        610        989        995  820.080996 541.849213
#> 5  1882         36        116        567       1163  250.970188 239.400313
#> 6   538          0          4         35        499   31.891121  43.586768
#> 7    70          0          2          1         67    1.388533   2.903588
#>   exp.freq.2 exp.freq.3 obs.prop.0  obs.prop.1 obs.prop.2  obs.prop.3
#> 1   75.99131   31.58038 0.95929593 0.035203520 0.00550055 0.000000000
#> 2  299.28184  170.77657 0.69193045 0.216228266 0.08559964 0.006241641
#> 3  635.31564  530.60222 0.48011546 0.215202053 0.19243105 0.112251443
#> 4  800.06054  957.00925 0.16832318 0.195575505 0.31708881 0.319012504
#> 5  510.32871  881.30079 0.01912859 0.061636557 0.30127524 0.617959617
#> 6  133.12581  329.39630 0.00000000 0.007434944 0.06505576 0.927509294
#> 7   13.62404   52.08384 0.00000000 0.028571429 0.01428571 0.957142857
#>   exp.prob.0 exp.prob.1 exp.prob.2 exp.prob.3   raw.rsd.0   raw.rsd.1
#> 1 0.71817762 0.16348168  0.0835988 0.03474190  0.24111831 -0.12827816
#> 2 0.60100173 0.18943144  0.1334293 0.07613757  0.09092872  0.02679683
#> 3 0.42842506 0.19764362  0.2037574 0.17017390  0.05169040  0.01755843
#> 4 0.26293075 0.17372530  0.2565119 0.30683208 -0.09460757  0.02185020
#> 5 0.13335292 0.12720527  0.2711630 0.46827885 -0.11422433 -0.06556871
#> 6 0.05927718 0.08101630  0.2474457 0.61226079 -0.05927718 -0.07358135
#> 7 0.01983618 0.04147983  0.1946291 0.74405490 -0.01983618 -0.01290840
#>     raw.rsd.2   raw.rsd.3
#> 1 -0.07809825 -0.03474190
#> 2 -0.04782962 -0.06989593
#> 3 -0.01132638 -0.05792246
#> 4  0.06057693  0.01218043
#> 5  0.03011227  0.14968077
#> 6 -0.18238998  0.31524851
#> 7 -0.18034337  0.21308796
#> 
#> $V220
#>    total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1    668        668          0          0          0  473.75960   73.39859
#> 2    814        706        108          0          0  529.13005   93.91865
#> 3    970        749        115        106          0  576.05411  114.63599
#> 4    954        625        150        110         69  506.09850  113.40446
#> 5    908        338        140        206        224  433.08278  106.76677
#> 6   1291        494        133        197        467  531.98578  146.95643
#> 7    967        130        144        193        500  335.09925  104.05904
#> 8    795          0         42        230        523  227.78915   79.23713
#> 9    630          0          0         68        562  145.34100   56.79509
#> 10   162          0          0          0        162   30.80343   13.24365
#>    exp.freq.2 exp.freq.3 obs.prop.0 obs.prop.1 obs.prop.2 obs.prop.3 exp.prob.0
#> 1    56.21308   64.62874  1.0000000 0.00000000  0.0000000 0.00000000  0.7092210
#> 2    82.40645  108.54485  0.8673219 0.13267813  0.0000000 0.00000000  0.6500369
#> 3   112.77131  166.53860  0.7721649 0.11855670  0.1092784 0.00000000  0.5938702
#> 4   125.61607  208.88097  0.6551363 0.15723270  0.1153040 0.07232704  0.5305016
#> 5   130.11316  238.03729  0.3722467 0.15418502  0.2268722 0.24669604  0.4769634
#> 6   200.67677  411.38102  0.3826491 0.10302091  0.1525949 0.36173509  0.4120726
#> 7   159.73723  368.10448  0.1344364 0.14891417  0.1995863 0.51706308  0.3465349
#> 8   136.25246  351.72126  0.0000000 0.05283019  0.2893082 0.65786164  0.2865272
#> 9   109.71183  318.15208  0.0000000 0.00000000  0.1079365 0.89206349  0.2307000
#> 10   28.14728   89.80564  0.0000000 0.00000000  0.0000000 1.00000000  0.1901446
#>    exp.prob.1 exp.prob.2 exp.prob.3   raw.rsd.0     raw.rsd.1    raw.rsd.2
#> 1  0.10987813 0.08415131 0.09674961  0.29077905 -0.1098781285 -0.084151314
#> 2  0.11537918 0.10123642 0.13334748  0.21728495  0.0172989565 -0.101236424
#> 3  0.11818143 0.11625908 0.17168928  0.17829473  0.0003752725 -0.006980729
#> 4  0.11887260 0.13167303 0.21895280  0.12463470  0.0383601029 -0.016369042
#> 5  0.11758455 0.14329643 0.26215561 -0.10471672  0.0366004766  0.083575815
#> 6  0.11383147 0.15544289 0.31865300 -0.02942354 -0.0108105553 -0.002847999
#> 7  0.10761018 0.16518844 0.38066647 -0.21209850  0.0413039899  0.034397905
#> 8  0.09966934 0.17138674 0.44241668 -0.28652724 -0.0468391524  0.117921436
#> 9  0.09015094 0.17414577 0.50500330 -0.23069999 -0.0901509432 -0.066209258
#> 10 0.08175093 0.17374861 0.55435583 -0.19014463 -0.0817509344 -0.173748613
#>      raw.rsd.3
#> 1  -0.09674961
#> 2  -0.13334748
#> 3  -0.17168928
#> 4  -0.14662576
#> 5  -0.01545957
#> 6   0.04308209
#> 7   0.13639661
#> 8   0.21544495
#> 9   0.38706020
#> 10  0.44564417
#> 
#> $V222
#>   total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3  exp.freq.0 exp.freq.1
#> 1   394        345         45          4          0  284.596081   93.20502
#> 2  1394        931        419         28         16  802.986135  439.68162
#> 3  3200       1563       1205        244        188 1333.490397 1119.03319
#> 4  3861        859       1429        588        985  951.026569 1229.78715
#> 5  3760        173        746        601       2240  395.801024  818.96704
#> 6  2455          0        108        281       2066   91.813756  296.06164
#> 7   934          0          0         60        874   11.677873   56.96378
#> 8   425          0          0          0        425    1.381017   10.86376
#>   exp.freq.2  exp.freq.3 obs.prop.0 obs.prop.1 obs.prop.2 obs.prop.3
#> 1   9.927696    6.271204 0.87563452 0.11421320 0.01015228 0.00000000
#> 2  75.661490   75.670754 0.66786227 0.30057389 0.02008608 0.01147776
#> 3 295.122521  452.353892 0.48843750 0.37656250 0.07625000 0.05875000
#> 4 499.773884 1180.412400 0.22248122 0.37011137 0.15229215 0.25511526
#> 5 532.552274 2012.679661 0.04601064 0.19840426 0.15984043 0.59574468
#> 6 300.028719 1767.095888 0.00000000 0.04399185 0.11446029 0.84154786
#> 7  87.325305  778.033045 0.00000000 0.00000000 0.06423983 0.93576017
#> 8  26.857641  385.897579 0.00000000 0.00000000 0.00000000 1.00000000
#>    exp.prob.0 exp.prob.1 exp.prob.2 exp.prob.3    raw.rsd.0   raw.rsd.1
#> 1 0.722325079 0.23656096 0.02519720 0.01591676  0.153309439 -0.12234776
#> 2 0.576030226 0.31541006 0.05427653 0.05428318  0.091832041 -0.01483617
#> 3 0.416715749 0.34969787 0.09222579 0.14136059  0.071721751  0.02686463
#> 4 0.246316128 0.31851519 0.12944157 0.30572712 -0.023834905  0.05159618
#> 5 0.105266230 0.21781038 0.14163624 0.53528714 -0.059255592 -0.01940613
#> 6 0.037398679 0.12059537 0.12221129 0.71979466 -0.037398679 -0.07660352
#> 7 0.012503076 0.06098905 0.09349604 0.83301183 -0.012503076 -0.06098905
#> 8 0.003249452 0.02556180 0.06319445 0.90799430 -0.003249452 -0.02556180
#>      raw.rsd.2   raw.rsd.3
#> 1 -0.015044915 -0.01591676
#> 2 -0.034190452 -0.04280542
#> 3 -0.015975788 -0.08261059
#> 4  0.022850587 -0.05061186
#> 5  0.018204183  0.06045754
#> 6 -0.007751006  0.12175320
#> 7 -0.029256215  0.10274835
#> 8 -0.063194450  0.09200570
#> 
#> $V223
#>   total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1   173        169          4          0          0  150.41445   15.45807
#> 2   415        385         25          5          0  334.82686   47.63071
#> 3  1676       1368        195        113          0 1198.01840  235.50067
#> 4  3762       2611        653        340        158 2167.05329  603.67504
#> 5  3255       1205        542        920        588 1285.76857  519.69396
#> 6  3713        612        470       1199       1432  867.84595  496.99478
#> 7  2516         96        210        743       1467  308.63792  245.60750
#> 8  1234         12         16        231        975   70.07068   77.83838
#> 9   435          0          2         69        364   10.29288   16.09384
#>    exp.freq.2  exp.freq.3  obs.prop.0  obs.prop.1 obs.prop.2 obs.prop.3
#> 1    5.819837    1.307636 0.976878613 0.023121387 0.00000000 0.00000000
#> 2   24.822381    7.720057 0.927710843 0.060240964 0.01204819 0.00000000
#> 3  169.594012   72.886910 0.816229117 0.116348449 0.06742243 0.00000000
#> 4  616.065011  375.206651 0.694045720 0.173577884 0.09037746 0.04199894
#> 5  769.524106  680.013368 0.370199693 0.166513057 0.28264209 0.18064516
#> 6 1042.679193 1305.480081 0.164826286 0.126582278 0.32291947 0.38567196
#> 7  716.017554 1245.737023 0.038155803 0.083465819 0.29531002 0.58306836
#> 8  316.767761  769.323182 0.009724473 0.012965964 0.18719611 0.79011345
#> 9   92.287390  316.325896 0.000000000 0.004597701 0.15862069 0.83678161
#>   exp.prob.0 exp.prob.1 exp.prob.2  exp.prob.3   raw.rsd.0    raw.rsd.1
#> 1 0.86944771 0.08935302 0.03364068 0.007558592  0.10743091 -0.066231637
#> 2 0.80681170 0.11477279 0.05981297 0.018602548  0.12089914 -0.054531821
#> 3 0.71480812 0.14051353 0.10118974 0.043488610  0.10142100 -0.024165080
#> 4 0.57603756 0.16046652 0.16375997 0.099735952  0.11800816  0.013111365
#> 5 0.39501338 0.15966020 0.23641294 0.208913477 -0.02481369  0.006852854
#> 6 0.23373174 0.13385262 0.28081853 0.351597113 -0.06890545 -0.007270342
#> 7 0.12267008 0.09761824 0.28458567 0.495126003 -0.08451428 -0.014152426
#> 8 0.05678337 0.06307810 0.25669997 0.623438559 -0.04705889 -0.050112140
#> 9 0.02366179 0.03699732 0.21215492 0.727185967 -0.02366179 -0.032399624
#>     raw.rsd.2    raw.rsd.3
#> 1 -0.03364068 -0.007558592
#> 2 -0.04776477 -0.018602548
#> 3 -0.03376731 -0.043488610
#> 4 -0.07338251 -0.057737015
#> 5  0.04622915 -0.028268316
#> 6  0.04210094  0.034074850
#> 7  0.01072434  0.087942360
#> 8 -0.06950386  0.166674893
#> 9 -0.05353423  0.109595642
#> 
#> $V226
#>   total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 exp.freq.0 exp.freq.1
#> 1   439        439          0          0          0  296.56451   73.90696
#> 2   771        591        153         27          0  460.26543  140.76126
#> 3  1250        727        229        231         63  630.83738  238.17158
#> 4  1585        690        368        360        167  642.20668  299.75893
#> 5  1939        600        369        581        389  613.21468  346.75748
#> 6  1486        295        226        528        437  339.05169  236.41948
#> 7  1203         67        143        358        635  188.92608  161.67764
#> 8   568          0          2        230        336   59.04750   61.75040
#> 9   246          0          0          1        245   15.61834   20.30524
#>   exp.freq.2 exp.freq.3 obs.prop.0  obs.prop.1  obs.prop.2 obs.prop.3
#> 1   50.06356   18.46497  1.0000000 0.000000000 0.000000000  0.0000000
#> 2  117.01139   52.96192  0.7665370 0.198443580 0.035019455  0.0000000
#> 3  244.41769  136.57335  0.5816000 0.183200000 0.184800000  0.0504000
#> 4  380.31162  262.72277  0.4353312 0.232176656 0.227129338  0.1053628
#> 5  532.97800  446.04984  0.3094379 0.190304281 0.299638989  0.2006189
#> 6  448.09570  462.43313  0.1985195 0.152086137 0.355316285  0.2940781
#> 7  376.07822  476.31805  0.0556941 0.118869493 0.297589360  0.5278470
#> 8  175.52869  271.67341  0.0000000 0.003521127 0.404929577  0.5915493
#> 9   71.79002  138.28640  0.0000000 0.000000000 0.004065041  0.9959350
#>   exp.prob.0 exp.prob.1 exp.prob.2 exp.prob.3    raw.rsd.0    raw.rsd.1
#> 1  0.6755456 0.16835299  0.1140400 0.04206144  0.324454423 -0.168352990
#> 2  0.5969720 0.18256973  0.1517657 0.06869250  0.169564942  0.015873847
#> 3  0.5046699 0.19053727  0.1955341 0.10925868  0.076930095 -0.007337267
#> 4  0.4051777 0.18912235  0.2399442 0.16575569  0.030153517  0.043054304
#> 5  0.3162531 0.17883315  0.2748726 0.23004117 -0.006815206  0.011471131
#> 6  0.2281640 0.15909790  0.3015449 0.31119322 -0.029644473 -0.007011765
#> 7  0.1570458 0.13439538  0.3126170 0.39594186 -0.101351690 -0.015525888
#> 8  0.1039569 0.10871549  0.3090294 0.47829825 -0.103956874 -0.105194360
#> 9  0.0634892 0.08254161  0.2918294 0.56213983 -0.063489200 -0.082541609
#>     raw.rsd.2   raw.rsd.3
#> 1 -0.11403999 -0.04206144
#> 2 -0.11674628 -0.06869250
#> 3 -0.01073415 -0.05885868
#> 4 -0.01281490 -0.06039292
#> 5  0.02476637 -0.02942230
#> 6  0.05377140 -0.01711516
#> 7 -0.01502761  0.13190519
#> 8  0.09590019  0.11325104
#> 9 -0.28776432  0.43379513
#> 


## Example 2
## Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select the first two dichotomous items and the last polytomous item
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[c(1:2, 55), ]

# Generate ability values from a standard normal distribution
set.seed(10)
score <- rnorm(1000, mean = 0, sd = 1)

# Simulate response data
data <- simdat(x = x, theta = score, D = 1)

# Compute item fit statistics
fit2 <- irtfit(
  x = x, score = score, data = data, group.method = "equal.freq",
  n.width = 11, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
)

# View the fit statistics
fit2$fit_stat
#>     id     X2     G2 df.X2 df.G2 crit.val.X2 crit.val.G2  p.X2  p.G2 outfit
#> 1 CMC1  8.541  8.684     8    11       15.51       19.68 0.382 0.651  1.007
#> 2 CMC2  8.812 10.095     7    10       14.07       18.31 0.266 0.432  0.843
#> 3 AFR3 42.396 43.817    39    44       54.57       60.48 0.327 0.479  0.950
#>   infit    N overSR.prop
#> 1 1.003 1000       0.000
#> 2 1.003 1000       0.000
#> 3 0.955 1000       0.036

# View the contingency tables used to compute fit statistics
fit2$contingency.fitstat
#> $CMC1
#>    total obs.freq.0 obs.freq.1 exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1     91         70         21   61.91419   29.08581  0.7692308  0.2307692
#> 2     91         56         35   58.97557   32.02443  0.6153846  0.3846154
#> 3     91         58         33   56.57872   34.42128  0.6373626  0.3626374
#> 4     91         55         36   54.61312   36.38688  0.6043956  0.3956044
#> 5     91         50         41   52.85925   38.14075  0.5494505  0.4505495
#> 6     90         53         37   49.96450   40.03550  0.5888889  0.4111111
#> 7     91         42         49   47.89655   43.10345  0.4615385  0.5384615
#> 8     91         44         47   45.00921   45.99079  0.4835165  0.5164835
#> 9     91         45         46   41.96040   49.03960  0.4945055  0.5054945
#> 10    91         36         55   38.18202   52.81798  0.3956044  0.6043956
#> 11    91         36         55   30.10521   60.89479  0.3956044  0.6043956
#>    exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1   0.6803758  0.3196242  0.088855008 -0.088855008
#> 2   0.6480831  0.3519169 -0.032698520  0.032698520
#> 3   0.6217442  0.3782558  0.015618466 -0.015618466
#> 4   0.6001442  0.3998558  0.004251382 -0.004251382
#> 5   0.5808709  0.4191291 -0.031420344  0.031420344
#> 6   0.5551611  0.4448389  0.033727786 -0.033727786
#> 7   0.5263357  0.4736643 -0.064797264  0.064797264
#> 8   0.4946067  0.5053933 -0.011090181  0.011090181
#> 9   0.4611033  0.5388967  0.033402221 -0.033402221
#> 10  0.4195826  0.5804174 -0.023978208  0.023978208
#> 11  0.3308265  0.6691735  0.064777861 -0.064777861
#> 
#> $CMC2
#>    total obs.freq.0 obs.freq.1 exp.freq.0 exp.freq.1 obs.prop.0 obs.prop.1
#> 1     91         56         35  59.847423   31.15258 0.61538462  0.3846154
#> 2     91         41         50  40.052623   50.94738 0.45054945  0.5494505
#> 3     91         28         63  26.361491   64.63851 0.30769231  0.6923077
#> 4     91         14         77  18.322089   72.67791 0.15384615  0.8461538
#> 5     91         17         74  13.239430   77.76057 0.18681319  0.8131868
#> 6     90         11         79   8.575110   81.42489 0.12222222  0.8777778
#> 7     91          4         87   5.498524   85.50148 0.04395604  0.9560440
#> 8     91          6         85   3.409008   87.59099 0.06593407  0.9340659
#> 9     91          1         90   2.103482   88.89652 0.01098901  0.9890110
#> 10   182          0        182   1.535797  180.46420 0.00000000  1.0000000
#>     exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1
#> 1  0.657663986  0.3423360 -0.042279370  0.042279370
#> 2  0.440138716  0.5598613  0.010410735 -0.010410735
#> 3  0.289686715  0.7103133  0.018005593 -0.018005593
#> 4  0.201341642  0.7986584 -0.047495489  0.047495489
#> 5  0.145488246  0.8545118  0.041324940 -0.041324940
#> 6  0.095278998  0.9047210  0.026943224 -0.026943224
#> 7  0.060423345  0.9395767 -0.016467301  0.016467301
#> 8  0.037461630  0.9625384  0.028472436 -0.028472436
#> 9  0.023115182  0.9768848 -0.012126171  0.012126171
#> 10 0.008438447  0.9915616 -0.008438447  0.008438447
#> 
#> $AFR3
#>    total obs.freq.0 obs.freq.1 obs.freq.2 obs.freq.3 obs.freq.4 exp.freq.0
#> 1     91         70          9          7          3          2  64.382932
#> 2     91         53         13          7          5         13  52.714982
#> 3     91         35         18         14          5         19  45.010677
#> 4     91         43         15         13          5         15  39.682474
#> 5     91         39         20          8         12         12  35.552292
#> 6     90         38         12         14          8         18  30.464603
#> 7     91         30         13         19         11         18  26.322302
#> 8     91         20         18         13          8         32  22.209861
#> 9     91         17         17         13         10         34  18.601260
#> 10    91         13         11         12          7         48  14.938697
#> 11    91          7          5          7         11         61   9.220059
#>    exp.freq.1 exp.freq.2 exp.freq.3 exp.freq.4 obs.prop.0 obs.prop.1 obs.prop.2
#> 1   10.550129   6.381702   3.612889   6.072348 0.76923077 0.09890110 0.07692308
#> 2   13.386387   9.152465   5.591455  10.154711 0.58241758 0.14285714 0.07692308
#> 3   14.469492  10.824328   7.029807  13.665695 0.38461538 0.19780220 0.15384615
#> 4   14.785266  11.830812   8.068613  16.632834 0.47252747 0.16483516 0.14285714
#> 5   14.754971  12.480281   8.881067  19.331389 0.42857143 0.21978022 0.08791209
#> 6   14.232327  12.883558   9.684825  22.734686 0.42222222 0.13333333 0.15555556
#> 7   13.688501  13.270488  10.580120  27.138589 0.32967033 0.14285714 0.20879121
#> 8   12.704352  13.174325  11.168933  31.742529 0.21978022 0.19780220 0.14285714
#> 9   11.544767  12.751245  11.491416  36.611312 0.18681319 0.18681319 0.14285714
#> 10  10.056405  11.893203  11.501006  42.610689 0.14285714 0.12087912 0.13186813
#> 11   7.029158   9.345416  10.316679  55.088688 0.07692308 0.05494505 0.07692308
#>    obs.prop.3 obs.prop.4 exp.prob.0 exp.prob.1 exp.prob.2 exp.prob.3 exp.prob.4
#> 1  0.03296703 0.02197802  0.7075048 0.11593548 0.07012859 0.03970208  0.0667291
#> 2  0.05494505 0.14285714  0.5792855 0.14710316 0.10057654 0.06144456  0.1115902
#> 3  0.05494505 0.20879121  0.4946228 0.15900541 0.11894866 0.07725063  0.1501725
#> 4  0.05494505 0.16483516  0.4360711 0.16247545 0.13000892 0.08866608  0.1827784
#> 5  0.13186813 0.13186813  0.3906845 0.16214254 0.13714595 0.09759414  0.2124328
#> 6  0.08888889 0.20000000  0.3384956 0.15813696 0.14315065 0.10760917  0.2526076
#> 7  0.12087912 0.19780220  0.2892561 0.15042309 0.14582953 0.11626506  0.2982263
#> 8  0.08791209 0.35164835  0.2440644 0.13960827 0.14477280 0.12273553  0.3488190
#> 9  0.10989011 0.37362637  0.2044095 0.12686557 0.14012357 0.12627930  0.4023221
#> 10 0.07692308 0.52747253  0.1641615 0.11050994 0.13069454 0.12638469  0.4682493
#> 11 0.12087912 0.67032967  0.1013193 0.07724349 0.10269687 0.11337010  0.6053702
#>       raw.rsd.0    raw.rsd.1    raw.rsd.2    raw.rsd.3    raw.rsd.4
#> 1   0.061726017 -0.017034382  0.006794485 -0.006735045 -0.044751075
#> 2   0.003132066 -0.004246016 -0.023653459 -0.006499501  0.031266910
#> 3  -0.110007440  0.038796789  0.034897491 -0.022305576  0.058618736
#> 4   0.036456328  0.002359713  0.012848220 -0.033721026 -0.017943235
#> 5   0.037886906  0.057637678 -0.049233859  0.034273988 -0.080564714
#> 6   0.083726632 -0.024803631  0.012404906 -0.018720281 -0.052607626
#> 7   0.040414266 -0.007565948  0.062961675  0.004614064 -0.100424057
#> 8  -0.024284190  0.058193931 -0.001915655 -0.034823441  0.002829356
#> 9  -0.017596263  0.059947617  0.002733575 -0.016389191 -0.028695738
#> 10 -0.021304362  0.010369178  0.001173592 -0.049461610  0.059223202
#> 11 -0.024396252 -0.022298440 -0.025773797  0.007509020  0.064959469
#> 

# Plot raw and standardized residuals for the first item (dichotomous)
plot(x = fit2, item.loc = 1, type = "both", ci.method = "wald",
     show.table = TRUE, ylim.sr.adjust = TRUE)

#>                   interval       point total obs.freq.0 obs.freq.1 obs.prop.0
#> 1    [-3.012164,-1.359216) -1.76327705    91         70         21  0.7692308
#> 2   [-1.359216,-0.8869352) -1.12022860    91         56         35  0.6153846
#> 3   [-0.8869352,-0.600789) -0.73065886    91         58         33  0.6373626
#> 4   [-0.600789,-0.3503181) -0.46176301    91         55         36  0.6043956
#> 5  [-0.3503181,-0.1182627) -0.24799010    91         50         41  0.5494505
#> 6   [-0.1182627,0.1417995)  0.00944368    90         53         37  0.5888889
#> 7    [0.1417995,0.4048853)  0.27078570    91         42         49  0.4615385
#> 8    [0.4048853,0.6625933)  0.53501489    91         44         47  0.4835165
#> 9    [0.6625933,0.9468582)  0.79571925    91         45         46  0.4945055
#> 10    [0.9468582,1.260075)  1.10229030    91         36         55  0.3956044
#> 11      [1.260075,3.54114]  1.73576474    91         36         55  0.3956044
#>    obs.prop.1 exp.prob.0 exp.prob.1    raw.rsd.0    raw.rsd.1       se.0
#> 1   0.2307692  0.6803758  0.3196242  0.088855008 -0.088855008 0.04888477
#> 2   0.3846154  0.6480831  0.3519169 -0.032698520  0.032698520 0.05006275
#> 3   0.3626374  0.6217442  0.3782558  0.015618466 -0.015618466 0.05083677
#> 4   0.3956044  0.6001442  0.3998558  0.004251382 -0.004251382 0.05135217
#> 5   0.4505495  0.5808709  0.4191291 -0.031420344  0.031420344 0.05172411
#> 6   0.4111111  0.5551611  0.4448389  0.033727786 -0.033727786 0.05238291
#> 7   0.5384615  0.5263357  0.4736643 -0.064797264  0.064797264 0.05234149
#> 8   0.5164835  0.4946067  0.5053933 -0.011090181  0.011090181 0.05241119
#> 9   0.5054945  0.4611033  0.5388967  0.033402221 -0.033402221 0.05225540
#> 10  0.6043956  0.4195826  0.5804174 -0.023978208  0.023978208 0.05173188
#> 11  0.6043956  0.3308265  0.6691735  0.064777861 -0.064777861 0.04932293
#>          se.1   std.rsd.0   std.rsd.1
#> 1  0.04888477  1.81764194 -1.81764194
#> 2  0.05006275 -0.65315069  0.65315069
#> 3  0.05083677  0.30722772 -0.30722772
#> 4  0.05135217  0.08278874 -0.08278874
#> 5  0.05172411 -0.60746032  0.60746032
#> 6  0.05238291  0.64387000 -0.64387000
#> 7  0.05234149 -1.23797144  1.23797144
#> 8  0.05241119 -0.21159947  0.21159947
#> 9  0.05225540  0.63921088 -0.63921088
#> 10 0.05173188 -0.46350931  0.46350931
#> 11 0.04932293  1.31334183 -1.31334183

# Plot raw and standardized residuals for the third item (polytomous)
plot(x = fit2, item.loc = 3, type = "both", ci.method = "wald",
     show.table = FALSE, ylim.sr.adjust = TRUE)

# }
```
