# Classical Test Theory (CTT) Item and Test Analysis

Computes classical test theory (CTT) statistics from a scored item
response data set (dichotomous or polytomous). For each item, computes
its difficulty, an item-total correlation reflecting how well the item
discriminates between high- and low-scoring examinees, and Cronbach's
alpha recomputed with that item removed, with optional flagging of items
whose difficulty or discrimination falls outside conventional quality
thresholds. At the test level, computes Cronbach's alpha as a
reliability estimate (in both raw and standardized forms), the
associated standard error of measurement (SEM), and the average item
difficulty and discrimination across the test. The distribution of total
(raw) scores across examinees is also tabulated. The result is returned
as an object of class `"ctt"`, with
[`print.ctt()`](https://hwangQ.github.io/irtQ/reference/print.ctt.md)
and
[`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md)
methods that display a condensed and a full report, respectively,
analogous to how
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) pairs
with
`print.est_irt()`/[`summary.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md).

## Usage

``` r
ctt(
  data,
  item.id = NULL,
  cats = NULL,
  correct = FALSE,
  missing = NA,
  flag = TRUE,
  crit.p = c(0.1, 0.95),
  crit.dis = 0.2
)
```

## Arguments

- data:

  A data frame or matrix of already-scored item responses, with
  examinees in rows and items in columns. Item scores must range from 0
  to `cats[j] - 1` for each item j (0/1 for a dichotomous item; 0, 1, 2,
  ... for a polytomous/partial-credit item).

- item.id:

  A character vector of item identifiers, in the same order as the
  columns of `data`. If `NULL` (default), item IDs are generated
  automatically as `paste0("V", 1:ncol(data))`, following the convention
  used elsewhere in irtQ (e.g.,
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)).
  Note that this generated ID does not fall back to `colnames(data)`;
  pass `item.id` explicitly to label items using the column names of
  `data` or any other identifier scheme.

- cats:

  A numeric vector giving the number of score categories for each item
  (e.g., 2 for a dichotomous item), following the `cats` convention used
  elsewhere in irtQ (see, e.g.,
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)).
  If `NULL` (default), the number of categories for each item is
  inferred from the observed maximum score in `data` (i.e.,
  `max(data[, j], na.rm = TRUE) + 1`); supply `cats` explicitly whenever
  the maximum possible score may not have been observed in the sample.

- correct:

  Logical. Both the raw (uncorrected) item-total correlation - where an
  item is correlated with the total score that includes its own
  contribution - and the corrected item-total correlation - excluding
  its own contribution - are always computed and reported as separate
  columns. This argument only controls which of the two is used for
  flagging: if `TRUE`, flagging is based on the corrected value; if
  `FALSE` (default), on the raw value.

- missing:

  A value indicating missing responses in `data`, analogous to the
  `missing` argument in
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) and
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
  Any cell equal to `missing` is recoded to `NA` before analysis.
  Default is `NA`. Examinees with any remaining missing item response
  are excluded listwise from all statistics computed by this function.

- flag:

  Logical. If `TRUE` (default), items are flagged when their difficulty
  or discrimination falls outside the thresholds given in `crit.p` and
  `crit.dis`.

- crit.p:

  A numeric vector of length two giving the lower and upper difficulty
  bounds used for flagging: difficulty below the first value is flagged
  as too difficult, and difficulty above the second value is flagged as
  too easy. Default is `c(0.10, 0.95)`.

- crit.dis:

  A single numeric value giving the minimum acceptable discrimination
  (item-total correlation); items strictly below this value are flagged
  as poorly discriminating. Default is `0.20`.

## Value

An object of class `"ctt"`, a list with elements:

- item:

  A data frame with one row per item, containing the item label, number
  of score categories, difficulty, the raw (uncorrected) item-total
  correlation (`discrimination_raw`), the corrected (item-excluded)
  item-total correlation (`discrimination_corrected`),
  alpha-with-item-removed, and, if `flag = TRUE`, a `flag` column.

- crit:

  A list echoing the `crit.p` and `crit.dis` thresholds used for
  flagging.

- alpha:

  A one-row test-level summary data frame containing `n_examinee`,
  `n_item`, `alpha`, `alpha_std`, `sem`, `mean_difficulty`,
  `mean_discrimination_raw`, and `mean_discrimination_corrected`. Both
  the raw (`alpha`) and standardized (`alpha_std`) forms of Cronbach's
  alpha are always included; standardized alpha standardizes every item
  to unit variance before combining them, which differs meaningfully
  from raw alpha mainly when items vary widely in scale or format.

- freq:

  The total-score frequency distribution table returned by
  [`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md)
  (`score`, `freq`, `pct`, `cum_pct`).

- call:

  The matched call, as used by the
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://hwangQ.github.io/irtQ/reference/summary.md)
  methods for `"ctt"` objects.

## Details

Difficulty for item j is the mean observed item score divided by the
item's maximum possible score, \$\$p_j = \frac{\bar{X}\_j}{m_j},\$\$
where \\m_j\\ is `cats[j] - 1`. For a dichotomous item (`cats[j] = 2`)
this reduces to the familiar proportion-correct difficulty index; for a
polytomous item it expresses the average score as a proportion of the
maximum attainable score, keeping difficulty on a common 0-1 scale
regardless of the number of score categories.

Discrimination for item j is the Pearson correlation between the item
score and the total score, always computed and reported in two forms:
the raw (uncorrected) item-total correlation, which correlates the item
with the total score that includes the item's own contribution, and the
corrected (item-excluded) item-total correlation, which excludes it. For
a dichotomous item, the raw item-total correlation is mathematically
equivalent to the point-biserial correlation. `correct` selects which of
the two feeds the discrimination flagging criterion (`crit.dis`); both
are always reported as separate columns regardless of `correct` (this
mirrors how, e.g., `psych::alpha()` reports the raw item-total
correlation as `raw.r` and the corrected version as `r.drop`).

Alpha-with-item-removed for item j is Cronbach's alpha (see below)
recomputed using only the remaining items, so that a low value flags an
item whose removal would increase the overall reliability of the test.

At the test level, two forms of Cronbach's alpha are always computed and
reported. Raw alpha uses the standard variance-based formula \$\$\alpha
= \frac{k}{k - 1}\left(1 - \frac{\sum\_{j}
\sigma_j^2}{\sigma_X^2}\right),\$\$ where \\k\\ is the number of items,
\\\sigma_j^2\\ is the variance of item j, and \\\sigma_X^2\\ is the
variance of the total score X. Raw alpha reflects the reliability of the
actual (unweighted) total score obtained by simply summing the item
scores - the score most tests actually use for reporting and decisions.
Standardized alpha instead first standardizes every item to unit
variance before combining them, \$\$\alpha\_{std} = \frac{k \bar{r}}{1 +
(k - 1)\bar{r}},\$\$ where \\\bar{r}\\ is the average pairwise
correlation among all items. Standardized alpha differs meaningfully
from raw alpha mainly when items vary substantially in scale or format
(e.g., a mix of dichotomous and polytomous items with very different
score ranges); when all items share the same scale and format, the two
are typically close, and raw alpha remains the more directly
interpretable of the two since it matches the reliability of the score
actually used in practice. See Cronbach (1951) for the original
derivation of coefficient alpha, and Osburn (2000) for a discussion
contrasting the raw (covariance-based) and standardized
(correlation-based) forms.

The standard error of measurement (SEM) is \\SEM = SD(X)\sqrt{1 -
\alpha}\\ (using raw alpha), the standard CTT relationship between test
reliability and measurement precision.

The total-score frequency distribution is computed internally by
[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md),
using the row sums of `data` after applying the same `missing` recoding
and listwise deletion of incomplete rows used throughout the rest of
this function, so that the frequency distribution reflects exactly the
same set of examinees and the same total-score definition used elsewhere
in the analysis.
[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md)
is also exported separately, for callers who only need a frequency
distribution for an arbitrary vector of integer total scores.

## References

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika*, *16*(3), 297-334.
https://doi.org/10.1007/BF02310555

Osburn, H. G. (2000). Coefficient alpha and related internal consistency
reliability coefficients. *Psychological Methods*, *5*(3), 343-355.
https://doi.org/10.1037/1082-989X.5.3.343

## See also

[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md),
[`print.ctt()`](https://hwangQ.github.io/irtQ/reference/print.ctt.md),
[`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# A small dichotomous example
set.seed(1)
dat <- data.frame(matrix(rbinom(300 * 8, 1, 0.6), nrow = 300))
out <- ctt(data = dat)
out
#> 
#> Call:
#> ctt(data = dat)
#> 
#> Classical Test Theory (CTT) Analysis
#>  Number of items: 8
#>  Number of examinees: 300
#>  Cronbach's alpha: 0.081
#>  SEM: 1.372
#>  Mean difficulty: 0.613
#>  Mean discrimination (raw / corrected): 0.367 / 0.028
#>  Flagged items: 0 of 8
#> 
#> Use summary() for the full item-level and frequency-distribution report.
summary(out)
#> 
#> Call:
#> ctt(data = dat)
#> 
#> Item-Level Statistics
#>   item  cats  difficulty  discrimination_raw  discrimination_corrected
#>     V1     2       0.630               0.344                     0.007
#>     V2     2       0.587               0.346                     0.001
#>     V3     2       0.603               0.403                     0.066
#>     V4     2       0.620               0.365                     0.027
#>     V5     2       0.620               0.389                     0.054
#>     V6     2       0.590               0.385                     0.044
#>     V7     2       0.623               0.337                    -0.002
#>     V8     2       0.630               0.364                     0.028
#>   alpha_removed  flag
#>           0.088      
#>           0.093      
#>           0.041      
#>           0.072      
#>           0.051      
#>           0.059      
#>           0.095      
#>           0.072      
#> 
#> Flagging thresholds: difficulty in [0.1, 0.95], discrimination >= 0.2
#> 0 of 8 item(s) flagged.
#> 
#> Test-Level Reliability Summary
#>   n_examinee  n_item  alpha  alpha_std    sem  mean_difficulty
#>          300       8  0.081      0.081  1.372            0.613
#>   mean_discrimination_raw  mean_discrimination_corrected
#>                     0.367                          0.028
#> 
#> Total-Score Frequency Distribution
#>   score  freq    pct  cum_pct
#>       1     4   1.33     1.33
#>       2    11   3.67     5.00
#>       3    34  11.33    16.33
#>       4    61  20.33    36.67
#>       5    87  29.00    65.67
#>       6    64  21.33    87.00
#>       7    32  10.67    97.67
#>       8     7   2.33   100.00
#> 

# A more realistic, mixed-format example: simulate response data for a
# 55-item test (50 dichotomous 3PLM items + 5 polytomous GRM items) from
# item parameters imported from a flexMIRT output file bundled with irtQ,
# then run ctt() on the simulated data
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
set.seed(2)
theta <- rnorm(1000)
dat_mixed <- simdat(x = x, theta = theta, D = 1)
out_mixed <- ctt(data = dat_mixed, item.id = x$id, cats = x$cats)
summary(out_mixed)
#> 
#> Call:
#> ctt(data = dat_mixed, item.id = x$id, cats = x$cats)
#> 
#> Item-Level Statistics
#>    item  cats  difficulty  discrimination_raw  discrimination_corrected
#>    CMC1     2       0.461               0.219                     0.183
#>    CMC2     2       0.828               0.443                     0.419
#>    CMC3     2       0.501               0.344                     0.310
#>    CMC4     2       0.674               0.423                     0.393
#>    CMC5     2       0.599               0.338                     0.304
#>    CMC6     2       0.407               0.499                     0.470
#>    CMC7     2       0.382               0.350                     0.317
#>    CMC8     2       0.438               0.312                     0.277
#>    CMC9     2       0.519               0.326                     0.291
#>   CMC10     2       0.569               0.487                     0.458
#>   CMC11     2       0.671               0.381                     0.350
#>   CMC12     2       0.368               0.289                     0.255
#>   CMC13     2       0.328               0.327                     0.294
#>   CMC14     2       0.586               0.399                     0.367
#>   CMC15     2       0.636               0.462                     0.432
#>   CMC16     2       0.563               0.598                     0.573
#>   CMC17     2       0.623               0.421                     0.390
#>   CMC18     2       0.393               0.260                     0.225
#>   CMC19     2       0.824               0.506                     0.484
#>   CMC20     2       0.886               0.336                     0.314
#>   CMC21     2       0.824               0.418                     0.393
#>   CMC22     2       0.705               0.269                     0.236
#>   CMC23     2       0.606               0.408                     0.376
#>   CMC24     2       0.392               0.228                     0.193
#>   CMC25     2       0.826               0.203                     0.175
#>   CMC26     2       0.885               0.308                     0.285
#>   CMC27     2       0.571               0.437                     0.406
#>   CMC28     2       0.642               0.503                     0.475
#>   CMC29     2       0.847               0.338                     0.313
#>   CMC30     2       0.533               0.284                     0.248
#>   CMC31     2       0.413               0.268                     0.233
#>   CMC32     2       0.810               0.435                     0.410
#>   CMC33     2       0.819               0.362                     0.336
#>   CMC34     2       0.548               0.471                     0.441
#>   CMC35     2       0.627               0.452                     0.422
#>   CMC36     2       0.427               0.306                     0.271
#>   CMC37     2       0.670               0.522                     0.495
#>   CMC38     2       0.689               0.286                     0.253
#>    CFR1     5       0.752               0.667                     0.603
#>    CFR2     5       0.478               0.651                     0.571
#>    AMC1     2       0.489               0.407                     0.374
#>    AMC2     2       0.889               0.416                     0.396
#>    AMC3     2       0.473               0.440                     0.408
#>    AMC4     2       0.634               0.396                     0.364
#>    AMC5     2       0.301               0.201                     0.167
#>    AMC6     2       0.255               0.256                     0.224
#>    AMC7     2       0.597               0.439                     0.408
#>    AMC8     2       0.582               0.418                     0.386
#>    AMC9     2       0.494               0.454                     0.423
#>   AMC10     2       0.275               0.351                     0.321
#>   AMC11     2       0.829               0.462                     0.439
#>   AMC12     2       0.738               0.310                     0.279
#>    AFR1     5       0.394               0.552                     0.459
#>    AFR2     5       0.724               0.588                     0.509
#>    AFR3     5       0.448               0.500                     0.393
#>   alpha_removed  flag
#>           0.885      
#>           0.883      
#>           0.883      
#>           0.883      
#>           0.884      
#>           0.882      
#>           0.883      
#>           0.884      
#>           0.884      
#>           0.882      
#>           0.883      
#>           0.884      
#>           0.884      
#>           0.883      
#>           0.882      
#>           0.881      
#>           0.883      
#>           0.884      
#>           0.882      
#>           0.884      
#>           0.883      
#>           0.884      
#>           0.883      
#>           0.885      
#>           0.885      
#>           0.884      
#>           0.883      
#>           0.882      
#>           0.884      
#>           0.884      
#>           0.884      
#>           0.883      
#>           0.883      
#>           0.882      
#>           0.882      
#>           0.884      
#>           0.882      
#>           0.884      
#>           0.878      
#>           0.881      
#>           0.883      
#>           0.883      
#>           0.882      
#>           0.883      
#>           0.885      
#>           0.884      
#>           0.883      
#>           0.883      
#>           0.882      
#>           0.883      
#>           0.883      
#>           0.884      
#>           0.884      
#>           0.882      
#>           0.888      
#> 
#> Flagging thresholds: difficulty in [0.1, 0.95], discrimination >= 0.2
#> 0 of 55 item(s) flagged.
#> 
#> Test-Level Reliability Summary
#>   n_examinee  n_item  alpha  alpha_std   sem  mean_difficulty
#>         1000      55  0.885      0.903  4.45             0.59
#>   mean_discrimination_raw  mean_discrimination_corrected
#>                     0.395                           0.36
#> 
#> Total-Score Frequency Distribution
#>   score  freq  pct  cum_pct
#>       8     1  0.1      0.1
#>       9     0  0.0      0.1
#>      10     3  0.3      0.4
#>      11     2  0.2      0.6
#>      12     2  0.2      0.8
#>      13     5  0.5      1.3
#>      14     3  0.3      1.6
#>      15    10  1.0      2.6
#>      16     8  0.8      3.4
#>      17     7  0.7      4.1
#>      18     8  0.8      4.9
#>      19    13  1.3      6.2
#>      20    10  1.0      7.2
#>      21    14  1.4      8.6
#>      22    15  1.5     10.1
#>      23    13  1.3     11.4
#>      24    14  1.4     12.8
#>      25    15  1.5     14.3
#>      26    11  1.1     15.4
#>      27    20  2.0     17.4
#>      28    14  1.4     18.8
#>      29    22  2.2     21.0
#>      30    20  2.0     23.0
#>      31    30  3.0     26.0
#>      32    29  2.9     28.9
#>      33    23  2.3     31.2
#>      34    15  1.5     32.7
#>      35    23  2.3     35.0
#>      36    32  3.2     38.2
#>      37    33  3.3     41.5
#>      38    26  2.6     44.1
#>      39    17  1.7     45.8
#>      40    29  2.9     48.7
#>      41    24  2.4     51.1
#>      42    21  2.1     53.2
#>      43    21  2.1     55.3
#>      44    21  2.1     57.4
#>      45    30  3.0     60.4
#>      46    31  3.1     63.5
#>      47    19  1.9     65.4
#>      48    21  2.1     67.5
#>      49    27  2.7     70.2
#>      50    30  3.0     73.2
#>      51    16  1.6     74.8
#>      52    26  2.6     77.4
#>      53    32  3.2     80.6
#>      54    25  2.5     83.1
#>      55    22  2.2     85.3
#>      56    25  2.5     87.8
#>      57    15  1.5     89.3
#>      58    19  1.9     91.2
#>      59    19  1.9     93.1
#>      60    12  1.2     94.3
#>      61    11  1.1     95.4
#>      62    10  1.0     96.4
#>      63     6  0.6     97.0
#>      64     8  0.8     97.8
#>      65    11  1.1     98.9
#>      66     8  0.8     99.7
#>      67     0  0.0     99.7
#>      68     1  0.1     99.8
#>      69     2  0.2    100.0
#> 
```
