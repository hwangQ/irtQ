# Frequency Distribution Table for Total Scores

Computes a frequency distribution table for a vector of total (raw)
scores: the frequency, percentage, and cumulative percentage of each
score value, as commonly reported with classical test theory (CTT) item
analysis results.
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) calls this
function to build the total-score frequency distribution in its output,
and it can also be used on its own for any vector of integer total
scores.

## Usage

``` r
freq_score(score, missing = NA)
```

## Arguments

- score:

  A numeric vector of integer total (raw) scores, one value per
  examinee. A factor is converted through its labels (for example, the
  level `"12"` becomes the score 12), not through its level codes.

- missing:

  A value indicating missing scores in `score`, analogous to the
  `missing` argument in
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) and
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
  Any element equal to `missing` is recoded to `NA` before tabulation.
  Default is `NA`. Examinees with a (remaining) missing score are
  excluded from the table, with a warning reporting how many were
  dropped.

## Value

A data frame with one row per score value from `min(score)` to
`max(score)`, containing:

- score:

  the score value.

- freq:

  the number of examinees with that score.

- pct:

  the percentage of examinees with that score.

- cum_pct:

  the cumulative percentage up to and including that score.

## Details

The table includes every integer score value spanning the observed
minimum to maximum score (not just the values that were actually
observed), so that a score with zero examinees still appears in the
table with a frequency of 0, matching how a raw-score frequency table is
conventionally reported. Percentages are computed relative to the number
of non-missing scores and rounded to two decimal places; cumulative
percentages are the running sum of the unrounded percentages, rounded to
two decimal places only in the final output, so rounding error does not
accumulate and the last cumulative percentage is 100.

This function requires integer scores, as is standard for a raw total
score. It does not bin or group values, and it stops with an error when
a score is not a whole number. A score within 1e-8 of a whole number,
such as the result of a floating-point sum, is counted as that whole
number.

The output table always spans `min(score)` to `max(score)`, so its size
scales with the observed score *range*, not the sample size; a single
unusually large or small outlier score will produce a correspondingly
large table.

## See also

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
set.seed(1)
score <- rbinom(500, size = 20, prob = 0.6)
freq_score(score)
#>    score freq  pct cum_pct
#> 1      6    2  0.4     0.4
#> 2      7   10  2.0     2.4
#> 3      8   18  3.6     6.0
#> 4      9   36  7.2    13.2
#> 5     10   53 10.6    23.8
#> 6     11   73 14.6    38.4
#> 7     12   94 18.8    57.2
#> 8     13   93 18.6    75.8
#> 9     14   67 13.4    89.2
#> 10    15   41  8.2    97.4
#> 11    16    8  1.6    99.0
#> 12    17    3  0.6    99.6
#> 13    18    2  0.4   100.0
```
