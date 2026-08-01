# Frequency Distribution Table for Total Scores

Computes a frequency distribution table for a vector of total (raw)
scores - frequency, percentage, and cumulative percentage for each score
value - commonly reported alongside classical test theory (CTT) item
analysis results. This is the same function
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) calls
internally to build the total-score frequency distribution included in
its output, but it is also exported and fully usable on its own for any
vector of integer total scores.

## Usage

``` r
freq_score(score, missing = NA)
```

## Arguments

- score:

  A numeric vector of total (raw) scores, one value per examinee.

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
two decimal places only in the final output, so that rounding error does
not accumulate across rows and the final cumulative percentage totals
almost exactly 100 (subject only to ordinary rounding).

This function assumes scores already lie on an integer (or otherwise
evenly spaced discrete) scale, as is standard for a raw total score; it
does not bin or group continuous values, and raises an error if any
non-integer score is supplied rather than silently dropping it.

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
