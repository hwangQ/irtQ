# Classical Test Theory (CTT) Analysis

## Overview

The rest of **irtQ** centers on item response theory (IRT): item
parameters are estimated on a common ability scale, and item/test
statistics generally depend on those estimated parameters. Many testing
programs, however, still report - or additionally report - classical
test theory (CTT) statistics: simple, total-score-based indices of item
difficulty, item discrimination, and test reliability that require no
IRT model fitting at all. irtQ provides a small set of CTT functions for
exactly this purpose, working directly from raw selected-response data
or an already-scored item response matrix. For a fuller treatment of
classical test theory, see Crocker & Algina (1986).

| Function | Purpose |
|----|----|
| [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md) | Score raw selected-response data (option choices) into a 0/1 matrix (covered in the *Utility Functions* article; briefly revisited below) |
| [`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md) | Frequency distribution table (frequency, percentage, cumulative percentage) for a vector of total scores |
| [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md) | Option/category response distribution and point-biserial correlations, including distractor flagging |
| [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) | Full item- and test-level CTT analysis: difficulty, discrimination, alpha-with-item-removed, reliability, SEM, and the total-score frequency distribution |

``` r

library(irtQ)
```

### From Raw Responses to a Scored Matrix

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) expects an
already-scored item response matrix (0/1 for a dichotomous item, or 0,
1, 2, … for a polytomous item), the same convention used throughout irtQ
(e.g., by
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)), and
[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md)
expects a vector of total scores, such as the
[`rowSums()`](https://rdrr.io/r/base/colSums.html) of a scored matrix.
Raw testing data, however, usually records the option each examinee
actually chose. The *Utility Functions* article covers
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
in full - including its support for numeric, Latin-letter, and
general-label (e.g., Korean-syllable) option-coding schemes - but a
short example is worth repeating here, since it is the natural first
step before calling
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) on raw data
(or computing the total scores that
[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md)
tabulates):

``` r

# Raw responses for 5 examinees on 3 letter-coded items; item 2 has one
# double-marked response ("B,D") and item 3 has one blank (NA) response
raw3 <- data.frame(
  V1 = c("A", "B", NA, "a", "C"),
  V2 = c("D", "D", "D", "B,D", "D"),
  V3 = c("B", "A", "C", "C", NA)
)
key3 <- c("A", "D", "C")

scored3 <- score_resp(data = raw3, key = key3)
scored3$scored          # 0/1 matrix - the input format ctt() expects
#>   V1 V2 V3
#> 1  1  1  0
#> 2  0  1  0
#> 3  0  1  1
#> 4  1  0  1
#> 5  0  1  0
scored3$resp_summary     # per-item blank/double-mark/invalid counts
#>   item key n n_correct n_wrong n_blank n_double n_invalid pct_blank pct_double
#> 1   V1   A 5         2       3       1        0         0        20          0
#> 2   V2   D 5         4       1       0        1         0         0         20
#> 3   V3   C 5         2       3       1        0         0        20          0
```

[`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md),
described next, is the one CTT function that instead accepts raw
selected-response data directly (together with an answer key): it calls
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
internally so that its distractor analysis and the rest of the CTT
pipeline always agree on how a given response is scored.

The two example data sets used throughout the remainder of this article
are defined once here and reused below:

``` r

# Data set 1: a realistic, already-scored mixed-format response matrix -
# 1,000 examinees on a 55-item test (50 dichotomous 3PLM items + 5 polytomous
# GRM items), simulated from item parameters imported from a flexMIRT output
# file bundled with irtQ. This mirrors the example used in `?ctt`
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
set.seed(2026)
theta <- rnorm(1000)
dat_mixed <- simdat(x = x, theta = theta, D = 1)

dim(dat_mixed)   # 1,000 examinees x 55 items
#> [1] 1000   55
table(x$cats)    # 50 items with cats = 2 (dichotomous), 5 items with cats = 5 (GRM)
#> 
#>  2  5 
#> 50  5

# Data set 2: a small dichotomous-only response matrix (300 examinees on 15
# simulated 3PLM items), used for a simpler, easier-to-read walk-through of
# ctt()'s item- and test-level output
set.seed(1)
x_bin <- shape_df(
  par.drm = list(a = rep(1.5, 15), b = seq(-1.5, 1.5, length.out = 15),
                 g = rep(0.2, 15)),
  cats = 2, model = "3PLM"
)
dat_bin <- simdat(x = x_bin, theta = rnorm(300), D = 1)
dim(dat_bin)     # 300 examinees x 15 items
#> [1] 300  15
```

------------------------------------------------------------------------

## `freq_score()`: Total-Score Frequency Distribution

[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md)
tabulates a vector of total (raw) scores into a frequency distribution
table: frequency, percentage, and cumulative percentage for every
integer score value spanning the observed minimum to maximum (including
any value with zero observations in between, so gaps in the distribution
are visible rather than silently skipped). This is the same table
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) builds
internally for its own `$freq` element, but it is also fully usable on
its own for any vector of integer total scores - for instance, the row
sums of `dat_mixed` computed above:

``` r

total_mixed <- rowSums(dat_mixed)
freq_out <- freq_score(score = total_mixed)

head(freq_out)
#>   score freq pct cum_pct
#> 1     9    2 0.2     0.2
#> 2    10    1 0.1     0.3
#> 3    11    5 0.5     0.8
#> 4    12    3 0.3     1.1
#> 5    13    3 0.3     1.4
#> 6    14    4 0.4     1.8
tail(freq_out)
#>    score freq pct cum_pct
#> 56    64    9 0.9    98.5
#> 57    65    9 0.9    99.4
#> 58    66    0 0.0    99.4
#> 59    67    3 0.3    99.7
#> 60    68    2 0.2    99.9
#> 61    69    1 0.1   100.0
sum(freq_out$freq) == nrow(dat_mixed)   # every examinee accounted for
#> [1] TRUE
```

------------------------------------------------------------------------

## `ctt_distr()`: Option/Category Response Distribution and Distractor Analysis

[`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md)
reports, for every item, how examinees distributed themselves across
response options (or, for already-scored polytomous data, score
categories), together with each option’s/category’s point-biserial
correlation with the total score. It operates in one of two modes,
chosen automatically by whether a `key` is supplied.

### Selected-response mode: raw option data plus an answer key

When `key` is supplied, `data` is treated as raw, unscored option
responses - exactly the format
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
expects - and
[`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md)
additionally reports the blank/double-marking rate for each item:

``` r

raw <- data.frame(
  V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
  V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
)
key <- c(1, 4)
out_distr <- ctt_distr(data = raw, key = key)

out_distr$distr   # one row per item x option: freq, pct, pb_raw, pb_corrected, flag
#>   item option is_key freq  pct pb_raw pb_corrected flag
#> 1   V1      1   TRUE    4 50.0  0.626       -0.258     
#> 2   V1      2  FALSE    2 25.0 -0.602       -0.149     
#> 3   V1      3  FALSE    1 12.5 -0.079        0.293     
#> 4   V1      4  FALSE    1 12.5 -0.079        0.293     
#> 5   V2      1  FALSE    1 12.5 -0.709       -0.378     
#> 6   V2      2  FALSE    1 12.5 -0.079        0.378     
#> 7   V2      3  FALSE    0  0.0     NA           NA     
#> 8   V2      4   TRUE    5 62.5  0.592       -0.258
out_distr$omit    # one row per item: pct_blank, pct_double
#>   item pct_blank pct_double
#> 1   V1         0        0.0
#> 2   V2         0       12.5
```

With only two items, the corrected total is the score on the other item
alone, so `pb_corrected` is unstable here; use a full-length test in
practice.

Each item’s option-coding scheme (numeric, Latin-letter, or general
label) is inferred independently from its own `key` value, using exactly
the same rule as
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md) -
so a letter-coded or Korean-syllable-coded test form is analyzed
identically to a numerically coded one, once options are mapped.

**A scope note on mixed-format data:** selected-response mode is
restricted to **dichotomous, single-key items only**, since it scores
`data` via
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
internally, and a selected-response item has exactly one correct option
and therefore cannot yield a partial-credit score by construction. A
mixed dichotomous/polytomous test form cannot be analyzed in this mode.
The scored-category mode described further below has no such restriction
and fully supports a mixed-format test form scored in a single call -
see the *Scored-category mode* section.

### Flagging an attractive distractor

A distractor (a non-key option) whose point-biserial correlation with
the total score is positive - meaning it is chosen disproportionately by
*higher*-scoring examinees rather than lower-scoring ones - is a common
symptom of a miskeyed or ambiguously worded item.
[`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md)
flags any distractor exceeding `crit.distractor` (default `0`):

``` r

# A single item where the wrong option "2" is chosen only by the 5
# highest-ranked examinees, and the keyed option "1" only by the 5 lowest-
# ranked examinees (total scores supplied externally here purely to make the
# pattern explicit; in practice `total` is usually the test's own total score)
raw_flag <- data.frame(
  target = c("2", "2", "2", "2", "2", "1", "1", "1", "1", "1")
)
total_rank <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
out_flag <- ctt_distr(data = raw_flag, key = c(1), total = total_rank)

out_flag$distr
#>   item option is_key freq pct pb_raw pb_corrected
#> 1   V1      1   TRUE    5  50  -0.87       -0.905
#> 2   V1      2  FALSE    5  50   0.87        0.905
#>                                         flag
#> 1                                           
#> 2 distractor more attractive to high scorers
```

Option `"2"`’s `pb_raw` is strongly positive, and it is flagged as
“distractor more attractive to high scorers” - exactly the pattern to
look for when reviewing a flagged item for a possible keying error.

### Scored-category mode: already-scored item response data, mixed formats included

When `key` is left at its default (`NULL`), `data` is instead treated as
an already-scored response matrix, and
[`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md)
reports each score category’s response proportion and category-total
correlation, with no notion of a “correct” option to flag distractors
against (`is_key` is `NA` and `omit` is `NULL` in this mode).

`cats` is always a *per-item* vector in this mode - there is no
assumption anywhere that every item shares the same number of score
categories - so a single call can score a genuinely **mixed
dichotomous/polytomous test form** together, exactly as
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) does further
below. The call below scores all 55 items of `dat_mixed` (the 50
dichotomous items and the 5 polytomous GRM items) at once, not just a
same-format subset:

``` r

out_cat <- ctt_distr(data = dat_mixed, item.id = x$id, cats = x$cats)

# one dichotomous item (2 category rows) and one polytomous GRM item (5
# category rows), both scored within this same mixed-format call
dichot_id <- x$id[x$cats == 2][1]
poly_id   <- x$id[x$cats > 2][1]
out_cat$distr[out_cat$distr$item %in% c(dichot_id, poly_id), ]
#>    item option is_key freq  pct pb_raw pb_corrected flag
#> 1  CMC1      0     NA  548 54.8 -0.219       -0.182     
#> 2  CMC1      1     NA  452 45.2  0.219        0.182     
#> 77 CFR1      0     NA   85  8.5 -0.400       -0.352     
#> 78 CFR1      1     NA   96  9.6 -0.301       -0.268     
#> 79 CFR1      2     NA  119 11.9 -0.221       -0.206     
#> 80 CFR1      3     NA  134 13.4 -0.073       -0.078     
#> 81 CFR1      4     NA  566 56.6  0.599        0.546
```

Each item’s own `cats` value determines how many category rows it
contributes to `out_cat$distr` - two rows for a dichotomous item, five
for a 5-category GRM item - regardless of what any other item in the
same call looks like.

For `dichot_id` above, notice that `pb_raw` for category 0 and category
1 are equal in magnitude but opposite in sign - this is expected, not a
computation error. With only two categories, choosing category 0 is the
exact complement of choosing category 1 (every scored response is one or
the other), so the “in category 0” indicator equals one minus the “in
category 1” indicator, and correlating a variable’s complement with any
reference variable simply flips the sign of the correlation while
leaving its magnitude unchanged. This exact symmetry is specific to
two-category (dichotomous) items; it does not hold for `poly_id` or any
other item with three or more categories, where no such simple
complementary relationship exists between any two categories.

------------------------------------------------------------------------

## `ctt()`: Full Item- and Test-Level CTT Analysis

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) is the main
entry point for a complete CTT analysis: for each item it computes
difficulty, an item-total correlation (reported in both raw and
corrected forms), and Cronbach’s alpha recomputed with that item
removed, with optional flagging of items whose difficulty or
discrimination falls outside conventional thresholds; at the test level
it computes Cronbach’s alpha (raw and standardized), the standard error
of measurement (SEM), and average item difficulty/discrimination; and it
tabulates the total-score frequency distribution via
[`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md).

### A simple dichotomous example

``` r

out_bin <- ctt(data = dat_bin)

out_bin$item    # difficulty, discrimination_raw/corrected, alpha-with-item-removed, flag
#>    item cats difficulty discrimination_raw discrimination_corrected
#> 1    V1    2      0.887              0.309                    0.205
#> 2    V2    2      0.837              0.369                    0.250
#> 3    V3    2      0.813              0.479                    0.365
#> 4    V4    2      0.763              0.489                    0.364
#> 5    V5    2      0.790              0.329                    0.195
#> 6    V6    2      0.697              0.403                    0.257
#> 7    V7    2      0.637              0.514                    0.375
#> 8    V8    2      0.630              0.498                    0.356
#> 9    V9    2      0.560              0.376                    0.215
#> 10  V10    2      0.500              0.497                    0.349
#> 11  V11    2      0.467              0.470                    0.319
#> 12  V12    2      0.457              0.459                    0.307
#> 13  V13    2      0.407              0.382                    0.223
#> 14  V14    2      0.313              0.395                    0.247
#> 15  V15    2      0.323              0.386                    0.236
#>    alpha_removed flag
#> 1          0.668     
#> 2          0.664     
#> 3          0.651     
#> 4          0.649     
#> 5          0.670     
#> 6          0.663     
#> 7          0.646     
#> 8          0.649     
#> 9          0.669     
#> 10         0.650     
#> 11         0.654     
#> 12         0.656     
#> 13         0.668     
#> 14         0.664     
#> 15         0.666
out_bin$alpha   # test-level alpha, alpha_std, sem, mean difficulty/discrimination
#>   n_examinee n_item alpha alpha_std   sem mean_difficulty
#> 1        300     15 0.675     0.675 1.648           0.605
#>   mean_discrimination_raw mean_discrimination_corrected
#> 1                   0.424                         0.284
```

Difficulty for item $`j`$ is the mean item score divided by its maximum
possible score, $`p_j = \bar{X}_j / m_j`$ - the familiar
proportion-correct index for a dichotomous item. Discrimination is the
Pearson correlation between the item score and the total score; the
**raw** (uncorrected) form correlates an item with the total score that
includes its own contribution, while the **corrected** (item-excluded)
form removes it. Both are always computed and reported as separate
columns; the `correct` argument only selects which of the two feeds the
discrimination-flagging criterion (`crit.dis`).

### A mixed-format example

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md)’s difficulty
definition, $`p_j = \bar{X}_j / m_j`$, keeps difficulty on a common 0-1
scale regardless of the number of score categories, so it generalizes
directly to a mixed dichotomous/polytomous test such as `dat_mixed`:

``` r

out_mixed <- ctt(data = dat_mixed, item.id = x$id, cats = x$cats)

# the 5 polytomous (GRM) items, for comparison with the ctt_distr() output above
out_mixed$item[out_mixed$item$cats > 2, ]
#>    item cats difficulty discrimination_raw discrimination_corrected
#> 39 CFR1    5      0.750              0.650                    0.584
#> 40 CFR2    5      0.452              0.623                    0.537
#> 53 AFR1    5      0.418              0.585                    0.494
#> 54 AFR2    5      0.712              0.563                    0.480
#> 55 AFR3    5      0.454              0.466                    0.357
#>    alpha_removed flag
#> 39         0.876     
#> 40         0.879     
#> 53         0.880     
#> 54         0.879     
#> 55         0.885
```

### Reliability: raw alpha and standardized alpha

Two forms of Cronbach’s alpha are always computed and returned in
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md)’s `$alpha`
element. Raw alpha uses the standard variance-based formula
``` math
\alpha = \frac{k}{k - 1}\left(1 - \frac{\sum_{j} \sigma_j^2}{\sigma_X^2}\right),
```
and reflects the reliability of the actual (unweighted) total score
obtained by simply summing the item scores - the score most tests
actually use for reporting and decisions (Cronbach, 1951). Standardized
alpha instead first standardizes every item to unit variance before
combining them,
``` math
\alpha_{std} = \frac{k \bar{r}}{1 + (k - 1)\bar{r}},
```
where $`\bar{r}`$ is the average pairwise inter-item correlation; it
differs meaningfully from raw alpha mainly when items vary substantially
in scale or format, such as the mix of dichotomous and polytomous items
in `dat_mixed` (Osburn, 2000):

``` r

out_bin$alpha[, c("alpha", "alpha_std")]     # dichotomous-only: close together
#>   alpha alpha_std
#> 1 0.675     0.675
out_mixed$alpha[, c("alpha", "alpha_std")]   # mixed-format: can differ more
#>   alpha alpha_std
#> 1 0.882     0.901
```

The standard error of measurement, $`SEM = SD(X)\sqrt{1 - \alpha}`$
(using raw alpha), is reported alongside alpha as the standard CTT
expression of measurement precision.

### `print()` and `summary()` methods

[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) returns an
object of class `"ctt"`.
[`print.ctt()`](https://hwangQ.github.io/irtQ/reference/print.ctt.md)
displays a condensed report;
[`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md)
(paired with
[`print.summary.ctt()`](https://hwangQ.github.io/irtQ/reference/print.summary.ctt.md))
prints the full item-level table, the test-level reliability summary,
and the total-score frequency distribution - analogous to how
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) pairs
with
`print.est_irt()`/[`summary.est_irt()`](https://hwangQ.github.io/irtQ/reference/summary.md)
elsewhere in irtQ:

``` r

out_bin              # print.ctt(): condensed report
#> 
#> Call:
#> ctt(data = dat_bin)
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
summary(out_bin)     # summary.ctt(): full report
#> 
#> Call:
#> ctt(data = dat_bin)
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
```

------------------------------------------------------------------------

## Summary

| Function | Input | Key output | Used by |
|----|----|----|----|
| [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md) | Raw option choices + answer key | 0/1 scored response matrix | Output feeds [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md), [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md), [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md); called by [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md) (selected-response mode) |
| [`freq_score()`](https://hwangQ.github.io/irtQ/reference/freq_score.md) | Vector of total (raw) scores | Frequency / percentage / cumulative-percentage table | [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) |
| [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md) | Raw option data + key, or a scored matrix | Option/category response distribution, point-biserial correlations, distractor flags | Item and distractor review |
| [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) | Scored item response matrix | Item-level statistics, test-level reliability (alpha/SEM), total-score frequency distribution | [`print.ctt()`](https://hwangQ.github.io/irtQ/reference/print.ctt.md), [`summary.ctt()`](https://hwangQ.github.io/irtQ/reference/summary.ctt.md) |

------------------------------------------------------------------------

## References

Crocker, L., & Algina, J. (1986). *Introduction to classical and modern
test theory*. Holt, Rinehart; Winston.

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika*, *16*(3), 297–334.
<https://doi.org/10.1007/BF02310555>

Osburn, H. G. (2000). Coefficient alpha and related internal consistency
reliability coefficients. *Psychological Methods*, *5*(3), 343–355.
<https://doi.org/10.1037/1082-989X.5.3.343>
