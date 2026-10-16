# Classical Test Theory Item Analysis

Computes classical test theory (CTT) item statistics for dichotomous or
polytomous item response data: difficulty, the item-total correlation
(discrimination), and Cronbach's alpha with the item removed. Items that
fall outside commonly used quality thresholds can optionally be flagged.

## Usage

``` r
ctt_item(
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
  ... for a polytomous/partial-credit item). An error is raised when a
  score is not a whole number between 0 and `cats[j] - 1`; a score
  within 1e-8 of a whole number is treated as that whole number.

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
  inferred from the observed maximum score in `data` plus one, with a
  minimum of two (i.e., `max(max(data[, j]) + 1, 2)`); supply `cats`
  explicitly whenever the maximum possible score may not have been
  observed in the sample. Each value must be a whole number of at least
  2.

- correct:

  Logical. Both the raw item-total correlation (the item score, or in
  [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md)
  an option or category indicator, correlated with a total score that
  includes the item) and the corrected item-total correlation
  (correlated with a total score that excludes the item) are always
  computed and reported as separate columns. This argument only selects
  which of the two is used for flagging: the corrected value if `TRUE`,
  the raw value if `FALSE` (default). It drives the `flag` column of
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md) and the
  `crit.distractor` check of
  [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md).

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
  `crit.dis`, or when their discrimination is undefined.

- crit.p:

  A numeric vector of length two giving the lower and upper difficulty
  bounds used for flagging: an item with difficulty below the first
  value is flagged as too difficult ("difficulty too low"), and an item
  with difficulty above the second value is flagged as too easy
  ("difficulty too high"). Default is `c(0.10, 0.95)`.

- crit.dis:

  A single numeric value giving the minimum acceptable discrimination
  (item-total correlation); an item whose discrimination is strictly
  below this value is flagged as poorly discriminating ("discrimination
  too low"). An item whose discrimination is undefined (a constant item,
  or a constant reference total score) is flagged "discrimination
  undefined". Default is `0.20`.

## Value

A list with two elements:

- item:

  A data frame with one row per item and the columns `item` (item ID),
  `cats` (number of score categories), `difficulty`,
  `discrimination_raw` (raw item-total correlation),
  `discrimination_corrected` (corrected item-total correlation),
  `alpha_removed` (alpha with the item removed), and, only when
  `flag = TRUE`, `flag` (the criteria the item meets, among "difficulty
  too low", "difficulty too high", "discrimination too low", and
  "discrimination undefined", separated by "; ", or `""` when none).
  Discrimination is flagged using the column selected by `correct`.
  Statistics are rounded to three decimal places.

- crit:

  A list with elements `crit.p`, `crit.dis`, and `correct`, the
  thresholds and the correlation used for flagging (`correct = TRUE` for
  the corrected item-total correlation, `FALSE` for the raw one).

## Details

Difficulty for item j is defined generally as the mean observed score
divided by the item's maximum possible score,
`mean(data[, j]) / (cats[j] - 1)`, computed after listwise deletion. For
a dichotomous item (`cats[j] = 2`), this reduces to the familiar
proportion-correct difficulty index. For a polytomous item, this
expresses the average score as a proportion of the maximum attainable
score, so that difficulty remains interpretable on the same 0-1 scale
regardless of the number of score categories.

Discrimination for item j is the Pearson correlation between the item
score and the total score, which for a dichotomous item is
mathematically equivalent to the point-biserial correlation. A
correlation is `NA` when the item score or the reference total score is
constant. Both the uncorrected (raw) item-total correlation and the
corrected (item-excluded) item-total correlation are always computed and
returned as separate columns; see, e.g., Crocker and Algina (1986) for
discussion of both conventions. This mirrors how some software reports
both side by side (e.g., `psych::alpha()` reports the raw item-total
correlation as `raw.r` and the corrected version as `r.drop`) rather
than defaulting to one or the other. The `correct` argument only selects
which of the two feeds the discrimination flagging criterion (see
`crit.dis`).

Alpha-with-item-removed for item j is Cronbach's alpha recomputed on the
remaining `ncol(data) - 1` items, using the same variance-based formula
used for an overall, test-level alpha
(`k / (k - 1) * (1 - sum(item variances) / total variance)`). The value
is `NA` when the total score of the remaining items has zero variance or
when fewer than two items remain (a two-item test), since alpha is then
undefined.

## References

Crocker, L., & Algina, J. (1986). *Introduction to classical and modern
test theory*. Holt, Rinehart and Winston.

## See also

[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
