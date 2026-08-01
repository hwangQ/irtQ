# Classical Test Theory Item Analysis

Computes traditional classical test theory (CTT) item statistics -
difficulty, an item-total correlation (discrimination), and Cronbach's
alpha with the item removed - for dichotomous or polytomous item
response data, along with optional flagging of items that fall outside
commonly used quality thresholds.

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
  item (or, in
  [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md),
  a response option/category) is correlated with the total score that
  includes its own contribution - and the corrected item-total
  correlation - excluding its own contribution - are always computed and
  reported as separate columns (see **Details**/**Value**). This
  argument only controls which of the two is used for flagging: if
  `TRUE`, flagging (in `ctt_item()`'s `flag` column,
  [`ctt_alpha()`](https://hwangQ.github.io/irtQ/reference/ctt_alpha.md)'s
  summary discrimination figure, and, in
  [`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md),
  the `crit.distractor` check) is based on the corrected value; if
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

A list with two elements:

- item:

  A data frame with one row per item, containing the item label, number
  of score categories, difficulty, the raw (uncorrected) item-total
  correlation (`discrimination_raw`), the corrected (item-excluded)
  item-total correlation (`discrimination_corrected`), and
  alpha-with-item-removed, plus a `flag` column (a character description
  of which criteria were triggered, or `""` when none were; present only
  when `flag = TRUE`). Flagging is based on whichever of
  `discrimination_raw`/`discrimination_corrected` is selected by
  `correct`.

- crit:

  A list echoing the `crit.p` and `crit.dis` thresholds used, for
  reference in downstream reporting.

## Details

Difficulty for item j is defined generally as the mean observed score
divided by the item's maximum possible score,
`mean(data[, j], na.rm = TRUE) / (cats[j] - 1)`. For a dichotomous item
(`cats[j] = 2`), this reduces to the familiar proportion-correct
difficulty index. For a polytomous item, this expresses the average
score as a proportion of the maximum attainable score, so that
difficulty remains interpretable on the same 0-1 scale regardless of the
number of score categories.

Discrimination for item j is the Pearson correlation between the item
score and the total score, which for a dichotomous item is
mathematically equivalent to the point-biserial correlation. Both the
uncorrected (raw) item-total correlation and the corrected
(item-excluded) item-total correlation are always computed and returned
as separate columns; see, e.g., Crocker and Algina (1986) for discussion
of both conventions. This mirrors how some software reports both side by
side (e.g., `psych::alpha()` reports the raw item-total correlation as
`raw.r` and the corrected version as `r.drop`) rather than defaulting to
one or the other. The `correct` argument only selects which of the two
feeds the discrimination flagging criterion (see `crit.dis`).

Alpha-with-item-removed for item j is Cronbach's alpha recomputed on the
remaining `ncol(data) - 1` items, using the same variance-based formula
used for an overall, test-level alpha
(`k / (k - 1) * (1 - sum(item variances) / total variance)`). A value of
`NA` is returned wherever a needed variance is zero (e.g., a constant
item, or fewer than two items remaining), since the relevant ratio is
then undefined.

## References

Crocker, L., & Algina, J. (1986). *Introduction to classical and modern
test theory*. Holt, Rinehart and Winston.

## See also

[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
