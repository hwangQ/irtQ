# Test-Level Reliability Summary (Cronbach's Alpha)

Computes a test-level classical test theory (CTT) reliability summary
from scored item response data: Cronbach's alpha (raw and standardized
forms), the standard error of measurement (SEM), and the average item
difficulty and discrimination.

## Usage

``` r
ctt_alpha(data, item.id = NULL, cats = NULL, correct = FALSE, missing = NA)
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

  Logical. Accepted for consistency with
  [`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md)
  and [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md). Both
  the raw and corrected mean item-total correlations are always reported
  and `ctt_alpha()` has no flagging step, so `correct` does not change
  the result.

- missing:

  A value indicating missing responses in `data`, analogous to the
  `missing` argument in
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) and
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
  Any cell equal to `missing` is recoded to `NA` before analysis.
  Default is `NA`. Examinees with any remaining missing item response
  are excluded listwise from all statistics computed by this function.

## Value

A one-row data frame with the following columns. Values other than the
counts are rounded to three decimal places.

- n_examinee:

  number of examinees included (after listwise deletion).

- n_item:

  number of items.

- alpha:

  Cronbach's alpha, raw (covariance-based) form; `NA` when the total
  score has zero variance.

- alpha_std:

  Cronbach's alpha, standardized (correlation-based) form: raw alpha
  computed after standardizing every item with nonzero variance to unit
  variance, with a constant item kept in the item count as in `alpha`;
  `NA` when it is undefined. It equals
  `k * r_bar / (1 + (k - 1) * r_bar)` when no item is constant. See
  **Details** for when it differs meaningfully from `alpha`.

- sem:

  the standard error of measurement (based on raw alpha); `NA` when
  `alpha` is `NA`.

- mean_difficulty:

  average item difficulty.

- mean_discrimination_raw:

  average raw (uncorrected) item-total correlation across items with a
  defined value; `NA` when there is none.

- mean_discrimination_corrected:

  average corrected (item-excluded) item-total correlation across items
  with a defined value; `NA` when there is none.

## Details

Two forms of Cronbach's alpha are always computed and reported as
separate columns:

Raw alpha (`alpha`) uses the standard variance-based formula,
`alpha = (k / (k - 1)) * (1 - sum(item variances) / total-score variance)`,
where k is the number of items. This formula is a general reliability
coefficient that applies unchanged to dichotomous and polytomous item
scores alike, and it reflects the reliability of the unweighted total
score obtained by summing the item scores, which is the score most tests
use for reporting and decisions. In everyday terms, raw alpha asks: "if
every item kept its own natural scale and spread, how consistently do
these items agree with each other?"

Standardized alpha (`alpha_std`) is raw alpha computed after
standardizing every item with nonzero variance to unit variance. When no
item is constant, it equals
`alpha_std = (k * r_bar) / (1 + (k - 1) * r_bar)`, where `r_bar` is the
average pairwise correlation among the items. A constant item has no
defined correlation and cannot be standardized, so it stays in the item
count k with zero variance, as in raw alpha, and `r_bar` is the average
correlation among the items that vary. When only one item varies,
`alpha_std` is 0, as is raw alpha; it is `NA` when no item varies or the
standardized total score is constant. In everyday terms, standardized
alpha asks: "if every item counted equally regardless of how much it
happens to vary in this particular sample, how consistently would these
items agree with each other?" Because it removes the influence of any
single item's variance, standardized alpha is most useful when items
differ substantially in scale or format (e.g., a mix of dichotomous and
polytomous items with very different score ranges); when all items share
the same scale and format (as with a dichotomous selected-response test
scored 0/1), raw and standardized alpha are typically close, and raw
alpha remains the more directly interpretable of the two since it
matches the reliability of the score actually used in practice. See
Cronbach (1951) for the original derivation of coefficient alpha, and
Osburn (2000) for a discussion contrasting the raw (covariance-based)
and standardized (correlation-based) forms.

The standard error of measurement (SEM) is computed as
`SEM = SD(total score) * sqrt(1 - alpha)` (using raw alpha), following
the standard CTT relationship between test reliability and measurement
precision.

Average difficulty and average discrimination are the simple means of
the per-item values computed by
[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md),
ignoring `NA` values (e.g., the discrimination of a constant item). As
in [`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md),
the raw and corrected item-total correlations are averaged separately
and both are reported. `ctt_alpha()` has no flagging step, so `correct`
does not change the result.

Average difficulty depends on `cats`. When `cats` is not supplied, it is
inferred for each item as the observed maximum score plus one, with a
minimum of two, as in
[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md). An
item that every examinee scores 0 on therefore enters the average with a
difficulty of 0. For a polytomous item whose highest score category is
not observed in the sample, `cats` is inferred too small and the item's
difficulty is biased upward, so supply `cats` explicitly whenever the
maximum possible score may not have been observed.

## References

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika, 16*(3), 297-334.
[doi:10.1007/BF02310555](https://doi.org/10.1007/BF02310555) .

Osburn, H. G. (2000). Coefficient alpha and related internal consistency
reliability coefficients. *Psychological Methods, 5*(3), 343-355.
[doi:10.1037/1082-989X.5.3.343](https://doi.org/10.1037/1082-989X.5.3.343)
.

## See also

[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md),
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
