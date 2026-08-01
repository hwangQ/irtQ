# Test-Level Reliability Summary (Cronbach's Alpha)

Computes a test-level classical test theory (CTT) reliability summary -
Cronbach's alpha (both the raw and standardized forms), the standard
error of measurement (SEM), and the average item difficulty and
discrimination - from scored item response data.

## Usage

``` r
ctt_alpha(data, item.id = NULL, cats = NULL, correct = FALSE, missing = NA)
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
  `TRUE`, flagging (in
  [`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md)'s
  `flag` column, `ctt_alpha()`'s summary discrimination figure, and, in
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

## Value

A one-row data frame containing:

- n_examinee:

  number of examinees included (after listwise deletion).

- n_item:

  number of items.

- alpha:

  Cronbach's alpha, raw (covariance-based) form.

- alpha_std:

  Cronbach's alpha, standardized (correlation-based) form - equivalent
  to computing raw alpha after first standardizing every item to unit
  variance. See **Details** for when this differs meaningfully from
  `alpha`.

- sem:

  the standard error of measurement (based on raw alpha).

- mean_difficulty:

  average item difficulty.

- mean_discrimination_raw:

  average raw (uncorrected) item-total correlation across items.

- mean_discrimination_corrected:

  average corrected (item-excluded) item-total correlation across items.

## Details

Two forms of Cronbach's alpha are always computed and reported as
separate columns:

Raw alpha (`alpha`) uses the standard variance-based formula,
`alpha = (k / (k - 1)) * (1 - sum(item variances) / total-score variance)`,
where k is the number of items. This formula is a general reliability
coefficient that applies unchanged to dichotomous and polytomous item
scores alike, and it reflects the reliability of the actual (unweighted)
total score obtained by simply summing the item scores - the score most
tests actually use for reporting and decisions. In everyday terms, raw
alpha asks: "if every item kept its own natural scale and spread, how
consistently do these items agree with each other?"

Standardized alpha (`alpha_std`) instead first standardizes every item
to the same scale (unit variance) before combining them, using the
equivalent formula `alpha_std = (k * r_bar) / (1 + (k - 1) * r_bar)`,
where `r_bar` is the average pairwise correlation among all items. In
everyday terms, standardized alpha asks: "if every item counted equally
regardless of how much it happens to vary in this particular sample, how
consistently would these items agree with each other?" Because it
removes the influence of any single item's variance, standardized alpha
is most useful when items differ substantially in scale or format (e.g.,
a mix of dichotomous and polytomous items with very different score
ranges); when all items share the same scale and format (as with a
dichotomous selected-response test scored 0/1), raw and standardized
alpha are typically close, and raw alpha remains the more directly
interpretable of the two since it matches the reliability of the score
actually used in practice. See Cronbach (1951) for the original
derivation of coefficient alpha, and Osburn (2000) for a discussion
contrasting the raw (covariance-based) and standardized
(correlation-based) forms.

The standard error of measurement (SEM) is computed as
`SEM = SD(total score) * sqrt(1 - alpha)` (using raw alpha), following
the standard CTT relationship between test reliability and measurement
precision.

Average difficulty and average discrimination are the simple means of
the per-item difficulty and discrimination values computed by
[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md) (any
`NA` per-item value, e.g. from a constant item, is excluded via
`na.rm = TRUE`). As in
[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md),
both the raw (uncorrected) and corrected (item-excluded) item-total
correlations are always averaged and reported as separate columns.
`ctt_alpha()` itself has no flagging step, so `correct` has no further
effect here beyond being passed through to
[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md) for
internal consistency.

Because `mean_difficulty`/`mean_discrimination_raw`/
`mean_discrimination_corrected` are averaged with `na.rm = TRUE`, they
inherit
[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md)'s
`cats`-auto-inference caveat: an item whose observed score range never
reaches its true maximum (most notably, an item that every examinee
scores 0 on) will have its difficulty silently excluded from the average
rather than contributing a `0`, which can bias the reported mean upward
in extreme/degenerate samples. This is unlikely to matter with a
reasonably large, non-degenerate sample, but supplying `cats` explicitly
avoids the issue entirely. If every item in `data` is degenerate in this
way, these means will be `NaN` (mean of an empty/all- `NA` vector)
rather than `NA`.

## References

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika*, *16*(3), 297-334.
https://doi.org/10.1007/BF02310555

Osburn, H. G. (2000). Coefficient alpha and related internal consistency
reliability coefficients. *Psychological Methods*, *5*(3), 343-355.
https://doi.org/10.1037/1082-989X.5.3.343

## See also

[`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md),
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
