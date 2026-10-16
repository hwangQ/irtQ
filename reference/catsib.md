# CATSIB DIF Detection Procedure

This function performs DIF analysis on items using the CATSIB procedure
(Nandakumar & Roussos, 2004), a modified version of SIBTEST (Shealy &
Stout, 1993). The CATSIB procedure is suitable for computerized adaptive
testing (CAT) environments. In CATSIB, examinees are matched on
IRT-based ability estimates that have been adjusted using a regression
correction method (Shealy & Stout, 1993) to reduce statistical bias in
the CATSIB statistic caused by impact.

## Usage

``` r
catsib(
  x = NULL,
  data,
  score = NULL,
  se = NULL,
  group,
  focal.name,
  item.skip = NULL,
  D = 1,
  n.bin = c(80, 10),
  min.binsize = 3,
  max.del = 0.075,
  weight.group = c("comb", "foc", "ref"),
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  max.iter = 10,
  min.resp = NULL,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  rho.screen = TRUE,
  ...
)

# S3 method for class 'catsib'
print(x, digits = max(2L, getOption("digits") - 5L), ...)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.). See
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function. `x` is needed to estimate the abilities internally, that is,
  when `score = NULL` or `purify = TRUE`. It can be `NULL` when `score`
  and `se` are supplied and `purify = FALSE`; the items are then labeled
  `item.1`, `item.2`, and so on. Default is `NULL`.

- data:

  A matrix of examinees' scored responses (0 or 1). Rows represent
  examinees and columns represent items, in the order of the rows of `x`
  when `x` is supplied.

- score:

  A numeric vector of examinees' ability estimates (theta values), in
  the same row order as `data`, with their standard errors in `se`. The
  estimates should be based on the item parameters in `x`. If `NULL`,
  abilities and standard errors are estimated internally with
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  using `method`, `range`, `norm.prior`, `nquad`, `weights`, and
  `ncore`, and `x` must be supplied. When `purify = TRUE`, abilities and
  standard errors are re-estimated internally at every purification
  iteration, so a supplied `score` is used only in the initial analysis.
  A missing value in `score` excludes the examinee from the analysis,
  also during purification, with a warning. Default is `NULL`.

- se:

  A numeric vector of the standard errors of the ability estimates in
  `score`, in the same order. The standard errors are used to estimate
  the reliability of the ability estimates in the regression correction
  (see **Details**). It must be supplied when `score` is supplied, and
  it is ignored when `score = NULL`. Default is `NULL`.

- group:

  A numeric or character vector indicating examinees' group membership.
  The length of the vector must match the number of rows in the response
  data matrix.

- focal.name:

  A single numeric or character value specifying the focal group. For
  instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating the
  focal group, set `focal.name = 1`.

- item.skip:

  A numeric vector of item positions (columns of `data`, or rows of
  `x`), or a logical vector with one value per item, identifying the
  items to exclude from the analysis. Skipped items still contribute to
  the ability estimates, but their results (`beta`, `se`, `z.beta`, `p`,
  and the contingency table) are set to `NA`, and they are never flagged
  or removed during purification; their sample sizes (`n.ref`, `n.foc`,
  and `n.total`) are still reported, as in
  [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) and
  [`crdif()`](https://hwangQ.github.io/irtQ/reference/crdif.md). If
  `NULL`, all items are analyzed. Default is `NULL`.

- D:

  A scaling factor in IRT models, used only when abilities are estimated
  internally. Default is `1`.

- n.bin:

  A numeric vector of two positive integers giving the maximum and
  minimum numbers of equal-width bins (intervals) of the corrected
  ability scale, in this order (the maximum must not be smaller than the
  minimum). For each item, the number of bins starts at the maximum and
  decreases by one until the condition set by `max.del` is met; if it is
  not met with the minimum number, the minimum number is used. Default
  is `c(80, 10)`. See **Details**.

- min.binsize:

  An integer of 2 or more specifying the minimum number of examinees of
  each group in a bin. A bin is used in the computation of
  \\\hat{\beta}\\ only if it contains at least `min.binsize` examinees
  of both the reference and focal groups who responded to the item.
  Default is `3`, the value used by Nandakumar and Roussos (2004). See
  **Details**.

- max.del:

  A numeric value specifying the maximum proportion of the examinees of
  each group who responded to the item that may be excluded because they
  fall in bins that are not used (see `min.binsize`). The number of bins
  is reduced until the excluded proportion is at most `max.del` in both
  groups. It must be at least 0 and less than 1. Default is `0.075`. See
  **Details**.

- weight.group:

  A character string specifying the target ability distribution used to
  weight the bins in \\\hat{\beta}\\ and its standard error: `"comb"`
  for the combined distribution of the reference and focal groups,
  `"foc"` for the focal group, or `"ref"` for the reference group. The
  weight of a bin is the proportion of the examinees of the target
  distribution in the bins used who fall in that bin. Default is
  `"comb"`.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) of the
  two-sided test of the CATSIB statistic. An item is flagged when its
  p-value is less than or equal to `alpha`. The comparison uses
  unrounded p-values, although the p-values reported in `dif_stat` are
  rounded to four decimal places. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- max.iter:

  A positive integer specifying the maximum number of iterations allowed
  for the purification process. Default is `10`. Each iteration removes
  one item, so `max.iter` should be at least as large as the number of
  items expected to be flagged. If the limit is reached while flagged
  items remain, a warning is issued, the `complete` element of the
  purification results is `FALSE`, and the items flagged in the last
  iteration are added to the flagged items. This argument is not passed
  to
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md).

- min.resp:

  A single nonnegative number, usually a positive integer, specifying
  the minimum number of item responses that an examinee must have to be
  used in the analysis. All responses of examinees with fewer than
  `min.resp` (but at least one) responses are set to `NA` before the
  ability estimation, in the initial analysis and at every purification
  iteration, also when `score` is supplied, and these examinees are
  excluded from the analysis. A warning reports the number of examinees
  excluded in the initial analysis, and one warning at the end of the
  purification reports the examinees excluded during it. A large
  `min.resp` can leave a group with fewer than two ability estimates
  within `range` in a later purification iteration, which stops the
  purification with a warning. If `NULL`, no minimum is applied. Default
  is `NULL`. See **Details** for more information.

- method:

  A character string indicating the scoring method used to estimate the
  abilities when `score = NULL` and in every purification iteration.
  Available options are:

  - `"ML"`: Maximum likelihood estimation

  - `"MLF"`: Maximum likelihood estimation with fences

  - `"WL"`: Weighted likelihood estimation (Warm, 1989)

  - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)

  - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)

  The methods `"EAP.SUM"` and `"INV.TCC"` of
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  are not supported, because they do not give an ability estimate for
  each response pattern. `"ML"` or `"WL"` is recommended, because the
  regression correction assumes unbiased ability estimates (see
  **Details**). Default is `"ML"`.

- range:

  A numeric vector of length two specifying the lower and upper bounds
  of the ability scale. It is used when abilities are estimated
  internally with the scoring methods `"ML"`, `"MLF"`, `"WL"`, and
  `"MAP"`. In addition, ability estimates at or beyond these bounds are
  excluded when the mean, the variance, and the error variance of the
  ability estimates are computed for the regression correction (see
  **Details**); this also applies when `score` is supplied, so `range`
  should match the bounds used to obtain `score`. Default is `c(-5, 5)`.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. These values are used to
  generate the Gaussian quadrature points and weights. Ignored if
  `method` is `"ML"` or `"WL"`. Default is `c(0, 1)`.

- nquad:

  An integer indicating the number of Gaussian quadrature points to be
  generated from the normal prior distribution. Used only when `method`
  is `"EAP"`. Ignored for `"ML"`, `"WL"`, and `"MAP"`. Default is 41.

- weights:

  A two-column matrix or data frame containing the quadrature points (in
  the first column) and their corresponding weights (in the second
  column) for the latent variable prior distribution. The weights and
  points can be conveniently generated using the function
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).

  If `NULL` and `method = "EAP"`, default quadrature values are
  generated based on the `norm.prior` and `nquad` arguments. Ignored if
  `method` is `"ML"`, `"WL"`, or `"MAP"`.

- ncore:

  An integer specifying the number of logical CPU cores used by
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  when abilities are estimated internally. The test statistics
  themselves are not computed in parallel. Default is `1`.

- verbose:

  Logical. If `TRUE`, the progress of the purification procedure is
  printed to the console. Default is `TRUE`.

- rho.screen:

  Logical. If `TRUE`, the ability estimates whose squared standard error
  is equal to or greater than the variance of the ability estimates of
  their group are left out of the computation of the mean, the variance,
  and the error variance that define \\\hat{\rho}^2_G\\ in the
  regression correction. These examinees are still matched and used to
  compute \\\hat{\beta}\\. If `FALSE`, \\\hat{\rho}^2_G\\ is computed
  from all ability estimates within `range`, as in Nandakumar and
  Roussos (2004, Equation 7). The screening is not part of the original
  procedure and is added in irtQ. Default is `TRUE`. See **Details**.

- ...:

  Additional arguments passed to
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  when abilities are estimated internally, for example `tol`. The
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  arguments `se`, `max.iter`, and `missing` cannot be passed this way,
  because `catsib()` uses arguments with these names for its own
  purposes.

- digits:

  Not used; kept for consistency with the other print methods.

## Value

This function returns an object of class `"catsib"`, which is a list
with the following five components:

- no_purify:

  A list containing the results of the DIF analysis without applying a
  purification procedure. This list includes:

  dif_stat

  :   A data frame with one row per item and the columns `id`, `beta`,
      `se`, `z.beta`, `p`, `n.ref`, `n.foc`, and `n.total`: the item ID,
      the CATSIB statistic \\\hat{\beta}\\, its standard error, the
      standardized statistic, the two-sided p-value, and the numbers of
      examinees of the reference and focal groups in the bins used to
      compute \\\hat{\beta}\\, and their sum. A positive \\\hat{\beta}\\
      indicates DIF in favor of the reference group. All values are
      rounded to four decimal places.

  dif_item

  :   A numeric vector of the positions (columns of `data`) of the
      flagged items, or `NULL` if no item is flagged.

  contingency

  :   A list with one element per item, named by item ID. Each element
      is a data frame with one row per bin used to compute
      \\\hat{\beta}\\ and the columns `bin` (the interval of corrected
      ability estimates), `n.ref`, `prop.ref`, `var.ref`, `n.foc`,
      `prop.foc`, and `var.foc` (the number of examinees, the proportion
      of correct responses, and the sample variance of the item
      responses of each group in the bin), `n.total`, `weight` (the
      weight of the bin), `beta` (the weighted difference in the
      proportions of correct responses), and `var.beta` (the
      contribution of the bin to the squared standard error), followed
      by a row whose `bin` is `Total`, with the totals of `n.ref`,
      `n.foc`, `n.total`, `beta`, and `var.beta`. The element is `NA`
      for a skipped item.

- purify:

  A logical value indicating whether purification was requested (the
  `purify` argument).

- with_purify:

  A list with the results of the purification procedure. All elements
  are `NULL` when `purify = FALSE`. The elements are:

  dif_stat

  :   A data frame with the same columns as `no_purify$dif_stat` and an
      additional column `n.iter`. For an item removed during
      purification, the row reports the statistics from the iteration in
      which the item was removed, and `n.iter` is that iteration minus 1
      (0 for an item removed on the basis of the initial analysis). For
      the other items, the row reports the statistics from the last
      iteration, and `n.iter` is the number of that iteration. If no
      item is flagged in the initial analysis, this is
      `no_purify$dif_stat` with `n.iter = 0`.

  dif_item

  :   A numeric vector of the positions (rows of `x`) of the items
      flagged by the CATSIB statistic, sorted in ascending order. It
      contains the items removed during purification and, if the
      purification stops before it is complete (see `complete`) with
      items left, the items flagged in the last iteration (iteration
      `n.iter`). It contains all items if every item is flagged. `NULL`
      if no item is flagged in the initial analysis.

  n.iter

  :   The number of purification iterations performed (0 if no item is
      flagged in the initial analysis). An iteration after which the
      statistics cannot be recomputed (see **Details**) is not counted.

  complete

  :   A logical value. `TRUE` if the procedure stopped because no
      remaining item was flagged (including the case in which no item is
      flagged in the initial analysis), and `FALSE` if it stopped
      because `max.iter` was reached, because every item was flagged,
      because a group had fewer than two ability estimates within
      `range`, or because no remaining item had a usable bin.

  contingency

  :   A list with the same structure as `no_purify$contingency`. For an
      item removed during purification, it holds the table from the
      iteration in which the item was removed; for the other items, it
      holds the table from the last iteration. If no item is flagged in
      the initial analysis, this is `no_purify$contingency`.

- alpha:

  A numeric value indicating the significance level (\\\alpha\\) used
  for the tests.

- call:

  The matched function call.

## Details

CATSIB (Nandakumar & Roussos, 2004) estimates the amount of DIF in an
item as \\\hat{\beta} = \sum\_{k} (\hat{P}\_{R,k} - \hat{P}\_{F,k})
\hat{p}\_k\\, where \\\hat{P}\_{R,k}\\ and \\\hat{P}\_{F,k}\\ are the
proportions of correct responses of the reference and focal groups in
bin k of the corrected ability scale, and \\\hat{p}\_k\\ is the weight
of bin k (see `weight.group`). A positive \\\hat{\beta}\\ indicates DIF
in favor of the reference group. The standard error of \\\hat{\beta}\\
is estimated from the sample variances of the item responses in each bin
and group (Nandakumar & Roussos, 2004, Equation 10), and the test
statistic \\\hat{\beta} / \mathrm{SE}(\hat{\beta})\\ (`z.beta`) is
compared with the standard normal distribution in a two-sided test.
Examinees are matched on \\\hat{\theta}^{\ast}\\, an estimate of the
expected value of \\\theta\\ given \\\hat{\theta}\\ obtained by the
regression correction described below. The range of
\\\hat{\theta}^{\ast}\\ is divided into K equal-width bins, and
examinees are classified into the bins by their \\\hat{\theta}^{\ast}\\.
A bin is used only if it contains at least `min.binsize` examinees of
both groups who responded to the item; Nandakumar and Roussos (2004)
required at least three, which is the default.

To choose the number of bins K, `catsib()` follows the rule of
Nandakumar and Roussos (2004). For each item, K starts at a large value
(80 by default), and if more than 7.5% of the examinees of either group
who responded to the item would be excluded because they fall in bins
that are not used, K is reduced by one, until at most 7.5% of each group
is excluded. Nandakumar and Roussos (2004) did not allow fewer than 20
bins, but they reported that a minimum of 10 bins reduced the estimation
bias caused by excluding examinees from sparse bins. `catsib()`
therefore uses a minimum of 10 bins by default; if the condition is not
met with the minimum number, the minimum number is used. The maximum and
minimum numbers of bins are set by `n.bin`, and the maximum proportion
of excluded examinees by `max.del`.

Li and Stout (1996) and Nandakumar and Roussos (2004) used the combined
distribution of the two groups as the target ability distribution, which
is the default of `weight.group`. See Nandakumar and Roussos (2004) for
further details about the CATSIB method.

Because the ability distributions of the reference and focal groups can
differ (impact), matching examinees directly on \\\hat{\theta}\\ would
bias \\\hat{\beta}\\ toward DIF in favor of the group with the higher
mean ability. CATSIB therefore applies a regression correction
separately in each group G (Nandakumar & Roussos, 2004, Equation 7):
\\\hat{\theta}^{\ast} = \bar{\hat{\theta}}\_G + \hat{\rho}^2_G
(\hat{\theta} - \bar{\hat{\theta}}\_G)\\, where
\\\bar{\hat{\theta}}\_G\\ is the mean of the ability estimates in group
G and \\\hat{\rho}^2_G = 1 - \hat{\sigma}^2\_{e,G} /
\hat{\sigma}^2\_{\hat{\theta},G}\\ is the estimated reliability of the
ability estimates in group G. Here \\\hat{\sigma}^2\_{e,G}\\ is the mean
of the squared standard errors (`se`) and
\\\hat{\sigma}^2\_{\hat{\theta},G}\\ is the variance of the ability
estimates in group G. Ability estimates at or beyond a bound of `range`
are excluded when these quantities are computed, because their standard
errors are not meaningful; this also applies when `score` is supplied,
so `range` should match the bounds used to obtain `score`.

Because \\\hat{\sigma}^2\_{e,G}\\ is a mean of squared standard errors,
a few examinees with very large standard errors, typically examinees
with very low or very high scores whose abilities are poorly determined,
can make \\\hat{\rho}^2_G\\ much smaller than the reliability of the
other ability estimates, or even negative. The corrected ability
estimates then shrink too much toward the group mean, which distorts the
matching and can inflate the Type I error rate, particularly during
purification. When `rho.screen = TRUE` (the default), the examinees
whose squared standard error is greater than or equal to
\\\hat{\sigma}^2\_{\hat{\theta},G}\\ of their group are left out of the
computation of the mean, the variance, and the error variance that
define \\\hat{\rho}^2_G\\, unless fewer than two estimates would be
left. For these examinees the individual reliability \\1 - SE^2 /
\hat{\sigma}^2\_{\hat{\theta},G}\\ is zero or negative, that is, their
ability estimates carry no information about their abilities beyond the
spread of the group. They still receive corrected ability estimates and
are used in the matching and in the computation of \\\hat{\beta}\\. This
screening is not part of the original procedure (Nandakumar & Roussos,
2004) and is added in irtQ; set `rho.screen = FALSE` to compute
\\\hat{\rho}^2_G\\ from all ability estimates within `range` as in the
original procedure. Even with the screening, the Type I error rate of
CATSIB can be higher than the nominal level, and in data without DIF the
purification can increase the number of false positives. In addition,
\\\hat{\rho}^2_G\\ is bounded to the interval \\\[0.05, 1\]\\, and a
warning is issued when its unbounded value is below 0.05. A negative
value would reverse the order of the ability estimates within a group,
and a value near zero would place almost all corrected ability estimates
of a group in one bin. This bound is not part of the original procedure
and is added in irtQ. Even with the bound, the corrected ability
estimates of the two groups may not overlap when the groups differ in
ability; then no bin contains enough examinees of both groups for an
item, \\\hat{\beta}\\, its standard error, `z.beta`, and the p-value are
`NA`, one warning names the items, and the items are not flagged.

The regression correction assumes that the ability estimates are
unbiased, that is, that \\\hat{\theta}\\ equals \\\theta\\ plus an error
with mean zero. The Bayesian estimates (`"MAP"` and `"EAP"`) are already
shrunk toward the prior mean, so the correction shrinks them a second
time. Use `"ML"` or `"WL"`, especially when the abilities are
re-estimated during the purification.

CATSIB, like SIBTEST (Shealy & Stout, 1993), was designed and evaluated
for uniform DIF (Nandakumar & Roussos, 2004). Because \\\hat{\beta}\\
sums signed differences across bins, differences of opposite sign can
cancel, so its sensitivity to nonuniform or mixed DIF is limited, and
results for items suspected of nonuniform DIF should be interpreted with
caution.

Nandakumar and Roussos (2004) did not propose a purification procedure
for CATSIB, but `catsib()` can apply an iterative purification procedure
in the manner of Lim et al. (2022). It requires `x`. When
`purify = TRUE`, at each iteration the flagged item with the smallest
p-value is removed (the p-values are compared on the log scale), the
abilities and their standard errors are re-estimated from the remaining
items with the scoring method specified in `method`, and the CATSIB
statistics of the remaining items are recomputed. Supplied `score` and
`se` are therefore used only in the initial analysis. The procedure
stops when no remaining item is flagged or when `max.iter` iterations
have been performed. If no item is flagged in the initial analysis, no
iteration is performed. If every item is flagged, the procedure stops
with a warning, because no item is left to estimate the abilities, and
`complete` is `FALSE`. The procedure also stops with a warning and
`complete = FALSE` when, after an item is removed, a group has fewer
than two ability estimates within `range`, which the regression
correction requires, or no bin contains enough examinees of both groups
for any remaining item; the remaining items then keep the results of the
previous iteration, and the items flagged in that iteration are added to
the flagged items.

Scoring based on a small number of item responses can lead to large
standard errors, potentially reducing the accuracy of DIF detection with
CATSIB. The `min.resp` argument excludes such examinees. For example, if
`min.resp = 5`, examinees who responded to fewer than five items (but at
least one) have all their responses treated as missing (i.e., `NA`)
before the abilities are estimated, in the initial analysis and at every
purification iteration, and a warning reports their number in the
initial analysis. As a result, their ability estimates are missing, and
they are not used in the computation of the CATSIB statistics. When
`score` is supplied, the scores of these examinees are also treated as
missing. During purification, the count is based on the items that
remain after the flagged items are removed, and an examinee excluded
once stays excluded. If `min.resp = NULL`, a score is computed for any
examinee with at least one valid item response.

Examinees whose ability estimate is `NA`, such as examinees without any
item response or examinees with `NA` in a supplied `score`, are excluded
from the analysis, and a warning reports the number of examinees with
item responses who are excluded because of a missing value in `score`;
these examinees stay excluded during purification. During purification,
an examinee whose remaining responses are all removed with the flagged
items, or who has fewer than `min.resp` responses to the remaining
items, has no ability estimate and is excluded, and one warning at the
end reports the number of such examinees. The sample sizes `n.ref` and
`n.foc` count only the examinees in the bins used to compute
\\\hat{\beta}\\.

The `group` argument should contain exactly two distinct values, one for
the reference group and one for the focal group identified by
`focal.name`, and its length must equal the number of rows in `data`.
Examinees with `NA` in `group` are excluded from the analysis. The
function stops when the reference group or the focal group has fewer
than two ability estimates within `range`, because the regression
correction cannot be computed.

The print method of this class has an argument `digits` that is not
used; the statistics are printed rounded to three decimal places.

## Methods (by generic)

- `print(catsib)`: Print method for an object of class `catsib` (given
  as `x`). It prints the flagged items and the CATSIB statistics of the
  analyses without and with purification, rounded to three decimal
  places.

## References

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
in a microcomputer environment. *Applied Psychological Measurement,
6*(4), 431-444.
[doi:10.1177/014662168200600405](https://doi.org/10.1177/014662168200600405)
.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals
of item response theory*. Newbury Park, CA: Sage.

Li, H.-H., & Stout, W. (1996). A new procedure for detection of crossing
DIF. *Psychometrika, 61*(4), 647-677.
[doi:10.1007/BF02294041](https://doi.org/10.1007/BF02294041) .

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Nandakumar, R., & Roussos, L. (2004). Evaluation of the CATSIB DIF
procedure in a pretest setting. *Journal of Educational and Behavioral
Statistics, 29*(2), 177-199.
[doi:10.3102/10769986029002177](https://doi.org/10.3102/10769986029002177)
.

Shealy, R., & Stout, W. (1993). A model-based standardization approach
that separates true bias/DIF from group ability differences and detects
test bias/DTF as well as item bias/DIF. *Psychometrika, 58*(2), 159-194.
[doi:10.1007/BF02294572](https://doi.org/10.1007/BF02294572) .

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627) .

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
[`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
# Load required package
library("dplyr")
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

## Uniform DIF Detection
###############################################
# (1) Simulate data with true uniform DIF
###############################################

# Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select 36 3PLM items that are non-DIF
par_nstd <-
  bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
  dplyr::filter(.data$model == "3PLM") %>%
  dplyr::filter(dplyr::row_number() %in% 1:36) %>%
  dplyr::select(1:6)
par_nstd$id <- paste0("nondif", 1:36)

# Generate four new items to contain uniform DIF
difpar_ref <-
  shape_df(
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = rep(0.15, 4)),
    item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
  )

# Introduce uniform DIF in the focal group by shifting b-parameters
difpar_foc <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))

# Combine the 4 DIF and 36 non-DIF items for both reference and focal groups
# Therefore, the first four items now exhibit uniform DIF
par_ref <- rbind(difpar_ref, par_nstd)
par_foc <- rbind(difpar_foc, par_nstd)

# Generate true theta values
set.seed(123)
theta_ref <- rnorm(500, 0.0, 1.0)
theta_foc <- rnorm(500, 0.0, 1.0)

# Simulate response data
resp_ref <- simdat(par_ref, theta = theta_ref, D = 1)
resp_foc <- simdat(par_foc, theta = theta_foc, D = 1)
data <- rbind(resp_ref, resp_foc)

###############################################
# (2) Estimate item and ability parameters
#     using the aggregated data
###############################################

# Estimate item parameters
est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -25472.9403, Max-Change: 1.888841 EM iteration: 2, Loglike: -22654.4441, Max-Change: 0.356802 EM iteration: 3, Loglike: -22635.6173, Max-Change: 0.157888 EM iteration: 4, Loglike: -22633.7926, Max-Change: 0.074945 EM iteration: 5, Loglike: -22633.4312, Max-Change: 0.039693 EM iteration: 6, Loglike: -22633.3151, Max-Change: 0.022643 EM iteration: 7, Loglike: -22633.2668, Max-Change: 0.013669 EM iteration: 8, Loglike: -22633.2439, Max-Change: 0.008635 EM iteration: 9, Loglike: -22633.2322, Max-Change: 0.005655 EM iteration: 10, Loglike: -22633.2259, Max-Change: 0.00381 EM iteration: 11, Loglike: -22633.2225, Max-Change: 0.002624 EM iteration: 12, Loglike: -22633.2205, Max-Change: 0.001837 EM iteration: 13, Loglike: -22633.2194, Max-Change: 0.001303 EM iteration: 14, Loglike: -22633.2188, Max-Change: 0.000932 EM iteration: 15, Loglike: -22633.2184, Max-Change: 0.000671 EM iteration: 16, Loglike: -22633.2182, Max-Change: 0.000486 EM iteration: 17, Loglike: -22633.2180, Max-Change: 0.000353 EM iteration: 18, Loglike: -22633.2180, Max-Change: 0.000258 EM iteration: 19, Loglike: -22633.2179, Max-Change: 0.000188 EM iteration: 20, Loglike: -22633.2179, Max-Change: 0.000144 EM iteration: 21, Loglike: -22633.2179, Max-Change: 0.000116 EM iteration: 22, Loglike: -22633.2179, Max-Change: 9.7e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 1.9 seconds. 
est_par <- est_mod$par.est

# Estimate ability parameters using ML
theta_est <- est_score(x = est_par, data = data, method = "ML")
score <- theta_est$est.theta
se <- theta_est$se.theta

###############################################
# (3) Conduct DIF analysis
###############################################
# Create a vector of group membership indicators
# where '1' indicates the focal group
group <- c(rep(0, 500), rep(1, 500))

# (1) Compute the CATSIB statistic using provided scores and standard errors,
#     without purification (x can be NULL in this case)
dif_1 <- catsib(
  x = NULL, data = data, D = 1, score = score, se = se, group = group, focal.name = 1,
  weight.group = "comb", alpha = 0.05, missing = NA, purify = FALSE
)
print(dif_1)
#> 
#> Call:
#> catsib(x = NULL, data = data, score = score, se = se, group = group, 
#>     focal.name = 1, D = 1, weight.group = "comb", alpha = 0.05, 
#>     missing = NA, purify = FALSE)
#> 
#> DIF analysis using CATSIB method 
#> 
#>  1. Without purification 
#> 
#>   - Potential DIF Items: 
#>     1, 2, 3, 4, 11, 14 
#>   - Test Statistic: 
#> 
#>         id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   item.1   474   465     939  0.104 0.032  3.268 0.001 ***
#> 2   item.2   474   465     939  0.178 0.030  5.935 0.000 ***
#> 3   item.3   474   465     939  0.068 0.031  2.195 0.028   *
#> 4   item.4   474   465     939  0.170 0.028  6.093 0.000 ***
#> 5   item.5   474   465     939  0.015 0.032  0.474 0.636    
#> 6   item.6   474   465     939 -0.027 0.019 -1.394 0.163    
#> 7   item.7   474   465     939 -0.019 0.031 -0.605 0.545    
#> 8   item.8   474   465     939 -0.034 0.029 -1.175 0.240    
#> 9   item.9   474   465     939  0.002 0.030  0.083 0.934    
#> 10 item.10   474   465     939 -0.039 0.027 -1.451 0.147    
#> 11 item.11   474   465     939 -0.063 0.031 -2.038 0.042   *
#> 12 item.12   474   465     939 -0.002 0.032 -0.064 0.949    
#> 13 item.13   474   465     939  0.011 0.032  0.331 0.741    
#> 14 item.14   474   465     939 -0.066 0.028 -2.345 0.019   *
#> 15 item.15   474   465     939 -0.026 0.029 -0.913 0.361    
#> 16 item.16   474   465     939 -0.004 0.031 -0.114 0.909    
#> 17 item.17   474   465     939 -0.018 0.029 -0.611 0.541    
#> 18 item.18   474   465     939  0.025 0.029  0.847 0.397    
#> 19 item.19   474   465     939 -0.053 0.028 -1.855 0.064   .
#> 20 item.20   474   465     939  0.003 0.024  0.126 0.900    
#> 21 item.21   474   465     939 -0.025 0.028 -0.878 0.380    
#> 22 item.22   474   465     939 -0.033 0.032 -1.040 0.298    
#> 23 item.23   474   465     939 -0.018 0.018 -1.004 0.316    
#> 24 item.24   474   465     939  0.009 0.018  0.523 0.601    
#> 25 item.25   474   465     939 -0.024 0.022 -1.113 0.266    
#> 26 item.26   474   465     939 -0.011 0.030 -0.355 0.723    
#> 27 item.27   474   465     939 -0.029 0.029 -1.012 0.312    
#> 28 item.28   474   465     939 -0.031 0.031 -0.992 0.321    
#> 29 item.29   474   465     939 -0.019 0.025 -0.781 0.435    
#> 30 item.30   474   465     939 -0.026 0.021 -1.224 0.221    
#> 31 item.31   474   465     939 -0.018 0.030 -0.584 0.559    
#> 32 item.32   474   465     939 -0.031 0.026 -1.220 0.222    
#> 33 item.33   474   465     939  0.025 0.022  1.143 0.253    
#> 34 item.34   474   465     939  0.015 0.032  0.472 0.637    
#> 35 item.35   474   465     939 -0.021 0.031 -0.674 0.500    
#> 36 item.36   474   465     939  0.016 0.023  0.685 0.493    
#> 37 item.37   474   465     939 -0.012 0.021 -0.558 0.577    
#> 38 item.38   474   465     939  0.054 0.031  1.721 0.085   .
#> 39 item.39   474   465     939 -0.030 0.029 -1.046 0.295    
#> 40 item.40   474   465     939  0.026 0.031  0.856 0.392    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (2) Compute the CATSIB statistic with purification (x is required; not
#     run: it takes about 7 seconds because the abilities of the examinees
#     are re-estimated at every iteration)
#
#   dif_2 <- catsib(
#     x = est_par, data = data, D = 1, score = score, se = se, group = group,
#     focal.name = 1, weight.group = "comb", alpha = 0.05, missing = NA,
#     purify = TRUE
#   )
#   print(dif_2)
# }
```
