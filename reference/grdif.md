# Generalized IRT residual-based DIF detection framework for multiple groups (GRDIF)

This function computes three GRDIF statistics, \\GRDIF\_{R}\\,
\\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\, for analyzing differential item
functioning (DIF) among multiple groups (Lim et al., 2024). They are
designed mainly to detect uniform DIF, nonuniform DIF, and mixed DIF,
respectively.

## Usage

``` r
grdif(x, ...)

# Default S3 method
grdif(
  x,
  data,
  score = NULL,
  group,
  focal.name,
  D = 1,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("grdifrs", "grdifr", "grdifs"),
  max.iter = 10,
  min.resp = NULL,
  post.hoc = TRUE,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'est_irt'
grdif(
  x,
  score = NULL,
  group,
  focal.name,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("grdifrs", "grdifr", "grdifs"),
  max.iter = 10,
  min.resp = NULL,
  post.hoc = TRUE,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'est_item'
grdif(
  x,
  group,
  focal.name,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("grdifrs", "grdifr", "grdifs"),
  max.iter = 10,
  min.resp = NULL,
  post.hoc = TRUE,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'grdif'
print(x, digits = max(2L, getOption("digits") - 5L), ...)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.), or an object of class `est_irt`
  or `est_item` obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).
  All items must be dichotomous (`"1PLM"`, `"2PLM"`, `"3PLM"`, or
  `"DRM"`); the function stops when any item is polytomous.

  See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
  or [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- ...:

  Additional arguments passed to
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  when abilities are estimated internally, for example `tol` or `se`.
  The
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  arguments `max.iter` and `missing` cannot be passed this way, because
  the function uses arguments with these names for its own purposes.

- data:

  A matrix of examinees' scored responses (0 or 1) to the dichotomous
  items specified in the `x` argument. Rows represent examinees and
  columns represent items.

- score:

  A numeric vector of examinees' ability estimates (theta values), in
  the same row order as `data`. The estimates should be based on the
  item parameters in `x` (see **Details**). If `NULL`, abilities are
  estimated internally with
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  using `method`, `range`, `norm.prior`, `nquad`, `weights`, and
  `ncore`. When `purify = TRUE`, abilities are re-estimated internally
  at every purification iteration, so a supplied `score` is used only in
  the initial analysis. A missing value in `score` excludes the examinee
  from the analysis, also during purification, with a warning. The
  `est_item` method has no `score` argument; it uses the abilities
  stored in the object. Default is `NULL`.

- group:

  A numeric or character vector of examinees' group membership, with one
  element per row of `data`. It should contain the value of the
  reference group and the values in `focal.name`. Examinees whose value
  is not in `focal.name` form the reference group, and the function
  stops when they have more than one value. Examinees with `NA` in
  `group` are excluded from the analysis.

- focal.name:

  A numeric or character vector of the values in `group` that identify
  the focal groups, with at least one value. For example, given
  `group = c(0, 0, 1, 2, 2, 3, 3)`, where '0' indicates the reference
  group and '1', '2', and '3' indicate three focal groups, set
  `focal.name = c(1, 2, 3)`. The order of `focal.name` determines the
  order of the focal groups in the results (for example, the columns
  `n.foc1`, `n.foc2`, and `n.foc3`); the GRDIF statistics do not depend
  on it. All examinees whose group value is not in `focal.name` form the
  reference group (see **Details**). Each value must be found in
  `group`, and the values must be distinct.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) of the
  tests. An item is flagged by a statistic when its p-value is less than
  or equal to `alpha`. The comparison uses unrounded p-values, although
  the p-values reported in `dif_stat` are rounded to four decimal
  places. The same level is used in the post-hoc pairwise analyses,
  without an adjustment for multiple comparisons. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying the GRDIF statistic used for
  purification: "grdifrs" for \\GRDIF\_{RS}\\, "grdifr" for
  \\GRDIF\_{R}\\, or "grdifs" for \\GRDIF\_{S}\\. Used only when
  `purify = TRUE`. Default is `"grdifrs"`.

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
  `min.resp` can leave a group without any ability estimate in a later
  purification iteration, which stops the purification with a warning.
  If `NULL`, no minimum is applied. Default is `NULL`. See **Details**
  for more information.

- post.hoc:

  Logical. If `TRUE`, post-hoc pairwise RDIF analyses are performed, for
  every pair of groups, on the items flagged by each GRDIF statistic.
  The p-values of the post-hoc analyses are not adjusted for multiple
  comparisons. See **Details**. Default is `TRUE`.

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
  each response pattern. Default is `"ML"`.

- range:

  A numeric vector of length two specifying the lower and upper bounds
  of the ability scale used when abilities are estimated internally with
  the scoring methods `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. Default is
  `c(-5, 5)`, the default of the other DIF functions.

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

- digits:

  Not used; kept for consistency with the other print methods.

## Value

This function returns an object of class `"grdif"`, which is a list with
the following five components:

- no_purify:

  A list of sub-objects presenting the results of DIF analysis without a
  purification procedure. These include:

  dif_stat

  :   A data frame with one row per item and the columns `id`, `grdifr`,
      `grdifs`, `grdifrs`, `p.grdifr`, `p.grdifs`, `p.grdifrs`, `n.ref`,
      `n.foc1`, `n.foc2`, ..., and `n.total`: the item ID, the three
      GRDIF statistics, their p-values, the numbers of examinees in the
      reference group and in each focal group (in the order of
      `focal.name`) who responded to the item, and their sum. Statistics
      and p-values are rounded to four decimal places.

  moments

  :   A list of three data frames, `mu`, `var`, and `covariance`, with
      the null moments of the MRR and MSR of each group, computed
      analytically under the null hypothesis of no DIF (Lim et al.,
      2024, Equations 3 to 5). `mu` and `var` have the columns `id`,
      `mrr.ref`, `msr.ref`, `mrr.foc1`, `msr.foc1`, and so on, with the
      null means and null variances of the MRR and MSR of the reference
      group and of each focal group. `covariance` has the columns `id`,
      `ref`, `foc1`, and so on, with the null covariance between the MRR
      and MSR of each group.

  dif_item

  :   A list with the elements `grdifr`, `grdifs`, and `grdifrs`, each
      giving the positions (rows of `x`) of the items flagged by the
      corresponding statistic, or `NULL` if no item is flagged.

  score

  :   A numeric vector of the ability estimates used in the initial
      analysis: the supplied `score`, or the internal estimates when
      `score = NULL`.

  post.hoc

  :   A list with the elements `by.grdifr`, `by.grdifs`, and
      `by.grdifrs`. Each element is a data frame with the post-hoc
      results for the items flagged by the corresponding GRDIF
      statistic, with one row per item and pair of groups and the
      columns `id`, `group.pair`, `rdifr`, `z.rdifr`, `rdifs`,
      `z.rdifs`, `rdifrs`, `p.rdifr`, `p.rdifs`, `p.rdifrs`, `n.ref`,
      `n.foc`, and `n.total`, defined as in the `dif_stat` element of
      [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md);
      `n.ref` and `n.foc` refer to the first and second groups in
      `group.pair`. An element is `NULL` if no item is flagged by the
      statistic or if `post.hoc = FALSE`.

- purify:

  A logical value indicating whether purification was requested (the
  `purify` argument).

- with_purify:

  A list with the results of the purification procedure. All elements
  are `NULL` when `purify = FALSE`. The elements are:

  purify.by

  :   A character string indicating which GRDIF statistic was used for
      purification: "grdifr", "grdifs", or "grdifrs", corresponding to
      \\GRDIF\_{R}\\, \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\, respectively.

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

  moments

  :   A list of three data frames with the same columns as in
      `no_purify$moments` and an additional column `n.iter`, defined in
      the same way as in `dif_stat`.

  dif_item

  :   A numeric vector of the positions (rows of `x`) of the items
      flagged by the `purify.by` statistic, sorted in ascending order.
      It contains the items removed during purification and, if
      `max.iter` is reached or the purification stops because a group
      has no ability estimate left, the items flagged in the last
      iteration (iteration `n.iter`). It contains all items if every
      item is flagged. `NULL` if no item is flagged in the initial
      analysis.

  n.iter

  :   The number of purification iterations performed (0 if no item is
      flagged in the initial analysis). The iteration that removes the
      last item, or after which a group has no ability estimate left, is
      not counted, because the statistics are not recomputed.

  score

  :   A numeric vector of the ability estimates from the last iteration
      counted in `n.iter`. `NULL` if no item is flagged in the initial
      analysis, because purification is not carried out; the initial
      estimates are then in `no_purify$score`.

  post.hoc

  :   A data frame with the same columns as the elements of
      `no_purify$post.hoc` and an additional column `n_iter`. It
      contains the post-hoc results of each item removed during
      purification, computed with the ability estimates of the iteration
      in which the item was removed; `n_iter` is defined as `n.iter` in
      `dif_stat`. `NULL` if `post.hoc = FALSE` or if no item is flagged
      in the initial analysis.

  complete

  :   A logical value. `TRUE` if the procedure stopped because no
      remaining item was flagged (including the case in which no item is
      flagged in the initial analysis), and `FALSE` if it stopped
      because `max.iter` was reached, because every item was flagged, or
      because a group had no examinee with an ability estimate left.

- alpha:

  A numeric value indicating the significance level (\\\alpha\\) used
  for the tests.

- call:

  The matched function call.

## Details

The GRDIF framework (Lim et al., 2024) generalizes the RDIF framework
(Lim & Choe, 2023; Lim et al., 2022) to more than two groups. For each
item and group, the mean raw residual (MRR) and the mean squared
residual (MSR) are computed, where a residual is the observed item score
minus the model-expected item score evaluated at the examinee's ability
estimate. Under the null hypothesis that the item has no DIF between any
pair of groups, the MRRs and MSRs asymptotically follow normal
distributions whose means, variances, and covariances are computed
analytically from the model-predicted probabilities (Lim et al., 2024,
Equations 2 to 6). \\GRDIF\_{R}\\, \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\
are quadratic forms of the contrasts between the focal groups and the
reference group of the MRRs, of the MSRs, and of both, centered at their
null means, in the inverse of their null covariance matrices (Lim et
al., 2024, Equation 10). Under the null hypothesis, they asymptotically
follow \\\chi^{2}\\ distributions with G-1, G-1, and 2(G-1) degrees of
freedom, respectively, where G is the number of groups. The statistics
do not depend on which group is the reference group. For two groups,
\\GRDIF\_{R}\\ and \\GRDIF\_{S}\\ equal the squares of the standardized
\\RDIF\_{R}\\ and \\RDIF\_{S}\\ statistics of
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
\\GRDIF\_{RS}\\ equals \\RDIF\_{RS}\\, and the p-values are the same.

`grdif()` supports dichotomous items only, because the GRDIF framework
has been developed and evaluated for dichotomous items (Lim et al.,
2024). The function stops with an error when `x` includes any item with
three or more score categories or any item of the "GRM" or "GPCM" model.

To compute the GRDIF statistics, `grdif()` requires: (1) item parameter
estimates obtained from aggregate data (regardless of group membership),
(2) examinees' ability estimates (e.g., ML), and (3) examinees' item
response data. Note that the ability estimates must be based on the
aggregate-data item parameters. The item parameter estimates should be
provided in the `x` argument, the ability estimates in the `score`
argument, and the response data in the `data` argument. If ability
estimates are not provided (i.e., `score = NULL`), `grdif()` estimates
them with the scoring method specified in the `method` argument (e.g.,
`method = "ML"`). When `x` is an object of class `est_irt`, the item
parameter estimates, the response data, and the scaling constant `D` are
taken from the object. When `x` is an object of class `est_item`, the
ability values stored in the object are also used as `score`.

In rare cases, the covariance matrix of the contrasts that define
\\GRDIF\_{R}\\, \\GRDIF\_{S}\\, or \\GRDIF\_{RS}\\ is singular or nearly
singular (its reciprocal condition number is 1e-10 or less), for example
when all examinees of a group who responded to an item have the same
probability of a correct response, so that the MSR is a linear function
of the MRR. Then the statistic is computed with the Moore-Penrose
generalized inverse of the covariance matrix and compared with a
chi-square distribution whose degrees of freedom equal the rank of the
matrix (Moore, 1977), and one warning names the items. This handling is
not part of the original method (Lim et al., 2024) and is added in irtQ.
In most data the matrix is not singular, and the degrees of freedom are
G-1, G-1, and 2(G-1).

The `group` argument contains the group membership of each examinee as
numeric or character values, and its length must equal the number of
rows in `data`. The values listed in `focal.name` identify the focal
groups, and all other examinees form the reference group. The `group`
vector should therefore contain one value for the reference group in
addition to the values in `focal.name`, and the function stops when it
contains more than one such value. Examinees with `NA` in `group` are
excluded from the analysis.

As in the RDIF framework, an iterative purification procedure can be
applied (Lim et al., 2022, 2024). When `purify = TRUE`, the statistic
specified in `purify.by` drives the procedure (e.g.,
`purify.by = "grdifrs"`). At each iteration, the flagged item with the
smallest p-value of the `purify.by` statistic is removed (the p-values
are compared on the log scale), the abilities are re-estimated from the
remaining items with the scoring method specified in `method`, and the
GRDIF statistics of the remaining items are recomputed. A supplied
`score` is therefore used only in the initial analysis. The procedure
stops when no remaining item is flagged or when `max.iter` iterations
have been performed. If no item is flagged in the initial analysis, no
iteration is performed. If every item is flagged, the procedure stops
with a warning, because no item is left to estimate the abilities, and
`complete` is `FALSE`. The procedure also stops with a warning and
`complete = FALSE` when, after an item is removed, no examinee of the
reference group or of a focal group has an ability estimate left (see
`min.resp` below); the remaining items then keep the results and the
scores of the previous iteration, and the items flagged in that
iteration are added to the flagged items.

Scoring based on a small number of item responses can lead to large
standard errors, potentially reducing the accuracy of DIF detection in
the GRDIF framework. The `min.resp` argument excludes such examinees.
For example, if `min.resp = 5`, examinees who responded to fewer than
five items (but at least one) have all their responses treated as
missing (i.e., `NA`) before the abilities are estimated, in the initial
analysis and at every purification iteration, and a warning reports
their number in the initial analysis. As a result, their ability
estimates are missing, and they are not used in the computation of the
GRDIF statistics. When `score` is supplied, the scores of these
examinees are also treated as missing. During purification, the count is
based on the items that remain after the flagged items are removed, and
an examinee excluded once stays excluded. If `min.resp = NULL`, a score
is computed for any examinee with at least one valid item response.

Examinees whose ability estimate is `NA`, such as examinees without any
item response or examinees with `NA` in a supplied `score`, are excluded
from the computation of the GRDIF statistics and from the post-hoc
analyses. A warning reports the number of examinees with item responses
who are excluded because of a missing value in `score`; these examinees
stay excluded during purification. During purification, an examinee
whose remaining responses are all removed with the flagged items, or who
has fewer than `min.resp` responses to the remaining items, has no
ability estimate and is excluded, and one warning at the end reports the
number of such examinees. The sample sizes `n.ref`, `n.foc1`, `n.foc2`,
and so on count only the examinees who are used.

If no examinee of one of the groups responded to an item, the GRDIF
statistics and p-values of the item are `NaN`, and the item is not
flagged.

The GRDIF statistics test only whether DIF is present between at least
one pair of groups (Lim et al., 2024). When `post.hoc = TRUE`, the three
RDIF statistics of
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) are computed
for every pair of distinct values in `group`, for each item flagged by
\\GRDIF\_{R}\\, \\GRDIF\_{S}\\, or \\GRDIF\_{RS}\\. For example, with
four groups A, B, C, and D, the pairs A-B, A-C, A-D, B-C, B-D, and C-D
are compared. Each pairwise analysis uses the responses of the examinees
in the two groups, their ability estimates from the GRDIF analysis, and
the item parameters in `x`; the abilities are not re-estimated for the
pair. In a pair labeled "a & b" in the column `group.pair`, a is the
reference group and b is the focal group of the pair, so \\RDIF\_{R}\\
and \\RDIF\_{S}\\ are the residual means of group b minus those of group
a. The reference group of the GRDIF analysis is the reference group of
every pair that includes it (Lim et al., 2024, Equation 14), and for a
pair of focal groups, the group listed earlier in `focal.name` is the
reference group. The p-values are not adjusted for multiple comparisons;
to control the familywise error rate, compare them with a corrected
level, for example `alpha` divided by the number of pairs (a Bonferroni
correction). When `purify = TRUE`, the post-hoc results in `with_purify`
are computed for each item when it is removed during purification, with
the ability estimates of that iteration; the items added to the flagged
items when the purification ends before it is complete (`max.iter`
reached or no ability estimate left in a group) are not analyzed.

The print method of this class has an argument `digits` that is not
used; the statistics are printed rounded to three decimal places.

## Methods (by class)

- `grdif(default)`: Default method for a data frame `x` containing item
  metadata.

- `grdif(est_irt)`: Method for an object of class `est_irt` created by
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md). The
  response data and the scaling constant `D` are taken from the object.

- `grdif(est_item)`: Method for an object of class `est_item` created by
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).
  The response data, the ability estimates, and the scaling constant `D`
  are taken from the object.

## Methods (by generic)

- `print(grdif)`: Print method for an object of class `grdif` (given as
  `x`). It prints the flagged items and the GRDIF statistics of the
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

Lim, H., & Choe, E. M. (2023). Detecting differential item functioning
in CAT using IRT residual DIF approach. *Journal of Educational
Measurement, 60*(4), 626-650.
[doi:10.1111/jedm.12366](https://doi.org/10.1111/jedm.12366) .

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Lim, H., Zhu, D., Choe, E. M., & Han, K. T. (2024). Detecting
differential item functioning among multiple groups using IRT residual
DIF framework. *Journal of Educational Measurement, 61*(4), 656-681.
[doi:10.1111/jedm.12415](https://doi.org/10.1111/jedm.12415) .

Moore, D. S. (1977). Generalized inverses, Wald's method, and the
construction of chi-squared tests of fit. *Journal of the American
Statistical Association, 72*(357), 131-137.
[doi:10.1080/01621459.1977.10479921](https://doi.org/10.1080/01621459.1977.10479921)
.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627) .

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
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
# Load required library
library("dplyr")

## Uniform DIF detection for four groups (1 reference, 3 focal)
########################################################
# (1) Manipulate uniform DIF for all three focal groups
########################################################

# Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select 36 non-DIF items modeled under 3PLM
par_nstd <-
  bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
  dplyr::filter(.data$model == "3PLM") %>%
  dplyr::filter(dplyr::row_number() %in% 1:36) %>%
  dplyr::select(1:6)
par_nstd$id <- paste0("nondif", 1:36)

# Generate four new items on which uniform DIF will be imposed
difpar_ref <-
  shape_df(
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = rep(0.15, 4)),
    item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
  )

# Introduce DIF by shifting b-parameters differently for each focal group
difpar_foc1 <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(0.7, 0.7, 0, 0))
difpar_foc2 <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(0, 0, 0.7, 0.7))
difpar_foc3 <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(-0.4, -0.4, -0.5, -0.5))

# Combine the 4 DIF and 36 non-DIF items for all four groups
# Therefore, each of the first four items has uniform DIF in two of the three
# focal groups
par_ref <- rbind(difpar_ref, par_nstd)
par_foc1 <- rbind(difpar_foc1, par_nstd)
par_foc2 <- rbind(difpar_foc2, par_nstd)
par_foc3 <- rbind(difpar_foc3, par_nstd)

# Generate true abilities from different distributions
set.seed(128)
theta_ref <- rnorm(500, 0.0, 1.0)
theta_foc1 <- rnorm(500, -1.0, 1.0)
theta_foc2 <- rnorm(500, 1.0, 1.0)
theta_foc3 <- rnorm(500, 0.5, 1.0)

# Simulate response data for each group
resp_ref <- irtQ::simdat(par_ref, theta = theta_ref, D = 1)
resp_foc1 <- irtQ::simdat(par_foc1, theta = theta_foc1, D = 1)
resp_foc2 <- irtQ::simdat(par_foc2, theta = theta_foc2, D = 1)
resp_foc3 <- irtQ::simdat(par_foc3, theta = theta_foc3, D = 1)
data <- rbind(resp_ref, resp_foc1, resp_foc2, resp_foc3)

########################################################
# (2) Estimate item and ability parameters
#     using aggregated data
########################################################

# Estimate item parameters
est_mod <- irtQ::est_irt(data = data, D = 1, model = "3PLM")
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -48572.5390, Max-Change: 1.424379 EM iteration: 2, Loglike: -43508.4186, Max-Change: 0.387077 EM iteration: 3, Loglike: -43457.1832, Max-Change: 0.190572 EM iteration: 4, Loglike: -43441.4336, Max-Change: 0.109565 EM iteration: 5, Loglike: -43433.3523, Max-Change: 0.073514 EM iteration: 6, Loglike: -43428.7906, Max-Change: 0.05436 EM iteration: 7, Loglike: -43426.1172, Max-Change: 0.041286 EM iteration: 8, Loglike: -43424.5103, Max-Change: 0.031746 EM iteration: 9, Loglike: -43423.5226, Max-Change: 0.024547 EM iteration: 10, Loglike: -43422.9015, Max-Change: 0.019027 EM iteration: 11, Loglike: -43422.5011, Max-Change: 0.014762 EM iteration: 12, Loglike: -43422.2358, Max-Change: 0.011457 EM iteration: 13, Loglike: -43422.0545, Max-Change: 0.008891 EM iteration: 14, Loglike: -43421.9266, Max-Change: 0.006897 EM iteration: 15, Loglike: -43421.8334, Max-Change: 0.005349 EM iteration: 16, Loglike: -43421.7633, Max-Change: 0.004145 EM iteration: 17, Loglike: -43421.7091, Max-Change: 0.00321 EM iteration: 18, Loglike: -43421.6661, Max-Change: 0.002491 EM iteration: 19, Loglike: -43421.6313, Max-Change: 0.00196 EM iteration: 20, Loglike: -43421.6027, Max-Change: 0.001545 EM iteration: 21, Loglike: -43421.5789, Max-Change: 0.001268 EM iteration: 22, Loglike: -43421.5588, Max-Change: 0.00111 EM iteration: 23, Loglike: -43421.5418, Max-Change: 0.000978 EM iteration: 24, Loglike: -43421.5273, Max-Change: 0.000866 EM iteration: 25, Loglike: -43421.5149, Max-Change: 0.000778 EM iteration: 26, Loglike: -43421.5043, Max-Change: 0.000702 EM iteration: 27, Loglike: -43421.4951, Max-Change: 0.000635 EM iteration: 28, Loglike: -43421.4872, Max-Change: 0.000576 EM iteration: 29, Loglike: -43421.4804, Max-Change: 0.000524 EM iteration: 30, Loglike: -43421.4745, Max-Change: 0.000477 EM iteration: 31, Loglike: -43421.4694, Max-Change: 0.000435 EM iteration: 32, Loglike: -43421.4649, Max-Change: 0.000397 EM iteration: 33, Loglike: -43421.4610, Max-Change: 0.000364 EM iteration: 34, Loglike: -43421.4577, Max-Change: 0.000334 EM iteration: 35, Loglike: -43421.4547, Max-Change: 0.000308 EM iteration: 36, Loglike: -43421.4522, Max-Change: 0.000283 EM iteration: 37, Loglike: -43421.4499, Max-Change: 0.000261 EM iteration: 38, Loglike: -43421.4479, Max-Change: 0.00024 EM iteration: 39, Loglike: -43421.4462, Max-Change: 0.000221 EM iteration: 40, Loglike: -43421.4447, Max-Change: 0.000204 EM iteration: 41, Loglike: -43421.4433, Max-Change: 0.000188 EM iteration: 42, Loglike: -43421.4422, Max-Change: 0.000173 EM iteration: 43, Loglike: -43421.4411, Max-Change: 0.000159 EM iteration: 44, Loglike: -43421.4402, Max-Change: 0.000147 EM iteration: 45, Loglike: -43421.4394, Max-Change: 0.000135 EM iteration: 46, Loglike: -43421.4387, Max-Change: 0.000124 EM iteration: 47, Loglike: -43421.4380, Max-Change: 0.000114 EM iteration: 48, Loglike: -43421.4374, Max-Change: 0.000105 EM iteration: 49, Loglike: -43421.4369, Max-Change: 9.7e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 2.3 seconds. 
est_par <- est_mod$par.est

# Estimate ability parameters using MLE
score <- irtQ::est_score(x = est_par, data = data, method = "ML")$est.theta

########################################################
# (3) Conduct DIF analysis
########################################################

# Create a group membership vector:
# 0 = reference group; 1, 2, 3 = focal groups
group <- c(rep(0, 500), rep(1, 500), rep(2, 500), rep(3, 500))

# (a) Compute GRDIF statistics without purification,
#     and perform post-hoc pairwise comparisons for flagged items
dif_nopuri <- grdif(
  x = est_par, data = data, score = score, group = group,
  focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
  purify = FALSE, post.hoc = TRUE
)
print(dif_nopuri)
#> 
#> Call:
#> grdif.default(x = est_par, data = data, score = score, group = group, 
#>     focal.name = c(1, 2, 3), D = 1, alpha = 0.05, purify = FALSE, 
#>     post.hoc = TRUE)
#> 
#> DIF analysis using three GRDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by GRDIF(R): 
#>     1, 2, 3, 4, 20 
#>   - DIF Items identified by GRDIF(S): 
#>     1, 3, 4 
#>   - DIF Items identified by GRDIF(RS): 
#>     1, 2, 3, 4, 14, 20 
#>   - GRDIF Statistics: 
#> 
#>     id n.ref n.foc1 n.foc2 n.foc3 grdifr p.grdifr     grdifs p.grdifs    
#> 1   V1   500    500    500    500 14.402    0.002  **  8.270    0.041   *
#> 2   V2   500    500    500    500 53.260    0.000 ***  7.702    0.053   .
#> 3   V3   500    500    500    500 21.906    0.000 *** 19.841    0.000 ***
#> 4   V4   500    500    500    500 59.872    0.000 *** 44.745    0.000 ***
#> 5   V5   500    500    500    500  5.315    0.150      3.689    0.297    
#> 6   V6   500    500    500    500  1.844    0.605      0.898    0.826    
#> 7   V7   500    500    500    500  1.028    0.795      2.784    0.426    
#> 8   V8   500    500    500    500  1.116    0.773      4.242    0.236    
#> 9   V9   500    500    500    500  0.468    0.926      0.417    0.937    
#> 10 V10   500    500    500    500  1.685    0.640      5.254    0.154    
#> 11 V11   500    500    500    500  3.806    0.283      1.177    0.758    
#> 12 V12   500    500    500    500  0.812    0.847      0.801    0.849    
#> 13 V13   500    500    500    500  0.804    0.849      2.776    0.427    
#> 14 V14   500    500    500    500  5.050    0.168      6.102    0.107    
#> 15 V15   500    500    500    500  4.212    0.239      5.028    0.170    
#> 16 V16   500    500    500    500  1.296    0.730      0.481    0.923    
#> 17 V17   500    500    500    500  2.820    0.420      4.018    0.260    
#> 18 V18   500    500    500    500  1.882    0.597      2.038    0.565    
#> 19 V19   500    500    500    500  1.156    0.764      0.184    0.980    
#> 20 V20   500    500    500    500 14.987    0.002  **  0.709    0.871    
#> 21 V21   500    500    500    500  5.964    0.113      6.248    0.100   .
#> 22 V22   500    500    500    500  7.566    0.056   .  6.701    0.082   .
#> 23 V23   500    500    500    500  0.404    0.939      4.028    0.258    
#> 24 V24   500    500    500    500  6.357    0.096   .  5.456    0.141    
#> 25 V25   500    500    500    500  3.843    0.279      4.487    0.214    
#> 26 V26   500    500    500    500  0.824    0.844      1.327    0.723    
#> 27 V27   500    500    500    500  0.932    0.818      1.772    0.621    
#> 28 V28   500    500    500    500  7.060    0.070   .  1.733    0.630    
#> 29 V29   500    500    500    500  0.513    0.916      0.494    0.920    
#> 30 V30   500    500    500    500  0.726    0.867      1.934    0.586    
#> 31 V31   500    500    500    500  1.415    0.702      4.600    0.204    
#> 32 V32   500    500    500    500  1.946    0.584      1.613    0.656    
#> 33 V33   500    500    500    500  1.566    0.667      5.352    0.148    
#> 34 V34   500    500    500    500  4.252    0.236      7.012    0.072   .
#> 35 V35   500    500    500    500  1.112    0.774      2.256    0.521    
#> 36 V36   500    500    500    500  2.443    0.486      0.420    0.936    
#> 37 V37   500    500    500    500  0.467    0.926      0.571    0.903    
#> 38 V38   500    500    500    500  6.266    0.099   .  1.726    0.631    
#> 39 V39   500    500    500    500  2.092    0.554      1.432    0.698    
#> 40 V40   500    500    500    500  1.303    0.728      5.156    0.161    
#>    grdifrs p.grdifrs    
#> 1    17.24     0.008  **
#> 2    76.62     0.000 ***
#> 3    24.45     0.000 ***
#> 4    64.43     0.000 ***
#> 5     5.98     0.425    
#> 6     1.92     0.927    
#> 7     3.58     0.733    
#> 8     4.69     0.585    
#> 9     1.04     0.984    
#> 10    5.81     0.445    
#> 11    4.17     0.653    
#> 12    2.01     0.919    
#> 13    3.08     0.799    
#> 14   14.82     0.022   *
#> 15    6.00     0.423    
#> 16    3.62     0.728    
#> 17    6.05     0.418    
#> 18    2.63     0.853    
#> 19    1.34     0.969    
#> 20   15.96     0.014   *
#> 21   10.05     0.123    
#> 22   11.32     0.079   .
#> 23    4.93     0.552    
#> 24    9.70     0.138    
#> 25    9.64     0.141    
#> 26    3.91     0.689    
#> 27    3.53     0.740    
#> 28    9.01     0.173    
#> 29    1.16     0.979    
#> 30    5.88     0.437    
#> 31    6.72     0.347    
#> 32    4.61     0.595    
#> 33    8.24     0.221    
#> 34    8.53     0.202    
#> 35    2.37     0.883    
#> 36    3.25     0.776    
#> 37    1.36     0.968    
#> 38   10.34     0.111    
#> 39    4.77     0.573    
#> 40    9.14     0.166    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# Display post-hoc pairwise comparison results
print(dif_nopuri$no_purify$post.hoc)
#> $by.grdifr
#>     id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1   V1      0 & 1 -0.0591 -2.0484 -0.0149 -1.6572  6.0615  0.0405  0.0975
#> 2   V1      0 & 2 -0.0453 -1.6517 -0.0222  1.3734  3.1023  0.0986  0.1696
#> 3   V1      0 & 3  0.0356  1.2562 -0.0160 -0.3421  1.6008  0.2090  0.7323
#> 4   V1      1 & 2  0.0138  0.5009 -0.0073  2.8499  8.7763  0.6164  0.0044
#> 5   V1      1 & 3  0.0947  3.3321 -0.0011  1.2514 12.4413  0.0009  0.2108
#> 6   V1      2 & 3  0.0809  2.9989  0.0062 -1.6333  9.0714  0.0027  0.1024
#> 7   V2      0 & 1 -0.0980 -3.6243  0.0231  0.9144 17.6473  0.0003  0.3605
#> 8   V2      0 & 2 -0.0273 -1.1118 -0.0553 -0.2014  1.5799  0.2662  0.8404
#> 9   V2      0 & 3  0.0895  3.4605 -0.0408 -1.8841 13.5694  0.0005  0.0596
#> 10  V2      1 & 2  0.0708  2.8292 -0.0784 -1.0990 10.0051  0.0047  0.2718
#> 11  V2      1 & 3  0.1875  7.1227 -0.0639 -2.6860 68.9692  0.0000  0.0072
#> 12  V2      2 & 3  0.1168  4.9212  0.0145 -1.6721 24.5433  0.0000  0.0945
#> 13  V3      0 & 1 -0.0202 -0.6886  0.0154 -0.5843  1.1812  0.4911  0.5590
#> 14  V3      0 & 2 -0.0739 -2.7289 -0.0142  2.2621  7.5237  0.0064  0.0237
#> 15  V3      0 & 3  0.0485  1.7294 -0.0429 -2.5155  6.3308  0.0837  0.0119
#> 16  V3      1 & 2 -0.0538 -1.9297 -0.0296  2.9332  9.1012  0.0536  0.0034
#> 17  V3      1 & 3  0.0687  2.3856 -0.0583 -2.1799  7.6235  0.0171  0.0293
#> 18  V3      2 & 3  0.1225  4.6129 -0.0288 -4.4157 22.5836  0.0000  0.0000
#> 19  V4      0 & 1 -0.0259 -0.9319  0.0433  0.8278  1.4762  0.3514  0.4078
#> 20  V4      0 & 2 -0.1206 -4.9869  0.0036  5.4702 34.3231  0.0000  0.0000
#> 21  V4      0 & 3  0.0521  2.0285 -0.0319 -0.7679  4.2249  0.0425  0.4426
#> 22  V4      1 & 2 -0.0947 -3.7033 -0.0397  4.8454 31.5761  0.0002  0.0000
#> 23  V4      1 & 3  0.0780  2.8913 -0.0752 -1.5948  9.9383  0.0038  0.1108
#> 24  V4      2 & 3  0.1727  7.4400 -0.0355 -6.0648 57.0393  0.0000  0.0000
#> 25 V20      0 & 1 -0.0365 -1.5211 -0.0123  0.0022  2.5292  0.1282  0.9982
#> 26 V20      0 & 2  0.0130  0.5761 -0.0518 -0.6875  0.6457  0.5645  0.4918
#> 27 V20      0 & 3 -0.0649 -2.7241 -0.0160 -0.0213  7.5662  0.0064  0.9830
#> 28 V20      1 & 2  0.0495  2.2427 -0.0396 -0.6506  5.8726  0.0249  0.5153
#> 29 V20      1 & 3 -0.0284 -1.2205 -0.0037 -0.0224  1.5446  0.2223  0.9821
#> 30 V20      2 & 3 -0.0780 -3.5569  0.0358  0.6700 13.1469  0.0004  0.5028
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0483   500   500    1000
#> 2    0.2120   500   500    1000
#> 3    0.4492   500   500    1000
#> 4    0.0124   500   500    1000
#> 5    0.0020   500   500    1000
#> 6    0.0107   500   500    1000
#> 7    0.0001   500   500    1000
#> 8    0.4539   500   500    1000
#> 9    0.0011   500   500    1000
#> 10   0.0067   500   500    1000
#> 11   0.0000   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.5540   500   500    1000
#> 14   0.0232   500   500    1000
#> 15   0.0422   500   500    1000
#> 16   0.0106   500   500    1000
#> 17   0.0221   500   500    1000
#> 18   0.0000   500   500    1000
#> 19   0.4780   500   500    1000
#> 20   0.0000   500   500    1000
#> 21   0.1209   500   500    1000
#> 22   0.0000   500   500    1000
#> 23   0.0069   500   500    1000
#> 24   0.0000   500   500    1000
#> 25   0.2824   500   500    1000
#> 26   0.7241   500   500    1000
#> 27   0.0228   500   500    1000
#> 28   0.0531   500   500    1000
#> 29   0.4620   500   500    1000
#> 30   0.0014   500   500    1000
#> 
#> $by.grdifs
#>    id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1  V1      0 & 1 -0.0591 -2.0484 -0.0149 -1.6572  6.0615  0.0405  0.0975
#> 2  V1      0 & 2 -0.0453 -1.6517 -0.0222  1.3734  3.1023  0.0986  0.1696
#> 3  V1      0 & 3  0.0356  1.2562 -0.0160 -0.3421  1.6008  0.2090  0.7323
#> 4  V1      1 & 2  0.0138  0.5009 -0.0073  2.8499  8.7763  0.6164  0.0044
#> 5  V1      1 & 3  0.0947  3.3321 -0.0011  1.2514 12.4413  0.0009  0.2108
#> 6  V1      2 & 3  0.0809  2.9989  0.0062 -1.6333  9.0714  0.0027  0.1024
#> 7  V3      0 & 1 -0.0202 -0.6886  0.0154 -0.5843  1.1812  0.4911  0.5590
#> 8  V3      0 & 2 -0.0739 -2.7289 -0.0142  2.2621  7.5237  0.0064  0.0237
#> 9  V3      0 & 3  0.0485  1.7294 -0.0429 -2.5155  6.3308  0.0837  0.0119
#> 10 V3      1 & 2 -0.0538 -1.9297 -0.0296  2.9332  9.1012  0.0536  0.0034
#> 11 V3      1 & 3  0.0687  2.3856 -0.0583 -2.1799  7.6235  0.0171  0.0293
#> 12 V3      2 & 3  0.1225  4.6129 -0.0288 -4.4157 22.5836  0.0000  0.0000
#> 13 V4      0 & 1 -0.0259 -0.9319  0.0433  0.8278  1.4762  0.3514  0.4078
#> 14 V4      0 & 2 -0.1206 -4.9869  0.0036  5.4702 34.3231  0.0000  0.0000
#> 15 V4      0 & 3  0.0521  2.0285 -0.0319 -0.7679  4.2249  0.0425  0.4426
#> 16 V4      1 & 2 -0.0947 -3.7033 -0.0397  4.8454 31.5761  0.0002  0.0000
#> 17 V4      1 & 3  0.0780  2.8913 -0.0752 -1.5948  9.9383  0.0038  0.1108
#> 18 V4      2 & 3  0.1727  7.4400 -0.0355 -6.0648 57.0393  0.0000  0.0000
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0483   500   500    1000
#> 2    0.2120   500   500    1000
#> 3    0.4492   500   500    1000
#> 4    0.0124   500   500    1000
#> 5    0.0020   500   500    1000
#> 6    0.0107   500   500    1000
#> 7    0.5540   500   500    1000
#> 8    0.0232   500   500    1000
#> 9    0.0422   500   500    1000
#> 10   0.0106   500   500    1000
#> 11   0.0221   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.4780   500   500    1000
#> 14   0.0000   500   500    1000
#> 15   0.1209   500   500    1000
#> 16   0.0000   500   500    1000
#> 17   0.0069   500   500    1000
#> 18   0.0000   500   500    1000
#> 
#> $by.grdifrs
#>     id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1   V1      0 & 1 -0.0591 -2.0484 -0.0149 -1.6572  6.0615  0.0405  0.0975
#> 2   V1      0 & 2 -0.0453 -1.6517 -0.0222  1.3734  3.1023  0.0986  0.1696
#> 3   V1      0 & 3  0.0356  1.2562 -0.0160 -0.3421  1.6008  0.2090  0.7323
#> 4   V1      1 & 2  0.0138  0.5009 -0.0073  2.8499  8.7763  0.6164  0.0044
#> 5   V1      1 & 3  0.0947  3.3321 -0.0011  1.2514 12.4413  0.0009  0.2108
#> 6   V1      2 & 3  0.0809  2.9989  0.0062 -1.6333  9.0714  0.0027  0.1024
#> 7   V2      0 & 1 -0.0980 -3.6243  0.0231  0.9144 17.6473  0.0003  0.3605
#> 8   V2      0 & 2 -0.0273 -1.1118 -0.0553 -0.2014  1.5799  0.2662  0.8404
#> 9   V2      0 & 3  0.0895  3.4605 -0.0408 -1.8841 13.5694  0.0005  0.0596
#> 10  V2      1 & 2  0.0708  2.8292 -0.0784 -1.0990 10.0051  0.0047  0.2718
#> 11  V2      1 & 3  0.1875  7.1227 -0.0639 -2.6860 68.9692  0.0000  0.0072
#> 12  V2      2 & 3  0.1168  4.9212  0.0145 -1.6721 24.5433  0.0000  0.0945
#> 13  V3      0 & 1 -0.0202 -0.6886  0.0154 -0.5843  1.1812  0.4911  0.5590
#> 14  V3      0 & 2 -0.0739 -2.7289 -0.0142  2.2621  7.5237  0.0064  0.0237
#> 15  V3      0 & 3  0.0485  1.7294 -0.0429 -2.5155  6.3308  0.0837  0.0119
#> 16  V3      1 & 2 -0.0538 -1.9297 -0.0296  2.9332  9.1012  0.0536  0.0034
#> 17  V3      1 & 3  0.0687  2.3856 -0.0583 -2.1799  7.6235  0.0171  0.0293
#> 18  V3      2 & 3  0.1225  4.6129 -0.0288 -4.4157 22.5836  0.0000  0.0000
#> 19  V4      0 & 1 -0.0259 -0.9319  0.0433  0.8278  1.4762  0.3514  0.4078
#> 20  V4      0 & 2 -0.1206 -4.9869  0.0036  5.4702 34.3231  0.0000  0.0000
#> 21  V4      0 & 3  0.0521  2.0285 -0.0319 -0.7679  4.2249  0.0425  0.4426
#> 22  V4      1 & 2 -0.0947 -3.7033 -0.0397  4.8454 31.5761  0.0002  0.0000
#> 23  V4      1 & 3  0.0780  2.8913 -0.0752 -1.5948  9.9383  0.0038  0.1108
#> 24  V4      2 & 3  0.1727  7.4400 -0.0355 -6.0648 57.0393  0.0000  0.0000
#> 25 V14      0 & 1 -0.0448 -1.6362 -0.0101 -1.3402  3.3918  0.1018  0.1802
#> 26 V14      0 & 2 -0.0136 -0.5357 -0.0425  0.4584  0.3820  0.5921  0.6467
#> 27 V14      0 & 3 -0.0501 -1.8867 -0.0354 -1.6599  7.5622  0.0592  0.0969
#> 28 V14      1 & 2  0.0312  1.2205 -0.0324  1.7459  4.1363  0.2223  0.0808
#> 29 V14      1 & 3 -0.0053 -0.1978 -0.0253 -0.2443  0.0814  0.8432  0.8070
#> 30 V14      2 & 3 -0.0365 -1.4790  0.0071 -2.0749 11.3496  0.1391  0.0380
#> 31 V20      0 & 1 -0.0365 -1.5211 -0.0123  0.0022  2.5292  0.1282  0.9982
#> 32 V20      0 & 2  0.0130  0.5761 -0.0518 -0.6875  0.6457  0.5645  0.4918
#> 33 V20      0 & 3 -0.0649 -2.7241 -0.0160 -0.0213  7.5662  0.0064  0.9830
#> 34 V20      1 & 2  0.0495  2.2427 -0.0396 -0.6506  5.8726  0.0249  0.5153
#> 35 V20      1 & 3 -0.0284 -1.2205 -0.0037 -0.0224  1.5446  0.2223  0.9821
#> 36 V20      2 & 3 -0.0780 -3.5569  0.0358  0.6700 13.1469  0.0004  0.5028
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0483   500   500    1000
#> 2    0.2120   500   500    1000
#> 3    0.4492   500   500    1000
#> 4    0.0124   500   500    1000
#> 5    0.0020   500   500    1000
#> 6    0.0107   500   500    1000
#> 7    0.0001   500   500    1000
#> 8    0.4539   500   500    1000
#> 9    0.0011   500   500    1000
#> 10   0.0067   500   500    1000
#> 11   0.0000   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.5540   500   500    1000
#> 14   0.0232   500   500    1000
#> 15   0.0422   500   500    1000
#> 16   0.0106   500   500    1000
#> 17   0.0221   500   500    1000
#> 18   0.0000   500   500    1000
#> 19   0.4780   500   500    1000
#> 20   0.0000   500   500    1000
#> 21   0.1209   500   500    1000
#> 22   0.0000   500   500    1000
#> 23   0.0069   500   500    1000
#> 24   0.0000   500   500    1000
#> 25   0.1834   500   500    1000
#> 26   0.8261   500   500    1000
#> 27   0.0228   500   500    1000
#> 28   0.1264   500   500    1000
#> 29   0.9601   500   500    1000
#> 30   0.0034   500   500    1000
#> 31   0.2824   500   500    1000
#> 32   0.7241   500   500    1000
#> 33   0.0228   500   500    1000
#> 34   0.0531   500   500    1000
#> 35   0.4620   500   500    1000
#> 36   0.0014   500   500    1000
#> 

# (b) Compute GRDIF statistics with purification based on GRDIF_R,
#     including post-hoc comparisons (not run: it takes about 10 seconds
#     because the abilities of the examinees are re-estimated at every
#     iteration)
#
#   dif_puri_r <- grdif(
#     x = est_par, data = data, score = score, group = group,
#     focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
#     purify = TRUE, purify.by = "grdifr", post.hoc = TRUE
#   )
#   print(dif_puri_r)
#
#   # Display post-hoc results after purification
#   print(dif_puri_r$with_purify$post.hoc)
#
#   purify.by = "grdifs" or "grdifrs" can be used in the same way.
# }
```
