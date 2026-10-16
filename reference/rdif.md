# IRT Residual-Based Differential Item Functioning (RDIF) Detection Framework

This function computes three RDIF statistics for each item:
\\RDIF\_{R}\\, \\RDIF\_{S}\\, and \\RDIF\_{RS}\\ (Lim & Choe, 2023; Lim
et al., 2022). \\RDIF\_{R}\\ primarily captures differences in raw
residuals between two groups, which are typically associated with
uniform DIF. \\RDIF\_{S}\\ primarily captures differences in squared
residuals, which are typically associated with nonuniform DIF.
\\RDIF\_{RS}\\ jointly considers both types of differences and is
capable of detecting both uniform and nonuniform DIF.

## Usage

``` r
rdif(x, ...)

# Default S3 method
rdif(
  x,
  data,
  score = NULL,
  group,
  focal.name,
  item.skip = NULL,
  D = 1,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("rdifrs", "rdifr", "rdifs"),
  max.iter = 10,
  min.resp = NULL,
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
rdif(
  x,
  score = NULL,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("rdifrs", "rdifr", "rdifs"),
  max.iter = 10,
  min.resp = NULL,
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
rdif(
  x,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("rdifrs", "rdifr", "rdifs"),
  max.iter = 10,
  min.resp = NULL,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'rdif'
print(x, digits = max(2L, getOption("digits") - 5L), ...)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.), or an object of class `est_irt`
  or `est_item` obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

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

  A matrix of examinees' item responses, with rows for examinees and
  columns for items in the same order as the items in `x`. Responses to
  an item with K score categories must be integers from 0 to K - 1.
  Missing responses are coded as `NA` or, in functions that have a
  `missing` argument, as the value given in `missing`.

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

  A numeric or character vector indicating examinees' group membership.
  The length of the vector must match the number of rows in the response
  data matrix.

- focal.name:

  A single numeric or character value specifying the focal group. For
  instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating the
  focal group, set `focal.name = 1`.

- item.skip:

  A numeric vector of item positions (row numbers of `x`), or a logical
  vector with one value per item, identifying the items to exclude from
  the analysis. Skipped items still contribute to the ability estimates,
  but their results (statistics, p-values, and moments) are set to `NA`,
  and they are never flagged or removed during purification. If `NULL`,
  all items are analyzed. Default is `NULL`.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) of the
  tests. An item is flagged by a statistic when its p-value is less than
  or equal to `alpha`. The comparison uses unrounded p-values, although
  the p-values reported in `dif_stat` are rounded to four decimal
  places. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying the RDIF statistic used for
  purification: "rdifrs" for \\RDIF\_{RS}\\, "rdifr" for \\RDIF\_{R}\\,
  or "rdifs" for \\RDIF\_{S}\\. Used only when `purify = TRUE`. Default
  is `"rdifrs"`.

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
  of the ability scale. This is used for the following scoring methods:
  `"ML"`, `"WL"`, and `"MAP"`. Default is `c(-5, 5)`.

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

This function returns an object of class `"rdif"`, which is a list with
the following five components:

- no_purify:

  A list of sub-objects containing the results of DIF analysis without
  applying a purification procedure. The sub-objects include:

  dif_stat

  :   A data frame with one row per item and the columns `id`, `rdifr`,
      `z.rdifr`, `rdifs`, `z.rdifs`, `rdifrs`, `p.rdifr`, `p.rdifs`,
      `p.rdifrs`, `n.ref`, `n.foc`, and `n.total`: the item ID, the
      \\RDIF\_{R}\\ statistic and its standardized value, the
      \\RDIF\_{S}\\ statistic and its standardized value, the
      \\RDIF\_{RS}\\ statistic, the p-values of the three statistics
      (two-sided for \\RDIF\_{R}\\ and \\RDIF\_{S}\\), the numbers of
      examinees in the reference and focal groups who responded to the
      item, and their sum. Statistics and p-values are rounded to four
      decimal places. \\RDIF\_{RS}\\ has no standardized value because
      it is a \\\chi^{2}\\-based statistic.

  moments

  :   A data frame with the columns `id`, `mu.rdifr`, `sigma.rdifr`,
      `mu.rdifs`, `sigma.rdifs`, and `covariance`: the item ID, the null
      mean and standard deviation of \\RDIF\_{R}\\, the null mean and
      standard deviation of \\RDIF\_{S}\\, and the null covariance
      between \\RDIF\_{R}\\ and \\RDIF\_{S}\\. These moments are
      computed analytically under the null hypothesis of no DIF and are
      used to standardize the statistics.

  dif_item

  :   A list with the elements `rdifr`, `rdifs`, and `rdifrs`, each
      giving the positions (rows of `x`) of the items flagged by the
      corresponding statistic, or `NULL` if no item is flagged.

  score

  :   A numeric vector of the ability estimates used in the initial
      analysis: the supplied `score`, or the internal estimates when
      `score = NULL`.

- purify:

  A logical value indicating whether purification was requested (the
  `purify` argument).

- with_purify:

  A list with the results of the purification procedure. All elements
  are `NULL` when `purify = FALSE`. The elements are:

  purify.by

  :   A character string indicating the RDIF statistic used for
      purification. Possible values are "rdifr", "rdifs", and "rdifrs",
      corresponding to \\RDIF\_{R}\\, \\RDIF\_{S}\\, and \\RDIF\_{RS}\\,
      respectively.

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

  :   A data frame with the same columns as `no_purify$moments` and an
      additional column `n.iter`, defined in the same way as in
      `dif_stat`.

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

The RDIF framework (Lim & Choe, 2023; Lim et al., 2022) consists of
three IRT residual-based statistics: \\RDIF\_{R}\\, \\RDIF\_{S}\\, and
\\RDIF\_{RS}\\. For each item, a residual is the observed item score
minus the model-expected item score evaluated at the examinee's ability
estimate. \\RDIF\_{R}\\ is the difference between the focal and
reference groups in the mean raw residual, and \\RDIF\_{S}\\ is the
difference in the mean squared residual. Under the null hypothesis that
the item has no DIF, \\RDIF\_{R}\\ and \\RDIF\_{S}\\ asymptotically
follow normal distributions whose means and variances are computed
analytically from the model-predicted probabilities (Lim et al., 2022,
Equations 9, 10, 12, and 13). The null mean of \\RDIF\_{R}\\ is zero,
whereas the null mean of \\RDIF\_{S}\\ is generally not zero. Each
statistic is standardized by subtracting its null mean and dividing by
its null standard deviation (`z.rdifr` and `z.rdifs`), and a two-sided
p-value is obtained from the standard normal distribution.
\\RDIF\_{RS}\\ is the quadratic form of the vector of \\RDIF\_{R}\\ and
\\RDIF\_{S}\\, centered at their null means, in the inverse of their
null covariance matrix. Under the null hypothesis, \\RDIF\_{RS}\\
asymptotically follows a \\\chi^{2}\\ distribution with 2 degrees of
freedom (Lim et al., 2022, Equation 20).

In rare cases, the covariance matrix of \\RDIF\_{R}\\ and \\RDIF\_{S}\\
is singular or nearly singular (its reciprocal condition number is 1e-10
or less), for example when all examinees who responded to an item have
the same probability of a correct response, so that \\RDIF\_{S}\\ is a
linear function of \\RDIF\_{R}\\. Then \\RDIF\_{RS}\\ is computed with
the Moore-Penrose generalized inverse of the covariance matrix and
compared with a chi-square distribution whose degrees of freedom equal
the rank of the matrix (Moore, 1977), and a warning names the item. This
handling is not part of the original method (Lim et al., 2022) and is
added in irtQ. In most data the matrix is not singular, and
\\RDIF\_{RS}\\ has two degrees of freedom.

`rdif()` accepts both dichotomous and polytomous items. For a polytomous
item, the residual is the observed item score minus the model-expected
item score, and the null means, variances, and covariance are computed
from the model-predicted category probabilities in the same way, so the
three statistics assess net DIF (Jung & Lim, 2026; Lim, Malatesta, &
Lee, 2024). To evaluate global DIF in polytomous items, use
[`crdif()`](https://hwangQ.github.io/irtQ/reference/crdif.md).

To compute the RDIF statistics, the `rdif()` function requires: (1) item
parameter estimates obtained from aggregate data (regardless of group
membership), (2) examinees' ability estimates (e.g., ML), and (3)
examinees' item response data. Note that the ability estimates must be
based on the aggregate-data item parameters. The item parameter
estimates should be provided in the `x` argument, the ability estimates
in the `score` argument, and the response data in the `data` argument.
If ability estimates are not provided (i.e., `score = NULL`), `rdif()`
will estimate them automatically using the scoring method specified via
the `method` argument (e.g., `method = "ML"`). When `x` is an object of
class `est_irt`, the item parameter estimates, the response data, and
the scaling constant `D` are taken from the object. When `x` is an
object of class `est_item`, the ability values stored in the object are
also used as `score`.

The `group` argument should be a vector containing exactly two distinct
values (either numeric or character), representing the reference and
focal groups, and the function stops when it contains more than two. Its
length must match the number of rows in the response data, where each
element corresponds to an examinee. Once `group` is specified, a single
numeric or character value must be provided in the `focal.name` argument
to indicate which level in `group` represents the focal group. Use
[`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md) to compare
more than two groups. Examinees with `NA` in `group` are excluded from
both groups.

Examinees whose ability estimate is `NA`, such as examinees without any
item response or examinees with `NA` in a supplied `score`, are excluded
from the computation of the RDIF statistics. A warning reports the
number of examinees with item responses who are excluded because of a
missing value in `score`; these examinees stay excluded during
purification. During purification, an examinee whose remaining responses
are all removed with the flagged items, or who has fewer than `min.resp`
responses to the remaining items, has no ability estimate and is
excluded, and one warning at the end reports the number of such
examinees. The sample sizes `n.ref` and `n.foc` count only the examinees
who are used.

If no examinee of one of the two groups responded to an item, the
statistics and p-values of the item are `NaN`, and the item is not
flagged.

Similar to other DIF detection approaches, the RDIF framework supports
an iterative purification procedure (Lim et al., 2022). When
`purify = TRUE`, the statistic specified in `purify.by` drives the
procedure (e.g., `purify.by = "rdifrs"`). At each iteration, the flagged
item with the smallest p-value of the `purify.by` statistic is removed
(the p-values are compared on the log scale), the abilities are
re-estimated from the remaining items with the scoring method specified
in `method`, and the RDIF statistics of the remaining items are
recomputed. A supplied `score` is therefore used only in the initial
analysis. The procedure stops when no remaining item is flagged or when
`max.iter` iterations have been performed, where `max.iter` must be a
single whole number of at least 1. If no item is flagged in the initial
analysis, no iteration is performed. If every item is flagged, the
procedure stops with a warning, because no item is left to estimate the
abilities, and `complete` is `FALSE`. The procedure also stops with a
warning and `complete = FALSE` when, after an item is removed, no
examinee of the reference group or of the focal group has an ability
estimate left (see `min.resp` below); the remaining items then keep the
results of the previous iteration, and the items flagged in that
iteration are added to the flagged items.

Scoring based on a small number of item responses can lead to large
standard errors, potentially reducing the accuracy of DIF detection in
the RDIF framework. The `min.resp` argument excludes such examinees. For
example, if `min.resp = 5`, examinees who responded to fewer than five
items (but at least one) have all their responses treated as missing
(i.e., `NA`) before the abilities are estimated, in the initial analysis
and at every purification iteration. As a result, their ability
estimates are missing, and they are not used in the computation of the
RDIF statistics. When `score` is supplied, the scores of these examinees
are also treated as missing in the initial analysis. During
purification, the count is based on the items that remain after the
flagged items are removed, and an examinee excluded once stays excluded.
Because the count is based on the remaining items, a large `min.resp`
can leave a group without any ability estimate in a later iteration, and
the purification then stops with a warning. If `min.resp = NULL`, a
score is computed for any examinee with at least one valid item
response.

The print method of this class has an argument `digits` that is not
used; the statistics are printed rounded to three decimal places.

## Methods (by class)

- `rdif(default)`: Default method for a data frame `x` containing item
  metadata.

- `rdif(est_irt)`: Method for an object of class `est_irt` created by
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md). The
  response data and the scaling factor `D` are taken from the object.

- `rdif(est_item)`: Method for an object of class `est_item` created by
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).
  The response data, ability estimates, and the scaling factor `D` are
  taken from the object.

## Methods (by generic)

- `print(rdif)`: Print method for an object of class `rdif` (given as
  `x`). It prints the flagged items and the RDIF statistics of the
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

Jung, H., & Lim, H. (2026, April). Detecting global and net DIF in
polytomous items using RDIF. Paper presented at the annual meeting of
the National Council on Measurement in Education, Los Angeles, CA.

Lim, H., & Choe, E. M. (2023). Detecting differential item functioning
in CAT using IRT residual DIF approach. *Journal of Educational
Measurement, 60*(4), 626-650.
[doi:10.1111/jedm.12366](https://doi.org/10.1111/jedm.12366) .

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Lim, H., Malatesta, J., & Lee, Y. (2024, July). Advancing polytomous DIF
detection with the residual DIF framework. Paper presented at the annual
International Meeting of the Psychometric Society, Prague, Czech
Republic.

Moore, D. S. (1977). Generalized inverses, Wald's method, and the
construction of chi-squared tests of fit. *Journal of the American
Statistical Association, 72*(357), 131-137.
[doi:10.1080/01621459.1977.10479921](https://doi.org/10.1080/01621459.1977.10479921)
.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627) .

## See also

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

## Uniform DIF detection
###############################################
# (1) Generate data with known uniform DIF
###############################################

# Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select 36 non-DIF items using the 3PLM model
par_nstd <-
  bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
  dplyr::filter(.data$model == "3PLM") %>%
  dplyr::filter(dplyr::row_number() %in% 1:36) %>%
  dplyr::select(1:6)
par_nstd$id <- paste0("nondif", 1:36)

# Generate 4 new DIF items for the reference group
difpar_ref <-
  shape_df(
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = rep(0.15, 4)),
    item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
  )

# Add uniform DIF by shifting the b-parameters for the focal group
difpar_foc <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))

# Combine the DIF and non-DIF items for both reference and focal groups
# Therefore, the first 4 items exhibit uniform DIF
par_ref <- rbind(difpar_ref, par_nstd)
par_foc <- rbind(difpar_foc, par_nstd)

# Generate true ability values
set.seed(123)
theta_ref <- rnorm(500, 0.0, 1.0)
theta_foc <- rnorm(500, 0.0, 1.0)

# Simulate response data
resp_ref <- simdat(par_ref, theta = theta_ref, D = 1)
resp_foc <- simdat(par_foc, theta = theta_foc, D = 1)
data <- rbind(resp_ref, resp_foc)

###############################################
# (2) Estimate item and ability parameters
#     from the combined response data
###############################################

# Estimate item parameters
est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -25472.9403, Max-Change: 1.888841 EM iteration: 2, Loglike: -22654.4441, Max-Change: 0.356802 EM iteration: 3, Loglike: -22635.6173, Max-Change: 0.157888 EM iteration: 4, Loglike: -22633.7926, Max-Change: 0.074945 EM iteration: 5, Loglike: -22633.4312, Max-Change: 0.039693 EM iteration: 6, Loglike: -22633.3151, Max-Change: 0.022643 EM iteration: 7, Loglike: -22633.2668, Max-Change: 0.013669 EM iteration: 8, Loglike: -22633.2439, Max-Change: 0.008635 EM iteration: 9, Loglike: -22633.2322, Max-Change: 0.005655 EM iteration: 10, Loglike: -22633.2259, Max-Change: 0.00381 EM iteration: 11, Loglike: -22633.2225, Max-Change: 0.002624 EM iteration: 12, Loglike: -22633.2205, Max-Change: 0.001837 EM iteration: 13, Loglike: -22633.2194, Max-Change: 0.001303 EM iteration: 14, Loglike: -22633.2188, Max-Change: 0.000932 EM iteration: 15, Loglike: -22633.2184, Max-Change: 0.000671 EM iteration: 16, Loglike: -22633.2182, Max-Change: 0.000486 EM iteration: 17, Loglike: -22633.2180, Max-Change: 0.000353 EM iteration: 18, Loglike: -22633.2180, Max-Change: 0.000258 EM iteration: 19, Loglike: -22633.2179, Max-Change: 0.000188 EM iteration: 20, Loglike: -22633.2179, Max-Change: 0.000144 EM iteration: 21, Loglike: -22633.2179, Max-Change: 0.000116 EM iteration: 22, Loglike: -22633.2179, Max-Change: 9.7e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 1.58 seconds. 
est_par <- est_mod$par.est

# Estimate ability parameters using ML
score <- est_score(x = est_par, data = data, method = "ML")$est.theta

###############################################
# (3) Perform DIF analysis
###############################################

# Define group membership: 1 = focal group
group <- c(rep(0, 500), rep(1, 500))

# (a)-1 Compute RDIF statistics with provided ability scores
#       (no purification)
dif_nopuri_1 <- rdif(
  x = est_par, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05
)
print(dif_nopuri_1)
#> 
#> Call:
#> rdif.default(x = est_par, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05)
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 3, 4, 11, 14 
#>   - DIF Items identified by RDIF(S): 
#>     4, 5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs    
#> 1   V1   500   500 -0.112   0.000 ***  0.009   0.170     14.359    0.001 ***
#> 2   V2   500   500 -0.162   0.000 *** -0.006   0.730     33.029    0.000 ***
#> 3   V3   500   500 -0.069   0.018   *  0.013   0.189      5.622    0.060   .
#> 4   V4   500   500 -0.161   0.000 ***  0.038   0.001 *** 38.669    0.000 ***
#> 5   V5   500   500 -0.009   0.778      0.020   0.003  **  9.902    0.007  **
#> 6   V6   500   500  0.011   0.593     -0.001   0.463      0.546    0.761    
#> 7   V7   500   500  0.016   0.596      0.005   0.427      0.858    0.651    
#> 8   V8   500   500  0.025   0.359      0.002   0.954      1.223    0.542    
#> 9   V9   500   500 -0.014   0.618      0.008   0.439      0.669    0.716    
#> 10 V10   500   500  0.027   0.304      0.010   0.416      1.173    0.556    
#> 11 V11   500   500  0.064   0.029   * -0.008   0.640      7.608    0.022   *
#> 12 V12   500   500  0.000   0.988      0.001   0.662      0.246    0.884    
#> 13 V13   500   500 -0.012   0.683     -0.003   0.739      0.352    0.839    
#> 14 V14   500   500  0.064   0.019   * -0.007   0.542      5.722    0.057   .
#> 15 V15   500   500  0.017   0.540     -0.006   0.430      0.651    0.722    
#> 16 V16   500   500  0.005   0.873     -0.016   0.280      2.498    0.287    
#> 17 V17   500   500  0.008   0.767     -0.008   0.663      1.020    0.601    
#> 18 V18   500   500 -0.021   0.438      0.007   0.551      0.749    0.688    
#> 19 V19   500   500  0.049   0.069   . -0.021   0.047   *  5.436    0.066   .
#> 20 V20   500   500 -0.006   0.789     -0.019   0.140      2.217    0.330    
#> 21 V21   500   500  0.015   0.578     -0.002   0.767      0.327    0.849    
#> 22 V22   500   500  0.025   0.407     -0.005   0.569      1.663    0.435    
#> 23 V23   500   500  0.014   0.458     -0.005   0.205      1.614    0.446    
#> 24 V24   500   500 -0.001   0.963      0.015   0.530      0.979    0.613    
#> 25 V25   500   500  0.018   0.386     -0.005   0.351      0.924    0.630    
#> 26 V26   500   500  0.022   0.423      0.002   0.935      1.550    0.461    
#> 27 V27   500   500  0.021   0.440     -0.002   0.779      0.597    0.742    
#> 28 V28   500   500  0.022   0.462     -0.004   0.854      1.257    0.533    
#> 29 V29   500   500  0.005   0.838      0.006   0.953      0.447    0.800    
#> 30 V30   500   500  0.027   0.196     -0.007   0.277      1.673    0.433    
#> 31 V31   500   500  0.010   0.714     -0.003   0.774      0.184    0.912    
#> 32 V32   500   500  0.024   0.336     -0.008   0.298      1.664    0.435    
#> 33 V33   500   500 -0.028   0.193      0.018   0.390      1.850    0.396    
#> 34 V34   500   500 -0.013   0.673     -0.009   0.325      1.319    0.517    
#> 35 V35   500   500  0.009   0.770     -0.017   0.084   .  3.969    0.137    
#> 36 V36   500   500 -0.021   0.353      0.022   0.224      1.477    0.478    
#> 37 V37   500   500  0.005   0.802      0.003   0.661      0.204    0.903    
#> 38 V38   500   500 -0.054   0.064   . -0.009   0.446      4.139    0.126    
#> 39 V39   500   500  0.022   0.426      0.006   0.799      0.974    0.614    
#> 40 V40   500   500 -0.016   0.578     -0.008   0.625      0.369    0.832    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (a)-2 Compute RDIF statistics without providing ability scores
#       (no purification)
dif_nopuri_2 <- rdif(
  x = est_par, data = data, score = NULL,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  method = "ML"
)
print(dif_nopuri_2)
#> 
#> Call:
#> rdif.default(x = est_par, data = data, score = NULL, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, method = "ML")
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 3, 4, 11, 14 
#>   - DIF Items identified by RDIF(S): 
#>     4, 5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs    
#> 1   V1   500   500 -0.112   0.000 ***  0.009   0.170     14.359    0.001 ***
#> 2   V2   500   500 -0.162   0.000 *** -0.006   0.730     33.029    0.000 ***
#> 3   V3   500   500 -0.069   0.018   *  0.013   0.189      5.622    0.060   .
#> 4   V4   500   500 -0.161   0.000 ***  0.038   0.001 *** 38.669    0.000 ***
#> 5   V5   500   500 -0.009   0.778      0.020   0.003  **  9.902    0.007  **
#> 6   V6   500   500  0.011   0.593     -0.001   0.463      0.546    0.761    
#> 7   V7   500   500  0.016   0.596      0.005   0.427      0.858    0.651    
#> 8   V8   500   500  0.025   0.359      0.002   0.954      1.223    0.542    
#> 9   V9   500   500 -0.014   0.618      0.008   0.439      0.669    0.716    
#> 10 V10   500   500  0.027   0.304      0.010   0.416      1.173    0.556    
#> 11 V11   500   500  0.064   0.029   * -0.008   0.640      7.608    0.022   *
#> 12 V12   500   500  0.000   0.988      0.001   0.662      0.246    0.884    
#> 13 V13   500   500 -0.012   0.683     -0.003   0.739      0.352    0.839    
#> 14 V14   500   500  0.064   0.019   * -0.007   0.542      5.722    0.057   .
#> 15 V15   500   500  0.017   0.540     -0.006   0.430      0.651    0.722    
#> 16 V16   500   500  0.005   0.873     -0.016   0.280      2.498    0.287    
#> 17 V17   500   500  0.008   0.767     -0.008   0.663      1.020    0.601    
#> 18 V18   500   500 -0.021   0.438      0.007   0.551      0.749    0.688    
#> 19 V19   500   500  0.049   0.069   . -0.021   0.047   *  5.436    0.066   .
#> 20 V20   500   500 -0.006   0.789     -0.019   0.140      2.217    0.330    
#> 21 V21   500   500  0.015   0.578     -0.002   0.767      0.327    0.849    
#> 22 V22   500   500  0.025   0.407     -0.005   0.569      1.663    0.435    
#> 23 V23   500   500  0.014   0.458     -0.005   0.205      1.614    0.446    
#> 24 V24   500   500 -0.001   0.963      0.015   0.530      0.979    0.613    
#> 25 V25   500   500  0.018   0.386     -0.005   0.351      0.924    0.630    
#> 26 V26   500   500  0.022   0.423      0.002   0.935      1.550    0.461    
#> 27 V27   500   500  0.021   0.440     -0.002   0.779      0.597    0.742    
#> 28 V28   500   500  0.022   0.462     -0.004   0.854      1.257    0.533    
#> 29 V29   500   500  0.005   0.838      0.006   0.953      0.447    0.800    
#> 30 V30   500   500  0.027   0.196     -0.007   0.277      1.673    0.433    
#> 31 V31   500   500  0.010   0.714     -0.003   0.774      0.184    0.912    
#> 32 V32   500   500  0.024   0.336     -0.008   0.298      1.664    0.435    
#> 33 V33   500   500 -0.028   0.193      0.018   0.390      1.850    0.396    
#> 34 V34   500   500 -0.013   0.673     -0.009   0.325      1.319    0.517    
#> 35 V35   500   500  0.009   0.770     -0.017   0.084   .  3.969    0.137    
#> 36 V36   500   500 -0.021   0.353      0.022   0.224      1.477    0.478    
#> 37 V37   500   500  0.005   0.802      0.003   0.661      0.204    0.903    
#> 38 V38   500   500 -0.054   0.064   . -0.009   0.446      4.139    0.126    
#> 39 V39   500   500  0.022   0.426      0.006   0.799      0.974    0.614    
#> 40 V40   500   500 -0.016   0.578     -0.008   0.625      0.369    0.832    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (b) Compute RDIF statistics with purification based on RDIF(RS)
dif_puri_rs <- rdif(
  x = est_par, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "rdifrs"
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 Iteration: 3 Iteration: 4 Iteration: 5 
#> Purification is finished. 
print(dif_puri_rs)
#> 
#> Call:
#> rdif.default(x = est_par, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, purify = TRUE, purify.by = "rdifrs")
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 3, 4, 11, 14 
#>   - DIF Items identified by RDIF(S): 
#>     4, 5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs    
#> 1   V1   500   500 -0.112   0.000 ***  0.009   0.170     14.359    0.001 ***
#> 2   V2   500   500 -0.162   0.000 *** -0.006   0.730     33.029    0.000 ***
#> 3   V3   500   500 -0.069   0.018   *  0.013   0.189      5.622    0.060   .
#> 4   V4   500   500 -0.161   0.000 ***  0.038   0.001 *** 38.669    0.000 ***
#> 5   V5   500   500 -0.009   0.778      0.020   0.003  **  9.902    0.007  **
#> 6   V6   500   500  0.011   0.593     -0.001   0.463      0.546    0.761    
#> 7   V7   500   500  0.016   0.596      0.005   0.427      0.858    0.651    
#> 8   V8   500   500  0.025   0.359      0.002   0.954      1.223    0.542    
#> 9   V9   500   500 -0.014   0.618      0.008   0.439      0.669    0.716    
#> 10 V10   500   500  0.027   0.304      0.010   0.416      1.173    0.556    
#> 11 V11   500   500  0.064   0.029   * -0.008   0.640      7.608    0.022   *
#> 12 V12   500   500  0.000   0.988      0.001   0.662      0.246    0.884    
#> 13 V13   500   500 -0.012   0.683     -0.003   0.739      0.352    0.839    
#> 14 V14   500   500  0.064   0.019   * -0.007   0.542      5.722    0.057   .
#> 15 V15   500   500  0.017   0.540     -0.006   0.430      0.651    0.722    
#> 16 V16   500   500  0.005   0.873     -0.016   0.280      2.498    0.287    
#> 17 V17   500   500  0.008   0.767     -0.008   0.663      1.020    0.601    
#> 18 V18   500   500 -0.021   0.438      0.007   0.551      0.749    0.688    
#> 19 V19   500   500  0.049   0.069   . -0.021   0.047   *  5.436    0.066   .
#> 20 V20   500   500 -0.006   0.789     -0.019   0.140      2.217    0.330    
#> 21 V21   500   500  0.015   0.578     -0.002   0.767      0.327    0.849    
#> 22 V22   500   500  0.025   0.407     -0.005   0.569      1.663    0.435    
#> 23 V23   500   500  0.014   0.458     -0.005   0.205      1.614    0.446    
#> 24 V24   500   500 -0.001   0.963      0.015   0.530      0.979    0.613    
#> 25 V25   500   500  0.018   0.386     -0.005   0.351      0.924    0.630    
#> 26 V26   500   500  0.022   0.423      0.002   0.935      1.550    0.461    
#> 27 V27   500   500  0.021   0.440     -0.002   0.779      0.597    0.742    
#> 28 V28   500   500  0.022   0.462     -0.004   0.854      1.257    0.533    
#> 29 V29   500   500  0.005   0.838      0.006   0.953      0.447    0.800    
#> 30 V30   500   500  0.027   0.196     -0.007   0.277      1.673    0.433    
#> 31 V31   500   500  0.010   0.714     -0.003   0.774      0.184    0.912    
#> 32 V32   500   500  0.024   0.336     -0.008   0.298      1.664    0.435    
#> 33 V33   500   500 -0.028   0.193      0.018   0.390      1.850    0.396    
#> 34 V34   500   500 -0.013   0.673     -0.009   0.325      1.319    0.517    
#> 35 V35   500   500  0.009   0.770     -0.017   0.084   .  3.969    0.137    
#> 36 V36   500   500 -0.021   0.353      0.022   0.224      1.477    0.478    
#> 37 V37   500   500  0.005   0.802      0.003   0.661      0.204    0.903    
#> 38 V38   500   500 -0.054   0.064   . -0.009   0.446      4.139    0.126    
#> 39 V39   500   500  0.022   0.426      0.006   0.799      0.974    0.614    
#> 40 V40   500   500 -0.016   0.578     -0.008   0.625      0.369    0.832    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 5
#>   - RDIF statistic used for purification: RDIF(RS)
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 3, 4, 5 
#>   - RDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs
#> 1   V1      2   500   500 -0.118   0.000 ***  0.010   0.115     16.279    0.000
#> 2   V2      1   500   500 -0.168   0.000 *** -0.006   0.769     35.752    0.000
#> 3   V3      4   500   500 -0.079   0.007  **  0.013   0.103      7.429    0.024
#> 4   V4      0   500   500 -0.161   0.000 ***  0.038   0.001 *** 38.669    0.000
#> 5   V5      3   500   500 -0.016   0.592      0.021   0.003  ** 10.881    0.004
#> 6   V6      5   500   500  0.001   0.956      0.000   0.867      0.032    0.984
#> 7   V7      5   500   500  0.002   0.932      0.006   0.301      1.070    0.586
#> 8   V8      5   500   500  0.014   0.609      0.001   0.832      0.681    0.712
#> 9   V9      5   500   500 -0.027   0.341      0.009   0.239      1.758    0.415
#> 10 V10      5   500   500  0.008   0.768      0.009   0.430      0.635    0.728
#> 11 V11      5   500   500  0.051   0.081   . -0.006   0.760      4.670    0.097
#> 12 V12      5   500   500 -0.010   0.741      0.001   0.654      0.549    0.760
#> 13 V13      5   500   500 -0.022   0.462     -0.003   0.983      0.575    0.750
#> 14 V14      5   500   500  0.048   0.079   . -0.006   0.859      3.081    0.214
#> 15 V15      5   500   500  0.006   0.836     -0.006   0.692      0.159    0.924
#> 16 V16      5   500   500 -0.007   0.794     -0.016   0.276      1.541    0.463
#> 17 V17      5   500   500 -0.005   0.864     -0.006   0.804      0.062    0.970
#> 18 V18      5   500   500 -0.036   0.190      0.010   0.240      2.396    0.302
#> 19 V19      5   500   500  0.034   0.204     -0.022   0.091   .  3.425    0.180
#> 20 V20      5   500   500 -0.027   0.271     -0.020   0.206      2.673    0.263
#> 21 V21      5   500   500  0.000   0.994     -0.006   0.825      0.052    0.974
#> 22 V22      5   500   500  0.013   0.662     -0.007   0.608      0.750    0.687
#> 23 V23      5   500   500  0.004   0.833     -0.004   0.455      0.605    0.739
#> 24 V24      5   500   500 -0.007   0.700      0.014   0.349      1.241    0.538
#> 25 V25      5   500   500  0.010   0.635     -0.005   0.549      0.360    0.835
#> 26 V26      5   500   500  0.014   0.614      0.003   0.785      1.552    0.460
#> 27 V27      5   500   500  0.008   0.773     -0.002   0.984      0.106    0.949
#> 28 V28      5   500   500  0.010   0.733     -0.005   0.824      0.438    0.803
#> 29 V29      5   500   500 -0.001   0.980      0.005   0.807      0.317    0.853
#> 30 V30      5   500   500  0.021   0.322     -0.005   0.569      1.298    0.523
#> 31 V31      5   500   500 -0.004   0.884     -0.003   0.976      0.021    0.989
#> 32 V32      5   500   500  0.008   0.757     -0.007   0.629      0.280    0.869
#> 33 V33      5   500   500 -0.035   0.100   .  0.017   0.252      2.862    0.239
#> 34 V34      5   500   500 -0.024   0.412     -0.011   0.349      1.851    0.396
#> 35 V35      5   500   500 -0.002   0.945     -0.018   0.089   .  3.226    0.199
#> 36 V36      5   500   500 -0.030   0.180      0.021   0.110      2.602    0.272
#> 37 V37      5   500   500 -0.004   0.865      0.002   0.981      0.076    0.963
#> 38 V38      5   500   500 -0.068   0.020   * -0.007   0.737      5.559    0.062
#> 39 V39      5   500   500  0.008   0.780      0.005   0.527      0.696    0.706
#> 40 V40      5   500   500 -0.029   0.314     -0.009   0.616      1.016    0.602
#>       
#> 1  ***
#> 2  ***
#> 3    *
#> 4  ***
#> 5   **
#> 6     
#> 7     
#> 8     
#> 9     
#> 10    
#> 11   .
#> 12    
#> 13    
#> 14    
#> 15    
#> 16    
#> 17    
#> 18    
#> 19    
#> 20    
#> 21    
#> 22    
#> 23    
#> 24    
#> 25    
#> 26    
#> 27    
#> 28    
#> 29    
#> 30    
#> 31    
#> 32    
#> 33    
#> 34    
#> 35    
#> 36    
#> 37    
#> 38   .
#> 39    
#> 40    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 

# purify.by = "rdifr" or "rdifs" can be used in the same way.

# }
```
