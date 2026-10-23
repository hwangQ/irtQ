# Residual-Based DIF Detection Framework Using Categorical Residuals (RDIF-CR)

This function computes three statistics of the residual-based DIF
detection framework using categorical residuals (RDIF-CR),
\\RDIF\_{R}-CR\\, \\RDIF\_{S}-CR\\, and \\RDIF\_{RS}-CR\\, for detecting
global differential item functioning (DIF), particularly in polytomously
scored items. The RDIF-CR framework evaluates DIF by comparing
categorical residual vectors, which are calculated as the difference
between a one-hot encoded response vector (with 1 for the selected
category and 0 for all others) and the vector of IRT model-predicted
probabilities across all score categories. This approach enables
fine-grained detection of global DIF patterns at the category level.

## Usage

``` r
crdif(x, ...)

# Default S3 method
crdif(
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
  purify.by = c("crdifrs", "crdifr", "crdifs"),
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
crdif(
  x,
  score = NULL,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("crdifrs", "crdifr", "crdifs"),
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
crdif(
  x,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("crdifrs", "crdifr", "crdifs"),
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

# S3 method for class 'crdif'
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

  A character string specifying the RDIF-CR statistic used for
  purification: `"crdifrs"` for \\RDIF\_{RS}-CR\\, `"crdifr"` for
  \\RDIF\_{R}-CR\\, or `"crdifs"` for \\RDIF\_{S}-CR\\. Used only when
  `purify = TRUE`. Default is `"crdifrs"`.

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

This function returns an object of class `"crdif"`, which is a list with
the following five components:

- no_purify:

  A list of sub-objects containing the results of DIF analysis without
  applying a purification procedure. The sub-objects include:

  dif_stat

  :   A data frame with one row per item and the columns `id`, `crdifr`,
      `df.crdifr`, `crdifs`, `df.crdifs`, `crdifrs`, `df.crdifrs`,
      `p.crdifr`, `p.crdifs`, `p.crdifrs`, `n.ref`, `n.foc`, and
      `n.total`: the item ID, each RDIF-CR statistic followed by the
      degrees of freedom used, the p-values of the three statistics, the
      numbers of examinees in the reference and focal groups who
      responded to the item, and their sum. Statistics and p-values are
      rounded to four decimal places.

  moments

  :   A list with the elements `mu.crdifr`, `mu.crdifs`, `mu.crdifrs`,
      `cov.crdifr`, `cov.crdifs`, and `cov.crdifrs`. Each element is a
      list named by item ID that holds, for each item, the null mean
      vector or the null covariance matrix of the category-level
      statistics. For an item with K score categories (K of 3 or more),
      the vectors have lengths K, K, and 2K, and the matrices have
      dimensions K by K, K by K, and 2K by 2K. For a dichotomous item,
      only score category 1 is used, so the lengths are 1, 1, and 2.

  dif_item

  :   A list with the elements `crdifr`, `crdifs`, and `crdifrs`, each
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

  :   A character string indicating the RDIF-CR statistic used for
      purification. Possible values are "crdifr", "crdifs", or
      "crdifrs".

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

  :   A list with the same structure as `no_purify$moments`. For an item
      removed during purification, it holds the moments from the
      iteration in which the item was removed; for the other items, it
      holds the moments from the last iteration. If no item is flagged
      in the initial analysis, this is `no_purify$moments`.

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

According to Penfield (2010), differential item functioning (DIF) in
polytomously scored items can be conceptualized in two forms: global DIF
and net DIF. Global DIF refers to differences between groups in the
conditional probabilities of responding in specific score categories,
thus offering a fine-grained view of DIF at the category level. In
contrast, net DIF summarizes these differences into a single value
representing the overall impact of DIF on the item's expected score.

The RDIF framework using categorical residuals (RDIF-CR), implemented in
`crdif()`, extends the original residual-based DIF framework proposed by
Lim et al. (2022) to detect global DIF in polytomous items (Jung & Lim,
2026; Lim, Malatesta, & Lee, 2024). This framework includes three
statistics: \\RDIF\_{R}-CR\\, \\RDIF\_{S}-CR\\, and \\RDIF\_{RS}-CR\\,
each designed to capture different aspects of group-level differences in
categorical response patterns.

To illustrate how the RDIF-CR framework operates, consider an item with
five ordered score categories (\\k \in \\0,1,2,3,4\\\\). Suppose an
examinee with latent ability \\\theta\\ responds with category 2. The
one-hot encoded response vector for this response is \\(0,0,1,0,0)^T\\.
Assume that the IRT model predicts the category probabilities as \\(0.1,
0.2, 0.4, 0.25, 0.05)^T\\, so that the expected item score is \\0(0.1) +
1(0.2) + 2(0.4) + 3(0.25) + 4(0.05) = 1.95\\. In the RDIF-CR framework,
the categorical residual vector is calculated by subtracting the
predicted probability vector from the one-hot response vector, resulting
in \\(-0.1, -0.2, 0.6, -0.25, -0.05)^T\\.

In contrast, net DIF is assessed with a single item score residual, the
observed item score minus the expected item score. In this example, the
residual would be \\2 - 1.95 = 0.05\\. To detect net DIF, use
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md).

\\RDIF\_{R}-CR\\ is the quadratic form of the vector of category-level
differences in mean raw residuals between the focal and reference
groups, centered at its null mean, in a generalized inverse of its null
covariance matrix. \\RDIF\_{S}-CR\\ is the corresponding statistic for
the mean squared residuals, and \\RDIF\_{RS}-CR\\ uses the two vectors
combined. Under the null hypothesis of no DIF, each statistic
asymptotically follows a chi-square distribution whose degrees of
freedom equal the rank of the covariance matrix (Moore, 1977). For an
item with K score categories (K of 3 or more), the K category residuals
of each examinee sum to zero, so the elements of the \\RDIF\_{R}-CR\\
vector also sum to zero and its covariance matrix has rank K - 1.
\\RDIF\_{R}-CR\\ therefore has K - 1 degrees of freedom. It is computed
after dropping one score category, which gives the same value as the
Moore-Penrose generalized inverse of the full covariance matrix, and the
same category is dropped from the raw residual part of
\\RDIF\_{RS}-CR\\. The squared residuals are not subject to a constraint
that holds for all examinees, so \\RDIF\_{S}-CR\\ has K degrees of
freedom, and \\RDIF\_{RS}-CR\\ has 2K - 1. For a dichotomous item, the
two category residuals are the negatives of each other, so only the
residual for score category 1 is used, and the three statistics have 1,
1, and 2 degrees of freedom; then \\RDIF\_{R}-CR\\ and \\RDIF\_{S}-CR\\
equal the squares of the standardized \\RDIF\_{R}\\ and \\RDIF\_{S}\\
statistics of
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
\\RDIF\_{RS}-CR\\ equals \\RDIF\_{RS}\\, and the p-values and flagged
items are identical to those of
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md).

In rare cases, for example when all examinees who responded to an item
have the same category probabilities, a covariance matrix has an even
lower rank. Then the statistic is computed with the Moore-Penrose
generalized inverse and compared with a chi-square distribution whose
degrees of freedom equal the rank, and a warning names the item. The
columns `df.crdifr`, `df.crdifs`, and `df.crdifrs` of `dif_stat` report
the degrees of freedom used.

The purification procedure, the `min.resp` argument, the handling of
examinees with `NA` ability estimates, the coding of `group`, and the
`NaN` results of an item to which no examinee of one group responded
follow [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) (see
its Details). The p-value, not the statistic, is used to choose the item
to remove because the degrees of freedom of the RDIF-CR statistics
differ between items with different numbers of score categories.

The print method of this class has an argument `digits` that is not
used; the statistics are printed rounded to three decimal places.

## Methods (by class)

- `crdif(default)`: Default method for a data frame `x` containing item
  metadata.

- `crdif(est_irt)`: Method for an object of class `est_irt` created by
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md). The
  response data and the scaling factor `D` are taken from the object.

- `crdif(est_item)`: Method for an object of class `est_item` created by
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).
  The response data, ability estimates, and the scaling factor `D` are
  taken from the object.

## Methods (by generic)

- `print(crdif)`: Print method for an object of class `crdif` (given as
  `x`). It prints the flagged items and the RDIF-CR statistics of the
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

Penfield, R. D. (2010). Distinguishing between net and global DIF in
polytomous items. *Journal of Educational Measurement, 47*(2), 129-149.
[doi:10.1111/j.1745-3984.2010.00105.x](https://doi.org/10.1111/j.1745-3984.2010.00105.x)
.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627) .

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{

############################################################################
# This example demonstrates how to detect global DIF in polytomous items
# using the RDIF-CR framework implemented in `irtQ::crdif()`.
# Simulated response data are generated from 5 GRM items with 4 score
# categories. DIF is introduced in the 1st and 5th items.
############################################################################

###############################################
# (1) Simulate response data with DIF
###############################################

set.seed(1)

# Generate ability parameters for 1000 examinees in each group
# Reference and focal groups follow N(0, 1.5^2)
theta_ref <- rnorm(1000, 0, 1.5)
theta_foc <- rnorm(1000, 0, 1.5)

# Combine abilities from both groups
theta_all <- c(theta_ref, theta_foc)

# Define item parameters using `irtQ::shape_df()`
# Items 1 and 5 are intentionally modified to exhibit DIF
par_ref <- irtQ::shape_df(
  par.prm = list(
    a = c(1, 1, 1, 2, 2),
    d = list(c(-2, 0, 1),
             c(-2, 0, 2),
             c(-2, 0, 1),
             c(-1, 0, 2),
             c(-2, 0, 0.5))
  ),
  cats = 4, model = "GRM"
)

par_foc <- irtQ::shape_df(
  par.prm = list(
    a = c(2, 1, 1, 2, 0.5),
    d = list(c(-0.5, 0, 0.5),
             c(-2, 0, 2),
             c(-2, 0, 1),
             c(-1, 0, 2),
             c(-1.5, -1, 0))
  ),
  cats = 4, model = "GRM"
)

# Generate response data
resp_ref <- irtQ::simdat(x = par_ref, theta = theta_ref, D = 1)
resp_foc <- irtQ::simdat(x = par_foc, theta = theta_foc, D = 1)

# Combine response data across groups
data <- rbind(resp_ref, resp_foc)

###############################################
# (2) Estimate item and ability parameters
###############################################

# Estimate GRM item parameters using `irtQ::est_irt()`
fit_mod <- irtQ::est_irt(data = data, D = 1, model = "GRM", cats = 4)
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -13079.1582, Max-Change: 0.716758 EM iteration: 2, Loglike: -12041.8257, Max-Change: 0.450889 EM iteration: 3, Loglike: -11967.1516, Max-Change: 0.304131 EM iteration: 4, Loglike: -11944.2648, Max-Change: 0.209059 EM iteration: 5, Loglike: -11935.3740, Max-Change: 0.146712 EM iteration: 6, Loglike: -11931.5055, Max-Change: 0.10458 EM iteration: 7, Loglike: -11929.6966, Max-Change: 0.075403 EM iteration: 8, Loglike: -11928.8066, Max-Change: 0.054821 EM iteration: 9, Loglike: -11928.3510, Max-Change: 0.040101 EM iteration: 10, Loglike: -11928.1100, Max-Change: 0.029467 EM iteration: 11, Loglike: -11927.9790, Max-Change: 0.021725 EM iteration: 12, Loglike: -11927.9059, Max-Change: 0.016058 EM iteration: 13, Loglike: -11927.8642, Max-Change: 0.011891 EM iteration: 14, Loglike: -11927.8399, Max-Change: 0.008818 EM iteration: 15, Loglike: -11927.8254, Max-Change: 0.006545 EM iteration: 16, Loglike: -11927.8167, Max-Change: 0.004862 EM iteration: 17, Loglike: -11927.8113, Max-Change: 0.003614 EM iteration: 18, Loglike: -11927.8080, Max-Change: 0.002687 EM iteration: 19, Loglike: -11927.8058, Max-Change: 0.001998 EM iteration: 20, Loglike: -11927.8045, Max-Change: 0.001486 EM iteration: 21, Loglike: -11927.8036, Max-Change: 0.001105 EM iteration: 22, Loglike: -11927.8030, Max-Change: 0.000822 EM iteration: 23, Loglike: -11927.8026, Max-Change: 0.000611 EM iteration: 24, Loglike: -11927.8024, Max-Change: 0.000454 EM iteration: 25, Loglike: -11927.8022, Max-Change: 0.000338 EM iteration: 26, Loglike: -11927.8021, Max-Change: 0.000251 EM iteration: 27, Loglike: -11927.8020, Max-Change: 0.000186 EM iteration: 28, Loglike: -11927.8019, Max-Change: 0.000138 EM iteration: 29, Loglike: -11927.8019, Max-Change: 0.000103 EM iteration: 30, Loglike: -11927.8019, Max-Change: 7.6e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 0.28 seconds. 

# Extract estimated item parameters
x <- fit_mod$par.est

# Estimate ability scores using ML method
score <- est_score(x = x, data = data, method = "ML")$est.theta

###############################################
# (3) Perform RDIF-CR DIF analysis
###############################################

# Define group membership: 1 = focal group
group <- c(rep(0, 1000), rep(1, 1000))

# (a) DIF detection without purification
dif_nopuri <- crdif(
  x = x, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05
)
print(dif_nopuri)
#> 
#> Call:
#> crdif.default(x = x, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05)
#> 
#> DIF analysis using three RDIF-CR statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R)-CR: 
#>     1, 5 
#>   - DIF Items identified by RDIF(S)-CR: 
#>     1, 5 
#>   - DIF Items identified by RDIF(RS)-CR: 
#>     1, 5 
#>   - RDIF-CR Statistics: 
#> 
#>   id n.ref n.foc crdifr df.crdifr p.crdifr     crdifs df.crdifs p.crdifs    
#> 1 V1  1000  1000 181.97         3    0.000 *** 156.06         4    0.000 ***
#> 2 V2  1000  1000   2.00         3    0.574       8.19         4    0.085   .
#> 3 V3  1000  1000   3.86         3    0.277       2.68         4    0.613    
#> 4 V4  1000  1000   5.57         3    0.134       5.29         4    0.258    
#> 5 V5  1000  1000 409.88         3    0.000 *** 606.02         4    0.000 ***
#>   crdifrs df.crdifrs p.crdifrs    
#> 1  194.66          7     0.000 ***
#> 2    9.25          7     0.235    
#> 3    8.57          7     0.285    
#> 4   12.15          7     0.096   .
#> 5  664.04          7     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (b) DIF detection with purification using RDIF_{RS}-CR (not run: it takes
#     about 15 seconds because the abilities of the examinees are
#     re-estimated at every iteration)
#
#   dif_puri <- crdif(
#     x = x, data = data, score = score,
#     group = group, focal.name = 1, D = 1, alpha = 0.05,
#     purify = TRUE, purify.by = "crdifrs"
#   )
#   print(dif_puri)
#
#   purify.by = "crdifr" or "crdifs" can be used in the same way.

# }
```
