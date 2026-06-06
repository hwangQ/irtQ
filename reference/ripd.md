# Residual-based Item Parameter Drift (RIPD) Detection Framework

This function computes three RIPD statistics—\\RIPD\_{R}\\,
\\RIPD\_{S}\\, and \\RIPD\_{RS}\\—for each item. \\RIPD\_{R}\\ captures
differences in mean raw residuals between groups, which is typically
indicative of uniform item parameter drift (IPD). \\RIPD\_{S}\\ captures
differences in mean squared residuals between groups, reflecting
nonuniform IPD. \\RIPD\_{RS}\\, a combined chi-square-based index, is
sensitive to both uniform and nonuniform IPD.

## Usage

``` r
ripd(x, ...)

# Default S3 method
ripd(
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
  purify.by = c("ripdrs", "ripdr", "ripds"),
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
ripd(
  x,
  score = NULL,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("ripdrs", "ripdr", "ripds"),
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
ripd(
  x,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("ripdrs", "ripdr", "ripds"),
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
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.); or an object of class `est_irt`
  obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md), or
  `est_item` from
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

  See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
  or [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- ...:

  Additional arguments passed to the
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  function.

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- score:

  A numeric vector containing examinees' ability estimates (theta
  values). If not provided, `ripd()` will estimate ability parameters
  internally before computing the RIPD statistics. See
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  for more information on scoring methods. Default is `NULL`.

- group:

  A numeric or character vector indicating group membership of
  examinees. The length of vector should the same with the number of
  rows in the response data matrix.

- focal.name:

  A single numeric or character value specifying the focal group. For
  instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating the
  focal group, set `focal.name = 1`.

- item.skip:

  A numeric vector of item indices to exclude from IPD analysis. If
  NULL, all items are included. Useful for omitting specific items based
  on prior insights.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for
  hypothesis testing using the RIPD statistics. Default is `0.05`.

- missing:

  A value indicating missing values in the response data set. Default is
  NA.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying which RIPD statistic is used to perform
  the purification. Available options are "ripdrs" for \\RIPD\_{RS}\\,
  "ripdr" for \\RIPD\_{R}\\, and "ripds" for \\RIPD\_{S}\\.

- max.iter:

  A positive integer specifying the maximum number of iterations allowed
  for the purification process. Default is `10`.

- min.resp:

  A positive integer specifying the minimum number of valid item
  responses required from an examinee in order to compute an ability
  estimate. Default is `NULL`.

- method:

  A character string indicating the scoring method to use. Available
  options are:

  - `"ML"`: Maximum likelihood estimation

  - `"WL"`: Weighted likelihood estimation (Warm, 1989)

  - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)

  - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)

  Default is `"ML"`.

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

  An integer specifying the number of logical CPU cores to use for
  parallel processing. Default is `1`. See
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  for details.

- verbose:

  Logical. If `TRUE`, progress messages from the purification procedure
  will be displayed; if `FALSE`, the messages will be suppressed.
  Default is `TRUE`.

## Value

This function returns a list containing four main components:

- no_purify:

  A list of sub-objects containing the results of IPD analysis without
  applying a purification procedure. The sub-objects include:

  ipd_stat

  :   A data frame summarizing the RIPD analysis results for all items.
      The columns include: item ID, \\RIPD\_{R}\\ statistic,
      standardized \\RIPD\_{R}\\, \\RIPD\_{S}\\ statistic, standardized
      \\RIPD\_{S}\\, \\RIPD\_{RS}\\ statistic, p-values for
      \\RIPD\_{R}\\, \\RIPD\_{S}\\, and \\RIPD\_{RS}\\, sample sizes for
      the reference and focal groups, and total sample size. Note that
      \\RIPD\_{RS}\\ does not have a standardized value because it is a
      \\\chi^{2}\\-based statistic.

  moments

  :   A data frame reporting the first and second moments of the RIPD
      statistics. The columns include: item ID, mean and standard
      deviation of \\RIPD\_{R}\\, mean and standard deviation of
      \\RIPD\_{S}\\, and the covariance between \\RIPD\_{R}\\ and
      \\RIPD\_{S}\\.

  ipd_item

  :   A list of three numeric vectors identifying items flagged as
      drifting by each RIPD statistic: \\RIPD\_{R}\\, \\RIPD\_{S}\\, and
      \\RIPD\_{RS}\\.

  score

  :   A numeric vector of ability estimates used to compute the RIPD
      statistics.

- purify:

  A logical value indicating whether the purification procedure was
  applied.

- with_purify:

  A list of sub-objects containing the results of IPD analysis with a
  purification procedure. The sub-objects include:

  purify.by

  :   A character string indicating the RIPD statistic used for
      purification. Possible values are "ripdr", "ripds", and "ripdrs",
      corresponding to \\RIPD\_{R}\\, \\RIPD\_{S}\\, and \\RIPD\_{RS}\\,
      respectively.

  ipd_stat

  :   A data frame reporting the RIPD analysis results for all items
      from the final iteration. Same structure as in `no_purify`, with
      one additional column indicating the iteration number in which
      each result was obtained.

  moments

  :   A data frame reporting the moments of RIPD statistics from the
      final iteration. Includes the same columns as in `no_purify`, with
      an additional column for the iteration number.

  ipd_item

  :   A list of numeric item indices identified as IPD items across
      iterations.

  n.iter

  :   An integer indicating the total number of iterations performed
      during the purification process.

  score

  :   A numeric vector of purified ability estimates used to compute the
      final RIPD statistics.

  complete

  :   A logical value indicating whether the purification process
      converged. If `FALSE`, the maximum number of iterations was
      reached without convergence.

- alpha:

  A numeric value indicating the significance level (\\\alpha\\) used in
  hypothesis testing for RIPD statistics.

## Methods (by class)

- `ripd(default)`: Default method for computing the three RIPD
  statistics using a data frame `x` that contains item metadata

- `ripd(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `ripd(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

## References

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Lim, H., & H, K. T. (2025, April). *IRT residual-based approach to
detecting item parameter drift in CAT*. Paper presented at the annual
conference of the National Council on Measurement in Education (NCME),
Denver, CO

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>
