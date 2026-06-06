# IRT Residual-Based Differential Item Functioning (RDIF) Detection Framework

This function computes three RDIF statistics for each item:
\\RDIF\_{R}\\, \\RDIF\_{S}\\, and \\RDIF\_{RS}\\ (Lim & Choe, 2023; Lim,
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
  values). If not provided, `rdif()` will estimate ability parameters
  internally before computing the RDIF statistics. See
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  for more information on scoring methods. Default is `NULL`.

- group:

  A numeric or character vector indicating examinees' group membership.
  The length of the vector must match the number of rows in the response
  data matrix.

- focal.name:

  A single numeric or character value specifying the focal group. For
  instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating the
  focal group, set `focal.name = 1`.

- item.skip:

  A numeric vector of item indices to exclude from DIF analysis. If
  `NULL`, all items are included. Useful for omitting specific items
  based on prior insights.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for
  hypothesis testing using the RDIF statistics. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying which RDIF statistic is used to perform
  the purification. Available options are "rdifrs" for \\RDIF\_{RS}\\,
  "rdifr" for \\RDIF\_{R}\\, and "rdifs" for \\RDIF\_{S}\\.

- max.iter:

  A positive integer specifying the maximum number of iterations allowed
  for the purification process. Default is `10`.

- min.resp:

  A positive integer specifying the minimum number of valid item
  responses required from an examinee in order to compute an ability
  estimate. Default is `NULL`. See **Details** for more information.

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

  A list of sub-objects containing the results of DIF analysis without
  applying a purification procedure. The sub-objects include:

  dif_stat

  :   A data frame summarizing the RDIF analysis results for all items.
      The columns include: item ID, \\RDIF\_{R}\\ statistic,
      standardized \\RDIF\_{R}\\, \\RDIF\_{S}\\ statistic, standardized
      \\RDIF\_{S}\\, \\RDIF\_{RS}\\ statistic, p-values for
      \\RDIF\_{R}\\, \\RDIF\_{S}\\, and \\RDIF\_{RS}\\, sample sizes for
      the reference and focal groups, and total sample size. Note that
      \\RDIF\_{RS}\\ does not have a standardized value because it is a
      \\\chi^{2}\\-based statistic.

  moments

  :   A data frame reporting the first and second moments of the RDIF
      statistics. The columns include: item ID, mean and standard
      deviation of \\RDIF\_{R}\\, mean and standard deviation of
      \\RDIF\_{S}\\, and the covariance between \\RDIF\_{R}\\ and
      \\RDIF\_{S}\\.

  dif_item

  :   A list of three numeric vectors identifying items flagged as DIF
      by each RDIF statistic: \\RDIF\_{R}\\, \\RDIF\_{S}\\, and
      \\RDIF\_{RS}\\.

  score

  :   A numeric vector of ability estimates used to compute the RDIF
      statistics.

- purify:

  A logical value indicating whether the purification procedure was
  applied.

- with_purify:

  A list of sub-objects containing the results of DIF analysis with a
  purification procedure. The sub-objects include:

  purify.by

  :   A character string indicating the RDIF statistic used for
      purification. Possible values are "rdifr", "rdifs", and "rdifrs",
      corresponding to \\RDIF\_{R}\\, \\RDIF\_{S}\\, and \\RDIF\_{RS}\\,
      respectively.

  dif_stat

  :   A data frame reporting the RDIF analysis results for all items
      across the final iteration. Same structure as in `no_purify`, with
      one additional column indicating the iteration number in which
      each result was obtained.

  moments

  :   A data frame reporting the moments of RDIF statistics across the
      final iteration. Includes the same columns as in `no_purify`, with
      an additional column for the iteration number.

  dif_item

  :   A list of three numeric vectors identifying DIF items flagged by
      each RDIF statistic.

  n.iter

  :   An integer indicating the total number of iterations performed
      during the purification process.

  score

  :   A numeric vector of purified ability estimates used to compute the
      final RDIF statistics.

  complete

  :   A logical value indicating whether the purification process
      converged. If FALSE, the maximum number of iterations was reached
      without convergence.

- alpha:

  A numeric value indicating the significance level (\\\alpha\\) used in
  hypothesis testing for RDIF statistics.

## Details

The RDIF framework (Lim & Choe, 2023; Lim et al., 2022) consists of
three IRT residual-based statistics: \\RDIF\_{R}\\, \\RDIF\_{S}\\, and
\\RDIF\_{RS}\\. Under the null hypothesis that a test contains no DIF
items, \\RDIF\_{R}\\ and \\RDIF\_{S}\\ asymptotically follow standard
normal distributions. \\RDIF\_{RS}\\ is based on a bivariate normal
distribution of the \\RDIF\_{R}\\ and \\RDIF\_{S}\\ statistics, and
under the null hypothesis, it asymptotically follows a \\\chi^{2}\\
distribution with 2 degrees of freedom. See Lim et al. (2022) for more
details about the RDIF framework.

The `rdif()` function computes all three RDIF statistics: \\RDIF\_{R}\\,
\\RDIF\_{S}\\, and \\RDIF\_{RS}\\. The current version of `rdif()`
supports both dichotomous and polytomous item response data. Note that
for polytomous items, net DIF are assessed. To evaluate global DIF for
polytomous items, use
[`crdif()`](https://hwangQ.github.io/irtQ/reference/crdif.md) function.

To compute the RDIF statistics, the `rdif()` function requires: (1) item
parameter estimates obtained from aggregate data (regardless of group
membership), (2) examinees' ability estimates (e.g., ML), and (3)
examinees' item response data. Note that the ability estimates must be
based on the aggregate-data item parameters. The item parameter
estimates should be provided in the `x` argument, the ability estimates
in the `score` argument, and the response data in the `data` argument.
If ability estimates are not provided (i.e., `score = NULL`), `rdif()`
will estimate them automatically using the scoring method specified via
the `method` argument (e.g., `method = "ML"`).

The `group` argument should be a vector containing exactly two distinct
values (either numeric or character), representing the reference and
focal groups. Its length must match the number of rows in the response
data, where each element corresponds to an examinee. Once `group` is
specified, a single numeric or character value must be provided in the
`focal.name` argument to indicate which level in `group` represents the
focal group.

Similar to other DIF detection approaches, the RDIF framework supports
an iterative purification process. When `purify = TRUE`, purification is
conducted using one of the RDIF statistics specified in the `purify.by`
argument (e.g., `purify.by = "rdifrs"`). At each iteration, examinees'
ability estimates are recalculated based on the set of purified items
using the scoring method specified in the `method` argument. The
purification process continues until no additional DIF items are
identified or the maximum number of iterations specified in `max.iter`
is reached. See Lim et al. (2022) for more details on the purification
procedure.

Scoring based on a small number of item responses can lead to large
standard errors, potentially reducing the accuracy of DIF detection in
the RDIF framework. The `min.resp` argument can be used to exclude
examinees with insufficient response data from scoring, especially
during the purification process. For example, if `min.resp` is not NULL
(e.g., `min.resp = 5`), examinees who responded to fewer than five items
will have all their responses treated as missing (i.e., NA). As a
result, their ability estimates will also be missing and will not be
used in the computation of RDIF statistics. If `min.resp = NULL`, a
score will be computed for any examinee with at least one valid item
response.

## Methods (by class)

- `rdif(default)`: Default method for computing the three RDIF
  statistics using a data frame `x` that contains item metadata

- `rdif(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `rdif(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

## References

Lim, H., & Choe, E. M. (2023). Detecting differential item functioning
in CAT using IRT residual DIF approach. *Journal of Educational
Measurement, 60*(4), 626-650.
[doi:10.1111/jedm.12366](https://doi.org/10.1111/jedm.12366) .

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

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
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = 0.15),
    item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
  )

# Add uniform DIF by shifting the b-parameters for the focal group
difpar_foc <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))

# Combine the DIF and non-DIF items for both reference and focal groups
# Therefor, the first 4 items exhibit uniform DIF
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
#>  EM iteration: 1, Loglike: -25510.5305, Max-Change: 1.315321 EM iteration: 2, Loglike: -22606.9434, Max-Change: 0.354576 EM iteration: 3, Loglike: -22588.3376, Max-Change: 0.156275 EM iteration: 4, Loglike: -22586.3606, Max-Change: 0.075037 EM iteration: 5, Loglike: -22585.8560, Max-Change: 0.039037 EM iteration: 6, Loglike: -22585.6468, Max-Change: 0.022879 EM iteration: 7, Loglike: -22585.5393, Max-Change: 0.014254 EM iteration: 8, Loglike: -22585.4778, Max-Change: 0.009304 EM iteration: 9, Loglike: -22585.4400, Max-Change: 0.006289 EM iteration: 10, Loglike: -22585.4156, Max-Change: 0.004361 EM iteration: 11, Loglike: -22585.3992, Max-Change: 0.00308 EM iteration: 12, Loglike: -22585.3877, Max-Change: 0.002204 EM iteration: 13, Loglike: -22585.3795, Max-Change: 0.001592 EM iteration: 14, Loglike: -22585.3736, Max-Change: 0.001277 EM iteration: 15, Loglike: -22585.3692, Max-Change: 0.001067 EM iteration: 16, Loglike: -22585.3659, Max-Change: 0.000905 EM iteration: 17, Loglike: -22585.3634, Max-Change: 0.000772 EM iteration: 18, Loglike: -22585.3616, Max-Change: 0.000663 EM iteration: 19, Loglike: -22585.3602, Max-Change: 0.000572 EM iteration: 20, Loglike: -22585.3592, Max-Change: 0.000496 EM iteration: 21, Loglike: -22585.3584, Max-Change: 0.000431 EM iteration: 22, Loglike: -22585.3579, Max-Change: 0.000377 EM iteration: 23, Loglike: -22585.3575, Max-Change: 0.00033 EM iteration: 24, Loglike: -22585.3572, Max-Change: 0.00029 EM iteration: 25, Loglike: -22585.3571, Max-Change: 0.000255 EM iteration: 26, Loglike: -22585.3570, Max-Change: 0.000225 EM iteration: 27, Loglike: -22585.3569, Max-Change: 0.000198 EM iteration: 28, Loglike: -22585.3569, Max-Change: 0.000175 EM iteration: 29, Loglike: -22585.3569, Max-Change: 0.000155 EM iteration: 30, Loglike: -22585.3569, Max-Change: 0.000137 EM iteration: 31, Loglike: -22585.3569, Max-Change: 0.000122 EM iteration: 32, Loglike: -22585.3569, Max-Change: 0.000108 EM iteration: 33, Loglike: -22585.3570, Max-Change: 9.6e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 1.99 seconds. 
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
#>     5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11, 14, 19 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs    
#> 1   V1   500   500 -0.110   0.000 ***  0.008   0.204    13.750    0.001 ***
#> 2   V2   500   500 -0.175   0.000 *** -0.020   0.163    42.916    0.000 ***
#> 3   V3   500   500 -0.059   0.044   * -0.003   0.688     4.261    0.119    
#> 4   V4   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000 ***
#> 5   V5   500   500 -0.006   0.830      0.021   0.003 **  9.607    0.008  **
#> 6   V6   500   500  0.012   0.534      0.000   0.463     0.579    0.749    
#> 7   V7   500   500  0.018   0.533      0.005   0.470     0.854    0.653    
#> 8   V8   500   500  0.028   0.317      0.002   0.887     1.329    0.514    
#> 9   V9   500   500 -0.012   0.689      0.008   0.499     0.495    0.781    
#> 10 V10   500   500  0.031   0.229      0.011   0.390     1.536    0.464    
#> 11 V11   500   500  0.066   0.023   * -0.008   0.638     8.193    0.017   *
#> 12 V12   500   500  0.002   0.961      0.002   0.613     0.296    0.863    
#> 13 V13   500   500 -0.010   0.737     -0.003   0.750     0.272    0.873    
#> 14 V14   500   500  0.068   0.013   * -0.008   0.427     6.525    0.038   *
#> 15 V15   500   500  0.019   0.488     -0.006   0.404     0.752    0.687    
#> 16 V16   500   500  0.007   0.809     -0.016   0.324     2.394    0.302    
#> 17 V17   500   500  0.011   0.696     -0.008   0.690     1.180    0.554    
#> 18 V18   500   500 -0.018   0.505      0.007   0.608     0.554    0.758    
#> 19 V19   500   500  0.052   0.055   . -0.021   0.034  *  6.103    0.047   *
#> 20 V20   500   500 -0.002   0.936     -0.019   0.113     2.510    0.285    
#> 21 V21   500   500  0.018   0.499     -0.003   0.655     0.520    0.771    
#> 22 V22   500   500  0.028   0.358     -0.005   0.578     1.837    0.399    
#> 23 V23   500   500  0.015   0.414     -0.005   0.172     1.879    0.391    
#> 24 V24   500   500  0.000   0.980      0.014   0.601     0.813    0.666    
#> 25 V25   500   500  0.020   0.350     -0.004   0.318     1.062    0.588    
#> 26 V26   500   500  0.024   0.383      0.002   0.851     1.482    0.476    
#> 27 V27   500   500  0.024   0.385     -0.001   0.753     0.758    0.685    
#> 28 V28   500   500  0.024   0.412     -0.004   0.864     1.477    0.478    
#> 29 V29   500   500  0.006   0.784      0.006   0.991     0.509    0.775    
#> 30 V30   500   500  0.028   0.176     -0.007   0.254     1.839    0.399    
#> 31 V31   500   500  0.014   0.635     -0.002   0.803     0.252    0.881    
#> 32 V32   500   500  0.028   0.270     -0.008   0.238     2.153    0.341    
#> 33 V33   500   500 -0.026   0.217      0.019   0.402     1.643    0.440    
#> 34 V34   500   500 -0.010   0.734     -0.009   0.290     1.392    0.499    
#> 35 V35   500   500  0.011   0.713     -0.017   0.093  .  3.930    0.140    
#> 36 V36   500   500 -0.019   0.395      0.021   0.268     1.227    0.541    
#> 37 V37   500   500  0.007   0.743      0.003   0.614     0.258    0.879    
#> 38 V38   500   500 -0.051   0.078   . -0.008   0.427     3.844    0.146    
#> 39 V39   500   500  0.025   0.366      0.006   0.933     1.022    0.600    
#> 40 V40   500   500 -0.013   0.645     -0.008   0.622     0.305    0.859    
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
#>     5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11, 14, 19 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs    
#> 1   V1   500   500 -0.110   0.000 ***  0.008   0.204    13.750    0.001 ***
#> 2   V2   500   500 -0.175   0.000 *** -0.020   0.163    42.916    0.000 ***
#> 3   V3   500   500 -0.059   0.044   * -0.003   0.688     4.261    0.119    
#> 4   V4   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000 ***
#> 5   V5   500   500 -0.006   0.830      0.021   0.003 **  9.607    0.008  **
#> 6   V6   500   500  0.012   0.534      0.000   0.463     0.579    0.749    
#> 7   V7   500   500  0.018   0.533      0.005   0.470     0.854    0.653    
#> 8   V8   500   500  0.028   0.317      0.002   0.887     1.329    0.514    
#> 9   V9   500   500 -0.012   0.689      0.008   0.499     0.495    0.781    
#> 10 V10   500   500  0.031   0.229      0.011   0.390     1.536    0.464    
#> 11 V11   500   500  0.066   0.023   * -0.008   0.638     8.193    0.017   *
#> 12 V12   500   500  0.002   0.961      0.002   0.613     0.296    0.863    
#> 13 V13   500   500 -0.010   0.737     -0.003   0.750     0.272    0.873    
#> 14 V14   500   500  0.068   0.013   * -0.008   0.427     6.525    0.038   *
#> 15 V15   500   500  0.019   0.488     -0.006   0.404     0.752    0.687    
#> 16 V16   500   500  0.007   0.809     -0.016   0.324     2.394    0.302    
#> 17 V17   500   500  0.011   0.696     -0.008   0.690     1.180    0.554    
#> 18 V18   500   500 -0.018   0.505      0.007   0.608     0.554    0.758    
#> 19 V19   500   500  0.052   0.055   . -0.021   0.034  *  6.103    0.047   *
#> 20 V20   500   500 -0.002   0.936     -0.019   0.113     2.510    0.285    
#> 21 V21   500   500  0.018   0.499     -0.003   0.655     0.520    0.771    
#> 22 V22   500   500  0.028   0.358     -0.005   0.578     1.837    0.399    
#> 23 V23   500   500  0.015   0.414     -0.005   0.172     1.879    0.391    
#> 24 V24   500   500  0.000   0.980      0.014   0.601     0.813    0.666    
#> 25 V25   500   500  0.020   0.350     -0.004   0.318     1.062    0.588    
#> 26 V26   500   500  0.024   0.383      0.002   0.851     1.482    0.476    
#> 27 V27   500   500  0.024   0.385     -0.001   0.753     0.758    0.685    
#> 28 V28   500   500  0.024   0.412     -0.004   0.864     1.477    0.478    
#> 29 V29   500   500  0.006   0.784      0.006   0.991     0.509    0.775    
#> 30 V30   500   500  0.028   0.176     -0.007   0.254     1.839    0.399    
#> 31 V31   500   500  0.014   0.635     -0.002   0.803     0.252    0.881    
#> 32 V32   500   500  0.028   0.270     -0.008   0.238     2.153    0.341    
#> 33 V33   500   500 -0.026   0.217      0.019   0.402     1.643    0.440    
#> 34 V34   500   500 -0.010   0.734     -0.009   0.290     1.392    0.499    
#> 35 V35   500   500  0.011   0.713     -0.017   0.093  .  3.930    0.140    
#> 36 V36   500   500 -0.019   0.395      0.021   0.268     1.227    0.541    
#> 37 V37   500   500  0.007   0.743      0.003   0.614     0.258    0.879    
#> 38 V38   500   500 -0.051   0.078   . -0.008   0.427     3.844    0.146    
#> 39 V39   500   500  0.025   0.366      0.006   0.933     1.022    0.600    
#> 40 V40   500   500 -0.013   0.645     -0.008   0.622     0.305    0.859    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (b)-1 Compute RDIF statistics with purification based on RDIF(R)
dif_puri_r <- rdif(
  x = est_par, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "rdifr"
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 Iteration: 3 Iteration: 4 Iteration: 5 
#> Purification is finished. 
print(dif_puri_r)
#> 
#> Call:
#> rdif.default(x = est_par, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, purify = TRUE, purify.by = "rdifr")
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 3, 4, 11, 14 
#>   - DIF Items identified by RDIF(S): 
#>     5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11, 14, 19 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs    
#> 1   V1   500   500 -0.110   0.000 ***  0.008   0.204    13.750    0.001 ***
#> 2   V2   500   500 -0.175   0.000 *** -0.020   0.163    42.916    0.000 ***
#> 3   V3   500   500 -0.059   0.044   * -0.003   0.688     4.261    0.119    
#> 4   V4   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000 ***
#> 5   V5   500   500 -0.006   0.830      0.021   0.003 **  9.607    0.008  **
#> 6   V6   500   500  0.012   0.534      0.000   0.463     0.579    0.749    
#> 7   V7   500   500  0.018   0.533      0.005   0.470     0.854    0.653    
#> 8   V8   500   500  0.028   0.317      0.002   0.887     1.329    0.514    
#> 9   V9   500   500 -0.012   0.689      0.008   0.499     0.495    0.781    
#> 10 V10   500   500  0.031   0.229      0.011   0.390     1.536    0.464    
#> 11 V11   500   500  0.066   0.023   * -0.008   0.638     8.193    0.017   *
#> 12 V12   500   500  0.002   0.961      0.002   0.613     0.296    0.863    
#> 13 V13   500   500 -0.010   0.737     -0.003   0.750     0.272    0.873    
#> 14 V14   500   500  0.068   0.013   * -0.008   0.427     6.525    0.038   *
#> 15 V15   500   500  0.019   0.488     -0.006   0.404     0.752    0.687    
#> 16 V16   500   500  0.007   0.809     -0.016   0.324     2.394    0.302    
#> 17 V17   500   500  0.011   0.696     -0.008   0.690     1.180    0.554    
#> 18 V18   500   500 -0.018   0.505      0.007   0.608     0.554    0.758    
#> 19 V19   500   500  0.052   0.055   . -0.021   0.034  *  6.103    0.047   *
#> 20 V20   500   500 -0.002   0.936     -0.019   0.113     2.510    0.285    
#> 21 V21   500   500  0.018   0.499     -0.003   0.655     0.520    0.771    
#> 22 V22   500   500  0.028   0.358     -0.005   0.578     1.837    0.399    
#> 23 V23   500   500  0.015   0.414     -0.005   0.172     1.879    0.391    
#> 24 V24   500   500  0.000   0.980      0.014   0.601     0.813    0.666    
#> 25 V25   500   500  0.020   0.350     -0.004   0.318     1.062    0.588    
#> 26 V26   500   500  0.024   0.383      0.002   0.851     1.482    0.476    
#> 27 V27   500   500  0.024   0.385     -0.001   0.753     0.758    0.685    
#> 28 V28   500   500  0.024   0.412     -0.004   0.864     1.477    0.478    
#> 29 V29   500   500  0.006   0.784      0.006   0.991     0.509    0.775    
#> 30 V30   500   500  0.028   0.176     -0.007   0.254     1.839    0.399    
#> 31 V31   500   500  0.014   0.635     -0.002   0.803     0.252    0.881    
#> 32 V32   500   500  0.028   0.270     -0.008   0.238     2.153    0.341    
#> 33 V33   500   500 -0.026   0.217      0.019   0.402     1.643    0.440    
#> 34 V34   500   500 -0.010   0.734     -0.009   0.290     1.392    0.499    
#> 35 V35   500   500  0.011   0.713     -0.017   0.093  .  3.930    0.140    
#> 36 V36   500   500 -0.019   0.395      0.021   0.268     1.227    0.541    
#> 37 V37   500   500  0.007   0.743      0.003   0.614     0.258    0.879    
#> 38 V38   500   500 -0.051   0.078   . -0.008   0.427     3.844    0.146    
#> 39 V39   500   500  0.025   0.366      0.006   0.933     1.022    0.600    
#> 40 V40   500   500 -0.013   0.645     -0.008   0.622     0.305    0.859    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 5
#>   - RDIF statistic used for purification: RDIF(R)
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 3, 4, 38 
#>   - RDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs
#> 1   V1      2   500   500 -0.118   0.000 ***  0.009   0.123    16.190    0.000
#> 2   V2      1   500   500 -0.184   0.000 *** -0.020   0.162    47.376    0.000
#> 3   V3      3   500   500 -0.073   0.013   * -0.002   0.985     6.178    0.046
#> 4   V4      0   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000
#> 5   V5      5   500   500 -0.018   0.546      0.020   0.003 ** 10.858    0.004
#> 6   V6      5   500   500 -0.001   0.963     -0.001   0.812     0.121    0.941
#> 7   V7      5   500   500  0.002   0.951      0.007   0.302     1.065    0.587
#> 8   V8      5   500   500  0.013   0.640      0.000   0.912     0.454    0.797
#> 9   V9      5   500   500 -0.028   0.324      0.009   0.270     1.669    0.434
#> 10 V10      5   500   500  0.008   0.761      0.010   0.416     0.673    0.714
#> 11 V11      5   500   500  0.051   0.081   . -0.006   0.747     4.697    0.096
#> 12 V12      5   500   500 -0.010   0.732      0.001   0.706     0.466    0.792
#> 13 V13      5   500   500 -0.023   0.450     -0.002   0.981     0.605    0.739
#> 14 V14      5   500   500  0.047   0.086   . -0.005   0.849     2.959    0.228
#> 15 V15      5   500   500  0.004   0.869     -0.006   0.690     0.168    0.919
#> 16 V16      5   500   500 -0.008   0.789     -0.015   0.264     1.632    0.442
#> 17 V17      5   500   500 -0.004   0.876     -0.006   0.754     0.110    0.946
#> 18 V18      5   500   500 -0.037   0.176      0.009   0.284     2.336    0.311
#> 19 V19      5   500   500  0.033   0.223     -0.022   0.086  .  3.409    0.182
#> 20 V20      5   500   500 -0.028   0.247     -0.022   0.153     3.213    0.201
#> 21 V21      5   500   500 -0.001   0.974     -0.006   0.803     0.075    0.963
#> 22 V22      5   500   500  0.014   0.653     -0.007   0.525     0.976    0.614
#> 23 V23      5   500   500  0.002   0.922     -0.004   0.495     0.567    0.753
#> 24 V24      5   500   500 -0.009   0.624      0.015   0.308     1.322    0.516
#> 25 V25      5   500   500  0.008   0.697     -0.005   0.589     0.294    0.863
#> 26 V26      5   500   500  0.013   0.636      0.002   0.823     1.217    0.544
#> 27 V27      5   500   500  0.007   0.804     -0.003   0.932     0.062    0.969
#> 28 V28      5   500   500  0.010   0.723     -0.004   0.821     0.462    0.793
#> 29 V29      5   500   500 -0.002   0.947      0.005   0.779     0.304    0.859
#> 30 V30      5   500   500  0.019   0.361     -0.004   0.615     1.153    0.562
#> 31 V31      5   500   500 -0.005   0.861     -0.004   0.934     0.044    0.978
#> 32 V32      5   500   500  0.007   0.795     -0.008   0.610     0.286    0.867
#> 33 V33      5   500   500 -0.037   0.085   .  0.018   0.222     3.099    0.212
#> 34 V34      5   500   500 -0.025   0.402     -0.011   0.318     2.026    0.363
#> 35 V35      5   500   500 -0.002   0.936     -0.018   0.086  .  3.279    0.194
#> 36 V36      5   500   500 -0.032   0.156      0.022   0.098  .  2.814    0.245
#> 37 V37      5   500   500 -0.005   0.800      0.002   0.990     0.128    0.938
#> 38 V38      4   500   500 -0.067   0.022   * -0.006   0.717     5.471    0.065
#> 39 V39      5   500   500  0.007   0.801      0.006   0.487     0.777    0.678
#> 40 V40      5   500   500 -0.029   0.310     -0.008   0.605     1.030    0.597
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

# (b)-2 Compute RDIF statistics with purification based on RDIF(S)
dif_puri_s <- rdif(
  x = est_par, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "rdifs"
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 
#> Purification is finished. 
print(dif_puri_s)
#> 
#> Call:
#> rdif.default(x = est_par, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, purify = TRUE, purify.by = "rdifs")
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 3, 4, 11, 14 
#>   - DIF Items identified by RDIF(S): 
#>     5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11, 14, 19 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs    
#> 1   V1   500   500 -0.110   0.000 ***  0.008   0.204    13.750    0.001 ***
#> 2   V2   500   500 -0.175   0.000 *** -0.020   0.163    42.916    0.000 ***
#> 3   V3   500   500 -0.059   0.044   * -0.003   0.688     4.261    0.119    
#> 4   V4   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000 ***
#> 5   V5   500   500 -0.006   0.830      0.021   0.003 **  9.607    0.008  **
#> 6   V6   500   500  0.012   0.534      0.000   0.463     0.579    0.749    
#> 7   V7   500   500  0.018   0.533      0.005   0.470     0.854    0.653    
#> 8   V8   500   500  0.028   0.317      0.002   0.887     1.329    0.514    
#> 9   V9   500   500 -0.012   0.689      0.008   0.499     0.495    0.781    
#> 10 V10   500   500  0.031   0.229      0.011   0.390     1.536    0.464    
#> 11 V11   500   500  0.066   0.023   * -0.008   0.638     8.193    0.017   *
#> 12 V12   500   500  0.002   0.961      0.002   0.613     0.296    0.863    
#> 13 V13   500   500 -0.010   0.737     -0.003   0.750     0.272    0.873    
#> 14 V14   500   500  0.068   0.013   * -0.008   0.427     6.525    0.038   *
#> 15 V15   500   500  0.019   0.488     -0.006   0.404     0.752    0.687    
#> 16 V16   500   500  0.007   0.809     -0.016   0.324     2.394    0.302    
#> 17 V17   500   500  0.011   0.696     -0.008   0.690     1.180    0.554    
#> 18 V18   500   500 -0.018   0.505      0.007   0.608     0.554    0.758    
#> 19 V19   500   500  0.052   0.055   . -0.021   0.034  *  6.103    0.047   *
#> 20 V20   500   500 -0.002   0.936     -0.019   0.113     2.510    0.285    
#> 21 V21   500   500  0.018   0.499     -0.003   0.655     0.520    0.771    
#> 22 V22   500   500  0.028   0.358     -0.005   0.578     1.837    0.399    
#> 23 V23   500   500  0.015   0.414     -0.005   0.172     1.879    0.391    
#> 24 V24   500   500  0.000   0.980      0.014   0.601     0.813    0.666    
#> 25 V25   500   500  0.020   0.350     -0.004   0.318     1.062    0.588    
#> 26 V26   500   500  0.024   0.383      0.002   0.851     1.482    0.476    
#> 27 V27   500   500  0.024   0.385     -0.001   0.753     0.758    0.685    
#> 28 V28   500   500  0.024   0.412     -0.004   0.864     1.477    0.478    
#> 29 V29   500   500  0.006   0.784      0.006   0.991     0.509    0.775    
#> 30 V30   500   500  0.028   0.176     -0.007   0.254     1.839    0.399    
#> 31 V31   500   500  0.014   0.635     -0.002   0.803     0.252    0.881    
#> 32 V32   500   500  0.028   0.270     -0.008   0.238     2.153    0.341    
#> 33 V33   500   500 -0.026   0.217      0.019   0.402     1.643    0.440    
#> 34 V34   500   500 -0.010   0.734     -0.009   0.290     1.392    0.499    
#> 35 V35   500   500  0.011   0.713     -0.017   0.093  .  3.930    0.140    
#> 36 V36   500   500 -0.019   0.395      0.021   0.268     1.227    0.541    
#> 37 V37   500   500  0.007   0.743      0.003   0.614     0.258    0.879    
#> 38 V38   500   500 -0.051   0.078   . -0.008   0.427     3.844    0.146    
#> 39 V39   500   500  0.025   0.366      0.006   0.933     1.022    0.600    
#> 40 V40   500   500 -0.013   0.645     -0.008   0.622     0.305    0.859    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 2
#>   - RDIF statistic used for purification: RDIF(S)
#>   - DIF Items identified by RDIF(S): 
#>     5, 19 
#>   - RDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs
#> 1   V1      2   500   500 -0.109   0.000 ***  0.009   0.181    13.707    0.001
#> 2   V2      2   500   500 -0.173   0.000 *** -0.022   0.147    41.923    0.000
#> 3   V3      2   500   500 -0.058   0.049   * -0.003   0.711     4.062    0.131
#> 4   V4      2   500   500 -0.187   0.000 ***  0.012   0.293    49.400    0.000
#> 5   V5      0   500   500 -0.006   0.830      0.021   0.003 **  9.607    0.008
#> 6   V6      2   500   500  0.014   0.483      0.000   0.471     0.622    0.733
#> 7   V7      2   500   500  0.020   0.506      0.004   0.503     0.836    0.658
#> 8   V8      2   500   500  0.029   0.294      0.002   0.852     1.408    0.494
#> 9   V9      2   500   500 -0.010   0.726      0.009   0.463     0.553    0.758
#> 10 V10      2   500   500  0.033   0.198      0.012   0.309     1.847    0.397
#> 11 V11      2   500   500  0.067   0.021   * -0.007   0.735     7.861    0.020
#> 12 V12      2   500   500  0.002   0.944      0.001   0.608     0.295    0.863
#> 13 V13      2   500   500 -0.009   0.759     -0.002   0.844     0.166    0.920
#> 14 V14      2   500   500  0.070   0.011   * -0.008   0.421     6.910    0.032
#> 15 V15      2   500   500  0.020   0.458     -0.005   0.401     0.791    0.673
#> 16 V16      2   500   500  0.008   0.789     -0.015   0.356     2.212    0.331
#> 17 V17      2   500   500  0.011   0.687     -0.006   0.812     0.780    0.677
#> 18 V18      2   500   500 -0.017   0.547      0.009   0.552     0.553    0.758
#> 19 V19      1   500   500  0.052   0.056   . -0.020   0.040  *  5.878    0.053
#> 20 V20      2   500   500  0.001   0.966     -0.020   0.089  .  2.901    0.234
#> 21 V21      2   500   500  0.020   0.457     -0.003   0.659     0.604    0.739
#> 22 V22      2   500   500  0.028   0.346     -0.004   0.686     1.590    0.452
#> 23 V23      2   500   500  0.017   0.368     -0.007   0.112     2.538    0.281
#> 24 V24      2   500   500  0.002   0.929      0.015   0.595     1.002    0.606
#> 25 V25      2   500   500  0.021   0.320     -0.005   0.267     1.282    0.527
#> 26 V26      2   500   500  0.025   0.364      0.001   0.793     1.414    0.493
#> 27 V27      2   500   500  0.026   0.354     -0.002   0.688     0.860    0.651
#> 28 V28      2   500   500  0.025   0.403     -0.003   0.908     1.389    0.499
#> 29 V29      2   500   500  0.007   0.758      0.005   0.937     0.354    0.838
#> 30 V30      2   500   500  0.029   0.161     -0.008   0.232     1.970    0.373
#> 31 V31      2   500   500  0.015   0.594      0.000   0.885     0.287    0.867
#> 32 V32      2   500   500  0.030   0.233     -0.009   0.176     2.690    0.261
#> 33 V33      2   500   500 -0.025   0.237      0.018   0.449     1.560    0.458
#> 34 V34      2   500   500 -0.009   0.762     -0.009   0.277     1.422    0.491
#> 35 V35      2   500   500  0.012   0.693     -0.018   0.082  .  4.234    0.120
#> 36 V36      2   500   500 -0.018   0.431      0.022   0.272     1.210    0.546
#> 37 V37      2   500   500  0.008   0.693      0.004   0.650     0.215    0.898
#> 38 V38      2   500   500 -0.050   0.086   . -0.008   0.442     3.644    0.162
#> 39 V39      2   500   500  0.027   0.330      0.006   0.963     1.142    0.565
#> 40 V40      2   500   500 -0.013   0.667     -0.008   0.642     0.268    0.875
#>       
#> 1  ***
#> 2  ***
#> 3     
#> 4  ***
#> 5   **
#> 6     
#> 7     
#> 8     
#> 9     
#> 10    
#> 11   *
#> 12    
#> 13    
#> 14   *
#> 15    
#> 16    
#> 17    
#> 18    
#> 19   .
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
#> 38    
#> 39    
#> 40    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 

# (b)-3 Compute RDIF statistics with purification based on RDIF(RS)
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
#>     5, 19 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 4, 5, 11, 14, 19 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs    
#> 1   V1   500   500 -0.110   0.000 ***  0.008   0.204    13.750    0.001 ***
#> 2   V2   500   500 -0.175   0.000 *** -0.020   0.163    42.916    0.000 ***
#> 3   V3   500   500 -0.059   0.044   * -0.003   0.688     4.261    0.119    
#> 4   V4   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000 ***
#> 5   V5   500   500 -0.006   0.830      0.021   0.003 **  9.607    0.008  **
#> 6   V6   500   500  0.012   0.534      0.000   0.463     0.579    0.749    
#> 7   V7   500   500  0.018   0.533      0.005   0.470     0.854    0.653    
#> 8   V8   500   500  0.028   0.317      0.002   0.887     1.329    0.514    
#> 9   V9   500   500 -0.012   0.689      0.008   0.499     0.495    0.781    
#> 10 V10   500   500  0.031   0.229      0.011   0.390     1.536    0.464    
#> 11 V11   500   500  0.066   0.023   * -0.008   0.638     8.193    0.017   *
#> 12 V12   500   500  0.002   0.961      0.002   0.613     0.296    0.863    
#> 13 V13   500   500 -0.010   0.737     -0.003   0.750     0.272    0.873    
#> 14 V14   500   500  0.068   0.013   * -0.008   0.427     6.525    0.038   *
#> 15 V15   500   500  0.019   0.488     -0.006   0.404     0.752    0.687    
#> 16 V16   500   500  0.007   0.809     -0.016   0.324     2.394    0.302    
#> 17 V17   500   500  0.011   0.696     -0.008   0.690     1.180    0.554    
#> 18 V18   500   500 -0.018   0.505      0.007   0.608     0.554    0.758    
#> 19 V19   500   500  0.052   0.055   . -0.021   0.034  *  6.103    0.047   *
#> 20 V20   500   500 -0.002   0.936     -0.019   0.113     2.510    0.285    
#> 21 V21   500   500  0.018   0.499     -0.003   0.655     0.520    0.771    
#> 22 V22   500   500  0.028   0.358     -0.005   0.578     1.837    0.399    
#> 23 V23   500   500  0.015   0.414     -0.005   0.172     1.879    0.391    
#> 24 V24   500   500  0.000   0.980      0.014   0.601     0.813    0.666    
#> 25 V25   500   500  0.020   0.350     -0.004   0.318     1.062    0.588    
#> 26 V26   500   500  0.024   0.383      0.002   0.851     1.482    0.476    
#> 27 V27   500   500  0.024   0.385     -0.001   0.753     0.758    0.685    
#> 28 V28   500   500  0.024   0.412     -0.004   0.864     1.477    0.478    
#> 29 V29   500   500  0.006   0.784      0.006   0.991     0.509    0.775    
#> 30 V30   500   500  0.028   0.176     -0.007   0.254     1.839    0.399    
#> 31 V31   500   500  0.014   0.635     -0.002   0.803     0.252    0.881    
#> 32 V32   500   500  0.028   0.270     -0.008   0.238     2.153    0.341    
#> 33 V33   500   500 -0.026   0.217      0.019   0.402     1.643    0.440    
#> 34 V34   500   500 -0.010   0.734     -0.009   0.290     1.392    0.499    
#> 35 V35   500   500  0.011   0.713     -0.017   0.093  .  3.930    0.140    
#> 36 V36   500   500 -0.019   0.395      0.021   0.268     1.227    0.541    
#> 37 V37   500   500  0.007   0.743      0.003   0.614     0.258    0.879    
#> 38 V38   500   500 -0.051   0.078   . -0.008   0.427     3.844    0.146    
#> 39 V39   500   500  0.025   0.366      0.006   0.933     1.022    0.600    
#> 40 V40   500   500 -0.013   0.645     -0.008   0.622     0.305    0.859    
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
#>     id n.iter n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs    rdifrs p.rdifrs
#> 1   V1      2   500   500 -0.118   0.000 ***  0.009   0.123    16.190    0.000
#> 2   V2      1   500   500 -0.184   0.000 *** -0.020   0.162    47.376    0.000
#> 3   V3      4   500   500 -0.074   0.012   * -0.002   0.937     6.296    0.043
#> 4   V4      0   500   500 -0.189   0.000 ***  0.013   0.255    50.605    0.000
#> 5   V5      3   500   500 -0.016   0.593      0.021   0.003 ** 10.940    0.004
#> 6   V6      5   500   500  0.001   0.954      0.000   0.867     0.032    0.984
#> 7   V7      5   500   500  0.002   0.931      0.006   0.300     1.075    0.584
#> 8   V8      5   500   500  0.014   0.610      0.001   0.833     0.674    0.714
#> 9   V9      5   500   500 -0.027   0.342      0.009   0.240     1.754    0.416
#> 10 V10      5   500   500  0.008   0.762      0.009   0.428     0.638    0.727
#> 11 V11      5   500   500  0.051   0.080   . -0.006   0.761     4.676    0.096
#> 12 V12      5   500   500 -0.010   0.742      0.001   0.657     0.542    0.762
#> 13 V13      5   500   500 -0.022   0.463     -0.002   0.989     0.570    0.752
#> 14 V14      5   500   500  0.048   0.079   . -0.006   0.852     3.082    0.214
#> 15 V15      5   500   500  0.006   0.837     -0.006   0.686     0.166    0.920
#> 16 V16      5   500   500 -0.007   0.793     -0.016   0.276     1.547    0.461
#> 17 V17      5   500   500 -0.005   0.868     -0.006   0.808     0.059    0.971
#> 18 V18      5   500   500 -0.036   0.190      0.010   0.244     2.389    0.303
#> 19 V19      5   500   500  0.034   0.206     -0.022   0.088  .  3.442    0.179
#> 20 V20      5   500   500 -0.027   0.273     -0.020   0.205     2.669    0.263
#> 21 V21      5   500   500  0.000   0.990     -0.006   0.824     0.053    0.974
#> 22 V22      5   500   500  0.013   0.653     -0.007   0.614     0.743    0.690
#> 23 V23      5   500   500  0.004   0.836     -0.004   0.465     0.582    0.748
#> 24 V24      5   500   500 -0.007   0.698      0.014   0.348     1.234    0.539
#> 25 V25      5   500   500  0.010   0.639     -0.005   0.547     0.362    0.834
#> 26 V26      5   500   500  0.014   0.614      0.003   0.794     1.471    0.479
#> 27 V27      5   500   500  0.008   0.771     -0.002   0.987     0.108    0.948
#> 28 V28      5   500   500  0.010   0.728     -0.004   0.834     0.425    0.808
#> 29 V29      5   500   500  0.000   0.983      0.005   0.808     0.316    0.854
#> 30 V30      5   500   500  0.020   0.324     -0.005   0.573     1.299    0.522
#> 31 V31      5   500   500 -0.004   0.886     -0.003   0.977     0.020    0.990
#> 32 V32      5   500   500  0.008   0.753     -0.008   0.621     0.292    0.864
#> 33 V33      5   500   500 -0.035   0.099   .  0.017   0.250     2.876    0.238
#> 34 V34      5   500   500 -0.024   0.413     -0.011   0.348     1.851    0.396
#> 35 V35      5   500   500 -0.002   0.946     -0.018   0.089  .  3.226    0.199
#> 36 V36      5   500   500 -0.030   0.179      0.021   0.112     2.585    0.274
#> 37 V37      5   500   500 -0.004   0.865      0.002   0.986     0.071    0.965
#> 38 V38      5   500   500 -0.068   0.020   * -0.007   0.726     5.567    0.062
#> 39 V39      5   500   500  0.008   0.782      0.005   0.530     0.688    0.709
#> 40 V40      5   500   500 -0.029   0.313     -0.009   0.609     1.017    0.602
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

# }
```
