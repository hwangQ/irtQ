# Generalized IRT residual-based DIF detection framework for multiple groups (GRDIF)

This function computes three GRDIF statistics, \\GRDIF\_{R}\\,
\\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\, for analyzing differential item
functioning (DIF) among multiple groups (Lim et al., 2024). They are
specialized to capture uniform DIF, nonuniform DIF, and mixed DIF,
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
  range = c(-4, 4),
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
  range = c(-4, 4),
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
  range = c(-4, 4),
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
  values). If not provided, `grdif()` will estimate ability parameters
  internally before computing the GRDIF statistics. See
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  for more information on scoring methods. Default is `NULL`.

- group:

  A numeric or character vector indicating examinees' group membership.
  The length of the vector must match the number of rows in the response
  data matrix.

- focal.name:

  A numeric or character vector specifying the levels associated with
  the focal groups. For example, consider
  `group = c(0, 0, 1, 2, 2, 3, 3)`, where '1', '2', and '3' indicate
  three distinct focal groups and '0' represents the reference group. In
  this case, set `focal.name = c(1, 2, 3)`.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for
  hypothesis testing using the GRDIF statistics. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying which GRDIF statistic is used to perform
  the purification. Available options are "grdifrs" for \\GRDIF\_{RS}\\,
  "grdifr" for \\GRDIF\_{R}\\, and "grdifs" for \\GRDIF\_{S}\\.

- max.iter:

  A positive integer specifying the maximum number of iterations allowed
  for the purification process. Default is `10`. If the limit is reached
  while flagged items remain, a warning is issued and the `complete`
  element of the purification results is `FALSE`; for tests with many
  items or many DIF items, consider increasing `max.iter`.

- min.resp:

  A positive integer specifying the minimum number of valid item
  responses required from an examinee in order to compute an ability
  estimate. Default is `NULL`. See **Details** for more information.

- post.hoc:

  A logical value indicating whether to perform post-hoc RDIF analyses
  for all possible pairwise group comparisons on items flagged as
  statistically significant. Default is TRUE. See details below.

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

  A list of sub-objects presenting the results of DIF analysis without a
  purification procedure. These include:

  dif_stat

  :   A data frame summarizing the results of the three GRDIF statistics
      for all evaluated items. Columns include item ID, \\GRDIF\_{R}\\,
      \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\ statistics; their
      corresponding p-values; the sample size of the reference group;
      the sample sizes of the focal groups; and the total sample size.

  moments

  :   A list of three data frames reporting the moments of mean raw
      residuals (MRRs) and mean squared residuals (MSRs) across all
      compared groups. The first contains means, the second variances,
      and the third covariances of MRRs and MSRs.

  dif_item

  :   A list of three numeric vectors indicating the items flagged as
      potential DIF items by each of the GRDIF statistics
      (\\GRDIF\_{R}\\, \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\).

  score

  :   A numeric vector of ability estimates used to compute the GRDIF
      statistics.

  post.hoc

  :   A list of three data frames containing post-hoc RDIF analysis
      results for all possible pairwise group comparisons. Each data
      frame corresponds to the results for items flagged by
      \\GRDIF\_{R}\\, \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\, respectively.

- purify:

  A logical value indicating whether the purification process was
  applied.

- with_purify:

  A list of sub-objects presenting the results of DIF analysis with a
  purification procedure. These include:

  purify.by

  :   A character string indicating which GRDIF statistic was used for
      purification: "grdifr", "grdifs", or "grdifrs", corresponding to
      \\GRDIF\_{R}\\, \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\, respectively.

  dif_stat

  :   A data frame summarizing the GRDIF results across iterations.
      Columns include item ID, the three GRDIF statistics and their
      p-values, sample size of the reference group, sample sizes of the
      focal groups, total sample size, and the iteration number at which
      each statistic was computed.

  moments

  :   A list of three data frames showing the MRR and MSR moments across
      iterations. The final column in each data frame indicates the
      iteration in which the statistics were computed.

  dif_item

  :   A numeric vector of the positions (rows of `x`) of the items
      flagged as DIF by the `purify.by` statistic across all
      purification iterations, sorted in ascending order.

  n.iter

  :   The total number of iterations executed during the purification
      process.

  score

  :   A numeric vector of the final purified ability estimates used for
      computing GRDIF statistics.

  post.hoc

  :   A data frame containing the post-hoc RDIF analysis results for
      flagged items across all possible pairwise group comparisons,
      updated at each iteration.

  complete

  :   A logical value indicating whether the purification process was
      completed. If `FALSE`, the process reached the maximum number of
      iterations without convergence.

- alpha:

  The significance level (\\\alpha\\) used for hypothesis testing of the
  GRDIF statistics.

## Details

The GRDIF framework (Lim et al., 2024) is a generalized version of the
RDIF detection framework, designed to assess DIF across multiple groups.
The framework includes three statistics: \\GRDIF\_{R}\\, \\GRDIF\_{S}\\,
and \\GRDIF\_{RS}\\, which are tailored to detect uniform, nonuniform,
and mixed DIF, respectively. Under the null hypothesis that the test
contains no DIF items, the statistics \\GRDIF\_{R}\\, \\GRDIF\_{S}\\,
and \\GRDIF\_{RS}\\ asymptotically follow \\\chi^{2}\\ distributions
with G-1, G-1, and 2(G-1) degrees of freedom, respectively, where *G*
represents the number of groups being compared. For more information on
the GRDIF framework, see Lim et al. (2024).

The `grdif()` function computes all three GRDIF statistics:
\\GRDIF\_{R}\\, \\GRDIF\_{S}\\, and \\GRDIF\_{RS}\\. It supports both
dichotomous and polytomous item response data. To compute these
statistics, `grdif()` requires: (1) item parameter estimates obtained
from aggregate data (regardless of group membership); (2) examinees'
ability estimates (e.g., ML); and (3) their item response data. Note
that ability estimates must be computed using the item parameters
estimated from the aggregate data. The item parameters should be
provided via the `x` argument, the ability estimates via the `score`
argument, and the response data via the `data` argument. If ability
estimates are not supplied (`score = NULL`), `grdif()` automatically
computes them using the scoring method specified in the `method`
argument (e.g., `method = "ML"`).

The `group` argument accepts a vector of numeric or character values,
indicating the group membership of examinees. The vector may contain
multiple distinct values, where one represents the reference group and
the others represent focal groups. Its length must match the number of
rows in the response data, with each value corresponding to an
examinee's group membership. Once `group` is specified, a numeric or
character vector must be supplied via the `focal.name` argument to
define which group(s) in `group` represent the focal groups. The
reference group is defined as the group not included in `focal.name`.

Similar to the original RDIF framework for two-group comparisons, the
GRDIF framework supports an iterative purification process. When
`purify = TRUE`, purification is conducted based on the GRDIF statistic
specified in the `purify.by` argument (e.g., `purify.by = "grdifrs"`).
During each iteration, examinees' latent abilities are re-estimated
using only the purified items, with the scoring method determined by the
`method` argument. The process continues until no additional DIF items
are flagged or until the number of iterations reaches the specified
`max.iter` limit. For details on the purification procedure, see Lim et
al. (2022).

Scoring based on a limited number of items can lead to large standard
errors, which may compromise the effectiveness of DIF detection within
the GRDIF framework. The `min.resp` argument can be used to exclude
ability estimates with substantial standard errors, especially during
the purification process. For example, if `min.resp` is not NULL (e.g.,
`min.resp = 5`), examinees whose total number of responses falls below
the specified threshold will have their responses treated as missing
values (i.e., NA). Consequently, their ability estimates will also be
missing and will not be used when computing the GRDIF statistics. If
`min.resp = NULL`, an examinee's score will be computed as long as at
least one response is available.

The `post.hoc` argument enables post-hoc RDIF analyses across all
possible pairwise group comparisons for items flagged as statistically
significant. For instance, consider four groups of examinees: A, B, C,
and D. If `post.hoc = TRUE`, the `grdif()` function will perform
pairwise RDIF analyses for each flagged item across all group pairs
(A-B, A-C, A-D, B-C, B-D, and C-D). This provides a more detailed
understanding of which specific group pairs exhibit DIF. Note that when
purification is enabled (i.e., `purify = TRUE`), post-hoc RDIF analyses
are conducted for each flagged item at each iteration of the
purification process.

## Methods (by class)

- `grdif(default)`: Default method to compute the three GRDIF statistics
  for multiple-group data using a data frame `x` that contains item
  metadata.

- `grdif(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `grdif(est_item)`: An object created by the function
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

Lim, H., Zhu, D., Choe, E. M., & Han, K. T. (2024). Detecting
differential item functioning among multiple groups using IRT residual
DIF framework. *Journal of Educational Measurement, 61*(4), 656-681.

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md)
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
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = .15),
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
# Therefore, the first four items contain uniform DIF across all focal groups
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
#>  EM iteration: 1, Loglike: -48577.1497, Max-Change: 1.427193 EM iteration: 2, Loglike: -43407.2535, Max-Change: 0.375334 EM iteration: 3, Loglike: -43351.5535, Max-Change: 0.185529 EM iteration: 4, Loglike: -43334.1667, Max-Change: 0.107977 EM iteration: 5, Loglike: -43325.2522, Max-Change: 0.075754 EM iteration: 6, Loglike: -43320.2571, Max-Change: 0.056555 EM iteration: 7, Loglike: -43317.3671, Max-Change: 0.043296 EM iteration: 8, Loglike: -43315.6636, Max-Change: 0.03352 EM iteration: 9, Loglike: -43314.6455, Max-Change: 0.026079 EM iteration: 10, Loglike: -43314.0296, Max-Change: 0.020332 EM iteration: 11, Loglike: -43313.6522, Max-Change: 0.015866 EM iteration: 12, Loglike: -43313.4177, Max-Change: 0.012384 EM iteration: 13, Loglike: -43313.2695, Max-Change: 0.009667 EM iteration: 14, Loglike: -43313.1737, Max-Change: 0.007545 EM iteration: 15, Loglike: -43313.1103, Max-Change: 0.005888 EM iteration: 16, Loglike: -43313.0669, Max-Change: 0.004593 EM iteration: 17, Loglike: -43313.0361, Max-Change: 0.003582 EM iteration: 18, Loglike: -43313.0135, Max-Change: 0.002793 EM iteration: 19, Loglike: -43312.9962, Max-Change: 0.002176 EM iteration: 20, Loglike: -43312.9826, Max-Change: 0.001694 EM iteration: 21, Loglike: -43312.9715, Max-Change: 0.001318 EM iteration: 22, Loglike: -43312.9623, Max-Change: 0.001025 EM iteration: 23, Loglike: -43312.9545, Max-Change: 0.000804 EM iteration: 24, Loglike: -43312.9478, Max-Change: 0.00065 EM iteration: 25, Loglike: -43312.9419, Max-Change: 0.00057 EM iteration: 26, Loglike: -43312.9368, Max-Change: 0.000502 EM iteration: 27, Loglike: -43312.9323, Max-Change: 0.000445 EM iteration: 28, Loglike: -43312.9283, Max-Change: 0.000399 EM iteration: 29, Loglike: -43312.9248, Max-Change: 0.000361 EM iteration: 30, Loglike: -43312.9217, Max-Change: 0.000327 EM iteration: 31, Loglike: -43312.9190, Max-Change: 0.000297 EM iteration: 32, Loglike: -43312.9165, Max-Change: 0.000271 EM iteration: 33, Loglike: -43312.9143, Max-Change: 0.000247 EM iteration: 34, Loglike: -43312.9124, Max-Change: 0.000226 EM iteration: 35, Loglike: -43312.9107, Max-Change: 0.000206 EM iteration: 36, Loglike: -43312.9091, Max-Change: 0.000189 EM iteration: 37, Loglike: -43312.9077, Max-Change: 0.000173 EM iteration: 38, Loglike: -43312.9065, Max-Change: 0.00016 EM iteration: 39, Loglike: -43312.9054, Max-Change: 0.000147 EM iteration: 40, Loglike: -43312.9044, Max-Change: 0.000136 EM iteration: 41, Loglike: -43312.9035, Max-Change: 0.000126 EM iteration: 42, Loglike: -43312.9027, Max-Change: 0.000116 EM iteration: 43, Loglike: -43312.9020, Max-Change: 0.000107 EM iteration: 44, Loglike: -43312.9013, Max-Change: 9.9e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 2.91 seconds. 
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
#> 1   V1   500    500    500    500 13.154    0.004  **  8.058    0.045   *
#> 2   V2   500    500    500    500 96.369    0.000 ***  2.000    0.572    
#> 3   V3   500    500    500    500 29.873    0.000 *** 19.246    0.000 ***
#> 4   V4   500    500    500    500 77.456    0.000 *** 41.597    0.000 ***
#> 5   V5   500    500    500    500  5.607    0.132      3.636    0.304    
#> 6   V6   500    500    500    500  1.658    0.646      0.928    0.819    
#> 7   V7   500    500    500    500  0.739    0.864      2.728    0.436    
#> 8   V8   500    500    500    500  1.220    0.748      4.916    0.178    
#> 9   V9   500    500    500    500  0.672    0.880      0.320    0.956    
#> 10 V10   500    500    500    500  1.779    0.620      5.262    0.154    
#> 11 V11   500    500    500    500  3.533    0.317      0.958    0.811    
#> 12 V12   500    500    500    500  1.025    0.795      0.630    0.889    
#> 13 V13   500    500    500    500  1.001    0.801      2.942    0.401    
#> 14 V14   500    500    500    500  5.495    0.139      5.752    0.124    
#> 15 V15   500    500    500    500  4.650    0.199      5.189    0.158    
#> 16 V16   500    500    500    500  1.504    0.681      0.346    0.951    
#> 17 V17   500    500    500    500  3.102    0.376      4.323    0.229    
#> 18 V18   500    500    500    500  2.234    0.525      2.210    0.530    
#> 19 V19   500    500    500    500  1.336    0.721      0.226    0.973    
#> 20 V20   500    500    500    500 16.500    0.001 ***  0.900    0.825    
#> 21 V21   500    500    500    500  5.680    0.128      6.536    0.088   .
#> 22 V22   500    500    500    500  7.513    0.057   .  6.500    0.090   .
#> 23 V23   500    500    500    500  0.392    0.942      4.418    0.220    
#> 24 V24   500    500    500    500  5.900    0.117      6.070    0.108    
#> 25 V25   500    500    500    500  4.882    0.181      4.858    0.182    
#> 26 V26   500    500    500    500  1.038    0.792      1.440    0.696    
#> 27 V27   500    500    500    500  1.213    0.750      1.710    0.635    
#> 28 V28   500    500    500    500  7.098    0.069   .  1.724    0.631    
#> 29 V29   500    500    500    500  0.735    0.865      0.500    0.919    
#> 30 V30   500    500    500    500  1.216    0.749      2.597    0.458    
#> 31 V31   500    500    500    500  1.558    0.669      4.717    0.194    
#> 32 V32   500    500    500    500  2.913    0.405      1.503    0.681    
#> 33 V33   500    500    500    500  2.173    0.537      5.840    0.120    
#> 34 V34   500    500    500    500  3.954    0.266      7.293    0.063   .
#> 35 V35   500    500    500    500  1.245    0.742      2.402    0.493    
#> 36 V36   500    500    500    500  3.275    0.351      0.668    0.881    
#> 37 V37   500    500    500    500  0.358    0.949      0.572    0.903    
#> 38 V38   500    500    500    500  6.205    0.102      1.738    0.628    
#> 39 V39   500    500    500    500  2.417    0.490      1.414    0.702    
#> 40 V40   500    500    500    500  1.538    0.674      4.804    0.187    
#>    grdifrs p.grdifrs    
#> 1    15.76     0.015   *
#> 2   102.89     0.000 ***
#> 3    36.03     0.000 ***
#> 4    83.28     0.000 ***
#> 5     6.24     0.396    
#> 6     1.77     0.939    
#> 7     3.43     0.753    
#> 8     5.51     0.480    
#> 9     1.17     0.979    
#> 10    5.81     0.444    
#> 11    3.81     0.702    
#> 12    1.97     0.922    
#> 13    3.31     0.769    
#> 14   15.32     0.018   *
#> 15    6.29     0.392    
#> 16    3.57     0.735    
#> 17    6.35     0.385    
#> 18    2.96     0.814    
#> 19    1.60     0.953    
#> 20   17.55     0.007  **
#> 21    9.66     0.140    
#> 22   11.12     0.085   .
#> 23    5.10     0.531    
#> 24    9.45     0.150    
#> 25   10.87     0.093   .
#> 26    4.17     0.654    
#> 27    3.64     0.725    
#> 28    9.28     0.158    
#> 29    1.24     0.975    
#> 30    6.31     0.390    
#> 31    7.46     0.281    
#> 32    5.48     0.484    
#> 33    8.61     0.196    
#> 34    8.55     0.200    
#> 35    2.59     0.859    
#> 36    4.09     0.664    
#> 37    1.42     0.965    
#> 38   10.45     0.107    
#> 39    5.48     0.484    
#> 40    8.63     0.195    
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
#> 1   V1      0 & 1 -0.0571 -1.9787 -0.0143 -1.5392  5.4726  0.0478  0.1238
#> 2   V1      0 & 2 -0.0453 -1.6508 -0.0214  1.4492  3.2189  0.0988  0.1473
#> 3   V1      0 & 3  0.0323  1.1415 -0.0163 -0.3473  1.3143  0.2536  0.7284
#> 4   V1      1 & 2  0.0118  0.4281 -0.0072  2.8172  8.4417  0.6686  0.0048
#> 5   V1      1 & 3  0.0894  3.1484 -0.0020  1.1352 11.0095  0.0016  0.2563
#> 6   V1      2 & 3  0.0776  2.8783  0.0051 -1.7120  8.2949  0.0040  0.0869
#> 7   V2      0 & 1 -0.1060 -4.5351 -0.0494 -1.3017 20.6424  0.0000  0.1930
#> 8   V2      0 & 2  0.0013  0.0559 -0.0332 -0.2033  0.0417  0.9554  0.8389
#> 9   V2      0 & 3  0.1196  4.9024 -0.0171 -0.6807 24.1845  0.0000  0.4961
#> 10  V2      1 & 2  0.1074  4.8742  0.0162  1.1210 24.0713  0.0000  0.2623
#> 11  V2      1 & 3  0.2256  9.8131  0.0323  0.6452 98.9891  0.0000  0.5188
#> 12  V2      2 & 3  0.1183  5.1209  0.0161 -0.4860 27.6335  0.0000  0.6270
#> 13  V3      0 & 1 -0.0269 -0.9246 -0.0225 -2.2768  5.7699  0.3552  0.0228
#> 14  V3      0 & 2 -0.0910 -3.2646 -0.0169  1.5449 10.8342  0.0011  0.1224
#> 15  V3      0 & 3  0.0565  1.9705 -0.0372 -2.4747  6.9738  0.0488  0.0133
#> 16  V3      1 & 2 -0.0641 -2.2990  0.0055  3.5583 15.5463  0.0215  0.0004
#> 17  V3      1 & 3  0.0834  2.9071 -0.0147 -0.2808  8.4681  0.0036  0.7789
#> 18  V3      2 & 3  0.1476  5.3711 -0.0202 -3.7051 28.8580  0.0000  0.0002
#> 19  V4      0 & 1 -0.0349 -1.3083  0.0130  0.4067  1.9900  0.1908  0.6842
#> 20  V4      0 & 2 -0.1382 -5.6951  0.0132  5.4802 41.4631  0.0000  0.0000
#> 21  V4      0 & 3  0.0610  2.3838 -0.0194  0.0301  6.8637  0.0171  0.9760
#> 22  V4      1 & 2 -0.1034 -4.2016  0.0002  4.9559 39.1282  0.0000  0.0000
#> 23  V4      1 & 3  0.0958  3.7010 -0.0324 -0.3708 13.8353  0.0002  0.7108
#> 24  V4      2 & 3  0.1992  8.4964 -0.0326 -5.3629 72.2070  0.0000  0.0000
#> 25 V20      0 & 1 -0.0327 -1.3665 -0.0125  0.1019  2.1325  0.1718  0.9188
#> 26 V20      0 & 2  0.0143  0.6303 -0.0520 -0.7158  0.7247  0.5285  0.4741
#> 27 V20      0 & 3 -0.0694 -2.9073 -0.0151  0.0406  8.6031  0.0036  0.9676
#> 28 V20      1 & 2  0.0470  2.1321 -0.0395 -0.7806  5.6204  0.0330  0.4350
#> 29 V20      1 & 3 -0.0367 -1.5751 -0.0026 -0.0642  2.5453  0.1152  0.9488
#> 30 V20      2 & 3 -0.0836 -3.8099  0.0368  0.7635 15.0786  0.0001  0.4452
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0648   500   500    1000
#> 2    0.2000   500   500    1000
#> 3    0.5183   500   500    1000
#> 4    0.0147   500   500    1000
#> 5    0.0041   500   500    1000
#> 6    0.0158   500   500    1000
#> 7    0.0000   500   500    1000
#> 8    0.9794   500   500    1000
#> 9    0.0000   500   500    1000
#> 10   0.0000   500   500    1000
#> 11   0.0000   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.0559   500   500    1000
#> 14   0.0044   500   500    1000
#> 15   0.0306   500   500    1000
#> 16   0.0004   500   500    1000
#> 17   0.0145   500   500    1000
#> 18   0.0000   500   500    1000
#> 19   0.3697   500   500    1000
#> 20   0.0000   500   500    1000
#> 21   0.0323   500   500    1000
#> 22   0.0000   500   500    1000
#> 23   0.0010   500   500    1000
#> 24   0.0000   500   500    1000
#> 25   0.3443   500   500    1000
#> 26   0.6960   500   500    1000
#> 27   0.0135   500   500    1000
#> 28   0.0602   500   500    1000
#> 29   0.2801   500   500    1000
#> 30   0.0005   500   500    1000
#> 
#> $by.grdifs
#>    id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1  V1      0 & 1 -0.0571 -1.9787 -0.0143 -1.5392  5.4726  0.0478  0.1238
#> 2  V1      0 & 2 -0.0453 -1.6508 -0.0214  1.4492  3.2189  0.0988  0.1473
#> 3  V1      0 & 3  0.0323  1.1415 -0.0163 -0.3473  1.3143  0.2536  0.7284
#> 4  V1      1 & 2  0.0118  0.4281 -0.0072  2.8172  8.4417  0.6686  0.0048
#> 5  V1      1 & 3  0.0894  3.1484 -0.0020  1.1352 11.0095  0.0016  0.2563
#> 6  V1      2 & 3  0.0776  2.8783  0.0051 -1.7120  8.2949  0.0040  0.0869
#> 7  V3      0 & 1 -0.0269 -0.9246 -0.0225 -2.2768  5.7699  0.3552  0.0228
#> 8  V3      0 & 2 -0.0910 -3.2646 -0.0169  1.5449 10.8342  0.0011  0.1224
#> 9  V3      0 & 3  0.0565  1.9705 -0.0372 -2.4747  6.9738  0.0488  0.0133
#> 10 V3      1 & 2 -0.0641 -2.2990  0.0055  3.5583 15.5463  0.0215  0.0004
#> 11 V3      1 & 3  0.0834  2.9071 -0.0147 -0.2808  8.4681  0.0036  0.7789
#> 12 V3      2 & 3  0.1476  5.3711 -0.0202 -3.7051 28.8580  0.0000  0.0002
#> 13 V4      0 & 1 -0.0349 -1.3083  0.0130  0.4067  1.9900  0.1908  0.6842
#> 14 V4      0 & 2 -0.1382 -5.6951  0.0132  5.4802 41.4631  0.0000  0.0000
#> 15 V4      0 & 3  0.0610  2.3838 -0.0194  0.0301  6.8637  0.0171  0.9760
#> 16 V4      1 & 2 -0.1034 -4.2016  0.0002  4.9559 39.1282  0.0000  0.0000
#> 17 V4      1 & 3  0.0958  3.7010 -0.0324 -0.3708 13.8353  0.0002  0.7108
#> 18 V4      2 & 3  0.1992  8.4964 -0.0326 -5.3629 72.2070  0.0000  0.0000
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0648   500   500    1000
#> 2    0.2000   500   500    1000
#> 3    0.5183   500   500    1000
#> 4    0.0147   500   500    1000
#> 5    0.0041   500   500    1000
#> 6    0.0158   500   500    1000
#> 7    0.0559   500   500    1000
#> 8    0.0044   500   500    1000
#> 9    0.0306   500   500    1000
#> 10   0.0004   500   500    1000
#> 11   0.0145   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.3697   500   500    1000
#> 14   0.0000   500   500    1000
#> 15   0.0323   500   500    1000
#> 16   0.0000   500   500    1000
#> 17   0.0010   500   500    1000
#> 18   0.0000   500   500    1000
#> 
#> $by.grdifrs
#>     id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1   V1      0 & 1 -0.0571 -1.9787 -0.0143 -1.5392  5.4726  0.0478  0.1238
#> 2   V1      0 & 2 -0.0453 -1.6508 -0.0214  1.4492  3.2189  0.0988  0.1473
#> 3   V1      0 & 3  0.0323  1.1415 -0.0163 -0.3473  1.3143  0.2536  0.7284
#> 4   V1      1 & 2  0.0118  0.4281 -0.0072  2.8172  8.4417  0.6686  0.0048
#> 5   V1      1 & 3  0.0894  3.1484 -0.0020  1.1352 11.0095  0.0016  0.2563
#> 6   V1      2 & 3  0.0776  2.8783  0.0051 -1.7120  8.2949  0.0040  0.0869
#> 7   V2      0 & 1 -0.1060 -4.5351 -0.0494 -1.3017 20.6424  0.0000  0.1930
#> 8   V2      0 & 2  0.0013  0.0559 -0.0332 -0.2033  0.0417  0.9554  0.8389
#> 9   V2      0 & 3  0.1196  4.9024 -0.0171 -0.6807 24.1845  0.0000  0.4961
#> 10  V2      1 & 2  0.1074  4.8742  0.0162  1.1210 24.0713  0.0000  0.2623
#> 11  V2      1 & 3  0.2256  9.8131  0.0323  0.6452 98.9891  0.0000  0.5188
#> 12  V2      2 & 3  0.1183  5.1209  0.0161 -0.4860 27.6335  0.0000  0.6270
#> 13  V3      0 & 1 -0.0269 -0.9246 -0.0225 -2.2768  5.7699  0.3552  0.0228
#> 14  V3      0 & 2 -0.0910 -3.2646 -0.0169  1.5449 10.8342  0.0011  0.1224
#> 15  V3      0 & 3  0.0565  1.9705 -0.0372 -2.4747  6.9738  0.0488  0.0133
#> 16  V3      1 & 2 -0.0641 -2.2990  0.0055  3.5583 15.5463  0.0215  0.0004
#> 17  V3      1 & 3  0.0834  2.9071 -0.0147 -0.2808  8.4681  0.0036  0.7789
#> 18  V3      2 & 3  0.1476  5.3711 -0.0202 -3.7051 28.8580  0.0000  0.0002
#> 19  V4      0 & 1 -0.0349 -1.3083  0.0130  0.4067  1.9900  0.1908  0.6842
#> 20  V4      0 & 2 -0.1382 -5.6951  0.0132  5.4802 41.4631  0.0000  0.0000
#> 21  V4      0 & 3  0.0610  2.3838 -0.0194  0.0301  6.8637  0.0171  0.9760
#> 22  V4      1 & 2 -0.1034 -4.2016  0.0002  4.9559 39.1282  0.0000  0.0000
#> 23  V4      1 & 3  0.0958  3.7010 -0.0324 -0.3708 13.8353  0.0002  0.7108
#> 24  V4      2 & 3  0.1992  8.4964 -0.0326 -5.3629 72.2070  0.0000  0.0000
#> 25 V14      0 & 1 -0.0431 -1.5723 -0.0106 -1.3392  3.2138  0.1159  0.1805
#> 26 V14      0 & 2 -0.0129 -0.5099 -0.0428  0.4440  0.3512  0.6101  0.6570
#> 27 V14      0 & 3 -0.0537 -2.0220 -0.0348 -1.5860  7.9560  0.0432  0.1127
#> 28 V14      1 & 2  0.0301  1.1775 -0.0322  1.7309  3.9750  0.2390  0.0835
#> 29 V14      1 & 3 -0.0107 -0.3982 -0.0242 -0.1730  0.1661  0.6905  0.8626
#> 30 V14      2 & 3 -0.0408 -1.6510  0.0079 -1.9889 11.9210  0.0987  0.0467
#> 31 V20      0 & 1 -0.0327 -1.3665 -0.0125  0.1019  2.1325  0.1718  0.9188
#> 32 V20      0 & 2  0.0143  0.6303 -0.0520 -0.7158  0.7247  0.5285  0.4741
#> 33 V20      0 & 3 -0.0694 -2.9073 -0.0151  0.0406  8.6031  0.0036  0.9676
#> 34 V20      1 & 2  0.0470  2.1321 -0.0395 -0.7806  5.6204  0.0330  0.4350
#> 35 V20      1 & 3 -0.0367 -1.5751 -0.0026 -0.0642  2.5453  0.1152  0.9488
#> 36 V20      2 & 3 -0.0836 -3.8099  0.0368  0.7635 15.0786  0.0001  0.4452
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0648   500   500    1000
#> 2    0.2000   500   500    1000
#> 3    0.5183   500   500    1000
#> 4    0.0147   500   500    1000
#> 5    0.0041   500   500    1000
#> 6    0.0158   500   500    1000
#> 7    0.0000   500   500    1000
#> 8    0.9794   500   500    1000
#> 9    0.0000   500   500    1000
#> 10   0.0000   500   500    1000
#> 11   0.0000   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.0559   500   500    1000
#> 14   0.0044   500   500    1000
#> 15   0.0306   500   500    1000
#> 16   0.0004   500   500    1000
#> 17   0.0145   500   500    1000
#> 18   0.0000   500   500    1000
#> 19   0.3697   500   500    1000
#> 20   0.0000   500   500    1000
#> 21   0.0323   500   500    1000
#> 22   0.0000   500   500    1000
#> 23   0.0010   500   500    1000
#> 24   0.0000   500   500    1000
#> 25   0.2005   500   500    1000
#> 26   0.8390   500   500    1000
#> 27   0.0187   500   500    1000
#> 28   0.1370   500   500    1000
#> 29   0.9203   500   500    1000
#> 30   0.0026   500   500    1000
#> 31   0.3443   500   500    1000
#> 32   0.6960   500   500    1000
#> 33   0.0135   500   500    1000
#> 34   0.0602   500   500    1000
#> 35   0.2801   500   500    1000
#> 36   0.0005   500   500    1000
#> 

# (b) Compute GRDIF statistics with purification
#     based on GRDIF_R, including post-hoc comparisons
dif_puri_r <- grdif(
  x = est_par, data = data, score = score, group = group,
  focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "grdifr", post.hoc = TRUE
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 Iteration: 3 Iteration: 4 Iteration: 5 Iteration: 6 
#> Purification is finished. 
print(dif_puri_r)
#> 
#> Call:
#> grdif.default(x = est_par, data = data, score = score, group = group, 
#>     focal.name = c(1, 2, 3), D = 1, alpha = 0.05, purify = TRUE, 
#>     purify.by = "grdifr", post.hoc = TRUE)
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
#> 1   V1   500    500    500    500 13.154    0.004  **  8.058    0.045   *
#> 2   V2   500    500    500    500 96.369    0.000 ***  2.000    0.572    
#> 3   V3   500    500    500    500 29.873    0.000 *** 19.246    0.000 ***
#> 4   V4   500    500    500    500 77.456    0.000 *** 41.597    0.000 ***
#> 5   V5   500    500    500    500  5.607    0.132      3.636    0.304    
#> 6   V6   500    500    500    500  1.658    0.646      0.928    0.819    
#> 7   V7   500    500    500    500  0.739    0.864      2.728    0.436    
#> 8   V8   500    500    500    500  1.220    0.748      4.916    0.178    
#> 9   V9   500    500    500    500  0.672    0.880      0.320    0.956    
#> 10 V10   500    500    500    500  1.779    0.620      5.262    0.154    
#> 11 V11   500    500    500    500  3.533    0.317      0.958    0.811    
#> 12 V12   500    500    500    500  1.025    0.795      0.630    0.889    
#> 13 V13   500    500    500    500  1.001    0.801      2.942    0.401    
#> 14 V14   500    500    500    500  5.495    0.139      5.752    0.124    
#> 15 V15   500    500    500    500  4.650    0.199      5.189    0.158    
#> 16 V16   500    500    500    500  1.504    0.681      0.346    0.951    
#> 17 V17   500    500    500    500  3.102    0.376      4.323    0.229    
#> 18 V18   500    500    500    500  2.234    0.525      2.210    0.530    
#> 19 V19   500    500    500    500  1.336    0.721      0.226    0.973    
#> 20 V20   500    500    500    500 16.500    0.001 ***  0.900    0.825    
#> 21 V21   500    500    500    500  5.680    0.128      6.536    0.088   .
#> 22 V22   500    500    500    500  7.513    0.057   .  6.500    0.090   .
#> 23 V23   500    500    500    500  0.392    0.942      4.418    0.220    
#> 24 V24   500    500    500    500  5.900    0.117      6.070    0.108    
#> 25 V25   500    500    500    500  4.882    0.181      4.858    0.182    
#> 26 V26   500    500    500    500  1.038    0.792      1.440    0.696    
#> 27 V27   500    500    500    500  1.213    0.750      1.710    0.635    
#> 28 V28   500    500    500    500  7.098    0.069   .  1.724    0.631    
#> 29 V29   500    500    500    500  0.735    0.865      0.500    0.919    
#> 30 V30   500    500    500    500  1.216    0.749      2.597    0.458    
#> 31 V31   500    500    500    500  1.558    0.669      4.717    0.194    
#> 32 V32   500    500    500    500  2.913    0.405      1.503    0.681    
#> 33 V33   500    500    500    500  2.173    0.537      5.840    0.120    
#> 34 V34   500    500    500    500  3.954    0.266      7.293    0.063   .
#> 35 V35   500    500    500    500  1.245    0.742      2.402    0.493    
#> 36 V36   500    500    500    500  3.275    0.351      0.668    0.881    
#> 37 V37   500    500    500    500  0.358    0.949      0.572    0.903    
#> 38 V38   500    500    500    500  6.205    0.102      1.738    0.628    
#> 39 V39   500    500    500    500  2.417    0.490      1.414    0.702    
#> 40 V40   500    500    500    500  1.538    0.674      4.804    0.187    
#>    grdifrs p.grdifrs    
#> 1    15.76     0.015   *
#> 2   102.89     0.000 ***
#> 3    36.03     0.000 ***
#> 4    83.28     0.000 ***
#> 5     6.24     0.396    
#> 6     1.77     0.939    
#> 7     3.43     0.753    
#> 8     5.51     0.480    
#> 9     1.17     0.979    
#> 10    5.81     0.444    
#> 11    3.81     0.702    
#> 12    1.97     0.922    
#> 13    3.31     0.769    
#> 14   15.32     0.018   *
#> 15    6.29     0.392    
#> 16    3.57     0.735    
#> 17    6.35     0.385    
#> 18    2.96     0.814    
#> 19    1.60     0.953    
#> 20   17.55     0.007  **
#> 21    9.66     0.140    
#> 22   11.12     0.085   .
#> 23    5.10     0.531    
#> 24    9.45     0.150    
#> 25   10.87     0.093   .
#> 26    4.17     0.654    
#> 27    3.64     0.725    
#> 28    9.28     0.158    
#> 29    1.24     0.975    
#> 30    6.31     0.390    
#> 31    7.46     0.281    
#> 32    5.48     0.484    
#> 33    8.61     0.196    
#> 34    8.55     0.200    
#> 35    2.59     0.859    
#> 36    4.09     0.664    
#> 37    1.42     0.965    
#> 38   10.45     0.107    
#> 39    5.48     0.484    
#> 40    8.63     0.195    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 6
#>   - GRDIF statistic used for purification: GRDIF(R)
#>   - DIF Items identified by GRDIF(R): 
#>     1, 2, 3, 4, 20, 24 
#>   - GRDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc1 n.foc2 n.foc3 grdifr p.grdifr     grdifs p.grdifs
#> 1   V1      3   500    500    500    500 17.959    0.000 ***  8.918    0.030
#> 2   V2      0   500    500    500    500 96.369    0.000 ***  2.000    0.572
#> 3   V3      2   500    500    500    500 34.519    0.000 *** 20.320    0.000
#> 4   V4      1   500    500    500    500 81.387    0.000 *** 41.038    0.000
#> 5   V5      6   500    500    500    500  4.006    0.261      2.663    0.447
#> 6   V6      6   500    500    500    500  2.274    0.517      0.840    0.840
#> 7   V7      6   500    500    500    500  2.158    0.540      2.902    0.407
#> 8   V8      6   500    500    500    500  0.654    0.884      4.189    0.242
#> 9   V9      6   500    500    500    500  0.328    0.955      0.481    0.923
#> 10 V10      6   500    500    500    500  0.979    0.806      4.924    0.178
#> 11 V11      6   500    500    500    500  4.117    0.249      0.939    0.816
#> 12 V12      6   500    500    500    500  0.426    0.935      1.419    0.701
#> 13 V13      6   500    500    500    500  0.247    0.970      1.294    0.731
#> 14 V14      6   500    500    500    500  5.364    0.147      9.885    0.020
#> 15 V15      6   500    500    500    500  3.167    0.367      4.340    0.227
#> 16 V16      6   500    500    500    500  1.131    0.770      0.308    0.959
#> 17 V17      6   500    500    500    500  1.771    0.621      3.608    0.307
#> 18 V18      6   500    500    500    500  1.280    0.734      1.711    0.634
#> 19 V19      6   500    500    500    500  1.089    0.780      0.039    0.998
#> 20 V20      4   500    500    500    500 10.774    0.013   *  0.261    0.967
#> 21 V21      6   500    500    500    500  7.617    0.055   .  5.653    0.130
#> 22 V22      6   500    500    500    500  6.856    0.077   .  5.938    0.115
#> 23 V23      6   500    500    500    500  1.226    0.747      2.218    0.528
#> 24 V24      5   500    500    500    500  8.239    0.041   *  7.726    0.052
#> 25 V25      6   500    500    500    500  2.629    0.452      2.330    0.507
#> 26 V26      6   500    500    500    500  0.357    0.949      1.785    0.618
#> 27 V27      6   500    500    500    500  0.375    0.945      2.320    0.509
#> 28 V28      6   500    500    500    500  5.696    0.127      2.241    0.524
#> 29 V29      6   500    500    500    500  0.309    0.958      0.246    0.970
#> 30 V30      6   500    500    500    500  0.382    0.944      0.748    0.862
#> 31 V31      6   500    500    500    500  2.025    0.567      4.930    0.177
#> 32 V32      6   500    500    500    500  0.319    0.957      0.932    0.818
#> 33 V33      6   500    500    500    500  0.581    0.901      4.083    0.253
#> 34 V34      6   500    500    500    500  6.327    0.097   .  8.231    0.042
#> 35 V35      6   500    500    500    500  0.671    0.880      2.488    0.477
#> 36 V36      6   500    500    500    500  1.290    0.731      0.374    0.946
#> 37 V37      6   500    500    500    500  1.462    0.691      1.537    0.674
#> 38 V38      6   500    500    500    500  5.422    0.143      2.220    0.528
#> 39 V39      6   500    500    500    500  0.942    0.815      1.878    0.598
#> 40 V40      6   500    500    500    500  0.764    0.858      5.295    0.151
#>        grdifrs p.grdifrs    
#> 1    *   20.32     0.002  **
#> 2       102.89     0.000 ***
#> 3  ***   40.17     0.000 ***
#> 4  ***   86.93     0.000 ***
#> 5         4.32     0.633    
#> 6         2.77     0.837    
#> 7         4.23     0.646    
#> 8         4.56     0.601    
#> 9         1.03     0.985    
#> 10        6.17     0.405    
#> 11        4.28     0.639    
#> 12        2.20     0.901    
#> 13        1.39     0.967    
#> 14   *   19.39     0.004  **
#> 15        5.65     0.463    
#> 16        2.87     0.825    
#> 17        4.20     0.649    
#> 18        2.43     0.876    
#> 19        1.46     0.962    
#> 20       12.04     0.061   .
#> 21       11.94     0.063   .
#> 22        9.06     0.170    
#> 23        3.33     0.766    
#> 24   .   12.21     0.058   .
#> 25        6.14     0.408    
#> 26        4.53     0.606    
#> 27        3.26     0.776    
#> 28        8.02     0.237    
#> 29        1.13     0.980    
#> 30        4.91     0.555    
#> 31        6.19     0.402    
#> 32        1.62     0.951    
#> 33        8.41     0.209    
#> 34   *   10.68     0.099   .
#> 35        2.53     0.865    
#> 36        1.65     0.949    
#> 37        3.14     0.791    
#> 38        9.81     0.133    
#> 39        3.24     0.778    
#> 40        7.71     0.260    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 

# Display post-hoc results before purification
print(dif_puri_r$no_purify$post.hoc)
#> $by.grdifr
#>     id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1   V1      0 & 1 -0.0571 -1.9787 -0.0143 -1.5392  5.4726  0.0478  0.1238
#> 2   V1      0 & 2 -0.0453 -1.6508 -0.0214  1.4492  3.2189  0.0988  0.1473
#> 3   V1      0 & 3  0.0323  1.1415 -0.0163 -0.3473  1.3143  0.2536  0.7284
#> 4   V1      1 & 2  0.0118  0.4281 -0.0072  2.8172  8.4417  0.6686  0.0048
#> 5   V1      1 & 3  0.0894  3.1484 -0.0020  1.1352 11.0095  0.0016  0.2563
#> 6   V1      2 & 3  0.0776  2.8783  0.0051 -1.7120  8.2949  0.0040  0.0869
#> 7   V2      0 & 1 -0.1060 -4.5351 -0.0494 -1.3017 20.6424  0.0000  0.1930
#> 8   V2      0 & 2  0.0013  0.0559 -0.0332 -0.2033  0.0417  0.9554  0.8389
#> 9   V2      0 & 3  0.1196  4.9024 -0.0171 -0.6807 24.1845  0.0000  0.4961
#> 10  V2      1 & 2  0.1074  4.8742  0.0162  1.1210 24.0713  0.0000  0.2623
#> 11  V2      1 & 3  0.2256  9.8131  0.0323  0.6452 98.9891  0.0000  0.5188
#> 12  V2      2 & 3  0.1183  5.1209  0.0161 -0.4860 27.6335  0.0000  0.6270
#> 13  V3      0 & 1 -0.0269 -0.9246 -0.0225 -2.2768  5.7699  0.3552  0.0228
#> 14  V3      0 & 2 -0.0910 -3.2646 -0.0169  1.5449 10.8342  0.0011  0.1224
#> 15  V3      0 & 3  0.0565  1.9705 -0.0372 -2.4747  6.9738  0.0488  0.0133
#> 16  V3      1 & 2 -0.0641 -2.2990  0.0055  3.5583 15.5463  0.0215  0.0004
#> 17  V3      1 & 3  0.0834  2.9071 -0.0147 -0.2808  8.4681  0.0036  0.7789
#> 18  V3      2 & 3  0.1476  5.3711 -0.0202 -3.7051 28.8580  0.0000  0.0002
#> 19  V4      0 & 1 -0.0349 -1.3083  0.0130  0.4067  1.9900  0.1908  0.6842
#> 20  V4      0 & 2 -0.1382 -5.6951  0.0132  5.4802 41.4631  0.0000  0.0000
#> 21  V4      0 & 3  0.0610  2.3838 -0.0194  0.0301  6.8637  0.0171  0.9760
#> 22  V4      1 & 2 -0.1034 -4.2016  0.0002  4.9559 39.1282  0.0000  0.0000
#> 23  V4      1 & 3  0.0958  3.7010 -0.0324 -0.3708 13.8353  0.0002  0.7108
#> 24  V4      2 & 3  0.1992  8.4964 -0.0326 -5.3629 72.2070  0.0000  0.0000
#> 25 V20      0 & 1 -0.0327 -1.3665 -0.0125  0.1019  2.1325  0.1718  0.9188
#> 26 V20      0 & 2  0.0143  0.6303 -0.0520 -0.7158  0.7247  0.5285  0.4741
#> 27 V20      0 & 3 -0.0694 -2.9073 -0.0151  0.0406  8.6031  0.0036  0.9676
#> 28 V20      1 & 2  0.0470  2.1321 -0.0395 -0.7806  5.6204  0.0330  0.4350
#> 29 V20      1 & 3 -0.0367 -1.5751 -0.0026 -0.0642  2.5453  0.1152  0.9488
#> 30 V20      2 & 3 -0.0836 -3.8099  0.0368  0.7635 15.0786  0.0001  0.4452
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0648   500   500    1000
#> 2    0.2000   500   500    1000
#> 3    0.5183   500   500    1000
#> 4    0.0147   500   500    1000
#> 5    0.0041   500   500    1000
#> 6    0.0158   500   500    1000
#> 7    0.0000   500   500    1000
#> 8    0.9794   500   500    1000
#> 9    0.0000   500   500    1000
#> 10   0.0000   500   500    1000
#> 11   0.0000   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.0559   500   500    1000
#> 14   0.0044   500   500    1000
#> 15   0.0306   500   500    1000
#> 16   0.0004   500   500    1000
#> 17   0.0145   500   500    1000
#> 18   0.0000   500   500    1000
#> 19   0.3697   500   500    1000
#> 20   0.0000   500   500    1000
#> 21   0.0323   500   500    1000
#> 22   0.0000   500   500    1000
#> 23   0.0010   500   500    1000
#> 24   0.0000   500   500    1000
#> 25   0.3443   500   500    1000
#> 26   0.6960   500   500    1000
#> 27   0.0135   500   500    1000
#> 28   0.0602   500   500    1000
#> 29   0.2801   500   500    1000
#> 30   0.0005   500   500    1000
#> 
#> $by.grdifs
#>    id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1  V1      0 & 1 -0.0571 -1.9787 -0.0143 -1.5392  5.4726  0.0478  0.1238
#> 2  V1      0 & 2 -0.0453 -1.6508 -0.0214  1.4492  3.2189  0.0988  0.1473
#> 3  V1      0 & 3  0.0323  1.1415 -0.0163 -0.3473  1.3143  0.2536  0.7284
#> 4  V1      1 & 2  0.0118  0.4281 -0.0072  2.8172  8.4417  0.6686  0.0048
#> 5  V1      1 & 3  0.0894  3.1484 -0.0020  1.1352 11.0095  0.0016  0.2563
#> 6  V1      2 & 3  0.0776  2.8783  0.0051 -1.7120  8.2949  0.0040  0.0869
#> 7  V3      0 & 1 -0.0269 -0.9246 -0.0225 -2.2768  5.7699  0.3552  0.0228
#> 8  V3      0 & 2 -0.0910 -3.2646 -0.0169  1.5449 10.8342  0.0011  0.1224
#> 9  V3      0 & 3  0.0565  1.9705 -0.0372 -2.4747  6.9738  0.0488  0.0133
#> 10 V3      1 & 2 -0.0641 -2.2990  0.0055  3.5583 15.5463  0.0215  0.0004
#> 11 V3      1 & 3  0.0834  2.9071 -0.0147 -0.2808  8.4681  0.0036  0.7789
#> 12 V3      2 & 3  0.1476  5.3711 -0.0202 -3.7051 28.8580  0.0000  0.0002
#> 13 V4      0 & 1 -0.0349 -1.3083  0.0130  0.4067  1.9900  0.1908  0.6842
#> 14 V4      0 & 2 -0.1382 -5.6951  0.0132  5.4802 41.4631  0.0000  0.0000
#> 15 V4      0 & 3  0.0610  2.3838 -0.0194  0.0301  6.8637  0.0171  0.9760
#> 16 V4      1 & 2 -0.1034 -4.2016  0.0002  4.9559 39.1282  0.0000  0.0000
#> 17 V4      1 & 3  0.0958  3.7010 -0.0324 -0.3708 13.8353  0.0002  0.7108
#> 18 V4      2 & 3  0.1992  8.4964 -0.0326 -5.3629 72.2070  0.0000  0.0000
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0648   500   500    1000
#> 2    0.2000   500   500    1000
#> 3    0.5183   500   500    1000
#> 4    0.0147   500   500    1000
#> 5    0.0041   500   500    1000
#> 6    0.0158   500   500    1000
#> 7    0.0559   500   500    1000
#> 8    0.0044   500   500    1000
#> 9    0.0306   500   500    1000
#> 10   0.0004   500   500    1000
#> 11   0.0145   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.3697   500   500    1000
#> 14   0.0000   500   500    1000
#> 15   0.0323   500   500    1000
#> 16   0.0000   500   500    1000
#> 17   0.0010   500   500    1000
#> 18   0.0000   500   500    1000
#> 
#> $by.grdifrs
#>     id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1   V1      0 & 1 -0.0571 -1.9787 -0.0143 -1.5392  5.4726  0.0478  0.1238
#> 2   V1      0 & 2 -0.0453 -1.6508 -0.0214  1.4492  3.2189  0.0988  0.1473
#> 3   V1      0 & 3  0.0323  1.1415 -0.0163 -0.3473  1.3143  0.2536  0.7284
#> 4   V1      1 & 2  0.0118  0.4281 -0.0072  2.8172  8.4417  0.6686  0.0048
#> 5   V1      1 & 3  0.0894  3.1484 -0.0020  1.1352 11.0095  0.0016  0.2563
#> 6   V1      2 & 3  0.0776  2.8783  0.0051 -1.7120  8.2949  0.0040  0.0869
#> 7   V2      0 & 1 -0.1060 -4.5351 -0.0494 -1.3017 20.6424  0.0000  0.1930
#> 8   V2      0 & 2  0.0013  0.0559 -0.0332 -0.2033  0.0417  0.9554  0.8389
#> 9   V2      0 & 3  0.1196  4.9024 -0.0171 -0.6807 24.1845  0.0000  0.4961
#> 10  V2      1 & 2  0.1074  4.8742  0.0162  1.1210 24.0713  0.0000  0.2623
#> 11  V2      1 & 3  0.2256  9.8131  0.0323  0.6452 98.9891  0.0000  0.5188
#> 12  V2      2 & 3  0.1183  5.1209  0.0161 -0.4860 27.6335  0.0000  0.6270
#> 13  V3      0 & 1 -0.0269 -0.9246 -0.0225 -2.2768  5.7699  0.3552  0.0228
#> 14  V3      0 & 2 -0.0910 -3.2646 -0.0169  1.5449 10.8342  0.0011  0.1224
#> 15  V3      0 & 3  0.0565  1.9705 -0.0372 -2.4747  6.9738  0.0488  0.0133
#> 16  V3      1 & 2 -0.0641 -2.2990  0.0055  3.5583 15.5463  0.0215  0.0004
#> 17  V3      1 & 3  0.0834  2.9071 -0.0147 -0.2808  8.4681  0.0036  0.7789
#> 18  V3      2 & 3  0.1476  5.3711 -0.0202 -3.7051 28.8580  0.0000  0.0002
#> 19  V4      0 & 1 -0.0349 -1.3083  0.0130  0.4067  1.9900  0.1908  0.6842
#> 20  V4      0 & 2 -0.1382 -5.6951  0.0132  5.4802 41.4631  0.0000  0.0000
#> 21  V4      0 & 3  0.0610  2.3838 -0.0194  0.0301  6.8637  0.0171  0.9760
#> 22  V4      1 & 2 -0.1034 -4.2016  0.0002  4.9559 39.1282  0.0000  0.0000
#> 23  V4      1 & 3  0.0958  3.7010 -0.0324 -0.3708 13.8353  0.0002  0.7108
#> 24  V4      2 & 3  0.1992  8.4964 -0.0326 -5.3629 72.2070  0.0000  0.0000
#> 25 V14      0 & 1 -0.0431 -1.5723 -0.0106 -1.3392  3.2138  0.1159  0.1805
#> 26 V14      0 & 2 -0.0129 -0.5099 -0.0428  0.4440  0.3512  0.6101  0.6570
#> 27 V14      0 & 3 -0.0537 -2.0220 -0.0348 -1.5860  7.9560  0.0432  0.1127
#> 28 V14      1 & 2  0.0301  1.1775 -0.0322  1.7309  3.9750  0.2390  0.0835
#> 29 V14      1 & 3 -0.0107 -0.3982 -0.0242 -0.1730  0.1661  0.6905  0.8626
#> 30 V14      2 & 3 -0.0408 -1.6510  0.0079 -1.9889 11.9210  0.0987  0.0467
#> 31 V20      0 & 1 -0.0327 -1.3665 -0.0125  0.1019  2.1325  0.1718  0.9188
#> 32 V20      0 & 2  0.0143  0.6303 -0.0520 -0.7158  0.7247  0.5285  0.4741
#> 33 V20      0 & 3 -0.0694 -2.9073 -0.0151  0.0406  8.6031  0.0036  0.9676
#> 34 V20      1 & 2  0.0470  2.1321 -0.0395 -0.7806  5.6204  0.0330  0.4350
#> 35 V20      1 & 3 -0.0367 -1.5751 -0.0026 -0.0642  2.5453  0.1152  0.9488
#> 36 V20      2 & 3 -0.0836 -3.8099  0.0368  0.7635 15.0786  0.0001  0.4452
#>    p.rdifrs n.ref n.foc n.total
#> 1    0.0648   500   500    1000
#> 2    0.2000   500   500    1000
#> 3    0.5183   500   500    1000
#> 4    0.0147   500   500    1000
#> 5    0.0041   500   500    1000
#> 6    0.0158   500   500    1000
#> 7    0.0000   500   500    1000
#> 8    0.9794   500   500    1000
#> 9    0.0000   500   500    1000
#> 10   0.0000   500   500    1000
#> 11   0.0000   500   500    1000
#> 12   0.0000   500   500    1000
#> 13   0.0559   500   500    1000
#> 14   0.0044   500   500    1000
#> 15   0.0306   500   500    1000
#> 16   0.0004   500   500    1000
#> 17   0.0145   500   500    1000
#> 18   0.0000   500   500    1000
#> 19   0.3697   500   500    1000
#> 20   0.0000   500   500    1000
#> 21   0.0323   500   500    1000
#> 22   0.0000   500   500    1000
#> 23   0.0010   500   500    1000
#> 24   0.0000   500   500    1000
#> 25   0.2005   500   500    1000
#> 26   0.8390   500   500    1000
#> 27   0.0187   500   500    1000
#> 28   0.1370   500   500    1000
#> 29   0.9203   500   500    1000
#> 30   0.0026   500   500    1000
#> 31   0.3443   500   500    1000
#> 32   0.6960   500   500    1000
#> 33   0.0135   500   500    1000
#> 34   0.0602   500   500    1000
#> 35   0.2801   500   500    1000
#> 36   0.0005   500   500    1000
#> 

# Display post-hoc results after purification
print(dif_puri_r$with_purify$post.hoc)
#>     id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1   V2      0 & 1 -0.1060 -4.5351 -0.0494 -1.3017 20.6424  0.0000  0.1930
#> 2   V2      0 & 2  0.0013  0.0559 -0.0332 -0.2033  0.0417  0.9554  0.8389
#> 3   V2      0 & 3  0.1196  4.9024 -0.0171 -0.6807 24.1845  0.0000  0.4961
#> 4   V2      1 & 2  0.1074  4.8742  0.0162  1.1210 24.0713  0.0000  0.2623
#> 5   V2      1 & 3  0.2256  9.8131  0.0323  0.6452 98.9891  0.0000  0.5188
#> 6   V2      2 & 3  0.1183  5.1209  0.0161 -0.4860 27.6335  0.0000  0.6270
#> 7   V4      0 & 1 -0.0408 -1.5337  0.0137  0.5547  2.8226  0.1251  0.5791
#> 8   V4      0 & 2 -0.1381 -5.6951  0.0134  5.5106 41.7333  0.0000  0.0000
#> 9   V4      0 & 3  0.0674  2.6303 -0.0170  0.1021  8.4403  0.0085  0.9187
#> 10  V4      1 & 2 -0.0974 -3.9671 -0.0003  4.8151 35.7952  0.0001  0.0000
#> 11  V4      1 & 3  0.1081  4.1774 -0.0306 -0.4468 17.6506  0.0000  0.6550
#> 12  V4      2 & 3  0.2055  8.7475 -0.0303 -5.3239 76.5180  0.0000  0.0000
#> 13  V3      0 & 1 -0.0319 -1.0951 -0.0227 -2.2444  5.9556  0.2735  0.0248
#> 14  V3      0 & 2 -0.0948 -3.4050 -0.0162  1.7262 11.6876  0.0007  0.0843
#> 15  V3      0 & 3  0.0636  2.2129 -0.0360 -2.4537  7.5664  0.0269  0.0141
#> 16  V3      1 & 2 -0.0630 -2.2617  0.0065  3.6955 16.1752  0.0237  0.0002
#> 17  V3      1 & 3  0.0954  3.3224 -0.0133 -0.2731 11.0527  0.0009  0.7848
#> 18  V3      2 & 3  0.1584  5.7670 -0.0198 -3.8541 33.2592  0.0000  0.0001
#> 19  V1      0 & 1 -0.0629 -2.1836 -0.0157 -1.6652  6.6177  0.0290  0.0959
#> 20  V1      0 & 2 -0.0514 -1.8796 -0.0231  1.4629  3.8687  0.0602  0.1435
#> 21  V1      0 & 3  0.0414  1.4585 -0.0170 -0.5477  2.1275  0.1447  0.5839
#> 22  V1      1 & 2  0.0115  0.4202 -0.0074  2.9357  9.2150  0.6743  0.0033
#> 23  V1      1 & 3  0.1043  3.6691 -0.0013  1.0676 14.3465  0.0002  0.2857
#> 24  V1      2 & 3  0.0928  3.4443  0.0061 -1.9103 11.9217  0.0006  0.0561
#> 25 V20      0 & 1 -0.0462 -1.9292 -0.0129 -0.0532  3.9816  0.0537  0.9576
#> 26 V20      0 & 2  0.0026  0.1149 -0.0511 -0.4672  0.2183  0.9085  0.6403
#> 27 V20      0 & 3 -0.0554 -2.3232 -0.0163 -0.2081  5.6225  0.0202  0.8352
#> 28 V20      1 & 2  0.0487  2.2224 -0.0382 -0.3874  5.2949  0.0263  0.6985
#> 29 V20      1 & 3 -0.0093 -0.3974 -0.0034 -0.1444  0.1625  0.6911  0.8852
#> 30 V20      2 & 3 -0.0580 -2.6539  0.0348  0.2525  7.6502  0.0080  0.8007
#> 31 V24      0 & 1 -0.0451 -2.0079  0.1025  2.4548  6.4872  0.0447  0.0141
#> 32 V24      0 & 2 -0.0392 -2.4982 -0.0290  2.4438  6.6578  0.0125  0.0145
#> 33 V24      0 & 3 -0.0156 -0.8792 -0.0027  1.6489  3.3539  0.3793  0.0992
#> 34 V24      1 & 2  0.0059  0.2951 -0.1315 -0.2623  0.1012  0.7679  0.7931
#> 35 V24      1 & 3  0.0295  1.3734 -0.1052 -0.8798  1.8975  0.1696  0.3790
#> 36 V24      2 & 3  0.0237  1.6545  0.0263 -0.7047  4.4341  0.0980  0.4810
#>    p.rdifrs n.ref n.foc n.total n_iter
#> 1    0.0000   500   500    1000      0
#> 2    0.9794   500   500    1000      0
#> 3    0.0000   500   500    1000      0
#> 4    0.0000   500   500    1000      0
#> 5    0.0000   500   500    1000      0
#> 6    0.0000   500   500    1000      0
#> 7    0.2438   500   500    1000      1
#> 8    0.0000   500   500    1000      1
#> 9    0.0147   500   500    1000      1
#> 10   0.0000   500   500    1000      1
#> 11   0.0001   500   500    1000      1
#> 12   0.0000   500   500    1000      1
#> 13   0.0509   500   500    1000      2
#> 14   0.0029   500   500    1000      2
#> 15   0.0228   500   500    1000      2
#> 16   0.0003   500   500    1000      2
#> 17   0.0040   500   500    1000      2
#> 18   0.0000   500   500    1000      2
#> 19   0.0366   500   500    1000      3
#> 20   0.1445   500   500    1000      3
#> 21   0.3452   500   500    1000      3
#> 22   0.0100   500   500    1000      3
#> 23   0.0008   500   500    1000      3
#> 24   0.0026   500   500    1000      3
#> 25   0.1366   500   500    1000      4
#> 26   0.8966   500   500    1000      4
#> 27   0.0601   500   500    1000      4
#> 28   0.0708   500   500    1000      4
#> 29   0.9220   500   500    1000      4
#> 30   0.0218   500   500    1000      4
#> 31   0.0390   500   500    1000      5
#> 32   0.0358   500   500    1000      5
#> 33   0.1869   500   500    1000      5
#> 34   0.9506   500   500    1000      5
#> 35   0.3872   500   500    1000      5
#> 36   0.1089   500   500    1000      5
# }
```
