# Residual-Based DIF Detection Framework Using Categorical Residuals (RDIF-CR)

This function computes three statistics of the residual-based DIF
detection framework using categorical residuals
(RDIF-CR)—\\RDIF\_{R}-CR\\, \\RDIF\_{S}-CR\\, and \\RDIF\_{RS}-CR\\—for
detecting global differential item functioning (DIF), particularly in
polytomously scored items. The RDIF-CR framework evaluates DIF by
comparing categorical residual vectors, which are calculated as the
difference between a one-hot encoded response vector (with 1 for the
selected category and 0 for all others) and the IRT model–predicted
probability vector across all score categories. This approach enables
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
  values). If not provided, `crdif()` will estimate ability parameters
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
  hypothesis testing using the CRDIF statistics. Default is `0.05`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying which RDIF statistic is used to perform
  the purification. Available options are `"crdifrs"` for
  \\RDIF\_{RS}-CR\\, `"crdifr"` for \\RDIF\_{R}-CR\\, and `"crdifs"` for
  \\RDIF\_{S}-CR\\.

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

  A list of sub-objects containing the results of DIF analysis without
  applying a purification procedure. The sub-objects include:

  dif_stat

  :   A data frame summarizing the RDIF-CR analysis results for all
      items. Columns include item ID, \\RDIF\_{R}-CR\\, degrees of
      freedom, \\RDIF\_{S}-CR\\, degrees of freedom, \\RDIF\_{RS}-CR\\,
      degrees of freedom, associated p-values, and sample sizes for the
      reference and focal groups.

  moments

  :   A list containing the first and second moments (means and
      covariance matrices) of the RDIF-CR statistics. The elements
      include: `mu.crdifr`, `mu.crdifs`, `mu.crdifrs` (means), and
      `cov.crdifr`, `cov.crdifs`, `cov.crdifrs` (covariances), each
      indexed by item ID.

  dif_item

  :   A list of three numeric vectors identifying items flagged as DIF
      based on each statistic: `crdifr`, `crdifs`, and `crdifrs`.

  score

  :   A numeric vector of ability estimates used to compute the RDIF-CR
      statistics. These may be user-supplied or internally estimated.

- purify:

  A logical value indicating whether a purification procedure was
  applied.

- with_purify:

  A list of sub-objects containing the results of DIF analysis after
  applying the purification procedure. The sub-objects include:

  purify.by

  :   A character string indicating the RDIF-CR statistic used for
      purification. Possible values are "crdifr", "crdifs", or
      "crdifrs".

  dif_stat

  :   A data frame summarizing the final RDIF-CR statistics after
      purification. Same structure as in `no_purify`, with an additional
      column indicating the iteration in which the result was obtained.

  moments

  :   A list of moments (means and covariance matrices) of the RDIF-CR
      statistics for all items, updated based on the final iteration.

  dif_item

  :   A list of three numeric vectors identifying items flagged as DIF
      at any iteration, by each statistic.

  n.iter

  :   An integer indicating the number of iterations performed during
      the purification procedure.

  score

  :   A numeric vector of updated ability estimates used in the final
      iteration.

  complete

  :   A logical value indicating whether the purification process
      converged. If `FALSE`, the maximum number of iterations was
      reached before convergence.

- alpha:

  A numeric value indicating the significance level (\\\alpha\\) used
  for hypothesis testing with RDIF-CR statistics.

## Details

According to Penfield (2010), differential item functioning (DIF) in
polytomously scored items can be conceptualized in two forms: global DIF
and net DIF. Global DIF refers to differences between groups in the
conditional probabilities of responding in specific score categories,
thus offering a fine-grained view of DIF at the category level. In
contrast, net DIF summarizes these differences into a single value
representing the overall impact of DIF on the item’s expected score.

The RDIF framework using categorical residuals (RDIF-CR), implemented in
`crdif()`, extends the original residual-based DIF framework proposed by
Lim et al. (2022) to detect global DIF in polytomous items. This
framework includes three statistics: \\RDIF\_{R}-CR\\, \\RDIF\_{S}-CR\\,
and \\RDIF\_{RS}-CR\\, each designed to capture different aspects of
group-level differences in categorical response patterns.

To illustrate how the RDIF-CR framework operates, consider an item with
five ordered score categories (\\k \in \\0,1,2,3,4\\\\). Suppose an
examinee with latent ability \\\theta\\ responds with category 2. The
one-hot encoded response vector for this response is \\(0,0,1,0,0)^T\\.
Assume that the IRT model estimates the examinee’s expected score as 2.5
and predicts the category probabilities as \\(0.1, 0.2, 0.4, 0.25,
0.05)^T\\. In the RDIF-CR framework, the categorical residual vector is
calculated by subtracting the predicted probability vector from the
one-hot response vector, resulting in \\(-0.1, -0.2, 0.6, -0.25,
-0.05)^T\\.

In contrast to the RDIF-CR framework, net DIF is assessed using a
unidimensional item score residual. In this example, the residual would
be \\2 - 2.5 = -0.5\\. For detecting net DIF, the
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) function
should be used instead.

Note that for dichotomous items, `crdif()` and
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) yield
identical results. This is because the categorical probability vector
for a binary item reduces to a scalar difference, making the global and
net DIF evaluations mathematically equivalent.

## Methods (by class)

- `crdif(default)`: Default method for computing the three RDIF-CR
  statistics using a data frame `x` that contains item metadata

- `crdif(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `crdif(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

## References

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Penfield, R. D. (2010). Distinguishing between net and global DIF in
polytomous items. *Journal of Educational Measurement, 47*(2), 129–149.

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
#> Estimation is finished in 0.44 seconds. 

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
#> 1 V1  1000  1000 181.97         4    0.000 *** 156.06         4    0.000 ***
#> 2 V2  1000  1000   2.00         4    0.737       8.19         4    0.085   .
#> 3 V3  1000  1000   3.86         4    0.425       2.68         4    0.613    
#> 4 V4  1000  1000   5.57         4    0.234       5.29         4    0.258    
#> 5 V5  1000  1000 409.88         4    0.000 *** 606.02         4    0.000 ***
#>   crdifrs df.crdifrs p.crdifrs    
#> 1  194.66          8     0.000 ***
#> 2    9.25          8     0.322    
#> 3    8.57          8     0.380    
#> 4   12.15          8     0.145    
#> 5  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (b) DIF detection with purification using RDIF_{R}-CR
dif_puri_1 <- crdif(
  x = x, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "crdifr"
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 
#> Purification is finished. 
print(dif_puri_1)
#> 
#> Call:
#> crdif.default(x = x, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, purify = TRUE, purify.by = "crdifr")
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
#> 1 V1  1000  1000 181.97         4    0.000 *** 156.06         4    0.000 ***
#> 2 V2  1000  1000   2.00         4    0.737       8.19         4    0.085   .
#> 3 V3  1000  1000   3.86         4    0.425       2.68         4    0.613    
#> 4 V4  1000  1000   5.57         4    0.234       5.29         4    0.258    
#> 5 V5  1000  1000 409.88         4    0.000 *** 606.02         4    0.000 ***
#>   crdifrs df.crdifrs p.crdifrs    
#> 1  194.66          8     0.000 ***
#> 2    9.25          8     0.322    
#> 3    8.57          8     0.380    
#> 4   12.15          8     0.145    
#> 5  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 2
#>   - RDIF-CR statistic used for purification: RDIF(R)-CR
#>   - DIF Items identified by RDIF(R)-CR: 
#>     1, 5 
#>   - RDIF-CR Statistics: 
#> 
#>   id n.iter n.ref n.foc  crdifr df.crdifr p.crdifr      crdifs df.crdifs
#> 1 V1      1  1000  1000 171.536         4    0.000 *** 144.633         4
#> 2 V2      2  1000  1000   0.115         4    0.998       3.229         4
#> 3 V3      2  1000  1000   2.976         4    0.562       2.493         4
#> 4 V4      2  1000  1000   1.096         4    0.895       0.102         4
#> 5 V5      0  1000  1000 409.884         4    0.000 *** 606.018         4
#>   p.crdifs     crdifrs df.crdifrs p.crdifrs    
#> 1    0.000 ***  182.31          8     0.000 ***
#> 2    0.520        3.41          8     0.906    
#> 3    0.646        8.08          8     0.425    
#> 4    0.999        2.06          8     0.979    
#> 5    0.000 ***  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 

# (c) DIF detection with purification using RDIF_{S}-CR
dif_puri_2 <- crdif(
  x = x, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "crdifs"
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 
#> Purification is finished. 
print(dif_puri_2)
#> 
#> Call:
#> crdif.default(x = x, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, purify = TRUE, purify.by = "crdifs")
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
#> 1 V1  1000  1000 181.97         4    0.000 *** 156.06         4    0.000 ***
#> 2 V2  1000  1000   2.00         4    0.737       8.19         4    0.085   .
#> 3 V3  1000  1000   3.86         4    0.425       2.68         4    0.613    
#> 4 V4  1000  1000   5.57         4    0.234       5.29         4    0.258    
#> 5 V5  1000  1000 409.88         4    0.000 *** 606.02         4    0.000 ***
#>   crdifrs df.crdifrs p.crdifrs    
#> 1  194.66          8     0.000 ***
#> 2    9.25          8     0.322    
#> 3    8.57          8     0.380    
#> 4   12.15          8     0.145    
#> 5  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 2
#>   - RDIF-CR statistic used for purification: RDIF(S)-CR
#>   - DIF Items identified by RDIF(S)-CR: 
#>     1, 5 
#>   - RDIF-CR Statistics: 
#> 
#>   id n.iter n.ref n.foc  crdifr df.crdifr p.crdifr      crdifs df.crdifs
#> 1 V1      1  1000  1000 171.536         4    0.000 *** 144.633         4
#> 2 V2      2  1000  1000   0.115         4    0.998       3.229         4
#> 3 V3      2  1000  1000   2.976         4    0.562       2.493         4
#> 4 V4      2  1000  1000   1.096         4    0.895       0.102         4
#> 5 V5      0  1000  1000 409.884         4    0.000 *** 606.018         4
#>   p.crdifs     crdifrs df.crdifrs p.crdifrs    
#> 1    0.000 ***  182.31          8     0.000 ***
#> 2    0.520        3.41          8     0.906    
#> 3    0.646        8.08          8     0.425    
#> 4    0.999        2.06          8     0.979    
#> 5    0.000 ***  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 

# (d) DIF detection with purification using RDIF_{RS}-CR
dif_puri_3 <- crdif(
  x = x, data = data, score = score,
  group = group, focal.name = 1, D = 1, alpha = 0.05,
  purify = TRUE, purify.by = "crdifrs"
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 
#> Purification is finished. 
print(dif_puri_3)
#> 
#> Call:
#> crdif.default(x = x, data = data, score = score, group = group, 
#>     focal.name = 1, D = 1, alpha = 0.05, purify = TRUE, purify.by = "crdifrs")
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
#> 1 V1  1000  1000 181.97         4    0.000 *** 156.06         4    0.000 ***
#> 2 V2  1000  1000   2.00         4    0.737       8.19         4    0.085   .
#> 3 V3  1000  1000   3.86         4    0.425       2.68         4    0.613    
#> 4 V4  1000  1000   5.57         4    0.234       5.29         4    0.258    
#> 5 V5  1000  1000 409.88         4    0.000 *** 606.02         4    0.000 ***
#>   crdifrs df.crdifrs p.crdifrs    
#> 1  194.66          8     0.000 ***
#> 2    9.25          8     0.322    
#> 3    8.57          8     0.380    
#> 4   12.15          8     0.145    
#> 5  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 2
#>   - RDIF-CR statistic used for purification: RDIF(RS)-CR
#>   - DIF Items identified by RDIF(RS)-CR: 
#>     1, 5 
#>   - RDIF-CR Statistics: 
#> 
#>   id n.iter n.ref n.foc  crdifr df.crdifr p.crdifr      crdifs df.crdifs
#> 1 V1      1  1000  1000 171.536         4    0.000 *** 144.633         4
#> 2 V2      2  1000  1000   0.115         4    0.998       3.229         4
#> 3 V3      2  1000  1000   2.976         4    0.562       2.493         4
#> 4 V4      2  1000  1000   1.096         4    0.895       0.102         4
#> 5 V5      0  1000  1000 409.884         4    0.000 *** 606.018         4
#>   p.crdifs     crdifrs df.crdifrs p.crdifrs    
#> 1    0.000 ***  182.31          8     0.000 ***
#> 2    0.520        3.41          8     0.906    
#> 3    0.646        8.08          8     0.425    
#> 4    0.999        2.06          8     0.979    
#> 5    0.000 ***  664.04          8     0.000 ***
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 

# }
```
