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
  ...
)
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
  function.

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- score:

  A numeric vector containing examinees' ability estimates (theta
  values). If not provided, `catsib()` will estimate ability parameters
  internally before computing the CATSIB statistics. See
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  for more information on scoring methods. Default is `NULL`.

- se:

  A vector of standard errors corresponding to the ability estimates.
  The order of the standard errors must match the order of the ability
  estimates provided in the `score` argument. Default is `NULL`.

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

- n.bin:

  A numeric vector of two positive integers specifying the maximum and
  minimum numbers of bins (or intervals) on the ability scale. The first
  and second values represent the maximum and minimum numbers of bins,
  respectively. Default is `c(80, 10)`. See the **Details** section
  below for more information.

- min.binsize:

  A positive integer specifying the minimum number of examinees required
  in each bin. To ensure stable statistical estimation, each bin must
  contain at least the specified number of examinees from both the
  reference and focal groups in order to be included in the calculation
  of \\\hat{\beta}\\. Bins that do not meet this minimum are excluded
  from the computation. Default is 3. See the **Details** section for
  further explanation.

- max.del:

  A numeric value specifying the maximum allowable proportion of
  examinees that may be excluded from either the reference or focal
  group during the binning process. This threshold is used when
  determining the number of bins on the ability scale automatically.
  Default is 0.075. See the **Details** section for more information.

- weight.group:

  A character string specifying the target ability distribution used to
  compute the expected DIF measure \\\hat{\beta}\\ and its corresponding
  standard error. Available options are: `"comb"` for the combined
  distribution of both the reference and focal groups, `"foc"` for the
  focal group's distribution, and `"ref"` for the reference group's
  distribution. Default is `"comb"`. See the **Details** section below
  for more information.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for the
  hypothesis test associated with the CATSIB (*beta*) statistic. Default
  is 0.05.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

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

- ...:

  Additional arguments passed to the
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  function.

## Value

This function returns a list consisting of four elements:

- no_purify:

  A list containing the results of the DIF analysis without applying a
  purification procedure. This list includes:

  dif_stat

  :   A data frame containing the results of the CATSIB statistics for
      all evaluated items. The columns include the item ID, CATSIB
      (*beta*) statistic, standard error of *beta*, standardized *beta*,
      p-value for *beta*, sample size of the reference group, sample
      size of the focal group, and total sample size.

  dif_item

  :   A numeric vector identifying items flagged as potential DIF items
      based on the CATSIB statistic.

  contingency

  :   A list of contingency tables used for computing the CATSIB
      statistics for each item.

- purify:

  A logical value indicating whether a purification procedure was
  applied.

- with_purify:

  A list containing the results of the DIF analysis with a purification
  procedure. This list includes:

  dif_stat

  :   A data frame containing the results of the CATSIB statistics for
      all evaluated items. The columns include the item ID, CATSIB
      (*beta*) statistic, standard error of *beta*, standardized *beta*,
      p-value for *beta*, sample size of the reference group, sample
      size of the focal group, total sample size, and the iteration
      number (*n*) in which the CATSIB statistics were computed.

  dif_item

  :   A numeric vector identifying items flagged as potential DIF items
      based on the CATSIB statistic.

  n.iter

  :   An integer indicating the total number of iterations performed
      during the purification process.

  complete

  :   A logical value indicating whether the purification process was
      completed. If FALSE, the process reached the maximum number of
      iterations without full convergence.

  contingency

  :   A list of contingency tables used for computing the CATSIB
      statistics for each item during the purification process.

- alpha:

  The significance level \\\alpha\\ used to compute the p-values of the
  CATSIB statistics.

## Details

In the CATSIB procedure (Nandakumar & Roussos, 2004),
\\\hat{\theta}^{\ast}\\— the expected value of \\\theta\\ regressed on
\\\hat{\theta}\\—is a continuous variable. The range of
\\\hat{\theta}^{\ast}\\ is divided into *K* equal-width intervals, and
examinees are classified into one of these *K* intervals based on their
\\\hat{\theta}^{\ast}\\ values. Any interval containing fewer than three
examinees from either the reference or focal group is excluded from the
computation of \\\hat{\beta}\\, the DIF effect size, to ensure
statistical stability. According to Nandakumar and Roussos (2004), the
default minimum bin size is 3, which can be controlled via the
`min.binsize` argument.

To determine an appropriate number of intervals (*K*), `catsib()`
automatically decreases *K* from a large starting value (e.g., 80) based
on the rule proposed by Nandakumar and Roussos (2004). Specifically, if
more than 7.5% of examinees in either the reference or focal group would
be excluded due to small bin sizes, the number of bins is reduced by one
and the process is repeated. This continues until the retained examinees
in each group comprise at least 92.5\\ few bins, they recommended a
minimum of *K* = 10. Therefore, the default maximum and minimum number
of bins are set to 80 and 10, respectively, via `n.bin`. Likewise, the
maximum allowable proportion of excluded examinees is set to 0.075 by
default through the `max.del` argument.

When it comes to the target ability distribution used to compute
\\\hat{\beta}\\, Li and Stout (1996) and Nandakumar and Roussos (2004)
employed the combined-group target ability distribution, which is the
default option in `weight.group`. See Nandakumar and Roussos (2004) for
further details about the CATSIB method.

Although Nandakumar and Roussos (2004) did not propose a purification
procedure for DIF analysis using CATSIB, `catsib()` can implement an
iterative purification process in a manner similar to that of Lim et al.
(2022). Specifically, at each iteration, examinees' latent abilities are
recalculated using the purified set of items and the scoring method
specified in the `method` argument. The iterative purification process
terminates either when no additional DIF items are detected or when the
number of iterations reaches the limit set by `max.iter`. See Lim et al.
(2022) for more details on the purification procedure.

Scoring based on a limited number of items may result in large standard
errors, which can negatively affect the effectiveness of DIF detection
using the CATSIB procedure. The `min.resp` argument can be used to
prevent the use of scores with large standard errors, particularly
during the purification process. For example, if `min.resp` is not NULL
(e.g., `min.resp = 5`), item responses from examinees whose total number
of valid responses is below the specified threshold are treated as
missing (i.e., NA). As a result, their ability estimates are also
treated as missing and are excluded from the CATSIB statistic
computation. If `min.resp = NULL`, a score will be computed for any
examinee with at least one valid item response.

Note that the regression correction (Eq. 7 in Nandakumar & Roussos,
2004) assumes \\\hat{\rho}^2\\ (the estimated reliability of ability
estimates) lies in \\\[0, 1\]\\. In practice, however, \\\hat{\rho}^2\\
can become negative when the mean squared standard error of ability
estimates exceeds the observed variance of ability estimates — a
situation that can arise when (a) the number of items is very small, (b)
a purification procedure removes many items, or (c) items exhibiting
nonuniform DIF inflate the standard errors of focal group examinees. A
negative \\\hat{\rho}^2\\ causes the regression correction to amplify
rather than attenuate group differences, leading to inflated Type I
error rates. Even a small positive \\\hat{\rho}^2\\ (e.g., 0.03) can
collapse the corrected ability scores so tightly around each group's
mean that, when ability impact exists between groups, the two groups'
corrected score distributions no longer overlap. This leaves no bins
containing examinees from both groups, resulting in \\\hat{\beta} = 0\\
and \\\text{SE}(\hat{\beta}) = 0\\ for every item, which causes the
purification loop to terminate early with invalid statistics. To prevent
this correction collapse, `catsib()` enforces a floor of 0.05 on
\\\hat{\rho}^2\\ — i.e., \\\hat{\rho}^2 = \max(0.05, \min(1, 1 -
\hat{\sigma}\_e^2 / \hat{\sigma}\_{\hat{\theta}}^2))\\ — so that a
minimum degree of score spread is always preserved. When the unclamped
\\\hat{\rho}^2\\ falls below 0.05 for either group, a warning is issued
and DIF results from that iteration should be interpreted with caution.
This situation typically arises during purification when too few items
remain to yield reliable ability estimates. Users should also be aware
that CATSIB, like its predecessor SIBTEST (Shealy & Stout, 1993), was
originally designed and validated for detecting uniform DIF. Its
statistical behavior under nonuniform or mixed DIF conditions has not
been formally evaluated, and caution is warranted when interpreting
results for items suspected of nonuniform DIF.

## References

Li, H. H., & Stout, W. (1996). A new procedure for detection of crossing
DIF. *Psychometrika, 61*(4), 647-677.

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement*.

Nandakumar, R., & Roussos, L. (2004). Evaluation of the CATSIB DIF
procedure in a pretest setting. *Journal of Educational and Behavioral
Statistics, 29*(2), 177-199.

Shealy, R. T., & Stout, W. F. (1993). A model-based standardization
approach that separates true bias/DIF from group ability differences and
detects test bias/DIF as well as item bias/DIF. *Psychometrika, 58*,
159–194.

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
[est_irt](https://hwangQ.github.io/irtQ/reference/est_irt.md),
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
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = 0.15),
    item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
  )

# Introduce uniform DIF in the focal group by shifting b-parameters
difpar_foc <-
  difpar_ref %>%
  dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))

# Combine the 4 DIF and 36 non-DIF items for both reference and focal groups
# Threfore, the first four items now exhibit uniform DIF
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
#>  EM iteration: 1, Loglike: -25510.5305, Max-Change: 1.315321 EM iteration: 2, Loglike: -22606.9434, Max-Change: 0.354576 EM iteration: 3, Loglike: -22588.3376, Max-Change: 0.156275 EM iteration: 4, Loglike: -22586.3606, Max-Change: 0.075037 EM iteration: 5, Loglike: -22585.8560, Max-Change: 0.039037 EM iteration: 6, Loglike: -22585.6468, Max-Change: 0.022879 EM iteration: 7, Loglike: -22585.5393, Max-Change: 0.014254 EM iteration: 8, Loglike: -22585.4778, Max-Change: 0.009304 EM iteration: 9, Loglike: -22585.4400, Max-Change: 0.006289 EM iteration: 10, Loglike: -22585.4156, Max-Change: 0.004361 EM iteration: 11, Loglike: -22585.3992, Max-Change: 0.00308 EM iteration: 12, Loglike: -22585.3877, Max-Change: 0.002204 EM iteration: 13, Loglike: -22585.3795, Max-Change: 0.001592 EM iteration: 14, Loglike: -22585.3736, Max-Change: 0.001277 EM iteration: 15, Loglike: -22585.3692, Max-Change: 0.001067 EM iteration: 16, Loglike: -22585.3659, Max-Change: 0.000905 EM iteration: 17, Loglike: -22585.3634, Max-Change: 0.000772 EM iteration: 18, Loglike: -22585.3616, Max-Change: 0.000663 EM iteration: 19, Loglike: -22585.3602, Max-Change: 0.000572 EM iteration: 20, Loglike: -22585.3592, Max-Change: 0.000496 EM iteration: 21, Loglike: -22585.3584, Max-Change: 0.000431 EM iteration: 22, Loglike: -22585.3579, Max-Change: 0.000377 EM iteration: 23, Loglike: -22585.3575, Max-Change: 0.00033 EM iteration: 24, Loglike: -22585.3572, Max-Change: 0.00029 EM iteration: 25, Loglike: -22585.3571, Max-Change: 0.000255 EM iteration: 26, Loglike: -22585.3570, Max-Change: 0.000225 EM iteration: 27, Loglike: -22585.3569, Max-Change: 0.000198 EM iteration: 28, Loglike: -22585.3569, Max-Change: 0.000175 EM iteration: 29, Loglike: -22585.3569, Max-Change: 0.000155 EM iteration: 30, Loglike: -22585.3569, Max-Change: 0.000137 EM iteration: 31, Loglike: -22585.3569, Max-Change: 0.000122 EM iteration: 32, Loglike: -22585.3569, Max-Change: 0.000108 EM iteration: 33, Loglike: -22585.3570, Max-Change: 9.6e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 2.41 seconds. 
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

# (a)-1 Compute the CATSIB statistic using provided scores,
#       without purification
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
#>     1, 2, 4, 11, 14, 19 
#>   - Test Statistic: 
#> 
#>         id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   item.1   470   466     936  0.110 0.032  3.462 0.000 ***
#> 2   item.2   470   466     936  0.177 0.028  6.203 0.000 ***
#> 3   item.3   470   466     936  0.054 0.031  1.730 0.084   .
#> 4   item.4   470   466     936  0.181 0.028  6.592 0.000 ***
#> 5   item.5   470   466     936  0.024 0.033  0.736 0.462    
#> 6   item.6   470   466     936 -0.013 0.020 -0.653 0.514    
#> 7   item.7   470   466     936 -0.030 0.032 -0.962 0.336    
#> 8   item.8   470   466     936 -0.036 0.029 -1.225 0.221    
#> 9   item.9   470   466     936  0.001 0.031  0.020 0.984    
#> 10 item.10   470   466     936 -0.052 0.027 -1.891 0.059   .
#> 11 item.11   470   466     936 -0.068 0.031 -2.201 0.028   *
#> 12 item.12   470   466     936 -0.018 0.032 -0.540 0.589    
#> 13 item.13   470   466     936  0.001 0.032  0.038 0.970    
#> 14 item.14   470   466     936 -0.066 0.029 -2.262 0.024   *
#> 15 item.15   470   466     936 -0.029 0.029 -1.021 0.307    
#> 16 item.16   470   466     936 -0.008 0.031 -0.251 0.802    
#> 17 item.17   470   466     936 -0.021 0.030 -0.712 0.477    
#> 18 item.18   470   466     936  0.005 0.029  0.172 0.864    
#> 19 item.19   470   466     936 -0.072 0.028 -2.530 0.011   *
#> 20 item.20   470   466     936 -0.010 0.025 -0.410 0.682    
#> 21 item.21   470   466     936 -0.037 0.029 -1.284 0.199    
#> 22 item.22   470   466     936 -0.051 0.031 -1.635 0.102    
#> 23 item.23   470   466     936 -0.009 0.018 -0.492 0.623    
#> 24 item.24   470   466     936  0.010 0.018  0.569 0.570    
#> 25 item.25   470   466     936 -0.026 0.023 -1.163 0.245    
#> 26 item.26   470   466     936 -0.018 0.029 -0.612 0.540    
#> 27 item.27   470   466     936 -0.030 0.030 -1.008 0.314    
#> 28 item.28   470   466     936 -0.037 0.031 -1.191 0.234    
#> 29 item.29   470   466     936 -0.013 0.025 -0.527 0.598    
#> 30 item.30   470   466     936 -0.022 0.021 -1.044 0.296    
#> 31 item.31   470   466     936 -0.035 0.030 -1.187 0.235    
#> 32 item.32   470   466     936 -0.050 0.026 -1.919 0.055   .
#> 33 item.33   470   466     936  0.025 0.022  1.134 0.257    
#> 34 item.34   470   466     936  0.004 0.031  0.138 0.890    
#> 35 item.35   470   466     936 -0.016 0.032 -0.524 0.600    
#> 36 item.36   470   466     936  0.011 0.023  0.468 0.640    
#> 37 item.37   470   466     936 -0.012 0.021 -0.584 0.559    
#> 38 item.38   470   466     936  0.049 0.032  1.565 0.117    
#> 39 item.39   470   466     936 -0.040 0.029 -1.409 0.159    
#> 40 item.40   470   466     936  0.029 0.031  0.933 0.351    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented. 
#> 

# (a)-2 Compute the CATSIB statistic using provided scores,
#       with purification
dif_2 <- catsib(
  x = est_par, data = data, D = 1, score = score, se = se, group = group, focal.name = 1,
  weight.group = "comb", alpha = 0.05, missing = NA, purify = TRUE
)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 Iteration: 3 Iteration: 4 Iteration: 5 Iteration: 6 
#> Purification is finished. 
print(dif_2)
#> 
#> Call:
#> catsib(x = est_par, data = data, score = score, se = se, group = group, 
#>     focal.name = 1, D = 1, weight.group = "comb", alpha = 0.05, 
#>     missing = NA, purify = TRUE)
#> 
#> DIF analysis using CATSIB method 
#> 
#>  1. Without purification 
#> 
#>   - Potential DIF Items: 
#>     1, 2, 4, 11, 14, 19 
#>   - Test Statistic: 
#> 
#>     id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   V1   470   466     936  0.110 0.032  3.462 0.000 ***
#> 2   V2   470   466     936  0.177 0.028  6.203 0.000 ***
#> 3   V3   470   466     936  0.054 0.031  1.730 0.084   .
#> 4   V4   470   466     936  0.181 0.028  6.592 0.000 ***
#> 5   V5   470   466     936  0.024 0.033  0.736 0.462    
#> 6   V6   470   466     936 -0.013 0.020 -0.653 0.514    
#> 7   V7   470   466     936 -0.030 0.032 -0.962 0.336    
#> 8   V8   470   466     936 -0.036 0.029 -1.225 0.221    
#> 9   V9   470   466     936  0.001 0.031  0.020 0.984    
#> 10 V10   470   466     936 -0.052 0.027 -1.891 0.059   .
#> 11 V11   470   466     936 -0.068 0.031 -2.201 0.028   *
#> 12 V12   470   466     936 -0.018 0.032 -0.540 0.589    
#> 13 V13   470   466     936  0.001 0.032  0.038 0.970    
#> 14 V14   470   466     936 -0.066 0.029 -2.262 0.024   *
#> 15 V15   470   466     936 -0.029 0.029 -1.021 0.307    
#> 16 V16   470   466     936 -0.008 0.031 -0.251 0.802    
#> 17 V17   470   466     936 -0.021 0.030 -0.712 0.477    
#> 18 V18   470   466     936  0.005 0.029  0.172 0.864    
#> 19 V19   470   466     936 -0.072 0.028 -2.530 0.011   *
#> 20 V20   470   466     936 -0.010 0.025 -0.410 0.682    
#> 21 V21   470   466     936 -0.037 0.029 -1.284 0.199    
#> 22 V22   470   466     936 -0.051 0.031 -1.635 0.102    
#> 23 V23   470   466     936 -0.009 0.018 -0.492 0.623    
#> 24 V24   470   466     936  0.010 0.018  0.569 0.570    
#> 25 V25   470   466     936 -0.026 0.023 -1.163 0.245    
#> 26 V26   470   466     936 -0.018 0.029 -0.612 0.540    
#> 27 V27   470   466     936 -0.030 0.030 -1.008 0.314    
#> 28 V28   470   466     936 -0.037 0.031 -1.191 0.234    
#> 29 V29   470   466     936 -0.013 0.025 -0.527 0.598    
#> 30 V30   470   466     936 -0.022 0.021 -1.044 0.296    
#> 31 V31   470   466     936 -0.035 0.030 -1.187 0.235    
#> 32 V32   470   466     936 -0.050 0.026 -1.919 0.055   .
#> 33 V33   470   466     936  0.025 0.022  1.134 0.257    
#> 34 V34   470   466     936  0.004 0.031  0.138 0.890    
#> 35 V35   470   466     936 -0.016 0.032 -0.524 0.600    
#> 36 V36   470   466     936  0.011 0.023  0.468 0.640    
#> 37 V37   470   466     936 -0.012 0.021 -0.584 0.559    
#> 38 V38   470   466     936  0.049 0.032  1.565 0.117    
#> 39 V39   470   466     936 -0.040 0.029 -1.409 0.159    
#> 40 V40   470   466     936  0.029 0.031  0.933 0.351    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 6
#>   - Potential DIF Items: 
#>     1, 2, 3, 4, 11, 38 
#>   - Test Statistic: 
#> 
#>     id n.iter n.ref n.foc n.total   beta    se z.beta     p    
#> 1   V1      2   467   467     934  0.126 0.032  3.881 0.000 ***
#> 2   V2      1   467   463     930  0.187 0.029  6.531 0.000 ***
#> 3   V3      3   476   463     939  0.079 0.031  2.536 0.011   *
#> 4   V4      0   470   466     936  0.181 0.028  6.592 0.000 ***
#> 5   V5      6   467   471     938  0.015 0.033  0.448 0.654    
#> 6   V6      6   467   471     938  0.015 0.019  0.756 0.450    
#> 7   V7      6   467   471     938 -0.003 0.031 -0.097 0.923    
#> 8   V8      6   467   471     938 -0.016 0.030 -0.552 0.581    
#> 9   V9      6   467   471     938  0.029 0.030  0.962 0.336    
#> 10 V10      6   467   471     938 -0.032 0.026 -1.187 0.235    
#> 11 V11      5   470   471     941 -0.063 0.030 -2.078 0.038   *
#> 12 V12      6   467   471     938  0.008 0.032  0.239 0.811    
#> 13 V13      6   467   471     938  0.011 0.032  0.338 0.735    
#> 14 V14      6   467   471     938 -0.054 0.028 -1.908 0.056   .
#> 15 V15      6   467   471     938 -0.017 0.029 -0.598 0.550    
#> 16 V16      6   467   471     938  0.003 0.030  0.102 0.919    
#> 17 V17      6   467   471     938 -0.023 0.029 -0.796 0.426    
#> 18 V18      6   467   471     938  0.031 0.029  1.071 0.284    
#> 19 V19      6   467   471     938 -0.047 0.029 -1.639 0.101    
#> 20 V20      6   467   471     938  0.025 0.025  1.018 0.309    
#> 21 V21      6   467   471     938 -0.004 0.029 -0.135 0.892    
#> 22 V22      6   467   471     938 -0.039 0.032 -1.237 0.216    
#> 23 V23      6   467   471     938  0.001 0.019  0.032 0.975    
#> 24 V24      6   467   471     938  0.013 0.018  0.730 0.465    
#> 25 V25      6   467   471     938 -0.009 0.023 -0.375 0.708    
#> 26 V26      6   467   471     938 -0.003 0.029 -0.102 0.919    
#> 27 V27      6   467   471     938 -0.001 0.029 -0.038 0.970    
#> 28 V28      6   467   471     938 -0.034 0.030 -1.100 0.271    
#> 29 V29      6   467   471     938  0.001 0.025  0.023 0.982    
#> 30 V30      6   467   471     938 -0.001 0.021 -0.068 0.946    
#> 31 V31      6   467   471     938 -0.003 0.030 -0.096 0.924    
#> 32 V32      6   467   471     938 -0.025 0.026 -0.963 0.336    
#> 33 V33      6   467   471     938  0.034 0.022  1.534 0.125    
#> 34 V34      6   467   471     938  0.011 0.032  0.350 0.726    
#> 35 V35      6   467   471     938 -0.012 0.032 -0.381 0.703    
#> 36 V36      6   467   471     938  0.038 0.023  1.619 0.105    
#> 37 V37      6   467   471     938  0.012 0.021  0.578 0.563    
#> 38 V38      4   473   464     937  0.073 0.031  2.328 0.020   *
#> 39 V39      6   467   471     938 -0.031 0.029 -1.060 0.289    
#> 40 V40      6   467   471     938  0.019 0.031  0.604 0.546    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
# }
```
