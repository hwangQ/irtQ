# Extract Components from 'est_irt', 'est_mg', or 'est_item' Objects

Extracts internal components from an object of class `est_irt` (from
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)),
`est_mg` (from
[`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md)), or
`est_item` (from
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md)).

## Usage

``` r
getirt(x, ...)

# S3 method for class 'est_irt'
getirt(x, what, ...)

# S3 method for class 'est_mg'
getirt(x, what, ...)

# S3 method for class 'est_item'
getirt(x, what, ...)
```

## Arguments

- x:

  An object of class `est_irt`, `est_mg`, or `est_item` as returned by
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
  [`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md), or
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md),
  respectively.

- ...:

  Additional arguments passed to or from other methods.

- what:

  A character string specifying the name of the internal component to
  extract.

## Value

The internal component extracted from an object of class `est_irt`,
`est_mg`, or `est_item`, depending on the input to the `x` argument.

## Details

The following components can be extracted from an object of class
`est_irt` created by
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md):

- estimates:

  A data frame containing both the item parameter estimates and their
  corresponding standard errors.

- par.est:

  A data frame containing only the item parameter estimates.

- se.est:

  A data frame containing the standard errors of the item parameter
  estimates, calculated using the cross-product approximation method
  (Meilijson, 1989).

- pos.par:

  A data frame indicating the position index of each estimated item
  parameter. This is useful when interpreting the variance-covariance
  matrix.

- covariance:

  A variance-covariance matrix of the item parameter estimates.

- loglikelihood:

  The total marginal log-likelihood value summed across all items.

- aic:

  Akaike Information Criterion (AIC) based on the marginal
  log-likelihood.

- bic:

  Bayesian Information Criterion (BIC) based on the marginal
  log-likelihood.

- group.par:

  A data frame containing the mean, variance, and standard deviation of
  the latent variable's prior distribution.

- weights:

  A two-column data frame containing quadrature points (first column)
  and corresponding weights (second column) of the (updated) latent
  trait prior.

- posterior.dist:

  A matrix of normalized posterior densities for all response patterns
  at each quadrature point. Rows represent examinees, and columns
  represent quadrature points.

- data:

  A data frame of the examinee response dataset used in estimation.

- scale.D:

  The scaling constant (usually 1 or 1.7) used in the IRT model.

- ncase:

  The number of unique response patterns.

- nitem:

  The number of items included in the dataset.

- Etol:

  The convergence criterion for the E-step of the EM algorithm: the
  largest absolute change in the item parameter estimates between
  consecutive cycles. For FIPC with all items fixed, it applies to the
  change in the mean and variance of the prior distribution.

- MaxE:

  The maximum number of E-steps allowed during EM estimation.

- aprior:

  A list describing the prior distribution for item slope parameters.

- bprior:

  A list describing the prior distribution for item difficulty (or
  threshold) parameters.

- gprior:

  A list describing the prior distribution for item guessing parameters.

- npar.est:

  The total number of parameters estimated.

- niter:

  The number of EM cycles completed.

- maxpar.diff:

  The largest absolute change in the estimates in the last EM cycle.

- EMtime:

  Computation time (in seconds) for the EM algorithm.

- SEtime:

  Computation time (in seconds) for estimating standard errors.

- TotalTime:

  Total computation time (in seconds) for model estimation.

- test.1:

  A message indicating whether the convergence criteria were met (M-step
  convergence and the EM criterion).

- test.2:

  Result of the second-order test indicating whether the information
  matrix was positive definite (a condition for maximum likelihood).

- var.note:

  A note indicating whether the variance-covariance matrix was
  successfully derived from the information matrix.

- fipc:

  Logical value indicating whether Fixed Item Parameter Calibration
  (FIPC) was applied.

- fipc.method:

  The specific method used for FIPC.

- fix.loc:

  An integer vector indicating the positions of fixed items used during
  FIPC.

Components that can be extracted from an object of class `est_mg`
created by
[`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md) include:

- estimates:

  A list with two components: `overall` and `group`.

  - `overall`: A data frame containing item parameter estimates and
    their standard errors, based on the combined data set across all
    groups.

  - `group`: A list of group-specific data frames containing item
    parameter estimates and standard errors for each group.

- par.est:

  Same structure as `estimates`, but containing only the item parameter
  estimates (without standard errors).

- se.est:

  Same structure as `estimates`, but containing only the standard errors
  of the item parameter estimates. The standard errors are computed
  using the cross-product approximation method (Meilijson, 1989).

- pos.par:

  A data frame indicating the position index of each estimated
  parameter. This index is based on the combined item set across all
  groups and is useful when interpreting the variance-covariance matrix.

- covariance:

  A variance-covariance matrix for the item parameter estimates based on
  the combined data from all groups.

- loglikelihood:

  A list with `overall` and `group` components:

  - `overall`: The marginal log-likelihood summed over all unique items
    across all groups.

  - `group`: Group-specific marginal log-likelihood values.

- aic:

  Akaike Information Criterion (AIC) computed from the overall
  log-likelihood.

- bic:

  Bayesian Information Criterion (BIC) computed from the overall
  log-likelihood.

- group.par:

  A list of group-specific summary statistics (mean, variance, and
  standard deviation) of the latent trait prior distribution.

- weights:

  A list of two-column data frames (one per group) containing the
  quadrature points (first column) and the corresponding weights (second
  column) for the updated prior distributions.

- posterior.dist:

  A matrix of normalized posterior densities for all response patterns
  at each quadrature point. Rows correspond to individuals, and columns
  to quadrature points.

- data:

  A list with `overall` and `group` components, each containing examinee
  response data.

- scale.D:

  The scaling constant used in the IRT model (typically 1 or 1.7).

- ncase:

  A list with `overall` and `group` components indicating the number of
  response patterns in each.

- nitem:

  A list with `overall` and `group` components indicating the number of
  items in the respective response sets.

- Etol:

  The convergence criterion for the E-step of the EM algorithm: the
  largest absolute change in the item parameter estimates between
  consecutive cycles. For FIPC with all items fixed, it applies to the
  change in the mean and variance of the prior distribution.

- MaxE:

  Maximum number of E-steps allowed in the EM algorithm.

- aprior:

  A list describing the prior distribution for item slope parameters.

- bprior:

  A list describing the prior distribution for item difficulty
  parameters.

- gprior:

  A list describing the prior distribution for item guessing parameters.

- npar.est:

  Total number of parameters estimated across all unique items.

- niter:

  Number of EM cycles completed.

- maxpar.diff:

  The largest absolute change in the estimates in the last EM cycle.

- EMtime:

  Computation time (in seconds) for EM estimation.

- SEtime:

  Computation time (in seconds) for estimating standard errors.

- TotalTime:

  Total computation time (in seconds) for model estimation.

- test.1:

  A message indicating whether the convergence criteria were met (M-step
  convergence and the EM criterion).

- test.2:

  Second-order condition test result indicating whether the information
  matrix is positive definite.

- var.note:

  A note indicating whether the variance-covariance matrix was
  successfully derived from the information matrix.

- fipc:

  Logical value indicating whether Fixed Item Parameter Calibration
  (FIPC) was used.

- fipc.method:

  The method used for FIPC.

- fix.loc:

  A list with `overall` and `group` components specifying the locations
  of fixed items when FIPC was applied.

Components that can be extracted from an object of class `est_item`
created by
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md)
include:

- estimates:

  A data frame containing both the item parameter estimates and their
  corresponding standard errors.

- par.est:

  A data frame containing only the item parameter estimates.

- se.est:

  A data frame containing the standard errors of the item parameter
  estimates, computed using observed information functions.

- pos.par:

  A data frame indicating the position index of each estimated item
  parameter. This is useful when interpreting the variance-covariance
  matrix.

- covariance:

  A variance-covariance matrix of the item parameter estimates.

- loglikelihood:

  The sum of log-likelihood values across all items in the complete data
  set.

- data:

  A data frame of examinee response data.

- score:

  A numeric vector of examinees' ability values used as fixed effects
  during estimation.

- scale.D:

  The scaling constant (typically 1 or 1.7) used in the IRT model.

- convergence:

  A character string indicating the convergence status of the item
  parameter estimation.

- nitem:

  The total number of items included in the response data.

- deleted.item:

  Items that contained no response data and were excluded from
  estimation.

- npar.est:

  The total number of estimated item parameters.

- n.response:

  An integer vector indicating the number of responses used to estimate
  parameters for each item.

- TotalTime:

  Total computation time (in seconds) for the estimation process.

See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md), and
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md) for
more details.

## Methods (by class)

- `getirt(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `getirt(est_mg)`: An object created by the function
  [`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md).

- `getirt(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
# Fit a 2PL model to the LSAT6 data
mod.2pl <- est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2)
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -3182.3860, Max-Change: 2.294233 EM iteration: 2, Loglike: -2561.0478, Max-Change: 0.58517 EM iteration: 3, Loglike: -2482.8950, Max-Change: 0.323232 EM iteration: 4, Loglike: -2469.4379, Max-Change: 0.180623 EM iteration: 5, Loglike: -2467.2706, Max-Change: 0.106034 EM iteration: 6, Loglike: -2466.8615, Max-Change: 0.066546 EM iteration: 7, Loglike: -2466.7487, Max-Change: 0.044361 EM iteration: 8, Loglike: -2466.7049, Max-Change: 0.030996 EM iteration: 9, Loglike: -2466.6843, Max-Change: 0.022467 EM iteration: 10, Loglike: -2466.6736, Max-Change: 0.016791 EM iteration: 11, Loglike: -2466.6674, Max-Change: 0.012897 EM iteration: 12, Loglike: -2466.6636, Max-Change: 0.010155 EM iteration: 13, Loglike: -2466.6611, Max-Change: 0.008179 EM iteration: 14, Loglike: -2466.6593, Max-Change: 0.00672 EM iteration: 15, Loglike: -2466.6580, Max-Change: 0.005616 EM iteration: 16, Loglike: -2466.6571, Max-Change: 0.00476 EM iteration: 17, Loglike: -2466.6563, Max-Change: 0.004082 EM iteration: 18, Loglike: -2466.6557, Max-Change: 0.003532 EM iteration: 19, Loglike: -2466.6553, Max-Change: 0.003079 EM iteration: 20, Loglike: -2466.6549, Max-Change: 0.00270 EM iteration: 21, Loglike: -2466.6546, Max-Change: 0.002378 EM iteration: 22, Loglike: -2466.6544, Max-Change: 0.002102 EM iteration: 23, Loglike: -2466.6542, Max-Change: 0.001863 EM iteration: 24, Loglike: -2466.6540, Max-Change: 0.001655 EM iteration: 25, Loglike: -2466.6539, Max-Change: 0.001473 EM iteration: 26, Loglike: -2466.6538, Max-Change: 0.001313 EM iteration: 27, Loglike: -2466.6537, Max-Change: 0.001172 EM iteration: 28, Loglike: -2466.6537, Max-Change: 0.001047 EM iteration: 29, Loglike: -2466.6536, Max-Change: 0.000936 EM iteration: 30, Loglike: -2466.6536, Max-Change: 0.000838 EM iteration: 31, Loglike: -2466.6535, Max-Change: 0.000751 EM iteration: 32, Loglike: -2466.6535, Max-Change: 0.000673 EM iteration: 33, Loglike: -2466.6535, Max-Change: 0.000603 EM iteration: 34, Loglike: -2466.6535, Max-Change: 0.000541 EM iteration: 35, Loglike: -2466.6534, Max-Change: 0.000486 EM iteration: 36, Loglike: -2466.6534, Max-Change: 0.000436 EM iteration: 37, Loglike: -2466.6534, Max-Change: 0.000392 EM iteration: 38, Loglike: -2466.6534, Max-Change: 0.000352 EM iteration: 39, Loglike: -2466.6534, Max-Change: 0.000316 EM iteration: 40, Loglike: -2466.6534, Max-Change: 0.000284 EM iteration: 41, Loglike: -2466.6534, Max-Change: 0.000256 EM iteration: 42, Loglike: -2466.6534, Max-Change: 0.00023 EM iteration: 43, Loglike: -2466.6534, Max-Change: 0.000207 EM iteration: 44, Loglike: -2466.6534, Max-Change: 0.000186 EM iteration: 45, Loglike: -2466.6534, Max-Change: 0.000167 EM iteration: 46, Loglike: -2466.6534, Max-Change: 0.000151 EM iteration: 47, Loglike: -2466.6534, Max-Change: 0.000136 EM iteration: 48, Loglike: -2466.6534, Max-Change: 0.000122 EM iteration: 49, Loglike: -2466.6534, Max-Change: 0.00011 EM iteration: 50, Loglike: -2466.6534, Max-Change: 9.9e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 0.36 seconds. 

# Extract item parameter estimates
(est.par <- getirt(mod.2pl, what = "par.est"))
#>   id cats model     par.1      par.2 par.3
#> 1 V1    2  2PLM 0.8255888 -3.3590452    NA
#> 2 V2    2  2PLM 0.7229138 -1.3697963    NA
#> 3 V3    2  2PLM 0.8904047 -0.2797759    NA
#> 4 V4    2  2PLM 0.6884844 -1.8661188    NA
#> 5 V5    2  2PLM 0.6570762 -3.1250079    NA

# Extract standard error estimates
(est.se <- getirt(mod.2pl, what = "se.est"))
#>   id cats model     par.1      par.2 par.3
#> 1 V1    2  2PLM 0.2538702 0.85483653    NA
#> 2 V2    2  2PLM 0.1879120 0.30935133    NA
#> 3 V3    2  2PLM 0.2298908 0.09951049    NA
#> 4 V4    2  2PLM 0.1862128 0.43815986    NA
#> 5 V5    2  2PLM 0.2027270 0.84271138    NA

# Extract the variance-covariance matrix of item parameter estimates
(cov.mat <- getirt(mod.2pl, what = "covariance"))
#>                [,1]         [,2]          [,3]          [,4]         [,5]
#>  [1,]  0.0644500663  0.213516110 -0.0059296609 -1.012578e-02 -0.004681038
#>  [2,]  0.2135161095  0.730745496 -0.0205129366 -3.393420e-02 -0.016087744
#>  [3,] -0.0059296609 -0.020512937  0.0353109379  5.459288e-02 -0.011084779
#>  [4,] -0.0101257773 -0.033934198  0.0545928849  9.569824e-02 -0.016831654
#>  [5,] -0.0046810378 -0.016087744 -0.0110847789 -1.683165e-02  0.052849776
#>  [6,] -0.0011271119 -0.003110929 -0.0023765393 -2.603282e-03  0.012397532
#>  [7,] -0.0031463866 -0.010012688 -0.0005355098  4.910666e-06 -0.011308992
#>  [8,] -0.0075724067 -0.022812226  0.0001639620  3.457859e-03 -0.026321355
#>  [9,]  0.0002965241  0.002494229 -0.0037729160 -5.976366e-03 -0.004409260
#> [10,]  0.0023222383  0.015568907 -0.0160475746 -2.471635e-02 -0.015974475
#>               [,6]          [,7]         [,8]          [,9]        [,10]
#>  [1,] -0.001127112 -3.146387e-03 -0.007572407  0.0002965241  0.002322238
#>  [2,] -0.003110929 -1.001269e-02 -0.022812226  0.0024942292  0.015568907
#>  [3,] -0.002376539 -5.355098e-04  0.000163962 -0.0037729160 -0.016047575
#>  [4,] -0.002603282  4.910666e-06  0.003457859 -0.0059763663 -0.024716348
#>  [5,]  0.012397532 -1.130899e-02 -0.026321355 -0.0044092598 -0.015974475
#>  [6,]  0.009902337 -2.893417e-03 -0.005765609 -0.0010248799 -0.002351144
#>  [7,] -0.002893417  3.467521e-02  0.078608493 -0.0039184840 -0.017434966
#>  [8,] -0.005765609  7.860849e-02  0.191984062 -0.0094264375 -0.041401179
#>  [9,] -0.001024880 -3.918484e-03 -0.009426438  0.0410982364  0.168071445
#> [10,] -0.002351144 -1.743497e-02 -0.041401179  0.1680714454  0.710162472
# }
```
