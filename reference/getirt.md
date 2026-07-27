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

  The convergence criterion used for the E-step in the EM algorithm.

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

  The maximum change in parameter estimates at convergence.

- EMtime:

  Computation time (in seconds) for the EM algorithm.

- SEtime:

  Computation time (in seconds) for estimating standard errors.

- TotalTime:

  Total computation time (in seconds) for model estimation.

- test.1:

  Result of the first-order test indicating whether the gradients were
  sufficiently close to zero.

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

  Convergence criterion used for the E-step in the EM algorithm.

- MaxE:

  Maximum number of E-steps allowed in the EM algorithm.

- aprior:

  A list describing the prior distribution for item slope parameters.

- gprior:

  A list describing the prior distribution for item guessing parameters.

- npar.est:

  Total number of parameters estimated across all unique items.

- niter:

  Number of EM cycles completed.

- maxpar.diff:

  Maximum change in item parameter estimates at convergence.

- EMtime:

  Computation time (in seconds) for EM estimation.

- SEtime:

  Computation time (in seconds) for estimating standard errors.

- TotalTime:

  Total computation time (in seconds) for model estimation.

- test.1:

  First-order condition test result indicating whether gradients
  converged sufficiently.

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
#>  EM iteration: 1, Loglike: -3182.3860, Max-Change: 0.390402 EM iteration: 2, Loglike: -2561.0478, Max-Change: 0.00000 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 0.04 seconds. 

# Extract item parameter estimates
(est.par <- getirt(mod.2pl, what = "par.est"))
#>   id cats model     par.1        par.2 par.3
#> 1 V1    2  2PLM 0.8851380 -2.872384950    NA
#> 2 V2    2  2PLM 0.8780243 -0.890284439    NA
#> 3 V3    2  2PLM 0.9305760  0.003668382    NA
#> 4 V4    2  2PLM 0.8617347 -1.269897767    NA
#> 5 V5    2  2PLM 0.8410096 -2.250964707    NA

# Extract standard error estimates
(est.se <- getirt(mod.2pl, what = "se.est"))
#>   id cats model     par.1      par.2 par.3
#> 1 V1    2  2PLM 0.2452172 0.67487124    NA
#> 2 V2    2  2PLM 0.1905223 0.18450934    NA
#> 3 V3    2  2PLM 0.2020198 0.08146577    NA
#> 4 V4    2  2PLM 0.1919434 0.25426752    NA
#> 5 V5    2  2PLM 0.2087045 0.48313552    NA

# Extract the variance-covariance matrix of item parameter estimates
(cov.mat <- getirt(mod.2pl, what = "covariance"))
#>                [,1]         [,2]          [,3]          [,4]         [,5]
#>  [1,]  0.0601314834  0.161626270 -0.0051257769 -0.0044224511 -0.002071661
#>  [2,]  0.1616262702  0.455451195 -0.0136565487 -0.0103462255 -0.005500742
#>  [3,] -0.0051257769 -0.013656549  0.0362987520  0.0305074017 -0.007660103
#>  [4,] -0.0044224511 -0.010346225  0.0305074017  0.0340436966 -0.005426263
#>  [5,] -0.0020716614 -0.005500742 -0.0076601034 -0.0054262631  0.040812001
#>  [6,]  0.0003284084  0.001904617  0.0005559771  0.0017297229  0.001821232
#>  [7,] -0.0019311592 -0.003705538 -0.0009465802  0.0006292671 -0.008775219
#>  [8,] -0.0020030971 -0.001814204  0.0008034353  0.0040386727 -0.010316758
#>  [9,]  0.0015223739  0.006454168 -0.0032010354 -0.0020718830 -0.003243842
#> [10,]  0.0049990824  0.020820838 -0.0064320092 -0.0026824699 -0.005159683
#>                [,6]          [,7]          [,8]          [,9]        [,10]
#>  [1,]  3.284084e-04 -1.931159e-03 -0.0020030971  0.0015223739  0.004999082
#>  [2,]  1.904617e-03 -3.705538e-03 -0.0018142036  0.0064541684  0.020820838
#>  [3,]  5.559771e-04 -9.465802e-04  0.0008034353 -0.0032010354 -0.006432009
#>  [4,]  1.729723e-03  6.292671e-04  0.0040386727 -0.0020718830 -0.002682470
#>  [5,]  1.821232e-03 -8.775219e-03 -0.0103167581 -0.0032438424 -0.005159683
#>  [6,]  6.636672e-03 -5.923697e-05  0.0011872961  0.0003003195  0.002347453
#>  [7,] -5.923697e-05  3.684226e-02  0.0450056386 -0.0031576436 -0.007049767
#>  [8,]  1.187296e-03  4.500564e-02  0.0646519739 -0.0034621648 -0.006508739
#>  [9,]  3.003195e-04 -3.157644e-03 -0.0034621648  0.0435575616  0.097528414
#> [10,]  2.347453e-03 -7.049767e-03 -0.0065087389  0.0975284136  0.233419935
# }
```
