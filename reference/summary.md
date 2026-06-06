# Summary of Item Calibration Results

This S3 method summarizes the IRT calibration results from an object of
class `est_irt`, `est_mg`, or `est_item`, which are returned by the
functions
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md), and
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md),
respectively.

## Usage

``` r
summary(object, ...)

# S3 method for class 'est_irt'
summary(object, ...)

# S3 method for class 'est_mg'
summary(object, ...)

# S3 method for class 'est_item'
summary(object, ...)
```

## Arguments

- object:

  An object of class `est_irt`, `est_mg`, or `est_item`.

- ...:

  Additional arguments passed to or from other methods (currently not
  used).

## Value

A list of internal components extracted from the given object. In
addition, the summary method prints an overview of the IRT calibration
results to the console.

## Methods (by class)

- `summary(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `summary(est_mg)`: An object created by the function
  [`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md).

- `summary(est_item)`: An object created by the function
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
# Fit the 1PL model to LSAT6 data and constrain the slope parameters to be equal
fit.1pl <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2, fix.a.1pl = FALSE)
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -3182.3860, Max-Change: 0.385069 EM iteration: 2, Loglike: -2561.3380, Max-Change: 0.00000 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 0.03 seconds. 

# Display the calibration summary
summary(fit.1pl)
#> 
#> Call:
#> est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2, fix.a.1pl = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items: 5
#>  Number of Cases: 1000
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 1e-04
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 6
#>  Number of fixed items: 0
#>  Number of E-step cycles completed: 2
#>  Maximum parameter change: 0
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 0.02
#>  Standard error computation: 0
#>  Total computation: 0.03
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 4966.362
#>  Akaike Information Criterion (AIC): 4978.362
#>  Bayesian Information Criterion (BIC): 5007.809
#>  Item Parameters: 
#>    id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3
#> 1  V1     2   1PLM   0.88  0.07  -2.88  0.25     NA    NA
#> 2  V2     2   1PLM   0.88    NA  -0.88  0.11     NA    NA
#> 3  V3     2   1PLM   0.88    NA  -0.01  0.08     NA    NA
#> 4  V4     2   1PLM   0.88    NA  -1.24  0.13     NA    NA
#> 5  V5     2   1PLM   0.88    NA  -2.15  0.20     NA    NA
#>  Group Parameters: 
#>            mu  sigma2  sigma
#> estimates   0       1      1
#> se         NA      NA     NA
#> 
# }
```
