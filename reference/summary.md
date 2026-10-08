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
#>  EM iteration: 1, Loglike: -3182.3860, Max-Change: 2.293929 EM iteration: 2, Loglike: -2561.3380, Max-Change: 0.58111 EM iteration: 3, Loglike: -2483.1811, Max-Change: 0.31473 EM iteration: 4, Loglike: -2469.6884, Max-Change: 0.171175 EM iteration: 5, Loglike: -2467.5148, Max-Change: 0.096225 EM iteration: 6, Loglike: -2467.1096, Max-Change: 0.056965 EM iteration: 7, Loglike: -2467.0029, Max-Change: 0.035311 EM iteration: 8, Loglike: -2466.9648, Max-Change: 0.022579 EM iteration: 9, Loglike: -2466.9493, Max-Change: 0.014699 EM iteration: 10, Loglike: -2466.9427, Max-Change: 0.009659 EM iteration: 11, Loglike: -2466.9398, Max-Change: 0.006377 EM iteration: 12, Loglike: -2466.9386, Max-Change: 0.004219 EM iteration: 13, Loglike: -2466.9380, Max-Change: 0.002795 EM iteration: 14, Loglike: -2466.9378, Max-Change: 0.001852 EM iteration: 15, Loglike: -2466.9377, Max-Change: 0.001228 EM iteration: 16, Loglike: -2466.9376, Max-Change: 0.000814 EM iteration: 17, Loglike: -2466.9376, Max-Change: 0.000539 EM iteration: 18, Loglike: -2466.9376, Max-Change: 0.000358 EM iteration: 19, Loglike: -2466.9376, Max-Change: 0.000237 EM iteration: 20, Loglike: -2466.9376, Max-Change: 0.000157 EM iteration: 21, Loglike: -2466.9376, Max-Change: 0.000104 EM iteration: 22, Loglike: -2466.9376, Max-Change: 6.9e-05 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 0.1 seconds. 

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
#>  Number of E-step cycles completed: 22
#>  Maximum parameter change: 6.905886e-05
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 0.08
#>  Standard error computation: 0
#>  Total computation: 0.1
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 4933.875
#>  Akaike Information Criterion (AIC): 4945.875
#>  Bayesian Information Criterion (BIC): 4975.322
#>  Item Parameters: 
#>    id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3
#> 1  V1     2   1PLM   0.76  0.07  -3.62  0.32     NA    NA
#> 2  V2     2   1PLM   0.76    NA  -1.32  0.14     NA    NA
#> 3  V3     2   1PLM   0.76    NA  -0.32  0.10     NA    NA
#> 4  V4     2   1PLM   0.76    NA  -1.73  0.17     NA    NA
#> 5  V5     2   1PLM   0.76    NA  -2.78  0.25     NA    NA
#>  Group Parameters: 
#>            mu  sigma2  sigma
#> estimates   0       1      1
#> se         NA      NA     NA
#> 
# }
```
