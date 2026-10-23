# Summary of Item Calibration Results

Summarizes the item calibration results in an object of class `est_irt`,
`est_mg`, or `est_item`, returned by
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

# S3 method for class 'est_irt'
print(x, digits = max(2L, getOption("digits") - 5L), ...)

# S3 method for class 'est_mg'
print(x, digits = max(2L, getOption("digits") - 5L), ...)

# S3 method for class 'est_item'
print(x, digits = max(2L, getOption("digits") - 5L), ...)

# S3 method for class 'summary.est_irt'
print(x, digits = max(2L, getOption("digits") - 5L), ...)

# S3 method for class 'summary.est_mg'
print(x, digits = max(2L, getOption("digits") - 5L), ...)

# S3 method for class 'summary.est_item'
print(x, digits = max(2L, getOption("digits") - 5L), ...)
```

## Arguments

- object:

  An object of class `est_irt`, `est_mg`, or `est_item`.

- ...:

  Additional arguments passed to or from other methods (currently not
  used).

- x:

  An object of class `est_irt`, `est_mg`, or `est_item`, or the summary
  of one of them, to be printed.

- digits:

  For a summary object, the number of decimal places used to round the
  tables of item and group parameter estimates. Not used when an
  `est_irt`, `est_mg`, or `est_item` object itself is printed.

## Value

An object of class `summary.est_irt`, `summary.est_mg`, or
`summary.est_item`: a list of selected components of `object` (the call,
the numbers of items and examinees, the estimation settings and
convergence messages, the log-likelihood, AIC and BIC (not for
`est_item`), and the item and group parameter estimates). Printing the
object displays an overview of the calibration results.

## Details

irtQ exports its own generic `summary()`, which calls the method of the
class of `object` and therefore works as
[`base::summary()`](https://rdrr.io/r/base/summary.html) for other
objects. Printing a summary object rounds the tables of item and group
parameter estimates to the number of decimal places given in the
`digits` argument of the print method (default
`max(2L, getOption("digits") - 5L)`).

## Methods (by class)

- `summary(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

- `summary(est_mg)`: An object created by the function
  [`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md).

- `summary(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

## Functions

- `print(est_irt)`: Print method for an object of class `est_irt` (given
  as `x`). It prints the number of EM cycles, the number of quadrature
  points, the two convergence checks, a note on the computation of the
  variance-covariance matrix, and the log-likelihood.

- `print(est_mg)`: Print method for an object of class `est_mg` (given
  as `x`). It prints the number of EM cycles, the number of quadrature
  points, the two convergence checks, a note on the computation of the
  variance-covariance matrix, and the log-likelihood.

- `print(est_item)`: Print method for an object of class `est_item`
  (given as `x`). It prints the convergence message of the fixed ability
  parameter calibration and the log-likelihood.

- `print(summary.est_irt)`: Print method for the summary of an `est_irt`
  object (given as `x`). It prints the full calibration summary with the
  tables rounded to `digits` decimal places.

- `print(summary.est_mg)`: Print method for the summary of an `est_mg`
  object (given as `x`). It prints the full calibration summary with the
  tables rounded to `digits` decimal places.

- `print(summary.est_item)`: Print method for the summary of an
  `est_item` object (given as `x`). It prints the full calibration
  summary with the tables rounded to `digits` decimal places.

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
# Fit the 1PL model to the LSAT6 data, constraining the slopes to be equal
fit.1pl <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2,
                   fix.a.1pl = FALSE, verbose = FALSE)

# Display the calibration summary
summary(fit.1pl)
#> 
#> Call:
#> est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2, fix.a.1pl = FALSE, 
#>     verbose = FALSE)
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
#>    - Item parameters: 6
#>    - Latent distribution parameters (mean and variance): 0
#>    - Latent distribution weights (empirical histogram): 0
#>  Number of fixed items: 0
#>  Number of E-step cycles completed: 22
#>  Maximum parameter change: 6.905886e-05
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 0.05
#>  Standard error computation: 0
#>  Total computation: 0.05
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
