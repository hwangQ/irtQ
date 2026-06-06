# Asymptotic Variance-Covariance Matrices of Item Parameter Estimates

This function computes the analytical asymptotic variance-covariance
matrices of item parameter estimates for dichotomous and polytomous IRT
models, without requiring examinee response data. Given a set of item
parameter estimates and the corresponding sample sizes, the function
derives the matrices using analytical formulas (e.g., Li & Lissitz,
2004; Thissen & Wainer, 1982). The square roots of the diagonal elements
(variances) provide the asymptotic standard errors of the maximum
likelihood estimates.

## Usage

``` r
covirt(
  x,
  D = 1,
  nstd = 1000,
  pcm.loc = NULL,
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL
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

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- nstd:

  An integer or a vector of integers indicating the sample size(s). If a
  vector is provided, its length must match the number of items in the
  `x` argument. Default is 1,000. See **Details**.

- pcm.loc:

  A vector of integers indicating the positions of items calibrated
  using the partial credit model (PCM). For PCM items, the
  variance-covariance matrices are computed only for the item category
  difficulty parameters. Default is `NULL`. See the **Details** for more
  information.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. These values are used to
  generate the Gaussian quadrature points and weights when
  `weights = NULL`. Default is `c(0, 1)`.

- nquad:

  An integer indicating the number of Gaussian quadrature points to be
  generated from the normal prior distribution. The specified value is
  used when `weights` is not `NULL`. Default is 41.

- weights:

  An optional two-column data frame or matrix where the first column is
  the quadrature points (nodes) and the second column is the
  corresponding weights. This is typically used in quadrature-based IRT
  analysis.

## Value

A named list with the following two components:

- `cov`: A named list of variance-covariance matrices for item parameter
  estimates. Each element corresponds to a single item and contains a
  square matrix whose dimensions match the number of estimated
  parameters for that item. For dichotomous items, this typically
  includes slopes and intercepts. For polytomous items, it includes
  category difficulty parameters (for PCM) or both slope and difficulty
  (or threshold) parameters (for GRM and GPCM).

- `se`: A named list of vectors containing the asymptotic standard
  errors (ASEs) of the item parameter estimates, computed as the square
  roots of the diagonal elements of each corresponding
  variance-covariance matrix in `cov`.

The names of the list elements in both `cov` and `se` correspond to the
item identifiers (e.g., item names or labels) as given in the first
column of the input `x`.

## Details

The standard errors obtained from this analytical approach are generally
considered lower bounds of the true standard errors (Thissen & Wainer,
1982). Thus, they may serve as useful approximations for evaluating the
precision of item parameter estimates when empirical standard errors are
not reported in the literature or research reports.

If the item parameters provided in the `x` argument were calibrated
using different sample sizes, a corresponding vector of sample sizes
must be specified via the `nstd` argument. For example, suppose you wish
to compute the variance-covariance matrices of five 3PLM items that were
calibrated using 500, 600, 1,000, 2,000, and 700 examinees,
respectively. In this case, set `nstd = c(500, 600, 1000, 2000, 700)`.

Since the item metadata allows only `"GPCM"` to denote both the partial
credit model (PCM) and the generalized partial credit model (GPCM), PCM
items must be explicitly identified using the `pcm.loc` argument. This
is necessary because the category difficulty parameters of PCM items
require special handling when computing variance-covariance matrices.
For instance, if you wish to compute the matrices for five polytomous
items and the last two were calibrated using PCM, then specify
`pcm.loc = c(4, 5)`.

## References

Li, Y. & Lissitz, R. (2004). Applications of the analytically derived
asymptotic standard errors of item response theory item parameter
estimates. *Journal of educational measurement, 41*(2), 85-117.

Thissen, D. & Wainer, H. (1982). Some standard errors in item response
theory. *Psychometrika, 47*, 397-412.

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# Example using a "-prm.txt" file exported from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select the first two dichotomous items and the last polytomous item
x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df[c(1:2, 55), ]

# Compute the variance-covariance matrices assuming a sample size of 2,000
covirt(x, D = 1, nstd = 2000, norm.prior = c(0, 1), nquad = 41)
#> $cov
#> $cov$CMC1
#>            par.1      par.2       par.3
#> par.1 0.04217319 0.02741566 0.012644892
#> par.2 0.02741566 0.06741859 0.014367975
#> par.3 0.01264489 0.01436798 0.004767954
#> 
#> $cov$CMC2
#>            par.1       par.2       par.3
#> par.1 0.02778076 0.016981307 0.008625830
#> par.2 0.01698131 0.015693052 0.008648866
#> par.3 0.00862583 0.008648866 0.005782880
#> 
#> $cov$AFR3
#>               par.1         par.2         par.3         par.4         par.5
#> par.1  1.453330e-03  0.0009384134 -6.872098e-06 -0.0008555713 -0.0015951007
#> par.2  9.384134e-04  0.0027747596  1.456908e-03  0.0005107048 -0.0002225561
#> par.3 -6.872098e-06  0.0014569084  2.206776e-03  0.0016094777  0.0012267762
#> par.4 -8.555713e-04  0.0005107048  1.609478e-03  0.0028580848  0.0027270730
#> par.5 -1.595101e-03 -0.0002225561  1.226776e-03  0.0027270730  0.0043192115
#> 
#> 
#> $se
#> $se$CMC1
#>      par.1      par.2      par.3 
#> 0.20536112 0.25965090 0.06905037 
#> 
#> $se$CMC2
#>      par.1      par.2      par.3 
#> 0.16667561 0.12527191 0.07604525 
#> 
#> $se$AFR3
#>      par.1      par.2      par.3      par.4      par.5 
#> 0.03812256 0.05267599 0.04697633 0.05346106 0.06572071 
#> 
#> 
```
