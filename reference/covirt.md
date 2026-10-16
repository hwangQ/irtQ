# Asymptotic Variance-Covariance Matrices of Item Parameter Estimates

This function computes the analytical asymptotic variance-covariance
matrices of item parameter estimates for dichotomous and polytomous IRT
models without examinee response data. For each item, the matrix is the
inverse of the expected (Fisher) information of the item parameters,
which is integrated over an assumed ability distribution and multiplied
by the calibration sample size (Thissen & Wainer, 1982; Li & Lissitz,
2004). The square roots of the diagonal elements are the asymptotic
standard errors of the maximum likelihood estimates.

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

  A data frame containing item metadata (item IDs, numbers of score
  categories, IRT models, and item parameters). The item parameters are
  treated as the true values at which the information is evaluated; for
  1PLM items and for PCM items specified in `pcm.loc`, the slope in `x`
  is held fixed. See
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  the format of the item metadata, which can be created with
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md).

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- nstd:

  An integer, or a vector of integers with one value per item in `x`,
  giving the number of examinees used to calibrate the items. The
  information of each item is multiplied by its sample size, so the
  standard errors are proportional to 1 / sqrt(`nstd`). Default is 1000.
  See **Details**.

- pcm.loc:

  A vector of integers giving the positions (rows of `x`) of GPCM items
  that were calibrated with the partial credit model (PCM). For these
  items, the slope is held at the value in `x` and the
  variance-covariance matrices are computed only for the threshold
  parameters. Positions of items that are not GPCM items are ignored.
  Default is `NULL`. See **Details**.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal ability distribution over which the
  information is integrated. Used to generate the Gaussian quadrature
  points and weights when `weights = NULL`. Default is `c(0, 1)`.

- nquad:

  An integer specifying the number of Gaussian quadrature points
  generated from the normal distribution given in `norm.prior`. Used
  only when `weights = NULL`. Default is 41.

- weights:

  An optional two-column data frame or matrix of quadrature points
  (first column) and weights (second column) that defines the ability
  distribution. The weights should sum to 1, because they are not
  normalized. If `weights` is supplied, `norm.prior` and `nquad` are
  ignored. Such a data frame can be created with
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).
  Default is `NULL`.

## Value

A named list with the following two components:

- `cov`: A named list of variance-covariance matrices, one per item.
  Each matrix contains only the estimated parameters of the item, and
  its row and column names are the parameter columns of `x`: `par.2`
  (difficulty) for the 1PLM; `par.1` and `par.2` (slope and difficulty)
  for the 2PLM; `par.1` to `par.3` (slope, difficulty, and guessing) for
  the 3PLM; `par.1` to `par.K` (slope and K - 1 thresholds) for GRM and
  GPCM items with K score categories; and `par.2` to `par.K` (thresholds
  only) for PCM items specified in `pcm.loc`. Covariances between
  parameters of different items are not computed. If the information
  matrix of an item is singular, a warning is issued and the matrix is
  filled with `NA`.

- `se`: A named list of vectors containing the asymptotic standard
  errors (ASEs) of the item parameter estimates, computed as the square
  roots of the diagonal elements of each corresponding
  variance-covariance matrix in `cov`.

The names of the list elements in both `cov` and `se` correspond to the
item identifiers (e.g., item names or labels) as given in the first
column of the input `x`.

## Details

For each item, the asymptotic variance-covariance matrix is the inverse
of the expected (Fisher) information matrix of the item parameters,
which is integrated over the ability distribution with quadrature and
multiplied by the sample size (Thissen & Wainer, 1982; Li & Lissitz,
2004). For an item with \\K\\ score categories, the element of the
information matrix for a pair of item parameters \\\xi_s\\ and \\\xi_t\\
is

\$\$I(\xi_s, \xi_t) = N \sum\_{q=1}^{Q} A(X_q) \sum\_{k=1}^{K}
\frac{1}{P_k(X_q)} \frac{\partial P_k(X_q)}{\partial \xi_s}
\frac{\partial P_k(X_q)}{\partial \xi_t},\$\$

where \\N\\ is the sample size (`nstd`), \\X_q\\ and \\A(X_q)\\ are the
\\Q\\ quadrature points and weights that represent the ability
distribution (`norm.prior` and `nquad`, or `weights`), and \\P_k(X_q)\\
is the probability of a response in category \\k\\ at \\X_q\\. This is
the information of a categorical (multinomial) response. For a
dichotomous item, it reduces to the familiar form in which each term is
divided by \\P(X_q) \[1 - P(X_q)\]\\ (Li & Lissitz, 2004, Equation A4).

For polytomous items, the information above differs from Equation A8 of
Li and Lissitz (2004). Their Equation A6 writes the likelihood of a
response (and Equation A7 its logarithm) as \\\prod\_{k} P_k^{u_k} (1 -
P_k)^{1 - u_k}\\, where \\u_k\\ is 1 for the observed category and 0
otherwise. The \\K\\ category indicators are thus treated as independent
Bernoulli variables, and Equation A8 accordingly divides each term by
\\P_k (1 - P_k)\\ instead of \\P_k\\. However, the indicators of one
response sum to 1, so the response follows a multinomial distribution
with likelihood \\\prod\_{k} P_k^{u_k}\\, whose information divides each
term by \\P_k\\ only. Because \\1 - P_k \< 1\\, Equation A8 overstates
the information and understates the standard errors. For example, for an
item with two categories, Equation A8 gives exactly twice the
information of their Equation A4 for dichotomous items. This function
uses the multinomial information.

The multinomial information was checked by a simulation that matches the
setting of this function: abilities drawn from N(0, 1) and treated as
known, maximum likelihood estimation, a sample size of 1,000, and 1,000
replications for GRM and GPCM items with three to five categories. The
standard errors from the multinomial information agreed with the
standard deviations of the estimates across replications (ratios of 0.96
to 1.08), whereas those from Equation A8 were about 10 to 23 percent
smaller (ratios of 0.77 to 0.90). For a dichotomous (2PLM) item, the two
forms give the same results.

The matrices are derived under the assumptions that the IRT model fits
the data, that the item parameters are estimated by maximum likelihood,
and that the ability distribution is known (Li & Lissitz, 2004). Because
the abilities are treated as known, the uncertainty in the abilities is
ignored, and prior distributions on the item parameters are not taken
into account. The resulting standard errors are therefore generally
considered lower bounds of the standard errors of an actual calibration
(Thissen & Wainer, 1982). They may serve as approximations for
evaluating the precision of item parameter estimates when empirical
standard errors are not reported.

If the item parameters provided in the `x` argument were calibrated
using different sample sizes, a corresponding vector of sample sizes
must be specified via the `nstd` argument. For example, suppose you wish
to compute the variance-covariance matrices of five 3PLM items that were
calibrated using 500, 600, 1,000, 2,000, and 700 examinees,
respectively. In this case, set `nstd = c(500, 600, 1000, 2000, 700)`.

Since the item metadata allows only `"GPCM"` to denote both the partial
credit model (PCM) and the generalized partial credit model (GPCM), PCM
items must be explicitly identified using the `pcm.loc` argument. For
these items, the slope is not treated as an estimated parameter: it is
held at the value in `x`, and only the threshold parameters enter the
information matrix. For instance, if the last two of five polytomous
items were calibrated using the PCM, specify `pcm.loc = c(4, 5)`. A GRM
or GPCM item with two score categories is treated as a 2PLM item, and a
PCM item with two score categories as a 1PLM item with a fixed slope.

## References

Li, Y. H., & Lissitz, R. W. (2004). Applications of the analytically
derived asymptotic standard errors of item response theory item
parameter estimates. *Journal of Educational Measurement, 41*(2),
85-117.
[doi:10.1111/j.1745-3984.2004.tb01109.x](https://doi.org/10.1111/j.1745-3984.2004.tb01109.x)
.

Thissen, D., & Wainer, H. (1982). Some standard errors in item response
theory. *Psychometrika, 47*(4), 397-412.
[doi:10.1007/BF02293705](https://doi.org/10.1007/BF02293705) .

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
#>               par.1         par.2        par.3         par.4         par.5
#> par.1  2.223338e-03  0.0016396610 4.901224e-05 -0.0013690094 -0.0026195815
#> par.2  1.639661e-03  0.0044055650 2.225489e-03  0.0006531171 -0.0005557895
#> par.3  4.901224e-05  0.0022254891 2.951196e-03  0.0022103306  0.0017965206
#> par.4 -1.369009e-03  0.0006531171 2.210331e-03  0.0039665170  0.0041980802
#> par.5 -2.619582e-03 -0.0005557895 1.796521e-03  0.0041980802  0.0067007153
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
#> 0.04715229 0.06637443 0.05432491 0.06298029 0.08185790 
#> 
#> 
```
