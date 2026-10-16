# Polytomous Response Model (PRM) Probabilities (GRM and GPCM)

This function computes the probability of selecting each response
category for an item, given a set of theta values, using the graded
response model (GRM) or the (generalized) partial credit model (GPCM).

## Usage

``` r
prm(theta, a, d, D = 1, pr.model = c("GRM", "GPCM"))
```

## Arguments

- theta:

  A numeric vector of ability values (latent traits).

- a:

  A numeric value of the item discrimination (slope) parameter.

- d:

  A numeric vector of the threshold parameters of the item. For the GRM,
  the thresholds must be in increasing order.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- pr.model:

  A character string specifying the polytomous IRT model: `"GRM"` for
  the graded response model or `"GPCM"` for the (generalized) partial
  credit model. Default is `"GRM"`.

## Value

A matrix of category response probabilities, where each row corresponds
to a theta value and each column represents a score category of the
item.

## Details

When computing category probabilities using the partial credit model
(PCM), set `a = 1`.

For `pr.model = "GPCM"`, `d` contains the threshold parameters \\b_v =
\beta - \tau_v\\, that is, the overall item location minus the threshold
of each score category. An item with *K* score categories requires *K -
1* threshold parameters. The term for the lowest score category is fixed
at 0 and is not supplied. For example, a GPCM item with five categories
requires four threshold parameters.

The GRM probabilities are bounded to the range from 1e-10 to 1 - 1e-10,
whereas the GPCM probabilities are bounded below only by the smallest
positive normalized double value (about 2.2e-308). The category
probabilities in
[`traceline()`](https://hwangQ.github.io/irtQ/reference/traceline.md)
and [`info()`](https://hwangQ.github.io/irtQ/reference/info.md) are
bounded below by 1e-10 for both models.

For more details on the parameterization of the (generalized) partial
credit model, refer to the *IRT Models* section in the
[irtQ-package](https://hwangQ.github.io/irtQ/reference/irtQ-package.md)
documentation.

## See also

[`drm()`](https://hwangQ.github.io/irtQ/reference/drm.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Category probabilities for an item with four categories
## using the generalized partial credit model (GPCM)
prm(theta = c(-0.2, 0, 0.5), a = 1.4, d = c(-0.2, 0, 0.5), D = 1, pr.model = "GPCM")
#>            [,1]      [,2]      [,3]       [,4]
#> [1,] 0.32900822 0.3290082 0.2486591 0.09332451
#> [2,] 0.23237945 0.3074682 0.3074682 0.15268418
#> [3,] 0.06946582 0.1850886 0.3727228 0.37272276

## Category probabilities for an item with five categories
## using the graded response model (GRM)
prm(theta = c(-0.2, 0, 0.5), a = 1.2, d = c(-0.4, -0.2, 0.4, 1.5), D = 1, pr.model = "GRM")
#>           [,1]       [,2]      [,3]      [,4]      [,5]
#> [1,] 0.4402864 0.05971365 0.1726070 0.2123263 0.1150667
#> [2,] 0.3822521 0.05803423 0.1774615 0.2404011 0.1418511
#> [3,] 0.2535060 0.04802877 0.1685012 0.2984888 0.2314752
```
