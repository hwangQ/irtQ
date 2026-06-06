# Dichotomous Response Model (DRM) Probabilities

This function computes the probability of a correct response for
multiple items given a set of theta values using the 1PL, 2PL, or 3PL
item response models.

## Usage

``` r
drm(theta, a, b, g = NULL, D = 1)
```

## Arguments

- theta:

  A numeric vector of ability values (latent traits).

- a:

  A numeric vector of item discrimination (slope) parameters.

- b:

  A numeric vector of item difficulty parameters.

- g:

  A numeric vector of item guessing parameters. Not required for 1PL or
  2PL models.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

## Value

A matrix of response probabilities, where rows represent ability values
(theta) and columns represent items.

## Details

If `g` is not specified, the function assumes a guessing parameter of 0
for all items, corresponding to the 1PL or 2PL model. The function
automatically adjusts the model form based on the presence of `g`.

## See also

[`prm()`](https://hwangQ.github.io/irtQ/reference/prm.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1: theta and item parameters for 3PL model
drm(c(-0.1, 0.0, 1.5), a = c(1, 2), b = c(0, 1), g = c(0.2, 0.1), D = 1)
#>           [,1]      [,2]
#> [1,] 0.5800167 0.1897754
#> [2,] 0.6000000 0.2072826
#> [3,] 0.8540596 0.7579527

## Example 2: single theta value with 2PL item parameters
drm(0.0, a = c(1, 2), b = c(0, 1), D = 1)
#>      [,1]      [,2]
#> [1,]  0.5 0.1192029

## Example 3: multiple theta values with a single item (3PL model)
drm(c(-0.1, 0.0, 1.5), a = 1, b = 1, g = 0.2, D = 1)
#>           [,1]
#> [1,] 0.3997919
#> [2,] 0.4151531
#> [3,] 0.6979675
```
