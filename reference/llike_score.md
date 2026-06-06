# Log-Likelihood of Ability Parameters

This function computes the log-likelihood values for a set of ability
parameters, given item parameters and response data

## Usage

``` r
llike_score(
  x,
  data,
  theta,
  D = 1,
  method = "ML",
  norm.prior = c(0, 1),
  fence.a = 3,
  fence.b = NULL,
  missing = NA
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

- theta:

  A numeric vector of ability values at which to evaluate the
  log-likelihood function.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- method:

  A character string specifying the estimation method. Available options
  include:

  - `"ML"`: Maximum Likelihood Estimation

  - `"MLF"`: Maximum Likelihood Estimation with Fences

  - `"MAP"`: Maximum A Posteriori Estimation Default is `"ML"`.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution (used only when
  `method = "MAP"`). Default is `c(0, 1)`. Ignored for `"ML"` and
  `"MLF"`.

- fence.a:

  A numeric value specifying the item slope parameter (i.e.,
  *a*-parameter) for the two imaginary items used in MLF. See
  **Details** below. Default is 3.0.

- fence.b:

  A numeric vector of length two specifying the lower and upper bounds
  of the item difficulty parameters (i.e., *b*-parameters) for the two
  imaginary items in MLF. If `fence.b = NULL`, the values specified in
  the `range` argument are used instead. Default is NULL.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

## Value

A data frame of log-likelihood values.

- Each **row** corresponds to an ability value (`theta`).

- Each **column** corresponds to an examinee’s response pattern.

## Details

This function evaluates the log-likelihood of a given ability (`theta`)
for one or more examinees, based on item parameters (`x`) and item
response data (`data`).

If `method = "MLF"` is selected, the function appends two virtual
"fence" items to the item pool with fixed parameters. These artificial
items help avoid unstable likelihood functions near the boundaries of
the ability scale.

For example, to compute the log-likelihood curves of two examinees'
responses to the same test items, supply a 2-row matrix to `data` and a
vector of ability values to `theta`.

## Examples

``` r
## Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Read item parameters and convert them to item metadata
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df

# Generate ability values from N(0, 1)
set.seed(10)
score <- rnorm(5, mean = 0, sd = 1)

# Simulate response data
data <- simdat(x = x, theta = score, D = 1)

# Specify ability values for log-likelihood evaluation
theta <- seq(-3, 3, 0.5)

# Compute log-likelihood values (using MLE)
llike_score(x = x, data = data, theta = theta, D = 1, method = "ML")
#>       Resp.1    Resp.2     Resp.3     Resp.4    Resp.5
#> 1  -66.50068 -67.04230  -42.73512  -52.76904 -72.98057
#> 2  -61.55760 -61.62838  -40.17094  -48.41493 -67.88062
#> 3  -55.54537 -55.33109  -37.55970  -43.38720 -61.81710
#> 4  -48.75428 -48.55454  -35.55436  -38.20796 -55.12921
#> 5  -42.16532 -42.30980  -35.37400  -34.02680 -48.71271
#> 6  -37.16077 -37.89791  -38.33277  -32.22230 -43.71779
#> 7  -34.86263 -36.36229  -45.20468  -33.79543 -41.03180
#> 8  -35.78900 -38.18251  -56.10789  -39.14549 -41.11082
#> 9  -39.85249 -43.25350  -70.66780  -48.09411 -44.00121
#> 10 -46.65728 -51.18881  -88.30519  -60.15318 -49.51090
#> 11 -55.69829 -61.49921 -108.39067  -74.72129 -57.28597
#> 12 -66.37010 -73.57568 -130.23480  -91.10918 -66.78270
#> 13 -78.09534 -86.81804 -153.20654 -108.67614 -77.42433
```
