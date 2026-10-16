# Log-Likelihood of Ability Parameters

This function computes the log-likelihood values for a set of ability
values, given item parameters and response data.

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

  A matrix of examinees' item responses, with rows for examinees and
  columns for items in the same order as the items in `x`. Responses to
  an item with K score categories must be integers from 0 to K - 1.
  Missing responses are coded as `NA` or, in functions that have a
  `missing` argument, as the value given in `missing`.

- theta:

  A numeric vector of ability values at which to evaluate the
  log-likelihood function.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- method:

  A character string specifying the estimation method. Available options
  are:

  - `"ML"`: Maximum likelihood estimation

  - `"MLF"`: Maximum likelihood estimation with fences

  - `"MAP"`: Maximum a posteriori estimation

  Default is `"ML"`.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. Used only when
  `method = "MAP"`. Default is `c(0, 1)`.

- fence.a:

  A numeric value specifying the slope (*a*) parameter of the two fence
  items used in MLF. See **Details** below. Default is 3.0.

- fence.b:

  A numeric vector of length two specifying the difficulty (*b*)
  parameters of the lower and upper fence items. If `NULL`, `c(-5, 5)`
  is used. Used only when `method = "MLF"`. Default is `NULL`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

## Value

A data frame of log-likelihood values (log-posterior values for
`method = "MAP"`).

- Each **row** corresponds to a value of `theta`.

- Each **column** corresponds to an examinee and is named `Resp.1`,
  `Resp.2`, and so on. Examinees without any observed response have `NA`
  in their column.

## Details

This function evaluates the log-likelihood at each value of `theta` for
one or more examinees, based on item parameters (`x`) and item response
data (`data`). For `method = "MAP"`, the log density of the normal prior
is added, so the values are log-posterior values.

If `method = "MLF"`, two fence items with slope `fence.a` and
difficulties `fence.b` are appended, with the lower fence answered
correctly and the upper fence answered incorrectly. The fences give the
likelihood a finite maximum for all-correct and all-incorrect patterns.
See
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
for details.

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

# Compute log-likelihood values (ML)
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
