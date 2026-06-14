# Simulated Response Data

This function generates simulated response data for single-format or
mixed-format test forms. For dichotomous item response data, the IRT
1PL, 2PL, and 3PL models are supported. For polytomous item response
data, the graded response model (GRM), the partial credit model (PCM),
and the generalized partial credit model (GPCM) are supported.

## Usage

``` r
simdat(
  x = NULL,
  theta,
  a.drm,
  b.drm,
  g.drm = NULL,
  a.prm,
  d.prm,
  cats,
  pr.model,
  D = 1
)
```

## Arguments

- x:

  A data frame containing item metadata. This metadata is required to
  retrieve essential information for each item (e.g., number of score
  categories, IRT model type, etc.) necessary for calibration. You can
  create an empty item metadata frame using the function
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md).
  See **below** for more details. Default is `NULL`.

- theta:

  A numeric vector of ability (theta) values.

- a.drm:

  A numeric vector of item discrimination (slope) parameters for
  dichotomous IRT models.

- b.drm:

  A numeric vector of item difficulty parameters for dichotomous IRT
  models.

- g.drm:

  A numeric vector of guessing parameters for dichotomous IRT models.

- a.prm:

  A numeric vector of item discrimination (slope) parameters for
  polytomous IRT models.

- d.prm:

  A list of numeric vectors, where each vector contains difficulty
  (threshold) parameters for a polytomous item.

- cats:

  A numeric vector indicating the number of score categories for each
  item.

- pr.model:

  A character vector specifying the polytomous IRT model used to
  simulate responses for each polytomous item. Each element should be
  either "GRM" (graded response model) or "GPCM" (generalized partial
  credit model).

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

## Value

A matrix or vector of simulated item responses. If a matrix is returned,
rows correspond to examinees (theta values) and columns to items.

## Details

There are two ways to generate simulated response data. The first is by
providing a data frame of item metadata using the argument `x`. This
data frame must follow a specific structure: the first column should
contain item IDs, the second column should contain the number of unique
score categories for each item, and the third column should specify the
IRT model to be fitted to each item. Available IRT models are:

- `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for dichotomous item data

- `"GRM"` and `"GPCM"` for polytomous item data

Note that `"DRM"` serves as a general label covering all dichotomous IRT
models (i.e., `"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and
`"GPCM"` represent the graded response model and (generalized) partial
credit model, respectively.

The subsequent columns should contain the item parameters for the
specified models. For dichotomous items, the fourth, fifth, and sixth
columns represent item discrimination (slope), item difficulty, and item
guessing parameters, respectively. When `"1PLM"` or `"2PLM"` is
specified in the third column, `NA`s must be entered in the sixth column
for the guessing parameters.

For polytomous items, the item discrimination (slope) parameter should
appear in the fourth column, and the item difficulty (or threshold)
parameters for category boundaries should occupy the fifth through the
last columns. When the number of unique score categories differs across
items, unused parameter cells should be filled with `NA`s.

In the irtQ package, the threshold parameters for GPCM items are
expressed as the item location (or overall difficulty) minus the
threshold values for each score category. Note that when a GPCM item has
*K* unique score categories, *K - 1* threshold parameters are required,
since the threshold for the first category boundary is always fixed at
0. For example, if a GPCM item has five score categories, four threshold
parameters must be provided.

An example of a data frame for a single-format test is shown below:

|       |       |        |       |       |       |        |       |
|-------|-------|--------|-------|-------|-------|--------|-------|
| ITEM1 | 2     | 1PLM   | 1.000 | 1.461 | NA    | ITEM2  | 2     |
| 2PLM  | 1.921 | -1.049 | NA    | ITEM3 | 2     | 3PLM   | 1.736 |
| 1.501 | 0.203 | ITEM4  | 2     | 3PLM  | 0.835 | -1.049 | 0.182 |

An example of a data frame for a mixed-format test is shown below:

|       |     |      |       |        |        |        |        |
|-------|-----|------|-------|--------|--------|--------|--------|
| ITEM1 | 2   | 1PLM | 1.000 | 1.461  | NA     | NA     | NA     |
| ITEM2 | 2   | 2PLM | 1.921 | -1.049 | NA     | NA     | NA     |
| ITEM3 | 2   | 3PLM | 0.926 | 0.394  | 0.099  | NA     | NA     |
| ITEM4 | 2   | DRM  | 1.052 | -0.407 | 0.201  | NA     | NA     |
| ITEM5 | 4   | GRM  | 1.913 | -1.869 | -1.238 | -0.714 | NA     |
| ITEM6 | 5   | GRM  | 1.278 | -0.724 | -0.068 | 0.568  | 1.072  |
| ITEM7 | 4   | GPCM | 1.137 | -0.374 | 0.215  | 0.848  | NA     |
| ITEM8 | 5   | GPCM | 1.233 | -2.078 | -1.347 | -0.705 | -0.116 |

See the *IRT Models* section in the
[irtQ-package](https://hwangQ.github.io/irtQ/reference/irtQ-package.md)
documentation for more details about the IRT models used in the irtQ
package. A convenient way to create a data frame for the argument `x` is
by using the function
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md).

The second approach is to simulate response data by directly specifying
item parameters, instead of providing a metadata data frame via the `x`
argument (see examples below). In this case, the following arguments
must also be specified: `theta`, `cats`, `pr.model`, and `D`.

The `g.drm` argument is only required when simulating dichotomous item
responses under the 3PL model. It can be omitted entirely if all
dichotomous items follow the 1PL or 2PL model. However, if the test
includes a mixture of 1PL, 2PL, and 3PL items, the `g.drm` vector must
be specified for all items, using `NA` for non-3PL items. For example,
if a test consists of four dichotomous items where the first two follow
the 3PL model and the third and fourth follow the 1PL and 2PL models
respectively, then `g.drm = c(0.2, 0.1, NA, NA)` should be used.

For dichotomous items, each element in `cats` should be set to 2. For
polytomous items, the number of unique score categories should be
specified in `cats`. When simulating data for a mixed-format test, it is
important to specify `cats` in the correct item order. For example,
suppose responses are simulated for 10 examinees across 5 items,
including 3 dichotomous items and 2 polytomous items (each with 3
categories), where the second and fourth items are polytomous. In this
case, `cats = c(2, 3, 2, 3, 2)` should be used.

Furthermore, if the two polytomous items are modeled using the graded
response model and the generalized partial credit model, respectively,
then `pr.model = c("GRM", "GPCM")`.

## See also

[`drm()`](https://hwangQ.github.io/irtQ/reference/drm.md),
[`prm()`](https://hwangQ.github.io/irtQ/reference/prm.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1:
## Simulate response data for a mixed-format test.
## The first two polytomous items use the generalized partial credit model (GPCM),
## and the last polytomous item uses the graded response model (GRM).
# Generate theta values for 100 examinees
theta <- rnorm(100)

# Set item parameters for three dichotomous items under the 3PL model
a.drm <- c(1, 1.2, 1.3)
b.drm <- c(-1, 0, 1)
g.drm <- rep(0.2, 3)

# Set item parameters for three polytomous items
# These items have 4, 4, and 5 response categories, respectively
a.prm <- c(1.3, 1.2, 1.7)
d.prm <- list(c(-1.2, -0.3, 0.4), c(-0.2, 0.5, 1.6), c(-1.7, 0.2, 1.1, 2.0))

# Specify the number of score categories for all items
# This vector also determines the location of polytomous items
cats <- c(2, 2, 4, 4, 5, 2)

# Specify the IRT models for the polytomous items
pr.model <- c("GPCM", "GPCM", "GRM")

# Simulate the response data
simdat(
  theta = theta, a.drm = a.drm, b.drm = b.drm, g.drm = NULL,
  a.prm = a.prm, d.prm = d.prm, cats = cats, pr.model = pr.model, D = 1
)
#>        [,1] [,2] [,3] [,4] [,5] [,6]
#>   [1,]    1    1    3    3    3    1
#>   [2,]    0    0    3    1    3    0
#>   [3,]    1    0    0    0    1    0
#>   [4,]    1    1    2    1    1    0
#>   [5,]    1    0    2    1    1    0
#>   [6,]    1    1    2    3    1    1
#>   [7,]    1    1    2    1    1    0
#>   [8,]    0    0    2    0    2    0
#>   [9,]    0    0    0    1    1    0
#>  [10,]    1    0    3    1    2    0
#>  [11,]    1    0    0    1    1    0
#>  [12,]    1    0    0    1    1    0
#>  [13,]    1    0    1    0    2    0
#>  [14,]    0    1    0    0    1    0
#>  [15,]    1    0    2    2    1    1
#>  [16,]    1    0    2    2    2    0
#>  [17,]    0    1    2    0    1    1
#>  [18,]    1    1    2    1    1    1
#>  [19,]    1    1    0    0    1    0
#>  [20,]    0    1    1    1    2    0
#>  [21,]    1    1    3    3    2    1
#>  [22,]    1    1    2    3    2    0
#>  [23,]    1    1    3    0    3    0
#>  [24,]    1    1    3    1    1    0
#>  [25,]    1    1    3    2    1    0
#>  [26,]    1    1    1    0    1    1
#>  [27,]    1    0    0    0    1    0
#>  [28,]    1    0    2    2    1    0
#>  [29,]    0    0    3    2    1    1
#>  [30,]    1    0    0    0    1    0
#>  [31,]    1    1    2    2    1    0
#>  [32,]    0    0    2    0    1    0
#>  [33,]    1    1    3    2    2    0
#>  [34,]    1    1    2    2    2    0
#>  [35,]    0    1    1    0    1    0
#>  [36,]    1    0    3    1    2    0
#>  [37,]    1    0    2    0    4    1
#>  [38,]    0    1    2    0    1    0
#>  [39,]    1    0    2    1    1    0
#>  [40,]    1    0    3    2    4    0
#>  [41,]    1    1    2    1    1    1
#>  [42,]    0    1    0    1    0    0
#>  [43,]    1    1    2    0    2    0
#>  [44,]    0    0    2    1    1    0
#>  [45,]    1    1    2    1    1    0
#>  [46,]    0    1    2    0    2    0
#>  [47,]    1    0    1    2    2    1
#>  [48,]    0    1    3    3    4    1
#>  [49,]    0    1    2    1    1    0
#>  [50,]    0    1    1    0    2    0
#>  [51,]    1    1    3    2    3    0
#>  [52,]    1    0    0    1    0    0
#>  [53,]    0    1    2    0    2    0
#>  [54,]    1    0    2    0    1    0
#>  [55,]    1    1    3    3    4    1
#>  [56,]    1    1    2    1    1    1
#>  [57,]    0    0    0    0    0    0
#>  [58,]    1    0    1    0    1    0
#>  [59,]    1    1    1    1    1    0
#>  [60,]    0    0    0    1    1    0
#>  [61,]    1    1    1    1    0    0
#>  [62,]    0    1    0    0    1    0
#>  [63,]    1    0    0    0    1    0
#>  [64,]    0    1    2    2    3    0
#>  [65,]    1    1    3    2    4    1
#>  [66,]    1    1    3    2    4    1
#>  [67,]    0    0    0    1    1    0
#>  [68,]    1    0    3    3    2    0
#>  [69,]    0    1    1    0    1    0
#>  [70,]    0    0    0    0    1    0
#>  [71,]    0    1    3    2    3    1
#>  [72,]    0    0    1    0    2    0
#>  [73,]    0    1    2    0    3    0
#>  [74,]    0    1    3    2    3    1
#>  [75,]    1    0    1    0    1    0
#>  [76,]    0    1    3    1    0    0
#>  [77,]    1    1    2    1    1    0
#>  [78,]    0    1    3    2    2    0
#>  [79,]    1    1    1    1    1    1
#>  [80,]    1    0    2    2    4    0
#>  [81,]    1    1    3    1    1    0
#>  [82,]    1    0    3    0    0    0
#>  [83,]    0    1    2    2    1    1
#>  [84,]    0    0    0    0    0    0
#>  [85,]    1    1    2    1    2    1
#>  [86,]    0    0    1    0    0    0
#>  [87,]    1    0    2    1    2    1
#>  [88,]    0    1    1    1    3    0
#>  [89,]    0    1    2    2    3    1
#>  [90,]    1    1    3    2    3    0
#>  [91,]    1    1    3    1    3    1
#>  [92,]    0    1    0    0    1    0
#>  [93,]    1    1    2    1    1    0
#>  [94,]    0    0    1    0    1    0
#>  [95,]    0    0    3    1    1    0
#>  [96,]    1    0    1    0    1    0
#>  [97,]    0    0    3    1    1    1
#>  [98,]    1    0    2    1    1    0
#>  [99,]    1    1    2    3    1    1
#> [100,]    1    0    2    2    1    0

## Example 2:
## Simulate response data for a single-format test using the 2PL model
# Specify score categories (2 for each dichotomous item)
cats <- rep(2, 3)

# Simulate the response data
simdat(theta = theta, a.drm = a.drm, b.drm = b.drm, cats = cats, D = 1)
#>        [,1] [,2] [,3]
#>   [1,]    1    1    0
#>   [2,]    0    0    1
#>   [3,]    0    0    0
#>   [4,]    0    1    0
#>   [5,]    1    0    0
#>   [6,]    1    1    1
#>   [7,]    1    0    0
#>   [8,]    1    1    0
#>   [9,]    0    0    0
#>  [10,]    1    0    0
#>  [11,]    0    0    0
#>  [12,]    1    1    0
#>  [13,]    1    1    0
#>  [14,]    1    0    0
#>  [15,]    0    0    0
#>  [16,]    1    1    0
#>  [17,]    1    1    0
#>  [18,]    1    1    0
#>  [19,]    1    1    0
#>  [20,]    1    1    0
#>  [21,]    1    1    1
#>  [22,]    1    1    1
#>  [23,]    1    1    0
#>  [24,]    1    1    0
#>  [25,]    1    0    1
#>  [26,]    1    0    0
#>  [27,]    0    0    0
#>  [28,]    0    0    0
#>  [29,]    1    1    1
#>  [30,]    0    1    0
#>  [31,]    1    1    0
#>  [32,]    1    0    0
#>  [33,]    1    1    1
#>  [34,]    1    1    0
#>  [35,]    1    1    0
#>  [36,]    1    1    1
#>  [37,]    1    1    0
#>  [38,]    1    1    0
#>  [39,]    1    1    0
#>  [40,]    1    1    0
#>  [41,]    1    1    0
#>  [42,]    0    0    0
#>  [43,]    0    0    0
#>  [44,]    1    1    0
#>  [45,]    1    1    0
#>  [46,]    1    0    1
#>  [47,]    0    0    0
#>  [48,]    1    1    1
#>  [49,]    1    1    0
#>  [50,]    0    1    0
#>  [51,]    1    1    1
#>  [52,]    1    1    0
#>  [53,]    0    0    0
#>  [54,]    0    0    0
#>  [55,]    1    0    0
#>  [56,]    0    0    0
#>  [57,]    0    1    0
#>  [58,]    0    0    0
#>  [59,]    1    0    0
#>  [60,]    0    0    0
#>  [61,]    0    1    0
#>  [62,]    0    1    0
#>  [63,]    0    0    0
#>  [64,]    0    1    0
#>  [65,]    1    1    1
#>  [66,]    0    1    0
#>  [67,]    0    0    0
#>  [68,]    1    0    1
#>  [69,]    1    0    1
#>  [70,]    0    0    0
#>  [71,]    1    1    1
#>  [72,]    0    0    0
#>  [73,]    1    0    0
#>  [74,]    1    1    0
#>  [75,]    1    0    0
#>  [76,]    1    1    0
#>  [77,]    1    1    0
#>  [78,]    0    1    0
#>  [79,]    1    0    0
#>  [80,]    1    1    0
#>  [81,]    1    0    1
#>  [82,]    1    0    0
#>  [83,]    1    0    1
#>  [84,]    0    0    0
#>  [85,]    1    0    1
#>  [86,]    0    0    0
#>  [87,]    1    1    1
#>  [88,]    1    1    1
#>  [89,]    1    1    0
#>  [90,]    1    1    1
#>  [91,]    1    1    0
#>  [92,]    1    0    0
#>  [93,]    1    0    0
#>  [94,]    0    0    0
#>  [95,]    0    0    0
#>  [96,]    0    1    0
#>  [97,]    0    0    0
#>  [98,]    1    1    0
#>  [99,]    1    1    0
#> [100,]    1    1    0

## Example 3:
## Simulate response data using a "-prm.txt" file exported from flexMIRT
# Load the flexMIRT parameter file
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Convert the flexMIRT parameters to item metadata
test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df

# Simulate the response data using the item metadata
simdat(x = test_flex, theta = theta, D = 1)
#>        [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>   [1,]    1    1    1    1    1    1    1    1    1     1     1     0     0
#>   [2,]    0    1    0    1    1    1    1    1    0     1     1     0     0
#>   [3,]    1    1    1    1    1    0    0    0    1     0     0     0     0
#>   [4,]    0    1    0    0    0    0    0    1    0     1     1     0     0
#>   [5,]    1    0    0    1    1    1    1    0    0     1     1     1     1
#>   [6,]    0    1    1    0    1    1    1    1    1     1     1     0     1
#>   [7,]    0    1    0    0    0    0    1    0    0     0     1     0     0
#>   [8,]    1    1    0    1    1    1    0    0    0     0     1     0     1
#>   [9,]    0    0    1    0    0    0    0    0    1     0     1     0     0
#>  [10,]    0    1    0    1    1    0    1    1    1     0     1     1     0
#>  [11,]    0    0    0    1    1    0    0    0    0     0     0     0     0
#>  [12,]    1    1    0    1    1    0    1    0    0     1     1     0     1
#>  [13,]    0    1    1    1    1    1    1    0    1     1     1     1     1
#>  [14,]    0    1    0    1    0    0    0    1    1     1     0     1     0
#>  [15,]    1    0    1    1    0    0    0    0    0     0     1     0     0
#>  [16,]    0    1    0    1    1    0    0    1    0     0     1     1     0
#>  [17,]    1    1    0    1    1    1    0    0    1     0     1     0     0
#>  [18,]    1    1    1    1    1    1    0    0    1     1     1     0     1
#>  [19,]    0    0    0    0    1    0    0    0    1     0     1     0     1
#>  [20,]    1    1    0    1    1    1    1    1    1     1     1     0     1
#>  [21,]    0    1    1    0    1    0    0    1    0     1     1     1     0
#>  [22,]    0    1    1    1    1    1    1    0    0     1     1     0     1
#>  [23,]    1    1    1    0    1    0    0    1    0     1     0     0     0
#>  [24,]    1    1    0    1    1    0    0    1    0     0     1     0     0
#>  [25,]    0    1    0    1    1    1    0    1    0     1     1     1     1
#>  [26,]    0    1    0    0    1    1    0    0    1     1     1     0     0
#>  [27,]    0    1    1    0    0    0    0    0    0     0     0     0     0
#>  [28,]    0    1    0    1    0    0    0    1    1     1     0     0     0
#>  [29,]    0    1    1    1    0    0    0    0    1     1     1     0     0
#>  [30,]    0    0    0    0    1    0    0    0    0     1     1     0     0
#>  [31,]    1    1    1    1    1    1    1    1    1     1     0     1     0
#>  [32,]    1    1    0    0    0    1    0    0    0     0     1     0     0
#>  [33,]    1    1    1    1    1    1    0    0    1     1     1     1     0
#>  [34,]    0    1    1    1    1    1    1    1    0     1     1     0     0
#>  [35,]    0    1    1    1    1    1    0    0    0     0     1     0     0
#>  [36,]    1    1    0    1    1    1    1    1    0     0     0     1     0
#>  [37,]    0    1    0    0    0    1    0    1    0     1     0     1     0
#>  [38,]    1    1    0    1    1    1    0    0    1     1     1     1     1
#>  [39,]    1    1    1    1    0    0    0    1    1     0     0     0     0
#>  [40,]    1    1    1    1    0    1    1    0    0     1     1     0     0
#>  [41,]    1    1    1    1    1    0    1    1    1     1     1     0     0
#>  [42,]    1    1    0    0    0    0    0    0    1     0     0     0     0
#>  [43,]    1    1    0    1    0    0    1    1    1     0     1     1     0
#>  [44,]    0    1    1    1    0    0    1    1    1     0     0     0     1
#>  [45,]    0    1    1    1    0    0    0    1    0     1     1     0     0
#>  [46,]    0    1    0    0    1    0    0    0    0     1     0     0     0
#>  [47,]    0    1    0    1    1    0    0    0    1     0     1     0     1
#>  [48,]    1    1    0    0    1    1    0    1    1     1     0     1     1
#>  [49,]    1    1    1    1    1    0    0    1    0     1     1     1     0
#>  [50,]    0    1    0    0    1    0    0    0    1     0     1     0     0
#>  [51,]    0    1    1    1    1    1    1    0    1     1     1     1     1
#>  [52,]    0    0    0    1    0    0    0    0    0     1     0     0     1
#>  [53,]    0    1    0    0    0    0    0    0    1     0     0     0     1
#>  [54,]    0    1    0    0    0    0    1    0    0     0     0     0     0
#>  [55,]    1    1    1    1    1    0    1    1    1     0     1     1     1
#>  [56,]    1    1    1    0    0    1    1    0    0     0     1     0     0
#>  [57,]    0    0    0    0    0    0    1    0    0     0     1     0     0
#>  [58,]    0    1    0    0    0    1    0    0    0     1     1     0     0
#>  [59,]    0    1    0    1    0    0    1    0    1     1     0     0     0
#>  [60,]    1    0    0    1    0    0    0    1    1     0     0     0     1
#>  [61,]    0    1    0    1    0    0    0    0    1     1     0     0     1
#>  [62,]    1    1    0    1    0    0    0    0    0     0     0     0     0
#>  [63,]    1    0    1    0    0    0    0    0    1     0     0     0     1
#>  [64,]    1    1    0    0    1    0    1    1    1     1     1     1     0
#>  [65,]    0    1    1    0    1    1    1    1    1     1     1     1     1
#>  [66,]    0    1    1    1    1    1    0    1    1     1     1     0     0
#>  [67,]    0    1    1    0    1    0    0    0    0     0     0     0     0
#>  [68,]    0    1    1    1    1    0    0    0    0     1     1     0     0
#>  [69,]    1    1    0    0    1    0    0    1    0     1     1     1     0
#>  [70,]    0    0    0    1    1    0    0    0    0     0     0     0     0
#>  [71,]    1    1    1    0    1    1    0    1    1     1     0     0     0
#>  [72,]    0    1    0    1    1    0    0    0    1     1     1     1     0
#>  [73,]    0    1    0    1    1    0    1    0    0     1     1     0     1
#>  [74,]    1    1    1    1    1    1    1    1    0     0     1     0     1
#>  [75,]    0    0    1    0    0    0    1    0    0     1     0     1     0
#>  [76,]    0    1    1    1    1    0    0    0    0     1     1     1     1
#>  [77,]    0    1    0    1    1    1    1    1    1     1     1     0     0
#>  [78,]    0    1    1    1    1    0    1    1    0     1     1     0     0
#>  [79,]    0    1    0    0    1    0    0    1    1     1     1     0     0
#>  [80,]    0    1    1    1    1    1    0    1    1     1     1     0     1
#>  [81,]    1    1    0    1    1    1    0    0    0     1     1     0     0
#>  [82,]    0    1    1    1    0    1    0    1    0     0     1     1     0
#>  [83,]    0    1    0    0    0    0    0    1    0     1     0     0     0
#>  [84,]    0    1    0    1    1    0    0    0    1     0     0     0     1
#>  [85,]    0    1    1    1    0    0    1    1    1     1     0     0     0
#>  [86,]    0    0    0    0    0    0    0    1    1     1     0     0     0
#>  [87,]    0    1    1    1    0    1    1    1    1     0     1     0     0
#>  [88,]    1    1    0    0    1    1    1    1    1     1     1     0     0
#>  [89,]    1    0    0    1    0    0    0    0    0     1     0     1     0
#>  [90,]    1    1    1    1    1    1    0    0    1     1     1     1     1
#>  [91,]    0    1    0    1    1    0    1    0    0     0     1     1     0
#>  [92,]    0    1    0    0    0    0    1    0    1     0     0     0     0
#>  [93,]    0    1    1    0    1    0    0    0    1     0     1     0     0
#>  [94,]    0    1    1    1    0    1    1    0    0     1     1     0     0
#>  [95,]    0    1    0    1    0    0    1    0    0     0     1     0     0
#>  [96,]    1    1    0    1    1    0    0    1    0     1     0     0     0
#>  [97,]    0    1    0    1    0    0    0    0    0     0     1     0     0
#>  [98,]    1    0    1    1    1    1    1    1    1     1     0     1     0
#>  [99,]    1    1    1    1    1    1    1    0    1     1     1     0     0
#> [100,]    0    1    0    0    0    1    1    0    0     1     1     0     0
#>        [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
#>   [1,]     1     1     1     1     1     1     1     1     1     1     0     1
#>   [2,]     0     0     1     0     0     1     1     1     1     1     0     0
#>   [3,]     0     1     0     0     1     0     1     0     1     0     0     1
#>   [4,]     0     0     0     0     0     1     1     1     1     1     0     1
#>   [5,]     1     0     0     0     0     0     1     0     1     1     1     1
#>   [6,]     1     1     1     1     0     1     1     1     1     1     1     1
#>   [7,]     0     0     1     1     0     1     1     1     1     1     1     1
#>   [8,]     1     1     1     1     0     1     1     1     1     1     0     1
#>   [9,]     1     1     0     1     0     1     1     1     1     1     1     1
#>  [10,]     0     0     1     0     0     1     1     1     1     1     0     1
#>  [11,]     0     0     0     1     0     0     1     0     0     1     0     0
#>  [12,]     0     0     0     1     1     1     1     1     1     0     1     1
#>  [13,]     0     1     1     1     0     1     1     1     1     1     0     1
#>  [14,]     0     1     0     0     0     0     1     1     1     0     1     1
#>  [15,]     1     0     1     1     1     1     1     0     1     1     0     1
#>  [16,]     1     1     1     1     0     1     1     1     1     0     0     1
#>  [17,]     1     1     1     0     1     1     1     1     1     0     1     1
#>  [18,]     1     1     0     1     1     1     1     1     0     1     1     1
#>  [19,]     0     1     1     0     0     0     1     1     1     0     1     1
#>  [20,]     1     1     1     1     0     1     1     1     1     1     1     1
#>  [21,]     0     1     1     1     1     1     1     1     1     1     0     1
#>  [22,]     1     1     1     0     0     1     1     1     1     1     0     1
#>  [23,]     0     0     1     0     0     1     1     1     0     1     0     1
#>  [24,]     1     1     1     1     1     1     1     1     1     1     0     1
#>  [25,]     1     1     0     1     0     1     1     1     1     1     1     1
#>  [26,]     1     1     1     1     0     1     1     1     0     0     0     1
#>  [27,]     0     0     0     0     0     1     1     1     0     0     0     0
#>  [28,]     0     0     0     0     0     0     0     0     0     0     0     1
#>  [29,]     0     1     0     1     0     1     1     1     1     0     1     1
#>  [30,]     0     0     1     1     1     1     0     1     0     0     0     1
#>  [31,]     1     1     1     1     0     1     1     1     1     1     1     1
#>  [32,]     0     1     1     0     0     1     1     1     1     0     1     1
#>  [33,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [34,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [35,]     1     1     0     1     0     1     1     0     1     1     0     1
#>  [36,]     1     1     1     0     1     1     1     1     1     1     1     1
#>  [37,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [38,]     0     1     0     1     1     1     1     1     1     1     1     1
#>  [39,]     1     1     1     1     0     1     1     1     0     1     0     0
#>  [40,]     1     0     0     1     1     1     1     0     0     1     1     1
#>  [41,]     1     1     1     0     0     1     1     1     1     1     0     1
#>  [42,]     0     0     0     0     0     0     0     1     1     0     0     0
#>  [43,]     0     0     0     1     0     1     1     1     1     0     0     1
#>  [44,]     1     0     1     0     0     1     1     1     1     1     0     0
#>  [45,]     1     1     0     0     1     1     1     1     1     0     0     1
#>  [46,]     0     1     0     1     1     1     1     1     0     0     0     1
#>  [47,]     1     0     1     1     0     1     1     1     1     1     0     1
#>  [48,]     1     1     1     1     1     1     1     1     0     1     1     1
#>  [49,]     0     1     1     0     1     1     1     1     0     1     0     1
#>  [50,]     1     0     0     0     0     1     1     1     1     1     0     0
#>  [51,]     1     1     1     0     1     1     1     1     1     1     1     1
#>  [52,]     1     0     0     0     0     1     0     1     1     1     0     1
#>  [53,]     0     0     0     0     0     1     1     0     0     1     1     1
#>  [54,]     1     1     0     0     1     0     1     0     0     0     0     0
#>  [55,]     1     0     0     0     0     1     1     1     1     0     1     1
#>  [56,]     1     0     1     1     0     1     1     1     0     1     1     0
#>  [57,]     0     0     0     1     0     0     0     0     1     0     0     0
#>  [58,]     0     0     0     1     0     1     1     1     1     1     0     1
#>  [59,]     1     1     0     1     1     1     1     1     1     1     0     0
#>  [60,]     1     0     0     0     0     1     1     0     1     0     1     1
#>  [61,]     1     1     0     0     0     0     0     0     0     1     0     1
#>  [62,]     1     0     0     1     1     1     1     1     0     0     0     1
#>  [63,]     1     1     0     0     0     1     0     0     0     0     0     1
#>  [64,]     0     0     1     1     1     1     1     1     1     1     0     1
#>  [65,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [66,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [67,]     0     1     0     0     0     1     1     0     0     0     0     1
#>  [68,]     1     1     1     1     0     1     1     1     1     0     1     1
#>  [69,]     1     0     1     1     0     1     1     1     1     1     0     1
#>  [70,]     1     0     1     0     0     1     1     1     1     1     0     1
#>  [71,]     1     1     1     1     0     1     1     1     1     1     1     1
#>  [72,]     1     0     0     0     0     1     1     1     0     0     0     1
#>  [73,]     0     0     0     0     0     1     1     1     0     1     1     0
#>  [74,]     0     1     1     1     0     1     1     1     0     1     1     1
#>  [75,]     1     0     1     1     1     1     1     1     0     0     0     1
#>  [76,]     1     1     1     0     0     1     1     1     1     0     0     1
#>  [77,]     1     0     0     1     1     1     1     1     1     1     0     1
#>  [78,]     1     1     1     0     0     1     1     1     1     1     0     1
#>  [79,]     0     1     1     0     0     1     1     1     1     1     0     1
#>  [80,]     1     1     1     1     1     1     1     1     1     0     1     1
#>  [81,]     1     0     1     1     0     1     1     1     0     1     1     1
#>  [82,]     1     1     0     1     0     1     1     1     1     0     0     1
#>  [83,]     1     1     0     0     0     1     1     1     1     1     1     1
#>  [84,]     1     0     0     1     1     0     0     1     0     0     0     1
#>  [85,]     0     1     1     1     1     1     1     1     0     1     0     1
#>  [86,]     0     0     0     1     0     1     0     0     0     0     0     0
#>  [87,]     1     0     0     1     0     1     1     1     0     1     0     1
#>  [88,]     1     1     1     1     1     1     1     0     1     1     1     1
#>  [89,]     0     1     0     1     0     1     1     0     0     1     0     1
#>  [90,]     1     0     1     0     1     1     1     1     1     1     1     1
#>  [91,]     1     1     1     1     0     1     1     1     0     0     1     1
#>  [92,]     0     0     0     1     1     1     1     1     1     1     0     1
#>  [93,]     1     0     0     0     0     0     1     1     1     1     1     0
#>  [94,]     1     1     1     0     0     1     1     1     1     1     0     1
#>  [95,]     1     1     0     0     1     1     0     1     1     1     1     1
#>  [96,]     0     0     0     0     1     1     1     0     0     0     1     0
#>  [97,]     0     0     0     0     0     1     1     1     1     0     0     1
#>  [98,]     1     1     1     0     0     1     1     1     0     0     1     1
#>  [99,]     0     1     1     1     1     1     1     1     1     1     1     1
#> [100,]     1     1     0     0     1     1     1     1     1     1     0     1
#>        [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37]
#>   [1,]     1     1     1     1     0     0     1     1     1     1     1     1
#>   [2,]     1     1     1     0     1     0     1     1     0     0     0     1
#>   [3,]     1     0     0     1     1     0     0     0     0     0     1     1
#>   [4,]     1     0     0     1     1     0     1     1     1     1     0     1
#>   [5,]     1     0     1     0     1     1     0     1     0     1     0     1
#>   [6,]     1     0     1     1     1     1     1     1     1     1     1     1
#>   [7,]     1     1     0     1     0     0     1     1     1     1     0     1
#>   [8,]     1     1     1     1     1     1     1     1     0     1     1     1
#>   [9,]     0     0     0     1     1     0     1     1     0     0     0     0
#>  [10,]     1     1     1     1     1     0     1     1     1     1     1     1
#>  [11,]     1     0     0     0     1     1     1     1     1     0     0     0
#>  [12,]     1     0     0     1     1     0     1     0     0     1     0     0
#>  [13,]     1     1     1     1     0     1     1     1     1     1     0     0
#>  [14,]     1     0     0     1     0     0     1     1     1     1     0     1
#>  [15,]     0     0     0     1     1     0     1     1     1     1     0     0
#>  [16,]     1     0     1     1     0     0     0     1     1     1     1     1
#>  [17,]     1     1     1     0     0     0     1     1     1     1     1     1
#>  [18,]     1     1     0     1     1     0     1     1     0     0     1     1
#>  [19,]     0     1     0     1     0     1     1     1     0     0     0     0
#>  [20,]     1     1     0     1     0     0     1     1     0     1     0     1
#>  [21,]     1     0     1     1     1     0     1     1     0     1     0     1
#>  [22,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [23,]     1     1     1     1     0     0     0     1     0     0     0     1
#>  [24,]     1     0     1     1     0     1     1     0     0     1     1     1
#>  [25,]     1     0     1     1     1     1     1     1     1     0     1     1
#>  [26,]     1     1     0     1     1     0     0     1     1     1     0     0
#>  [27,]     1     0     1     0     0     0     1     0     1     0     0     0
#>  [28,]     1     0     0     1     0     0     1     1     0     0     1     1
#>  [29,]     1     0     1     1     0     1     1     1     1     1     0     0
#>  [30,]     1     0     0     1     0     0     1     1     0     1     1     1
#>  [31,]     1     1     1     1     0     1     1     1     1     1     0     1
#>  [32,]     1     1     0     0     0     0     1     1     1     0     0     0
#>  [33,]     1     1     1     1     1     0     1     1     0     1     1     1
#>  [34,]     1     1     1     1     0     0     1     1     0     1     0     1
#>  [35,]     1     0     1     1     1     1     1     1     1     1     0     0
#>  [36,]     1     0     1     1     1     0     1     1     0     1     0     1
#>  [37,]     1     0     1     1     1     0     1     1     1     1     0     1
#>  [38,]     1     1     1     1     0     0     1     1     1     1     1     1
#>  [39,]     1     0     1     1     0     0     1     1     0     1     0     0
#>  [40,]     1     1     1     1     1     0     1     1     0     1     0     1
#>  [41,]     1     1     1     1     1     0     1     1     0     1     1     1
#>  [42,]     1     0     1     1     0     0     0     0     1     0     0     1
#>  [43,]     0     0     0     1     0     0     1     1     1     1     1     1
#>  [44,]     1     1     1     0     0     0     1     1     1     0     1     0
#>  [45,]     1     0     1     1     0     1     1     1     1     0     0     1
#>  [46,]     1     0     0     1     1     0     1     1     1     0     1     1
#>  [47,]     1     1     1     1     0     1     1     1     1     1     0     0
#>  [48,]     1     0     1     1     0     0     1     1     1     1     1     1
#>  [49,]     1     1     1     1     1     0     1     1     1     1     0     1
#>  [50,]     1     1     0     1     0     0     0     1     0     0     0     0
#>  [51,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [52,]     1     1     1     1     1     0     1     0     1     1     0     1
#>  [53,]     0     1     0     0     0     0     0     1     0     0     0     1
#>  [54,]     1     0     1     1     1     0     0     0     1     0     1     0
#>  [55,]     1     1     1     1     0     1     1     1     0     1     0     1
#>  [56,]     1     0     1     1     0     0     1     1     0     1     0     1
#>  [57,]     1     1     0     0     0     0     1     1     0     1     0     0
#>  [58,]     0     0     0     1     0     0     0     1     0     0     0     0
#>  [59,]     1     1     0     1     0     0     1     1     1     1     1     0
#>  [60,]     1     0     1     0     0     0     0     1     0     0     0     0
#>  [61,]     1     0     0     1     1     0     0     0     1     1     0     0
#>  [62,]     1     0     0     0     0     0     0     0     0     1     1     0
#>  [63,]     1     0     1     1     0     0     0     0     1     1     0     1
#>  [64,]     1     1     1     1     1     1     1     0     1     1     0     1
#>  [65,]     1     1     1     1     1     1     1     1     1     1     0     1
#>  [66,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [67,]     0     0     0     1     1     0     0     1     0     0     0     0
#>  [68,]     1     1     1     1     1     0     1     1     1     1     0     1
#>  [69,]     1     1     1     1     1     0     1     1     0     1     1     0
#>  [70,]     1     0     0     0     0     0     1     0     0     0     0     1
#>  [71,]     1     1     1     1     1     1     1     0     0     1     1     1
#>  [72,]     1     1     1     1     0     1     1     1     1     1     0     0
#>  [73,]     0     0     0     1     1     0     1     1     0     1     0     1
#>  [74,]     1     1     1     1     1     1     1     0     1     0     0     1
#>  [75,]     1     1     0     1     0     0     1     0     0     0     0     0
#>  [76,]     1     0     0     0     0     0     0     0     0     1     0     0
#>  [77,]     1     1     0     1     1     0     1     1     1     0     1     0
#>  [78,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [79,]     1     1     0     1     0     0     1     0     0     1     0     1
#>  [80,]     1     1     1     1     0     1     1     1     1     0     1     1
#>  [81,]     1     1     1     0     1     0     1     1     1     1     1     1
#>  [82,]     0     0     0     0     0     1     1     1     1     1     1     0
#>  [83,]     0     1     0     1     1     1     1     0     1     1     0     1
#>  [84,]     1     0     0     0     1     0     0     1     0     1     0     1
#>  [85,]     1     0     0     0     1     0     1     1     1     1     0     1
#>  [86,]     0     0     1     0     1     0     1     0     0     0     0     1
#>  [87,]     1     1     1     1     0     0     1     1     1     1     0     1
#>  [88,]     1     1     1     1     1     1     1     1     0     1     0     1
#>  [89,]     1     1     1     1     1     1     1     1     1     1     1     1
#>  [90,]     1     1     1     1     0     1     1     1     1     1     1     1
#>  [91,]     1     1     1     0     1     1     1     1     1     1     0     0
#>  [92,]     1     1     0     1     1     1     1     1     1     1     0     1
#>  [93,]     1     1     0     1     1     0     1     1     0     1     1     1
#>  [94,]     1     0     0     1     0     1     1     1     1     0     1     1
#>  [95,]     1     1     0     1     0     1     0     1     0     1     0     1
#>  [96,]     1     0     0     1     0     1     1     1     0     0     1     0
#>  [97,]     1     0     0     1     1     1     0     1     0     0     0     1
#>  [98,]     1     1     1     1     0     0     1     1     1     1     0     1
#>  [99,]     1     1     1     1     1     0     1     1     0     1     0     1
#> [100,]     1     0     1     1     0     1     1     1     1     1     1     0
#>        [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49]
#>   [1,]     0     4     4     1     1     0     1     0     1     1     1     1
#>   [2,]     1     4     0     1     1     0     1     1     1     0     0     0
#>   [3,]     1     3     0     0     0     0     1     0     0     0     0     0
#>   [4,]     0     1     2     1     1     1     1     0     1     1     0     0
#>   [5,]     1     3     3     1     1     0     1     1     0     0     0     1
#>   [6,]     1     4     2     1     1     1     1     1     1     1     1     0
#>   [7,]     1     2     0     0     1     0     0     0     0     1     0     1
#>   [8,]     1     2     4     0     1     0     1     0     0     1     0     1
#>   [9,]     1     1     2     1     1     0     0     0     0     0     1     0
#>  [10,]     1     4     2     0     1     1     1     0     0     1     1     1
#>  [11,]     1     3     0     0     0     0     1     1     0     1     0     0
#>  [12,]     0     3     4     1     1     1     0     1     0     1     0     1
#>  [13,]     0     4     4     1     1     0     0     1     1     1     1     1
#>  [14,]     1     1     4     0     1     1     1     0     1     0     1     0
#>  [15,]     1     4     2     1     1     0     0     0     0     0     1     1
#>  [16,]     1     4     3     1     1     0     0     1     0     1     0     0
#>  [17,]     1     3     3     0     1     1     0     0     0     1     1     0
#>  [18,]     0     4     2     1     1     1     1     0     0     1     1     1
#>  [19,]     1     3     4     0     1     0     1     0     0     0     1     0
#>  [20,]     1     4     4     1     1     0     0     0     0     1     1     1
#>  [21,]     0     3     2     1     1     1     1     0     1     1     1     1
#>  [22,]     1     4     4     1     1     1     1     0     0     1     1     1
#>  [23,]     0     4     1     0     1     1     1     0     0     0     1     0
#>  [24,]     1     4     4     0     1     1     1     0     0     1     0     0
#>  [25,]     1     4     2     0     1     1     1     1     0     1     1     1
#>  [26,]     1     1     3     1     1     0     1     0     0     0     1     0
#>  [27,]     1     0     0     0     1     0     0     0     0     0     0     0
#>  [28,]     1     1     0     0     1     0     1     0     0     1     1     0
#>  [29,]     1     3     0     1     1     0     1     0     0     0     1     0
#>  [30,]     0     2     1     1     1     0     0     0     0     0     0     0
#>  [31,]     1     4     3     1     1     0     1     0     1     1     1     1
#>  [32,]     1     4     0     0     1     0     1     0     1     1     0     0
#>  [33,]     0     3     4     1     1     1     1     1     0     1     1     1
#>  [34,]     0     4     4     1     1     1     1     0     0     1     1     1
#>  [35,]     1     3     0     0     1     0     1     0     1     1     1     0
#>  [36,]     0     3     2     0     1     1     1     0     1     0     0     1
#>  [37,]     1     4     4     0     1     0     1     0     0     1     1     1
#>  [38,]     1     4     4     1     1     1     1     0     1     1     1     1
#>  [39,]     1     4     0     1     1     0     0     0     0     0     1     1
#>  [40,]     1     4     0     0     1     0     1     0     0     1     0     0
#>  [41,]     0     4     1     1     1     0     0     1     0     1     0     0
#>  [42,]     1     0     4     0     0     0     0     0     0     0     1     0
#>  [43,]     0     3     3     0     1     0     0     1     1     0     0     1
#>  [44,]     0     3     0     0     1     0     0     0     0     0     1     0
#>  [45,]     1     3     2     1     1     0     1     0     0     1     0     1
#>  [46,]     0     4     4     0     1     0     1     0     0     1     1     1
#>  [47,]     1     4     2     1     1     0     0     1     1     0     0     0
#>  [48,]     1     1     4     1     1     1     1     0     0     1     1     1
#>  [49,]     0     2     1     0     1     1     1     0     0     0     0     1
#>  [50,]     0     4     3     0     1     0     0     0     0     0     1     1
#>  [51,]     1     4     0     1     1     1     1     0     0     1     1     1
#>  [52,]     1     1     0     1     1     0     0     0     0     0     0     1
#>  [53,]     0     2     0     0     1     1     1     0     0     0     0     1
#>  [54,]     0     0     0     0     1     0     0     0     0     0     0     0
#>  [55,]     0     4     4     1     1     0     1     0     0     0     1     1
#>  [56,]     1     2     0     1     1     0     1     0     1     0     0     1
#>  [57,]     0     0     0     1     1     0     1     0     1     0     0     0
#>  [58,]     0     4     0     0     1     0     0     0     1     1     1     0
#>  [59,]     0     4     3     1     1     1     1     0     0     1     0     0
#>  [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [61,]     0     4     0     0     1     0     0     1     0     0     1     0
#>  [62,]     1     4     3     1     1     0     1     0     0     1     1     1
#>  [63,]     0     0     0     0     1     0     0     0     0     0     0     0
#>  [64,]     1     4     0     1     1     0     1     0     0     1     1     0
#>  [65,]     1     4     2     1     1     1     1     1     1     1     1     1
#>  [66,]     1     4     4     0     1     1     1     0     0     1     1     1
#>  [67,]     1     2     0     1     0     1     0     0     0     1     0     0
#>  [68,]     0     4     2     0     1     0     0     0     1     1     0     1
#>  [69,]     1     4     4     1     1     1     0     0     0     0     1     1
#>  [70,]     0     1     4     0     0     0     1     0     1     0     0     1
#>  [71,]     0     4     4     1     1     1     1     0     0     0     1     0
#>  [72,]     1     2     4     1     1     0     0     0     0     1     1     0
#>  [73,]     1     3     3     0     1     1     1     0     0     1     1     0
#>  [74,]     0     4     1     0     1     1     1     0     0     1     1     1
#>  [75,]     0     0     4     1     1     0     0     1     0     1     0     1
#>  [76,]     0     4     2     0     1     1     1     0     0     1     1     1
#>  [77,]     1     4     1     0     1     1     1     1     0     1     0     1
#>  [78,]     0     4     4     1     1     1     1     0     1     0     0     1
#>  [79,]     0     4     2     0     1     0     0     0     0     0     0     0
#>  [80,]     0     4     4     1     1     0     1     0     0     1     1     1
#>  [81,]     1     4     1     1     1     1     0     0     0     0     1     0
#>  [82,]     0     4     1     0     1     1     0     0     1     0     1     0
#>  [83,]     0     4     0     0     1     0     1     0     0     1     1     1
#>  [84,]     1     1     0     0     1     0     1     0     0     0     0     0
#>  [85,]     1     4     3     0     1     0     1     0     0     1     0     1
#>  [86,]     1     0     0     1     0     0     0     1     1     0     0     0
#>  [87,]     1     4     1     1     1     1     1     0     0     0     1     1
#>  [88,]     1     3     1     0     1     0     1     0     0     1     1     1
#>  [89,]     1     4     2     1     1     1     1     0     0     1     0     0
#>  [90,]     0     4     1     1     1     1     0     0     1     1     1     1
#>  [91,]     1     4     4     1     1     1     1     0     0     1     0     1
#>  [92,]     1     1     0     0     1     0     0     0     1     0     1     0
#>  [93,]     1     3     4     0     1     0     0     1     0     0     0     0
#>  [94,]     0     0     0     0     0     0     0     0     0     1     0     0
#>  [95,]     0     3     0     0     1     0     1     0     0     0     0     1
#>  [96,]     1     0     1     0     0     0     0     0     0     0     0     0
#>  [97,]     0     3     0     0     1     0     0     0     0     1     1     1
#>  [98,]     1     4     4     1     1     1     1     1     1     1     1     1
#>  [99,]     1     4     4     1     1     1     1     0     1     1     0     1
#> [100,]     1     3     4     1     1     1     1     0     1     0     0     0
#>        [,50] [,51] [,52] [,53] [,54] [,55]
#>   [1,]     0     1     1     0     2     4
#>   [2,]     1     1     1     2     4     3
#>   [3,]     0     1     1     3     1     4
#>   [4,]     0     1     1     1     4     0
#>   [5,]     0     1     1     3     1     1
#>   [6,]     0     1     1     4     4     2
#>   [7,]     0     0     0     0     2     3
#>   [8,]     0     1     1     0     4     3
#>   [9,]     0     0     1     1     0     4
#>  [10,]     1     0     1     3     2     3
#>  [11,]     0     0     0     1     0     0
#>  [12,]     1     0     1     0     2     1
#>  [13,]     1     1     1     1     4     2
#>  [14,]     0     0     1     0     3     2
#>  [15,]     0     1     0     0     3     0
#>  [16,]     0     1     1     2     2     4
#>  [17,]     0     0     1     0     2     3
#>  [18,]     1     1     1     2     3     4
#>  [19,]     0     0     1     0     1     0
#>  [20,]     1     1     1     0     4     1
#>  [21,]     0     1     0     1     4     1
#>  [22,]     0     1     1     0     4     0
#>  [23,]     1     0     1     0     3     0
#>  [24,]     0     1     1     3     4     3
#>  [25,]     0     1     1     4     4     3
#>  [26,]     1     1     1     0     1     4
#>  [27,]     1     0     0     3     0     2
#>  [28,]     0     0     1     4     1     0
#>  [29,]     0     0     0     4     3     3
#>  [30,]     0     0     1     0     0     2
#>  [31,]     1     1     1     4     0     4
#>  [32,]     0     0     0     0     0     3
#>  [33,]     0     1     1     0     4     4
#>  [34,]     0     1     0     2     2     4
#>  [35,]     0     1     0     1     4     4
#>  [36,]     0     1     1     0     4     0
#>  [37,]     0     1     0     4     4     4
#>  [38,]     1     1     1     4     4     2
#>  [39,]     0     1     1     0     4     0
#>  [40,]     0     1     1     4     1     4
#>  [41,]     0     1     1     2     2     4
#>  [42,]     0     0     0     0     1     0
#>  [43,]     0     0     1     0     2     4
#>  [44,]     0     1     0     0     4     1
#>  [45,]     0     1     1     0     4     4
#>  [46,]     0     1     1     4     1     0
#>  [47,]     0     1     1     3     4     0
#>  [48,]     0     1     1     4     0     3
#>  [49,]     0     1     1     4     4     0
#>  [50,]     0     1     0     0     2     0
#>  [51,]     0     1     0     0     4     4
#>  [52,]     0     0     1     0     4     2
#>  [53,]     0     1     1     0     4     1
#>  [54,]     0     0     0     0     1     0
#>  [55,]     0     1     1     4     4     3
#>  [56,]     1     1     1     4     4     1
#>  [57,]     0     0     0     2     0     0
#>  [58,]     0     1     1     2     0     4
#>  [59,]     1     1     1     1     1     0
#>  [60,]     0     1     0     0     4     2
#>  [61,]     0     1     0     0     3     1
#>  [62,]     0     1     0     4     3     1
#>  [63,]     0     0     1     0     1     0
#>  [64,]     1     1     0     2     1     4
#>  [65,]     0     1     1     0     4     4
#>  [66,]     0     1     1     4     4     4
#>  [67,]     1     0     1     0     4     0
#>  [68,]     1     1     0     4     2     0
#>  [69,]     0     1     0     0     4     4
#>  [70,]     0     0     0     0     2     1
#>  [71,]     1     1     1     0     4     1
#>  [72,]     0     1     1     2     1     0
#>  [73,]     0     1     1     1     1     4
#>  [74,]     0     1     1     0     4     4
#>  [75,]     0     1     1     0     1     2
#>  [76,]     0     1     1     1     2     0
#>  [77,]     1     1     1     4     0     3
#>  [78,]     0     1     1     0     4     3
#>  [79,]     0     1     1     0     2     0
#>  [80,]     0     1     0     0     2     4
#>  [81,]     1     0     1     0     2     1
#>  [82,]     1     0     1     2     1     4
#>  [83,]     0     1     1     0     4     0
#>  [84,]     0     0     0     0     1     0
#>  [85,]     1     1     1     4     4     0
#>  [86,]     0     0     0     0     1     0
#>  [87,]     0     1     1     1     2     1
#>  [88,]     0     1     1     4     2     4
#>  [89,]     0     1     0     1     4     0
#>  [90,]     1     1     1     4     4     4
#>  [91,]     1     1     1     2     3     3
#>  [92,]     0     0     0     0     1     0
#>  [93,]     0     1     1     1     4     0
#>  [94,]     1     1     1     4     4     3
#>  [95,]     0     0     1     0     0     0
#>  [96,]     0     1     1     0     4     0
#>  [97,]     1     1     1     1     4     4
#>  [98,]     0     0     1     1     4     2
#>  [99,]     0     1     1     4     4     4
#> [100,]     0     1     0     3     3     2
```
