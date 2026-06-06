# Create a Data Frame of Item Metadata

This function creates a data frame of item metadata—including item
parameters, the number of score categories, and IRT model
specifications—to be used in various IRT-related analyses within the
irtQ package.

## Usage

``` r
shape_df(
  par.drm = list(a = NULL, b = NULL, g = NULL),
  par.prm = list(a = NULL, d = NULL),
  item.id = NULL,
  cats,
  model,
  default.par = FALSE
)
```

## Arguments

- par.drm:

  A list containing three numeric vectors for dichotomous item
  parameters: item discrimination (`a`), item difficulty (`b`), and
  guessing parameters (`g`).

- par.prm:

  A list containing polytomous item parameters. The list must include a
  numeric vector `a` for item discrimination (slope) parameters, and a
  list `d` of numeric vectors specifying difficulty (or threshold)
  parameters for each item. See the **Details** section for more
  information.

- item.id:

  A character vector of item IDs. If `NULL`, default IDs (e.g., "V1",
  "V2", ...) are assigned automatically.

- cats:

  A numeric vector indicating the number of score categories for each
  item.

- model:

  A character vector specifying the IRT model for each item. Available
  options are `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for dichotomous
  items, and `"GRM"` and `"GPCM"` for polytomous items. The label
  `"DRM"` serves as a general category that encompasses all dichotomous
  models (`"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and `"GPCM"`
  refer to the graded response model and (generalized) partial credit
  model, respectively.

- default.par:

  Logical. If `TRUE`, default item parameters are generated based on the
  specified `cats` and `model`. In this case, the slope parameter is set
  to 1, all difficulty (or threshold) parameters are set to 0, and the
  guessing parameter is set to 0.2 for `"3PLM"` or `"DRM"` items. The
  default is `FALSE`.

## Value

A data frame containing item metadata, including item IDs, number of
score categories, IRT model types, and associated item parameters. This
data frame can be used as input for other functions in the irtQ package,
such as
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md).

## Details

For any item where `"1PLM"` or `"2PLM"` is specified in `model`, the
guessing parameter will be set to `NA`. If `model` is a vector of length
1, the specified model will be replicated across all items.

As in the
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md)
function, when constructing a mixed-format test form, it is important to
specify the `cats` argument to reflect the correct number of score
categories for each item, in the exact order that the items appear. See
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
further guidance on how to specify `cats`.

When specifying item parameters using `par.drm` and/or `par.prm`, the
internal structure and ordering of elements must be followed.

- `par.drm` should be a list with three components:

  - `a`: a numeric vector of slope parameters

  - `b`: a numeric vector of difficulty parameters

  - `g`: a numeric vector of guessing parameters

- `par.prm` should be a list with two components:

  - `a`: a numeric vector of slope parameters for polytomous items

  - `d`: a list of numeric vectors specifying threshold (or step)
    parameters for each polytomous item

For items following the (generalized) partial credit model (`"GPCM"`),
the threshold (or step) parameters are computed as the overall item
difficulty (location) minus the category-specific thresholds. Therefore,
for an item with `m` score categories, `m - 1` step parameters must be
provided, since the first category threshold is fixed and does not
contribute to category probabilities.

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df_fipc()`](https://hwangQ.github.io/irtQ/reference/shape_df_fipc.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## A mixed-format test form
## containing five dichotomous items and two polytomous items
# Create a list of dichotomous item parameters
par.drm <- list(
  a = c(1.1, 1.2, 0.9, 1.8, 1.4),
  b = c(0.1, -1.6, -0.2, 1.0, 1.2),
  g = rep(0.2, 5)
)

# Create a list of polytomous item parameters
par.prm <- list(
  a = c(1.4, 0.6),
  d = list(
    c(0.0, -1.9, 1.2),
    c(0.4, -1.1, 1.5, 0.2)
  )
)

# Create a numeric vector indicating the number of score categories for each item
cats <- c(2, 4, 2, 2, 5, 2, 2)

# Create a character vector specifying the IRT model for each item
model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")

# Generate an item metadata set using the specified parameters
shape_df(par.drm = par.drm, par.prm = par.prm, cats = cats, model = model)
#>   id cats model par.1 par.2 par.3 par.4 par.5
#> 1 V1    2   DRM   1.1   0.1   0.2    NA    NA
#> 2 V2    4   GRM   1.4   0.0  -1.9   1.2    NA
#> 3 V3    2   DRM   1.2  -1.6   0.2    NA    NA
#> 4 V4    2   DRM   0.9  -0.2   0.2    NA    NA
#> 5 V5    5  GPCM   0.6   0.4  -1.1   1.5   0.2
#> 6 V6    2   DRM   1.8   1.0   0.2    NA    NA
#> 7 V7    2   DRM   1.4   1.2   0.2    NA    NA

## An empty item metadata frame with five dichotomous items and two polytomous items
# Create a numeric vector indicating the number of score categories for each item
cats <- c(2, 4, 3, 2, 5, 2, 2)

# Create a character vector specifying the IRT model for each item
model <- c("1PLM", "GRM", "GRM", "2PLM", "GPCM", "DRM", "3PLM")

# Generate an item metadata frame with default parameters
shape_df(cats = cats, model = model, default.par = TRUE)
#>   id cats model par.1 par.2 par.3 par.4 par.5
#> 1 V1    2  1PLM     1     0    NA    NA    NA
#> 2 V2    4   GRM     1     0   0.0     0    NA
#> 3 V3    3   GRM     1     0   0.0    NA    NA
#> 4 V4    2  2PLM     1     0    NA    NA    NA
#> 5 V5    5  GPCM     1     0   0.0     0     0
#> 6 V6    2   DRM     1     0   0.2    NA    NA
#> 7 V7    2  3PLM     1     0   0.2    NA    NA

## A single-format test form consisting of five dichotomous items
# Generate the item metadata
shape_df(par.drm = par.drm, cats = rep(2, 5), model = "DRM")
#>   id cats model par.1 par.2 par.3
#> 1 V1    2   DRM   1.1   0.1   0.2
#> 2 V2    2   DRM   1.2  -1.6   0.2
#> 3 V3    2   DRM   0.9  -0.2   0.2
#> 4 V4    2   DRM   1.8   1.0   0.2
#> 5 V5    2   DRM   1.4   1.2   0.2
```
