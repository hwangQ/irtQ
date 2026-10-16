# Create a Data Frame of Item Metadata

This function creates a data frame of item metadata (item parameters,
numbers of score categories, and IRT models) for use in IRT analyses
with the irtQ package.

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
  guessing parameters (`g`), each with one value per item with
  `cats = 2`. A single value of `g` is used for every dichotomous item.
  If `g` is omitted or `NULL`, the guessing parameters are set to 0.

- par.prm:

  A list containing polytomous item parameters. The list must include a
  numeric vector `a` of slope parameters and a list `d` of numeric
  vectors of threshold parameters, with one slope and one vector of
  `cats - 1` thresholds per item with `cats > 2`. See the **Details**
  section for more information.

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

A data frame of item metadata with columns `id`, `cats`, `model`, and
`par.1`, `par.2`, ..., one row per item. The number of parameter columns
is the larger of 3 and the maximum of `cats`; unused cells are `NA`.
This data frame can be used as input for other functions in the irtQ
package, such as
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md).

## Details

For any item where `"1PLM"` or `"2PLM"` is specified in `model`, the
guessing parameter is set to `NA`. If `cats` or `model` has length 1, it
is recycled across all items.

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

  - `d`: a list of numeric vectors specifying the threshold parameters
    for each polytomous item

For GPCM items, the threshold parameters are \\b_v = \beta - \tau_v\\,
the overall item location minus the threshold of each score category. An
item with *K* score categories requires *K - 1* threshold parameters,
because the term for the lowest score category is fixed at 0 and is not
supplied. For GRM items, the thresholds must be in increasing order.

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
    c(-1.9, 0.0, 1.2),
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
#> 2 V2    4   GRM   1.4  -1.9   0.0   1.2    NA
#> 3 V3    2   DRM   1.2  -1.6   0.2    NA    NA
#> 4 V4    2   DRM   0.9  -0.2   0.2    NA    NA
#> 5 V5    5  GPCM   0.6   0.4  -1.1   1.5   0.2
#> 6 V6    2   DRM   1.8   1.0   0.2    NA    NA
#> 7 V7    2   DRM   1.4   1.2   0.2    NA    NA

## An item metadata frame with default parameters for five dichotomous and two polytomous items
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
