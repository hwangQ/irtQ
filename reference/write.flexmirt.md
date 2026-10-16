# Write a "-prm.txt" File for flexMIRT

This function writes a flexMIRT-compatible "-prm.txt" file (Cai, 2017).
It supports unidimensional IRT models only. This function was developed
by modifying `read.flexmirt()` from Pritikin & Falk (2020).

## Usage

``` r
write.flexmirt(
  x,
  file = NULL,
  norm.pop = c(0, 1),
  rePar = TRUE,
  mgroup = FALSE,
  group.name = NULL
)
```

## Arguments

- x:

  A data frame of item metadata (e.g., item parameters, number of
  categories, model types) for a single group, or a list of such data
  frames for multiple groups. See
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  item metadata format. You can also create metadata using
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md).

- file:

  A character string giving the path of the file to be written. It must
  be supplied. flexMIRT parameter files usually end in "-prm.txt".

- norm.pop:

  A numeric vector of length two giving the mean and the variance of the
  normal ability distribution of a single group, or a list of such
  vectors, one per group (e.g.,
  `norm.pop = list(c(0, 1), c(0, 0.64), c(0.5, 1.44))` for three
  groups). The second value is written to the variance field of the
  file. If `mgroup = TRUE` and a single vector is given, it is used for
  all groups. Default is `c(0, 1)`.

- rePar:

  Logical. Indicates how the parameters in `x` are given. If `TRUE`
  (default), `x` is taken to hold flexMIRT slope-intercept values:
  slopes, intercepts, and logits of the guessing parameters, and for
  GPCM items the scoring function and intercept values of the nominal
  model; they are written as they are. If `FALSE`, `x` is taken to hold
  irtQ item metadata (slopes, difficulty or step parameters, and
  guessing parameters), which is converted to the flexMIRT
  parameterization. Item metadata from irtQ (for example, from
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md) or
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md))
  must be written with `rePar = FALSE`. Note that the meaning is the
  reverse of the `rePar` argument of
  [`bring.flexmirt()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md),
  where `TRUE` converts the values of the file to irtQ metadata.

- mgroup:

  A logical value indicating whether the file includes multiple groups.
  Default is `FALSE`.

- group.name:

  A character vector of group names. If `NULL`, group names are
  automatically generated (e.g., "Group1", "Group2", ...).

## Value

Called for its side effect of writing the file to `file`. The value is
returned invisibly and is not meant to be used.

## Details

The function writes items of the models `"1PLM"`, `"2PLM"`, `"3PLM"`,
`"DRM"`, `"GRM"`, and `"GPCM"`. Dichotomous items are written as
flexMIRT 3PL (`"3PLM"` and `"DRM"`) or graded (`"1PLM"` and `"2PLM"`)
items, so
[`bring.flexmirt()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md)
reads them back as `"3PLM"` or `"2PLM"` items. Every item and every
group is written as one line, whatever the width of the console. The
parameters are written without rescaling. flexMIRT uses the logistic
metric (D = 1); when `x` was calibrated with another scaling constant,
multiply the slopes by that constant before writing.

## References

Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item
analysis and test scoring (Computer software). Chapel Hill, NC: Vector
Psychometric Group.

Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research
environment for item response theory method development. *Applied
Psychological Measurement, 44*(7-8), 561-562.
[doi:10.1177/0146621620929431](https://doi.org/10.1177/0146621620929431)
.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
## 1. Create a "-prm.txt" file for a single group
##    using the simulated CAT data
# 1-(1) Extract the item metadata
x <- simCAT_MX$item.prm

# 1-(2) Set the name of the "-prm.txt" file
temp_prm <- file.path(tempdir(), "single_group_temp-prm.txt")

# 1-(3) Write the "-prm.txt" file
write.flexmirt(x, file = temp_prm, norm.pop = c(0, 1), rePar = FALSE)

# 1-(4) Read the file back and compare the item parameters
x_back <- bring.flexmirt(file = temp_prm, type = "par")$Group1$full_df
head(x_back)
#>   id cats model     par.1       par.2 par.3 par.4
#> 1 V1    2  2PLM 0.9664715 -0.84085552    NA    NA
#> 2 V2    2  2PLM 0.9152754  1.38435940    NA    NA
#> 3 V3    2  2PLM 1.3454796 -1.25549187    NA    NA
#> 4 V4    2  2PLM 0.7925184  0.07014276    NA    NA
#> 5 V5    2  2PLM 1.0862914  1.71144078    NA    NA
#> 6 V6    2  2PLM 1.1311496 -0.60290796    NA    NA
max(abs(as.matrix(x[, 4:7]) - as.matrix(x_back[, 4:7])), na.rm = TRUE)
#> [1] 1.70888e-07

## 2. Create a "-prm.txt" file for multiple groups
##    using simulated multi-group data
# 2-(1) Extract the item metadata
x <- simMG$item.prm

# Set the name of the "-prm.txt" file
temp_prm <- file.path(tempdir(), "mg_group_temp-prm1.txt")

# Write the "-prm.txt" file, giving the mean and variance of each group
write.flexmirt(x,
  file = temp_prm, norm.pop = list(c(0, 1), c(0.5, 0.64), c(-0.3, 1.69)),
  rePar = FALSE, mgroup = TRUE, group.name = c("GR1", "GR2", "GR3")
)

# Or write the "-prm.txt" file so that
# all groups share the same ability distribution
# and group names are generated automatically
temp_prm <- file.path(tempdir(), "mg_group_temp-prm2.txt")
write.flexmirt(x,
  file = temp_prm, norm.pop = c(0, 1),
  rePar = FALSE, mgroup = TRUE, group.name = NULL
)
# }
```
