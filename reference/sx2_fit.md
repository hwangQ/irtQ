# S-X2 Fit Statistic

Computes the \\S\text{-}X^2\\ item fit statistic proposed by Orlando and
Thissen (2000, 2003). This statistic evaluates the fit of IRT models by
comparing observed and expected item response frequencies across summed
score groups.

## Usage

``` r
sx2_fit(x, ...)

# Default S3 method
sx2_fit(
  x,
  data,
  D = 1,
  alpha = 0.05,
  min.collapse = 1,
  norm.prior = c(0, 1),
  nquad = 30,
  weights,
  pcm.loc = NULL,
  ...
)

# S3 method for class 'est_item'
sx2_fit(
  x,
  alpha = 0.05,
  min.collapse = 1,
  norm.prior = c(0, 1),
  nquad = 30,
  weights,
  pcm.loc = NULL,
  ...
)

# S3 method for class 'est_irt'
sx2_fit(
  x,
  alpha = 0.05,
  min.collapse = 1,
  norm.prior = c(0, 1),
  nquad = 30,
  weights,
  pcm.loc = NULL,
  ...
)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.); or an object of class `est_irt`
  obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md), or
  `est_item` from
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

  See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
  or [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- ...:

  Additional arguments passed to or from other methods.

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for the
  hypothesis test associated with the \\S\text{-}X^2\\ statistic.
  Default is 0.05.

- min.collapse:

  An integer specifying the minimum expected frequency required per cell
  before adjacent cells are collapsed. Default is 1. See **Details**.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. These values are used to
  generate the Gaussian quadrature points and weights. Ignored if
  `method` is `"ML"`, `"MLF"`, `"WL"`, or `"INV.TCC"`. Default is
  `c(0, 1)`.

- nquad:

  An integer specifying the number of Gaussian quadrature points used to
  approximate the normal prior distribution. Default is 30.

- weights:

  A two-column matrix or data frame containing the quadrature points
  (first column) and their corresponding weights (second column) for the
  latent ability distribution. If omitted, default values are generated
  using
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)
  according to the `norm.prior` and `nquad` arguments.

- pcm.loc:

  An optional integer vector indicating the row indices of items that
  follow the partial credit model (PCM), where slope parameters are
  fixed. Default is `NULL`.

## Value

A list containing the following components:

- fit_stat:

  A data frame summarizing the \\S\text{-}X^2\\ fit statistics for all
  items, including the chi-square value, degrees of freedom, critical
  value, and p-value.

- item_df:

  A data frame containing the item metadata as specified in the input
  argument `x`.

- exp_freq:

  A list of collapsed expected frequency tables for all items.

- obs_freq:

  A list of collapsed observed frequency tables for all items.

- exp_prob:

  A list of collapsed expected probability tables for all items.

- obs_prop:

  A list of collapsed observed proportion tables for all items.

## Details

The accuracy of the \\\chi^{2}\\ approximation in item fit statistics
can be compromised when expected cell frequencies in contingency tables
are too small (Orlando & Thissen, 2000). To address this issue, Orlando
and Thissen (2000) proposed collapsing adjacent summed score groups to
ensure a minimum expected frequency of at least 1.

However, applying this collapsing approach directly to polytomous item
data can result in excessive information loss (Kang & Chen, 2008). To
mitigate this, Kang and Chen (2008) instead collapsed adjacent response
categories *within* each summed score group, maintaining a minimum
expected frequency of 1 per category. The same collapsing strategies are
implemented in `sx2_fit()`. If a different minimum expected frequency is
desired, it can be specified via the `min.collapse` argument.

When an item is labeled as "DRM" in the item metadata, it is treated as
a 3PLM item when computing the degrees of freedom for the
\\S\text{-}X^2\\ statistic.

Additionally, any missing responses in the `data` are automatically
replaced with incorrect responses (i.e., 0s).

## Methods (by class)

- `sx2_fit(default)`: Default method for computing \\S\text{-}X^{2}\\
  fit statistics from a data frame `x` containing item metadata.

- `sx2_fit(est_item)`: An object created by the function
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

- `sx2_fit(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

## References

Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X2 item
fit index for polytomous IRT models. *Journal of Educational
Measurement, 45*(4), 391-406.

Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for
dichotomous item response theory models. *Applied Psychological
Measurement, 24*(1), 50-64.

Orlando, M., & Thissen, D. (2003). Further investigation of the
performance of S-X2: An item fit index for use with dichotomous item
response theory models. *Applied Psychological Measurement, 27*(4),
289-298.

## See also

[`irtfit()`](https://hwangQ.github.io/irtQ/reference/irtfit.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1: All five polytomous IRT items follow the GRM
## Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select the item metadata
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df

# Generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(500, mean = 0, sd = 1)

# Simulate response data
data <- simdat(x = x, theta = score, D = 1)

# \donttest{
# Compute fit statistics
fit1 <- sx2_fit(x = x, data = data, nquad = 30)

# Display fit statistics
fit1$fit_stat
#>       id   chisq  df crit.val     p
#> 1   CMC1  41.820  46   62.830 0.648
#> 2   CMC2  33.471  30   43.773 0.302
#> 3   CMC3  44.763  44   60.481 0.440
#> 4   CMC4  38.579  43   59.304 0.663
#> 5   CMC5  37.867  44   60.481 0.731
#> 6   CMC6  30.877  37   52.192 0.751
#> 7   CMC7  35.587  44   60.481 0.813
#> 8   CMC8  41.343  44   60.481 0.586
#> 9   CMC9  34.567  45   61.656 0.870
#> 10 CMC10  45.242  40   55.758 0.262
#> 11 CMC11  38.145  43   59.304 0.682
#> 12 CMC12  49.043  44   60.481 0.278
#> 13 CMC13  49.212  42   58.124 0.207
#> 14 CMC14  42.408  42   58.124 0.453
#> 15 CMC15  37.455  42   58.124 0.671
#> 16 CMC16  27.217  32   46.194 0.707
#> 17 CMC17  45.594  40   55.758 0.251
#> 18 CMC18  58.619  44   60.481 0.069
#> 19 CMC19  20.195  27   40.113 0.822
#> 20 CMC20  27.962  27   40.113 0.413
#> 21 CMC21  52.977  32   46.194 0.011
#> 22 CMC22  34.690  44   60.481 0.841
#> 23 CMC23  41.056  43   59.304 0.556
#> 24 CMC24  61.967  46   62.830 0.058
#> 25 CMC25  43.106  40   55.758 0.340
#> 26 CMC26  31.319  32   46.194 0.501
#> 27 CMC27  56.343  42   58.124 0.069
#> 28 CMC28  25.142  35   49.802 0.891
#> 29 CMC29  45.948  32   46.194 0.053
#> 30 CMC30  37.810  43   59.304 0.695
#> 31 CMC31  54.912  43   59.304 0.105
#> 32 CMC32  37.874  35   49.802 0.340
#> 33 CMC33  38.690  32   46.194 0.193
#> 34 CMC34  45.301  41   56.942 0.297
#> 35 CMC35  45.151  39   54.572 0.230
#> 36 CMC36  48.687  44   60.481 0.290
#> 37 CMC37  40.247  36   50.998 0.288
#> 38 CMC38  44.503  45   61.656 0.493
#> 39  CFR1  97.733  88  110.898 0.224
#> 40  CFR2 108.445 120  146.567 0.767
#> 41  AMC1  46.444  42   58.124 0.294
#> 42  AMC2  27.622  24   36.415 0.276
#> 43  AMC3  41.031  40   55.758 0.425
#> 44  AMC4  34.027  42   58.124 0.805
#> 45  AMC5  45.964  45   61.656 0.432
#> 46  AMC6  40.593  43   59.304 0.576
#> 47  AMC7  60.867  42   58.124 0.030
#> 48  AMC8  54.508  41   56.942 0.077
#> 49  AMC9  30.830  40   55.758 0.851
#> 50 AMC10  42.111  41   56.942 0.423
#> 51 AMC11  29.044  32   46.194 0.617
#> 52 AMC12  51.016  41   56.942 0.136
#> 53  AFR1 135.326 119  145.461 0.145
#> 54  AFR2 121.217 110  135.480 0.219
#> 55  AFR3 141.069 123  149.885 0.127
# }

## Example 2: Items 39 and 40 follow the GRM, and items 53, 54, and 55
##            follow the PCM (with slope parameters fixed to 1)
# Replace the model names with "GPCM" and
# set the slope parameters of items 53–55 to 1
x[53:55, 3] <- "GPCM"
x[53:55, 4] <- 1

# Generate examinees' abilities from N(0, 1)
set.seed(25)
score <- rnorm(1000, mean = 0, sd = 1)

# Simulate response data
data <- simdat(x = x, theta = score, D = 1)

# \donttest{
# Compute fit statistics
fit2 <- sx2_fit(x = x, data = data, nquad = 30, pcm.loc = 53:55)

# Display fit statistics
fit2$fit_stat
#>       id   chisq  df crit.val     p
#> 1   CMC1  38.993  52   69.832 0.909
#> 2   CMC2  40.598  36   50.998 0.275
#> 3   CMC3  54.997  48   65.171 0.227
#> 4   CMC4  62.705  47   64.001 0.062
#> 5   CMC5  43.025  48   65.171 0.676
#> 6   CMC6  39.764  44   60.481 0.654
#> 7   CMC7  27.195  49   66.339 0.995
#> 8   CMC8  53.367  49   66.339 0.310
#> 9   CMC9  47.774  51   68.669 0.603
#> 10 CMC10  46.878  45   61.656 0.395
#> 11 CMC11  69.239  47   64.001 0.019
#> 12 CMC12  60.003  50   67.505 0.157
#> 13 CMC13  43.799  49   66.339 0.683
#> 14 CMC14  51.288  47   64.001 0.309
#> 15 CMC15  40.811  46   62.830 0.689
#> 16 CMC16  43.177  40   55.758 0.337
#> 17 CMC17  47.159  46   62.830 0.425
#> 18 CMC18  72.673  50   67.505 0.020
#> 19 CMC19  36.399  33   47.400 0.313
#> 20 CMC20  37.697  38   53.384 0.483
#> 21 CMC21  35.660  38   53.384 0.578
#> 22 CMC22  49.990  49   66.339 0.434
#> 23 CMC23  71.263  48   65.171 0.016
#> 24 CMC24  53.683  52   69.832 0.410
#> 25 CMC25  47.714  48   65.171 0.484
#> 26 CMC26  45.599  42   58.124 0.325
#> 27 CMC27  36.997  46   62.830 0.826
#> 28 CMC28  35.186  42   58.124 0.762
#> 29 CMC29  57.828  41   56.942 0.042
#> 30 CMC30  39.178  49   66.339 0.841
#> 31 CMC31  57.011  49   66.339 0.202
#> 32 CMC32  36.683  40   55.758 0.620
#> 33 CMC33  37.643  42   58.124 0.663
#> 34 CMC34  44.357  46   62.830 0.541
#> 35 CMC35  64.767  46   62.830 0.035
#> 36 CMC36  62.038  50   67.505 0.118
#> 37 CMC37  41.595  41   56.942 0.445
#> 38 CMC38  53.261  50   67.505 0.350
#> 39  CFR1 167.078 116  142.138 0.001
#> 40  CFR2 151.980 153  182.865 0.508
#> 41  AMC1  64.537  48   65.171 0.056
#> 42  AMC2  30.042  33   47.400 0.615
#> 43  AMC3  38.483  46   62.830 0.777
#> 44  AMC4  51.250  47   64.001 0.311
#> 45  AMC5  66.037  51   68.669 0.077
#> 46  AMC6  47.302  49   66.339 0.542
#> 47  AMC7  48.906  46   62.830 0.357
#> 48  AMC8  40.019  45   61.656 0.683
#> 49  AMC9  46.916  46   62.830 0.435
#> 50 AMC10  64.741  47   64.001 0.044
#> 51 AMC11  45.410  38   53.384 0.191
#> 52 AMC12  59.689  48   65.171 0.120
#> 53  AFR1 110.521 122  148.779 0.763
#> 54  AFR2 103.503 115  141.030 0.771
#> 55  AFR3 126.169 126  153.198 0.479
# }
```
