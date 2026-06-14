# Fixed ability parameter calibration

This function performs fixed ability parameter calibration (FAPC), often
called Stocking's (1988) Method A, which is the maximum likelihood
estimation of item parameters given ability estimates (Baker & Kim,
2004; Ban et al., 2001; Stocking, 1988). It can also be considered a
special case of joint maximum likelihood estimation in which only one
cycle of item parameter estimation is conducted, conditioned on the
given ability estimates (Birnbaum, 1968). FAPC is a potentially useful
method for calibrating pretest (or newly developed) items in
computerized adaptive testing (CAT), as it enables placing their
parameter estimates on the same scale as operational items. In addition,
it can be used to recalibrate operational items in the item bank to
evaluate potential parameter drift (Chen & Wang, 2016; Stocking, 1988).

## Usage

``` r
est_item(
  x = NULL,
  data,
  score,
  D = 1,
  model = NULL,
  cats = NULL,
  item.id = NULL,
  fix.a.1pl = FALSE,
  fix.a.gpcm = FALSE,
  fix.g = FALSE,
  a.val.1pl = 1,
  a.val.gpcm = 1,
  g.val = 0.2,
  use.aprior = FALSE,
  use.bprior = FALSE,
  use.gprior = TRUE,
  aprior = list(dist = "lnorm", params = c(0, 0.5)),
  bprior = list(dist = "norm", params = c(0, 1)),
  gprior = list(dist = "beta", params = c(5, 17)),
  missing = NA,
  use.startval = FALSE,
  control = list(eval.max = 500, iter.max = 200, x.tol = 1e-04),
  verbose = TRUE
)
```

## Arguments

- x:

  A data frame containing item metadata. This metadata is required to
  retrieve essential information for each item (e.g., number of score
  categories, IRT model type, etc.) necessary for calibration. You can
  create an empty item metadata frame using the function
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md).

  When `use.startval = TRUE`, the item parameters specified in the
  metadata will be used as starting values for parameter estimation. If
  `x = NULL`, both `model` and `cats` arguments must be specified. See
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. Default is `NULL`.

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- score:

  A numeric vector of examinees' ability estimates (theta values). The
  length of this vector must match the number of rows in the response
  data.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- model:

  A character vector specifying the IRT model to fit each item.
  Available values are:

  - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items

  - `"GRM"`, `"GPCM"` for polytomous items

  Here, `"GRM"` denotes the graded response model and `"GPCM"` the
  (generalized) partial credit model. Note that `"DRM"` serves as a
  general label covering all three dichotomous IRT models. If a single
  model name is provided, it is recycled for all items. This argument is
  only used when `x = NULL`. Default is `NULL`.

- cats:

  Numeric vector specifying the number of score categories per item. For
  dichotomous items, this should be 2. If a single value is supplied, it
  will be recycled across all items. When `cats = NULL` and all models
  specified in the `model` argument are dichotomous (`"1PLM"`, `"2PLM"`,
  `"3PLM"`, or `"DRM"`), the function defaults to 2 categories per item.
  This argument is used only when `x = NULL`. Default is `NULL`.

- item.id:

  Character vector of item identifiers. If `NULL`, IDs are generated
  automatically. When `fipc = TRUE`, a provided `item.id` will override
  any IDs present in `x`. Default is `NULL`.

- fix.a.1pl:

  Logical. If `TRUE`, the slope parameters of all 1PLM items are fixed
  to `a.val.1pl`; otherwise, they are constrained to be equal and
  estimated. Default is `FALSE`.

- fix.a.gpcm:

  Logical. If `TRUE`, GPCM items are calibrated as PCM with slopes fixed
  to `a.val.gpcm`; otherwise, each item's slope is estimated. Default is
  `FALSE`.

- fix.g:

  Logical. If `TRUE`, all 3PLM guessing parameters are fixed to `g.val`;
  otherwise, each guessing parameter is estimated. Default is `FALSE`.

- a.val.1pl:

  Numeric. Value to which the slope parameters of 1PLM items are fixed
  when `fix.a.1pl = TRUE`. Default is 1.

- a.val.gpcm:

  Numeric. Value to which the slope parameters of GPCM items are fixed
  when `fix.a.gpcm = TRUE`. Default is 1.

- g.val:

  Numeric. Value to which the guessing parameters of 3PLM items are
  fixed when `fix.g = TRUE`. Default is 0.2.

- use.aprior:

  Logical. If `TRUE`, applies a prior distribution to all item
  discrimination (slope) parameters during calibration. Default is
  `FALSE`.

- use.bprior:

  Logical. If `TRUE`, applies a prior distribution to all item
  difficulty (or threshold) parameters during calibration. Default is
  `FALSE`.

- use.gprior:

  Logical. If `TRUE`, applies a prior distribution to all 3PLM guessing
  parameters during calibration. Default is `TRUE`.

- aprior, bprior, gprior:

  A list specifying the prior distribution for all item discrimination
  (slope), difficulty (or threshold), guessing parameters. Three
  distributions are supported: Beta, Log-normal, and Normal. The list
  must have two elements:

  - `dist`: A character string, one of `"beta"`, `"lnorm"`, or `"norm"`.

  - `params`: A numeric vector of length two giving the distribution’s
    parameters. For details on each parameterization, see
    [`stats::dbeta()`](https://rdrr.io/r/stats/Beta.html),
    [`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html), and
    [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html).

  Defaults are:

  - `aprior = list(dist = "lnorm", params = c(0.0, 0.5))`

  - `bprior = list(dist = "norm", params = c(0.0, 1.0))`

  - `gprior = list(dist = "beta", params = c(5, 16))`

  for discrimination, difficulty, and guessing parameters, respectively.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- use.startval:

  Logical. If `TRUE`, the item parameters provided in the item metadata
  (i.e., the `x` argument) are used as starting values for item
  parameter estimation. Otherwise, internally generated starting values
  are used. Default is `FALSE`.

- control:

  A named list of options passed directly to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). These
  parameters define settings for the item parameter estimation process,
  such as the maximum number of iterations. By default:
  `control = list(eval.max = 500, iter.max = 200, x.tol = 1e-4)`, where

  - `eval.max` = 500 limits the number of function evaluations

  - `iter.max` = 200 caps the number of internal optimizer iterations

  - `x.tol` = 1e‑4 sets the absolute change threshold in parameter
    values below which
    [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) considers
    the solution to have converged Users may additionally supply other
    [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) control options
    (such as `abs.tol`, `rel.tol`, `trace`, etc.) as needed.

- verbose:

  Logical. If `FALSE`, all progress messages are suppressed. Default is
  `TRUE`.

## Value

This function returns an object of class `est_item`. The returned object
contains the following components:

- estimates:

  A data frame containing both the item parameter estimates and their
  corresponding standard errors.

- par.est:

  A data frame of item parameter estimates, structured according to the
  item metadata format.

- se.est:

  A data frame of standard errors for the item parameter estimates,
  computed based on the observed information functions

- pos.par:

  A data frame indicating the position of each item parameter within the
  estimation vector. Useful for interpreting the variance-covariance
  matrix.

- covariance:

  A variance-covariance matrix of the item parameter estimates.

- loglikelihood:

  The total log-likelihood value computed across all estimated items
  based on the complete response data.

- data:

  A data frame of examinees' response data.

- score:

  A vector of examinees' ability estimates used as fixed values during
  item parameter estimation.

- scale.D:

  The scaling factor used in the IRT model.

- convergence:

  A message indicating whether item parameter estimation successfully
  converged.

- nitem:

  The total number of items in the response data.

- deleted.item:

  Items with no response data. These items are excluded from the item
  parameter estimation.

- npar.est:

  The total number of parameters estimated.

- n.response:

  An integer vector indicating the number of valid responses for each
  item used in the item parameter estimation.

- TotalTime:

  Total computation time in seconds.

Note that you can easily extract components from the output using the
[`getirt()`](https://hwangQ.github.io/irtQ/reference/getirt.md)
function.

## Details

In most cases, the function `est_item()` returns successfully converged
item parameter estimates using its default internal starting values.
However, if convergence issues arise during calibration, one possible
solution is to use alternative starting values. If item parameter values
are already specified in the item metadata (i.e., the `x` argument),
they can be used as starting values for item parameter calibration by
setting `use.startval = TRUE`.

## References

Baker, F. B., & Kim, S. H. (2004). *Item response theory: Parameter
estimation techniques.* CRC Press.

Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A
comparative study of on-line pretest item calibration/scaling methods in
computerized adaptive testing. *Journal of Educational Measurement,
38*(3), 191-212.

Birnbaum, A. (1968). Some latent trait models and their use in inferring
an examinee's ability. In F. M. Lord & M. R. Novick (Eds.), *Statistical
theories of mental test scores* (pp. 397-479). Reading, MA:
Addison-Wesley.

Chen, P., & Wang, C. (2016). A new online calibration method for
multidimensional computerized adaptive testing. *Psychometrika, 81*(3),
674-701.

Stocking, M. L. (1988). *Scale drift in on-line calibration* (Research
Rep. 88-28). Princeton, NJ: ETS.

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`getirt()`](https://hwangQ.github.io/irtQ/reference/getirt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Extract the item metadata
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df

# Modify the item metadata so that some items follow 1PLM, 2PLM, and GPCM
x[c(1:3, 5), 3] <- "1PLM"
x[c(1:3, 5), 4] <- 1
x[c(1:3, 5), 6] <- 0
x[c(4, 8:12), 3] <- "2PLM"
x[c(4, 8:12), 6] <- 0
x[54:55, 3] <- "GPCM"

# Generate examinees' abilities from N(0, 1)
set.seed(23)
score <- rnorm(500, mean = 0, sd = 1)

# Simulate response data based on the item metadata and ability values
data <- simdat(x = x, theta = score, D = 1)

# \donttest{
# 1) Estimate item parameters: constrain the slope parameters of 1PLM items
#    to be equal
(mod1 <- est_item(x, data, score,
  D = 1, fix.a.1pl = FALSE, use.gprior = TRUE,
  gprior = list(dist = "beta", params = c(5, 17)), use.startval = FALSE
))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished. 
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = FALSE, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         17)), use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15830.66
#> 
summary(mod1)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = FALSE, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         17)), use.startval = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items in Response Data: 55
#>  Number of Excluded Items: 0
#>  Number of free parameters: 162
#>  Number of Responses for Each Item: 
#>        id    n
#> 1    CMC1  500
#> 2    CMC2  500
#> 3    CMC3  500
#> 4    CMC4  500
#> 5    CMC5  500
#> 6    CMC6  500
#> 7    CMC7  500
#> 8    CMC8  500
#> 9    CMC9  500
#> 10  CMC10  500
#> 11  CMC11  500
#> 12  CMC12  500
#> 13  CMC13  500
#> 14  CMC14  500
#> 15  CMC15  500
#> 16  CMC16  500
#> 17  CMC17  500
#> 18  CMC18  500
#> 19  CMC19  500
#> 20  CMC20  500
#> 21  CMC21  500
#> 22  CMC22  500
#> 23  CMC23  500
#> 24  CMC24  500
#> 25  CMC25  500
#> 26  CMC26  500
#> 27  CMC27  500
#> 28  CMC28  500
#> 29  CMC29  500
#> 30  CMC30  500
#> 31  CMC31  500
#> 32  CMC32  500
#> 33  CMC33  500
#> 34  CMC34  500
#> 35  CMC35  500
#> 36  CMC36  500
#> 37  CMC37  500
#> 38  CMC38  500
#> 39   CFR1  500
#> 40   CFR2  500
#> 41   AMC1  500
#> 42   AMC2  500
#> 43   AMC3  500
#> 44   AMC4  500
#> 45   AMC5  500
#> 46   AMC6  500
#> 47   AMC7  500
#> 48   AMC8  500
#> 49   AMC9  500
#> 50  AMC10  500
#> 51  AMC11  500
#> 52  AMC12  500
#> 53   AFR1  500
#> 54   AFR2  500
#> 55   AFR3  500
#> 
#> Processing time (in seconds) 
#>  Total computation: 0.44
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31661.31
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.02  0.06   1.60  0.13     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.02    NA  -1.06  0.12     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.02    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.02    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   1.88  0.27   0.67  0.09   0.10  0.03     NA    NA
#> 7    CMC7     2   3PLM   0.88  0.17   1.03  0.23   0.13  0.05     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.35  0.34   1.31  0.17   0.17  0.04     NA    NA
#> 14  CMC14     2   3PLM   1.36  0.31   0.15  0.24   0.24  0.08     NA    NA
#> 15  CMC15     2   3PLM   1.53  0.27   0.01  0.17   0.20  0.07     NA    NA
#> 16  CMC16     2   3PLM   2.10  0.25   0.04  0.08   0.10  0.04     NA    NA
#> 17  CMC17     2   3PLM   1.02  0.15  -0.41  0.22   0.16  0.07     NA    NA
#> 18  CMC18     2   3PLM   1.27  0.38   1.42  0.20   0.22  0.05     NA    NA
#> 19  CMC19     2   3PLM   2.25  0.32  -1.11  0.14   0.17  0.07     NA    NA
#> 20  CMC20     2   3PLM   1.47  0.22  -1.74  0.22   0.18  0.08     NA    NA
#> 21  CMC21     2   3PLM   1.38  0.21  -1.25  0.23   0.20  0.08     NA    NA
#> 22  CMC22     2   3PLM   0.92  0.16  -0.55  0.28   0.19  0.08     NA    NA
#> 23  CMC23     2   3PLM   1.10  0.22  -0.12  0.27   0.22  0.09     NA    NA
#> 24  CMC24     2   3PLM   1.21  0.34   1.43  0.21   0.22  0.05     NA    NA
#> 25  CMC25     2   3PLM   0.83  0.16  -1.51  0.40   0.21  0.09     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.18  -2.16  0.35   0.19  0.08     NA    NA
#> 27  CMC27     2   3PLM   1.18  0.18   0.09  0.17   0.14  0.06     NA    NA
#> 28  CMC28     2   3PLM   2.19  0.31  -0.17  0.11   0.19  0.05     NA    NA
#> 29  CMC29     2   3PLM   2.48  0.54  -0.81  0.20   0.38  0.09     NA    NA
#> 30  CMC30     2   3PLM   1.88  0.45   0.69  0.15   0.34  0.05     NA    NA
#> 31  CMC31     2   3PLM   0.70  0.16   1.00  0.32   0.16  0.07     NA    NA
#> 32  CMC32     2   3PLM   1.73  0.30  -0.78  0.21   0.26  0.09     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.17  -1.45  0.28   0.19  0.08     NA    NA
#> 34  CMC34     2   3PLM   1.04  0.17   0.21  0.20   0.16  0.06     NA    NA
#> 35  CMC35     2   3PLM   1.36  0.19  -0.45  0.17   0.16  0.06     NA    NA
#> 36  CMC36     2   3PLM   0.88  0.17   0.97  0.23   0.14  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.13  0.26  -0.25  0.09   0.13  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.87  0.17  -0.31  0.32   0.20  0.09     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.81  0.39   0.74  0.14   0.28  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.70  0.25  -1.59  0.20   0.19  0.08     NA    NA
#> 43   AMC3     2   3PLM   1.30  0.25   0.68  0.16   0.16  0.05     NA    NA
#> 44   AMC4     2   3PLM   0.94  0.17  -0.18  0.26   0.18  0.07     NA    NA
#> 45   AMC5     2   3PLM   1.69  0.65   2.11  0.26   0.19  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.83  0.64   1.44  0.10   0.15  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.69  0.41   0.37  0.18   0.25  0.07     NA    NA
#> 48   AMC8     2   3PLM   1.65  0.29   0.39  0.14   0.20  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.55  0.26   0.49  0.13   0.15  0.05     NA    NA
#> 50  AMC10     2   3PLM   2.48  0.51   1.31  0.10   0.13  0.02     NA    NA
#> 51  AMC11     2   3PLM   1.73  0.23  -1.02  0.15   0.16  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.95  0.20  -0.83  0.38   0.24  0.10     NA    NA
#> 53   AFR1     5    GRM   1.14  0.10  -0.30  0.09   0.30  0.09   0.92  0.11
#> 54   AFR2     5   GPCM   1.33  0.11  -1.99  0.21  -1.31  0.15  -0.72  0.12
#> 55   AFR3     5   GPCM   0.89  0.07  -0.80  0.15   0.15  0.15   0.46  0.16
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.23  0.06
#> 40   1.09  0.10
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.35  0.13
#> 54  -0.21  0.10
#> 55   1.35  0.19
#> 
#>  Group Parameters: 
#>    mu  sigma  
#>  0.03   1.02  
#> 

# Extract the item parameter estimates
getirt(mod1, what = "par.est")
#>       id cats model     par.1        par.2      par.3      par.4      par.5
#> 1   CMC1    2  1PLM 1.0185279  1.601332202         NA         NA         NA
#> 2   CMC2    2  1PLM 1.0185279 -1.063754095         NA         NA         NA
#> 3   CMC3    2  1PLM 1.0185279  0.396517677         NA         NA         NA
#> 4   CMC4    2  2PLM 0.9627020 -0.425593384         NA         NA         NA
#> 5   CMC5    2  1PLM 1.0185279 -0.250363868         NA         NA         NA
#> 6   CMC6    2  3PLM 1.8799184  0.671716274  0.1030566         NA         NA
#> 7   CMC7    2  3PLM 0.8814907  1.027532298  0.1344986         NA         NA
#> 8   CMC8    2  2PLM 0.9232153  0.865152105         NA         NA         NA
#> 9   CMC9    2  2PLM 1.0043493  0.892276674         NA         NA         NA
#> 10 CMC10    2  2PLM 1.6128253  0.092530217         NA         NA         NA
#> 11 CMC11    2  2PLM 1.0684709 -0.373791037         NA         NA         NA
#> 12 CMC12    2  2PLM 0.9377669  1.096465720         NA         NA         NA
#> 13 CMC13    2  3PLM 1.3517210  1.310901882  0.1679386         NA         NA
#> 14 CMC14    2  3PLM 1.3633891  0.148678889  0.2381220         NA         NA
#> 15 CMC15    2  3PLM 1.5260724  0.006963704  0.2014372         NA         NA
#> 16 CMC16    2  3PLM 2.0983557  0.039632481  0.1024467         NA         NA
#> 17 CMC17    2  3PLM 1.0182338 -0.410450918  0.1596037         NA         NA
#> 18 CMC18    2  3PLM 1.2670142  1.423951501  0.2216976         NA         NA
#> 19 CMC19    2  3PLM 2.2480817 -1.110500323  0.1713341         NA         NA
#> 20 CMC20    2  3PLM 1.4661194 -1.736037297  0.1808714         NA         NA
#> 21 CMC21    2  3PLM 1.3829158 -1.251268609  0.2003236         NA         NA
#> 22 CMC22    2  3PLM 0.9246240 -0.550518227  0.1901104         NA         NA
#> 23 CMC23    2  3PLM 1.1035700 -0.119211336  0.2161302         NA         NA
#> 24 CMC24    2  3PLM 1.2149156  1.433332519  0.2175768         NA         NA
#> 25 CMC25    2  3PLM 0.8260590 -1.506818565  0.2126189         NA         NA
#> 26 CMC26    2  3PLM 1.0651606 -2.157381967  0.1905922         NA         NA
#> 27 CMC27    2  3PLM 1.1784410  0.087234898  0.1431743         NA         NA
#> 28 CMC28    2  3PLM 2.1935809 -0.168868628  0.1910359         NA         NA
#> 29 CMC29    2  3PLM 2.4827570 -0.810059744  0.3832586         NA         NA
#> 30 CMC30    2  3PLM 1.8822455  0.693193770  0.3367413         NA         NA
#> 31 CMC31    2  3PLM 0.6992865  1.001048202  0.1588988         NA         NA
#> 32 CMC32    2  3PLM 1.7343873 -0.783189247  0.2575467         NA         NA
#> 33 CMC33    2  3PLM 1.0695257 -1.448547191  0.1902857         NA         NA
#> 34 CMC34    2  3PLM 1.0448545  0.210734195  0.1566960         NA         NA
#> 35 CMC35    2  3PLM 1.3588210 -0.447313405  0.1579927         NA         NA
#> 36 CMC36    2  3PLM 0.8813308  0.967252482  0.1390619         NA         NA
#> 37 CMC37    2  3PLM 2.1252131 -0.246234819  0.1285889         NA         NA
#> 38 CMC38    2  3PLM 0.8710093 -0.314317203  0.2036931         NA         NA
#> 39  CFR1    5   GRM 2.0000867 -1.882481406 -1.2545794 -0.7031275 -0.2324418
#> 40  CFR2    5   GRM 1.3885998 -0.796583085 -0.1288726  0.6012371  1.0859619
#> 41  AMC1    2  3PLM 1.8103262  0.735519975  0.2756669         NA         NA
#> 42  AMC2    2  3PLM 1.7047488 -1.593681876  0.1864312         NA         NA
#> 43  AMC3    2  3PLM 1.3037867  0.677356694  0.1580301         NA         NA
#> 44  AMC4    2  3PLM 0.9414845 -0.175796355  0.1844536         NA         NA
#> 45  AMC5    2  3PLM 1.6874360  2.107853698  0.1899503         NA         NA
#> 46  AMC6    2  3PLM 2.8301910  1.444309538  0.1536947         NA         NA
#> 47  AMC7    2  3PLM 1.6902164  0.371603972  0.2524593         NA         NA
#> 48  AMC8    2  3PLM 1.6456775  0.389081754  0.1999853         NA         NA
#> 49  AMC9    2  3PLM 1.5475943  0.487734230  0.1516041         NA         NA
#> 50 AMC10    2  3PLM 2.4831968  1.310828955  0.1287528         NA         NA
#> 51 AMC11    2  3PLM 1.7344800 -1.020301934  0.1636095         NA         NA
#> 52 AMC12    2  3PLM 0.9513776 -0.828671559  0.2401560         NA         NA
#> 53  AFR1    5   GRM 1.1355459 -0.300165699  0.3028394  0.9182795  1.3525597
#> 54  AFR2    5  GPCM 1.3267609 -1.994449888 -1.3127161 -0.7196856 -0.2068705
#> 55  AFR3    5  GPCM 0.8941751 -0.795427791  0.1540945  0.4638485  1.3511384

# 2) Estimate item parameters: fix the slope parameters of 1PLM items to 1
(mod2 <- est_item(x, data, score,
  D = 1, fix.a.1pl = TRUE, a.val.1pl = 1, use.gprior = TRUE,
  gprior = list(dist = "beta", params = c(5, 17)), use.startval = FALSE
))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished. 
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>         params = c(5, 17)), use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15830.7
#> 
summary(mod2)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>         params = c(5, 17)), use.startval = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items in Response Data: 55
#>  Number of Excluded Items: 0
#>  Number of free parameters: 161
#>  Number of Responses for Each Item: 
#>        id    n
#> 1    CMC1  500
#> 2    CMC2  500
#> 3    CMC3  500
#> 4    CMC4  500
#> 5    CMC5  500
#> 6    CMC6  500
#> 7    CMC7  500
#> 8    CMC8  500
#> 9    CMC9  500
#> 10  CMC10  500
#> 11  CMC11  500
#> 12  CMC12  500
#> 13  CMC13  500
#> 14  CMC14  500
#> 15  CMC15  500
#> 16  CMC16  500
#> 17  CMC17  500
#> 18  CMC18  500
#> 19  CMC19  500
#> 20  CMC20  500
#> 21  CMC21  500
#> 22  CMC22  500
#> 23  CMC23  500
#> 24  CMC24  500
#> 25  CMC25  500
#> 26  CMC26  500
#> 27  CMC27  500
#> 28  CMC28  500
#> 29  CMC29  500
#> 30  CMC30  500
#> 31  CMC31  500
#> 32  CMC32  500
#> 33  CMC33  500
#> 34  CMC34  500
#> 35  CMC35  500
#> 36  CMC36  500
#> 37  CMC37  500
#> 38  CMC38  500
#> 39   CFR1  500
#> 40   CFR2  500
#> 41   AMC1  500
#> 42   AMC2  500
#> 43   AMC3  500
#> 44   AMC4  500
#> 45   AMC5  500
#> 46   AMC6  500
#> 47   AMC7  500
#> 48   AMC8  500
#> 49   AMC9  500
#> 50  AMC10  500
#> 51  AMC11  500
#> 52  AMC12  500
#> 53   AFR1  500
#> 54   AFR2  500
#> 55   AFR3  500
#> 
#> Processing time (in seconds) 
#>  Total computation: 0.44
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31661.4
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.00    NA   1.62  0.12     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.00    NA  -1.08  0.11     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.00    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.00    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   1.88  0.27   0.67  0.09   0.10  0.03     NA    NA
#> 7    CMC7     2   3PLM   0.88  0.17   1.03  0.23   0.13  0.05     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.35  0.34   1.31  0.17   0.17  0.04     NA    NA
#> 14  CMC14     2   3PLM   1.36  0.31   0.15  0.24   0.24  0.08     NA    NA
#> 15  CMC15     2   3PLM   1.53  0.27   0.01  0.17   0.20  0.07     NA    NA
#> 16  CMC16     2   3PLM   2.10  0.25   0.04  0.08   0.10  0.04     NA    NA
#> 17  CMC17     2   3PLM   1.02  0.15  -0.41  0.22   0.16  0.07     NA    NA
#> 18  CMC18     2   3PLM   1.27  0.38   1.42  0.20   0.22  0.05     NA    NA
#> 19  CMC19     2   3PLM   2.25  0.32  -1.11  0.14   0.17  0.07     NA    NA
#> 20  CMC20     2   3PLM   1.47  0.22  -1.74  0.22   0.18  0.08     NA    NA
#> 21  CMC21     2   3PLM   1.38  0.21  -1.25  0.23   0.20  0.08     NA    NA
#> 22  CMC22     2   3PLM   0.92  0.16  -0.55  0.28   0.19  0.08     NA    NA
#> 23  CMC23     2   3PLM   1.10  0.22  -0.12  0.27   0.22  0.09     NA    NA
#> 24  CMC24     2   3PLM   1.21  0.34   1.43  0.21   0.22  0.05     NA    NA
#> 25  CMC25     2   3PLM   0.83  0.16  -1.51  0.40   0.21  0.09     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.18  -2.16  0.35   0.19  0.08     NA    NA
#> 27  CMC27     2   3PLM   1.18  0.18   0.09  0.17   0.14  0.06     NA    NA
#> 28  CMC28     2   3PLM   2.19  0.31  -0.17  0.11   0.19  0.05     NA    NA
#> 29  CMC29     2   3PLM   2.48  0.54  -0.81  0.20   0.38  0.09     NA    NA
#> 30  CMC30     2   3PLM   1.88  0.45   0.69  0.15   0.34  0.05     NA    NA
#> 31  CMC31     2   3PLM   0.70  0.16   1.00  0.32   0.16  0.07     NA    NA
#> 32  CMC32     2   3PLM   1.73  0.30  -0.78  0.21   0.26  0.09     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.17  -1.45  0.28   0.19  0.08     NA    NA
#> 34  CMC34     2   3PLM   1.04  0.17   0.21  0.20   0.16  0.06     NA    NA
#> 35  CMC35     2   3PLM   1.36  0.19  -0.45  0.17   0.16  0.06     NA    NA
#> 36  CMC36     2   3PLM   0.88  0.17   0.97  0.23   0.14  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.13  0.26  -0.25  0.09   0.13  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.87  0.17  -0.31  0.32   0.20  0.09     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.81  0.39   0.74  0.14   0.28  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.70  0.25  -1.59  0.20   0.19  0.08     NA    NA
#> 43   AMC3     2   3PLM   1.30  0.25   0.68  0.16   0.16  0.05     NA    NA
#> 44   AMC4     2   3PLM   0.94  0.17  -0.18  0.26   0.18  0.07     NA    NA
#> 45   AMC5     2   3PLM   1.69  0.65   2.11  0.26   0.19  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.83  0.64   1.44  0.10   0.15  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.69  0.41   0.37  0.18   0.25  0.07     NA    NA
#> 48   AMC8     2   3PLM   1.65  0.29   0.39  0.14   0.20  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.55  0.26   0.49  0.13   0.15  0.05     NA    NA
#> 50  AMC10     2   3PLM   2.48  0.51   1.31  0.10   0.13  0.02     NA    NA
#> 51  AMC11     2   3PLM   1.73  0.23  -1.02  0.15   0.16  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.95  0.20  -0.83  0.38   0.24  0.10     NA    NA
#> 53   AFR1     5    GRM   1.14  0.10  -0.30  0.09   0.30  0.09   0.92  0.11
#> 54   AFR2     5   GPCM   1.33  0.11  -1.99  0.21  -1.31  0.15  -0.72  0.12
#> 55   AFR3     5   GPCM   0.89  0.07  -0.80  0.15   0.15  0.15   0.46  0.16
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.23  0.06
#> 40   1.09  0.10
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.35  0.13
#> 54  -0.21  0.10
#> 55   1.35  0.19
#> 
#>  Group Parameters: 
#>    mu  sigma  
#>  0.03   1.02  
#> 

# Extract the standard error estimates
getirt(mod2, what = "se.est")
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   CMC1    2  1PLM         NA 0.11819685         NA         NA         NA
#> 2   CMC2    2  1PLM         NA 0.10825953         NA         NA         NA
#> 3   CMC3    2  1PLM         NA 0.09953834         NA         NA         NA
#> 4   CMC4    2  2PLM 0.11608725 0.11102431         NA         NA         NA
#> 5   CMC5    2  1PLM         NA 0.09928963         NA         NA         NA
#> 6   CMC6    2  3PLM 0.26592875 0.09223933 0.03126001         NA         NA
#> 7   CMC7    2  3PLM 0.17007216 0.22868364 0.05043985         NA         NA
#> 8   CMC8    2  2PLM 0.11779263 0.13411409         NA         NA         NA
#> 9   CMC9    2  2PLM 0.12280227 0.12623451         NA         NA         NA
#> 10 CMC10    2  2PLM 0.15351879 0.06743190         NA         NA         NA
#> 11 CMC11    2  2PLM 0.12147067 0.09986851         NA         NA         NA
#> 12 CMC12    2  2PLM 0.12170130 0.14960502         NA         NA         NA
#> 13 CMC13    2  3PLM 0.33543956 0.16922355 0.04476102         NA         NA
#> 14 CMC14    2  3PLM 0.31025933 0.23821549 0.08293693         NA         NA
#> 15 CMC15    2  3PLM 0.26615652 0.17322453 0.06684887         NA         NA
#> 16 CMC16    2  3PLM 0.24782824 0.08378912 0.03671343         NA         NA
#> 17 CMC17    2  3PLM 0.15336572 0.21797789 0.06861050         NA         NA
#> 18 CMC18    2  3PLM 0.37951683 0.20475893 0.05152830         NA         NA
#> 19 CMC19    2  3PLM 0.31575318 0.13854698 0.07248859         NA         NA
#> 20 CMC20    2  3PLM 0.21769602 0.22444579 0.07986432         NA         NA
#> 21 CMC21    2  3PLM 0.20918112 0.22563675 0.08356158         NA         NA
#> 22 CMC22    2  3PLM 0.15632987 0.27810767 0.07985413         NA         NA
#> 23 CMC23    2  3PLM 0.21812942 0.27395715 0.08641957         NA         NA
#> 24 CMC24    2  3PLM 0.33644706 0.21099243 0.04799483         NA         NA
#> 25 CMC25    2  3PLM 0.15706343 0.40435671 0.09407396         NA         NA
#> 26 CMC26    2  3PLM 0.18446052 0.34745444 0.08484835         NA         NA
#> 27 CMC27    2  3PLM 0.17645756 0.16891064 0.05657179         NA         NA
#> 28 CMC28    2  3PLM 0.31487828 0.10951574 0.05247206         NA         NA
#> 29 CMC29    2  3PLM 0.54086125 0.19720209 0.08789680         NA         NA
#> 30 CMC30    2  3PLM 0.44551859 0.15180550 0.04975419         NA         NA
#> 31 CMC31    2  3PLM 0.15959759 0.32384772 0.06560236         NA         NA
#> 32 CMC32    2  3PLM 0.30476277 0.21153602 0.08816621         NA         NA
#> 33 CMC33    2  3PLM 0.16888358 0.27816088 0.08298992         NA         NA
#> 34 CMC34    2  3PLM 0.16969347 0.19516848 0.05936051         NA         NA
#> 35 CMC35    2  3PLM 0.18844135 0.16666456 0.06387891         NA         NA
#> 36 CMC36    2  3PLM 0.17198087 0.23078825 0.05312022         NA         NA
#> 37 CMC37    2  3PLM 0.25992108 0.09456053 0.04537614         NA         NA
#> 38 CMC38    2  3PLM 0.16843257 0.31736545 0.08539489         NA         NA
#> 39  CFR1    5   GRM 0.14447773 0.11978029 0.08307182 0.06307867 0.05703100
#> 40  CFR2    5   GRM 0.10798168 0.09234859 0.07429798 0.08110680 0.09995812
#> 41  AMC1    2  3PLM 0.38772085 0.14321180 0.04820526         NA         NA
#> 42  AMC2    2  3PLM 0.25187572 0.19813671 0.08045500         NA         NA
#> 43  AMC3    2  3PLM 0.24782407 0.16413179 0.05395429         NA         NA
#> 44  AMC4    2  3PLM 0.16675130 0.25922208 0.07498779         NA         NA
#> 45  AMC5    2  3PLM 0.64524089 0.26375779 0.03021145         NA         NA
#> 46  AMC6    2  3PLM 0.63645971 0.09751124 0.02218408         NA         NA
#> 47  AMC7    2  3PLM 0.40524702 0.18202698 0.06837224         NA         NA
#> 48  AMC8    2  3PLM 0.29111839 0.13833114 0.05205826         NA         NA
#> 49  AMC9    2  3PLM 0.25875454 0.13232181 0.04830486         NA         NA
#> 50 AMC10    2  3PLM 0.50849974 0.09515461 0.02399167         NA         NA
#> 51 AMC11    2  3PLM 0.22749532 0.15400610 0.06996146         NA         NA
#> 52 AMC12    2  3PLM 0.19683738 0.38464459 0.10437009         NA         NA
#> 53  AFR1    5   GRM 0.10188574 0.09221751 0.08931676 0.10902552 0.13282858
#> 54  AFR2    5  GPCM 0.10884827 0.21458427 0.15227381 0.11752786 0.10041751
#> 55  AFR3    5  GPCM 0.07361344 0.15373354 0.15429938 0.16257565 0.18586806

# 3) Estimate item parameters: fix the guessing parameters of 3PLM items to 0.2
(mod3 <- est_item(x, data, score,
  D = 1, fix.a.1pl = TRUE, fix.g = TRUE, a.val.1pl = 1, g.val = .2,
  use.startval = FALSE
))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished. 
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     fix.g = TRUE, a.val.1pl = 1, g.val = 0.2, use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15916.26
#> 
summary(mod3)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     fix.g = TRUE, a.val.1pl = 1, g.val = 0.2, use.startval = FALSE)
#> 
#> Summary of the Data 
#>  Number of Items in Response Data: 55
#>  Number of Excluded Items: 0
#>  Number of free parameters: 121
#>  Number of Responses for Each Item: 
#>        id    n
#> 1    CMC1  500
#> 2    CMC2  500
#> 3    CMC3  500
#> 4    CMC4  500
#> 5    CMC5  500
#> 6    CMC6  500
#> 7    CMC7  500
#> 8    CMC8  500
#> 9    CMC9  500
#> 10  CMC10  500
#> 11  CMC11  500
#> 12  CMC12  500
#> 13  CMC13  500
#> 14  CMC14  500
#> 15  CMC15  500
#> 16  CMC16  500
#> 17  CMC17  500
#> 18  CMC18  500
#> 19  CMC19  500
#> 20  CMC20  500
#> 21  CMC21  500
#> 22  CMC22  500
#> 23  CMC23  500
#> 24  CMC24  500
#> 25  CMC25  500
#> 26  CMC26  500
#> 27  CMC27  500
#> 28  CMC28  500
#> 29  CMC29  500
#> 30  CMC30  500
#> 31  CMC31  500
#> 32  CMC32  500
#> 33  CMC33  500
#> 34  CMC34  500
#> 35  CMC35  500
#> 36  CMC36  500
#> 37  CMC37  500
#> 38  CMC38  500
#> 39   CFR1  500
#> 40   CFR2  500
#> 41   AMC1  500
#> 42   AMC2  500
#> 43   AMC3  500
#> 44   AMC4  500
#> 45   AMC5  500
#> 46   AMC6  500
#> 47   AMC7  500
#> 48   AMC8  500
#> 49   AMC9  500
#> 50  AMC10  500
#> 51  AMC11  500
#> 52  AMC12  500
#> 53   AFR1  500
#> 54   AFR2  500
#> 55   AFR3  500
#> 
#> Processing time (in seconds) 
#>  Total computation: 0.26
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31832.52
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.00    NA   1.62  0.12     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.00    NA  -1.08  0.11     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.00    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.00    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   2.25  0.29   0.83  0.08   0.20    NA     NA    NA
#> 7    CMC7     2   3PLM   1.00  0.17   1.22  0.19   0.20    NA     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.53  0.28   1.37  0.15   0.20    NA     NA    NA
#> 14  CMC14     2   3PLM   1.26  0.18   0.05  0.10   0.20    NA     NA    NA
#> 15  CMC15     2   3PLM   1.52  0.20   0.00  0.09   0.20    NA     NA    NA
#> 16  CMC16     2   3PLM   2.36  0.27   0.18  0.07   0.20    NA     NA    NA
#> 17  CMC17     2   3PLM   1.06  0.15  -0.30  0.12   0.20    NA     NA    NA
#> 18  CMC18     2   3PLM   1.16  0.23   1.38  0.19   0.20    NA     NA    NA
#> 19  CMC19     2   3PLM   2.30  0.30  -1.07  0.10   0.20    NA     NA    NA
#> 20  CMC20     2   3PLM   1.48  0.22  -1.71  0.19   0.20    NA     NA    NA
#> 21  CMC21     2   3PLM   1.38  0.19  -1.25  0.16   0.20    NA     NA    NA
#> 22  CMC22     2   3PLM   0.93  0.14  -0.52  0.15   0.20    NA     NA    NA
#> 23  CMC23     2   3PLM   1.08  0.16  -0.17  0.12   0.20    NA     NA    NA
#> 24  CMC24     2   3PLM   1.13  0.23   1.40  0.19   0.20    NA     NA    NA
#> 25  CMC25     2   3PLM   0.82  0.15  -1.55  0.28   0.20    NA     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.18  -2.14  0.31   0.20    NA     NA    NA
#> 27  CMC27     2   3PLM   1.27  0.17   0.22  0.10   0.20    NA     NA    NA
#> 28  CMC28     2   3PLM   2.22  0.27  -0.15  0.07   0.20    NA     NA    NA
#> 29  CMC29     2   3PLM   1.84  0.28  -1.19  0.14   0.20    NA     NA    NA
#> 30  CMC30     2   3PLM   1.19  0.20   0.35  0.11   0.20    NA     NA    NA
#> 31  CMC31     2   3PLM   0.76  0.15   1.15  0.23   0.20    NA     NA    NA
#> 32  CMC32     2   3PLM   1.62  0.22  -0.90  0.12   0.20    NA     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.16  -1.43  0.21   0.20    NA     NA    NA
#> 34  CMC34     2   3PLM   1.11  0.16   0.33  0.12   0.20    NA     NA    NA
#> 35  CMC35     2   3PLM   1.42  0.18  -0.36  0.10   0.20    NA     NA    NA
#> 36  CMC36     2   3PLM   1.00  0.17   1.15  0.18   0.20    NA     NA    NA
#> 37  CMC37     2   3PLM   2.30  0.27  -0.15  0.07   0.20    NA     NA    NA
#> 38  CMC38     2   3PLM   0.87  0.14  -0.33  0.15   0.20    NA     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.45  0.22   0.57  0.10   0.20    NA     NA    NA
#> 42   AMC2     2   3PLM   1.72  0.25  -1.57  0.16   0.20    NA     NA    NA
#> 43   AMC3     2   3PLM   1.44  0.21   0.77  0.11   0.20    NA     NA    NA
#> 44   AMC4     2   3PLM   0.96  0.15  -0.13  0.13   0.20    NA     NA    NA
#> 45   AMC5     2   3PLM   1.83  0.53   2.10  0.25   0.20    NA     NA    NA
#> 46   AMC6     2   3PLM   3.32  0.68   1.50  0.09   0.20    NA     NA    NA
#> 47   AMC7     2   3PLM   1.48  0.21   0.25  0.09   0.20    NA     NA    NA
#> 48   AMC8     2   3PLM   1.65  0.23   0.39  0.09   0.20    NA     NA    NA
#> 49   AMC9     2   3PLM   1.72  0.23   0.59  0.09   0.20    NA     NA    NA
#> 50  AMC10     2   3PLM   3.21  0.62   1.40  0.09   0.20    NA     NA    NA
#> 51  AMC11     2   3PLM   1.78  0.22  -0.96  0.11   0.20    NA     NA    NA
#> 52  AMC12     2   3PLM   0.91  0.15  -0.96  0.19   0.20    NA     NA    NA
#> 53   AFR1     5    GRM   1.14  0.10  -0.30  0.09   0.30  0.09   0.92  0.11
#> 54   AFR2     5   GPCM   1.33  0.11  -1.99  0.21  -1.31  0.15  -0.72  0.12
#> 55   AFR3     5   GPCM   0.89  0.07  -0.80  0.15   0.15  0.15   0.46  0.16
#>     par.5  se.5
#> 1      NA    NA
#> 2      NA    NA
#> 3      NA    NA
#> 4      NA    NA
#> 5      NA    NA
#> 6      NA    NA
#> 7      NA    NA
#> 8      NA    NA
#> 9      NA    NA
#> 10     NA    NA
#> 11     NA    NA
#> 12     NA    NA
#> 13     NA    NA
#> 14     NA    NA
#> 15     NA    NA
#> 16     NA    NA
#> 17     NA    NA
#> 18     NA    NA
#> 19     NA    NA
#> 20     NA    NA
#> 21     NA    NA
#> 22     NA    NA
#> 23     NA    NA
#> 24     NA    NA
#> 25     NA    NA
#> 26     NA    NA
#> 27     NA    NA
#> 28     NA    NA
#> 29     NA    NA
#> 30     NA    NA
#> 31     NA    NA
#> 32     NA    NA
#> 33     NA    NA
#> 34     NA    NA
#> 35     NA    NA
#> 36     NA    NA
#> 37     NA    NA
#> 38     NA    NA
#> 39  -0.23  0.06
#> 40   1.09  0.10
#> 41     NA    NA
#> 42     NA    NA
#> 43     NA    NA
#> 44     NA    NA
#> 45     NA    NA
#> 46     NA    NA
#> 47     NA    NA
#> 48     NA    NA
#> 49     NA    NA
#> 50     NA    NA
#> 51     NA    NA
#> 52     NA    NA
#> 53   1.35  0.13
#> 54  -0.21  0.10
#> 55   1.35  0.19
#> 
#>  Group Parameters: 
#>    mu  sigma  
#>  0.03   1.02  
#> 

# Extract both item parameter and standard error estimates
getirt(mod2, what = "estimates")
#>       id cats model     par.1       se.1        par.2       se.2      par.3
#> 1   CMC1    2  1PLM 1.0000000         NA  1.621648179 0.11819685         NA
#> 2   CMC2    2  1PLM 1.0000000         NA -1.077997364 0.10825953         NA
#> 3   CMC3    2  1PLM 1.0000000         NA  0.400985408 0.09953834         NA
#> 4   CMC4    2  2PLM 0.9627020 0.11608725 -0.425593384 0.11102431         NA
#> 5   CMC5    2  1PLM 1.0000000         NA -0.254140164 0.09928963         NA
#> 6   CMC6    2  3PLM 1.8799184 0.26592875  0.671716274 0.09223933  0.1030566
#> 7   CMC7    2  3PLM 0.8814907 0.17007216  1.027532298 0.22868364  0.1344986
#> 8   CMC8    2  2PLM 0.9232153 0.11779263  0.865152105 0.13411409         NA
#> 9   CMC9    2  2PLM 1.0043493 0.12280227  0.892276674 0.12623451         NA
#> 10 CMC10    2  2PLM 1.6128253 0.15351879  0.092530217 0.06743190         NA
#> 11 CMC11    2  2PLM 1.0684709 0.12147067 -0.373791037 0.09986851         NA
#> 12 CMC12    2  2PLM 0.9377669 0.12170130  1.096465720 0.14960502         NA
#> 13 CMC13    2  3PLM 1.3517210 0.33543956  1.310901882 0.16922355  0.1679386
#> 14 CMC14    2  3PLM 1.3633891 0.31025933  0.148678889 0.23821549  0.2381220
#> 15 CMC15    2  3PLM 1.5260724 0.26615652  0.006963704 0.17322453  0.2014372
#> 16 CMC16    2  3PLM 2.0983557 0.24782824  0.039632481 0.08378912  0.1024467
#> 17 CMC17    2  3PLM 1.0182338 0.15336572 -0.410450918 0.21797789  0.1596037
#> 18 CMC18    2  3PLM 1.2670142 0.37951683  1.423951501 0.20475893  0.2216976
#> 19 CMC19    2  3PLM 2.2480817 0.31575318 -1.110500323 0.13854698  0.1713341
#> 20 CMC20    2  3PLM 1.4661194 0.21769602 -1.736037297 0.22444579  0.1808714
#> 21 CMC21    2  3PLM 1.3829158 0.20918112 -1.251268609 0.22563675  0.2003236
#> 22 CMC22    2  3PLM 0.9246240 0.15632987 -0.550518227 0.27810767  0.1901104
#> 23 CMC23    2  3PLM 1.1035700 0.21812942 -0.119211336 0.27395715  0.2161302
#> 24 CMC24    2  3PLM 1.2149156 0.33644706  1.433332519 0.21099243  0.2175768
#> 25 CMC25    2  3PLM 0.8260590 0.15706343 -1.506818565 0.40435671  0.2126189
#> 26 CMC26    2  3PLM 1.0651606 0.18446052 -2.157381967 0.34745444  0.1905922
#> 27 CMC27    2  3PLM 1.1784410 0.17645756  0.087234898 0.16891064  0.1431743
#> 28 CMC28    2  3PLM 2.1935809 0.31487828 -0.168868628 0.10951574  0.1910359
#> 29 CMC29    2  3PLM 2.4827570 0.54086125 -0.810059744 0.19720209  0.3832586
#> 30 CMC30    2  3PLM 1.8822455 0.44551859  0.693193770 0.15180550  0.3367413
#> 31 CMC31    2  3PLM 0.6992865 0.15959759  1.001048202 0.32384772  0.1588988
#> 32 CMC32    2  3PLM 1.7343873 0.30476277 -0.783189247 0.21153602  0.2575467
#> 33 CMC33    2  3PLM 1.0695257 0.16888358 -1.448547191 0.27816088  0.1902857
#> 34 CMC34    2  3PLM 1.0448545 0.16969347  0.210734195 0.19516848  0.1566960
#> 35 CMC35    2  3PLM 1.3588210 0.18844135 -0.447313405 0.16666456  0.1579927
#> 36 CMC36    2  3PLM 0.8813308 0.17198087  0.967252482 0.23078825  0.1390619
#> 37 CMC37    2  3PLM 2.1252131 0.25992108 -0.246234819 0.09456053  0.1285889
#> 38 CMC38    2  3PLM 0.8710093 0.16843257 -0.314317203 0.31736545  0.2036931
#> 39  CFR1    5   GRM 2.0000867 0.14447773 -1.882481406 0.11978029 -1.2545794
#> 40  CFR2    5   GRM 1.3885998 0.10798168 -0.796583085 0.09234859 -0.1288726
#> 41  AMC1    2  3PLM 1.8103262 0.38772085  0.735519975 0.14321180  0.2756669
#> 42  AMC2    2  3PLM 1.7047488 0.25187572 -1.593681876 0.19813671  0.1864312
#> 43  AMC3    2  3PLM 1.3037867 0.24782407  0.677356694 0.16413179  0.1580301
#> 44  AMC4    2  3PLM 0.9414845 0.16675130 -0.175796355 0.25922208  0.1844536
#> 45  AMC5    2  3PLM 1.6874360 0.64524089  2.107853698 0.26375779  0.1899503
#> 46  AMC6    2  3PLM 2.8301910 0.63645971  1.444309538 0.09751124  0.1536947
#> 47  AMC7    2  3PLM 1.6902164 0.40524702  0.371603972 0.18202698  0.2524593
#> 48  AMC8    2  3PLM 1.6456775 0.29111839  0.389081754 0.13833114  0.1999853
#> 49  AMC9    2  3PLM 1.5475943 0.25875454  0.487734230 0.13232181  0.1516041
#> 50 AMC10    2  3PLM 2.4831968 0.50849974  1.310828955 0.09515461  0.1287528
#> 51 AMC11    2  3PLM 1.7344800 0.22749532 -1.020301934 0.15400610  0.1636095
#> 52 AMC12    2  3PLM 0.9513776 0.19683738 -0.828671559 0.38464459  0.2401560
#> 53  AFR1    5   GRM 1.1355459 0.10188574 -0.300165699 0.09221751  0.3028394
#> 54  AFR2    5  GPCM 1.3267609 0.10884827 -1.994449888 0.21458427 -1.3127161
#> 55  AFR3    5  GPCM 0.8941751 0.07361344 -0.795427791 0.15373354  0.1540945
#>          se.3      par.4       se.4      par.5       se.5
#> 1          NA         NA         NA         NA         NA
#> 2          NA         NA         NA         NA         NA
#> 3          NA         NA         NA         NA         NA
#> 4          NA         NA         NA         NA         NA
#> 5          NA         NA         NA         NA         NA
#> 6  0.03126001         NA         NA         NA         NA
#> 7  0.05043985         NA         NA         NA         NA
#> 8          NA         NA         NA         NA         NA
#> 9          NA         NA         NA         NA         NA
#> 10         NA         NA         NA         NA         NA
#> 11         NA         NA         NA         NA         NA
#> 12         NA         NA         NA         NA         NA
#> 13 0.04476102         NA         NA         NA         NA
#> 14 0.08293693         NA         NA         NA         NA
#> 15 0.06684887         NA         NA         NA         NA
#> 16 0.03671343         NA         NA         NA         NA
#> 17 0.06861050         NA         NA         NA         NA
#> 18 0.05152830         NA         NA         NA         NA
#> 19 0.07248859         NA         NA         NA         NA
#> 20 0.07986432         NA         NA         NA         NA
#> 21 0.08356158         NA         NA         NA         NA
#> 22 0.07985413         NA         NA         NA         NA
#> 23 0.08641957         NA         NA         NA         NA
#> 24 0.04799483         NA         NA         NA         NA
#> 25 0.09407396         NA         NA         NA         NA
#> 26 0.08484835         NA         NA         NA         NA
#> 27 0.05657179         NA         NA         NA         NA
#> 28 0.05247206         NA         NA         NA         NA
#> 29 0.08789680         NA         NA         NA         NA
#> 30 0.04975419         NA         NA         NA         NA
#> 31 0.06560236         NA         NA         NA         NA
#> 32 0.08816621         NA         NA         NA         NA
#> 33 0.08298992         NA         NA         NA         NA
#> 34 0.05936051         NA         NA         NA         NA
#> 35 0.06387891         NA         NA         NA         NA
#> 36 0.05312022         NA         NA         NA         NA
#> 37 0.04537614         NA         NA         NA         NA
#> 38 0.08539489         NA         NA         NA         NA
#> 39 0.08307182 -0.7031275 0.06307867 -0.2324418 0.05703100
#> 40 0.07429798  0.6012371 0.08110680  1.0859619 0.09995812
#> 41 0.04820526         NA         NA         NA         NA
#> 42 0.08045500         NA         NA         NA         NA
#> 43 0.05395429         NA         NA         NA         NA
#> 44 0.07498779         NA         NA         NA         NA
#> 45 0.03021145         NA         NA         NA         NA
#> 46 0.02218408         NA         NA         NA         NA
#> 47 0.06837224         NA         NA         NA         NA
#> 48 0.05205826         NA         NA         NA         NA
#> 49 0.04830486         NA         NA         NA         NA
#> 50 0.02399167         NA         NA         NA         NA
#> 51 0.06996146         NA         NA         NA         NA
#> 52 0.10437009         NA         NA         NA         NA
#> 53 0.08931676  0.9182795 0.10902552  1.3525597 0.13282858
#> 54 0.15227381 -0.7196856 0.11752786 -0.2068705 0.10041751
#> 55 0.15429938  0.4638485 0.16257565  1.3511384 0.18586806

# }
```
