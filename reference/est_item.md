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
  gprior = list(dist = "beta", params = c(5, 16)),
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

  - `params`: A numeric vector of length two giving the distribution's
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

  - `x.tol` = 1e-4 sets the absolute change threshold in parameter
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

Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D. J. (2001). A
comparative study of on-line pretest item calibration/scaling methods in
computerized adaptive testing. *Journal of Educational Measurement,
38*(3), 191-212.
[doi:10.1111/j.1745-3984.2001.tb01123.x](https://doi.org/10.1111/j.1745-3984.2001.tb01123.x)
.

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
  gprior = list(dist = "beta", params = c(5, 16)), use.startval = FALSE
))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished. 
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = FALSE, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15832.81
#> 
summary(mod1)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = FALSE, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), use.startval = FALSE)
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
#>  Total computation: 0.47
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31665.62
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.02  0.06   1.60  0.13     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.02    NA  -1.06  0.12     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.02    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.02    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   1.88  0.27   0.67  0.09   0.10  0.03     NA    NA
#> 7    CMC7     2   3PLM   0.89  0.17   1.04  0.23   0.14  0.05     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.36  0.34   1.32  0.17   0.17  0.04     NA    NA
#> 14  CMC14     2   3PLM   1.39  0.32   0.17  0.24   0.25  0.08     NA    NA
#> 15  CMC15     2   3PLM   1.54  0.27   0.02  0.17   0.21  0.07     NA    NA
#> 16  CMC16     2   3PLM   2.10  0.25   0.04  0.08   0.10  0.04     NA    NA
#> 17  CMC17     2   3PLM   1.02  0.16  -0.40  0.22   0.17  0.07     NA    NA
#> 18  CMC18     2   3PLM   1.29  0.38   1.43  0.20   0.23  0.05     NA    NA
#> 19  CMC19     2   3PLM   2.26  0.32  -1.10  0.14   0.18  0.07     NA    NA
#> 20  CMC20     2   3PLM   1.47  0.22  -1.72  0.23   0.19  0.08     NA    NA
#> 21  CMC21     2   3PLM   1.39  0.21  -1.23  0.23   0.21  0.09     NA    NA
#> 22  CMC22     2   3PLM   0.93  0.16  -0.53  0.28   0.20  0.08     NA    NA
#> 23  CMC23     2   3PLM   1.12  0.23  -0.09  0.28   0.23  0.09     NA    NA
#> 24  CMC24     2   3PLM   1.23  0.34   1.44  0.21   0.22  0.05     NA    NA
#> 25  CMC25     2   3PLM   0.83  0.16  -1.47  0.41   0.22  0.10     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.19  -2.14  0.35   0.20  0.09     NA    NA
#> 27  CMC27     2   3PLM   1.18  0.18   0.10  0.17   0.15  0.06     NA    NA
#> 28  CMC28     2   3PLM   2.20  0.32  -0.16  0.11   0.19  0.05     NA    NA
#> 29  CMC29     2   3PLM   2.53  0.55  -0.79  0.19   0.40  0.09     NA    NA
#> 30  CMC30     2   3PLM   1.90  0.45   0.70  0.15   0.34  0.05     NA    NA
#> 31  CMC31     2   3PLM   0.71  0.16   1.02  0.33   0.16  0.07     NA    NA
#> 32  CMC32     2   3PLM   1.76  0.31  -0.76  0.21   0.27  0.09     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.17  -1.43  0.28   0.20  0.09     NA    NA
#> 34  CMC34     2   3PLM   1.05  0.17   0.22  0.20   0.16  0.06     NA    NA
#> 35  CMC35     2   3PLM   1.37  0.19  -0.44  0.17   0.16  0.07     NA    NA
#> 36  CMC36     2   3PLM   0.89  0.17   0.98  0.23   0.14  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.13  0.26  -0.24  0.09   0.13  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.88  0.17  -0.28  0.32   0.21  0.09     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.83  0.39   0.74  0.14   0.28  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.71  0.25  -1.58  0.20   0.19  0.08     NA    NA
#> 43   AMC3     2   3PLM   1.31  0.25   0.69  0.16   0.16  0.05     NA    NA
#> 44   AMC4     2   3PLM   0.95  0.17  -0.15  0.26   0.19  0.08     NA    NA
#> 45   AMC5     2   3PLM   1.70  0.65   2.11  0.26   0.19  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.84  0.64   1.44  0.10   0.15  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.72  0.41   0.39  0.18   0.26  0.07     NA    NA
#> 48   AMC8     2   3PLM   1.66  0.29   0.40  0.14   0.20  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.56  0.26   0.49  0.13   0.15  0.05     NA    NA
#> 50  AMC10     2   3PLM   2.49  0.51   1.31  0.10   0.13  0.02     NA    NA
#> 51  AMC11     2   3PLM   1.74  0.23  -1.01  0.16   0.17  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.97  0.21  -0.78  0.40   0.26  0.11     NA    NA
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
#>       id cats model     par.1       par.2      par.3      par.4      par.5
#> 1   CMC1    2  1PLM 1.0185279  1.60133220         NA         NA         NA
#> 2   CMC2    2  1PLM 1.0185279 -1.06375410         NA         NA         NA
#> 3   CMC3    2  1PLM 1.0185279  0.39651768         NA         NA         NA
#> 4   CMC4    2  2PLM 0.9627020 -0.42559338         NA         NA         NA
#> 5   CMC5    2  1PLM 1.0185279 -0.25036387         NA         NA         NA
#> 6   CMC6    2  3PLM 1.8846619  0.67350821  0.1041510         NA         NA
#> 7   CMC7    2  3PLM 0.8866580  1.03619252  0.1374716         NA         NA
#> 8   CMC8    2  2PLM 0.9232153  0.86515210         NA         NA         NA
#> 9   CMC9    2  2PLM 1.0043493  0.89227667         NA         NA         NA
#> 10 CMC10    2  2PLM 1.6128253  0.09253022         NA         NA         NA
#> 11 CMC11    2  2PLM 1.0684709 -0.37379104         NA         NA         NA
#> 12 CMC12    2  2PLM 0.9377669  1.09646572         NA         NA         NA
#> 13 CMC13    2  3PLM 1.3645089  1.31503341  0.1703375         NA         NA
#> 14 CMC14    2  3PLM 1.3898925  0.17217473  0.2472123         NA         NA
#> 15 CMC15    2  3PLM 1.5408435  0.01937278  0.2070791         NA         NA
#> 16 CMC16    2  3PLM 2.1026954  0.04186170  0.1039632         NA         NA
#> 17 CMC17    2  3PLM 1.0234932 -0.39553054  0.1653727         NA         NA
#> 18 CMC18    2  3PLM 1.2856593  1.43007415  0.2250732         NA         NA
#> 19 CMC19    2  3PLM 2.2594756 -1.10181795  0.1778769         NA         NA
#> 20 CMC20    2  3PLM 1.4700031 -1.72399703  0.1889865         NA         NA
#> 21 CMC21    2  3PLM 1.3912594 -1.23367452  0.2093723         NA         NA
#> 22 CMC22    2  3PLM 0.9317523 -0.52661586  0.1982627         NA         NA
#> 23 CMC23    2  3PLM 1.1204103 -0.09122518  0.2259506         NA         NA
#> 24 CMC24    2  3PLM 1.2290146  1.43876383  0.2204891         NA         NA
#> 25 CMC25    2  3PLM 0.8330125 -1.46925302  0.2244715         NA         NA
#> 26 CMC26    2  3PLM 1.0676986 -2.13985943  0.1999090         NA         NA
#> 27 CMC27    2  3PLM 1.1843144  0.09628328  0.1469717         NA         NA
#> 28 CMC28    2  3PLM 2.2044237 -0.16348433  0.1944618         NA         NA
#> 29 CMC29    2  3PLM 2.5307240 -0.78717016  0.3956548         NA         NA
#> 30 CMC30    2  3PLM 1.9041795  0.70161781  0.3404332         NA         NA
#> 31 CMC31    2  3PLM 0.7062845  1.02018205  0.1641299         NA         NA
#> 32 CMC32    2  3PLM 1.7564751 -0.76184344  0.2681714         NA         NA
#> 33 CMC33    2  3PLM 1.0744198 -1.42868292  0.1991516         NA         NA
#> 34 CMC34    2  3PLM 1.0509163  0.22205598  0.1609429         NA         NA
#> 35 CMC35    2  3PLM 1.3654004 -0.43712558  0.1629476         NA         NA
#> 36 CMC36    2  3PLM 0.8871087  0.97700579  0.1423811         NA         NA
#> 37 CMC37    2  3PLM 2.1311743 -0.24278137  0.1309809         NA         NA
#> 38 CMC38    2  3PLM 0.8813430 -0.28321005  0.2131754         NA         NA
#> 39  CFR1    5   GRM 2.0000867 -1.88248141 -1.2545794 -0.7031275 -0.2324418
#> 40  CFR2    5   GRM 1.3885998 -0.79658308 -0.1288726  0.6012371  1.0859619
#> 41  AMC1    2  3PLM 1.8272686  0.74245570  0.2788644         NA         NA
#> 42  AMC2    2  3PLM 1.7110855 -1.58183908  0.1947024         NA         NA
#> 43  AMC3    2  3PLM 1.3145385  0.68535647  0.1615132         NA         NA
#> 44  AMC4    2  3PLM 0.9496703 -0.15473831  0.1915267         NA         NA
#> 45  AMC5    2  3PLM 1.7029567  2.10717153  0.1910672         NA         NA
#> 46  AMC6    2  3PLM 2.8374610  1.44493420  0.1542768         NA         NA
#> 47  AMC7    2  3PLM 1.7188353  0.38574744  0.2586265         NA         NA
#> 48  AMC8    2  3PLM 1.6576325  0.39605794  0.2033722         NA         NA
#> 49  AMC9    2  3PLM 1.5568580  0.49337767  0.1543679         NA         NA
#> 50 AMC10    2  3PLM 2.4905705  1.31162574  0.1294143         NA         NA
#> 51 AMC11    2  3PLM 1.7412715 -1.01102653  0.1696460         NA         NA
#> 52 AMC12    2  3PLM 0.9688530 -0.77937288  0.2552984         NA         NA
#> 53  AFR1    5   GRM 1.1355459 -0.30016570  0.3028394  0.9182795  1.3525597
#> 54  AFR2    5  GPCM 1.3267609 -1.99444989 -1.3127161 -0.7196856 -0.2068705
#> 55  AFR3    5  GPCM 0.8941751 -0.79542779  0.1540945  0.4638485  1.3511384

# 2) Estimate item parameters: fix the slope parameters of 1PLM items to 1
(mod2 <- est_item(x, data, score,
  D = 1, fix.a.1pl = TRUE, a.val.1pl = 1, use.gprior = TRUE,
  gprior = list(dist = "beta", params = c(5, 16)), use.startval = FALSE
))
#> Starting... 
#> Parsing input... 
#> Estimating item parameters... 
#> Estimation is finished. 
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>         params = c(5, 16)), use.startval = FALSE)
#> 
#> Fixed ability parameter calibration (Stocking's Method A). 
#> All item parameters were successfully converged. 
#> 
#> Log-likelihood: -15832.85
#> 
summary(mod2)
#> 
#> Call:
#> est_item(x = x, data = data, score = score, D = 1, fix.a.1pl = TRUE, 
#>     a.val.1pl = 1, use.gprior = TRUE, gprior = list(dist = "beta", 
#>         params = c(5, 16)), use.startval = FALSE)
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
#>  Total computation: 0.47
#> 
#> Convergence of Solution 
#>  All item parameters were successfully converged.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 31665.71
#>  Item Parameters: 
#>        id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1    CMC1     2   1PLM   1.00    NA   1.62  0.12     NA    NA     NA    NA
#> 2    CMC2     2   1PLM   1.00    NA  -1.08  0.11     NA    NA     NA    NA
#> 3    CMC3     2   1PLM   1.00    NA   0.40  0.10     NA    NA     NA    NA
#> 4    CMC4     2   2PLM   0.96  0.12  -0.43  0.11     NA    NA     NA    NA
#> 5    CMC5     2   1PLM   1.00    NA  -0.25  0.10     NA    NA     NA    NA
#> 6    CMC6     2   3PLM   1.88  0.27   0.67  0.09   0.10  0.03     NA    NA
#> 7    CMC7     2   3PLM   0.89  0.17   1.04  0.23   0.14  0.05     NA    NA
#> 8    CMC8     2   2PLM   0.92  0.12   0.87  0.13     NA    NA     NA    NA
#> 9    CMC9     2   2PLM   1.00  0.12   0.89  0.13     NA    NA     NA    NA
#> 10  CMC10     2   2PLM   1.61  0.15   0.09  0.07     NA    NA     NA    NA
#> 11  CMC11     2   2PLM   1.07  0.12  -0.37  0.10     NA    NA     NA    NA
#> 12  CMC12     2   2PLM   0.94  0.12   1.10  0.15     NA    NA     NA    NA
#> 13  CMC13     2   3PLM   1.36  0.34   1.32  0.17   0.17  0.04     NA    NA
#> 14  CMC14     2   3PLM   1.39  0.32   0.17  0.24   0.25  0.08     NA    NA
#> 15  CMC15     2   3PLM   1.54  0.27   0.02  0.17   0.21  0.07     NA    NA
#> 16  CMC16     2   3PLM   2.10  0.25   0.04  0.08   0.10  0.04     NA    NA
#> 17  CMC17     2   3PLM   1.02  0.16  -0.40  0.22   0.17  0.07     NA    NA
#> 18  CMC18     2   3PLM   1.29  0.38   1.43  0.20   0.23  0.05     NA    NA
#> 19  CMC19     2   3PLM   2.26  0.32  -1.10  0.14   0.18  0.07     NA    NA
#> 20  CMC20     2   3PLM   1.47  0.22  -1.72  0.23   0.19  0.08     NA    NA
#> 21  CMC21     2   3PLM   1.39  0.21  -1.23  0.23   0.21  0.09     NA    NA
#> 22  CMC22     2   3PLM   0.93  0.16  -0.53  0.28   0.20  0.08     NA    NA
#> 23  CMC23     2   3PLM   1.12  0.23  -0.09  0.28   0.23  0.09     NA    NA
#> 24  CMC24     2   3PLM   1.23  0.34   1.44  0.21   0.22  0.05     NA    NA
#> 25  CMC25     2   3PLM   0.83  0.16  -1.47  0.41   0.22  0.10     NA    NA
#> 26  CMC26     2   3PLM   1.07  0.19  -2.14  0.35   0.20  0.09     NA    NA
#> 27  CMC27     2   3PLM   1.18  0.18   0.10  0.17   0.15  0.06     NA    NA
#> 28  CMC28     2   3PLM   2.20  0.32  -0.16  0.11   0.19  0.05     NA    NA
#> 29  CMC29     2   3PLM   2.53  0.55  -0.79  0.19   0.40  0.09     NA    NA
#> 30  CMC30     2   3PLM   1.90  0.45   0.70  0.15   0.34  0.05     NA    NA
#> 31  CMC31     2   3PLM   0.71  0.16   1.02  0.33   0.16  0.07     NA    NA
#> 32  CMC32     2   3PLM   1.76  0.31  -0.76  0.21   0.27  0.09     NA    NA
#> 33  CMC33     2   3PLM   1.07  0.17  -1.43  0.28   0.20  0.09     NA    NA
#> 34  CMC34     2   3PLM   1.05  0.17   0.22  0.20   0.16  0.06     NA    NA
#> 35  CMC35     2   3PLM   1.37  0.19  -0.44  0.17   0.16  0.07     NA    NA
#> 36  CMC36     2   3PLM   0.89  0.17   0.98  0.23   0.14  0.05     NA    NA
#> 37  CMC37     2   3PLM   2.13  0.26  -0.24  0.09   0.13  0.05     NA    NA
#> 38  CMC38     2   3PLM   0.88  0.17  -0.28  0.32   0.21  0.09     NA    NA
#> 39   CFR1     5    GRM   2.00  0.14  -1.88  0.12  -1.25  0.08  -0.70  0.06
#> 40   CFR2     5    GRM   1.39  0.11  -0.80  0.09  -0.13  0.07   0.60  0.08
#> 41   AMC1     2   3PLM   1.83  0.39   0.74  0.14   0.28  0.05     NA    NA
#> 42   AMC2     2   3PLM   1.71  0.25  -1.58  0.20   0.19  0.08     NA    NA
#> 43   AMC3     2   3PLM   1.31  0.25   0.69  0.16   0.16  0.05     NA    NA
#> 44   AMC4     2   3PLM   0.95  0.17  -0.15  0.26   0.19  0.08     NA    NA
#> 45   AMC5     2   3PLM   1.70  0.65   2.11  0.26   0.19  0.03     NA    NA
#> 46   AMC6     2   3PLM   2.84  0.64   1.44  0.10   0.15  0.02     NA    NA
#> 47   AMC7     2   3PLM   1.72  0.41   0.39  0.18   0.26  0.07     NA    NA
#> 48   AMC8     2   3PLM   1.66  0.29   0.40  0.14   0.20  0.05     NA    NA
#> 49   AMC9     2   3PLM   1.56  0.26   0.49  0.13   0.15  0.05     NA    NA
#> 50  AMC10     2   3PLM   2.49  0.51   1.31  0.10   0.13  0.02     NA    NA
#> 51  AMC11     2   3PLM   1.74  0.23  -1.01  0.16   0.17  0.07     NA    NA
#> 52  AMC12     2   3PLM   0.97  0.21  -0.78  0.40   0.26  0.11     NA    NA
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
#> 6   CMC6    2  3PLM 0.26669060 0.09232536 0.03138242         NA         NA
#> 7   CMC7    2  3PLM 0.17171857 0.22999854 0.05092018         NA         NA
#> 8   CMC8    2  2PLM 0.11779263 0.13411409         NA         NA         NA
#> 9   CMC9    2  2PLM 0.12280227 0.12623451         NA         NA         NA
#> 10 CMC10    2  2PLM 0.15351879 0.06743190         NA         NA         NA
#> 11 CMC11    2  2PLM 0.12147067 0.09986851         NA         NA         NA
#> 12 CMC12    2  2PLM 0.12170130 0.14960502         NA         NA         NA
#> 13 CMC13    2  3PLM 0.33830056 0.16875186 0.04452707         NA         NA
#> 14 CMC14    2  3PLM 0.32184562 0.23754432 0.08294538         NA         NA
#> 15 CMC15    2  3PLM 0.27056505 0.17343671 0.06714419         NA         NA
#> 16 CMC16    2  3PLM 0.24879579 0.08406154 0.03704144         NA         NA
#> 17 CMC17    2  3PLM 0.15550550 0.22209049 0.07041155         NA         NA
#> 18 CMC18    2  3PLM 0.38393374 0.20356474 0.05087416         NA         NA
#> 19 CMC19    2  3PLM 0.32044597 0.14051114 0.07448965         NA         NA
#> 20 CMC20    2  3PLM 0.21920822 0.22689431 0.08281388         NA         NA
#> 21 CMC21    2  3PLM 0.21215536 0.22895293 0.08608492         NA         NA
#> 22 CMC22    2  3PLM 0.15932023 0.28436101 0.08224986         NA         NA
#> 23 CMC23    2  3PLM 0.22727469 0.27935617 0.08847770         NA         NA
#> 24 CMC24    2  3PLM 0.33771111 0.21032077 0.04739131         NA         NA
#> 25 CMC25    2  3PLM 0.16035815 0.41455949 0.09844718         NA         NA
#> 26 CMC26    2  3PLM 0.18549120 0.35014363 0.08835993         NA         NA
#> 27 CMC27    2  3PLM 0.17838470 0.17043885 0.05738007         NA         NA
#> 28 CMC28    2  3PLM 0.31674844 0.10953967 0.05270285         NA         NA
#> 29 CMC29    2  3PLM 0.54619578 0.19184713 0.08610236         NA         NA
#> 30 CMC30    2  3PLM 0.44751927 0.15026890 0.04908313         NA         NA
#> 31 CMC31    2  3PLM 0.16316408 0.32764850 0.06684690         NA         NA
#> 32 CMC32    2  3PLM 0.31082457 0.21160840 0.08879304         NA         NA
#> 33 CMC33    2  3PLM 0.17077250 0.28271978 0.08601695         NA         NA
#> 34 CMC34    2  3PLM 0.17161904 0.19700301 0.06017396         NA         NA
#> 35 CMC35    2  3PLM 0.19069948 0.16865435 0.06511108         NA         NA
#> 36 CMC36    2  3PLM 0.17401848 0.23224001 0.05368437         NA         NA
#> 37 CMC37    2  3PLM 0.26144273 0.09498708 0.04587033         NA         NA
#> 38 CMC38    2  3PLM 0.17350226 0.32455119 0.08786132         NA         NA
#> 39  CFR1    5   GRM 0.14447773 0.11978029 0.08307182 0.06307867 0.05703100
#> 40  CFR2    5   GRM 0.10798168 0.09234859 0.07429798 0.08110680 0.09995812
#> 41  AMC1    2  3PLM 0.39127793 0.14263810 0.04793931         NA         NA
#> 42  AMC2    2  3PLM 0.25426145 0.20042801 0.08320344         NA         NA
#> 43  AMC3    2  3PLM 0.25146602 0.16446443 0.05423568         NA         NA
#> 44  AMC4    2  3PLM 0.17005564 0.26323511 0.07656965         NA         NA
#> 45  AMC5    2  3PLM 0.64633195 0.26189935 0.02992857         NA         NA
#> 46  AMC6    2  3PLM 0.63744855 0.09747301 0.02219971         NA         NA
#> 47  AMC7    2  3PLM 0.41094702 0.17861264 0.06713238         NA         NA
#> 48  AMC8    2  3PLM 0.29286642 0.13785201 0.05193618         NA         NA
#> 49  AMC9    2  3PLM 0.26091218 0.13241328 0.04845912         NA         NA
#> 50 AMC10    2  3PLM 0.50988663 0.09509815 0.02401244         NA         NA
#> 51 AMC11    2  3PLM 0.22999194 0.15624576 0.07189865         NA         NA
#> 52 AMC12    2  3PLM 0.20770051 0.39869458 0.10907235         NA         NA
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
#>  Total computation: 0.27
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
#>       id cats model     par.1       se.1       par.2       se.2      par.3
#> 1   CMC1    2  1PLM 1.0000000         NA  1.62164818 0.11819685         NA
#> 2   CMC2    2  1PLM 1.0000000         NA -1.07799736 0.10825953         NA
#> 3   CMC3    2  1PLM 1.0000000         NA  0.40098541 0.09953834         NA
#> 4   CMC4    2  2PLM 0.9627020 0.11608725 -0.42559338 0.11102431         NA
#> 5   CMC5    2  1PLM 1.0000000         NA -0.25414016 0.09928963         NA
#> 6   CMC6    2  3PLM 1.8846619 0.26669060  0.67350821 0.09232536  0.1041510
#> 7   CMC7    2  3PLM 0.8866580 0.17171857  1.03619252 0.22999854  0.1374716
#> 8   CMC8    2  2PLM 0.9232153 0.11779263  0.86515210 0.13411409         NA
#> 9   CMC9    2  2PLM 1.0043493 0.12280227  0.89227667 0.12623451         NA
#> 10 CMC10    2  2PLM 1.6128253 0.15351879  0.09253022 0.06743190         NA
#> 11 CMC11    2  2PLM 1.0684709 0.12147067 -0.37379104 0.09986851         NA
#> 12 CMC12    2  2PLM 0.9377669 0.12170130  1.09646572 0.14960502         NA
#> 13 CMC13    2  3PLM 1.3645089 0.33830056  1.31503341 0.16875186  0.1703375
#> 14 CMC14    2  3PLM 1.3898925 0.32184562  0.17217473 0.23754432  0.2472123
#> 15 CMC15    2  3PLM 1.5408435 0.27056505  0.01937278 0.17343671  0.2070791
#> 16 CMC16    2  3PLM 2.1026954 0.24879579  0.04186170 0.08406154  0.1039632
#> 17 CMC17    2  3PLM 1.0234932 0.15550550 -0.39553054 0.22209049  0.1653727
#> 18 CMC18    2  3PLM 1.2856593 0.38393374  1.43007415 0.20356474  0.2250732
#> 19 CMC19    2  3PLM 2.2594756 0.32044597 -1.10181795 0.14051114  0.1778769
#> 20 CMC20    2  3PLM 1.4700031 0.21920822 -1.72399703 0.22689431  0.1889865
#> 21 CMC21    2  3PLM 1.3912594 0.21215536 -1.23367452 0.22895293  0.2093723
#> 22 CMC22    2  3PLM 0.9317523 0.15932023 -0.52661586 0.28436101  0.1982627
#> 23 CMC23    2  3PLM 1.1204103 0.22727469 -0.09122518 0.27935617  0.2259506
#> 24 CMC24    2  3PLM 1.2290146 0.33771111  1.43876383 0.21032077  0.2204891
#> 25 CMC25    2  3PLM 0.8330125 0.16035815 -1.46925302 0.41455949  0.2244715
#> 26 CMC26    2  3PLM 1.0676986 0.18549120 -2.13985943 0.35014363  0.1999090
#> 27 CMC27    2  3PLM 1.1843144 0.17838470  0.09628328 0.17043885  0.1469717
#> 28 CMC28    2  3PLM 2.2044237 0.31674844 -0.16348433 0.10953967  0.1944618
#> 29 CMC29    2  3PLM 2.5307240 0.54619578 -0.78717016 0.19184713  0.3956548
#> 30 CMC30    2  3PLM 1.9041795 0.44751927  0.70161781 0.15026890  0.3404332
#> 31 CMC31    2  3PLM 0.7062845 0.16316408  1.02018205 0.32764850  0.1641299
#> 32 CMC32    2  3PLM 1.7564751 0.31082457 -0.76184344 0.21160840  0.2681714
#> 33 CMC33    2  3PLM 1.0744198 0.17077250 -1.42868292 0.28271978  0.1991516
#> 34 CMC34    2  3PLM 1.0509163 0.17161904  0.22205598 0.19700301  0.1609429
#> 35 CMC35    2  3PLM 1.3654004 0.19069948 -0.43712558 0.16865435  0.1629476
#> 36 CMC36    2  3PLM 0.8871087 0.17401848  0.97700579 0.23224001  0.1423811
#> 37 CMC37    2  3PLM 2.1311743 0.26144273 -0.24278137 0.09498708  0.1309809
#> 38 CMC38    2  3PLM 0.8813430 0.17350226 -0.28321005 0.32455119  0.2131754
#> 39  CFR1    5   GRM 2.0000867 0.14447773 -1.88248141 0.11978029 -1.2545794
#> 40  CFR2    5   GRM 1.3885998 0.10798168 -0.79658308 0.09234859 -0.1288726
#> 41  AMC1    2  3PLM 1.8272686 0.39127793  0.74245570 0.14263810  0.2788644
#> 42  AMC2    2  3PLM 1.7110855 0.25426145 -1.58183908 0.20042801  0.1947024
#> 43  AMC3    2  3PLM 1.3145385 0.25146602  0.68535647 0.16446443  0.1615132
#> 44  AMC4    2  3PLM 0.9496703 0.17005564 -0.15473831 0.26323511  0.1915267
#> 45  AMC5    2  3PLM 1.7029567 0.64633195  2.10717153 0.26189935  0.1910672
#> 46  AMC6    2  3PLM 2.8374610 0.63744855  1.44493420 0.09747301  0.1542768
#> 47  AMC7    2  3PLM 1.7188353 0.41094702  0.38574744 0.17861264  0.2586265
#> 48  AMC8    2  3PLM 1.6576325 0.29286642  0.39605794 0.13785201  0.2033722
#> 49  AMC9    2  3PLM 1.5568580 0.26091218  0.49337767 0.13241328  0.1543679
#> 50 AMC10    2  3PLM 2.4905705 0.50988663  1.31162574 0.09509815  0.1294143
#> 51 AMC11    2  3PLM 1.7412715 0.22999194 -1.01102653 0.15624576  0.1696460
#> 52 AMC12    2  3PLM 0.9688530 0.20770051 -0.77937288 0.39869458  0.2552984
#> 53  AFR1    5   GRM 1.1355459 0.10188574 -0.30016570 0.09221751  0.3028394
#> 54  AFR2    5  GPCM 1.3267609 0.10884827 -1.99444989 0.21458427 -1.3127161
#> 55  AFR3    5  GPCM 0.8941751 0.07361344 -0.79542779 0.15373354  0.1540945
#>          se.3      par.4       se.4      par.5       se.5
#> 1          NA         NA         NA         NA         NA
#> 2          NA         NA         NA         NA         NA
#> 3          NA         NA         NA         NA         NA
#> 4          NA         NA         NA         NA         NA
#> 5          NA         NA         NA         NA         NA
#> 6  0.03138242         NA         NA         NA         NA
#> 7  0.05092018         NA         NA         NA         NA
#> 8          NA         NA         NA         NA         NA
#> 9          NA         NA         NA         NA         NA
#> 10         NA         NA         NA         NA         NA
#> 11         NA         NA         NA         NA         NA
#> 12         NA         NA         NA         NA         NA
#> 13 0.04452707         NA         NA         NA         NA
#> 14 0.08294538         NA         NA         NA         NA
#> 15 0.06714419         NA         NA         NA         NA
#> 16 0.03704144         NA         NA         NA         NA
#> 17 0.07041155         NA         NA         NA         NA
#> 18 0.05087416         NA         NA         NA         NA
#> 19 0.07448965         NA         NA         NA         NA
#> 20 0.08281388         NA         NA         NA         NA
#> 21 0.08608492         NA         NA         NA         NA
#> 22 0.08224986         NA         NA         NA         NA
#> 23 0.08847770         NA         NA         NA         NA
#> 24 0.04739131         NA         NA         NA         NA
#> 25 0.09844718         NA         NA         NA         NA
#> 26 0.08835993         NA         NA         NA         NA
#> 27 0.05738007         NA         NA         NA         NA
#> 28 0.05270285         NA         NA         NA         NA
#> 29 0.08610236         NA         NA         NA         NA
#> 30 0.04908313         NA         NA         NA         NA
#> 31 0.06684690         NA         NA         NA         NA
#> 32 0.08879304         NA         NA         NA         NA
#> 33 0.08601695         NA         NA         NA         NA
#> 34 0.06017396         NA         NA         NA         NA
#> 35 0.06511108         NA         NA         NA         NA
#> 36 0.05368437         NA         NA         NA         NA
#> 37 0.04587033         NA         NA         NA         NA
#> 38 0.08786132         NA         NA         NA         NA
#> 39 0.08307182 -0.7031275 0.06307867 -0.2324418 0.05703100
#> 40 0.07429798  0.6012371 0.08110680  1.0859619 0.09995812
#> 41 0.04793931         NA         NA         NA         NA
#> 42 0.08320344         NA         NA         NA         NA
#> 43 0.05423568         NA         NA         NA         NA
#> 44 0.07656965         NA         NA         NA         NA
#> 45 0.02992857         NA         NA         NA         NA
#> 46 0.02219971         NA         NA         NA         NA
#> 47 0.06713238         NA         NA         NA         NA
#> 48 0.05193618         NA         NA         NA         NA
#> 49 0.04845912         NA         NA         NA         NA
#> 50 0.02401244         NA         NA         NA         NA
#> 51 0.07189865         NA         NA         NA         NA
#> 52 0.10907235         NA         NA         NA         NA
#> 53 0.08931676  0.9182795 0.10902552  1.3525597 0.13282858
#> 54 0.15227381 -0.7196856 0.11752786 -0.2068705 0.10041751
#> 55 0.15429938  0.4638485 0.16257565  1.3511384 0.18586806

# }
```
