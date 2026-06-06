# Multiple-group item calibration using MMLE-EM algorithm

This function performs multiple-group item calibration (Bock & Zimowski,
1997) using marginal maximum likelihood estimation via the
expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin, 1981). It
also supports multiple-group fixed item parameter calibration (MG-FIPC;
e.g., Kim & Kolen, 2016), which extends the single-group FIPC method
(Kim, 2006) to multiple-group settings. For dichotomous items, the
function supports one-, two-, and three-parameter logistic IRT models.
For polytomous items, the graded response model (GRM) and the
(generalized) partial credit model (GPCM) are available.

## Usage

``` r
est_mg(
  x = NULL,
  data,
  group.name = NULL,
  D = 1,
  model = NULL,
  cats = NULL,
  item.id = NULL,
  free.group = NULL,
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
  Quadrature = c(49, 6),
  weights = NULL,
  group.mean = 0,
  group.var = 1,
  EmpHist = FALSE,
  use.startval = FALSE,
  Etol = 0.001,
  MaxE = 500,
  control = list(eval.max = 500, iter.max = 200, x.tol = 1e-04),
  fipc = FALSE,
  fipc.method = "MEM",
  fix.loc = NULL,
  fix.id = NULL,
  se = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  A list containing item metadata for all groups to be analyzed. For
  example, if five groups are analyzed, the list should contain five
  elements, each representing the item metadata for one group. The order
  of the elements in the list must match the order of group names
  specified in the `group.name` argument.

  Each group's item metadata includes essential information for each
  item (e.g., number of score categories, IRT model type, etc.) required
  for calibration. See
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata.

  When `use.startval = TRUE`, the item parameters specified in the
  metadata will be used as starting values for parameter estimation. If
  `x = NULL`, both `model` and `cats` arguments must be specified. Note
  that when `fipc = TRUE` to implement MG-FIPC, the `x` argument must be
  specified and cannot be NULL. Default is `NULL`.

- data:

  A list containing item response matrices for all groups to be
  analyzed. For example, if five groups are analyzed, the list should
  include five elements, each representing the response data matrix for
  one group. The elements in the list must be ordered to match the group
  names specified in the `group.name` argument. Each matrix contains
  examinees' item responses corresponding to the item metadata for that
  group. In each matrix, rows represent examinees and columns represent
  items.

- group.name:

  A character vector indicating the names of the groups. For example, if
  five groups are analyzed, use
  `group.name = c("G1", "G2", "G3", "G4", "G5")`. Group names can be any
  valid character strings.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- model:

  A list containing character vectors specifying the IRT models used to
  calibrate items across all groups. For example, if five groups are
  analyzed, the list should contain five elements, each being a
  character vector of IRT model names for one group. The elements in the
  list must be ordered to match the group names specified in the
  `group.name` argument.

  Available IRT models include:

  - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items

  - `"GRM"`, `"GPCM"` for polytomous items

  Here, `"GRM"` denotes the graded response model and `"GPCM"` the
  (generalized) partial credit model. Note that `"DRM"` serves as a
  general label covering all three dichotomous IRT models.If a single
  model name is provided in any element of the list, it will be recycled
  across all items within that group.

  This argument is used only when `x = NULL` and `fipc = FALSE`. Default
  is `NULL`.

- cats:

  A list containing numeric vectors specifying the number of score
  categories for items in each group. For example, if five groups are
  analyzed, the list should contain five numeric vectors corresponding
  to the five groups. The elements in the list must be ordered
  consistently with the group names specified in the `group.name`
  argument.

  If a single numeric value is specified in any element of the list, it
  will be recycled across all items in the corresponding group. If
  `cats = NULL` and all models specified in the `model` argument are
  dichotomous (i.e., "1PLM", "2PLM", "3PLM", or "DRM"), the function
  assumes that all items have two score categories across all groups.

  This argument is used only when `x = NULL` and `fipc = FALSE`. Default
  is `NULL`.

- item.id:

  A list containing character vectors of item IDs for each group to be
  analyzed. For example, if five groups are analyzed, the list should
  contain five character vectors corresponding to the five groups. The
  elements in the list must be ordered consistently with the group names
  specified in the `group.name` argument.

  When `fipc = TRUE` and item IDs are provided via the `item.id`
  argument, the item IDs in the `x` argument will be overridden. Default
  is `NULL`.

- free.group:

  A numeric or character vector indicating the groups for which the
  scales (i.e., mean and standard deviation) of the latent ability
  distributions are freely estimated. The scales of the remaining groups
  (those not specified in this argument) are fixed using the values
  provided in the `group.mean` and `group.var` arguments, or from the
  `weights` argument.

  For example, suppose that five groups are analyzed with group names
  "G1", "G2", "G3", "G4", and "G5". To freely estimate the scales for
  groups 2 through 5, set `free.group = c(2, 3, 4, 5)` or
  `free.group = c("G2", "G3", "G4", "G5")`. In this case, the first
  group ("G1") will have a fixed scale (e.g., a mean of 0 and variance
  of 1 when `group.mean = 0`, `group.var = 1`, and `weights = NULL`).

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

- Quadrature:

  A numeric vector of length two:

  - first element: number of quadrature points

  - second element: symmetric bound (absolute value) for those points
    For example, `c(49, 6)` specifies 49 evenly spaced points from –6
    to 6. These points are used in the E-step of the EM algorithm.
    Default is `c(49, 6)`.

- weights:

  A two-column matrix or data frame containing the quadrature points (in
  the first column) and the corresponding weights (in the second column)
  for the latent ability prior distribution. If not `NULL`, the latent
  ability distributions for the groups not specified in the `free.group`
  argument are fixed to match the scale defined by the provided
  quadrature points and weights. The weights and points can be
  conveniently generated using the function
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).

  If `NULL`, a normal prior density is used instead, based on the
  information provided in the `Quadrature`, `group.mean`, and
  `group.var` arguments. Default is `NULL`.

- group.mean:

  A numeric value specifying the mean of the latent variable prior
  distribution when `weights = NULL`. Default is 0. For groups not
  specified in the `free.group` argument, their distribution means are
  fixed to this value in order to resolve the indeterminacy of the item
  parameter scale.

- group.var:

  A positive numeric value specifying the variance of the latent
  variable prior distribution when `weights = NULL`. Default is 1. For
  groups not specified in the `free.group` argument, their distribution
  variances are fixed to this value in order to resolve the
  indeterminacy of the item parameter scale.

- EmpHist:

  Logical. If `TRUE`, the empirical histograms of the latent ability
  prior distributions across all groups are estimated simultaneously
  with the item parameters using the approach proposed by Woods (2007).
  Item calibration is then performed relative to the estimated empirical
  priors.

- use.startval:

  Logical. If `TRUE`, the item parameters provided in the item metadata
  (i.e., the `x` argument) are used as starting values for item
  parameter estimation. Otherwise, internally generated starting values
  are used. Default is `FALSE`.

- Etol:

  A positive numeric value specifying the convergence criterion for the
  E-step of the EM algorithm. Default is 1e-3. Specifically, the EM
  algorithm terminates when the largest absolute difference in item
  parameter estimates between consecutive iterations is smaller than
  this value.

- MaxE:

  A positive integer specifying the maximum number of iterations for the
  E-step in the EM algorithm. Default is `500`.

- control:

  A named list of options passed directly to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) in each
  M‑step optimization of the EM algorithm. By default:
  `control = list(eval.max = 500, iter.max = 200, x.tol = 1e-4)`, where

  - `eval.max` = 500 limits the number of function evaluations

  - `iter.max` = 200 caps the number of internal optimizer iterations

  - `x.tol` = 1e‑4 sets the absolute change threshold in parameter
    values below which
    [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) considers
    the solution to have converged Users may additionally supply other
    [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) control options
    (such as `abs.tol`, `rel.tol`, `trace`, etc.) as needed.

- fipc:

  Logical. If `TRUE`, multiple-group fixed item parameter calibration
  (MG-FIPC) is applied during item parameter estimation. When
  `fipc = TRUE`, the information on which items are fixed must be
  provided via either `fix.loc` or `fix.id`. See below for details.

- fipc.method:

  A character string specifying the FIPC method. Available options are:

  - `"OEM"`: No Prior Weights Updating and One EM Cycle (NWU-OEM; Wainer
    & Mislevy, 1990)

  - `"MEM"`: Multiple Prior Weights Updating and Multiple EM Cycles
    (MWU-MEM; Kim, 2006) When `fipc.method = "OEM"`, the maximum number
    of E-steps is automatically set to 1, regardless of the value
    specified in `MaxE`.

- fix.loc:

  A list of positive integer vectors. Each internal vector specifies the
  positions of the items to be fixed in the item metadata (i.e., `x`)
  for each group when MG-FIPC is implemented (i.e., `fipc = TRUE`). The
  internal objects in the list must follow the same order as the group
  names provided in the `group.name` argument.

  For example, suppose three groups are analyzed. In the first group,
  the 1st, 3rd, and 5th items are fixed; in the second group, the 2nd,
  3rd, 4th, and 7th items are fixed; and in the third group, the 1st,
  2nd, and 6th items are fixed. Then
  `fix.loc = list(c(1, 3, 5), c(2, 3, 4, 7), c(1, 2, 6))`. Note that if
  the `fix.id` argument is not NULL, the information in `fix.loc` will
  be ignored. See below for details.

- fix.id:

  A vector of character strings specifying the IDs of items to be fixed
  when MG-FIPC is implemented (i.e., `fipc = TRUE`).

  For example, suppose that three groups are analyzed. In the first
  group, three items with IDs G1I1, C1I1, and C1I2 are fixed. In the
  second group, four items with IDs C1I1, C1I2, C2I1, and C2I2 are
  fixed. In the third group, three items with IDs C2I1, C2I2, and G3I1
  are fixed.

  In this case, there are six unique items fixed across the
  groups—namely, G1I1, C1I1, C1I2, C2I1, C2I2, and G3I1, because C1I1
  and C1I2 appear in both the first and second groups, while C2I1 and
  C2I2 appear in both the second and third groups. Thus, you should
  specify `fix.id = c("G1I1", "C1I1", "C1I2", "C2I1", "C2I2", "G3I1")`.
  Note that if the `fix.id` argument is not NULL, the information
  provided in `fix.loc` is ignored. See below for details.

- se:

  Logical. If `FALSE`, standard errors of the item parameter estimates
  are not computed. Default is `TRUE`.

- verbose:

  Logical. If `FALSE`, all progress messages, including information
  about the EM algorithm process, are suppressed. Default is `TRUE`.

## Value

This function returns an object of class `est_irt`. The returned object
contains the following components:

- estimates:

  A list containing two internal elements: `overall` and `group`. The
  `overall` element is a data frame with item parameter estimates and
  their standard errors, computed from the combined data across all
  groups. This data frame includes only the unique items across all
  groups. The `group` element is a list of group-specific data frames,
  each containing item parameter estimates and standard errors for that
  particular group.

- par.est:

  A list with the same structure as `estimates`, containing only the
  item parameter estimates (excluding standard errors), formatted
  according to the item metadata structure.

- se.est:

  A list with the same structure as `estimates`, but containing only the
  standard errors of the item parameter estimates. Note that the
  standard errors are calculated using the cross-product approximation
  method (Meilijson, 1989).

- pos.par:

  A data frame indicating the position index of each estimated item
  parameter. This index is based on the combined data set across all
  groups (i.e., the first internal object of `estimates`). The position
  information is useful for interpreting the variance-covariance matrix
  of item parameter estimates.

- covariance:

  A variance-covariance matrix of the item parameter estimates, based on
  the combined data set across all groups (i.e., the first internal
  object of `estimates`).

- loglikelihood:

  A list containing two internal objects (i.e., overall and group) of
  marginal log-likelihood values based on the observed data. The
  structure of the list matches that of `estimates`. Specifically, the
  `overall` component contains the total log-likelihood summed across
  all unique items from all groups, while the `group` component provides
  group-specific log-likelihood values.

- aic:

  A model fit statistic based on the Akaike Information Criterion (AIC),
  calculated from the log-likelihood of all unique items.

- bic:

  A model fit statistic based on the Bayesian Information Criterion
  (BIC), calculated from the log-likelihood of all unique items.

- group.par:

  A list containing summary statistics (i.e., mean, variance, and
  standard deviation) of the latent variable prior distributions across
  all groups.

- weights:

  A list of two-column data frames, where the first column contains
  quadrature points and the second column contains the corresponding
  weights of the (updated) latent variable prior distributions for each
  group.

- posterior.dist:

  A matrix of normalized posterior densities for all response patterns
  at each quadrature point. Rows and columns represent response patterns
  and quadrature points, respectively.

- data:

  A list containing two internal objects (i.e., overall and group)
  representing the examinees' response data sets. The structure of this
  list matches that of the `estimates` component.

- scale.D:

  The scaling factor used in the IRT model.

- ncase:

  A list containing two internal objects (i.e., overall and group)
  representing the total number of response patterns. The structure of
  this list matches that of the `estimates` component.

- nitem:

  A list containing two internal objects (i.e., overall and group)
  representing the total number of items included in the response data.
  The structure of this list matches that of the `estimates` component.

- Etol:

  The convergence criterion for the E-step of the EM algorithm.

- MaxE:

  The maximum number of E-steps allowed in the EM algorithm.

- aprior:

  A list describing the prior distribution used for discrimination
  parameters.

- bprior:

  A list describing the prior distribution used for difficulty
  parameters.

- gprior:

  A list describing the prior distribution used for guessing parameters.

- npar.est:

  The total number of parameters estimated across all unique items.

- niter:

  The number of completed EM cycles.

- maxpar.diff:

  The maximum absolute change in parameter estimates at convergence.

- EMtime:

  Time (in seconds) spent on EM cycles.

- SEtime:

  Time (in seconds) spent computing standard errors.

- TotalTime:

  Total computation time (in seconds).

- test.1:

  First-order test result indicating whether the gradient sufficiently
  vanished for solution stability.

- test.2:

  Second-order test result indicating whether the information matrix is
  positive definite, a necessary condition for identifying a local
  maximum.

- var.note:

  A note indicating whether the variance-covariance matrix was
  successfully obtained from the information matrix.

- fipc:

  Logical. Indicates whether FIPC was used.

- fipc.method:

  The method used for FIPC.

- fix.loc:

  A list containing two internal objects (i.e., overall and group)
  indicating the locations of fixed items when FIPC is applied. The
  structure of the list matches that of the 'estimates' component.

Note that you can easily extract components from the output using the
[`getirt()`](https://hwangQ.github.io/irtQ/reference/getirt.md)
function.

## Details

Multiple-group (MG) item calibration (Bock & Zimowski, 1996) provides a
unified framework for handling testing scenarios involving multiple
groups, such as nonequivalent groups equating, vertical scaling, and the
identification of differential item functioning (DIF). In such
applications, examinees from different groups typically respond to
either the same test form or to different forms that share common
(anchor) items.

The goal of MG item calibration is to estimate both item parameters and
latent ability distributions for all groups simultaneously (Bock &
Zimowski, 1996). The irtQ package implements MG calibration via the
`est_mg()` function, which uses marginal maximum likelihood estimation
through the expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin,
1981). In addition, the function supports multiple-group fixed item
parameter calibration (MG-FIPC; e.g., Kim & Kolen, 2016), which allows
the parameters of specific items to be fixed across groups.

In MG IRT analyses, it is common for multiple groups' test forms to
share some common (anchor) items. By default, the `est_mg()` function
automatically constrains items with identical item IDs across groups to
share the same parameter estimates.

Most of the features of the `est_mg()` function are similar to those of
the [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
function. The main difference is that several arguments in `est_mg()`
accept list objects containing elements for each group to be analyzed.
These arguments include `x`, `data`, `model`, `cats`, `item.id`,
`fix.loc` and `fix.id`.

Additionally, `est_mg()` introduces two new arguments: `group.name` and
`free.group`. The `group.name` argument is required to assign a unique
identifier to each group. The order of the list elements provided in
`x`, `data`, `model`, `cats`, `item.id`, `fix.loc` and `fix.id` must
match the order of group names specified in the `group.name` argument.

The `free.group` argument is required to indicate which groups have
their latent ability distribution scales (i.e., mean and standard
deviation) freely estimated. When no item parameters are fixed (i.e.,
`fipc = FALSE`), at least one group must have a fixed latent ability
scale (e.g., mean = 0 and variance = 1) among the multiple groups
sharing common items, in order to resolve the scale indeterminacy
inherent in IRT estimation. By specifying the groups in the `free.group`
argument, the scales for those groups will be freely estimated, while
the scales for all other groups not included in `free.group` will be
fixed using the values provided in the `group.mean` and `group.var`
arguments or from the `weights` argument.

Situations requiring the implementation of MG-FIPC typically arise when
new latent ability scales from multiple-group (MG) test data need to be
linked to an established scale (e.g., that of an existing item bank). In
a single run of the MG-FIPC procedure, the parameters of non-fixed
(freed) items across multiple test forms, as well as the latent ability
distributions for multiple groups, can be estimated on the same scale as
the fixed items (Kim & Kolen, 2016).

For example, suppose that three different test forms—Form 1, Form 2, and
Form 3—are administered to three nonequivalent groups: Group1, Group2,
and Group3. Form 1 and Form 2 share 12 common items (C1I1 to C1I12),
while Form 2 and Form 3 share 10 common items (C2I1 to C2I10). There are
no common items between Form 1 and Form 3. Also, assume that all unique
items in Form 1 are from an existing item bank and have already been
calibrated on the item bank's scale.

In this case, the goal of MG-FIPC is to estimate the parameters of all
items across the three test forms—except the unique items in Form 1— and
the latent ability distributions of the three groups, all on the same
scale as the item bank. To achieve this, the unique items in Form 1 must
be fixed during MG-FIPC to link the current MG test data to the item
bank scale.

The `est_mg()` function can implement MG-FIPC by setting `fipc = TRUE`.
In this case, the information on which items to fix must be provided
through either the `fix.loc` or `fix.id` argument. When using `fix.loc`,
you must supply a list of item positions (locations) to be fixed in each
group’s test form. For example, suppose that the test data from the
three groups above are analyzed. In the first group, the 1st, 3rd, and
5th items are fixed; in the second group, the 2nd, 3rd, 4th, and 7th
items are fixed; and in the third group, the 1st, 2nd, and 6th items are
fixed. In this case, you would specify:
`fix.loc = list(c(1, 3, 5), c(2, 3, 4, 7), c(1, 2, 6))`.

Alternatively, you can use `fix.id` to specify a character vector of
item IDs to be fixed across groups. In the first group, the items with
IDs G1I1, C1I1, and C1I2 are fixed; in the second group, the items with
IDs C1I1, C1I2, C2I1, and C2I2 are fixed; and in the third group, the
items with IDs C2I1, C2I2, and G3I1 are fixed. In this case, there are
six unique items to be fixed across all groups: G1I1, C1I1, C1I2, C2I1,
C2I2, and G3I1. You would then specify:
`fix.id = c("G1I1", "C1I1", "C1I2", "C2I1", "C2I2", "G3I1")`.

Note that when both `fix.loc` and `fix.id` are provided, the information
in `fix.id` takes precedence and overrides `fix.loc`.

## References

Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation
of item parameters: Application of an EM algorithm. *Psychometrika, 46*,
443-459.

Bock, R. D., & Zimowski, M. F. (1997). Multiple group IRT. In W. J. van
der Linden & R. K. Hambleton (Eds.), *Handbook of modern item response
theory* (pp. 433-448). New York: Springer.

Kim, S. (2006). A comparative study of IRT fixed parameter calibration
methods. *Journal of Educational Measurement, 43*(4), 355-381.

Kim, S., & Kolen, M. J. (2016). Multiple group IRT fixed-parameter
estimation for maintaining an established ability scale. *Center for
Advanced Studies in Measurement and Assessment Report, 49.*

Meilijson, I. (1989). A fast improvement to the EM algorithm on its own
terms. *Journal of the Royal Statistical Society: Series B
(Methodological), 51*, 127-138.

Woods, C. M. (2007). Empirical histograms in item response theory with
ordinal data. *Educational and Psychological Measurement, 67*(1), 73-87.

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`shape_df_fipc()`](https://hwangQ.github.io/irtQ/reference/shape_df_fipc.md),
[`getirt()`](https://hwangQ.github.io/irtQ/reference/getirt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
## ------------------------------------------------------------------------------
# 1. MG calibration using the simMG data
#  - Details:
#    (a) Constrain common items between groups to have
#        identical item parameters (i.e., items C1I1–C1I12 between
#        Groups 1 and 2, and items C2I1–C2I10 between Groups 2 and 3).
#    (b) Freely estimate the means and variances of the ability
#        distributions for all groups except the reference group,
#        where the mean and variance are fixed to 0 and 1, respectively.
## ------------------------------------------------------------------------------
# 1-(1). Freely estimate the means and variances of Groups 2 and 3
# Import the true item metadata for the three groups
x <- simMG$item.prm

# Extract model, score category, and item ID information
# from the item metadata for the three groups
model <- list(x$Group1$model, x$Group2$model, x$Group3$model)
cats <- list(x$Group1$cats, x$Group2$cats, x$Group3$cats)
item.id <- list(x$Group1$id, x$Group2$id, x$Group3$id)

# Import the simulated response data sets for the three groups
data <- simMG$res.dat

# Import the group names for the three groups
group.name <- simMG$group.name

# Specify Groups 2 and 3 as the free groups where the scale
# of the ability distributions will be freely estimated.
# Group 1 will serve as the reference group, where the scale
# of the ability distribution is fixed to the values specified
# via the 'group.mean' and 'group.var' arguments
free.group <- c(2, 3) # or use 'free.group <- group.name[2:3]'

# Estimate IRT parameters:
# As long as common items across groups share the same item IDs,
# their item parameters will be constrained to be equal across groups
# unless FIPC is implemented.
fit.1 <-
  est_mg(
    data = data, group.name = group.name, model = model,
    cats = cats, item.id = item.id, D = 1, free.group = free.group,
    use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16)),
    group.mean = 0, group.var = 1, EmpHist = TRUE, Etol = 0.001, MaxE = 500
  )
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -184619.3757, Max-Change: 2.188782 EM iteration: 2, Loglike: -159724.2438, Max-Change: 0.42374 EM iteration: 3, Loglike: -159631.0877, Max-Change: 0.166346 EM iteration: 4, Loglike: -159617.3743, Max-Change: 0.083614 EM iteration: 5, Loglike: -159611.4702, Max-Change: 0.052115 EM iteration: 6, Loglike: -159608.0480, Max-Change: 0.038966 EM iteration: 7, Loglike: -159605.7836, Max-Change: 0.032285 EM iteration: 8, Loglike: -159604.1351, Max-Change: 0.027248 EM iteration: 9, Loglike: -159602.8431, Max-Change: 0.023166 EM iteration: 10, Loglike: -159601.7733, Max-Change: 0.01977 EM iteration: 11, Loglike: -159600.8519, Max-Change: 0.016922 EM iteration: 12, Loglike: -159600.0358, Max-Change: 0.014527 EM iteration: 13, Loglike: -159599.2989, Max-Change: 0.012509 EM iteration: 14, Loglike: -159598.6241, Max-Change: 0.010869 EM iteration: 15, Loglike: -159597.9998, Max-Change: 0.009797 EM iteration: 16, Loglike: -159597.4177, Max-Change: 0.00890 EM iteration: 17, Loglike: -159596.8716, Max-Change: 0.008145 EM iteration: 18, Loglike: -159596.3569, Max-Change: 0.007856 EM iteration: 19, Loglike: -159595.8700, Max-Change: 0.007658 EM iteration: 20, Loglike: -159595.4080, Max-Change: 0.00744 EM iteration: 21, Loglike: -159594.9687, Max-Change: 0.007208 EM iteration: 22, Loglike: -159594.5501, Max-Change: 0.006968 EM iteration: 23, Loglike: -159594.1507, Max-Change: 0.006723 EM iteration: 24, Loglike: -159593.7692, Max-Change: 0.006477 EM iteration: 25, Loglike: -159593.4043, Max-Change: 0.006232 EM iteration: 26, Loglike: -159593.0551, Max-Change: 0.00599 EM iteration: 27, Loglike: -159592.7206, Max-Change: 0.005753 EM iteration: 28, Loglike: -159592.4000, Max-Change: 0.005522 EM iteration: 29, Loglike: -159592.0925, Max-Change: 0.005298 EM iteration: 30, Loglike: -159591.7975, Max-Change: 0.005081 EM iteration: 31, Loglike: -159591.5143, Max-Change: 0.004871 EM iteration: 32, Loglike: -159591.2423, Max-Change: 0.00467 EM iteration: 33, Loglike: -159590.9808, Max-Change: 0.004476 EM iteration: 34, Loglike: -159590.7295, Max-Change: 0.00429 EM iteration: 35, Loglike: -159590.4879, Max-Change: 0.004112 EM iteration: 36, Loglike: -159590.2553, Max-Change: 0.003942 EM iteration: 37, Loglike: -159590.0316, Max-Change: 0.003779 EM iteration: 38, Loglike: -159589.8161, Max-Change: 0.003624 EM iteration: 39, Loglike: -159589.6086, Max-Change: 0.003475 EM iteration: 40, Loglike: -159589.4087, Max-Change: 0.003333 EM iteration: 41, Loglike: -159589.2160, Max-Change: 0.003198 EM iteration: 42, Loglike: -159589.0303, Max-Change: 0.003068 EM iteration: 43, Loglike: -159588.8512, Max-Change: 0.002945 EM iteration: 44, Loglike: -159588.6785, Max-Change: 0.002827 EM iteration: 45, Loglike: -159588.5119, Max-Change: 0.002715 EM iteration: 46, Loglike: -159588.3511, Max-Change: 0.002607 EM iteration: 47, Loglike: -159588.1959, Max-Change: 0.002505 EM iteration: 48, Loglike: -159588.0460, Max-Change: 0.002407 EM iteration: 49, Loglike: -159587.9013, Max-Change: 0.002313 EM iteration: 50, Loglike: -159587.7616, Max-Change: 0.002228 EM iteration: 51, Loglike: -159587.6266, Max-Change: 0.002156 EM iteration: 52, Loglike: -159587.4961, Max-Change: 0.002086 EM iteration: 53, Loglike: -159587.3700, Max-Change: 0.002019 EM iteration: 54, Loglike: -159587.2481, Max-Change: 0.001955 EM iteration: 55, Loglike: -159587.1302, Max-Change: 0.001893 EM iteration: 56, Loglike: -159587.0163, Max-Change: 0.001833 EM iteration: 57, Loglike: -159586.9060, Max-Change: 0.001776 EM iteration: 58, Loglike: -159586.7993, Max-Change: 0.001721 EM iteration: 59, Loglike: -159586.6961, Max-Change: 0.001667 EM iteration: 60, Loglike: -159586.5962, Max-Change: 0.001616 EM iteration: 61, Loglike: -159586.4995, Max-Change: 0.001566 EM iteration: 62, Loglike: -159586.4058, Max-Change: 0.001519 EM iteration: 63, Loglike: -159586.3151, Max-Change: 0.001473 EM iteration: 64, Loglike: -159586.2272, Max-Change: 0.001429 EM iteration: 65, Loglike: -159586.1421, Max-Change: 0.001386 EM iteration: 66, Loglike: -159586.0596, Max-Change: 0.001358 EM iteration: 67, Loglike: -159585.9796, Max-Change: 0.001334 EM iteration: 68, Loglike: -159585.9021, Max-Change: 0.00131 EM iteration: 69, Loglike: -159585.8269, Max-Change: 0.001287 EM iteration: 70, Loglike: -159585.7540, Max-Change: 0.001264 EM iteration: 71, Loglike: -159585.6832, Max-Change: 0.001241 EM iteration: 72, Loglike: -159585.6146, Max-Change: 0.001219 EM iteration: 73, Loglike: -159585.5479, Max-Change: 0.001197 EM iteration: 74, Loglike: -159585.4832, Max-Change: 0.001175 EM iteration: 75, Loglike: -159585.4204, Max-Change: 0.001153 EM iteration: 76, Loglike: -159585.3593, Max-Change: 0.001132 EM iteration: 77, Loglike: -159585.3000, Max-Change: 0.001111 EM iteration: 78, Loglike: -159585.2424, Max-Change: 0.001091 EM iteration: 79, Loglike: -159585.1864, Max-Change: 0.00107 EM iteration: 80, Loglike: -159585.1319, Max-Change: 0.00105 EM iteration: 81, Loglike: -159585.0790, Max-Change: 0.001031 EM iteration: 82, Loglike: -159585.0274, Max-Change: 0.001011 EM iteration: 83, Loglike: -159584.9773, Max-Change: 0.000992 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 16.07 seconds. 

# Summary of the estimation
summary(fit.1)
#> 
#> Call:
#> est_mg(data = data, group.name = group.name, D = 1, model = model, 
#>     cats = cats, item.id = item.id, free.group = free.group, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), group.mean = 0, group.var = 1, EmpHist = TRUE, 
#>     Etol = 0.001, MaxE = 500)
#> 
#> Summary of the Data 
#>  Number of Items: 
#>   Overall: 116 unique items 
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of Cases: 
#>   Overall: 6000
#>   By group: 2000(Group1), 2000(Group2), 2000(Group3)
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 362
#>  Number of fixed items: 
#>   Overall: 0
#>   By group: 0(Group1), 0(Group2), 0(Group3)
#>  Number of E-step cycles completed: 83
#>  Maximum parameter change: 0.0009920582
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 15.39
#>  Standard error computation: 0.2
#>  Total computation: 16.07
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 
#>   Overall: 319169.9
#>   By group: 120347.958(Group1), 113944.899(Group2), 84877.056(Group3)
#> 
#>  Akaike Information Criterion (AIC): 319893.9
#>  Bayesian Information Criterion (BIC): 322319.1
#>  Item Parameters (Overall): 
#>         id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1     C1I1     2   3PLM   0.90  0.18   1.36  0.16   0.28  0.05     NA    NA
#> 2     C1I2     2   3PLM   2.12  0.15  -0.99  0.10   0.18  0.06     NA    NA
#> 3     C1I3     2   3PLM   1.05  0.12   0.61  0.13   0.18  0.05     NA    NA
#> 4     C1I4     2   3PLM   1.08  0.13  -0.17  0.22   0.28  0.07     NA    NA
#> 5     C1I5     2   3PLM   0.87  0.09  -0.15  0.22   0.18  0.07     NA    NA
#> 6     C1I6     2   3PLM   1.91  0.13   0.60  0.04   0.09  0.02     NA    NA
#> 7     C1I7     2   3PLM   1.10  0.13   1.10  0.09   0.15  0.03     NA    NA
#> 8     C1I8     2   3PLM   0.94  0.12   0.89  0.14   0.16  0.05     NA    NA
#> 9     C1I9     2   3PLM   0.89  0.12   0.63  0.18   0.21  0.06     NA    NA
#> 10   C1I10     2   3PLM   1.49  0.12   0.13  0.08   0.15  0.04     NA    NA
#> 11    G1I1     2   3PLM   0.96  0.11  -0.46  0.20   0.16  0.07     NA    NA
#> 12    G1I2     2   3PLM   0.90  0.14   1.20  0.14   0.11  0.04     NA    NA
#> 13    G1I3     2   3PLM   1.54  0.26   1.30  0.09   0.19  0.03     NA    NA
#> 14    G1I4     2   3PLM   1.57  0.22   0.28  0.12   0.30  0.04     NA    NA
#> 15    G1I5     2   3PLM   1.35  0.14  -0.17  0.13   0.16  0.05     NA    NA
#> 16    G1I6     2   3PLM   2.17  0.18   0.02  0.06   0.09  0.03     NA    NA
#> 17    G1I7     2   3PLM   1.47  0.16  -0.06  0.12   0.20  0.05     NA    NA
#> 18    G1I8     2   3PLM   2.53  0.47   1.17  0.07   0.32  0.02     NA    NA
#> 19    G1I9     2   3PLM   2.40  0.27  -0.93  0.11   0.24  0.06     NA    NA
#> 20   G1I10     2   3PLM   1.26  0.13  -1.76  0.23   0.23  0.10     NA    NA
#> 21   G1I11     2   3PLM   1.55  0.15  -1.11  0.16   0.21  0.08     NA    NA
#> 22   G1I12     2   3PLM   0.76  0.09  -0.80  0.31   0.20  0.08     NA    NA
#> 23   G1I13     2   3PLM   1.05  0.14  -0.12  0.21   0.21  0.07     NA    NA
#> 24   G1I14     2   3PLM   1.51  0.38   1.74  0.15   0.30  0.02     NA    NA
#> 25   G1I15     2   3PLM   0.87  0.10  -1.43  0.30   0.21  0.09     NA    NA
#> 26   G1I16     2   3PLM   1.05  0.11  -1.94  0.27   0.23  0.10     NA    NA
#> 27   G1I17     2   3PLM   1.06  0.13   0.26  0.15   0.16  0.05     NA    NA
#> 28   G1I18     2   3PLM   2.11  0.21  -0.08  0.08   0.24  0.04     NA    NA
#> 29   G1I19     2   3PLM   1.32  0.13  -1.39  0.19   0.20  0.08     NA    NA
#> 30   G1I20     2   3PLM   1.03  0.18   0.52  0.19   0.23  0.06     NA    NA
#> 31   G1I21     2   3PLM   0.92  0.13   0.76  0.15   0.13  0.05     NA    NA
#> 32   G1I22     2   3PLM   1.79  0.22  -0.63  0.16   0.36  0.06     NA    NA
#> 33   G1I23     2   3PLM   1.31  0.15  -1.16  0.23   0.26  0.09     NA    NA
#> 34   G1I24     2   3PLM   1.66  0.20   0.32  0.09   0.23  0.04     NA    NA
#> 35   G1I25     2   3PLM   1.60  0.18  -0.12  0.12   0.25  0.05     NA    NA
#> 36   G1I26     2   3PLM   1.91  0.26   0.66  0.08   0.25  0.03     NA    NA
#> 37   G1I27     2   3PLM   1.62  0.18  -1.55  0.20   0.27  0.10     NA    NA
#> 38   G1I28     2   3PLM   1.35  0.17   0.57  0.10   0.15  0.04     NA    NA
#> 39   G1I29     2   3PLM   0.92  0.09  -0.37  0.17   0.12  0.06     NA    NA
#> 40   G1I30     2   3PLM   1.06  0.29   2.24  0.25   0.17  0.03     NA    NA
#> 41   G1I31     2   3PLM   2.45  0.49   1.62  0.09   0.18  0.01     NA    NA
#> 42   G1I32     2   3PLM   1.11  0.12  -0.07  0.15   0.15  0.05     NA    NA
#> 43   G1I33     2   3PLM   1.62  0.18   0.18  0.09   0.16  0.04     NA    NA
#> 44   G1I34     2   3PLM   1.34  0.14   0.24  0.09   0.11  0.04     NA    NA
#> 45   G1I35     2   3PLM   1.35  0.17   1.28  0.08   0.07  0.02     NA    NA
#> 46   G1I36     2   3PLM   1.44  0.15  -1.23  0.19   0.22  0.09     NA    NA
#> 47   G1I37     2   3PLM   1.06  0.15  -0.58  0.27   0.27  0.09     NA    NA
#> 48   G1I38     5    GRM   1.06  0.06  -0.37  0.05   0.21  0.05   0.86  0.06
#> 49   C1I11     5    GRM   1.19  0.05  -2.21  0.09  -1.45  0.07  -0.75  0.05
#> 50   C1I12     5    GRM   0.91  0.04  -0.69  0.06   0.03  0.04   0.68  0.04
#> 51    G2I1     2   3PLM   1.76  0.19  -0.86  0.18   0.24  0.10     NA    NA
#> 52    G2I2     2   3PLM   0.84  0.11  -0.44  0.31   0.21  0.09     NA    NA
#> 53    G2I3     2   3PLM   1.09  0.14   0.09  0.20   0.18  0.08     NA    NA
#> 54    G2I4     2   3PLM   1.51  0.29   1.45  0.09   0.18  0.04     NA    NA
#> 55    G2I5     2   3PLM   0.71  0.10  -1.66  0.40   0.20  0.09     NA    NA
#> 56    G2I6     2   3PLM   1.09  0.15  -1.60  0.30   0.22  0.10     NA    NA
#> 57    G2I7     2   3PLM   1.33  0.14   0.26  0.12   0.13  0.05     NA    NA
#> 58    G2I8     2   3PLM   2.22  0.19  -0.09  0.08   0.14  0.05     NA    NA
#> 59    G2I9     2   3PLM   1.08  0.13  -1.56  0.28   0.21  0.09     NA    NA
#> 60   G2I10     2   3PLM   1.65  0.28   0.88  0.11   0.29  0.05     NA    NA
#> 61   G2I11     2   3PLM   0.96  0.13   0.86  0.14   0.11  0.05     NA    NA
#> 62   G2I12     2   3PLM   1.68  0.19  -0.68  0.19   0.26  0.10     NA    NA
#> 63   G2I13     2   3PLM   1.17  0.13  -1.30  0.24   0.21  0.09     NA    NA
#> 64   G2I14     2   3PLM   1.43  0.14   0.22  0.11   0.13  0.05     NA    NA
#> 65   G2I15     2   3PLM   1.58  0.21   0.05  0.17   0.27  0.08     NA    NA
#> 66   G2I16     2   3PLM   1.77  0.24   0.74  0.09   0.22  0.05     NA    NA
#> 67   G2I17     2   3PLM   2.11  0.22  -1.21  0.15   0.20  0.09     NA    NA
#> 68   G2I18     2   3PLM   1.79  0.21   0.63  0.08   0.14  0.04     NA    NA
#> 69   G2I19     2   3PLM   1.12  0.15   0.05  0.22   0.21  0.08     NA    NA
#> 70   G2I20     2   3PLM   1.70  0.47   2.16  0.18   0.16  0.02     NA    NA
#> 71   G2I21     2   3PLM   2.49  0.38   1.53  0.06   0.12  0.02     NA    NA
#> 72   G2I22     2   3PLM   1.50  0.19   0.24  0.15   0.22  0.07     NA    NA
#> 73   G2I23     2   3PLM   2.19  0.28   0.45  0.09   0.28  0.05     NA    NA
#> 74   G2I24     2   3PLM   1.43  0.12   0.36  0.07   0.07  0.03     NA    NA
#> 75   G2I25     2   3PLM   1.80  0.23   1.37  0.06   0.07  0.02     NA    NA
#> 76   G2I26     2   3PLM   1.90  0.20  -0.88  0.17   0.25  0.10     NA    NA
#> 77   G2I27     2   3PLM   0.99  0.12  -0.58  0.25   0.20  0.09     NA    NA
#> 78   G2I28     5    GRM   1.14  0.07  -0.37  0.07   0.14  0.05   0.78  0.05
#> 79    C2I1     2   3PLM   1.07  0.08  -0.29  0.11   0.13  0.04     NA    NA
#> 80    C2I2     2   3PLM   1.00  0.09   1.23  0.07   0.06  0.02     NA    NA
#> 81    C2I3     2   3PLM   1.63  0.16   1.37  0.05   0.13  0.01     NA    NA
#> 82    C2I4     2   3PLM   1.55  0.11   0.24  0.06   0.18  0.02     NA    NA
#> 83    C2I5     2   3PLM   1.43  0.10  -0.07  0.07   0.12  0.03     NA    NA
#> 84    C2I6     2   3PLM   2.04  0.10  -0.06  0.04   0.05  0.01     NA    NA
#> 85    C2I7     2   3PLM   1.56  0.11   0.00  0.07   0.18  0.02     NA    NA
#> 86    C2I8     2   3PLM   1.53  0.15   1.27  0.06   0.19  0.02     NA    NA
#> 87    C2I9     2   3PLM   2.17  0.13  -1.01  0.07   0.14  0.03     NA    NA
#> 88   C2I10     2   3PLM   1.58  0.12  -1.35  0.14   0.26  0.05     NA    NA
#> 89    G3I1     2   3PLM   1.52  0.14  -1.10  0.13   0.16  0.05     NA    NA
#> 90    G3I2     2   3PLM   0.77  0.10  -0.60  0.27   0.19  0.07     NA    NA
#> 91    G3I3     2   3PLM   1.15  0.12   0.08  0.12   0.17  0.03     NA    NA
#> 92    G3I4     2   3PLM   1.31  0.24   1.71  0.12   0.22  0.02     NA    NA
#> 93    G3I5     2   3PLM   0.79  0.11  -1.00  0.36   0.29  0.08     NA    NA
#> 94    G3I6     2   3PLM   1.18  0.15  -1.30  0.25   0.33  0.07     NA    NA
#> 95    G3I7     2   3PLM   1.37  0.13   0.26  0.09   0.14  0.03     NA    NA
#> 96    G3I8     2   3PLM   1.92  0.18  -0.08  0.07   0.15  0.02     NA    NA
#> 97    G3I9     2   3PLM   1.13  0.12  -1.20  0.21   0.21  0.06     NA    NA
#> 98   G3I10     2   3PLM   1.51  0.21   0.97  0.09   0.29  0.02     NA    NA
#> 99   G3I11     2   3PLM   0.90  0.10   0.80  0.11   0.08  0.03     NA    NA
#> 100  G3I12     2   3PLM   1.42  0.17  -0.78  0.15   0.26  0.05     NA    NA
#> 101  G3I13     2   3PLM   1.11  0.13  -1.10  0.22   0.24  0.06     NA    NA
#> 102  G3I14     2   3PLM   1.38  0.14   0.30  0.09   0.17  0.03     NA    NA
#> 103  G3I15     2   3PLM   1.29  0.13  -0.06  0.11   0.19  0.03     NA    NA
#> 104  G3I16     2   3PLM   1.58  0.18   0.74  0.08   0.20  0.02     NA    NA
#> 105  G3I17     2   3PLM   1.61  0.16  -1.55  0.15   0.19  0.06     NA    NA
#> 106  G3I18     2   3PLM   1.38  0.14   0.71  0.08   0.12  0.02     NA    NA
#> 107  G3I19     2   3PLM   0.99  0.11   0.09  0.14   0.17  0.04     NA    NA
#> 108  G3I20     2   3PLM   1.16  0.24   2.27  0.20   0.12  0.02     NA    NA
#> 109  G3I21     2   3PLM   3.06  0.57   1.63  0.06   0.14  0.01     NA    NA
#> 110  G3I22     2   3PLM   1.17  0.11   0.09  0.10   0.10  0.03     NA    NA
#> 111  G3I23     2   3PLM   1.79  0.17   0.32  0.06   0.15  0.02     NA    NA
#> 112  G3I24     2   3PLM   1.13  0.09   0.32  0.08   0.05  0.02     NA    NA
#> 113  G3I25     2   3PLM   1.38  0.15   1.35  0.08   0.05  0.01     NA    NA
#> 114  G3I26     2   3PLM   1.66  0.17  -1.00  0.12   0.22  0.04     NA    NA
#> 115  G3I27     2   3PLM   0.92  0.10  -0.62  0.21   0.18  0.06     NA    NA
#> 116  G3I28     5    GRM   0.94  0.05  -0.38  0.06   0.15  0.06   0.81  0.07
#>      par.5  se.5
#> 1       NA    NA
#> 2       NA    NA
#> 3       NA    NA
#> 4       NA    NA
#> 5       NA    NA
#> 6       NA    NA
#> 7       NA    NA
#> 8       NA    NA
#> 9       NA    NA
#> 10      NA    NA
#> 11      NA    NA
#> 12      NA    NA
#> 13      NA    NA
#> 14      NA    NA
#> 15      NA    NA
#> 16      NA    NA
#> 17      NA    NA
#> 18      NA    NA
#> 19      NA    NA
#> 20      NA    NA
#> 21      NA    NA
#> 22      NA    NA
#> 23      NA    NA
#> 24      NA    NA
#> 25      NA    NA
#> 26      NA    NA
#> 27      NA    NA
#> 28      NA    NA
#> 29      NA    NA
#> 30      NA    NA
#> 31      NA    NA
#> 32      NA    NA
#> 33      NA    NA
#> 34      NA    NA
#> 35      NA    NA
#> 36      NA    NA
#> 37      NA    NA
#> 38      NA    NA
#> 39      NA    NA
#> 40      NA    NA
#> 41      NA    NA
#> 42      NA    NA
#> 43      NA    NA
#> 44      NA    NA
#> 45      NA    NA
#> 46      NA    NA
#> 47      NA    NA
#> 48    1.42  0.08
#> 49   -0.12  0.03
#> 50    1.25  0.06
#> 51      NA    NA
#> 52      NA    NA
#> 53      NA    NA
#> 54      NA    NA
#> 55      NA    NA
#> 56      NA    NA
#> 57      NA    NA
#> 58      NA    NA
#> 59      NA    NA
#> 60      NA    NA
#> 61      NA    NA
#> 62      NA    NA
#> 63      NA    NA
#> 64      NA    NA
#> 65      NA    NA
#> 66      NA    NA
#> 67      NA    NA
#> 68      NA    NA
#> 69      NA    NA
#> 70      NA    NA
#> 71      NA    NA
#> 72      NA    NA
#> 73      NA    NA
#> 74      NA    NA
#> 75      NA    NA
#> 76      NA    NA
#> 77      NA    NA
#> 78    1.44  0.07
#> 79      NA    NA
#> 80      NA    NA
#> 81      NA    NA
#> 82      NA    NA
#> 83      NA    NA
#> 84      NA    NA
#> 85      NA    NA
#> 86      NA    NA
#> 87      NA    NA
#> 88      NA    NA
#> 89      NA    NA
#> 90      NA    NA
#> 91      NA    NA
#> 92      NA    NA
#> 93      NA    NA
#> 94      NA    NA
#> 95      NA    NA
#> 96      NA    NA
#> 97      NA    NA
#> 98      NA    NA
#> 99      NA    NA
#> 100     NA    NA
#> 101     NA    NA
#> 102     NA    NA
#> 103     NA    NA
#> 104     NA    NA
#> 105     NA    NA
#> 106     NA    NA
#> 107     NA    NA
#> 108     NA    NA
#> 109     NA    NA
#> 110     NA    NA
#> 111     NA    NA
#> 112     NA    NA
#> 113     NA    NA
#> 114     NA    NA
#> 115     NA    NA
#> 116   1.55  0.09
#>  Group Parameters: 
#>                      mu  sigma2  sigma
#> estimate(Group1)   0.00    1.00   1.00
#> se(Group1)           NA      NA     NA
#> estimate(Group2)   0.49    0.56   0.75
#> se(Group2)         0.02    0.02   0.01
#> estimate(Group3)  -0.37    1.95   1.40
#> se(Group3)         0.03    0.06   0.02
#> 

# Extract the item parameter estimates
getirt(fit.1, what = "par.est")
#> $overall
#>        id cats model     par.1         par.2       par.3      par.4      par.5
#> 1    C1I1    2  3PLM 0.9002993  1.3630996576  0.27971247         NA         NA
#> 2    C1I2    2  3PLM 2.1246417 -0.9940070191  0.18249617         NA         NA
#> 3    C1I3    2  3PLM 1.0541278  0.6086911591  0.18187156         NA         NA
#> 4    C1I4    2  3PLM 1.0789640 -0.1731865978  0.28172552         NA         NA
#> 5    C1I5    2  3PLM 0.8731636 -0.1539491078  0.18023951         NA         NA
#> 6    C1I6    2  3PLM 1.9100015  0.5956052553  0.08680994         NA         NA
#> 7    C1I7    2  3PLM 1.1042207  1.1038327947  0.15045515         NA         NA
#> 8    C1I8    2  3PLM 0.9407054  0.8862776984  0.15890193         NA         NA
#> 9    C1I9    2  3PLM 0.8926878  0.6255808127  0.20799520         NA         NA
#> 10  C1I10    2  3PLM 1.4872042  0.1323954152  0.15372172         NA         NA
#> 11   G1I1    2  3PLM 0.9648995 -0.4558711256  0.15992372         NA         NA
#> 12   G1I2    2  3PLM 0.8982317  1.2036385734  0.10668307         NA         NA
#> 13   G1I3    2  3PLM 1.5352805  1.2976450244  0.18624929         NA         NA
#> 14   G1I4    2  3PLM 1.5653523  0.2810572501  0.30326526         NA         NA
#> 15   G1I5    2  3PLM 1.3511034 -0.1735122324  0.16083579         NA         NA
#> 16   G1I6    2  3PLM 2.1741319  0.0170068932  0.08620387         NA         NA
#> 17   G1I7    2  3PLM 1.4655950 -0.0639529914  0.20160863         NA         NA
#> 18   G1I8    2  3PLM 2.5274624  1.1707732480  0.32453925         NA         NA
#> 19   G1I9    2  3PLM 2.3986360 -0.9296897943  0.24418284         NA         NA
#> 20  G1I10    2  3PLM 1.2585300 -1.7572134104  0.23391662         NA         NA
#> 21  G1I11    2  3PLM 1.5470937 -1.1118172631  0.20751568         NA         NA
#> 22  G1I12    2  3PLM 0.7556135 -0.8038513481  0.19656034         NA         NA
#> 23  G1I13    2  3PLM 1.0527888 -0.1185472914  0.20738601         NA         NA
#> 24  G1I14    2  3PLM 1.5102932  1.7406214138  0.29709350         NA         NA
#> 25  G1I15    2  3PLM 0.8664971 -1.4323893987  0.21465862         NA         NA
#> 26  G1I16    2  3PLM 1.0542506 -1.9446624563  0.22812622         NA         NA
#> 27  G1I17    2  3PLM 1.0571687  0.2574260233  0.15670798         NA         NA
#> 28  G1I18    2  3PLM 2.1077524 -0.0802270386  0.24053595         NA         NA
#> 29  G1I19    2  3PLM 1.3243225 -1.3873897669  0.20073972         NA         NA
#> 30  G1I20    2  3PLM 1.0303144  0.5198099412  0.22883325         NA         NA
#> 31  G1I21    2  3PLM 0.9230319  0.7599364214  0.13215608         NA         NA
#> 32  G1I22    2  3PLM 1.7870120 -0.6306797460  0.36010326         NA         NA
#> 33  G1I23    2  3PLM 1.3094006 -1.1637343305  0.25562388         NA         NA
#> 34  G1I24    2  3PLM 1.6645514  0.3233691790  0.23193826         NA         NA
#> 35  G1I25    2  3PLM 1.5989408 -0.1175965323  0.25385326         NA         NA
#> 36  G1I26    2  3PLM 1.9085326  0.6574101909  0.25475783         NA         NA
#> 37  G1I27    2  3PLM 1.6232889 -1.5482047548  0.26787427         NA         NA
#> 38  G1I28    2  3PLM 1.3497743  0.5653525105  0.15089732         NA         NA
#> 39  G1I29    2  3PLM 0.9229732 -0.3681170374  0.12451485         NA         NA
#> 40  G1I30    2  3PLM 1.0636199  2.2379302348  0.17259660         NA         NA
#> 41  G1I31    2  3PLM 2.4508248  1.6191425518  0.18202915         NA         NA
#> 42  G1I32    2  3PLM 1.1137113 -0.0722062408  0.14770787         NA         NA
#> 43  G1I33    2  3PLM 1.6242308  0.1841468639  0.16019712         NA         NA
#> 44  G1I34    2  3PLM 1.3396394  0.2438713584  0.10954369         NA         NA
#> 45  G1I35    2  3PLM 1.3498016  1.2826823306  0.06760955         NA         NA
#> 46  G1I36    2  3PLM 1.4377572 -1.2263600168  0.22297804         NA         NA
#> 47  G1I37    2  3PLM 1.0646367 -0.5799937856  0.27097889         NA         NA
#> 48  G1I38    5   GRM 1.0649768 -0.3653605941  0.21365349  0.8571501  1.4220840
#> 49  C1I11    5   GRM 1.1892683 -2.2080814019 -1.44963215 -0.7471002 -0.1188008
#> 50  C1I12    5   GRM 0.9141155 -0.6885203551  0.02631147  0.6789189  1.2487151
#> 51   G2I1    2  3PLM 1.7609541 -0.8578505772  0.24035894         NA         NA
#> 52   G2I2    2  3PLM 0.8436960 -0.4382789597  0.21438678         NA         NA
#> 53   G2I3    2  3PLM 1.0886052  0.0911696948  0.18458253         NA         NA
#> 54   G2I4    2  3PLM 1.5129142  1.4502930939  0.18027120         NA         NA
#> 55   G2I5    2  3PLM 0.7139700 -1.6581797024  0.19781812         NA         NA
#> 56   G2I6    2  3PLM 1.0919249 -1.5999033276  0.21671254         NA         NA
#> 57   G2I7    2  3PLM 1.3325688  0.2633001080  0.13036039         NA         NA
#> 58   G2I8    2  3PLM 2.2199875 -0.0891740329  0.13780823         NA         NA
#> 59   G2I9    2  3PLM 1.0832071 -1.5619697271  0.20682223         NA         NA
#> 60  G2I10    2  3PLM 1.6525972  0.8838581722  0.29075920         NA         NA
#> 61  G2I11    2  3PLM 0.9630202  0.8571579301  0.11294805         NA         NA
#> 62  G2I12    2  3PLM 1.6845080 -0.6793797278  0.25993819         NA         NA
#> 63  G2I13    2  3PLM 1.1744864 -1.3046157584  0.20979661         NA         NA
#> 64  G2I14    2  3PLM 1.4329862  0.2158522635  0.12986846         NA         NA
#> 65  G2I15    2  3PLM 1.5763766  0.0515443553  0.27114790         NA         NA
#> 66  G2I16    2  3PLM 1.7692947  0.7418363902  0.21979684         NA         NA
#> 67  G2I17    2  3PLM 2.1059175 -1.2085830261  0.20089939         NA         NA
#> 68  G2I18    2  3PLM 1.7861465  0.6306602725  0.14248762         NA         NA
#> 69  G2I19    2  3PLM 1.1236963  0.0462753631  0.21071941         NA         NA
#> 70  G2I20    2  3PLM 1.6957930  2.1576947264  0.16005055         NA         NA
#> 71  G2I21    2  3PLM 2.4882719  1.5329597738  0.11756098         NA         NA
#> 72  G2I22    2  3PLM 1.4964966  0.2397808194  0.21725633         NA         NA
#> 73  G2I23    2  3PLM 2.1906673  0.4452744274  0.28388427         NA         NA
#> 74  G2I24    2  3PLM 1.4315517  0.3567207426  0.07277537         NA         NA
#> 75  G2I25    2  3PLM 1.8049099  1.3706377058  0.07249339         NA         NA
#> 76  G2I26    2  3PLM 1.8993430 -0.8763222592  0.24753546         NA         NA
#> 77  G2I27    2  3PLM 0.9922602 -0.5809395819  0.20058057         NA         NA
#> 78  G2I28    5   GRM 1.1428941 -0.3726286563  0.13986303  0.7760909  1.4399545
#> 79   C2I1    2  3PLM 1.0724788 -0.2925185566  0.12706645         NA         NA
#> 80   C2I2    2  3PLM 1.0026458  1.2329410898  0.06038371         NA         NA
#> 81   C2I3    2  3PLM 1.6298126  1.3717691197  0.12593191         NA         NA
#> 82   C2I4    2  3PLM 1.5539481  0.2386833730  0.18401402         NA         NA
#> 83   C2I5    2  3PLM 1.4284843 -0.0690446133  0.11886473         NA         NA
#> 84   C2I6    2  3PLM 2.0356386 -0.0609918393  0.04627366         NA         NA
#> 85   C2I7    2  3PLM 1.5601571  0.0009650233  0.17678863         NA         NA
#> 86   C2I8    2  3PLM 1.5324423  1.2675649742  0.18822170         NA         NA
#> 87   C2I9    2  3PLM 2.1652174 -1.0147828257  0.14241930         NA         NA
#> 88  C2I10    2  3PLM 1.5756487 -1.3546338577  0.26239714         NA         NA
#> 89   G3I1    2  3PLM 1.5177028 -1.1012930873  0.16029848         NA         NA
#> 90   G3I2    2  3PLM 0.7684745 -0.5958229704  0.19239503         NA         NA
#> 91   G3I3    2  3PLM 1.1511024  0.0760256821  0.16953254         NA         NA
#> 92   G3I4    2  3PLM 1.3078975  1.7075770992  0.22441977         NA         NA
#> 93   G3I5    2  3PLM 0.7934450 -1.0028686375  0.28801113         NA         NA
#> 94   G3I6    2  3PLM 1.1835586 -1.2995355425  0.32753494         NA         NA
#> 95   G3I7    2  3PLM 1.3696476  0.2557793125  0.14270858         NA         NA
#> 96   G3I8    2  3PLM 1.9236418 -0.0843458754  0.14793716         NA         NA
#> 97   G3I9    2  3PLM 1.1300723 -1.1978738763  0.21206463         NA         NA
#> 98  G3I10    2  3PLM 1.5126694  0.9703312346  0.29033962         NA         NA
#> 99  G3I11    2  3PLM 0.9021794  0.7968831684  0.08388012         NA         NA
#> 100 G3I12    2  3PLM 1.4214420 -0.7766765567  0.25975672         NA         NA
#> 101 G3I13    2  3PLM 1.1082496 -1.0983147685  0.23504140         NA         NA
#> 102 G3I14    2  3PLM 1.3797649  0.3022772312  0.16690037         NA         NA
#> 103 G3I15    2  3PLM 1.2876135 -0.0573943526  0.19245588         NA         NA
#> 104 G3I16    2  3PLM 1.5824294  0.7372872340  0.20217816         NA         NA
#> 105 G3I17    2  3PLM 1.6058953 -1.5522762073  0.18645088         NA         NA
#> 106 G3I18    2  3PLM 1.3800007  0.7080448815  0.11881653         NA         NA
#> 107 G3I19    2  3PLM 0.9933705  0.0911167343  0.16557109         NA         NA
#> 108 G3I20    2  3PLM 1.1645170  2.2678742087  0.11821851         NA         NA
#> 109 G3I21    2  3PLM 3.0616274  1.6301700618  0.14040136         NA         NA
#> 110 G3I22    2  3PLM 1.1712500  0.0856353360  0.10276551         NA         NA
#> 111 G3I23    2  3PLM 1.7915979  0.3172563764  0.15045425         NA         NA
#> 112 G3I24    2  3PLM 1.1259988  0.3187456984  0.05486984         NA         NA
#> 113 G3I25    2  3PLM 1.3848088  1.3476897615  0.05145709         NA         NA
#> 114 G3I26    2  3PLM 1.6627860 -0.9983631680  0.21786458         NA         NA
#> 115 G3I27    2  3PLM 0.9242740 -0.6199576215  0.18433956         NA         NA
#> 116 G3I28    5   GRM 0.9426521 -0.3795334512  0.14801722  0.8117698  1.5515925
#> 
#> $group
#> $group$Group1
#>       id cats model     par.1       par.2       par.3      par.4      par.5
#> 1   C1I1    2  3PLM 0.9002993  1.36309966  0.27971247         NA         NA
#> 2   C1I2    2  3PLM 2.1246417 -0.99400702  0.18249617         NA         NA
#> 3   C1I3    2  3PLM 1.0541278  0.60869116  0.18187156         NA         NA
#> 4   C1I4    2  3PLM 1.0789640 -0.17318660  0.28172552         NA         NA
#> 5   C1I5    2  3PLM 0.8731636 -0.15394911  0.18023951         NA         NA
#> 6   C1I6    2  3PLM 1.9100015  0.59560526  0.08680994         NA         NA
#> 7   C1I7    2  3PLM 1.1042207  1.10383279  0.15045515         NA         NA
#> 8   C1I8    2  3PLM 0.9407054  0.88627770  0.15890193         NA         NA
#> 9   C1I9    2  3PLM 0.8926878  0.62558081  0.20799520         NA         NA
#> 10 C1I10    2  3PLM 1.4872042  0.13239542  0.15372172         NA         NA
#> 11  G1I1    2  3PLM 0.9648995 -0.45587113  0.15992372         NA         NA
#> 12  G1I2    2  3PLM 0.8982317  1.20363857  0.10668307         NA         NA
#> 13  G1I3    2  3PLM 1.5352805  1.29764502  0.18624929         NA         NA
#> 14  G1I4    2  3PLM 1.5653523  0.28105725  0.30326526         NA         NA
#> 15  G1I5    2  3PLM 1.3511034 -0.17351223  0.16083579         NA         NA
#> 16  G1I6    2  3PLM 2.1741319  0.01700689  0.08620387         NA         NA
#> 17  G1I7    2  3PLM 1.4655950 -0.06395299  0.20160863         NA         NA
#> 18  G1I8    2  3PLM 2.5274624  1.17077325  0.32453925         NA         NA
#> 19  G1I9    2  3PLM 2.3986360 -0.92968979  0.24418284         NA         NA
#> 20 G1I10    2  3PLM 1.2585300 -1.75721341  0.23391662         NA         NA
#> 21 G1I11    2  3PLM 1.5470937 -1.11181726  0.20751568         NA         NA
#> 22 G1I12    2  3PLM 0.7556135 -0.80385135  0.19656034         NA         NA
#> 23 G1I13    2  3PLM 1.0527888 -0.11854729  0.20738601         NA         NA
#> 24 G1I14    2  3PLM 1.5102932  1.74062141  0.29709350         NA         NA
#> 25 G1I15    2  3PLM 0.8664971 -1.43238940  0.21465862         NA         NA
#> 26 G1I16    2  3PLM 1.0542506 -1.94466246  0.22812622         NA         NA
#> 27 G1I17    2  3PLM 1.0571687  0.25742602  0.15670798         NA         NA
#> 28 G1I18    2  3PLM 2.1077524 -0.08022704  0.24053595         NA         NA
#> 29 G1I19    2  3PLM 1.3243225 -1.38738977  0.20073972         NA         NA
#> 30 G1I20    2  3PLM 1.0303144  0.51980994  0.22883325         NA         NA
#> 31 G1I21    2  3PLM 0.9230319  0.75993642  0.13215608         NA         NA
#> 32 G1I22    2  3PLM 1.7870120 -0.63067975  0.36010326         NA         NA
#> 33 G1I23    2  3PLM 1.3094006 -1.16373433  0.25562388         NA         NA
#> 34 G1I24    2  3PLM 1.6645514  0.32336918  0.23193826         NA         NA
#> 35 G1I25    2  3PLM 1.5989408 -0.11759653  0.25385326         NA         NA
#> 36 G1I26    2  3PLM 1.9085326  0.65741019  0.25475783         NA         NA
#> 37 G1I27    2  3PLM 1.6232889 -1.54820475  0.26787427         NA         NA
#> 38 G1I28    2  3PLM 1.3497743  0.56535251  0.15089732         NA         NA
#> 39 G1I29    2  3PLM 0.9229732 -0.36811704  0.12451485         NA         NA
#> 40 G1I30    2  3PLM 1.0636199  2.23793023  0.17259660         NA         NA
#> 41 G1I31    2  3PLM 2.4508248  1.61914255  0.18202915         NA         NA
#> 42 G1I32    2  3PLM 1.1137113 -0.07220624  0.14770787         NA         NA
#> 43 G1I33    2  3PLM 1.6242308  0.18414686  0.16019712         NA         NA
#> 44 G1I34    2  3PLM 1.3396394  0.24387136  0.10954369         NA         NA
#> 45 G1I35    2  3PLM 1.3498016  1.28268233  0.06760955         NA         NA
#> 46 G1I36    2  3PLM 1.4377572 -1.22636002  0.22297804         NA         NA
#> 47 G1I37    2  3PLM 1.0646367 -0.57999379  0.27097889         NA         NA
#> 48 G1I38    5   GRM 1.0649768 -0.36536059  0.21365349  0.8571501  1.4220840
#> 49 C1I11    5   GRM 1.1892683 -2.20808140 -1.44963215 -0.7471002 -0.1188008
#> 50 C1I12    5   GRM 0.9141155 -0.68852036  0.02631147  0.6789189  1.2487151
#> 
#> $group$Group2
#>       id cats model     par.1         par.2       par.3      par.4      par.5
#> 1   C1I1    2  3PLM 0.9002993  1.3630996576  0.27971247         NA         NA
#> 2   C1I2    2  3PLM 2.1246417 -0.9940070191  0.18249617         NA         NA
#> 3   C1I3    2  3PLM 1.0541278  0.6086911591  0.18187156         NA         NA
#> 4   C1I4    2  3PLM 1.0789640 -0.1731865978  0.28172552         NA         NA
#> 5   C1I5    2  3PLM 0.8731636 -0.1539491078  0.18023951         NA         NA
#> 6   C1I6    2  3PLM 1.9100015  0.5956052553  0.08680994         NA         NA
#> 7   C1I7    2  3PLM 1.1042207  1.1038327947  0.15045515         NA         NA
#> 8   C1I8    2  3PLM 0.9407054  0.8862776984  0.15890193         NA         NA
#> 9   C1I9    2  3PLM 0.8926878  0.6255808127  0.20799520         NA         NA
#> 10 C1I10    2  3PLM 1.4872042  0.1323954152  0.15372172         NA         NA
#> 11 C1I11    5   GRM 1.1892683 -2.2080814019 -1.44963215 -0.7471002 -0.1188008
#> 12 C1I12    5   GRM 0.9141155 -0.6885203551  0.02631147  0.6789189  1.2487151
#> 13  G2I1    2  3PLM 1.7609541 -0.8578505772  0.24035894         NA         NA
#> 14  G2I2    2  3PLM 0.8436960 -0.4382789597  0.21438678         NA         NA
#> 15  G2I3    2  3PLM 1.0886052  0.0911696948  0.18458253         NA         NA
#> 16  G2I4    2  3PLM 1.5129142  1.4502930939  0.18027120         NA         NA
#> 17  G2I5    2  3PLM 0.7139700 -1.6581797024  0.19781812         NA         NA
#> 18  G2I6    2  3PLM 1.0919249 -1.5999033276  0.21671254         NA         NA
#> 19  G2I7    2  3PLM 1.3325688  0.2633001080  0.13036039         NA         NA
#> 20  G2I8    2  3PLM 2.2199875 -0.0891740329  0.13780823         NA         NA
#> 21  G2I9    2  3PLM 1.0832071 -1.5619697271  0.20682223         NA         NA
#> 22 G2I10    2  3PLM 1.6525972  0.8838581722  0.29075920         NA         NA
#> 23 G2I11    2  3PLM 0.9630202  0.8571579301  0.11294805         NA         NA
#> 24 G2I12    2  3PLM 1.6845080 -0.6793797278  0.25993819         NA         NA
#> 25 G2I13    2  3PLM 1.1744864 -1.3046157584  0.20979661         NA         NA
#> 26 G2I14    2  3PLM 1.4329862  0.2158522635  0.12986846         NA         NA
#> 27 G2I15    2  3PLM 1.5763766  0.0515443553  0.27114790         NA         NA
#> 28 G2I16    2  3PLM 1.7692947  0.7418363902  0.21979684         NA         NA
#> 29 G2I17    2  3PLM 2.1059175 -1.2085830261  0.20089939         NA         NA
#> 30 G2I18    2  3PLM 1.7861465  0.6306602725  0.14248762         NA         NA
#> 31 G2I19    2  3PLM 1.1236963  0.0462753631  0.21071941         NA         NA
#> 32 G2I20    2  3PLM 1.6957930  2.1576947264  0.16005055         NA         NA
#> 33 G2I21    2  3PLM 2.4882719  1.5329597738  0.11756098         NA         NA
#> 34 G2I22    2  3PLM 1.4964966  0.2397808194  0.21725633         NA         NA
#> 35 G2I23    2  3PLM 2.1906673  0.4452744274  0.28388427         NA         NA
#> 36 G2I24    2  3PLM 1.4315517  0.3567207426  0.07277537         NA         NA
#> 37 G2I25    2  3PLM 1.8049099  1.3706377058  0.07249339         NA         NA
#> 38 G2I26    2  3PLM 1.8993430 -0.8763222592  0.24753546         NA         NA
#> 39 G2I27    2  3PLM 0.9922602 -0.5809395819  0.20058057         NA         NA
#> 40 G2I28    5   GRM 1.1428941 -0.3726286563  0.13986303  0.7760909  1.4399545
#> 41  C2I1    2  3PLM 1.0724788 -0.2925185566  0.12706645         NA         NA
#> 42  C2I2    2  3PLM 1.0026458  1.2329410898  0.06038371         NA         NA
#> 43  C2I3    2  3PLM 1.6298126  1.3717691197  0.12593191         NA         NA
#> 44  C2I4    2  3PLM 1.5539481  0.2386833730  0.18401402         NA         NA
#> 45  C2I5    2  3PLM 1.4284843 -0.0690446133  0.11886473         NA         NA
#> 46  C2I6    2  3PLM 2.0356386 -0.0609918393  0.04627366         NA         NA
#> 47  C2I7    2  3PLM 1.5601571  0.0009650233  0.17678863         NA         NA
#> 48  C2I8    2  3PLM 1.5324423  1.2675649742  0.18822170         NA         NA
#> 49  C2I9    2  3PLM 2.1652174 -1.0147828257  0.14241930         NA         NA
#> 50 C2I10    2  3PLM 1.5756487 -1.3546338577  0.26239714         NA         NA
#> 
#> $group$Group3
#>       id cats model     par.1         par.2      par.3     par.4    par.5
#> 1   C2I1    2  3PLM 1.0724788 -0.2925185566 0.12706645        NA       NA
#> 2   C2I2    2  3PLM 1.0026458  1.2329410898 0.06038371        NA       NA
#> 3   C2I3    2  3PLM 1.6298126  1.3717691197 0.12593191        NA       NA
#> 4   C2I4    2  3PLM 1.5539481  0.2386833730 0.18401402        NA       NA
#> 5   C2I5    2  3PLM 1.4284843 -0.0690446133 0.11886473        NA       NA
#> 6   C2I6    2  3PLM 2.0356386 -0.0609918393 0.04627366        NA       NA
#> 7   C2I7    2  3PLM 1.5601571  0.0009650233 0.17678863        NA       NA
#> 8   C2I8    2  3PLM 1.5324423  1.2675649742 0.18822170        NA       NA
#> 9   C2I9    2  3PLM 2.1652174 -1.0147828257 0.14241930        NA       NA
#> 10 C2I10    2  3PLM 1.5756487 -1.3546338577 0.26239714        NA       NA
#> 11  G3I1    2  3PLM 1.5177028 -1.1012930873 0.16029848        NA       NA
#> 12  G3I2    2  3PLM 0.7684745 -0.5958229704 0.19239503        NA       NA
#> 13  G3I3    2  3PLM 1.1511024  0.0760256821 0.16953254        NA       NA
#> 14  G3I4    2  3PLM 1.3078975  1.7075770992 0.22441977        NA       NA
#> 15  G3I5    2  3PLM 0.7934450 -1.0028686375 0.28801113        NA       NA
#> 16  G3I6    2  3PLM 1.1835586 -1.2995355425 0.32753494        NA       NA
#> 17  G3I7    2  3PLM 1.3696476  0.2557793125 0.14270858        NA       NA
#> 18  G3I8    2  3PLM 1.9236418 -0.0843458754 0.14793716        NA       NA
#> 19  G3I9    2  3PLM 1.1300723 -1.1978738763 0.21206463        NA       NA
#> 20 G3I10    2  3PLM 1.5126694  0.9703312346 0.29033962        NA       NA
#> 21 G3I11    2  3PLM 0.9021794  0.7968831684 0.08388012        NA       NA
#> 22 G3I12    2  3PLM 1.4214420 -0.7766765567 0.25975672        NA       NA
#> 23 G3I13    2  3PLM 1.1082496 -1.0983147685 0.23504140        NA       NA
#> 24 G3I14    2  3PLM 1.3797649  0.3022772312 0.16690037        NA       NA
#> 25 G3I15    2  3PLM 1.2876135 -0.0573943526 0.19245588        NA       NA
#> 26 G3I16    2  3PLM 1.5824294  0.7372872340 0.20217816        NA       NA
#> 27 G3I17    2  3PLM 1.6058953 -1.5522762073 0.18645088        NA       NA
#> 28 G3I18    2  3PLM 1.3800007  0.7080448815 0.11881653        NA       NA
#> 29 G3I19    2  3PLM 0.9933705  0.0911167343 0.16557109        NA       NA
#> 30 G3I20    2  3PLM 1.1645170  2.2678742087 0.11821851        NA       NA
#> 31 G3I21    2  3PLM 3.0616274  1.6301700618 0.14040136        NA       NA
#> 32 G3I22    2  3PLM 1.1712500  0.0856353360 0.10276551        NA       NA
#> 33 G3I23    2  3PLM 1.7915979  0.3172563764 0.15045425        NA       NA
#> 34 G3I24    2  3PLM 1.1259988  0.3187456984 0.05486984        NA       NA
#> 35 G3I25    2  3PLM 1.3848088  1.3476897615 0.05145709        NA       NA
#> 36 G3I26    2  3PLM 1.6627860 -0.9983631680 0.21786458        NA       NA
#> 37 G3I27    2  3PLM 0.9242740 -0.6199576215 0.18433956        NA       NA
#> 38 G3I28    5   GRM 0.9426521 -0.3795334512 0.14801722 0.8117698 1.551592
#> 
#> 

# Extract the standard error estimates
getirt(fit.1, what = "se.est")
#> $overall
#>        id cats model      par.1      par.2      par.3      par.4      par.5
#> 1    C1I1    2  3PLM 0.17989519 0.16035224 0.04893744         NA         NA
#> 2    C1I2    2  3PLM 0.14509575 0.09689007 0.05855244         NA         NA
#> 3    C1I3    2  3PLM 0.12333025 0.13013772 0.04568686         NA         NA
#> 4    C1I4    2  3PLM 0.12539952 0.21610221 0.06984577         NA         NA
#> 5    C1I5    2  3PLM 0.09159813 0.22015093 0.06797074         NA         NA
#> 6    C1I6    2  3PLM 0.13215821 0.04030740 0.01860230         NA         NA
#> 7    C1I7    2  3PLM 0.13398772 0.08891913 0.03307573         NA         NA
#> 8    C1I8    2  3PLM 0.12334991 0.13614025 0.04646945         NA         NA
#> 9    C1I9    2  3PLM 0.12236592 0.18495675 0.05639681         NA         NA
#> 10  C1I10    2  3PLM 0.11760294 0.08417251 0.03717477         NA         NA
#> 11   G1I1    2  3PLM 0.10655826 0.20217054 0.06672572         NA         NA
#> 12   G1I2    2  3PLM 0.13694638 0.13626565 0.03909387         NA         NA
#> 13   G1I3    2  3PLM 0.26266869 0.09156215 0.02556423         NA         NA
#> 14   G1I4    2  3PLM 0.21896209 0.11821220 0.04253288         NA         NA
#> 15   G1I5    2  3PLM 0.14227648 0.12915817 0.05186212         NA         NA
#> 16   G1I6    2  3PLM 0.18452247 0.05538028 0.02621436         NA         NA
#> 17   G1I7    2  3PLM 0.16049350 0.12171798 0.04882272         NA         NA
#> 18   G1I8    2  3PLM 0.46550649 0.06951993 0.01895467         NA         NA
#> 19   G1I9    2  3PLM 0.26711886 0.11407326 0.06322653         NA         NA
#> 20  G1I10    2  3PLM 0.12545580 0.22750745 0.09762660         NA         NA
#> 21  G1I11    2  3PLM 0.15057675 0.16486816 0.07867466         NA         NA
#> 22  G1I12    2  3PLM 0.09367336 0.31303395 0.08422275         NA         NA
#> 23  G1I13    2  3PLM 0.14115709 0.20788070 0.06872976         NA         NA
#> 24  G1I14    2  3PLM 0.38235203 0.15053141 0.02485445         NA         NA
#> 25  G1I15    2  3PLM 0.09621277 0.29800051 0.09181670         NA         NA
#> 26  G1I16    2  3PLM 0.11062833 0.26571952 0.09735725         NA         NA
#> 27  G1I17    2  3PLM 0.13151860 0.15421189 0.05197027         NA         NA
#> 28  G1I18    2  3PLM 0.20611213 0.08024248 0.03657769         NA         NA
#> 29  G1I19    2  3PLM 0.12685072 0.19375886 0.08485688         NA         NA
#> 30  G1I20    2  3PLM 0.18089581 0.18707909 0.05966072         NA         NA
#> 31  G1I21    2  3PLM 0.13455458 0.14996962 0.04703101         NA         NA
#> 32  G1I22    2  3PLM 0.21911500 0.15737913 0.06424983         NA         NA
#> 33  G1I23    2  3PLM 0.15053755 0.22982042 0.09322036         NA         NA
#> 34  G1I24    2  3PLM 0.19679685 0.09411360 0.03689301         NA         NA
#> 35  G1I25    2  3PLM 0.18243227 0.12257416 0.04913825         NA         NA
#> 36  G1I26    2  3PLM 0.25542584 0.07703704 0.02893084         NA         NA
#> 37  G1I27    2  3PLM 0.17691915 0.20386267 0.10127163         NA         NA
#> 38  G1I28    2  3PLM 0.16553187 0.09754762 0.03666743         NA         NA
#> 39  G1I29    2  3PLM 0.09087413 0.17161992 0.05512983         NA         NA
#> 40  G1I30    2  3PLM 0.28799644 0.24787264 0.02955840         NA         NA
#> 41  G1I31    2  3PLM 0.48866407 0.09220125 0.01439037         NA         NA
#> 42  G1I32    2  3PLM 0.12263290 0.15178218 0.05400259         NA         NA
#> 43  G1I33    2  3PLM 0.17530464 0.08877908 0.03732072         NA         NA
#> 44  G1I34    2  3PLM 0.13610866 0.09412723 0.03674769         NA         NA
#> 45  G1I35    2  3PLM 0.17137910 0.08344182 0.02039465         NA         NA
#> 46  G1I36    2  3PLM 0.14523548 0.19046168 0.08649362         NA         NA
#> 47  G1I37    2  3PLM 0.14678991 0.27080023 0.08750606         NA         NA
#> 48  G1I38    5   GRM 0.05824336 0.05464783 0.05036679 0.06337710 0.08439014
#> 49  C1I11    5   GRM 0.04768216 0.09429091 0.06737581 0.04635165 0.03448717
#> 50  C1I12    5   GRM 0.04132681 0.05629429 0.04082063 0.04311500 0.05653026
#> 51   G2I1    2  3PLM 0.19350221 0.17643375 0.09757730         NA         NA
#> 52   G2I2    2  3PLM 0.11489091 0.30524136 0.09207477         NA         NA
#> 53   G2I3    2  3PLM 0.13653157 0.20117207 0.07586750         NA         NA
#> 54   G2I4    2  3PLM 0.29245623 0.09333564 0.03808894         NA         NA
#> 55   G2I5    2  3PLM 0.10481807 0.39787167 0.08903051         NA         NA
#> 56   G2I6    2  3PLM 0.14518031 0.29678044 0.09513691         NA         NA
#> 57   G2I7    2  3PLM 0.13785758 0.12130037 0.05478446         NA         NA
#> 58   G2I8    2  3PLM 0.19269491 0.08072746 0.05145801         NA         NA
#> 59   G2I9    2  3PLM 0.13251099 0.27721615 0.09193314         NA         NA
#> 60  G2I10    2  3PLM 0.28338383 0.11486665 0.05004309         NA         NA
#> 61  G2I11    2  3PLM 0.12959812 0.14260521 0.05017878         NA         NA
#> 62  G2I12    2  3PLM 0.19165147 0.18551525 0.09840785         NA         NA
#> 63  G2I13    2  3PLM 0.13361413 0.24090086 0.09247457         NA         NA
#> 64  G2I14    2  3PLM 0.14344036 0.11234266 0.05306700         NA         NA
#> 65  G2I15    2  3PLM 0.20783033 0.17092977 0.07863958         NA         NA
#> 66  G2I16    2  3PLM 0.24141540 0.09487270 0.04614521         NA         NA
#> 67  G2I17    2  3PLM 0.21812356 0.14792918 0.08902957         NA         NA
#> 68  G2I18    2  3PLM 0.20825056 0.08176892 0.04319252         NA         NA
#> 69  G2I19    2  3PLM 0.14974281 0.21807451 0.08279427         NA         NA
#> 70  G2I20    2  3PLM 0.47112164 0.18433168 0.02351512         NA         NA
#> 71  G2I21    2  3PLM 0.37699918 0.05953101 0.01720186         NA         NA
#> 72  G2I22    2  3PLM 0.18888770 0.14944483 0.06872351         NA         NA
#> 73  G2I23    2  3PLM 0.28077865 0.09153703 0.04827635         NA         NA
#> 74  G2I24    2  3PLM 0.12130175 0.07457793 0.03386516         NA         NA
#> 75  G2I25    2  3PLM 0.22582407 0.06157614 0.02227878         NA         NA
#> 76  G2I26    2  3PLM 0.20134542 0.16552101 0.09828502         NA         NA
#> 77  G2I27    2  3PLM 0.11980495 0.24841774 0.08769248         NA         NA
#> 78  G2I28    5   GRM 0.07259152 0.06782688 0.04972658 0.04738506 0.07009125
#> 79   C2I1    2  3PLM 0.07967461 0.10986098 0.03524967         NA         NA
#> 80   C2I2    2  3PLM 0.08692171 0.06717781 0.01809171         NA         NA
#> 81   C2I3    2  3PLM 0.15623585 0.05282304 0.01408261         NA         NA
#> 82   C2I4    2  3PLM 0.11189913 0.06186328 0.02306786         NA         NA
#> 83   C2I5    2  3PLM 0.09565798 0.06810561 0.02521438         NA         NA
#> 84   C2I6    2  3PLM 0.10326422 0.03696659 0.01351802         NA         NA
#> 85   C2I7    2  3PLM 0.10574626 0.06606876 0.02489720         NA         NA
#> 86   C2I8    2  3PLM 0.15231896 0.05723105 0.01673332         NA         NA
#> 87   C2I9    2  3PLM 0.13437617 0.07069516 0.03095926         NA         NA
#> 88  C2I10    2  3PLM 0.12209895 0.14135303 0.05284333         NA         NA
#> 89   G3I1    2  3PLM 0.14295950 0.12757467 0.04598806         NA         NA
#> 90   G3I2    2  3PLM 0.09806822 0.27491495 0.06605199         NA         NA
#> 91   G3I3    2  3PLM 0.12223986 0.11810373 0.03455361         NA         NA
#> 92   G3I4    2  3PLM 0.23674338 0.12489657 0.02042690         NA         NA
#> 93   G3I5    2  3PLM 0.11027073 0.35844706 0.08124875         NA         NA
#> 94   G3I6    2  3PLM 0.15212744 0.25219492 0.07173387         NA         NA
#> 95   G3I7    2  3PLM 0.13192537 0.08582984 0.02599794         NA         NA
#> 96   G3I8    2  3PLM 0.17753632 0.06568100 0.02302363         NA         NA
#> 97   G3I9    2  3PLM 0.12073992 0.20628990 0.06412558         NA         NA
#> 98  G3I10    2  3PLM 0.21325449 0.09460878 0.02268460         NA         NA
#> 99  G3I11    2  3PLM 0.10223560 0.11256805 0.02810817         NA         NA
#> 100 G3I12    2  3PLM 0.16809784 0.14982785 0.04805076         NA         NA
#> 101 G3I13    2  3PLM 0.12969261 0.21859906 0.06476116         NA         NA
#> 102 G3I14    2  3PLM 0.14184270 0.08769788 0.02633227         NA         NA
#> 103 G3I15    2  3PLM 0.13462872 0.11103389 0.03394777         NA         NA
#> 104 G3I16    2  3PLM 0.17691571 0.07758799 0.02083063         NA         NA
#> 105 G3I17    2  3PLM 0.15688465 0.15428884 0.05926460         NA         NA
#> 106 G3I18    2  3PLM 0.14382506 0.07697407 0.02111022         NA         NA
#> 107 G3I19    2  3PLM 0.11393846 0.14464817 0.04000362         NA         NA
#> 108 G3I20    2  3PLM 0.23541976 0.19741351 0.01684581         NA         NA
#> 109 G3I21    2  3PLM 0.57479199 0.06270311 0.01022201         NA         NA
#> 110 G3I22    2  3PLM 0.11463786 0.09814469 0.02978440         NA         NA
#> 111 G3I23    2  3PLM 0.16694784 0.06445158 0.01994484         NA         NA
#> 112 G3I24    2  3PLM 0.09417036 0.07931539 0.02063766         NA         NA
#> 113 G3I25    2  3PLM 0.14627533 0.07996861 0.01264067         NA         NA
#> 114 G3I26    2  3PLM 0.16594755 0.12414916 0.04434197         NA         NA
#> 115 G3I27    2  3PLM 0.10437934 0.20932778 0.05657496         NA         NA
#> 116 G3I28    5   GRM 0.04771674 0.05965450 0.05957007 0.07120741 0.09365933
#> 
#> $group
#> $group$Group1
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   C1I1    2  3PLM 0.17989519 0.16035224 0.04893744         NA         NA
#> 2   C1I2    2  3PLM 0.14509575 0.09689007 0.05855244         NA         NA
#> 3   C1I3    2  3PLM 0.12333025 0.13013772 0.04568686         NA         NA
#> 4   C1I4    2  3PLM 0.12539952 0.21610221 0.06984577         NA         NA
#> 5   C1I5    2  3PLM 0.09159813 0.22015093 0.06797074         NA         NA
#> 6   C1I6    2  3PLM 0.13215821 0.04030740 0.01860230         NA         NA
#> 7   C1I7    2  3PLM 0.13398772 0.08891913 0.03307573         NA         NA
#> 8   C1I8    2  3PLM 0.12334991 0.13614025 0.04646945         NA         NA
#> 9   C1I9    2  3PLM 0.12236592 0.18495675 0.05639681         NA         NA
#> 10 C1I10    2  3PLM 0.11760294 0.08417251 0.03717477         NA         NA
#> 11  G1I1    2  3PLM 0.10655826 0.20217054 0.06672572         NA         NA
#> 12  G1I2    2  3PLM 0.13694638 0.13626565 0.03909387         NA         NA
#> 13  G1I3    2  3PLM 0.26266869 0.09156215 0.02556423         NA         NA
#> 14  G1I4    2  3PLM 0.21896209 0.11821220 0.04253288         NA         NA
#> 15  G1I5    2  3PLM 0.14227648 0.12915817 0.05186212         NA         NA
#> 16  G1I6    2  3PLM 0.18452247 0.05538028 0.02621436         NA         NA
#> 17  G1I7    2  3PLM 0.16049350 0.12171798 0.04882272         NA         NA
#> 18  G1I8    2  3PLM 0.46550649 0.06951993 0.01895467         NA         NA
#> 19  G1I9    2  3PLM 0.26711886 0.11407326 0.06322653         NA         NA
#> 20 G1I10    2  3PLM 0.12545580 0.22750745 0.09762660         NA         NA
#> 21 G1I11    2  3PLM 0.15057675 0.16486816 0.07867466         NA         NA
#> 22 G1I12    2  3PLM 0.09367336 0.31303395 0.08422275         NA         NA
#> 23 G1I13    2  3PLM 0.14115709 0.20788070 0.06872976         NA         NA
#> 24 G1I14    2  3PLM 0.38235203 0.15053141 0.02485445         NA         NA
#> 25 G1I15    2  3PLM 0.09621277 0.29800051 0.09181670         NA         NA
#> 26 G1I16    2  3PLM 0.11062833 0.26571952 0.09735725         NA         NA
#> 27 G1I17    2  3PLM 0.13151860 0.15421189 0.05197027         NA         NA
#> 28 G1I18    2  3PLM 0.20611213 0.08024248 0.03657769         NA         NA
#> 29 G1I19    2  3PLM 0.12685072 0.19375886 0.08485688         NA         NA
#> 30 G1I20    2  3PLM 0.18089581 0.18707909 0.05966072         NA         NA
#> 31 G1I21    2  3PLM 0.13455458 0.14996962 0.04703101         NA         NA
#> 32 G1I22    2  3PLM 0.21911500 0.15737913 0.06424983         NA         NA
#> 33 G1I23    2  3PLM 0.15053755 0.22982042 0.09322036         NA         NA
#> 34 G1I24    2  3PLM 0.19679685 0.09411360 0.03689301         NA         NA
#> 35 G1I25    2  3PLM 0.18243227 0.12257416 0.04913825         NA         NA
#> 36 G1I26    2  3PLM 0.25542584 0.07703704 0.02893084         NA         NA
#> 37 G1I27    2  3PLM 0.17691915 0.20386267 0.10127163         NA         NA
#> 38 G1I28    2  3PLM 0.16553187 0.09754762 0.03666743         NA         NA
#> 39 G1I29    2  3PLM 0.09087413 0.17161992 0.05512983         NA         NA
#> 40 G1I30    2  3PLM 0.28799644 0.24787264 0.02955840         NA         NA
#> 41 G1I31    2  3PLM 0.48866407 0.09220125 0.01439037         NA         NA
#> 42 G1I32    2  3PLM 0.12263290 0.15178218 0.05400259         NA         NA
#> 43 G1I33    2  3PLM 0.17530464 0.08877908 0.03732072         NA         NA
#> 44 G1I34    2  3PLM 0.13610866 0.09412723 0.03674769         NA         NA
#> 45 G1I35    2  3PLM 0.17137910 0.08344182 0.02039465         NA         NA
#> 46 G1I36    2  3PLM 0.14523548 0.19046168 0.08649362         NA         NA
#> 47 G1I37    2  3PLM 0.14678991 0.27080023 0.08750606         NA         NA
#> 48 G1I38    5   GRM 0.05824336 0.05464783 0.05036679 0.06337710 0.08439014
#> 49 C1I11    5   GRM 0.04768216 0.09429091 0.06737581 0.04635165 0.03448717
#> 50 C1I12    5   GRM 0.04132681 0.05629429 0.04082063 0.04311500 0.05653026
#> 
#> $group$Group2
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   C1I1    2  3PLM 0.17989519 0.16035224 0.04893744         NA         NA
#> 2   C1I2    2  3PLM 0.14509575 0.09689007 0.05855244         NA         NA
#> 3   C1I3    2  3PLM 0.12333025 0.13013772 0.04568686         NA         NA
#> 4   C1I4    2  3PLM 0.12539952 0.21610221 0.06984577         NA         NA
#> 5   C1I5    2  3PLM 0.09159813 0.22015093 0.06797074         NA         NA
#> 6   C1I6    2  3PLM 0.13215821 0.04030740 0.01860230         NA         NA
#> 7   C1I7    2  3PLM 0.13398772 0.08891913 0.03307573         NA         NA
#> 8   C1I8    2  3PLM 0.12334991 0.13614025 0.04646945         NA         NA
#> 9   C1I9    2  3PLM 0.12236592 0.18495675 0.05639681         NA         NA
#> 10 C1I10    2  3PLM 0.11760294 0.08417251 0.03717477         NA         NA
#> 11 C1I11    5   GRM 0.04768216 0.09429091 0.06737581 0.04635165 0.03448717
#> 12 C1I12    5   GRM 0.04132681 0.05629429 0.04082063 0.04311500 0.05653026
#> 13  G2I1    2  3PLM 0.19350221 0.17643375 0.09757730         NA         NA
#> 14  G2I2    2  3PLM 0.11489091 0.30524136 0.09207477         NA         NA
#> 15  G2I3    2  3PLM 0.13653157 0.20117207 0.07586750         NA         NA
#> 16  G2I4    2  3PLM 0.29245623 0.09333564 0.03808894         NA         NA
#> 17  G2I5    2  3PLM 0.10481807 0.39787167 0.08903051         NA         NA
#> 18  G2I6    2  3PLM 0.14518031 0.29678044 0.09513691         NA         NA
#> 19  G2I7    2  3PLM 0.13785758 0.12130037 0.05478446         NA         NA
#> 20  G2I8    2  3PLM 0.19269491 0.08072746 0.05145801         NA         NA
#> 21  G2I9    2  3PLM 0.13251099 0.27721615 0.09193314         NA         NA
#> 22 G2I10    2  3PLM 0.28338383 0.11486665 0.05004309         NA         NA
#> 23 G2I11    2  3PLM 0.12959812 0.14260521 0.05017878         NA         NA
#> 24 G2I12    2  3PLM 0.19165147 0.18551525 0.09840785         NA         NA
#> 25 G2I13    2  3PLM 0.13361413 0.24090086 0.09247457         NA         NA
#> 26 G2I14    2  3PLM 0.14344036 0.11234266 0.05306700         NA         NA
#> 27 G2I15    2  3PLM 0.20783033 0.17092977 0.07863958         NA         NA
#> 28 G2I16    2  3PLM 0.24141540 0.09487270 0.04614521         NA         NA
#> 29 G2I17    2  3PLM 0.21812356 0.14792918 0.08902957         NA         NA
#> 30 G2I18    2  3PLM 0.20825056 0.08176892 0.04319252         NA         NA
#> 31 G2I19    2  3PLM 0.14974281 0.21807451 0.08279427         NA         NA
#> 32 G2I20    2  3PLM 0.47112164 0.18433168 0.02351512         NA         NA
#> 33 G2I21    2  3PLM 0.37699918 0.05953101 0.01720186         NA         NA
#> 34 G2I22    2  3PLM 0.18888770 0.14944483 0.06872351         NA         NA
#> 35 G2I23    2  3PLM 0.28077865 0.09153703 0.04827635         NA         NA
#> 36 G2I24    2  3PLM 0.12130175 0.07457793 0.03386516         NA         NA
#> 37 G2I25    2  3PLM 0.22582407 0.06157614 0.02227878         NA         NA
#> 38 G2I26    2  3PLM 0.20134542 0.16552101 0.09828502         NA         NA
#> 39 G2I27    2  3PLM 0.11980495 0.24841774 0.08769248         NA         NA
#> 40 G2I28    5   GRM 0.07259152 0.06782688 0.04972658 0.04738506 0.07009125
#> 41  C2I1    2  3PLM 0.07967461 0.10986098 0.03524967         NA         NA
#> 42  C2I2    2  3PLM 0.08692171 0.06717781 0.01809171         NA         NA
#> 43  C2I3    2  3PLM 0.15623585 0.05282304 0.01408261         NA         NA
#> 44  C2I4    2  3PLM 0.11189913 0.06186328 0.02306786         NA         NA
#> 45  C2I5    2  3PLM 0.09565798 0.06810561 0.02521438         NA         NA
#> 46  C2I6    2  3PLM 0.10326422 0.03696659 0.01351802         NA         NA
#> 47  C2I7    2  3PLM 0.10574626 0.06606876 0.02489720         NA         NA
#> 48  C2I8    2  3PLM 0.15231896 0.05723105 0.01673332         NA         NA
#> 49  C2I9    2  3PLM 0.13437617 0.07069516 0.03095926         NA         NA
#> 50 C2I10    2  3PLM 0.12209895 0.14135303 0.05284333         NA         NA
#> 
#> $group$Group3
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   C2I1    2  3PLM 0.07967461 0.10986098 0.03524967         NA         NA
#> 2   C2I2    2  3PLM 0.08692171 0.06717781 0.01809171         NA         NA
#> 3   C2I3    2  3PLM 0.15623585 0.05282304 0.01408261         NA         NA
#> 4   C2I4    2  3PLM 0.11189913 0.06186328 0.02306786         NA         NA
#> 5   C2I5    2  3PLM 0.09565798 0.06810561 0.02521438         NA         NA
#> 6   C2I6    2  3PLM 0.10326422 0.03696659 0.01351802         NA         NA
#> 7   C2I7    2  3PLM 0.10574626 0.06606876 0.02489720         NA         NA
#> 8   C2I8    2  3PLM 0.15231896 0.05723105 0.01673332         NA         NA
#> 9   C2I9    2  3PLM 0.13437617 0.07069516 0.03095926         NA         NA
#> 10 C2I10    2  3PLM 0.12209895 0.14135303 0.05284333         NA         NA
#> 11  G3I1    2  3PLM 0.14295950 0.12757467 0.04598806         NA         NA
#> 12  G3I2    2  3PLM 0.09806822 0.27491495 0.06605199         NA         NA
#> 13  G3I3    2  3PLM 0.12223986 0.11810373 0.03455361         NA         NA
#> 14  G3I4    2  3PLM 0.23674338 0.12489657 0.02042690         NA         NA
#> 15  G3I5    2  3PLM 0.11027073 0.35844706 0.08124875         NA         NA
#> 16  G3I6    2  3PLM 0.15212744 0.25219492 0.07173387         NA         NA
#> 17  G3I7    2  3PLM 0.13192537 0.08582984 0.02599794         NA         NA
#> 18  G3I8    2  3PLM 0.17753632 0.06568100 0.02302363         NA         NA
#> 19  G3I9    2  3PLM 0.12073992 0.20628990 0.06412558         NA         NA
#> 20 G3I10    2  3PLM 0.21325449 0.09460878 0.02268460         NA         NA
#> 21 G3I11    2  3PLM 0.10223560 0.11256805 0.02810817         NA         NA
#> 22 G3I12    2  3PLM 0.16809784 0.14982785 0.04805076         NA         NA
#> 23 G3I13    2  3PLM 0.12969261 0.21859906 0.06476116         NA         NA
#> 24 G3I14    2  3PLM 0.14184270 0.08769788 0.02633227         NA         NA
#> 25 G3I15    2  3PLM 0.13462872 0.11103389 0.03394777         NA         NA
#> 26 G3I16    2  3PLM 0.17691571 0.07758799 0.02083063         NA         NA
#> 27 G3I17    2  3PLM 0.15688465 0.15428884 0.05926460         NA         NA
#> 28 G3I18    2  3PLM 0.14382506 0.07697407 0.02111022         NA         NA
#> 29 G3I19    2  3PLM 0.11393846 0.14464817 0.04000362         NA         NA
#> 30 G3I20    2  3PLM 0.23541976 0.19741351 0.01684581         NA         NA
#> 31 G3I21    2  3PLM 0.57479199 0.06270311 0.01022201         NA         NA
#> 32 G3I22    2  3PLM 0.11463786 0.09814469 0.02978440         NA         NA
#> 33 G3I23    2  3PLM 0.16694784 0.06445158 0.01994484         NA         NA
#> 34 G3I24    2  3PLM 0.09417036 0.07931539 0.02063766         NA         NA
#> 35 G3I25    2  3PLM 0.14627533 0.07996861 0.01264067         NA         NA
#> 36 G3I26    2  3PLM 0.16594755 0.12414916 0.04434197         NA         NA
#> 37 G3I27    2  3PLM 0.10437934 0.20932778 0.05657496         NA         NA
#> 38 G3I28    5   GRM 0.04771674 0.05965450 0.05957007 0.07120741 0.09365933
#> 
#> 

# Extract the group-level parameter estimates (i.e., scale parameters)
getirt(fit.1, what = "group.par")
#> $Group1
#>           mu sigma2 sigma
#> estimates  0      1     1
#> se        NA     NA    NA
#> 
#> $Group2
#>                   mu    sigma2      sigma
#> estimates 0.48545465 0.5575252 0.74667607
#> se        0.01669618 0.0176349 0.01180894
#> 
#> $Group3
#>                    mu     sigma2      sigma
#> estimates -0.36960845 1.95177274 1.39705860
#> se         0.03123918 0.06173591 0.02209496
#> 

# Extract the posterior latent ability distributions for each group
getirt(fit.1, what = "weights")
#> $Group1
#>    theta       weight
#> 1  -6.00 2.596404e-24
#> 2  -5.75 6.661462e-12
#> 3  -5.50 6.687304e-11
#> 4  -5.25 3.455143e-10
#> 5  -5.00 1.240179e-09
#> 6  -4.75 3.547688e-09
#> 7  -4.50 8.892961e-09
#> 8  -4.25 2.132090e-08
#> 9  -4.00 5.399524e-08
#> 10 -3.75 1.635844e-07
#> 11 -3.50 6.878012e-07
#> 12 -3.25 4.538677e-06
#> 13 -3.00 4.536928e-05
#> 14 -2.75 4.969458e-04
#> 15 -2.50 3.578851e-03
#> 16 -2.25 1.233388e-02
#> 17 -2.00 2.221829e-02
#> 18 -1.75 2.756719e-02
#> 19 -1.50 3.067494e-02
#> 20 -1.25 3.699025e-02
#> 21 -1.00 5.137755e-02
#> 22 -0.75 7.176369e-02
#> 23 -0.50 8.884787e-02
#> 24 -0.25 1.006945e-01
#> 25  0.00 1.009843e-01
#> 26  0.25 9.093334e-02
#> 27  0.50 8.369774e-02
#> 28  0.75 8.097105e-02
#> 29  1.00 7.169684e-02
#> 30  1.25 5.091446e-02
#> 31  1.50 3.038547e-02
#> 32  1.75 1.784466e-02
#> 33  2.00 1.099557e-02
#> 34  2.25 6.559709e-03
#> 35  2.50 3.586983e-03
#> 36  2.75 1.890639e-03
#> 37  3.00 1.057103e-03
#> 38  3.25 6.683050e-04
#> 39  3.50 4.688643e-04
#> 40  3.75 3.327412e-04
#> 41  4.00 2.150724e-04
#> 42  4.25 1.183608e-04
#> 43  4.50 5.413395e-05
#> 44  4.75 2.065077e-05
#> 45  5.00 6.683841e-06
#> 46  5.25 1.876663e-06
#> 47  5.50 4.678457e-07
#> 48  5.75 1.058922e-07
#> 49  6.00 2.221172e-08
#> 
#> $Group2
#>    theta       weight
#> 1  -6.00 4.462701e-59
#> 2  -5.75 5.478200e-57
#> 3  -5.50 1.250041e-54
#> 4  -5.25 5.639732e-52
#> 5  -5.00 5.211160e-49
#> 6  -4.75 9.838327e-46
#> 7  -4.50 3.615813e-42
#> 8  -4.25 2.347660e-38
#> 9  -4.00 2.341231e-34
#> 10 -3.75 3.033389e-30
#> 11 -3.50 4.296681e-26
#> 12 -3.25 5.685725e-22
#> 13 -3.00 6.113539e-18
#> 14 -2.75 4.561112e-14
#> 15 -2.50 1.824133e-10
#> 16 -2.25 2.490398e-07
#> 17 -2.00 6.146622e-05
#> 18 -1.75 1.684078e-03
#> 19 -1.50 7.306548e-03
#> 20 -1.25 1.587374e-02
#> 21 -1.00 2.333185e-02
#> 22 -0.75 1.961168e-02
#> 23 -0.50 3.881370e-02
#> 24 -0.25 1.156738e-01
#> 25  0.00 8.052857e-02
#> 26  0.25 8.668950e-02
#> 27  0.50 1.901804e-01
#> 28  0.75 1.359231e-01
#> 29  1.00 8.624615e-02
#> 30  1.25 8.381817e-02
#> 31  1.50 5.935414e-02
#> 32  1.75 3.231143e-02
#> 33  2.00 1.436083e-02
#> 34  2.25 3.834694e-03
#> 35  2.50 1.089214e-03
#> 36  2.75 7.241950e-04
#> 37  3.00 9.441983e-04
#> 38  3.25 1.046937e-03
#> 39  3.50 5.018765e-04
#> 40  3.75 8.379032e-05
#> 41  4.00 5.539575e-06
#> 42  4.25 1.915245e-07
#> 43  4.50 4.643638e-09
#> 44  4.75 1.007048e-10
#> 45  5.00 2.322071e-12
#> 46  5.25 6.327049e-14
#> 47  5.50 2.144754e-15
#> 48  5.75 9.160692e-17
#> 49  6.00 4.870649e-18
#> 
#> $Group3
#>    theta       weight
#> 1  -6.00 7.731298e-07
#> 2  -5.75 3.145368e-06
#> 3  -5.50 1.179526e-05
#> 4  -5.25 4.054836e-05
#> 5  -5.00 1.268999e-04
#> 6  -4.75 3.584830e-04
#> 7  -4.50 9.049299e-04
#> 8  -4.25 2.018717e-03
#> 9  -4.00 3.937687e-03
#> 10 -3.75 6.669793e-03
#> 11 -3.50 9.833099e-03
#> 12 -3.25 1.285721e-02
#> 13 -3.00 1.552609e-02
#> 14 -2.75 1.830649e-02
#> 15 -2.50 2.198421e-02
#> 16 -2.25 2.658226e-02
#> 17 -2.00 3.040018e-02
#> 18 -1.75 3.230747e-02
#> 19 -1.50 3.635519e-02
#> 20 -1.25 5.077580e-02
#> 21 -1.00 7.272466e-02
#> 22 -0.75 6.898185e-02
#> 23 -0.50 5.513930e-02
#> 24 -0.25 6.620643e-02
#> 25  0.00 8.981610e-02
#> 26  0.25 8.120609e-02
#> 27  0.50 5.730756e-02
#> 28  0.75 4.971849e-02
#> 29  1.00 4.623885e-02
#> 30  1.25 3.634235e-02
#> 31  1.50 3.452349e-02
#> 32  1.75 3.537759e-02
#> 33  2.00 2.023739e-02
#> 34  2.25 7.344692e-03
#> 35  2.50 2.991667e-03
#> 36  2.75 1.765572e-03
#> 37  3.00 1.415897e-03
#> 38  3.25 1.272697e-03
#> 39  3.50 1.056410e-03
#> 40  3.75 7.149887e-04
#> 41  4.00 3.785780e-04
#> 42  4.25 1.595304e-04
#> 43  4.50 5.584698e-05
#> 44  4.75 1.699530e-05
#> 45  5.00 4.664197e-06
#> 46  5.25 1.184031e-06
#> 47  5.50 2.823534e-07
#> 48  5.75 6.377004e-08
#> 49  6.00 1.368772e-08
#> 

# 1-(2). Alternatively, the same parameter estimation can be performed by
# inserting a list of item metadata for the groups into the 'x' argument.
# If the item metadata contains item parameters to be used as starting values,
# set 'use.startval = TRUE'.
# Also, specify the groups in which the ability distribution scales
# will be freely estimated using their group names.
free.group <- group.name[2:3]
fit.2 <-
  est_mg(
    x = x, data = data, group.name = group.name, D = 1,
    free.group = free.group, use.gprior = TRUE,
    gprior = list(dist = "beta", params = c(5, 16)),
    group.mean = 0, group.var = 1, EmpHist = TRUE, use.startval = TRUE,
    Etol = 0.001, MaxE = 500
  )
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -159736.5073, Max-Change: 0.400714 EM iteration: 2, Loglike: -159610.2190, Max-Change: 0.165418 EM iteration: 3, Loglike: -159601.3164, Max-Change: 0.088375 EM iteration: 4, Loglike: -159596.5635, Max-Change: 0.054733 EM iteration: 5, Loglike: -159593.6124, Max-Change: 0.040346 EM iteration: 6, Loglike: -159591.6852, Max-Change: 0.035828 EM iteration: 7, Loglike: -159590.3654, Max-Change: 0.030999 EM iteration: 8, Loglike: -159589.4215, Max-Change: 0.026538 EM iteration: 9, Loglike: -159588.7193, Max-Change: 0.022633 EM iteration: 10, Loglike: -159588.1782, Max-Change: 0.019291 EM iteration: 11, Loglike: -159587.7475, Max-Change: 0.016548 EM iteration: 12, Loglike: -159587.3947, Max-Change: 0.014582 EM iteration: 13, Loglike: -159587.0984, Max-Change: 0.012898 EM iteration: 14, Loglike: -159586.8441, Max-Change: 0.011459 EM iteration: 15, Loglike: -159586.6218, Max-Change: 0.010231 EM iteration: 16, Loglike: -159586.4248, Max-Change: 0.009182 EM iteration: 17, Loglike: -159586.2479, Max-Change: 0.008283 EM iteration: 18, Loglike: -159586.0876, Max-Change: 0.007509 EM iteration: 19, Loglike: -159585.9412, Max-Change: 0.006838 EM iteration: 20, Loglike: -159585.8067, Max-Change: 0.006252 EM iteration: 21, Loglike: -159585.6824, Max-Change: 0.005735 EM iteration: 22, Loglike: -159585.5673, Max-Change: 0.005275 EM iteration: 23, Loglike: -159585.4607, Max-Change: 0.004862 EM iteration: 24, Loglike: -159585.3616, Max-Change: 0.004488 EM iteration: 25, Loglike: -159585.2689, Max-Change: 0.004146 EM iteration: 26, Loglike: -159585.1819, Max-Change: 0.003832 EM iteration: 27, Loglike: -159585.0999, Max-Change: 0.003542 EM iteration: 28, Loglike: -159585.0226, Max-Change: 0.003272 EM iteration: 29, Loglike: -159584.9497, Max-Change: 0.00302 EM iteration: 30, Loglike: -159584.8808, Max-Change: 0.002784 EM iteration: 31, Loglike: -159584.8156, Max-Change: 0.002563 EM iteration: 32, Loglike: -159584.7538, Max-Change: 0.002355 EM iteration: 33, Loglike: -159584.6951, Max-Change: 0.00216 EM iteration: 34, Loglike: -159584.6391, Max-Change: 0.001976 EM iteration: 35, Loglike: -159584.5858, Max-Change: 0.001804 EM iteration: 36, Loglike: -159584.5348, Max-Change: 0.001641 EM iteration: 37, Loglike: -159584.4860, Max-Change: 0.001489 EM iteration: 38, Loglike: -159584.4392, Max-Change: 0.001351 EM iteration: 39, Loglike: -159584.3942, Max-Change: 0.00126 EM iteration: 40, Loglike: -159584.3510, Max-Change: 0.001174 EM iteration: 41, Loglike: -159584.3092, Max-Change: 0.001092 EM iteration: 42, Loglike: -159584.2690, Max-Change: 0.001015 EM iteration: 43, Loglike: -159584.2300, Max-Change: 0.000942 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 8.77 seconds. 

# Summary of the estimation
summary(fit.2)
#> 
#> Call:
#> est_mg(x = x, data = data, group.name = group.name, D = 1, free.group = free.group, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), group.mean = 0, group.var = 1, EmpHist = TRUE, 
#>     use.startval = TRUE, Etol = 0.001, MaxE = 500)
#> 
#> Summary of the Data 
#>  Number of Items: 
#>   Overall: 116 unique items 
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of Cases: 
#>   Overall: 6000
#>   By group: 2000(Group1), 2000(Group2), 2000(Group3)
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 362
#>  Number of fixed items: 
#>   Overall: 0
#>   By group: 0(Group1), 0(Group2), 0(Group3)
#>  Number of E-step cycles completed: 43
#>  Maximum parameter change: 0.0009420884
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 8.11
#>  Standard error computation: 0.43
#>  Total computation: 8.77
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 
#>   Overall: 319168.5
#>   By group: 120345.77(Group1), 113945.748(Group2), 84876.941(Group3)
#> 
#>  Akaike Information Criterion (AIC): 319892.5
#>  Bayesian Information Criterion (BIC): 322317.7
#>  Item Parameters (Overall): 
#>         id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1     C1I1     2   3PLM   0.91  0.18   1.38  0.16   0.28  0.05     NA    NA
#> 2     C1I2     2   3PLM   2.14  0.15  -0.95  0.10   0.19  0.06     NA    NA
#> 3     C1I3     2   3PLM   1.07  0.12   0.63  0.13   0.18  0.04     NA    NA
#> 4     C1I4     2   3PLM   1.09  0.13  -0.14  0.21   0.29  0.07     NA    NA
#> 5     C1I5     2   3PLM   0.88  0.09  -0.12  0.22   0.18  0.07     NA    NA
#> 6     C1I6     2   3PLM   1.93  0.13   0.61  0.04   0.09  0.02     NA    NA
#> 7     C1I7     2   3PLM   1.12  0.14   1.12  0.09   0.15  0.03     NA    NA
#> 8     C1I8     2   3PLM   0.95  0.12   0.90  0.13   0.16  0.05     NA    NA
#> 9     C1I9     2   3PLM   0.91  0.13   0.66  0.18   0.21  0.06     NA    NA
#> 10   C1I10     2   3PLM   1.50  0.12   0.16  0.08   0.16  0.04     NA    NA
#> 11    G1I1     2   3PLM   0.97  0.11  -0.43  0.20   0.16  0.07     NA    NA
#> 12    G1I2     2   3PLM   0.90  0.14   1.21  0.13   0.11  0.04     NA    NA
#> 13    G1I3     2   3PLM   1.55  0.27   1.30  0.09   0.19  0.03     NA    NA
#> 14    G1I4     2   3PLM   1.58  0.22   0.30  0.12   0.30  0.04     NA    NA
#> 15    G1I5     2   3PLM   1.36  0.14  -0.15  0.13   0.16  0.05     NA    NA
#> 16    G1I6     2   3PLM   2.18  0.18   0.03  0.05   0.09  0.03     NA    NA
#> 17    G1I7     2   3PLM   1.48  0.16  -0.04  0.12   0.20  0.05     NA    NA
#> 18    G1I8     2   3PLM   2.57  0.47   1.18  0.07   0.32  0.02     NA    NA
#> 19    G1I9     2   3PLM   2.43  0.27  -0.90  0.11   0.25  0.06     NA    NA
#> 20   G1I10     2   3PLM   1.26  0.13  -1.73  0.23   0.24  0.10     NA    NA
#> 21   G1I11     2   3PLM   1.55  0.15  -1.08  0.17   0.21  0.08     NA    NA
#> 22   G1I12     2   3PLM   0.76  0.10  -0.78  0.32   0.20  0.08     NA    NA
#> 23   G1I13     2   3PLM   1.06  0.14  -0.10  0.20   0.21  0.07     NA    NA
#> 24   G1I14     2   3PLM   1.52  0.38   1.74  0.15   0.30  0.02     NA    NA
#> 25   G1I15     2   3PLM   0.87  0.10  -1.41  0.30   0.22  0.09     NA    NA
#> 26   G1I16     2   3PLM   1.05  0.11  -1.93  0.27   0.23  0.10     NA    NA
#> 27   G1I17     2   3PLM   1.07  0.13   0.28  0.15   0.16  0.05     NA    NA
#> 28   G1I18     2   3PLM   2.13  0.21  -0.06  0.08   0.24  0.04     NA    NA
#> 29   G1I19     2   3PLM   1.33  0.13  -1.36  0.20   0.21  0.09     NA    NA
#> 30   G1I20     2   3PLM   1.04  0.18   0.54  0.18   0.23  0.06     NA    NA
#> 31   G1I21     2   3PLM   0.93  0.14   0.78  0.15   0.13  0.05     NA    NA
#> 32   G1I22     2   3PLM   1.82  0.22  -0.59  0.15   0.37  0.06     NA    NA
#> 33   G1I23     2   3PLM   1.32  0.15  -1.13  0.23   0.26  0.09     NA    NA
#> 34   G1I24     2   3PLM   1.68  0.20   0.34  0.09   0.23  0.04     NA    NA
#> 35   G1I25     2   3PLM   1.62  0.18  -0.09  0.12   0.26  0.05     NA    NA
#> 36   G1I26     2   3PLM   1.94  0.26   0.67  0.08   0.26  0.03     NA    NA
#> 37   G1I27     2   3PLM   1.63  0.19  -1.50  0.21   0.28  0.10     NA    NA
#> 38   G1I28     2   3PLM   1.36  0.17   0.58  0.10   0.15  0.04     NA    NA
#> 39   G1I29     2   3PLM   0.93  0.09  -0.35  0.17   0.13  0.06     NA    NA
#> 40   G1I30     2   3PLM   1.06  0.29   2.24  0.25   0.17  0.03     NA    NA
#> 41   G1I31     2   3PLM   2.47  0.49   1.62  0.09   0.18  0.01     NA    NA
#> 42   G1I32     2   3PLM   1.12  0.12  -0.05  0.15   0.15  0.05     NA    NA
#> 43   G1I33     2   3PLM   1.64  0.18   0.20  0.09   0.16  0.04     NA    NA
#> 44   G1I34     2   3PLM   1.35  0.14   0.26  0.09   0.11  0.04     NA    NA
#> 45   G1I35     2   3PLM   1.36  0.17   1.29  0.08   0.07  0.02     NA    NA
#> 46   G1I36     2   3PLM   1.45  0.15  -1.20  0.19   0.23  0.09     NA    NA
#> 47   G1I37     2   3PLM   1.08  0.15  -0.54  0.27   0.28  0.09     NA    NA
#> 48   G1I38     5    GRM   1.07  0.06  -0.35  0.05   0.23  0.05   0.87  0.06
#> 49   C1I11     5    GRM   1.19  0.05  -2.18  0.10  -1.43  0.07  -0.72  0.05
#> 50   C1I12     5    GRM   0.92  0.04  -0.66  0.06   0.05  0.04   0.70  0.04
#> 51    G2I1     2   3PLM   1.79  0.20  -0.81  0.17   0.24  0.10     NA    NA
#> 52    G2I2     2   3PLM   0.86  0.12  -0.40  0.30   0.21  0.09     NA    NA
#> 53    G2I3     2   3PLM   1.11  0.14   0.12  0.20   0.18  0.08     NA    NA
#> 54    G2I4     2   3PLM   1.53  0.30   1.46  0.09   0.18  0.04     NA    NA
#> 55    G2I5     2   3PLM   0.73  0.11  -1.60  0.39   0.20  0.09     NA    NA
#> 56    G2I6     2   3PLM   1.11  0.15  -1.54  0.29   0.22  0.09     NA    NA
#> 57    G2I7     2   3PLM   1.35  0.14   0.29  0.12   0.13  0.05     NA    NA
#> 58    G2I8     2   3PLM   2.26  0.20  -0.05  0.08   0.14  0.05     NA    NA
#> 59    G2I9     2   3PLM   1.10  0.13  -1.50  0.27   0.21  0.09     NA    NA
#> 60   G2I10     2   3PLM   1.68  0.29   0.90  0.11   0.29  0.05     NA    NA
#> 61   G2I11     2   3PLM   0.98  0.13   0.88  0.14   0.11  0.05     NA    NA
#> 62   G2I12     2   3PLM   1.71  0.19  -0.64  0.18   0.26  0.10     NA    NA
#> 63   G2I13     2   3PLM   1.19  0.14  -1.25  0.24   0.21  0.09     NA    NA
#> 64   G2I14     2   3PLM   1.46  0.15   0.25  0.11   0.13  0.05     NA    NA
#> 65   G2I15     2   3PLM   1.60  0.21   0.08  0.17   0.27  0.08     NA    NA
#> 66   G2I16     2   3PLM   1.80  0.25   0.76  0.09   0.22  0.05     NA    NA
#> 67   G2I17     2   3PLM   2.14  0.22  -1.16  0.15   0.20  0.09     NA    NA
#> 68   G2I18     2   3PLM   1.82  0.21   0.65  0.08   0.14  0.04     NA    NA
#> 69   G2I19     2   3PLM   1.15  0.15   0.08  0.21   0.21  0.08     NA    NA
#> 70   G2I20     2   3PLM   1.72  0.48   2.16  0.18   0.16  0.02     NA    NA
#> 71   G2I21     2   3PLM   2.52  0.38   1.54  0.06   0.12  0.02     NA    NA
#> 72   G2I22     2   3PLM   1.52  0.19   0.27  0.15   0.22  0.07     NA    NA
#> 73   G2I23     2   3PLM   2.22  0.28   0.47  0.09   0.28  0.05     NA    NA
#> 74   G2I24     2   3PLM   1.46  0.12   0.38  0.07   0.07  0.03     NA    NA
#> 75   G2I25     2   3PLM   1.84  0.23   1.38  0.06   0.07  0.02     NA    NA
#> 76   G2I26     2   3PLM   1.93  0.20  -0.83  0.16   0.25  0.10     NA    NA
#> 77   G2I27     2   3PLM   1.01  0.12  -0.54  0.24   0.20  0.09     NA    NA
#> 78   G2I28     5    GRM   1.16  0.07  -0.33  0.07   0.17  0.05   0.80  0.05
#> 79    C2I1     2   3PLM   1.10  0.08  -0.26  0.11   0.12  0.04     NA    NA
#> 80    C2I2     2   3PLM   1.02  0.09   1.24  0.07   0.06  0.02     NA    NA
#> 81    C2I3     2   3PLM   1.66  0.16   1.38  0.05   0.13  0.01     NA    NA
#> 82    C2I4     2   3PLM   1.59  0.11   0.27  0.06   0.18  0.02     NA    NA
#> 83    C2I5     2   3PLM   1.46  0.10  -0.03  0.07   0.12  0.03     NA    NA
#> 84    C2I6     2   3PLM   2.08  0.11  -0.02  0.04   0.05  0.01     NA    NA
#> 85    C2I7     2   3PLM   1.59  0.11   0.03  0.07   0.18  0.03     NA    NA
#> 86    C2I8     2   3PLM   1.56  0.16   1.28  0.06   0.19  0.02     NA    NA
#> 87    C2I9     2   3PLM   2.22  0.14  -0.95  0.07   0.14  0.03     NA    NA
#> 88   C2I10     2   3PLM   1.61  0.12  -1.30  0.14   0.25  0.06     NA    NA
#> 89    G3I1     2   3PLM   1.56  0.15  -1.04  0.12   0.16  0.05     NA    NA
#> 90    G3I2     2   3PLM   0.79  0.10  -0.55  0.27   0.19  0.07     NA    NA
#> 91    G3I3     2   3PLM   1.18  0.13   0.11  0.12   0.17  0.04     NA    NA
#> 92    G3I4     2   3PLM   1.34  0.24   1.71  0.12   0.22  0.02     NA    NA
#> 93    G3I5     2   3PLM   0.81  0.11  -0.97  0.35   0.28  0.08     NA    NA
#> 94    G3I6     2   3PLM   1.20  0.15  -1.26  0.25   0.31  0.08     NA    NA
#> 95    G3I7     2   3PLM   1.40  0.14   0.29  0.08   0.14  0.03     NA    NA
#> 96    G3I8     2   3PLM   1.97  0.18  -0.04  0.06   0.15  0.02     NA    NA
#> 97    G3I9     2   3PLM   1.16  0.12  -1.14  0.20   0.20  0.07     NA    NA
#> 98   G3I10     2   3PLM   1.55  0.22   0.98  0.09   0.29  0.02     NA    NA
#> 99   G3I11     2   3PLM   0.93  0.10   0.81  0.11   0.08  0.03     NA    NA
#> 100  G3I12     2   3PLM   1.46  0.17  -0.72  0.15   0.26  0.05     NA    NA
#> 101  G3I13     2   3PLM   1.14  0.13  -1.04  0.21   0.23  0.07     NA    NA
#> 102  G3I14     2   3PLM   1.41  0.15   0.33  0.09   0.17  0.03     NA    NA
#> 103  G3I15     2   3PLM   1.32  0.14  -0.02  0.11   0.19  0.03     NA    NA
#> 104  G3I16     2   3PLM   1.62  0.18   0.76  0.08   0.20  0.02     NA    NA
#> 105  G3I17     2   3PLM   1.66  0.16  -1.47  0.15   0.18  0.06     NA    NA
#> 106  G3I18     2   3PLM   1.41  0.15   0.73  0.08   0.12  0.02     NA    NA
#> 107  G3I19     2   3PLM   1.02  0.12   0.12  0.14   0.16  0.04     NA    NA
#> 108  G3I20     2   3PLM   1.19  0.24   2.25  0.19   0.12  0.02     NA    NA
#> 109  G3I21     2   3PLM   3.09  0.57   1.63  0.06   0.14  0.01     NA    NA
#> 110  G3I22     2   3PLM   1.20  0.12   0.12  0.10   0.10  0.03     NA    NA
#> 111  G3I23     2   3PLM   1.84  0.17   0.35  0.06   0.15  0.02     NA    NA
#> 112  G3I24     2   3PLM   1.16  0.10   0.35  0.08   0.05  0.02     NA    NA
#> 113  G3I25     2   3PLM   1.42  0.15   1.35  0.08   0.05  0.01     NA    NA
#> 114  G3I26     2   3PLM   1.70  0.17  -0.94  0.12   0.21  0.05     NA    NA
#> 115  G3I27     2   3PLM   0.95  0.11  -0.57  0.20   0.18  0.06     NA    NA
#> 116  G3I28     5    GRM   0.97  0.05  -0.32  0.06   0.19  0.06   0.83  0.07
#>      par.5  se.5
#> 1       NA    NA
#> 2       NA    NA
#> 3       NA    NA
#> 4       NA    NA
#> 5       NA    NA
#> 6       NA    NA
#> 7       NA    NA
#> 8       NA    NA
#> 9       NA    NA
#> 10      NA    NA
#> 11      NA    NA
#> 12      NA    NA
#> 13      NA    NA
#> 14      NA    NA
#> 15      NA    NA
#> 16      NA    NA
#> 17      NA    NA
#> 18      NA    NA
#> 19      NA    NA
#> 20      NA    NA
#> 21      NA    NA
#> 22      NA    NA
#> 23      NA    NA
#> 24      NA    NA
#> 25      NA    NA
#> 26      NA    NA
#> 27      NA    NA
#> 28      NA    NA
#> 29      NA    NA
#> 30      NA    NA
#> 31      NA    NA
#> 32      NA    NA
#> 33      NA    NA
#> 34      NA    NA
#> 35      NA    NA
#> 36      NA    NA
#> 37      NA    NA
#> 38      NA    NA
#> 39      NA    NA
#> 40      NA    NA
#> 41      NA    NA
#> 42      NA    NA
#> 43      NA    NA
#> 44      NA    NA
#> 45      NA    NA
#> 46      NA    NA
#> 47      NA    NA
#> 48    1.43  0.08
#> 49   -0.10  0.03
#> 50    1.26  0.06
#> 51      NA    NA
#> 52      NA    NA
#> 53      NA    NA
#> 54      NA    NA
#> 55      NA    NA
#> 56      NA    NA
#> 57      NA    NA
#> 58      NA    NA
#> 59      NA    NA
#> 60      NA    NA
#> 61      NA    NA
#> 62      NA    NA
#> 63      NA    NA
#> 64      NA    NA
#> 65      NA    NA
#> 66      NA    NA
#> 67      NA    NA
#> 68      NA    NA
#> 69      NA    NA
#> 70      NA    NA
#> 71      NA    NA
#> 72      NA    NA
#> 73      NA    NA
#> 74      NA    NA
#> 75      NA    NA
#> 76      NA    NA
#> 77      NA    NA
#> 78    1.45  0.07
#> 79      NA    NA
#> 80      NA    NA
#> 81      NA    NA
#> 82      NA    NA
#> 83      NA    NA
#> 84      NA    NA
#> 85      NA    NA
#> 86      NA    NA
#> 87      NA    NA
#> 88      NA    NA
#> 89      NA    NA
#> 90      NA    NA
#> 91      NA    NA
#> 92      NA    NA
#> 93      NA    NA
#> 94      NA    NA
#> 95      NA    NA
#> 96      NA    NA
#> 97      NA    NA
#> 98      NA    NA
#> 99      NA    NA
#> 100     NA    NA
#> 101     NA    NA
#> 102     NA    NA
#> 103     NA    NA
#> 104     NA    NA
#> 105     NA    NA
#> 106     NA    NA
#> 107     NA    NA
#> 108     NA    NA
#> 109     NA    NA
#> 110     NA    NA
#> 111     NA    NA
#> 112     NA    NA
#> 113     NA    NA
#> 114     NA    NA
#> 115     NA    NA
#> 116   1.55  0.09
#>  Group Parameters: 
#>                      mu  sigma2  sigma
#> estimate(Group1)   0.00    1.00   1.00
#> se(Group1)           NA      NA     NA
#> estimate(Group2)   0.51    0.54   0.73
#> se(Group2)         0.02    0.02   0.01
#> estimate(Group3)  -0.31    1.80   1.34
#> se(Group3)         0.03    0.06   0.02
#> 

## ------------------------------------------------------------------------------
# 2. MG calibration with FIPC using simMG data
#  - Details:
#    (a) Fix the parameters of the common items between the groups
#        (i.e., items C1I1–C1I12 between Groups 1 and 2, and
#        items C2I1–C2I10 between Groups 2 and 3)
#    (b) Freely estimate the means and variances of the ability
#        distributions for all three groups
## ------------------------------------------------------------------------------
# 2-(1). Freely estimate the means and variances for all three groups
# Set all three groups as free groups in which the scales
# of the ability distributions are to be freely estimated
free.group <- 1:3 # or use 'free.group <- group.name'

# Specify the locations of items to be fixed in each group's metadata
# For Group 1: C1I1–C1I12 are located in rows 1–10 and 49–50
# For Group 2: C1I1–C1I12 are in rows 1–12, and
#              C2I1–C2I10 are in rows 41–50
# For Group 3: C2I1–C2I10 are in rows 1–10
fix.loc <- list(
  c(1:10, 49:50),
  c(1:12, 41:50),
  c(1:10)
)

# Estimate IRT parameters using MG-FIPC:
# When FIPC is implemented, item metadata for all groups
# must be provided via the 'x' argument.
# For fixed items, their item parameters must be specified
# in the metadata. For non-fixed items, any placeholder values
# can be used in the metadata.
# Also set fipc = TRUE and fipc.method = "MEM"
fit.3 <-
  est_mg(
    x = x, data = data, group.name = group.name, D = 1,
    free.group = free.group, use.gprior = TRUE,
    gprior = list(dist = "beta", params = c(5, 16)),
    EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
    fipc.method = "MEM", fix.loc = fix.loc
  )
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -56756.7548, Max-Change: 3.554166 EM iteration: 2, Loglike: -160424.0878, Max-Change: 0.780368 EM iteration: 3, Loglike: -159755.3559, Max-Change: 0.276071 EM iteration: 4, Loglike: -159686.2158, Max-Change: 0.151073 EM iteration: 5, Loglike: -159657.5520, Max-Change: 0.102849 EM iteration: 6, Loglike: -159642.1127, Max-Change: 0.080257 EM iteration: 7, Loglike: -159633.1411, Max-Change: 0.066343 EM iteration: 8, Loglike: -159627.6844, Max-Change: 0.055853 EM iteration: 9, Loglike: -159624.2459, Max-Change: 0.047172 EM iteration: 10, Loglike: -159622.0120, Max-Change: 0.039781 EM iteration: 11, Loglike: -159620.5195, Max-Change: 0.033468 EM iteration: 12, Loglike: -159619.4950, Max-Change: 0.028101 EM iteration: 13, Loglike: -159618.7723, Max-Change: 0.023565 EM iteration: 14, Loglike: -159618.2481, Max-Change: 0.019752 EM iteration: 15, Loglike: -159617.8570, Max-Change: 0.016845 EM iteration: 16, Loglike: -159617.5566, Max-Change: 0.014864 EM iteration: 17, Loglike: -159617.3191, Max-Change: 0.01318 EM iteration: 18, Loglike: -159617.1261, Max-Change: 0.011751 EM iteration: 19, Loglike: -159616.9650, Max-Change: 0.010538 EM iteration: 20, Loglike: -159616.8273, Max-Change: 0.009507 EM iteration: 21, Loglike: -159616.7072, Max-Change: 0.008626 EM iteration: 22, Loglike: -159616.6005, Max-Change: 0.007869 EM iteration: 23, Loglike: -159616.5042, Max-Change: 0.007214 EM iteration: 24, Loglike: -159616.4163, Max-Change: 0.006643 EM iteration: 25, Loglike: -159616.3352, Max-Change: 0.006141 EM iteration: 26, Loglike: -159616.2596, Max-Change: 0.005695 EM iteration: 27, Loglike: -159616.1888, Max-Change: 0.005295 EM iteration: 28, Loglike: -159616.1221, Max-Change: 0.004934 EM iteration: 29, Loglike: -159616.0589, Max-Change: 0.004606 EM iteration: 30, Loglike: -159615.9989, Max-Change: 0.004306 EM iteration: 31, Loglike: -159615.9416, Max-Change: 0.00403 EM iteration: 32, Loglike: -159615.8868, Max-Change: 0.003774 EM iteration: 33, Loglike: -159615.8343, Max-Change: 0.003537 EM iteration: 34, Loglike: -159615.7838, Max-Change: 0.003317 EM iteration: 35, Loglike: -159615.7352, Max-Change: 0.003111 EM iteration: 36, Loglike: -159615.6884, Max-Change: 0.002919 EM iteration: 37, Loglike: -159615.6431, Max-Change: 0.00274 EM iteration: 38, Loglike: -159615.5994, Max-Change: 0.002572 EM iteration: 39, Loglike: -159615.5570, Max-Change: 0.002416 EM iteration: 40, Loglike: -159615.5159, Max-Change: 0.002269 EM iteration: 41, Loglike: -159615.4761, Max-Change: 0.002131 EM iteration: 42, Loglike: -159615.4373, Max-Change: 0.002003 EM iteration: 43, Loglike: -159615.3997, Max-Change: 0.001882 EM iteration: 44, Loglike: -159615.3630, Max-Change: 0.00177 EM iteration: 45, Loglike: -159615.3274, Max-Change: 0.001665 EM iteration: 46, Loglike: -159615.2926, Max-Change: 0.001566 EM iteration: 47, Loglike: -159615.2587, Max-Change: 0.001474 EM iteration: 48, Loglike: -159615.2256, Max-Change: 0.001388 EM iteration: 49, Loglike: -159615.1933, Max-Change: 0.001308 EM iteration: 50, Loglike: -159615.1617, Max-Change: 0.001232 EM iteration: 51, Loglike: -159615.1309, Max-Change: 0.001162 EM iteration: 52, Loglike: -159615.1007, Max-Change: 0.001097 EM iteration: 53, Loglike: -159615.0711, Max-Change: 0.001035 EM iteration: 54, Loglike: -159615.0422, Max-Change: 0.000978 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 8.91 seconds. 

# Summary of the estimation
summary(fit.3)
#> 
#> Call:
#> est_mg(x = x, data = data, group.name = group.name, D = 1, free.group = free.group, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE, 
#>     fipc.method = "MEM", fix.loc = fix.loc)
#> 
#> Summary of the Data 
#>  Number of Items: 
#>   Overall: 116 unique items 
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of Cases: 
#>   Overall: 6000
#>   By group: 2000(Group1), 2000(Group2), 2000(Group3)
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 294
#>  Number of fixed items: 
#>   Overall: 22
#>   By group: 12(Group1), 22(Group2), 10(Group3)
#>  Number of E-step cycles completed: 54
#>  Maximum parameter change: 0.0009779691
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 8.22
#>  Standard error computation: 0.15
#>  Total computation: 8.91
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 
#>   Overall: 319230.1
#>   By group: 120354.311(Group1), 113985.392(Group2), 84890.382(Group3)
#> 
#>  Akaike Information Criterion (AIC): 319818.1
#>  Bayesian Information Criterion (BIC): 321787.7
#>  Item Parameters (Overall): 
#>         id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1     C1I1     2   3PLM   0.76    NA   1.46    NA   0.26    NA     NA    NA
#> 2     C1I2     2   3PLM   1.92    NA  -1.05    NA   0.18    NA     NA    NA
#> 3     C1I3     2   3PLM   0.93    NA   0.39    NA   0.10    NA     NA    NA
#> 4     C1I4     2   3PLM   1.05    NA  -0.41    NA   0.20    NA     NA    NA
#> 5     C1I5     2   3PLM   0.87    NA  -0.12    NA   0.16    NA     NA    NA
#> 6     C1I6     2   3PLM   1.70    NA   0.63    NA   0.07    NA     NA    NA
#> 7     C1I7     2   3PLM   0.91    NA   1.02    NA   0.12    NA     NA    NA
#> 8     C1I8     2   3PLM   0.84    NA   0.80    NA   0.11    NA     NA    NA
#> 9     C1I9     2   3PLM   0.85    NA   0.85    NA   0.26    NA     NA    NA
#> 10   C1I10     2   3PLM   1.53    NA   0.09    NA   0.14    NA     NA    NA
#> 11    G1I1     2   3PLM   0.96  0.10  -0.46  0.20   0.16  0.07     NA    NA
#> 12    G1I2     2   3PLM   0.88  0.13   1.21  0.14   0.10  0.04     NA    NA
#> 13    G1I3     2   3PLM   1.49  0.26   1.31  0.09   0.18  0.03     NA    NA
#> 14    G1I4     2   3PLM   1.52  0.21   0.27  0.12   0.30  0.04     NA    NA
#> 15    G1I5     2   3PLM   1.33  0.14  -0.18  0.13   0.16  0.05     NA    NA
#> 16    G1I6     2   3PLM   2.15  0.18   0.01  0.05   0.08  0.03     NA    NA
#> 17    G1I7     2   3PLM   1.45  0.16  -0.07  0.12   0.20  0.05     NA    NA
#> 18    G1I8     2   3PLM   2.50  0.46   1.18  0.07   0.32  0.02     NA    NA
#> 19    G1I9     2   3PLM   2.40  0.26  -0.93  0.11   0.24  0.06     NA    NA
#> 20   G1I10     2   3PLM   1.26  0.13  -1.75  0.23   0.24  0.10     NA    NA
#> 21   G1I11     2   3PLM   1.55  0.15  -1.11  0.16   0.21  0.08     NA    NA
#> 22   G1I12     2   3PLM   0.75  0.09  -0.81  0.31   0.19  0.08     NA    NA
#> 23   G1I13     2   3PLM   1.04  0.14  -0.13  0.21   0.20  0.07     NA    NA
#> 24   G1I14     2   3PLM   1.46  0.37   1.77  0.15   0.30  0.03     NA    NA
#> 25   G1I15     2   3PLM   0.86  0.09  -1.44  0.29   0.21  0.09     NA    NA
#> 26   G1I16     2   3PLM   1.05  0.11  -1.95  0.27   0.23  0.10     NA    NA
#> 27   G1I17     2   3PLM   1.04  0.13   0.25  0.15   0.15  0.05     NA    NA
#> 28   G1I18     2   3PLM   2.08  0.20  -0.09  0.08   0.24  0.04     NA    NA
#> 29   G1I19     2   3PLM   1.32  0.13  -1.38  0.19   0.20  0.08     NA    NA
#> 30   G1I20     2   3PLM   1.01  0.17   0.51  0.19   0.22  0.06     NA    NA
#> 31   G1I21     2   3PLM   0.92  0.13   0.77  0.15   0.13  0.05     NA    NA
#> 32   G1I22     2   3PLM   1.77  0.21  -0.65  0.16   0.35  0.06     NA    NA
#> 33   G1I23     2   3PLM   1.31  0.15  -1.17  0.23   0.25  0.09     NA    NA
#> 34   G1I24     2   3PLM   1.64  0.19   0.32  0.10   0.23  0.04     NA    NA
#> 35   G1I25     2   3PLM   1.58  0.18  -0.13  0.12   0.25  0.05     NA    NA
#> 36   G1I26     2   3PLM   1.88  0.25   0.66  0.08   0.25  0.03     NA    NA
#> 37   G1I27     2   3PLM   1.63  0.18  -1.54  0.20   0.27  0.10     NA    NA
#> 38   G1I28     2   3PLM   1.33  0.16   0.57  0.10   0.15  0.04     NA    NA
#> 39   G1I29     2   3PLM   0.92  0.09  -0.37  0.17   0.12  0.05     NA    NA
#> 40   G1I30     2   3PLM   1.00  0.27   2.30  0.25   0.17  0.03     NA    NA
#> 41   G1I31     2   3PLM   2.39  0.46   1.64  0.09   0.18  0.01     NA    NA
#> 42   G1I32     2   3PLM   1.10  0.12  -0.08  0.15   0.15  0.05     NA    NA
#> 43   G1I33     2   3PLM   1.60  0.17   0.18  0.09   0.16  0.04     NA    NA
#> 44   G1I34     2   3PLM   1.33  0.13   0.24  0.09   0.11  0.04     NA    NA
#> 45   G1I35     2   3PLM   1.33  0.17   1.30  0.08   0.07  0.02     NA    NA
#> 46   G1I36     2   3PLM   1.44  0.14  -1.23  0.19   0.22  0.09     NA    NA
#> 47   G1I37     2   3PLM   1.05  0.14  -0.60  0.27   0.27  0.09     NA    NA
#> 48   G1I38     5    GRM   1.06  0.06  -0.37  0.05   0.22  0.05   0.87  0.06
#> 49   C1I11     5    GRM   1.23    NA  -2.08    NA  -1.35    NA  -0.71    NA
#> 50   C1I12     5    GRM   0.88    NA  -0.76    NA  -0.01    NA   0.67    NA
#> 51    G2I1     2   3PLM   1.74  0.19  -0.87  0.17   0.23  0.10     NA    NA
#> 52    G2I2     2   3PLM   0.83  0.11  -0.45  0.31   0.21  0.09     NA    NA
#> 53    G2I3     2   3PLM   1.06  0.13   0.08  0.20   0.18  0.07     NA    NA
#> 54    G2I4     2   3PLM   1.45  0.28   1.50  0.10   0.18  0.04     NA    NA
#> 55    G2I5     2   3PLM   0.70  0.10  -1.68  0.40   0.20  0.09     NA    NA
#> 56    G2I6     2   3PLM   1.08  0.14  -1.60  0.29   0.22  0.09     NA    NA
#> 57    G2I7     2   3PLM   1.30  0.13   0.27  0.12   0.13  0.05     NA    NA
#> 58    G2I8     2   3PLM   2.18  0.19  -0.09  0.08   0.14  0.05     NA    NA
#> 59    G2I9     2   3PLM   1.07  0.13  -1.57  0.28   0.21  0.09     NA    NA
#> 60   G2I10     2   3PLM   1.59  0.27   0.90  0.12   0.29  0.05     NA    NA
#> 61   G2I11     2   3PLM   0.93  0.12   0.88  0.15   0.11  0.05     NA    NA
#> 62   G2I12     2   3PLM   1.66  0.18  -0.70  0.18   0.25  0.10     NA    NA
#> 63   G2I13     2   3PLM   1.16  0.13  -1.31  0.24   0.21  0.09     NA    NA
#> 64   G2I14     2   3PLM   1.39  0.14   0.22  0.11   0.13  0.05     NA    NA
#> 65   G2I15     2   3PLM   1.53  0.20   0.03  0.17   0.26  0.08     NA    NA
#> 66   G2I16     2   3PLM   1.70  0.23   0.75  0.10   0.21  0.05     NA    NA
#> 67   G2I17     2   3PLM   2.11  0.22  -1.20  0.14   0.20  0.09     NA    NA
#> 68   G2I18     2   3PLM   1.72  0.20   0.64  0.08   0.14  0.04     NA    NA
#> 69   G2I19     2   3PLM   1.10  0.14   0.05  0.22   0.21  0.08     NA    NA
#> 70   G2I20     2   3PLM   1.62  0.45   2.24  0.19   0.16  0.02     NA    NA
#> 71   G2I21     2   3PLM   2.39  0.36   1.59  0.06   0.12  0.02     NA    NA
#> 72   G2I22     2   3PLM   1.43  0.17   0.22  0.15   0.20  0.07     NA    NA
#> 73   G2I23     2   3PLM   2.09  0.26   0.44  0.10   0.28  0.05     NA    NA
#> 74   G2I24     2   3PLM   1.40  0.12   0.37  0.08   0.07  0.03     NA    NA
#> 75   G2I25     2   3PLM   1.74  0.22   1.42  0.06   0.07  0.02     NA    NA
#> 76   G2I26     2   3PLM   1.89  0.20  -0.88  0.16   0.24  0.10     NA    NA
#> 77   G2I27     2   3PLM   0.98  0.12  -0.59  0.25   0.20  0.09     NA    NA
#> 78   G2I28     5    GRM   1.12  0.07  -0.37  0.07   0.15  0.05   0.80  0.05
#> 79    C2I1     2   3PLM   0.97    NA  -0.46    NA   0.05    NA     NA    NA
#> 80    C2I2     2   3PLM   0.85    NA   1.18    NA   0.01    NA     NA    NA
#> 81    C2I3     2   3PLM   1.43    NA   1.41    NA   0.10    NA     NA    NA
#> 82    C2I4     2   3PLM   1.48    NA   0.18    NA   0.17    NA     NA    NA
#> 83    C2I5     2   3PLM   1.27    NA  -0.23    NA   0.03    NA     NA    NA
#> 84    C2I6     2   3PLM   2.02    NA  -0.09    NA   0.01    NA     NA    NA
#> 85    C2I7     2   3PLM   1.37    NA  -0.13    NA   0.10    NA     NA    NA
#> 86    C2I8     2   3PLM   1.67    NA   1.25    NA   0.19    NA     NA    NA
#> 87    C2I9     2   3PLM   2.28    NA  -1.01    NA   0.10    NA     NA    NA
#> 88   C2I10     2   3PLM   1.42    NA  -1.65    NA   0.11    NA     NA    NA
#> 89    G3I1     2   3PLM   1.60  0.14  -1.04  0.12   0.14  0.05     NA    NA
#> 90    G3I2     2   3PLM   0.79  0.09  -0.62  0.25   0.16  0.07     NA    NA
#> 91    G3I3     2   3PLM   1.16  0.12   0.07  0.12   0.16  0.04     NA    NA
#> 92    G3I4     2   3PLM   1.30  0.24   1.72  0.12   0.22  0.02     NA    NA
#> 93    G3I5     2   3PLM   0.80  0.10  -1.09  0.33   0.24  0.09     NA    NA
#> 94    G3I6     2   3PLM   1.18  0.13  -1.38  0.24   0.26  0.08     NA    NA
#> 95    G3I7     2   3PLM   1.37  0.13   0.26  0.09   0.13  0.03     NA    NA
#> 96    G3I8     2   3PLM   1.95  0.18  -0.07  0.06   0.14  0.02     NA    NA
#> 97    G3I9     2   3PLM   1.18  0.11  -1.16  0.19   0.18  0.07     NA    NA
#> 98   G3I10     2   3PLM   1.50  0.21   0.98  0.10   0.29  0.02     NA    NA
#> 99   G3I11     2   3PLM   0.91  0.10   0.80  0.11   0.08  0.03     NA    NA
#> 100  G3I12     2   3PLM   1.44  0.17  -0.78  0.15   0.23  0.06     NA    NA
#> 101  G3I13     2   3PLM   1.15  0.13  -1.06  0.20   0.21  0.07     NA    NA
#> 102  G3I14     2   3PLM   1.38  0.14   0.31  0.09   0.16  0.03     NA    NA
#> 103  G3I15     2   3PLM   1.29  0.13  -0.06  0.11   0.18  0.04     NA    NA
#> 104  G3I16     2   3PLM   1.58  0.18   0.75  0.08   0.20  0.02     NA    NA
#> 105  G3I17     2   3PLM   1.72  0.15  -1.46  0.13   0.16  0.06     NA    NA
#> 106  G3I18     2   3PLM   1.37  0.14   0.71  0.08   0.11  0.02     NA    NA
#> 107  G3I19     2   3PLM   0.99  0.11   0.07  0.15   0.15  0.04     NA    NA
#> 108  G3I20     2   3PLM   1.14  0.23   2.30  0.20   0.12  0.02     NA    NA
#> 109  G3I21     2   3PLM   2.97  0.55   1.66  0.06   0.14  0.01     NA    NA
#> 110  G3I22     2   3PLM   1.19  0.11   0.10  0.10   0.10  0.03     NA    NA
#> 111  G3I23     2   3PLM   1.79  0.17   0.33  0.06   0.15  0.02     NA    NA
#> 112  G3I24     2   3PLM   1.15  0.09   0.34  0.08   0.05  0.02     NA    NA
#> 113  G3I25     2   3PLM   1.38  0.15   1.36  0.08   0.05  0.01     NA    NA
#> 114  G3I26     2   3PLM   1.72  0.17  -0.96  0.12   0.20  0.05     NA    NA
#> 115  G3I27     2   3PLM   0.95  0.10  -0.61  0.20   0.16  0.06     NA    NA
#> 116  G3I28     5    GRM   1.00  0.05  -0.31  0.05   0.19  0.05   0.82  0.07
#>      par.5  se.5
#> 1       NA    NA
#> 2       NA    NA
#> 3       NA    NA
#> 4       NA    NA
#> 5       NA    NA
#> 6       NA    NA
#> 7       NA    NA
#> 8       NA    NA
#> 9       NA    NA
#> 10      NA    NA
#> 11      NA    NA
#> 12      NA    NA
#> 13      NA    NA
#> 14      NA    NA
#> 15      NA    NA
#> 16      NA    NA
#> 17      NA    NA
#> 18      NA    NA
#> 19      NA    NA
#> 20      NA    NA
#> 21      NA    NA
#> 22      NA    NA
#> 23      NA    NA
#> 24      NA    NA
#> 25      NA    NA
#> 26      NA    NA
#> 27      NA    NA
#> 28      NA    NA
#> 29      NA    NA
#> 30      NA    NA
#> 31      NA    NA
#> 32      NA    NA
#> 33      NA    NA
#> 34      NA    NA
#> 35      NA    NA
#> 36      NA    NA
#> 37      NA    NA
#> 38      NA    NA
#> 39      NA    NA
#> 40      NA    NA
#> 41      NA    NA
#> 42      NA    NA
#> 43      NA    NA
#> 44      NA    NA
#> 45      NA    NA
#> 46      NA    NA
#> 47      NA    NA
#> 48    1.44  0.08
#> 49   -0.12    NA
#> 50    1.25    NA
#> 51      NA    NA
#> 52      NA    NA
#> 53      NA    NA
#> 54      NA    NA
#> 55      NA    NA
#> 56      NA    NA
#> 57      NA    NA
#> 58      NA    NA
#> 59      NA    NA
#> 60      NA    NA
#> 61      NA    NA
#> 62      NA    NA
#> 63      NA    NA
#> 64      NA    NA
#> 65      NA    NA
#> 66      NA    NA
#> 67      NA    NA
#> 68      NA    NA
#> 69      NA    NA
#> 70      NA    NA
#> 71      NA    NA
#> 72      NA    NA
#> 73      NA    NA
#> 74      NA    NA
#> 75      NA    NA
#> 76      NA    NA
#> 77      NA    NA
#> 78    1.48  0.07
#> 79      NA    NA
#> 80      NA    NA
#> 81      NA    NA
#> 82      NA    NA
#> 83      NA    NA
#> 84      NA    NA
#> 85      NA    NA
#> 86      NA    NA
#> 87      NA    NA
#> 88      NA    NA
#> 89      NA    NA
#> 90      NA    NA
#> 91      NA    NA
#> 92      NA    NA
#> 93      NA    NA
#> 94      NA    NA
#> 95      NA    NA
#> 96      NA    NA
#> 97      NA    NA
#> 98      NA    NA
#> 99      NA    NA
#> 100     NA    NA
#> 101     NA    NA
#> 102     NA    NA
#> 103     NA    NA
#> 104     NA    NA
#> 105     NA    NA
#> 106     NA    NA
#> 107     NA    NA
#> 108     NA    NA
#> 109     NA    NA
#> 110     NA    NA
#> 111     NA    NA
#> 112     NA    NA
#> 113     NA    NA
#> 114     NA    NA
#> 115     NA    NA
#> 116   1.52  0.09
#>  Group Parameters: 
#>                      mu  sigma2  sigma
#> estimate(Group1)  -0.01    1.01   1.01
#> se(Group1)         0.02    0.03   0.02
#> estimate(Group2)   0.50    0.58   0.76
#> se(Group2)         0.02    0.02   0.01
#> estimate(Group3)  -0.28    1.66   1.29
#> se(Group3)         0.03    0.05   0.02
#> 

# Extract the item parameter estimates
getirt(fit.3, what = "par.est")
#> $overall
#>        id cats model     par.1       par.2       par.3      par.4     par.5
#> 1    C1I1    2  3PLM 0.7600000  1.46000000  0.26000000         NA        NA
#> 2    C1I2    2  3PLM 1.9200000 -1.05000000  0.18000000         NA        NA
#> 3    C1I3    2  3PLM 0.9300000  0.39000000  0.10000000         NA        NA
#> 4    C1I4    2  3PLM 1.0500000 -0.41000000  0.20000000         NA        NA
#> 5    C1I5    2  3PLM 0.8700000 -0.12000000  0.16000000         NA        NA
#> 6    C1I6    2  3PLM 1.7000000  0.63000000  0.07000000         NA        NA
#> 7    C1I7    2  3PLM 0.9100000  1.02000000  0.12000000         NA        NA
#> 8    C1I8    2  3PLM 0.8400000  0.80000000  0.11000000         NA        NA
#> 9    C1I9    2  3PLM 0.8500000  0.85000000  0.26000000         NA        NA
#> 10  C1I10    2  3PLM 1.5300000  0.09000000  0.14000000         NA        NA
#> 11   G1I1    2  3PLM 0.9601094 -0.46323480  0.15743926         NA        NA
#> 12   G1I2    2  3PLM 0.8792351  1.21120501  0.10365721         NA        NA
#> 13   G1I3    2  3PLM 1.4922037  1.31254694  0.18433436         NA        NA
#> 14   G1I4    2  3PLM 1.5246246  0.26709705  0.29751014         NA        NA
#> 15   G1I5    2  3PLM 1.3343885 -0.18476888  0.15642097         NA        NA
#> 16   G1I6    2  3PLM 2.1457296  0.01122821  0.08373745         NA        NA
#> 17   G1I7    2  3PLM 1.4492421 -0.07232467  0.19831474         NA        NA
#> 18   G1I8    2  3PLM 2.4981366  1.18142115  0.32459740         NA        NA
#> 19   G1I9    2  3PLM 2.4036619 -0.93389393  0.24090792         NA        NA
#> 20  G1I10    2  3PLM 1.2603565 -1.75076139  0.23712596         NA        NA
#> 21  G1I11    2  3PLM 1.5482444 -1.11447546  0.20564573         NA        NA
#> 22  G1I12    2  3PLM 0.7501881 -0.81274221  0.19483606         NA        NA
#> 23  G1I13    2  3PLM 1.0391850 -0.13129452  0.20305864         NA        NA
#> 24  G1I14    2  3PLM 1.4641086  1.76817734  0.29587746         NA        NA
#> 25  G1I15    2  3PLM 0.8642507 -1.44005202  0.21246953         NA        NA
#> 26  G1I16    2  3PLM 1.0481967 -1.94936219  0.22996558         NA        NA
#> 27  G1I17    2  3PLM 1.0413440  0.24963719  0.15297903         NA        NA
#> 28  G1I18    2  3PLM 2.0821164 -0.08782853  0.23783714         NA        NA
#> 29  G1I19    2  3PLM 1.3246385 -1.38457105  0.20248019         NA        NA
#> 30  G1I20    2  3PLM 1.0078270  0.50901778  0.22362359         NA        NA
#> 31  G1I21    2  3PLM 0.9153677  0.76786397  0.13229024         NA        NA
#> 32  G1I22    2  3PLM 1.7709516 -0.64567759  0.35461144         NA        NA
#> 33  G1I23    2  3PLM 1.3100376 -1.16619310  0.25418714         NA        NA
#> 34  G1I24    2  3PLM 1.6392138  0.31576653  0.22824910         NA        NA
#> 35  G1I25    2  3PLM 1.5820253 -0.12733170  0.25023679         NA        NA
#> 36  G1I26    2  3PLM 1.8799489  0.66003503  0.25368918         NA        NA
#> 37  G1I27    2  3PLM 1.6294671 -1.53964203  0.27165399         NA        NA
#> 38  G1I28    2  3PLM 1.3331814  0.56546916  0.14913136         NA        NA
#> 39  G1I29    2  3PLM 0.9154613 -0.37004169  0.12429485         NA        NA
#> 40  G1I30    2  3PLM 1.0026860  2.29722730  0.16842602         NA        NA
#> 41  G1I31    2  3PLM 2.3881153  1.64223591  0.18163836         NA        NA
#> 42  G1I32    2  3PLM 1.1031428 -0.07825039  0.14531435         NA        NA
#> 43  G1I33    2  3PLM 1.5998724  0.17672610  0.15654349         NA        NA
#> 44  G1I34    2  3PLM 1.3266408  0.24185668  0.10791955         NA        NA
#> 45  G1I35    2  3PLM 1.3307122  1.29621629  0.06712955         NA        NA
#> 46  G1I36    2  3PLM 1.4375780 -1.23256008  0.21920779         NA        NA
#> 47  G1I37    2  3PLM 1.0523257 -0.59860569  0.26551312         NA        NA
#> 48  G1I38    5   GRM 1.0557817 -0.36687770  0.21653854  0.8654794  1.435633
#> 49  C1I11    5   GRM 1.2300000 -2.08000000 -1.35000000 -0.7100000 -0.120000
#> 50  C1I12    5   GRM 0.8800000 -0.76000000 -0.01000000  0.6700000  1.250000
#> 51   G2I1    2  3PLM 1.7426027 -0.86693468  0.23241989         NA        NA
#> 52   G2I2    2  3PLM 0.8271656 -0.44864390  0.21072500         NA        NA
#> 53   G2I3    2  3PLM 1.0594812  0.08411544  0.17868700         NA        NA
#> 54   G2I4    2  3PLM 1.4496890  1.49538158  0.17845021         NA        NA
#> 55   G2I5    2  3PLM 0.7000088 -1.68320349  0.19753124         NA        NA
#> 56   G2I6    2  3PLM 1.0845772 -1.60171524  0.21553624         NA        NA
#> 57   G2I7    2  3PLM 1.2984963  0.26873160  0.12783879         NA        NA
#> 58   G2I8    2  3PLM 2.1768681 -0.08931817  0.13528274         NA        NA
#> 59   G2I9    2  3PLM 1.0731850 -1.56763064  0.20570964         NA        NA
#> 60  G2I10    2  3PLM 1.5904952  0.90286013  0.28755659         NA        NA
#> 61  G2I11    2  3PLM 0.9314018  0.88217244  0.11179900         NA        NA
#> 62  G2I12    2  3PLM 1.6623481 -0.69600414  0.24746346         NA        NA
#> 63  G2I13    2  3PLM 1.1643479 -1.30658561  0.20913476         NA        NA
#> 64  G2I14    2  3PLM 1.3949743  0.22018889  0.12742097         NA        NA
#> 65  G2I15    2  3PLM 1.5256012  0.03274459  0.25930675         NA        NA
#> 66  G2I16    2  3PLM 1.6958062  0.75377859  0.21478247         NA        NA
#> 67  G2I17    2  3PLM 2.1138538 -1.19916943  0.19757122         NA        NA
#> 68  G2I18    2  3PLM 1.7193451  0.64117712  0.13803059         NA        NA
#> 69  G2I19    2  3PLM 1.0971580  0.04569867  0.20768248         NA        NA
#> 70  G2I20    2  3PLM 1.6230735  2.23612968  0.15952981         NA        NA
#> 71  G2I21    2  3PLM 2.3884669  1.58525725  0.11741983         NA        NA
#> 72  G2I22    2  3PLM 1.4332531  0.22129279  0.20450739         NA        NA
#> 73  G2I23    2  3PLM 2.0941037  0.44296609  0.27618451         NA        NA
#> 74  G2I24    2  3PLM 1.3973606  0.36942946  0.07288178         NA        NA
#> 75  G2I25    2  3PLM 1.7363207  1.41508375  0.07191836         NA        NA
#> 76  G2I26    2  3PLM 1.8858048 -0.88197041  0.23993266         NA        NA
#> 77  G2I27    2  3PLM 0.9754090 -0.59038894  0.19672791         NA        NA
#> 78  G2I28    5   GRM 1.1216203 -0.37176001  0.15040063  0.7990619  1.476485
#> 79   C2I1    2  3PLM 0.9700000 -0.46000000  0.05000000         NA        NA
#> 80   C2I2    2  3PLM 0.8500000  1.18000000  0.01000000         NA        NA
#> 81   C2I3    2  3PLM 1.4300000  1.41000000  0.10000000         NA        NA
#> 82   C2I4    2  3PLM 1.4800000  0.18000000  0.17000000         NA        NA
#> 83   C2I5    2  3PLM 1.2700000 -0.23000000  0.03000000         NA        NA
#> 84   C2I6    2  3PLM 2.0200000 -0.09000000  0.01000000         NA        NA
#> 85   C2I7    2  3PLM 1.3700000 -0.13000000  0.10000000         NA        NA
#> 86   C2I8    2  3PLM 1.6700000  1.25000000  0.19000000         NA        NA
#> 87   C2I9    2  3PLM 2.2800000 -1.01000000  0.10000000         NA        NA
#> 88  C2I10    2  3PLM 1.4200000 -1.65000000  0.11000000         NA        NA
#> 89   G3I1    2  3PLM 1.5973952 -1.04293518  0.14221272         NA        NA
#> 90   G3I2    2  3PLM 0.7882034 -0.61710502  0.16466619         NA        NA
#> 91   G3I3    2  3PLM 1.1574588  0.07488161  0.15704670         NA        NA
#> 92   G3I4    2  3PLM 1.3010656  1.72470571  0.22305233         NA        NA
#> 93   G3I5    2  3PLM 0.8030039 -1.08728779  0.23797751         NA        NA
#> 94   G3I6    2  3PLM 1.1828479 -1.37600012  0.26017944         NA        NA
#> 95   G3I7    2  3PLM 1.3732288  0.26144127  0.13453648         NA        NA
#> 96   G3I8    2  3PLM 1.9483776 -0.06609774  0.14303723         NA        NA
#> 97   G3I9    2  3PLM 1.1769474 -1.16424114  0.18147943         NA        NA
#> 98  G3I10    2  3PLM 1.4966087  0.98263824  0.28748102         NA        NA
#> 99  G3I11    2  3PLM 0.9143147  0.79757776  0.07731952         NA        NA
#> 100 G3I12    2  3PLM 1.4413693 -0.77635289  0.23363873         NA        NA
#> 101 G3I13    2  3PLM 1.1547205 -1.06358212  0.20951365         NA        NA
#> 102 G3I14    2  3PLM 1.3766468  0.30600869  0.15853392         NA        NA
#> 103 G3I15    2  3PLM 1.2913722 -0.05844768  0.17903882         NA        NA
#> 104 G3I16    2  3PLM 1.5776924  0.75034401  0.19938501         NA        NA
#> 105 G3I17    2  3PLM 1.7191587 -1.46022828  0.15855847         NA        NA
#> 106 G3I18    2  3PLM 1.3690901  0.71412800  0.11276074         NA        NA
#> 107 G3I19    2  3PLM 0.9937729  0.07222966  0.14756849         NA        NA
#> 108 G3I20    2  3PLM 1.1444553  2.29572750  0.11629624         NA        NA
#> 109 G3I21    2  3PLM 2.9714475  1.65685452  0.14005235         NA        NA
#> 110 G3I22    2  3PLM 1.1948944  0.10009253  0.09517841         NA        NA
#> 111 G3I23    2  3PLM 1.7915379  0.32938893  0.14617082         NA        NA
#> 112 G3I24    2  3PLM 1.1522431  0.33724431  0.05080719         NA        NA
#> 113 G3I25    2  3PLM 1.3842436  1.36390607  0.04998934         NA        NA
#> 114 G3I26    2  3PLM 1.7163634 -0.96107355  0.19794829         NA        NA
#> 115 G3I27    2  3PLM 0.9521484 -0.60866897  0.16408973         NA        NA
#> 116 G3I28    5   GRM 0.9955303 -0.31405937  0.18617231  0.8178505  1.523770
#> 
#> $group
#> $group$Group1
#>       id cats model     par.1       par.2       par.3      par.4     par.5
#> 1   C1I1    2  3PLM 0.7600000  1.46000000  0.26000000         NA        NA
#> 2   C1I2    2  3PLM 1.9200000 -1.05000000  0.18000000         NA        NA
#> 3   C1I3    2  3PLM 0.9300000  0.39000000  0.10000000         NA        NA
#> 4   C1I4    2  3PLM 1.0500000 -0.41000000  0.20000000         NA        NA
#> 5   C1I5    2  3PLM 0.8700000 -0.12000000  0.16000000         NA        NA
#> 6   C1I6    2  3PLM 1.7000000  0.63000000  0.07000000         NA        NA
#> 7   C1I7    2  3PLM 0.9100000  1.02000000  0.12000000         NA        NA
#> 8   C1I8    2  3PLM 0.8400000  0.80000000  0.11000000         NA        NA
#> 9   C1I9    2  3PLM 0.8500000  0.85000000  0.26000000         NA        NA
#> 10 C1I10    2  3PLM 1.5300000  0.09000000  0.14000000         NA        NA
#> 11  G1I1    2  3PLM 0.9601094 -0.46323480  0.15743926         NA        NA
#> 12  G1I2    2  3PLM 0.8792351  1.21120501  0.10365721         NA        NA
#> 13  G1I3    2  3PLM 1.4922037  1.31254694  0.18433436         NA        NA
#> 14  G1I4    2  3PLM 1.5246246  0.26709705  0.29751014         NA        NA
#> 15  G1I5    2  3PLM 1.3343885 -0.18476888  0.15642097         NA        NA
#> 16  G1I6    2  3PLM 2.1457296  0.01122821  0.08373745         NA        NA
#> 17  G1I7    2  3PLM 1.4492421 -0.07232467  0.19831474         NA        NA
#> 18  G1I8    2  3PLM 2.4981366  1.18142115  0.32459740         NA        NA
#> 19  G1I9    2  3PLM 2.4036619 -0.93389393  0.24090792         NA        NA
#> 20 G1I10    2  3PLM 1.2603565 -1.75076139  0.23712596         NA        NA
#> 21 G1I11    2  3PLM 1.5482444 -1.11447546  0.20564573         NA        NA
#> 22 G1I12    2  3PLM 0.7501881 -0.81274221  0.19483606         NA        NA
#> 23 G1I13    2  3PLM 1.0391850 -0.13129452  0.20305864         NA        NA
#> 24 G1I14    2  3PLM 1.4641086  1.76817734  0.29587746         NA        NA
#> 25 G1I15    2  3PLM 0.8642507 -1.44005202  0.21246953         NA        NA
#> 26 G1I16    2  3PLM 1.0481967 -1.94936219  0.22996558         NA        NA
#> 27 G1I17    2  3PLM 1.0413440  0.24963719  0.15297903         NA        NA
#> 28 G1I18    2  3PLM 2.0821164 -0.08782853  0.23783714         NA        NA
#> 29 G1I19    2  3PLM 1.3246385 -1.38457105  0.20248019         NA        NA
#> 30 G1I20    2  3PLM 1.0078270  0.50901778  0.22362359         NA        NA
#> 31 G1I21    2  3PLM 0.9153677  0.76786397  0.13229024         NA        NA
#> 32 G1I22    2  3PLM 1.7709516 -0.64567759  0.35461144         NA        NA
#> 33 G1I23    2  3PLM 1.3100376 -1.16619310  0.25418714         NA        NA
#> 34 G1I24    2  3PLM 1.6392138  0.31576653  0.22824910         NA        NA
#> 35 G1I25    2  3PLM 1.5820253 -0.12733170  0.25023679         NA        NA
#> 36 G1I26    2  3PLM 1.8799489  0.66003503  0.25368918         NA        NA
#> 37 G1I27    2  3PLM 1.6294671 -1.53964203  0.27165399         NA        NA
#> 38 G1I28    2  3PLM 1.3331814  0.56546916  0.14913136         NA        NA
#> 39 G1I29    2  3PLM 0.9154613 -0.37004169  0.12429485         NA        NA
#> 40 G1I30    2  3PLM 1.0026860  2.29722730  0.16842602         NA        NA
#> 41 G1I31    2  3PLM 2.3881153  1.64223591  0.18163836         NA        NA
#> 42 G1I32    2  3PLM 1.1031428 -0.07825039  0.14531435         NA        NA
#> 43 G1I33    2  3PLM 1.5998724  0.17672610  0.15654349         NA        NA
#> 44 G1I34    2  3PLM 1.3266408  0.24185668  0.10791955         NA        NA
#> 45 G1I35    2  3PLM 1.3307122  1.29621629  0.06712955         NA        NA
#> 46 G1I36    2  3PLM 1.4375780 -1.23256008  0.21920779         NA        NA
#> 47 G1I37    2  3PLM 1.0523257 -0.59860569  0.26551312         NA        NA
#> 48 G1I38    5   GRM 1.0557817 -0.36687770  0.21653854  0.8654794  1.435633
#> 49 C1I11    5   GRM 1.2300000 -2.08000000 -1.35000000 -0.7100000 -0.120000
#> 50 C1I12    5   GRM 0.8800000 -0.76000000 -0.01000000  0.6700000  1.250000
#> 
#> $group$Group2
#>       id cats model     par.1       par.2       par.3      par.4     par.5
#> 1   C1I1    2  3PLM 0.7600000  1.46000000  0.26000000         NA        NA
#> 2   C1I2    2  3PLM 1.9200000 -1.05000000  0.18000000         NA        NA
#> 3   C1I3    2  3PLM 0.9300000  0.39000000  0.10000000         NA        NA
#> 4   C1I4    2  3PLM 1.0500000 -0.41000000  0.20000000         NA        NA
#> 5   C1I5    2  3PLM 0.8700000 -0.12000000  0.16000000         NA        NA
#> 6   C1I6    2  3PLM 1.7000000  0.63000000  0.07000000         NA        NA
#> 7   C1I7    2  3PLM 0.9100000  1.02000000  0.12000000         NA        NA
#> 8   C1I8    2  3PLM 0.8400000  0.80000000  0.11000000         NA        NA
#> 9   C1I9    2  3PLM 0.8500000  0.85000000  0.26000000         NA        NA
#> 10 C1I10    2  3PLM 1.5300000  0.09000000  0.14000000         NA        NA
#> 11 C1I11    5   GRM 1.2300000 -2.08000000 -1.35000000 -0.7100000 -0.120000
#> 12 C1I12    5   GRM 0.8800000 -0.76000000 -0.01000000  0.6700000  1.250000
#> 13  G2I1    2  3PLM 1.7426027 -0.86693468  0.23241989         NA        NA
#> 14  G2I2    2  3PLM 0.8271656 -0.44864390  0.21072500         NA        NA
#> 15  G2I3    2  3PLM 1.0594812  0.08411544  0.17868700         NA        NA
#> 16  G2I4    2  3PLM 1.4496890  1.49538158  0.17845021         NA        NA
#> 17  G2I5    2  3PLM 0.7000088 -1.68320349  0.19753124         NA        NA
#> 18  G2I6    2  3PLM 1.0845772 -1.60171524  0.21553624         NA        NA
#> 19  G2I7    2  3PLM 1.2984963  0.26873160  0.12783879         NA        NA
#> 20  G2I8    2  3PLM 2.1768681 -0.08931817  0.13528274         NA        NA
#> 21  G2I9    2  3PLM 1.0731850 -1.56763064  0.20570964         NA        NA
#> 22 G2I10    2  3PLM 1.5904952  0.90286013  0.28755659         NA        NA
#> 23 G2I11    2  3PLM 0.9314018  0.88217244  0.11179900         NA        NA
#> 24 G2I12    2  3PLM 1.6623481 -0.69600414  0.24746346         NA        NA
#> 25 G2I13    2  3PLM 1.1643479 -1.30658561  0.20913476         NA        NA
#> 26 G2I14    2  3PLM 1.3949743  0.22018889  0.12742097         NA        NA
#> 27 G2I15    2  3PLM 1.5256012  0.03274459  0.25930675         NA        NA
#> 28 G2I16    2  3PLM 1.6958062  0.75377859  0.21478247         NA        NA
#> 29 G2I17    2  3PLM 2.1138538 -1.19916943  0.19757122         NA        NA
#> 30 G2I18    2  3PLM 1.7193451  0.64117712  0.13803059         NA        NA
#> 31 G2I19    2  3PLM 1.0971580  0.04569867  0.20768248         NA        NA
#> 32 G2I20    2  3PLM 1.6230735  2.23612968  0.15952981         NA        NA
#> 33 G2I21    2  3PLM 2.3884669  1.58525725  0.11741983         NA        NA
#> 34 G2I22    2  3PLM 1.4332531  0.22129279  0.20450739         NA        NA
#> 35 G2I23    2  3PLM 2.0941037  0.44296609  0.27618451         NA        NA
#> 36 G2I24    2  3PLM 1.3973606  0.36942946  0.07288178         NA        NA
#> 37 G2I25    2  3PLM 1.7363207  1.41508375  0.07191836         NA        NA
#> 38 G2I26    2  3PLM 1.8858048 -0.88197041  0.23993266         NA        NA
#> 39 G2I27    2  3PLM 0.9754090 -0.59038894  0.19672791         NA        NA
#> 40 G2I28    5   GRM 1.1216203 -0.37176001  0.15040063  0.7990619  1.476485
#> 41  C2I1    2  3PLM 0.9700000 -0.46000000  0.05000000         NA        NA
#> 42  C2I2    2  3PLM 0.8500000  1.18000000  0.01000000         NA        NA
#> 43  C2I3    2  3PLM 1.4300000  1.41000000  0.10000000         NA        NA
#> 44  C2I4    2  3PLM 1.4800000  0.18000000  0.17000000         NA        NA
#> 45  C2I5    2  3PLM 1.2700000 -0.23000000  0.03000000         NA        NA
#> 46  C2I6    2  3PLM 2.0200000 -0.09000000  0.01000000         NA        NA
#> 47  C2I7    2  3PLM 1.3700000 -0.13000000  0.10000000         NA        NA
#> 48  C2I8    2  3PLM 1.6700000  1.25000000  0.19000000         NA        NA
#> 49  C2I9    2  3PLM 2.2800000 -1.01000000  0.10000000         NA        NA
#> 50 C2I10    2  3PLM 1.4200000 -1.65000000  0.11000000         NA        NA
#> 
#> $group$Group3
#>       id cats model     par.1       par.2      par.3     par.4   par.5
#> 1   C2I1    2  3PLM 0.9700000 -0.46000000 0.05000000        NA      NA
#> 2   C2I2    2  3PLM 0.8500000  1.18000000 0.01000000        NA      NA
#> 3   C2I3    2  3PLM 1.4300000  1.41000000 0.10000000        NA      NA
#> 4   C2I4    2  3PLM 1.4800000  0.18000000 0.17000000        NA      NA
#> 5   C2I5    2  3PLM 1.2700000 -0.23000000 0.03000000        NA      NA
#> 6   C2I6    2  3PLM 2.0200000 -0.09000000 0.01000000        NA      NA
#> 7   C2I7    2  3PLM 1.3700000 -0.13000000 0.10000000        NA      NA
#> 8   C2I8    2  3PLM 1.6700000  1.25000000 0.19000000        NA      NA
#> 9   C2I9    2  3PLM 2.2800000 -1.01000000 0.10000000        NA      NA
#> 10 C2I10    2  3PLM 1.4200000 -1.65000000 0.11000000        NA      NA
#> 11  G3I1    2  3PLM 1.5973952 -1.04293518 0.14221272        NA      NA
#> 12  G3I2    2  3PLM 0.7882034 -0.61710502 0.16466619        NA      NA
#> 13  G3I3    2  3PLM 1.1574588  0.07488161 0.15704670        NA      NA
#> 14  G3I4    2  3PLM 1.3010656  1.72470571 0.22305233        NA      NA
#> 15  G3I5    2  3PLM 0.8030039 -1.08728779 0.23797751        NA      NA
#> 16  G3I6    2  3PLM 1.1828479 -1.37600012 0.26017944        NA      NA
#> 17  G3I7    2  3PLM 1.3732288  0.26144127 0.13453648        NA      NA
#> 18  G3I8    2  3PLM 1.9483776 -0.06609774 0.14303723        NA      NA
#> 19  G3I9    2  3PLM 1.1769474 -1.16424114 0.18147943        NA      NA
#> 20 G3I10    2  3PLM 1.4966087  0.98263824 0.28748102        NA      NA
#> 21 G3I11    2  3PLM 0.9143147  0.79757776 0.07731952        NA      NA
#> 22 G3I12    2  3PLM 1.4413693 -0.77635289 0.23363873        NA      NA
#> 23 G3I13    2  3PLM 1.1547205 -1.06358212 0.20951365        NA      NA
#> 24 G3I14    2  3PLM 1.3766468  0.30600869 0.15853392        NA      NA
#> 25 G3I15    2  3PLM 1.2913722 -0.05844768 0.17903882        NA      NA
#> 26 G3I16    2  3PLM 1.5776924  0.75034401 0.19938501        NA      NA
#> 27 G3I17    2  3PLM 1.7191587 -1.46022828 0.15855847        NA      NA
#> 28 G3I18    2  3PLM 1.3690901  0.71412800 0.11276074        NA      NA
#> 29 G3I19    2  3PLM 0.9937729  0.07222966 0.14756849        NA      NA
#> 30 G3I20    2  3PLM 1.1444553  2.29572750 0.11629624        NA      NA
#> 31 G3I21    2  3PLM 2.9714475  1.65685452 0.14005235        NA      NA
#> 32 G3I22    2  3PLM 1.1948944  0.10009253 0.09517841        NA      NA
#> 33 G3I23    2  3PLM 1.7915379  0.32938893 0.14617082        NA      NA
#> 34 G3I24    2  3PLM 1.1522431  0.33724431 0.05080719        NA      NA
#> 35 G3I25    2  3PLM 1.3842436  1.36390607 0.04998934        NA      NA
#> 36 G3I26    2  3PLM 1.7163634 -0.96107355 0.19794829        NA      NA
#> 37 G3I27    2  3PLM 0.9521484 -0.60866897 0.16408973        NA      NA
#> 38 G3I28    5   GRM 0.9955303 -0.31405937 0.18617231 0.8178505 1.52377
#> 
#> 

# Extract the standard error estimates
getirt(fit.3, what = "se.est")
#> $overall
#>        id cats model      par.1      par.2      par.3      par.4      par.5
#> 1    C1I1    2  3PLM         NA         NA         NA         NA         NA
#> 2    C1I2    2  3PLM         NA         NA         NA         NA         NA
#> 3    C1I3    2  3PLM         NA         NA         NA         NA         NA
#> 4    C1I4    2  3PLM         NA         NA         NA         NA         NA
#> 5    C1I5    2  3PLM         NA         NA         NA         NA         NA
#> 6    C1I6    2  3PLM         NA         NA         NA         NA         NA
#> 7    C1I7    2  3PLM         NA         NA         NA         NA         NA
#> 8    C1I8    2  3PLM         NA         NA         NA         NA         NA
#> 9    C1I9    2  3PLM         NA         NA         NA         NA         NA
#> 10  C1I10    2  3PLM         NA         NA         NA         NA         NA
#> 11   G1I1    2  3PLM 0.10426168 0.19914118 0.06583145         NA         NA
#> 12   G1I2    2  3PLM 0.13011827 0.13758723 0.03866179         NA         NA
#> 13   G1I3    2  3PLM 0.25521536 0.09288148 0.02598262         NA         NA
#> 14   G1I4    2  3PLM 0.21216024 0.12139692 0.04369469         NA         NA
#> 15   G1I5    2  3PLM 0.13737157 0.12762451 0.05144114         NA         NA
#> 16   G1I6    2  3PLM 0.17837536 0.05370898 0.02614114         NA         NA
#> 17   G1I7    2  3PLM 0.15675133 0.12128691 0.04842894         NA         NA
#> 18   G1I8    2  3PLM 0.46148303 0.06925484 0.01882411         NA         NA
#> 19   G1I9    2  3PLM 0.25992134 0.10854364 0.06176901         NA         NA
#> 20  G1I10    2  3PLM 0.12525002 0.22773298 0.09804119         NA         NA
#> 21  G1I11    2  3PLM 0.14906877 0.16092483 0.07743255         NA         NA
#> 22  G1I12    2  3PLM 0.09144445 0.31119356 0.08351947         NA         NA
#> 23  G1I13    2  3PLM 0.13678756 0.20744666 0.06854745         NA         NA
#> 24  G1I14    2  3PLM 0.36846246 0.15351046 0.02514203         NA         NA
#> 25  G1I15    2  3PLM 0.09478648 0.29483636 0.09102878         NA         NA
#> 26  G1I16    2  3PLM 0.10986557 0.26721069 0.09771458         NA         NA
#> 27  G1I17    2  3PLM 0.12646970 0.15451738 0.05176849         NA         NA
#> 28  G1I18    2  3PLM 0.20058125 0.08013120 0.03656393         NA         NA
#> 29  G1I19    2  3PLM 0.12663826 0.19277821 0.08479494         NA         NA
#> 30  G1I20    2  3PLM 0.17409361 0.19113321 0.06033487         NA         NA
#> 31  G1I21    2  3PLM 0.13256306 0.15087630 0.04704321         NA         NA
#> 32  G1I22    2  3PLM 0.21340356 0.15629299 0.06435608         NA         NA
#> 33  G1I23    2  3PLM 0.14898851 0.22567950 0.09214904         NA         NA
#> 34  G1I24    2  3PLM 0.19188388 0.09525298 0.03753622         NA         NA
#> 35  G1I25    2  3PLM 0.17877681 0.12332508 0.04964582         NA         NA
#> 36  G1I26    2  3PLM 0.25183655 0.07738128 0.02912676         NA         NA
#> 37  G1I27    2  3PLM 0.17771713 0.20258618 0.10060823         NA         NA
#> 38  G1I28    2  3PLM 0.16142785 0.09763951 0.03660579         NA         NA
#> 39  G1I29    2  3PLM 0.08893330 0.17133724 0.05492922         NA         NA
#> 40  G1I30    2  3PLM 0.26897282 0.25352870 0.03089404         NA         NA
#> 41  G1I31    2  3PLM 0.46169957 0.09359258 0.01431405         NA         NA
#> 42  G1I32    2  3PLM 0.11905116 0.15068461 0.05352540         NA         NA
#> 43  G1I33    2  3PLM 0.17012521 0.08891040 0.03755055         NA         NA
#> 44  G1I34    2  3PLM 0.13288626 0.09359115 0.03655936         NA         NA
#> 45  G1I35    2  3PLM 0.16791769 0.08361518 0.02044036         NA         NA
#> 46  G1I36    2  3PLM 0.14309878 0.18590021 0.08516936         NA         NA
#> 47  G1I37    2  3PLM 0.14194404 0.26869625 0.08697510         NA         NA
#> 48  G1I38    5   GRM 0.05705443 0.05322210 0.04926777 0.06290244 0.08428075
#> 49  C1I11    5   GRM         NA         NA         NA         NA         NA
#> 50  C1I12    5   GRM         NA         NA         NA         NA         NA
#> 51   G2I1    2  3PLM 0.18701977 0.17226079 0.09591121         NA         NA
#> 52   G2I2    2  3PLM 0.11055188 0.30549386 0.09111811         NA         NA
#> 53   G2I3    2  3PLM 0.12900982 0.20135385 0.07488982         NA         NA
#> 54   G2I4    2  3PLM 0.28006193 0.09631836 0.03881020         NA         NA
#> 55   G2I5    2  3PLM 0.10228622 0.40358987 0.08893723         NA         NA
#> 56   G2I6    2  3PLM 0.14204281 0.29399636 0.09476517         NA         NA
#> 57   G2I7    2  3PLM 0.13098448 0.12235808 0.05447273         NA         NA
#> 58   G2I8    2  3PLM 0.18554592 0.08066166 0.05149244         NA         NA
#> 59   G2I9    2  3PLM 0.13004997 0.27643093 0.09157031         NA         NA
#> 60  G2I10    2  3PLM 0.27183179 0.12027438 0.05121721         NA         NA
#> 61  G2I11    2  3PLM 0.12329519 0.14644859 0.05005698         NA         NA
#> 62  G2I12    2  3PLM 0.18329666 0.18051784 0.09670036         NA         NA
#> 63  G2I13    2  3PLM 0.13036275 0.23937423 0.09227547         NA         NA
#> 64  G2I14    2  3PLM 0.13691011 0.11354387 0.05286438         NA         NA
#> 65  G2I15    2  3PLM 0.19565359 0.17450556 0.08012343         NA         NA
#> 66  G2I16    2  3PLM 0.22938735 0.09966642 0.04736991         NA         NA
#> 67  G2I17    2  3PLM 0.21703924 0.14435249 0.08787073         NA         NA
#> 68  G2I18    2  3PLM 0.19775206 0.08387898 0.04368134         NA         NA
#> 69  G2I19    2  3PLM 0.14342944 0.22082957 0.08261270         NA         NA
#> 70  G2I20    2  3PLM 0.45488499 0.19197090 0.02387855         NA         NA
#> 71  G2I21    2  3PLM 0.35929005 0.05999541 0.01724135         NA         NA
#> 72  G2I22    2  3PLM 0.17462991 0.15301616 0.06937934         NA         NA
#> 73  G2I23    2  3PLM 0.26460348 0.09594883 0.04984224         NA         NA
#> 74  G2I24    2  3PLM 0.11632156 0.07551221 0.03397862         NA         NA
#> 75  G2I25    2  3PLM 0.21505543 0.06239024 0.02239864         NA         NA
#> 76  G2I26    2  3PLM 0.19522559 0.16121604 0.09680686         NA         NA
#> 77  G2I27    2  3PLM 0.11542196 0.24707543 0.08656505         NA         NA
#> 78  G2I28    5   GRM 0.07020432 0.06705796 0.04891269 0.04705132 0.07019674
#> 79   C2I1    2  3PLM         NA         NA         NA         NA         NA
#> 80   C2I2    2  3PLM         NA         NA         NA         NA         NA
#> 81   C2I3    2  3PLM         NA         NA         NA         NA         NA
#> 82   C2I4    2  3PLM         NA         NA         NA         NA         NA
#> 83   C2I5    2  3PLM         NA         NA         NA         NA         NA
#> 84   C2I6    2  3PLM         NA         NA         NA         NA         NA
#> 85   C2I7    2  3PLM         NA         NA         NA         NA         NA
#> 86   C2I8    2  3PLM         NA         NA         NA         NA         NA
#> 87   C2I9    2  3PLM         NA         NA         NA         NA         NA
#> 88  C2I10    2  3PLM         NA         NA         NA         NA         NA
#> 89   G3I1    2  3PLM 0.13986318 0.11556266 0.04828949         NA         NA
#> 90   G3I2    2  3PLM 0.08945397 0.25055681 0.06577738         NA         NA
#> 91   G3I3    2  3PLM 0.12166153 0.12145104 0.03777657         NA         NA
#> 92   G3I4    2  3PLM 0.23742580 0.12346854 0.02115715         NA         NA
#> 93   G3I5    2  3PLM 0.09747072 0.33207410 0.08571279         NA         NA
#> 94   G3I6    2  3PLM 0.13415716 0.24170387 0.08325843         NA         NA
#> 95   G3I7    2  3PLM 0.13182139 0.08662754 0.02782581         NA         NA
#> 96   G3I8    2  3PLM 0.17855395 0.06374914 0.02406454         NA         NA
#> 97   G3I9    2  3PLM 0.11323201 0.18515456 0.06633175         NA         NA
#> 98  G3I10    2  3PLM 0.21267233 0.09577230 0.02369318         NA         NA
#> 99  G3I11    2  3PLM 0.09929092 0.11096928 0.02888083         NA         NA
#> 100 G3I12    2  3PLM 0.16578402 0.15209966 0.05508812         NA         NA
#> 101 G3I13    2  3PLM 0.12510673 0.20472281 0.06955329         NA         NA
#> 102 G3I14    2  3PLM 0.14193927 0.08926242 0.02823279         NA         NA
#> 103 G3I15    2  3PLM 0.13431761 0.11420069 0.03726663         NA         NA
#> 104 G3I16    2  3PLM 0.17708800 0.07754598 0.02161000         NA         NA
#> 105 G3I17    2  3PLM 0.14820688 0.12802058 0.05977877         NA         NA
#> 106 G3I18    2  3PLM 0.14295569 0.07697678 0.02235150         NA         NA
#> 107 G3I19    2  3PLM 0.11089570 0.14924845 0.04377532         NA         NA
#> 108 G3I20    2  3PLM 0.23479290 0.19832103 0.01759904         NA         NA
#> 109 G3I21    2  3PLM 0.54769558 0.06208590 0.01022861         NA         NA
#> 110 G3I22    2  3PLM 0.11351503 0.09624776 0.03123112         NA         NA
#> 111 G3I23    2  3PLM 0.16621787 0.06369671 0.02083676         NA         NA
#> 112 G3I24    2  3PLM 0.09162888 0.07576239 0.02065526         NA         NA
#> 113 G3I25    2  3PLM 0.14521791 0.07831223 0.01301757         NA         NA
#> 114 G3I26    2  3PLM 0.16527458 0.12057791 0.04950403         NA         NA
#> 115 G3I27    2  3PLM 0.09992189 0.19902502 0.05904246         NA         NA
#> 116 G3I28    5   GRM 0.04870285 0.05360023 0.05471689 0.06669018 0.08789166
#> 
#> $group
#> $group$Group1
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   C1I1    2  3PLM         NA         NA         NA         NA         NA
#> 2   C1I2    2  3PLM         NA         NA         NA         NA         NA
#> 3   C1I3    2  3PLM         NA         NA         NA         NA         NA
#> 4   C1I4    2  3PLM         NA         NA         NA         NA         NA
#> 5   C1I5    2  3PLM         NA         NA         NA         NA         NA
#> 6   C1I6    2  3PLM         NA         NA         NA         NA         NA
#> 7   C1I7    2  3PLM         NA         NA         NA         NA         NA
#> 8   C1I8    2  3PLM         NA         NA         NA         NA         NA
#> 9   C1I9    2  3PLM         NA         NA         NA         NA         NA
#> 10 C1I10    2  3PLM         NA         NA         NA         NA         NA
#> 11  G1I1    2  3PLM 0.10426168 0.19914118 0.06583145         NA         NA
#> 12  G1I2    2  3PLM 0.13011827 0.13758723 0.03866179         NA         NA
#> 13  G1I3    2  3PLM 0.25521536 0.09288148 0.02598262         NA         NA
#> 14  G1I4    2  3PLM 0.21216024 0.12139692 0.04369469         NA         NA
#> 15  G1I5    2  3PLM 0.13737157 0.12762451 0.05144114         NA         NA
#> 16  G1I6    2  3PLM 0.17837536 0.05370898 0.02614114         NA         NA
#> 17  G1I7    2  3PLM 0.15675133 0.12128691 0.04842894         NA         NA
#> 18  G1I8    2  3PLM 0.46148303 0.06925484 0.01882411         NA         NA
#> 19  G1I9    2  3PLM 0.25992134 0.10854364 0.06176901         NA         NA
#> 20 G1I10    2  3PLM 0.12525002 0.22773298 0.09804119         NA         NA
#> 21 G1I11    2  3PLM 0.14906877 0.16092483 0.07743255         NA         NA
#> 22 G1I12    2  3PLM 0.09144445 0.31119356 0.08351947         NA         NA
#> 23 G1I13    2  3PLM 0.13678756 0.20744666 0.06854745         NA         NA
#> 24 G1I14    2  3PLM 0.36846246 0.15351046 0.02514203         NA         NA
#> 25 G1I15    2  3PLM 0.09478648 0.29483636 0.09102878         NA         NA
#> 26 G1I16    2  3PLM 0.10986557 0.26721069 0.09771458         NA         NA
#> 27 G1I17    2  3PLM 0.12646970 0.15451738 0.05176849         NA         NA
#> 28 G1I18    2  3PLM 0.20058125 0.08013120 0.03656393         NA         NA
#> 29 G1I19    2  3PLM 0.12663826 0.19277821 0.08479494         NA         NA
#> 30 G1I20    2  3PLM 0.17409361 0.19113321 0.06033487         NA         NA
#> 31 G1I21    2  3PLM 0.13256306 0.15087630 0.04704321         NA         NA
#> 32 G1I22    2  3PLM 0.21340356 0.15629299 0.06435608         NA         NA
#> 33 G1I23    2  3PLM 0.14898851 0.22567950 0.09214904         NA         NA
#> 34 G1I24    2  3PLM 0.19188388 0.09525298 0.03753622         NA         NA
#> 35 G1I25    2  3PLM 0.17877681 0.12332508 0.04964582         NA         NA
#> 36 G1I26    2  3PLM 0.25183655 0.07738128 0.02912676         NA         NA
#> 37 G1I27    2  3PLM 0.17771713 0.20258618 0.10060823         NA         NA
#> 38 G1I28    2  3PLM 0.16142785 0.09763951 0.03660579         NA         NA
#> 39 G1I29    2  3PLM 0.08893330 0.17133724 0.05492922         NA         NA
#> 40 G1I30    2  3PLM 0.26897282 0.25352870 0.03089404         NA         NA
#> 41 G1I31    2  3PLM 0.46169957 0.09359258 0.01431405         NA         NA
#> 42 G1I32    2  3PLM 0.11905116 0.15068461 0.05352540         NA         NA
#> 43 G1I33    2  3PLM 0.17012521 0.08891040 0.03755055         NA         NA
#> 44 G1I34    2  3PLM 0.13288626 0.09359115 0.03655936         NA         NA
#> 45 G1I35    2  3PLM 0.16791769 0.08361518 0.02044036         NA         NA
#> 46 G1I36    2  3PLM 0.14309878 0.18590021 0.08516936         NA         NA
#> 47 G1I37    2  3PLM 0.14194404 0.26869625 0.08697510         NA         NA
#> 48 G1I38    5   GRM 0.05705443 0.05322210 0.04926777 0.06290244 0.08428075
#> 49 C1I11    5   GRM         NA         NA         NA         NA         NA
#> 50 C1I12    5   GRM         NA         NA         NA         NA         NA
#> 
#> $group$Group2
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   C1I1    2  3PLM         NA         NA         NA         NA         NA
#> 2   C1I2    2  3PLM         NA         NA         NA         NA         NA
#> 3   C1I3    2  3PLM         NA         NA         NA         NA         NA
#> 4   C1I4    2  3PLM         NA         NA         NA         NA         NA
#> 5   C1I5    2  3PLM         NA         NA         NA         NA         NA
#> 6   C1I6    2  3PLM         NA         NA         NA         NA         NA
#> 7   C1I7    2  3PLM         NA         NA         NA         NA         NA
#> 8   C1I8    2  3PLM         NA         NA         NA         NA         NA
#> 9   C1I9    2  3PLM         NA         NA         NA         NA         NA
#> 10 C1I10    2  3PLM         NA         NA         NA         NA         NA
#> 11 C1I11    5   GRM         NA         NA         NA         NA         NA
#> 12 C1I12    5   GRM         NA         NA         NA         NA         NA
#> 13  G2I1    2  3PLM 0.18701977 0.17226079 0.09591121         NA         NA
#> 14  G2I2    2  3PLM 0.11055188 0.30549386 0.09111811         NA         NA
#> 15  G2I3    2  3PLM 0.12900982 0.20135385 0.07488982         NA         NA
#> 16  G2I4    2  3PLM 0.28006193 0.09631836 0.03881020         NA         NA
#> 17  G2I5    2  3PLM 0.10228622 0.40358987 0.08893723         NA         NA
#> 18  G2I6    2  3PLM 0.14204281 0.29399636 0.09476517         NA         NA
#> 19  G2I7    2  3PLM 0.13098448 0.12235808 0.05447273         NA         NA
#> 20  G2I8    2  3PLM 0.18554592 0.08066166 0.05149244         NA         NA
#> 21  G2I9    2  3PLM 0.13004997 0.27643093 0.09157031         NA         NA
#> 22 G2I10    2  3PLM 0.27183179 0.12027438 0.05121721         NA         NA
#> 23 G2I11    2  3PLM 0.12329519 0.14644859 0.05005698         NA         NA
#> 24 G2I12    2  3PLM 0.18329666 0.18051784 0.09670036         NA         NA
#> 25 G2I13    2  3PLM 0.13036275 0.23937423 0.09227547         NA         NA
#> 26 G2I14    2  3PLM 0.13691011 0.11354387 0.05286438         NA         NA
#> 27 G2I15    2  3PLM 0.19565359 0.17450556 0.08012343         NA         NA
#> 28 G2I16    2  3PLM 0.22938735 0.09966642 0.04736991         NA         NA
#> 29 G2I17    2  3PLM 0.21703924 0.14435249 0.08787073         NA         NA
#> 30 G2I18    2  3PLM 0.19775206 0.08387898 0.04368134         NA         NA
#> 31 G2I19    2  3PLM 0.14342944 0.22082957 0.08261270         NA         NA
#> 32 G2I20    2  3PLM 0.45488499 0.19197090 0.02387855         NA         NA
#> 33 G2I21    2  3PLM 0.35929005 0.05999541 0.01724135         NA         NA
#> 34 G2I22    2  3PLM 0.17462991 0.15301616 0.06937934         NA         NA
#> 35 G2I23    2  3PLM 0.26460348 0.09594883 0.04984224         NA         NA
#> 36 G2I24    2  3PLM 0.11632156 0.07551221 0.03397862         NA         NA
#> 37 G2I25    2  3PLM 0.21505543 0.06239024 0.02239864         NA         NA
#> 38 G2I26    2  3PLM 0.19522559 0.16121604 0.09680686         NA         NA
#> 39 G2I27    2  3PLM 0.11542196 0.24707543 0.08656505         NA         NA
#> 40 G2I28    5   GRM 0.07020432 0.06705796 0.04891269 0.04705132 0.07019674
#> 41  C2I1    2  3PLM         NA         NA         NA         NA         NA
#> 42  C2I2    2  3PLM         NA         NA         NA         NA         NA
#> 43  C2I3    2  3PLM         NA         NA         NA         NA         NA
#> 44  C2I4    2  3PLM         NA         NA         NA         NA         NA
#> 45  C2I5    2  3PLM         NA         NA         NA         NA         NA
#> 46  C2I6    2  3PLM         NA         NA         NA         NA         NA
#> 47  C2I7    2  3PLM         NA         NA         NA         NA         NA
#> 48  C2I8    2  3PLM         NA         NA         NA         NA         NA
#> 49  C2I9    2  3PLM         NA         NA         NA         NA         NA
#> 50 C2I10    2  3PLM         NA         NA         NA         NA         NA
#> 
#> $group$Group3
#>       id cats model      par.1      par.2      par.3      par.4      par.5
#> 1   C2I1    2  3PLM         NA         NA         NA         NA         NA
#> 2   C2I2    2  3PLM         NA         NA         NA         NA         NA
#> 3   C2I3    2  3PLM         NA         NA         NA         NA         NA
#> 4   C2I4    2  3PLM         NA         NA         NA         NA         NA
#> 5   C2I5    2  3PLM         NA         NA         NA         NA         NA
#> 6   C2I6    2  3PLM         NA         NA         NA         NA         NA
#> 7   C2I7    2  3PLM         NA         NA         NA         NA         NA
#> 8   C2I8    2  3PLM         NA         NA         NA         NA         NA
#> 9   C2I9    2  3PLM         NA         NA         NA         NA         NA
#> 10 C2I10    2  3PLM         NA         NA         NA         NA         NA
#> 11  G3I1    2  3PLM 0.13986318 0.11556266 0.04828949         NA         NA
#> 12  G3I2    2  3PLM 0.08945397 0.25055681 0.06577738         NA         NA
#> 13  G3I3    2  3PLM 0.12166153 0.12145104 0.03777657         NA         NA
#> 14  G3I4    2  3PLM 0.23742580 0.12346854 0.02115715         NA         NA
#> 15  G3I5    2  3PLM 0.09747072 0.33207410 0.08571279         NA         NA
#> 16  G3I6    2  3PLM 0.13415716 0.24170387 0.08325843         NA         NA
#> 17  G3I7    2  3PLM 0.13182139 0.08662754 0.02782581         NA         NA
#> 18  G3I8    2  3PLM 0.17855395 0.06374914 0.02406454         NA         NA
#> 19  G3I9    2  3PLM 0.11323201 0.18515456 0.06633175         NA         NA
#> 20 G3I10    2  3PLM 0.21267233 0.09577230 0.02369318         NA         NA
#> 21 G3I11    2  3PLM 0.09929092 0.11096928 0.02888083         NA         NA
#> 22 G3I12    2  3PLM 0.16578402 0.15209966 0.05508812         NA         NA
#> 23 G3I13    2  3PLM 0.12510673 0.20472281 0.06955329         NA         NA
#> 24 G3I14    2  3PLM 0.14193927 0.08926242 0.02823279         NA         NA
#> 25 G3I15    2  3PLM 0.13431761 0.11420069 0.03726663         NA         NA
#> 26 G3I16    2  3PLM 0.17708800 0.07754598 0.02161000         NA         NA
#> 27 G3I17    2  3PLM 0.14820688 0.12802058 0.05977877         NA         NA
#> 28 G3I18    2  3PLM 0.14295569 0.07697678 0.02235150         NA         NA
#> 29 G3I19    2  3PLM 0.11089570 0.14924845 0.04377532         NA         NA
#> 30 G3I20    2  3PLM 0.23479290 0.19832103 0.01759904         NA         NA
#> 31 G3I21    2  3PLM 0.54769558 0.06208590 0.01022861         NA         NA
#> 32 G3I22    2  3PLM 0.11351503 0.09624776 0.03123112         NA         NA
#> 33 G3I23    2  3PLM 0.16621787 0.06369671 0.02083676         NA         NA
#> 34 G3I24    2  3PLM 0.09162888 0.07576239 0.02065526         NA         NA
#> 35 G3I25    2  3PLM 0.14521791 0.07831223 0.01301757         NA         NA
#> 36 G3I26    2  3PLM 0.16527458 0.12057791 0.04950403         NA         NA
#> 37 G3I27    2  3PLM 0.09992189 0.19902502 0.05904246         NA         NA
#> 38 G3I28    5   GRM 0.04870285 0.05360023 0.05471689 0.06669018 0.08789166
#> 
#> 

# Extract the group parameter estimates (i.e., scale parameters)
getirt(fit.3, what = "group.par")
#> $Group1
#>                    mu     sigma2     sigma
#> estimates -0.01450724 1.01249024 1.0062257
#> se         0.02249989 0.03202576 0.0159138
#> 
#> $Group2
#>                   mu     sigma2      sigma
#> estimates 0.50413571 0.58171603 0.76270311
#> se        0.01705456 0.01840008 0.01206241
#> 
#> $Group3
#>                    mu     sigma2      sigma
#> estimates -0.28118672 1.66447266 1.29014443
#> se         0.02884851 0.05264841 0.02040408
#> 

# Extract the posterior latent ability distributions of the groups
getirt(fit.3, what = "weights")
#> $Group1
#>    theta       weight
#> 1  -6.00 7.662190e-12
#> 2  -5.75 2.451684e-11
#> 3  -5.50 7.119151e-11
#> 4  -5.25 1.891374e-10
#> 5  -5.00 4.676083e-10
#> 6  -4.75 1.108652e-09
#> 7  -4.50 2.644856e-09
#> 8  -4.25 6.812933e-09
#> 9  -4.00 2.081206e-08
#> 10 -3.75 8.401450e-08
#> 11 -3.50 4.930076e-07
#> 12 -3.25 4.326576e-06
#> 13 -3.00 5.074488e-05
#> 14 -2.75 5.872308e-04
#> 15 -2.50 4.332804e-03
#> 16 -2.25 1.444351e-02
#> 17 -2.00 2.242166e-02
#> 18 -1.75 2.444472e-02
#> 19 -1.50 2.775462e-02
#> 20 -1.25 3.781543e-02
#> 21 -1.00 5.911397e-02
#> 22 -0.75 7.686054e-02
#> 23 -0.50 8.228516e-02
#> 24 -0.25 1.049100e-01
#> 25  0.00 1.078960e-01
#> 26  0.25 8.138445e-02
#> 27  0.50 7.575408e-02
#> 28  0.75 8.485785e-02
#> 29  1.00 8.007142e-02
#> 30  1.25 4.680550e-02
#> 31  1.50 2.165959e-02
#> 32  1.75 1.619189e-02
#> 33  2.00 1.539816e-02
#> 34  2.25 8.170377e-03
#> 35  2.50 2.384944e-03
#> 36  2.75 7.954930e-04
#> 37  3.00 5.389228e-04
#> 38  3.25 6.906744e-04
#> 39  3.50 9.912521e-04
#> 40  3.75 9.167267e-04
#> 41  4.00 3.910055e-04
#> 42  4.25 7.029895e-05
#> 43  4.50 5.775249e-06
#> 44  4.75 2.563307e-07
#> 45  5.00 7.455139e-09
#> 46  5.25 1.703871e-10
#> 47  5.50 3.561747e-12
#> 48  5.75 7.640892e-14
#> 49  6.00 1.820327e-15
#> 
#> $Group2
#>    theta       weight
#> 1  -6.00 3.894817e-89
#> 2  -5.75 3.951307e-84
#> 3  -5.50 7.764026e-79
#> 4  -5.25 2.603004e-73
#> 5  -5.00 1.291279e-67
#> 6  -4.75 8.143648e-62
#> 7  -4.50 5.597431e-56
#> 8  -4.25 3.620629e-50
#> 9  -4.00 1.945094e-44
#> 10 -3.75 7.996024e-39
#> 11 -3.50 2.475576e-33
#> 12 -3.25 6.096104e-28
#> 13 -3.00 1.277570e-22
#> 14 -2.75 2.100737e-17
#> 15 -2.50 1.722718e-12
#> 16 -2.25 2.851606e-08
#> 17 -2.00 3.328725e-05
#> 18 -1.75 1.652181e-03
#> 19 -1.50 7.124061e-03
#> 20 -1.25 1.323762e-02
#> 21 -1.00 2.233184e-02
#> 22 -0.75 2.551256e-02
#> 23 -0.50 4.082653e-02
#> 24 -0.25 1.029622e-01
#> 25  0.00 9.278880e-02
#> 26  0.25 8.719857e-02
#> 27  0.50 1.665702e-01
#> 28  0.75 1.473383e-01
#> 29  1.00 8.777312e-02
#> 30  1.25 7.526865e-02
#> 31  1.50 6.301574e-02
#> 32  1.75 3.752019e-02
#> 33  2.00 1.769110e-02
#> 34  2.25 6.091144e-03
#> 35  2.50 1.740485e-03
#> 36  2.75 7.205189e-04
#> 37  3.00 5.733371e-04
#> 38  3.25 6.627171e-04
#> 39  3.50 6.883348e-04
#> 40  3.75 4.588606e-04
#> 41  4.00 1.743136e-04
#> 42  4.25 3.899078e-05
#> 43  4.50 5.702020e-06
#> 44  4.75 6.165462e-07
#> 45  5.00 5.519581e-08
#> 46  5.25 4.482665e-09
#> 47  5.50 3.531305e-10
#> 48  5.75 2.821138e-11
#> 49  6.00 2.345691e-12
#> 
#> $Group3
#>    theta       weight
#> 1  -6.00 6.924473e-08
#> 2  -5.75 2.943789e-07
#> 3  -5.50 1.166017e-06
#> 4  -5.25 4.289972e-06
#> 5  -5.00 1.460129e-05
#> 6  -4.75 4.573324e-05
#> 7  -4.50 1.309612e-04
#> 8  -4.25 3.402905e-04
#> 9  -4.00 7.962462e-04
#> 10 -3.75 1.668503e-03
#> 11 -3.50 3.132224e-03
#> 12 -3.25 5.329452e-03
#> 13 -3.00 8.445236e-03
#> 14 -2.75 1.295753e-02
#> 15 -2.50 1.981722e-02
#> 16 -2.25 2.965113e-02
#> 17 -2.00 3.982099e-02
#> 18 -1.75 4.375602e-02
#> 19 -1.50 4.160532e-02
#> 20 -1.25 4.561268e-02
#> 21 -1.00 6.646806e-02
#> 22 -0.75 8.026361e-02
#> 23 -0.50 6.427858e-02
#> 24 -0.25 6.407623e-02
#> 25  0.00 8.542032e-02
#> 26  0.25 8.548725e-02
#> 27  0.50 6.053382e-02
#> 28  0.75 4.722574e-02
#> 29  1.00 4.420511e-02
#> 30  1.25 3.945090e-02
#> 31  1.50 3.541132e-02
#> 32  1.75 3.265182e-02
#> 33  2.00 2.150663e-02
#> 34  2.25 9.341872e-03
#> 35  2.50 3.769655e-03
#> 36  2.75 1.886082e-03
#> 37  3.00 1.256293e-03
#> 38  3.25 1.030125e-03
#> 39  3.50 8.997076e-04
#> 40  3.75 7.280932e-04
#> 41  4.00 4.987633e-04
#> 42  4.25 2.792318e-04
#> 43  4.50 1.282712e-04
#> 44  4.75 4.948786e-05
#> 45  5.00 1.649488e-05
#> 46  5.25 4.873378e-06
#> 47  5.50 1.302434e-06
#> 48  5.75 3.195413e-07
#> 49  6.00 7.269746e-08
#> 

# 2-(2). Alternatively, MG-FIPC can be implemented by specifying the
# IDs of the items to be fixed using the 'fix.id' argument.
# Provide a character vector of fixed item IDs to 'fix.id'
fix.id <- c(paste0("C1I", 1:12), paste0("C2I", 1:10))
fit.4 <-
  est_mg(
    x = x, data = data, group.name = group.name, D = 1,
    free.group = free.group, use.gprior = TRUE,
    gprior = list(dist = "beta", params = c(5, 16)),
    EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
    fipc.method = "MEM", fix.id = fix.id
  )
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -56756.7548, Max-Change: 3.554166 EM iteration: 2, Loglike: -160424.0878, Max-Change: 0.780368 EM iteration: 3, Loglike: -159755.3559, Max-Change: 0.276071 EM iteration: 4, Loglike: -159686.2158, Max-Change: 0.151073 EM iteration: 5, Loglike: -159657.5520, Max-Change: 0.102849 EM iteration: 6, Loglike: -159642.1127, Max-Change: 0.080257 EM iteration: 7, Loglike: -159633.1411, Max-Change: 0.066343 EM iteration: 8, Loglike: -159627.6844, Max-Change: 0.055853 EM iteration: 9, Loglike: -159624.2459, Max-Change: 0.047172 EM iteration: 10, Loglike: -159622.0120, Max-Change: 0.039781 EM iteration: 11, Loglike: -159620.5195, Max-Change: 0.033468 EM iteration: 12, Loglike: -159619.4950, Max-Change: 0.028101 EM iteration: 13, Loglike: -159618.7723, Max-Change: 0.023565 EM iteration: 14, Loglike: -159618.2481, Max-Change: 0.019752 EM iteration: 15, Loglike: -159617.8570, Max-Change: 0.016845 EM iteration: 16, Loglike: -159617.5566, Max-Change: 0.014864 EM iteration: 17, Loglike: -159617.3191, Max-Change: 0.01318 EM iteration: 18, Loglike: -159617.1261, Max-Change: 0.011751 EM iteration: 19, Loglike: -159616.9650, Max-Change: 0.010538 EM iteration: 20, Loglike: -159616.8273, Max-Change: 0.009507 EM iteration: 21, Loglike: -159616.7072, Max-Change: 0.008626 EM iteration: 22, Loglike: -159616.6005, Max-Change: 0.007869 EM iteration: 23, Loglike: -159616.5042, Max-Change: 0.007214 EM iteration: 24, Loglike: -159616.4163, Max-Change: 0.006643 EM iteration: 25, Loglike: -159616.3352, Max-Change: 0.006141 EM iteration: 26, Loglike: -159616.2596, Max-Change: 0.005695 EM iteration: 27, Loglike: -159616.1888, Max-Change: 0.005295 EM iteration: 28, Loglike: -159616.1221, Max-Change: 0.004934 EM iteration: 29, Loglike: -159616.0589, Max-Change: 0.004606 EM iteration: 30, Loglike: -159615.9989, Max-Change: 0.004306 EM iteration: 31, Loglike: -159615.9416, Max-Change: 0.00403 EM iteration: 32, Loglike: -159615.8868, Max-Change: 0.003774 EM iteration: 33, Loglike: -159615.8343, Max-Change: 0.003537 EM iteration: 34, Loglike: -159615.7838, Max-Change: 0.003317 EM iteration: 35, Loglike: -159615.7352, Max-Change: 0.003111 EM iteration: 36, Loglike: -159615.6884, Max-Change: 0.002919 EM iteration: 37, Loglike: -159615.6431, Max-Change: 0.00274 EM iteration: 38, Loglike: -159615.5994, Max-Change: 0.002572 EM iteration: 39, Loglike: -159615.5570, Max-Change: 0.002416 EM iteration: 40, Loglike: -159615.5159, Max-Change: 0.002269 EM iteration: 41, Loglike: -159615.4761, Max-Change: 0.002131 EM iteration: 42, Loglike: -159615.4373, Max-Change: 0.002003 EM iteration: 43, Loglike: -159615.3997, Max-Change: 0.001882 EM iteration: 44, Loglike: -159615.3630, Max-Change: 0.00177 EM iteration: 45, Loglike: -159615.3274, Max-Change: 0.001665 EM iteration: 46, Loglike: -159615.2926, Max-Change: 0.001566 EM iteration: 47, Loglike: -159615.2587, Max-Change: 0.001474 EM iteration: 48, Loglike: -159615.2256, Max-Change: 0.001388 EM iteration: 49, Loglike: -159615.1933, Max-Change: 0.001308 EM iteration: 50, Loglike: -159615.1617, Max-Change: 0.001232 EM iteration: 51, Loglike: -159615.1309, Max-Change: 0.001162 EM iteration: 52, Loglike: -159615.1007, Max-Change: 0.001097 EM iteration: 53, Loglike: -159615.0711, Max-Change: 0.001035 EM iteration: 54, Loglike: -159615.0422, Max-Change: 0.000978 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 8.93 seconds. 

# Summary of the estimation
summary(fit.4)
#> 
#> Call:
#> est_mg(x = x, data = data, group.name = group.name, D = 1, free.group = free.group, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE, 
#>     fipc.method = "MEM", fix.id = fix.id)
#> 
#> Summary of the Data 
#>  Number of Items: 
#>   Overall: 116 unique items 
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of Cases: 
#>   Overall: 6000
#>   By group: 2000(Group1), 2000(Group2), 2000(Group3)
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 294
#>  Number of fixed items: 
#>   Overall: 22
#>   By group: 12(Group1), 22(Group2), 10(Group3)
#>  Number of E-step cycles completed: 54
#>  Maximum parameter change: 0.0009779691
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 8.26
#>  Standard error computation: 0.15
#>  Total computation: 8.93
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 
#>   Overall: 319230.1
#>   By group: 120354.311(Group1), 113985.392(Group2), 84890.382(Group3)
#> 
#>  Akaike Information Criterion (AIC): 319818.1
#>  Bayesian Information Criterion (BIC): 321787.7
#>  Item Parameters (Overall): 
#>         id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1     C1I1     2   3PLM   0.76    NA   1.46    NA   0.26    NA     NA    NA
#> 2     C1I2     2   3PLM   1.92    NA  -1.05    NA   0.18    NA     NA    NA
#> 3     C1I3     2   3PLM   0.93    NA   0.39    NA   0.10    NA     NA    NA
#> 4     C1I4     2   3PLM   1.05    NA  -0.41    NA   0.20    NA     NA    NA
#> 5     C1I5     2   3PLM   0.87    NA  -0.12    NA   0.16    NA     NA    NA
#> 6     C1I6     2   3PLM   1.70    NA   0.63    NA   0.07    NA     NA    NA
#> 7     C1I7     2   3PLM   0.91    NA   1.02    NA   0.12    NA     NA    NA
#> 8     C1I8     2   3PLM   0.84    NA   0.80    NA   0.11    NA     NA    NA
#> 9     C1I9     2   3PLM   0.85    NA   0.85    NA   0.26    NA     NA    NA
#> 10   C1I10     2   3PLM   1.53    NA   0.09    NA   0.14    NA     NA    NA
#> 11    G1I1     2   3PLM   0.96  0.10  -0.46  0.20   0.16  0.07     NA    NA
#> 12    G1I2     2   3PLM   0.88  0.13   1.21  0.14   0.10  0.04     NA    NA
#> 13    G1I3     2   3PLM   1.49  0.26   1.31  0.09   0.18  0.03     NA    NA
#> 14    G1I4     2   3PLM   1.52  0.21   0.27  0.12   0.30  0.04     NA    NA
#> 15    G1I5     2   3PLM   1.33  0.14  -0.18  0.13   0.16  0.05     NA    NA
#> 16    G1I6     2   3PLM   2.15  0.18   0.01  0.05   0.08  0.03     NA    NA
#> 17    G1I7     2   3PLM   1.45  0.16  -0.07  0.12   0.20  0.05     NA    NA
#> 18    G1I8     2   3PLM   2.50  0.46   1.18  0.07   0.32  0.02     NA    NA
#> 19    G1I9     2   3PLM   2.40  0.26  -0.93  0.11   0.24  0.06     NA    NA
#> 20   G1I10     2   3PLM   1.26  0.13  -1.75  0.23   0.24  0.10     NA    NA
#> 21   G1I11     2   3PLM   1.55  0.15  -1.11  0.16   0.21  0.08     NA    NA
#> 22   G1I12     2   3PLM   0.75  0.09  -0.81  0.31   0.19  0.08     NA    NA
#> 23   G1I13     2   3PLM   1.04  0.14  -0.13  0.21   0.20  0.07     NA    NA
#> 24   G1I14     2   3PLM   1.46  0.37   1.77  0.15   0.30  0.03     NA    NA
#> 25   G1I15     2   3PLM   0.86  0.09  -1.44  0.29   0.21  0.09     NA    NA
#> 26   G1I16     2   3PLM   1.05  0.11  -1.95  0.27   0.23  0.10     NA    NA
#> 27   G1I17     2   3PLM   1.04  0.13   0.25  0.15   0.15  0.05     NA    NA
#> 28   G1I18     2   3PLM   2.08  0.20  -0.09  0.08   0.24  0.04     NA    NA
#> 29   G1I19     2   3PLM   1.32  0.13  -1.38  0.19   0.20  0.08     NA    NA
#> 30   G1I20     2   3PLM   1.01  0.17   0.51  0.19   0.22  0.06     NA    NA
#> 31   G1I21     2   3PLM   0.92  0.13   0.77  0.15   0.13  0.05     NA    NA
#> 32   G1I22     2   3PLM   1.77  0.21  -0.65  0.16   0.35  0.06     NA    NA
#> 33   G1I23     2   3PLM   1.31  0.15  -1.17  0.23   0.25  0.09     NA    NA
#> 34   G1I24     2   3PLM   1.64  0.19   0.32  0.10   0.23  0.04     NA    NA
#> 35   G1I25     2   3PLM   1.58  0.18  -0.13  0.12   0.25  0.05     NA    NA
#> 36   G1I26     2   3PLM   1.88  0.25   0.66  0.08   0.25  0.03     NA    NA
#> 37   G1I27     2   3PLM   1.63  0.18  -1.54  0.20   0.27  0.10     NA    NA
#> 38   G1I28     2   3PLM   1.33  0.16   0.57  0.10   0.15  0.04     NA    NA
#> 39   G1I29     2   3PLM   0.92  0.09  -0.37  0.17   0.12  0.05     NA    NA
#> 40   G1I30     2   3PLM   1.00  0.27   2.30  0.25   0.17  0.03     NA    NA
#> 41   G1I31     2   3PLM   2.39  0.46   1.64  0.09   0.18  0.01     NA    NA
#> 42   G1I32     2   3PLM   1.10  0.12  -0.08  0.15   0.15  0.05     NA    NA
#> 43   G1I33     2   3PLM   1.60  0.17   0.18  0.09   0.16  0.04     NA    NA
#> 44   G1I34     2   3PLM   1.33  0.13   0.24  0.09   0.11  0.04     NA    NA
#> 45   G1I35     2   3PLM   1.33  0.17   1.30  0.08   0.07  0.02     NA    NA
#> 46   G1I36     2   3PLM   1.44  0.14  -1.23  0.19   0.22  0.09     NA    NA
#> 47   G1I37     2   3PLM   1.05  0.14  -0.60  0.27   0.27  0.09     NA    NA
#> 48   G1I38     5    GRM   1.06  0.06  -0.37  0.05   0.22  0.05   0.87  0.06
#> 49   C1I11     5    GRM   1.23    NA  -2.08    NA  -1.35    NA  -0.71    NA
#> 50   C1I12     5    GRM   0.88    NA  -0.76    NA  -0.01    NA   0.67    NA
#> 51    G2I1     2   3PLM   1.74  0.19  -0.87  0.17   0.23  0.10     NA    NA
#> 52    G2I2     2   3PLM   0.83  0.11  -0.45  0.31   0.21  0.09     NA    NA
#> 53    G2I3     2   3PLM   1.06  0.13   0.08  0.20   0.18  0.07     NA    NA
#> 54    G2I4     2   3PLM   1.45  0.28   1.50  0.10   0.18  0.04     NA    NA
#> 55    G2I5     2   3PLM   0.70  0.10  -1.68  0.40   0.20  0.09     NA    NA
#> 56    G2I6     2   3PLM   1.08  0.14  -1.60  0.29   0.22  0.09     NA    NA
#> 57    G2I7     2   3PLM   1.30  0.13   0.27  0.12   0.13  0.05     NA    NA
#> 58    G2I8     2   3PLM   2.18  0.19  -0.09  0.08   0.14  0.05     NA    NA
#> 59    G2I9     2   3PLM   1.07  0.13  -1.57  0.28   0.21  0.09     NA    NA
#> 60   G2I10     2   3PLM   1.59  0.27   0.90  0.12   0.29  0.05     NA    NA
#> 61   G2I11     2   3PLM   0.93  0.12   0.88  0.15   0.11  0.05     NA    NA
#> 62   G2I12     2   3PLM   1.66  0.18  -0.70  0.18   0.25  0.10     NA    NA
#> 63   G2I13     2   3PLM   1.16  0.13  -1.31  0.24   0.21  0.09     NA    NA
#> 64   G2I14     2   3PLM   1.39  0.14   0.22  0.11   0.13  0.05     NA    NA
#> 65   G2I15     2   3PLM   1.53  0.20   0.03  0.17   0.26  0.08     NA    NA
#> 66   G2I16     2   3PLM   1.70  0.23   0.75  0.10   0.21  0.05     NA    NA
#> 67   G2I17     2   3PLM   2.11  0.22  -1.20  0.14   0.20  0.09     NA    NA
#> 68   G2I18     2   3PLM   1.72  0.20   0.64  0.08   0.14  0.04     NA    NA
#> 69   G2I19     2   3PLM   1.10  0.14   0.05  0.22   0.21  0.08     NA    NA
#> 70   G2I20     2   3PLM   1.62  0.45   2.24  0.19   0.16  0.02     NA    NA
#> 71   G2I21     2   3PLM   2.39  0.36   1.59  0.06   0.12  0.02     NA    NA
#> 72   G2I22     2   3PLM   1.43  0.17   0.22  0.15   0.20  0.07     NA    NA
#> 73   G2I23     2   3PLM   2.09  0.26   0.44  0.10   0.28  0.05     NA    NA
#> 74   G2I24     2   3PLM   1.40  0.12   0.37  0.08   0.07  0.03     NA    NA
#> 75   G2I25     2   3PLM   1.74  0.22   1.42  0.06   0.07  0.02     NA    NA
#> 76   G2I26     2   3PLM   1.89  0.20  -0.88  0.16   0.24  0.10     NA    NA
#> 77   G2I27     2   3PLM   0.98  0.12  -0.59  0.25   0.20  0.09     NA    NA
#> 78   G2I28     5    GRM   1.12  0.07  -0.37  0.07   0.15  0.05   0.80  0.05
#> 79    C2I1     2   3PLM   0.97    NA  -0.46    NA   0.05    NA     NA    NA
#> 80    C2I2     2   3PLM   0.85    NA   1.18    NA   0.01    NA     NA    NA
#> 81    C2I3     2   3PLM   1.43    NA   1.41    NA   0.10    NA     NA    NA
#> 82    C2I4     2   3PLM   1.48    NA   0.18    NA   0.17    NA     NA    NA
#> 83    C2I5     2   3PLM   1.27    NA  -0.23    NA   0.03    NA     NA    NA
#> 84    C2I6     2   3PLM   2.02    NA  -0.09    NA   0.01    NA     NA    NA
#> 85    C2I7     2   3PLM   1.37    NA  -0.13    NA   0.10    NA     NA    NA
#> 86    C2I8     2   3PLM   1.67    NA   1.25    NA   0.19    NA     NA    NA
#> 87    C2I9     2   3PLM   2.28    NA  -1.01    NA   0.10    NA     NA    NA
#> 88   C2I10     2   3PLM   1.42    NA  -1.65    NA   0.11    NA     NA    NA
#> 89    G3I1     2   3PLM   1.60  0.14  -1.04  0.12   0.14  0.05     NA    NA
#> 90    G3I2     2   3PLM   0.79  0.09  -0.62  0.25   0.16  0.07     NA    NA
#> 91    G3I3     2   3PLM   1.16  0.12   0.07  0.12   0.16  0.04     NA    NA
#> 92    G3I4     2   3PLM   1.30  0.24   1.72  0.12   0.22  0.02     NA    NA
#> 93    G3I5     2   3PLM   0.80  0.10  -1.09  0.33   0.24  0.09     NA    NA
#> 94    G3I6     2   3PLM   1.18  0.13  -1.38  0.24   0.26  0.08     NA    NA
#> 95    G3I7     2   3PLM   1.37  0.13   0.26  0.09   0.13  0.03     NA    NA
#> 96    G3I8     2   3PLM   1.95  0.18  -0.07  0.06   0.14  0.02     NA    NA
#> 97    G3I9     2   3PLM   1.18  0.11  -1.16  0.19   0.18  0.07     NA    NA
#> 98   G3I10     2   3PLM   1.50  0.21   0.98  0.10   0.29  0.02     NA    NA
#> 99   G3I11     2   3PLM   0.91  0.10   0.80  0.11   0.08  0.03     NA    NA
#> 100  G3I12     2   3PLM   1.44  0.17  -0.78  0.15   0.23  0.06     NA    NA
#> 101  G3I13     2   3PLM   1.15  0.13  -1.06  0.20   0.21  0.07     NA    NA
#> 102  G3I14     2   3PLM   1.38  0.14   0.31  0.09   0.16  0.03     NA    NA
#> 103  G3I15     2   3PLM   1.29  0.13  -0.06  0.11   0.18  0.04     NA    NA
#> 104  G3I16     2   3PLM   1.58  0.18   0.75  0.08   0.20  0.02     NA    NA
#> 105  G3I17     2   3PLM   1.72  0.15  -1.46  0.13   0.16  0.06     NA    NA
#> 106  G3I18     2   3PLM   1.37  0.14   0.71  0.08   0.11  0.02     NA    NA
#> 107  G3I19     2   3PLM   0.99  0.11   0.07  0.15   0.15  0.04     NA    NA
#> 108  G3I20     2   3PLM   1.14  0.23   2.30  0.20   0.12  0.02     NA    NA
#> 109  G3I21     2   3PLM   2.97  0.55   1.66  0.06   0.14  0.01     NA    NA
#> 110  G3I22     2   3PLM   1.19  0.11   0.10  0.10   0.10  0.03     NA    NA
#> 111  G3I23     2   3PLM   1.79  0.17   0.33  0.06   0.15  0.02     NA    NA
#> 112  G3I24     2   3PLM   1.15  0.09   0.34  0.08   0.05  0.02     NA    NA
#> 113  G3I25     2   3PLM   1.38  0.15   1.36  0.08   0.05  0.01     NA    NA
#> 114  G3I26     2   3PLM   1.72  0.17  -0.96  0.12   0.20  0.05     NA    NA
#> 115  G3I27     2   3PLM   0.95  0.10  -0.61  0.20   0.16  0.06     NA    NA
#> 116  G3I28     5    GRM   1.00  0.05  -0.31  0.05   0.19  0.05   0.82  0.07
#>      par.5  se.5
#> 1       NA    NA
#> 2       NA    NA
#> 3       NA    NA
#> 4       NA    NA
#> 5       NA    NA
#> 6       NA    NA
#> 7       NA    NA
#> 8       NA    NA
#> 9       NA    NA
#> 10      NA    NA
#> 11      NA    NA
#> 12      NA    NA
#> 13      NA    NA
#> 14      NA    NA
#> 15      NA    NA
#> 16      NA    NA
#> 17      NA    NA
#> 18      NA    NA
#> 19      NA    NA
#> 20      NA    NA
#> 21      NA    NA
#> 22      NA    NA
#> 23      NA    NA
#> 24      NA    NA
#> 25      NA    NA
#> 26      NA    NA
#> 27      NA    NA
#> 28      NA    NA
#> 29      NA    NA
#> 30      NA    NA
#> 31      NA    NA
#> 32      NA    NA
#> 33      NA    NA
#> 34      NA    NA
#> 35      NA    NA
#> 36      NA    NA
#> 37      NA    NA
#> 38      NA    NA
#> 39      NA    NA
#> 40      NA    NA
#> 41      NA    NA
#> 42      NA    NA
#> 43      NA    NA
#> 44      NA    NA
#> 45      NA    NA
#> 46      NA    NA
#> 47      NA    NA
#> 48    1.44  0.08
#> 49   -0.12    NA
#> 50    1.25    NA
#> 51      NA    NA
#> 52      NA    NA
#> 53      NA    NA
#> 54      NA    NA
#> 55      NA    NA
#> 56      NA    NA
#> 57      NA    NA
#> 58      NA    NA
#> 59      NA    NA
#> 60      NA    NA
#> 61      NA    NA
#> 62      NA    NA
#> 63      NA    NA
#> 64      NA    NA
#> 65      NA    NA
#> 66      NA    NA
#> 67      NA    NA
#> 68      NA    NA
#> 69      NA    NA
#> 70      NA    NA
#> 71      NA    NA
#> 72      NA    NA
#> 73      NA    NA
#> 74      NA    NA
#> 75      NA    NA
#> 76      NA    NA
#> 77      NA    NA
#> 78    1.48  0.07
#> 79      NA    NA
#> 80      NA    NA
#> 81      NA    NA
#> 82      NA    NA
#> 83      NA    NA
#> 84      NA    NA
#> 85      NA    NA
#> 86      NA    NA
#> 87      NA    NA
#> 88      NA    NA
#> 89      NA    NA
#> 90      NA    NA
#> 91      NA    NA
#> 92      NA    NA
#> 93      NA    NA
#> 94      NA    NA
#> 95      NA    NA
#> 96      NA    NA
#> 97      NA    NA
#> 98      NA    NA
#> 99      NA    NA
#> 100     NA    NA
#> 101     NA    NA
#> 102     NA    NA
#> 103     NA    NA
#> 104     NA    NA
#> 105     NA    NA
#> 106     NA    NA
#> 107     NA    NA
#> 108     NA    NA
#> 109     NA    NA
#> 110     NA    NA
#> 111     NA    NA
#> 112     NA    NA
#> 113     NA    NA
#> 114     NA    NA
#> 115     NA    NA
#> 116   1.52  0.09
#>  Group Parameters: 
#>                      mu  sigma2  sigma
#> estimate(Group1)  -0.01    1.01   1.01
#> se(Group1)         0.02    0.03   0.02
#> estimate(Group2)   0.50    0.58   0.76
#> se(Group2)         0.02    0.02   0.01
#> estimate(Group3)  -0.28    1.66   1.29
#> se(Group3)         0.03    0.05   0.02
#> 

## ------------------------------------------------------------------------------
# 3. MG calibration with FIPC using simMG data
#    (Estimate group parameters only)
#  - Details:
#    (a) Fix all item parameters across all three groups
#    (b) Freely estimate the means and variances of the ability
#        distributions for all three groups
## ------------------------------------------------------------------------------
# 3-(1). Freely estimate the means and variances for all three groups
# Set all three groups as free groups in which the scales
# of the ability distributions will be freely estimated
free.group <- 1:3 # or use 'free.group <- group.name'

# Specify the locations of all fixed items in each group's metadata
fix.loc <- list(1:50, 1:50, 1:38)

# Estimate group parameters only using MG-FIPC
fit.5 <-
  est_mg(
    x = x, data = data, group.name = group.name, D = 1,
    free.group = free.group, use.gprior = TRUE,
    gprior = list(dist = "beta", params = c(5, 16)),
    EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
    fipc.method = "MEM", fix.loc = fix.loc
  )
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -159735.7617, Max-Change: 0.476403 EM iteration: 2, Loglike: -159719.0182, Max-Change: 0.138031 EM iteration: 3, Loglike: -159716.3699, Max-Change: 0.04617 EM iteration: 4, Loglike: -159715.2802, Max-Change: 0.018371 EM iteration: 5, Loglike: -159714.6363, Max-Change: 0.009115 EM iteration: 6, Loglike: -159714.1843, Max-Change: 0.005329 EM iteration: 7, Loglike: -159713.8343, Max-Change: 0.003463 EM iteration: 8, Loglike: -159713.5471, Max-Change: 0.002427 EM iteration: 9, Loglike: -159713.3033, Max-Change: 0.00181 EM iteration: 10, Loglike: -159713.0914, Max-Change: 0.001421 EM iteration: 11, Loglike: -159712.9045, Max-Change: 0.001165 EM iteration: 12, Loglike: -159712.7377, Max-Change: 0.000987 
#> Estimation is finished in 0.88 seconds. 

# Summary of the estimation
summary(fit.5)
#> 
#> Call:
#> est_mg(x = x, data = data, group.name = group.name, D = 1, free.group = free.group, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE, 
#>     fipc.method = "MEM", fix.loc = fix.loc)
#> 
#> Summary of the Data 
#>  Number of Items: 
#>   Overall: 116 unique items 
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of Cases: 
#>   Overall: 6000
#>   By group: 2000(Group1), 2000(Group2), 2000(Group3)
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 6
#>  Number of fixed items: 
#>   Overall: 116
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of E-step cycles completed: 12
#>  Maximum parameter change: 0.0009867474
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 0.65
#>  Standard error computation: 
#>  Total computation: 0.88
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are satisfied.
#>  Second-order test: Solution is a possible local maximum.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates was not estimated.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 
#>   Overall: 319425.5
#>   By group: 120445.354(Group1), 114033.356(Group2), 84946.765(Group3)
#> 
#>  Akaike Information Criterion (AIC): 319437.5
#>  Bayesian Information Criterion (BIC): 319477.7
#>  Item Parameters (Overall): 
#>         id  cats  model  par.1  se.1  par.2  se.2  par.3  se.3  par.4  se.4
#> 1     C1I1     2   3PLM   0.76    NA   1.46    NA   0.26    NA     NA    NA
#> 2     C1I2     2   3PLM   1.92    NA  -1.05    NA   0.18    NA     NA    NA
#> 3     C1I3     2   3PLM   0.93    NA   0.39    NA   0.10    NA     NA    NA
#> 4     C1I4     2   3PLM   1.05    NA  -0.41    NA   0.20    NA     NA    NA
#> 5     C1I5     2   3PLM   0.87    NA  -0.12    NA   0.16    NA     NA    NA
#> 6     C1I6     2   3PLM   1.70    NA   0.63    NA   0.07    NA     NA    NA
#> 7     C1I7     2   3PLM   0.91    NA   1.02    NA   0.12    NA     NA    NA
#> 8     C1I8     2   3PLM   0.84    NA   0.80    NA   0.11    NA     NA    NA
#> 9     C1I9     2   3PLM   0.85    NA   0.85    NA   0.26    NA     NA    NA
#> 10   C1I10     2   3PLM   1.53    NA   0.09    NA   0.14    NA     NA    NA
#> 11    G1I1     2   3PLM   1.00    NA  -0.46    NA   0.13    NA     NA    NA
#> 12    G1I2     2   3PLM   0.88    NA   1.18    NA   0.09    NA     NA    NA
#> 13    G1I3     2   3PLM   1.46    NA   1.41    NA   0.18    NA     NA    NA
#> 14    G1I4     2   3PLM   1.51    NA   0.18    NA   0.25    NA     NA    NA
#> 15    G1I5     2   3PLM   1.30    NA  -0.23    NA   0.11    NA     NA    NA
#> 16    G1I6     2   3PLM   2.05    NA  -0.09    NA   0.05    NA     NA    NA
#> 17    G1I7     2   3PLM   1.40    NA  -0.13    NA   0.18    NA     NA    NA
#> 18    G1I8     2   3PLM   1.70    NA   1.25    NA   0.27    NA     NA    NA
#> 19    G1I9     2   3PLM   2.31    NA  -1.01    NA   0.18    NA     NA    NA
#> 20   G1I10     2   3PLM   1.45    NA  -1.65    NA   0.19    NA     NA    NA
#> 21   G1I11     2   3PLM   1.63    NA  -1.19    NA   0.12    NA     NA    NA
#> 22   G1I12     2   3PLM   0.83    NA  -0.68    NA   0.20    NA     NA    NA
#> 23   G1I13     2   3PLM   0.98    NA  -0.26    NA   0.13    NA     NA    NA
#> 24   G1I14     2   3PLM   1.14    NA   1.68    NA   0.25    NA     NA    NA
#> 25   G1I15     2   3PLM   0.79    NA  -1.39    NA   0.26    NA     NA    NA
#> 26   G1I16     2   3PLM   1.09    NA  -1.85    NA   0.17    NA     NA    NA
#> 27   G1I17     2   3PLM   1.17    NA   0.07    NA   0.13    NA     NA    NA
#> 28   G1I18     2   3PLM   2.15    NA  -0.09    NA   0.21    NA     NA    NA
#> 29   G1I19     2   3PLM   1.28    NA  -1.38    NA   0.20    NA     NA    NA
#> 30   G1I20     2   3PLM   1.35    NA   0.82    NA   0.32    NA     NA    NA
#> 31   G1I21     2   3PLM   0.82    NA   0.71    NA   0.08    NA     NA    NA
#> 32   G1I22     2   3PLM   1.52    NA  -0.89    NA   0.26    NA     NA    NA
#> 33   G1I23     2   3PLM   1.27    NA  -1.31    NA   0.19    NA     NA    NA
#> 34   G1I24     2   3PLM   1.31    NA   0.19    NA   0.16    NA     NA    NA
#> 35   G1I25     2   3PLM   1.47    NA  -0.14    NA   0.23    NA     NA    NA
#> 36   G1I26     2   3PLM   1.47    NA   0.64    NA   0.23    NA     NA    NA
#> 37   G1I27     2   3PLM   1.76    NA  -1.53    NA   0.16    NA     NA    NA
#> 38   G1I28     2   3PLM   1.44    NA   0.54    NA   0.14    NA     NA    NA
#> 39   G1I29     2   3PLM   0.98    NA  -0.37    NA   0.13    NA     NA    NA
#> 40   G1I30     2   3PLM   0.99    NA   2.37    NA   0.16    NA     NA    NA
#> 41   G1I31     2   3PLM   2.27    NA   1.62    NA   0.18    NA     NA    NA
#> 42   G1I32     2   3PLM   1.23    NA  -0.07    NA   0.13    NA     NA    NA
#> 43   G1I33     2   3PLM   1.64    NA   0.17    NA   0.18    NA     NA    NA
#> 44   G1I34     2   3PLM   1.21    NA   0.24    NA   0.08    NA     NA    NA
#> 45   G1I35     2   3PLM   1.32    NA   1.34    NA   0.08    NA     NA    NA
#> 46   G1I36     2   3PLM   1.74    NA  -1.00    NA   0.25    NA     NA    NA
#> 47   G1I37     2   3PLM   0.97    NA  -0.73    NA   0.22    NA     NA    NA
#> 48   G1I38     5    GRM   1.14    NA  -0.37    NA   0.22    NA   0.85    NA
#> 49   C1I11     5    GRM   1.23    NA  -2.08    NA  -1.35    NA  -0.71    NA
#> 50   C1I12     5    GRM   0.88    NA  -0.76    NA  -0.01    NA   0.67    NA
#> 51    G2I1     2   3PLM   1.60    NA  -1.19    NA   0.04    NA     NA    NA
#> 52    G2I2     2   3PLM   0.80    NA  -0.68    NA   0.12    NA     NA    NA
#> 53    G2I3     2   3PLM   0.95    NA  -0.26    NA   0.05    NA     NA    NA
#> 54    G2I4     2   3PLM   1.11    NA   1.68    NA   0.17    NA     NA    NA
#> 55    G2I5     2   3PLM   0.76    NA  -1.39    NA   0.18    NA     NA    NA
#> 56    G2I6     2   3PLM   1.06    NA  -1.85    NA   0.09    NA     NA    NA
#> 57    G2I7     2   3PLM   1.14    NA   0.07    NA   0.05    NA     NA    NA
#> 58    G2I8     2   3PLM   2.12    NA  -0.09    NA   0.13    NA     NA    NA
#> 59    G2I9     2   3PLM   1.25    NA  -1.38    NA   0.12    NA     NA    NA
#> 60   G2I10     2   3PLM   1.32    NA   0.82    NA   0.24    NA     NA    NA
#> 61   G2I11     2   3PLM   0.79    NA   0.71    NA   0.01    NA     NA    NA
#> 62   G2I12     2   3PLM   1.49    NA  -0.89    NA   0.18    NA     NA    NA
#> 63   G2I13     2   3PLM   1.24    NA  -1.31    NA   0.11    NA     NA    NA
#> 64   G2I14     2   3PLM   1.28    NA   0.19    NA   0.08    NA     NA    NA
#> 65   G2I15     2   3PLM   1.44    NA  -0.14    NA   0.15    NA     NA    NA
#> 66   G2I16     2   3PLM   1.44    NA   0.64    NA   0.15    NA     NA    NA
#> 67   G2I17     2   3PLM   1.73    NA  -1.53    NA   0.08    NA     NA    NA
#> 68   G2I18     2   3PLM   1.41    NA   0.54    NA   0.06    NA     NA    NA
#> 69   G2I19     2   3PLM   0.95    NA  -0.37    NA   0.05    NA     NA    NA
#> 70   G2I20     2   3PLM   0.96    NA   2.37    NA   0.08    NA     NA    NA
#> 71   G2I21     2   3PLM   2.24    NA   1.62    NA   0.10    NA     NA    NA
#> 72   G2I22     2   3PLM   1.20    NA  -0.07    NA   0.05    NA     NA    NA
#> 73   G2I23     2   3PLM   1.61    NA   0.17    NA   0.10    NA     NA    NA
#> 74   G2I24     2   3PLM   1.18    NA   0.24    NA   0.01    NA     NA    NA
#> 75   G2I25     2   3PLM   1.29    NA   1.34    NA   0.01    NA     NA    NA
#> 76   G2I26     2   3PLM   1.71    NA  -1.00    NA   0.17    NA     NA    NA
#> 77   G2I27     2   3PLM   0.94    NA  -0.73    NA   0.14    NA     NA    NA
#> 78   G2I28     5    GRM   1.11    NA  -0.37    NA   0.14    NA   0.78    NA
#> 79    C2I1     2   3PLM   0.97    NA  -0.46    NA   0.05    NA     NA    NA
#> 80    C2I2     2   3PLM   0.85    NA   1.18    NA   0.01    NA     NA    NA
#> 81    C2I3     2   3PLM   1.43    NA   1.41    NA   0.10    NA     NA    NA
#> 82    C2I4     2   3PLM   1.48    NA   0.18    NA   0.17    NA     NA    NA
#> 83    C2I5     2   3PLM   1.27    NA  -0.23    NA   0.03    NA     NA    NA
#> 84    C2I6     2   3PLM   2.02    NA  -0.09    NA   0.01    NA     NA    NA
#> 85    C2I7     2   3PLM   1.37    NA  -0.13    NA   0.10    NA     NA    NA
#> 86    C2I8     2   3PLM   1.67    NA   1.25    NA   0.19    NA     NA    NA
#> 87    C2I9     2   3PLM   2.28    NA  -1.01    NA   0.10    NA     NA    NA
#> 88   C2I10     2   3PLM   1.42    NA  -1.65    NA   0.11    NA     NA    NA
#> 89    G3I1     2   3PLM   1.55    NA  -1.12    NA   0.07    NA     NA    NA
#> 90    G3I2     2   3PLM   0.75    NA  -0.61    NA   0.15    NA     NA    NA
#> 91    G3I3     2   3PLM   0.90    NA  -0.19    NA   0.08    NA     NA    NA
#> 92    G3I4     2   3PLM   1.06    NA   1.75    NA   0.20    NA     NA    NA
#> 93    G3I5     2   3PLM   0.71    NA  -1.32    NA   0.21    NA     NA    NA
#> 94    G3I6     2   3PLM   1.01    NA  -1.78    NA   0.12    NA     NA    NA
#> 95    G3I7     2   3PLM   1.09    NA   0.14    NA   0.08    NA     NA    NA
#> 96    G3I8     2   3PLM   2.07    NA  -0.02    NA   0.16    NA     NA    NA
#> 97    G3I9     2   3PLM   1.20    NA  -1.31    NA   0.15    NA     NA    NA
#> 98   G3I10     2   3PLM   1.27    NA   0.89    NA   0.27    NA     NA    NA
#> 99   G3I11     2   3PLM   0.74    NA   0.78    NA   0.04    NA     NA    NA
#> 100  G3I12     2   3PLM   1.44    NA  -0.82    NA   0.21    NA     NA    NA
#> 101  G3I13     2   3PLM   1.19    NA  -1.24    NA   0.14    NA     NA    NA
#> 102  G3I14     2   3PLM   1.23    NA   0.26    NA   0.11    NA     NA    NA
#> 103  G3I15     2   3PLM   1.39    NA  -0.07    NA   0.18    NA     NA    NA
#> 104  G3I16     2   3PLM   1.39    NA   0.71    NA   0.18    NA     NA    NA
#> 105  G3I17     2   3PLM   1.68    NA  -1.46    NA   0.11    NA     NA    NA
#> 106  G3I18     2   3PLM   1.36    NA   0.61    NA   0.09    NA     NA    NA
#> 107  G3I19     2   3PLM   0.90    NA  -0.30    NA   0.08    NA     NA    NA
#> 108  G3I20     2   3PLM   0.91    NA   2.44    NA   0.11    NA     NA    NA
#> 109  G3I21     2   3PLM   2.19    NA   1.69    NA   0.13    NA     NA    NA
#> 110  G3I22     2   3PLM   1.15    NA   0.00    NA   0.08    NA     NA    NA
#> 111  G3I23     2   3PLM   1.56    NA   0.24    NA   0.13    NA     NA    NA
#> 112  G3I24     2   3PLM   1.13    NA   0.31    NA   0.04    NA     NA    NA
#> 113  G3I25     2   3PLM   1.24    NA   1.41    NA   0.04    NA     NA    NA
#> 114  G3I26     2   3PLM   1.66    NA  -0.93    NA   0.20    NA     NA    NA
#> 115  G3I27     2   3PLM   0.89    NA  -0.66    NA   0.17    NA     NA    NA
#> 116  G3I28     5    GRM   1.06    NA  -0.30    NA   0.17    NA   0.78    NA
#>      par.5  se.5
#> 1       NA    NA
#> 2       NA    NA
#> 3       NA    NA
#> 4       NA    NA
#> 5       NA    NA
#> 6       NA    NA
#> 7       NA    NA
#> 8       NA    NA
#> 9       NA    NA
#> 10      NA    NA
#> 11      NA    NA
#> 12      NA    NA
#> 13      NA    NA
#> 14      NA    NA
#> 15      NA    NA
#> 16      NA    NA
#> 17      NA    NA
#> 18      NA    NA
#> 19      NA    NA
#> 20      NA    NA
#> 21      NA    NA
#> 22      NA    NA
#> 23      NA    NA
#> 24      NA    NA
#> 25      NA    NA
#> 26      NA    NA
#> 27      NA    NA
#> 28      NA    NA
#> 29      NA    NA
#> 30      NA    NA
#> 31      NA    NA
#> 32      NA    NA
#> 33      NA    NA
#> 34      NA    NA
#> 35      NA    NA
#> 36      NA    NA
#> 37      NA    NA
#> 38      NA    NA
#> 39      NA    NA
#> 40      NA    NA
#> 41      NA    NA
#> 42      NA    NA
#> 43      NA    NA
#> 44      NA    NA
#> 45      NA    NA
#> 46      NA    NA
#> 47      NA    NA
#> 48    1.38    NA
#> 49   -0.12    NA
#> 50    1.25    NA
#> 51      NA    NA
#> 52      NA    NA
#> 53      NA    NA
#> 54      NA    NA
#> 55      NA    NA
#> 56      NA    NA
#> 57      NA    NA
#> 58      NA    NA
#> 59      NA    NA
#> 60      NA    NA
#> 61      NA    NA
#> 62      NA    NA
#> 63      NA    NA
#> 64      NA    NA
#> 65      NA    NA
#> 66      NA    NA
#> 67      NA    NA
#> 68      NA    NA
#> 69      NA    NA
#> 70      NA    NA
#> 71      NA    NA
#> 72      NA    NA
#> 73      NA    NA
#> 74      NA    NA
#> 75      NA    NA
#> 76      NA    NA
#> 77      NA    NA
#> 78    1.44    NA
#> 79      NA    NA
#> 80      NA    NA
#> 81      NA    NA
#> 82      NA    NA
#> 83      NA    NA
#> 84      NA    NA
#> 85      NA    NA
#> 86      NA    NA
#> 87      NA    NA
#> 88      NA    NA
#> 89      NA    NA
#> 90      NA    NA
#> 91      NA    NA
#> 92      NA    NA
#> 93      NA    NA
#> 94      NA    NA
#> 95      NA    NA
#> 96      NA    NA
#> 97      NA    NA
#> 98      NA    NA
#> 99      NA    NA
#> 100     NA    NA
#> 101     NA    NA
#> 102     NA    NA
#> 103     NA    NA
#> 104     NA    NA
#> 105     NA    NA
#> 106     NA    NA
#> 107     NA    NA
#> 108     NA    NA
#> 109     NA    NA
#> 110     NA    NA
#> 111     NA    NA
#> 112     NA    NA
#> 113     NA    NA
#> 114     NA    NA
#> 115     NA    NA
#> 116   1.39    NA
#>  Group Parameters: 
#>                      mu  sigma2  sigma
#> estimate(Group1)   0.01    0.96   0.98
#> se(Group1)         0.02    0.03   0.02
#> estimate(Group2)   0.52    0.60   0.78
#> se(Group2)         0.02    0.02   0.01
#> estimate(Group3)  -0.28    1.63   1.28
#> se(Group3)         0.03    0.05   0.02
#> 

# Extract the group parameter estimates (i.e., scale parameters)
getirt(fit.5, what = "group.par")
#> $Group1
#>                    mu     sigma2      sigma
#> estimates 0.005021632 0.96136149 0.98049044
#> se        0.021924433 0.03040852 0.01550679
#> 
#> $Group2
#>                   mu     sigma2      sigma
#> estimates 0.51580977 0.60138258 0.77548861
#> se        0.01734045 0.01902214 0.01226462
#> 
#> $Group3
#>                    mu     sigma2      sigma
#> estimates -0.27507220 1.63376235 1.27818714
#> se         0.02858113 0.05167702 0.02021497
#> 

## ------------------------------------------------------------------------------
# 4. MG calibration with FIPC using simMG data
#    (Fix only the unique items of Group 1)
#  - Details:
#    (a) Fix item parameters of the unique items in Group 1 only
#    (b) Constrain the common items across groups to have
#        the same item parameters (i.e., C1I1–C1I12 between
#        Groups 1 and 2, and C2I1–C2I10 between Groups 2 and 3)
#    (c) Freely estimate the means and variances of the ability
#        distributions for all three groups
## ------------------------------------------------------------------------------
# 4-(1). Freely estimate the means and variances for all three groups
# Set all three groups as free groups in which the scales
# of the ability distributions will be freely estimated
free.group <- group.name # or use 'free.group <- 1:3'

# Specify the item IDs of the unique items in Group 1 to be fixed using
# the `fix.id` argument.
fix.id <- paste0("G1I", 1:38)

# Alternatively, use the 'fix.loc' argument as
# 'fix.loc = list(11:48, NULL, NULL)'

# Estimate IRT parameters using MG-FIPC
fit.6 <-
  est_mg(
    x = x, data = data, group.name = group.name, D = 1,
    free.group = free.group, use.gprior = TRUE,
    gprior = list(dist = "beta", params = c(5, 16)),
    EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
    fipc.method = "MEM", fix.loc = NULL, fix.id = fix.id
  )
#> Parsing input... 
#> Estimating item parameters... 
#>  EM iteration: 1, Loglike: -43183.7384, Max-Change: 432.57651 EM iteration: 2, Loglike: -171691.1579, Max-Change: 4501.79291 EM iteration: 3, Loglike: -167102.6836, Max-Change: 21514.29256 EM iteration: 4, Loglike: -161125.0804, Max-Change: 0.614658 EM iteration: 5, Loglike: -160697.8752, Max-Change: 15745.96473 EM iteration: 6, Loglike: -160607.1903, Max-Change: 0.223705 EM iteration: 7, Loglike: -160564.3088, Max-Change: 44701.60439 EM iteration: 8, Loglike: -160537.2418, Max-Change: 0.114086 EM iteration: 9, Loglike: -160522.5507, Max-Change: 0.085219 EM iteration: 10, Loglike: -160512.8183, Max-Change: 0.065365 EM iteration: 11, Loglike: -160505.9064, Max-Change: 0.051414 EM iteration: 12, Loglike: -160500.7067, Max-Change: 0.041546 EM iteration: 13, Loglike: -160496.6073, Max-Change: 0.037669 EM iteration: 14, Loglike: -160493.2514, Max-Change: 0.034155 EM iteration: 15, Loglike: -160490.4210, Max-Change: 0.031016 EM iteration: 16, Loglike: -160487.9770, Max-Change: 0.028238 EM iteration: 17, Loglike: -160485.8270, Max-Change: 0.025795 EM iteration: 18, Loglike: -160483.9076, Max-Change: 0.024255 EM iteration: 19, Loglike: -160482.1744, Max-Change: 0.023508 EM iteration: 20, Loglike: -160480.5948, Max-Change: 0.022739 EM iteration: 21, Loglike: -160479.1451, Max-Change: 0.021973 EM iteration: 22, Loglike: -160477.8067, Max-Change: 0.021227 EM iteration: 23, Loglike: -160476.5657, Max-Change: 0.02051 EM iteration: 24, Loglike: -160475.4104, Max-Change: 0.019825 EM iteration: 25, Loglike: -160474.3324, Max-Change: 0.019174 EM iteration: 26, Loglike: -160473.3235, Max-Change: 0.018558 EM iteration: 27, Loglike: -160472.3778, Max-Change: 0.017975 EM iteration: 28, Loglike: -160471.4899, Max-Change: 0.017422 EM iteration: 29, Loglike: -160470.6550, Max-Change: 0.016898 EM iteration: 30, Loglike: -160469.8690, Max-Change: 0.016399 EM iteration: 31, Loglike: -160469.1285, Max-Change: 0.015925 EM iteration: 32, Loglike: -160468.4300, Max-Change: 0.015472 EM iteration: 33, Loglike: -160467.7707, Max-Change: 0.015038 EM iteration: 34, Loglike: -160467.1479, Max-Change: 0.014622 EM iteration: 35, Loglike: -160466.5592, Max-Change: 0.014223 EM iteration: 36, Loglike: -160466.0024, Max-Change: 0.013838 EM iteration: 37, Loglike: -160465.4753, Max-Change: 0.013467 EM iteration: 38, Loglike: -160464.9761, Max-Change: 0.013108 EM iteration: 39, Loglike: -160464.5030, Max-Change: 0.012761 EM iteration: 40, Loglike: -160464.0543, Max-Change: 0.012424 EM iteration: 41, Loglike: -160463.6286, Max-Change: 0.012098 EM iteration: 42, Loglike: -160463.2244, Max-Change: 0.011782 EM iteration: 43, Loglike: -160462.8404, Max-Change: 0.011475 EM iteration: 44, Loglike: -160462.4753, Max-Change: 0.011176 EM iteration: 45, Loglike: -160462.1280, Max-Change: 0.010886 EM iteration: 46, Loglike: -160461.7973, Max-Change: 0.010604 EM iteration: 47, Loglike: -160461.4822, Max-Change: 0.01033 EM iteration: 48, Loglike: -160461.1818, Max-Change: 0.010063 EM iteration: 49, Loglike: -160460.8951, Max-Change: 0.009804 EM iteration: 50, Loglike: -160460.6215, Max-Change: 0.009551 EM iteration: 51, Loglike: -160460.3601, Max-Change: 0.009306 EM iteration: 52, Loglike: -160460.1103, Max-Change: 0.009067 EM iteration: 53, Loglike: -160459.8716, Max-Change: 0.008836 EM iteration: 54, Loglike: -160459.6434, Max-Change: 0.00861 EM iteration: 55, Loglike: -160459.4254, Max-Change: 0.008391 EM iteration: 56, Loglike: -160459.2169, Max-Change: 0.008178 EM iteration: 57, Loglike: -160459.0178, Max-Change: 0.00797 EM iteration: 58, Loglike: -160458.8274, Max-Change: 0.007769 EM iteration: 59, Loglike: -160458.6455, Max-Change: 0.007573 EM iteration: 60, Loglike: -160458.4716, Max-Change: 0.007383 EM iteration: 61, Loglike: -160458.3053, Max-Change: 0.007199 EM iteration: 62, Loglike: -160458.1463, Max-Change: 0.007019 EM iteration: 63, Loglike: -160457.9940, Max-Change: 0.006845 EM iteration: 64, Loglike: -160457.8483, Max-Change: 0.006676 EM iteration: 65, Loglike: -160457.7087, Max-Change: 0.006512 EM iteration: 66, Loglike: -160457.5748, Max-Change: 0.006352 EM iteration: 67, Loglike: -160457.4465, Max-Change: 0.006198 EM iteration: 68, Loglike: -160457.3234, Max-Change: 0.006047 EM iteration: 69, Loglike: -160457.2053, Max-Change: 0.005901 EM iteration: 70, Loglike: -160457.0919, Max-Change: 0.00576 EM iteration: 71, Loglike: -160456.9830, Max-Change: 0.005623 EM iteration: 72, Loglike: -160456.8784, Max-Change: 0.005489 EM iteration: 73, Loglike: -160456.7779, Max-Change: 0.00536 EM iteration: 74, Loglike: -160456.6812, Max-Change: 0.005234 EM iteration: 75, Loglike: -160456.5883, Max-Change: 0.005112 EM iteration: 76, Loglike: -160456.4989, Max-Change: 0.004994 EM iteration: 77, Loglike: -160456.4129, Max-Change: 0.004879 EM iteration: 78, Loglike: -160456.3300, Max-Change: 0.004767 EM iteration: 79, Loglike: -160456.2503, Max-Change: 0.004658 EM iteration: 80, Loglike: -160456.1735, Max-Change: 0.004553 EM iteration: 81, Loglike: -160456.0995, Max-Change: 0.00445 EM iteration: 82, Loglike: -160456.0281, Max-Change: 0.004351 EM iteration: 83, Loglike: -160455.9593, Max-Change: 0.004254 EM iteration: 84, Loglike: -160455.8929, Max-Change: 0.00416 EM iteration: 85, Loglike: -160455.8289, Max-Change: 0.004069 EM iteration: 86, Loglike: -160455.7671, Max-Change: 0.00398 EM iteration: 87, Loglike: -160455.7074, Max-Change: 0.003893 EM iteration: 88, Loglike: -160455.6498, Max-Change: 0.003809 EM iteration: 89, Loglike: -160455.5941, Max-Change: 0.003727 EM iteration: 90, Loglike: -160455.5403, Max-Change: 0.003648 EM iteration: 91, Loglike: -160455.4883, Max-Change: 0.00357 EM iteration: 92, Loglike: -160455.4381, Max-Change: 0.003495 EM iteration: 93, Loglike: -160455.3894, Max-Change: 0.003421 EM iteration: 94, Loglike: -160455.3423, Max-Change: 0.00335 EM iteration: 95, Loglike: -160455.2968, Max-Change: 0.00328 EM iteration: 96, Loglike: -160455.2527, Max-Change: 0.003212 EM iteration: 97, Loglike: -160455.2100, Max-Change: 0.003146 EM iteration: 98, Loglike: -160455.1686, Max-Change: 0.003082 EM iteration: 99, Loglike: -160455.1286, Max-Change: 0.003019 EM iteration: 100, Loglike: -160455.0897, Max-Change: 0.002958 EM iteration: 101, Loglike: -160455.0520, Max-Change: 0.002898 EM iteration: 102, Loglike: -160455.0155, Max-Change: 0.00284 EM iteration: 103, Loglike: -160454.9800, Max-Change: 0.002783 EM iteration: 104, Loglike: -160454.9456, Max-Change: 0.002728 EM iteration: 105, Loglike: -160454.9121, Max-Change: 0.002674 EM iteration: 106, Loglike: -160454.8797, Max-Change: 0.002622 EM iteration: 107, Loglike: -160454.8482, Max-Change: 0.00257 EM iteration: 108, Loglike: -160454.8176, Max-Change: 0.00252 EM iteration: 109, Loglike: -160454.7878, Max-Change: 0.002471 EM iteration: 110, Loglike: -160454.7588, Max-Change: 0.002424 EM iteration: 111, Loglike: -160454.7307, Max-Change: 0.002377 EM iteration: 112, Loglike: -160454.7033, Max-Change: 0.002332 EM iteration: 113, Loglike: -160454.6767, Max-Change: 0.002287 EM iteration: 114, Loglike: -160454.6507, Max-Change: 0.002244 EM iteration: 115, Loglike: -160454.6255, Max-Change: 0.002201 EM iteration: 116, Loglike: -160454.6009, Max-Change: 0.00216 EM iteration: 117, Loglike: -160454.5769, Max-Change: 0.002119 EM iteration: 118, Loglike: -160454.5536, Max-Change: 0.00208 EM iteration: 119, Loglike: -160454.5308, Max-Change: 0.002041 EM iteration: 120, Loglike: -160454.5086, Max-Change: 0.002003 EM iteration: 121, Loglike: -160454.4869, Max-Change: 0.001967 EM iteration: 122, Loglike: -160454.4658, Max-Change: 0.001931 EM iteration: 123, Loglike: -160454.4452, Max-Change: 0.001895 EM iteration: 124, Loglike: -160454.4251, Max-Change: 0.001861 EM iteration: 125, Loglike: -160454.4054, Max-Change: 0.001827 EM iteration: 126, Loglike: -160454.3862, Max-Change: 0.001794 EM iteration: 127, Loglike: -160454.3674, Max-Change: 0.001762 EM iteration: 128, Loglike: -160454.3491, Max-Change: 0.00173 EM iteration: 129, Loglike: -160454.3312, Max-Change: 0.001699 EM iteration: 130, Loglike: -160454.3137, Max-Change: 0.001669 EM iteration: 131, Loglike: -160454.2966, Max-Change: 0.00164 EM iteration: 132, Loglike: -160454.2798, Max-Change: 0.001611 EM iteration: 133, Loglike: -160454.2634, Max-Change: 0.001582 EM iteration: 134, Loglike: -160454.2473, Max-Change: 0.001555 EM iteration: 135, Loglike: -160454.2316, Max-Change: 0.001528 EM iteration: 136, Loglike: -160454.2162, Max-Change: 0.001501 EM iteration: 137, Loglike: -160454.2011, Max-Change: 0.001475 EM iteration: 138, Loglike: -160454.1863, Max-Change: 0.00145 EM iteration: 139, Loglike: -160454.1718, Max-Change: 0.001425 EM iteration: 140, Loglike: -160454.1577, Max-Change: 0.00140 EM iteration: 141, Loglike: -160454.1437, Max-Change: 0.001377 EM iteration: 142, Loglike: -160454.1301, Max-Change: 0.001353 EM iteration: 143, Loglike: -160454.1166, Max-Change: 0.00133 EM iteration: 144, Loglike: -160454.1035, Max-Change: 0.001308 EM iteration: 145, Loglike: -160454.0906, Max-Change: 0.001286 EM iteration: 146, Loglike: -160454.0779, Max-Change: 0.001265 EM iteration: 147, Loglike: -160454.0655, Max-Change: 0.001243 EM iteration: 148, Loglike: -160454.0533, Max-Change: 0.001223 EM iteration: 149, Loglike: -160454.0413, Max-Change: 0.001203 EM iteration: 150, Loglike: -160454.0295, Max-Change: 0.001183 EM iteration: 151, Loglike: -160454.0179, Max-Change: 0.001164 EM iteration: 152, Loglike: -160454.0065, Max-Change: 0.001145 EM iteration: 153, Loglike: -160453.9953, Max-Change: 0.001126 EM iteration: 154, Loglike: -160453.9842, Max-Change: 0.001108 EM iteration: 155, Loglike: -160453.9734, Max-Change: 0.00109 EM iteration: 156, Loglike: -160453.9627, Max-Change: 0.001072 EM iteration: 157, Loglike: -160453.9522, Max-Change: 0.001055 EM iteration: 158, Loglike: -160453.9419, Max-Change: 0.001038 EM iteration: 159, Loglike: -160453.9317, Max-Change: 0.001022 EM iteration: 160, Loglike: -160453.9217, Max-Change: 0.001006 EM iteration: 161, Loglike: -160453.9118, Max-Change: 0.00099 
#> Warning: Convergence criteria are not satisfied. 
#> Computing item parameter var-covariance matrix... 
#> Estimation is finished in 59.13 seconds. 

# Summary of the estimation
summary(fit.6)
#> 
#> Call:
#> est_mg(x = x, data = data, group.name = group.name, D = 1, free.group = free.group, 
#>     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 
#>         16)), EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE, 
#>     fipc.method = "MEM", fix.loc = NULL, fix.id = fix.id)
#> 
#> Summary of the Data 
#>  Number of Items: 
#>   Overall: 116 unique items 
#>   By group: 50(Group1), 50(Group2), 38(Group3)
#>  Number of Cases: 
#>   Overall: 6000
#>   By group: 2000(Group1), 2000(Group2), 2000(Group3)
#> 
#> Summary of Estimation Process 
#>  Maximum number of EM cycles: 500
#>  Convergence criterion of E-step: 0.001
#>  Number of rectangular quadrature points: 49
#>  Minimum & Maximum quadrature points: -6, 6
#>  Number of free parameters: 248
#>  Number of fixed items: 
#>   Overall: 38
#>   By group: 38(Group1), 0(Group2), 0(Group3)
#>  Number of E-step cycles completed: 161
#>  Maximum parameter change: 0.0009898914
#> 
#> Processing time (in seconds) 
#>  EM algorithm: 58.47
#>  Standard error computation: 0.13
#>  Total computation: 59.13
#> 
#> Convergence and Stability of Solution 
#>  First-order test: Convergence criteria are not satisfied.
#>  Second-order test: Information matrix of item parameter estimates is positive definite.
#>  Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Summary of Estimation Results 
#>  -2loglikelihood: 
#>   Overall: 320907.8
#>   By group: 120433.219(Group1), 114532.034(Group2), 85942.567(Group3)
#> 
#>  Akaike Information Criterion (AIC): 321403.8
#>  Bayesian Information Criterion (BIC): 323065.3
#>  Item Parameters (Overall): 
#>         id  cats  model  par.1  se.1     par.2      se.2  par.3  se.3  par.4
#> 1     C1I1     2   3PLM   0.88  0.17      1.37      0.17   0.27  0.05     NA
#> 2     C1I2     2   3PLM   2.17  0.14     -0.97      0.09   0.17  0.06     NA
#> 3     C1I3     2   3PLM   1.05  0.12      0.63      0.13   0.18  0.05     NA
#> 4     C1I4     2   3PLM   1.05  0.12     -0.24      0.22   0.25  0.07     NA
#> 5     C1I5     2   3PLM   0.87  0.09     -0.16      0.21   0.17  0.07     NA
#> 6     C1I6     2   3PLM   1.87  0.13      0.61      0.04   0.08  0.02     NA
#> 7     C1I7     2   3PLM   1.07  0.13      1.12      0.09   0.14  0.03     NA
#> 8     C1I8     2   3PLM   0.93  0.12      0.90      0.14   0.16  0.05     NA
#> 9     C1I9     2   3PLM   0.89  0.12      0.63      0.19   0.20  0.06     NA
#> 10   C1I10     2   3PLM   1.46  0.11      0.13      0.09   0.14  0.04     NA
#> 11    G1I1     2   3PLM   1.00    NA     -0.46        NA   0.13    NA     NA
#> 12    G1I2     2   3PLM   0.88    NA      1.18        NA   0.09    NA     NA
#> 13    G1I3     2   3PLM   1.46    NA      1.41        NA   0.18    NA     NA
#> 14    G1I4     2   3PLM   1.51    NA      0.18        NA   0.25    NA     NA
#> 15    G1I5     2   3PLM   1.30    NA     -0.23        NA   0.11    NA     NA
#> 16    G1I6     2   3PLM   2.05    NA     -0.09        NA   0.05    NA     NA
#> 17    G1I7     2   3PLM   1.40    NA     -0.13        NA   0.18    NA     NA
#> 18    G1I8     2   3PLM   1.70    NA      1.25        NA   0.27    NA     NA
#> 19    G1I9     2   3PLM   2.31    NA     -1.01        NA   0.18    NA     NA
#> 20   G1I10     2   3PLM   1.45    NA     -1.65        NA   0.19    NA     NA
#> 21   G1I11     2   3PLM   1.63    NA     -1.19        NA   0.12    NA     NA
#> 22   G1I12     2   3PLM   0.83    NA     -0.68        NA   0.20    NA     NA
#> 23   G1I13     2   3PLM   0.98    NA     -0.26        NA   0.13    NA     NA
#> 24   G1I14     2   3PLM   1.14    NA      1.68        NA   0.25    NA     NA
#> 25   G1I15     2   3PLM   0.79    NA     -1.39        NA   0.26    NA     NA
#> 26   G1I16     2   3PLM   1.09    NA     -1.85        NA   0.17    NA     NA
#> 27   G1I17     2   3PLM   1.17    NA      0.07        NA   0.13    NA     NA
#> 28   G1I18     2   3PLM   2.15    NA     -0.09        NA   0.21    NA     NA
#> 29   G1I19     2   3PLM   1.28    NA     -1.38        NA   0.20    NA     NA
#> 30   G1I20     2   3PLM   1.35    NA      0.82        NA   0.32    NA     NA
#> 31   G1I21     2   3PLM   0.82    NA      0.71        NA   0.08    NA     NA
#> 32   G1I22     2   3PLM   1.52    NA     -0.89        NA   0.26    NA     NA
#> 33   G1I23     2   3PLM   1.27    NA     -1.31        NA   0.19    NA     NA
#> 34   G1I24     2   3PLM   1.31    NA      0.19        NA   0.16    NA     NA
#> 35   G1I25     2   3PLM   1.47    NA     -0.14        NA   0.23    NA     NA
#> 36   G1I26     2   3PLM   1.47    NA      0.64        NA   0.23    NA     NA
#> 37   G1I27     2   3PLM   1.76    NA     -1.53        NA   0.16    NA     NA
#> 38   G1I28     2   3PLM   1.44    NA      0.54        NA   0.14    NA     NA
#> 39   G1I29     2   3PLM   0.98    NA     -0.37        NA   0.13    NA     NA
#> 40   G1I30     2   3PLM   0.99    NA      2.37        NA   0.16    NA     NA
#> 41   G1I31     2   3PLM   2.27    NA      1.62        NA   0.18    NA     NA
#> 42   G1I32     2   3PLM   1.23    NA     -0.07        NA   0.13    NA     NA
#> 43   G1I33     2   3PLM   1.64    NA      0.17        NA   0.18    NA     NA
#> 44   G1I34     2   3PLM   1.21    NA      0.24        NA   0.08    NA     NA
#> 45   G1I35     2   3PLM   1.32    NA      1.34        NA   0.08    NA     NA
#> 46   G1I36     2   3PLM   1.74    NA     -1.00        NA   0.25    NA     NA
#> 47   G1I37     2   3PLM   0.97    NA     -0.73        NA   0.22    NA     NA
#> 48   G1I38     5    GRM   1.14    NA     -0.37        NA   0.22    NA   0.85
#> 49   C1I11     5    GRM   1.21  0.05     -2.13      0.09  -1.39  0.06  -0.70
#> 50   C1I12     5    GRM   0.92  0.04     -0.66      0.05   0.05  0.04   0.70
#> 51    G2I1     2   3PLM   1.76  0.19     -0.85      0.17   0.22  0.09     NA
#> 52    G2I2     2   3PLM   0.00  0.09  21591.11  99999.00   0.21  0.09     NA
#> 53    G2I3     2   3PLM   1.08  0.14      0.12      0.20   0.18  0.08     NA
#> 54    G2I4     2   3PLM   1.53  0.30      1.47      0.09   0.18  0.04     NA
#> 55    G2I5     2   3PLM   0.00  0.09  16153.77  99999.00   0.21  0.09     NA
#> 56    G2I6     2   3PLM   1.07  0.15     -1.61      0.31   0.22  0.10     NA
#> 57    G2I7     2   3PLM   1.33  0.14      0.29      0.12   0.13  0.05     NA
#> 58    G2I8     2   3PLM   2.24  0.20     -0.06      0.08   0.14  0.05     NA
#> 59    G2I9     2   3PLM   1.08  0.13     -1.54      0.28   0.21  0.09     NA
#> 60   G2I10     2   3PLM   1.77  0.30      0.95      0.10   0.31  0.05     NA
#> 61   G2I11     2   3PLM   0.96  0.13      0.89      0.15   0.12  0.05     NA
#> 62   G2I12     2   3PLM   0.00  0.11  21352.80  99999.00   0.21  0.09     NA
#> 63   G2I13     2   3PLM   1.18  0.13     -1.27      0.24   0.21  0.09     NA
#> 64   G2I14     2   3PLM   1.41  0.14      0.23      0.11   0.12  0.05     NA
#> 65   G2I15     2   3PLM   1.58  0.21      0.09      0.17   0.28  0.08     NA
#> 66   G2I16     2   3PLM   1.79  0.25      0.78      0.09   0.22  0.05     NA
#> 67   G2I17     2   3PLM   2.07  0.22     -1.20      0.15   0.20  0.09     NA
#> 68   G2I18     2   3PLM   1.76  0.21      0.66      0.08   0.14  0.04     NA
#> 69   G2I19     2   3PLM   1.11  0.14      0.05      0.22   0.20  0.08     NA
#> 70   G2I20     2   3PLM   1.75  0.48      2.17      0.18   0.16  0.02     NA
#> 71   G2I21     2   3PLM   2.52  0.38      1.56      0.06   0.12  0.02     NA
#> 72   G2I22     2   3PLM   1.52  0.20      0.29      0.15   0.23  0.07     NA
#> 73   G2I23     2   3PLM   2.20  0.29      0.49      0.09   0.29  0.05     NA
#> 74   G2I24     2   3PLM   1.42  0.12      0.38      0.07   0.07  0.03     NA
#> 75   G2I25     2   3PLM   1.80  0.23      1.40      0.06   0.07  0.02     NA
#> 76   G2I26     2   3PLM   1.91  0.20     -0.84      0.17   0.25  0.10     NA
#> 77   G2I27     2   3PLM   1.00  0.12     -0.53      0.26   0.21  0.09     NA
#> 78   G2I28     5    GRM   1.14  0.07     -0.35      0.07   0.17  0.05   0.81
#> 79    C2I1     2   3PLM   1.08  0.08     -0.25      0.11   0.13  0.03     NA
#> 80    C2I2     2   3PLM   0.99  0.09      1.27      0.07   0.06  0.02     NA
#> 81    C2I3     2   3PLM   1.61  0.16      1.41      0.05   0.13  0.01     NA
#> 82    C2I4     2   3PLM   1.55  0.11      0.28      0.06   0.19  0.02     NA
#> 83    C2I5     2   3PLM   1.42  0.10     -0.03      0.07   0.12  0.03     NA
#> 84    C2I6     2   3PLM   2.06  0.11     -0.02      0.04   0.05  0.01     NA
#> 85    C2I7     2   3PLM   1.56  0.11      0.04      0.07   0.18  0.02     NA
#> 86    C2I8     2   3PLM   1.51  0.15      1.30      0.06   0.19  0.02     NA
#> 87    C2I9     2   3PLM   2.10  0.13     -1.03      0.07   0.13  0.03     NA
#> 88   C2I10     2   3PLM   0.00  0.04  45122.03  99999.00   0.21  0.09     NA
#> 89    G3I1     2   3PLM   1.50  0.14     -1.07      0.13   0.17  0.05     NA
#> 90    G3I2     2   3PLM   0.00  0.05    299.41  99999.00   0.24  0.10     NA
#> 91    G3I3     2   3PLM   1.10  0.12      0.07      0.12   0.16  0.03     NA
#> 92    G3I4     2   3PLM   1.32  0.24      1.75      0.12   0.23  0.02     NA
#> 93    G3I5     2   3PLM   0.00  0.05   4605.96  99999.00   0.21  0.09     NA
#> 94    G3I6     2   3PLM   1.17  0.15     -1.27      0.25   0.33  0.07     NA
#> 95    G3I7     2   3PLM   1.36  0.13      0.29      0.09   0.14  0.03     NA
#> 96    G3I8     2   3PLM   1.85  0.17     -0.07      0.07   0.14  0.02     NA
#> 97    G3I9     2   3PLM   1.12  0.12     -1.18      0.20   0.22  0.06     NA
#> 98   G3I10     2   3PLM   1.54  0.22      1.02      0.09   0.29  0.02     NA
#> 99   G3I11     2   3PLM   0.89  0.10      0.84      0.11   0.09  0.03     NA
#> 100  G3I12     2   3PLM   1.42  0.17     -0.73      0.15   0.27  0.05     NA
#> 101  G3I13     2   3PLM   1.13  0.13     -1.01      0.21   0.26  0.06     NA
#> 102  G3I14     2   3PLM   1.34  0.14      0.34      0.09   0.17  0.03     NA
#> 103  G3I15     2   3PLM   1.30  0.14      0.01      0.11   0.20  0.03     NA
#> 104  G3I16     2   3PLM   1.59  0.18      0.79      0.08   0.21  0.02     NA
#> 105  G3I17     2   3PLM   1.55  0.15     -1.59      0.15   0.17  0.06     NA
#> 106  G3I18     2   3PLM   1.36  0.14      0.74      0.08   0.12  0.02     NA
#> 107  G3I19     2   3PLM   0.98  0.11      0.13      0.15   0.17  0.04     NA
#> 108  G3I20     2   3PLM   1.15  0.23      2.32      0.20   0.12  0.02     NA
#> 109  G3I21     2   3PLM   3.02  0.56      1.68      0.06   0.14  0.01     NA
#> 110  G3I22     2   3PLM   1.13  0.11      0.09      0.10   0.09  0.03     NA
#> 111  G3I23     2   3PLM   1.77  0.17      0.35      0.07   0.15  0.02     NA
#> 112  G3I24     2   3PLM   1.12  0.09      0.35      0.08   0.06  0.02     NA
#> 113  G3I25     2   3PLM   1.36  0.15      1.40      0.08   0.05  0.01     NA
#> 114  G3I26     2   3PLM   1.61  0.16     -1.01      0.13   0.21  0.05     NA
#> 115  G3I27     2   3PLM   0.91  0.10     -0.62      0.21   0.18  0.06     NA
#> 116  G3I28     5    GRM   0.93  0.05     -0.36      0.06   0.18  0.06   0.86
#>      se.4  par.5  se.5
#> 1      NA     NA    NA
#> 2      NA     NA    NA
#> 3      NA     NA    NA
#> 4      NA     NA    NA
#> 5      NA     NA    NA
#> 6      NA     NA    NA
#> 7      NA     NA    NA
#> 8      NA     NA    NA
#> 9      NA     NA    NA
#> 10     NA     NA    NA
#> 11     NA     NA    NA
#> 12     NA     NA    NA
#> 13     NA     NA    NA
#> 14     NA     NA    NA
#> 15     NA     NA    NA
#> 16     NA     NA    NA
#> 17     NA     NA    NA
#> 18     NA     NA    NA
#> 19     NA     NA    NA
#> 20     NA     NA    NA
#> 21     NA     NA    NA
#> 22     NA     NA    NA
#> 23     NA     NA    NA
#> 24     NA     NA    NA
#> 25     NA     NA    NA
#> 26     NA     NA    NA
#> 27     NA     NA    NA
#> 28     NA     NA    NA
#> 29     NA     NA    NA
#> 30     NA     NA    NA
#> 31     NA     NA    NA
#> 32     NA     NA    NA
#> 33     NA     NA    NA
#> 34     NA     NA    NA
#> 35     NA     NA    NA
#> 36     NA     NA    NA
#> 37     NA     NA    NA
#> 38     NA     NA    NA
#> 39     NA     NA    NA
#> 40     NA     NA    NA
#> 41     NA     NA    NA
#> 42     NA     NA    NA
#> 43     NA     NA    NA
#> 44     NA     NA    NA
#> 45     NA     NA    NA
#> 46     NA     NA    NA
#> 47     NA     NA    NA
#> 48     NA   1.38    NA
#> 49   0.04  -0.09  0.03
#> 50   0.04   1.27  0.06
#> 51     NA     NA    NA
#> 52     NA     NA    NA
#> 53     NA     NA    NA
#> 54     NA     NA    NA
#> 55     NA     NA    NA
#> 56     NA     NA    NA
#> 57     NA     NA    NA
#> 58     NA     NA    NA
#> 59     NA     NA    NA
#> 60     NA     NA    NA
#> 61     NA     NA    NA
#> 62     NA     NA    NA
#> 63     NA     NA    NA
#> 64     NA     NA    NA
#> 65     NA     NA    NA
#> 66     NA     NA    NA
#> 67     NA     NA    NA
#> 68     NA     NA    NA
#> 69     NA     NA    NA
#> 70     NA     NA    NA
#> 71     NA     NA    NA
#> 72     NA     NA    NA
#> 73     NA     NA    NA
#> 74     NA     NA    NA
#> 75     NA     NA    NA
#> 76     NA     NA    NA
#> 77     NA     NA    NA
#> 78   0.05   1.47  0.07
#> 79     NA     NA    NA
#> 80     NA     NA    NA
#> 81     NA     NA    NA
#> 82     NA     NA    NA
#> 83     NA     NA    NA
#> 84     NA     NA    NA
#> 85     NA     NA    NA
#> 86     NA     NA    NA
#> 87     NA     NA    NA
#> 88     NA     NA    NA
#> 89     NA     NA    NA
#> 90     NA     NA    NA
#> 91     NA     NA    NA
#> 92     NA     NA    NA
#> 93     NA     NA    NA
#> 94     NA     NA    NA
#> 95     NA     NA    NA
#> 96     NA     NA    NA
#> 97     NA     NA    NA
#> 98     NA     NA    NA
#> 99     NA     NA    NA
#> 100    NA     NA    NA
#> 101    NA     NA    NA
#> 102    NA     NA    NA
#> 103    NA     NA    NA
#> 104    NA     NA    NA
#> 105    NA     NA    NA
#> 106    NA     NA    NA
#> 107    NA     NA    NA
#> 108    NA     NA    NA
#> 109    NA     NA    NA
#> 110    NA     NA    NA
#> 111    NA     NA    NA
#> 112    NA     NA    NA
#> 113    NA     NA    NA
#> 114    NA     NA    NA
#> 115    NA     NA    NA
#> 116  0.07   1.61  0.10
#>  Group Parameters: 
#>                      mu  sigma2  sigma
#> estimate(Group1)   0.01    0.96   0.98
#> se(Group1)         0.02    0.03   0.02
#> estimate(Group2)   0.51    0.56   0.75
#> se(Group2)         0.02    0.02   0.01
#> estimate(Group3)  -0.36    2.11   1.45
#> se(Group3)         0.03    0.07   0.02
#> 

# Extract the group parameter estimates (i.e., scale parameters)
getirt(fit.6, what = "group.par")
#> $Group1
#>                    mu     sigma2     sigma
#> estimates 0.008398825 0.95861133 0.9790870
#> se        0.021893051 0.03032153 0.0154846
#> 
#> $Group2
#>                   mu     sigma2      sigma
#> estimates 0.51380432 0.55799843 0.74699293
#> se        0.01670327 0.01764987 0.01181395
#> 
#> $Group3
#>                    mu     sigma2      sigma
#> estimates -0.36201704 2.10782685 1.45183568
#> se         0.03246403 0.06667201 0.02296128
#> 

# }

```
