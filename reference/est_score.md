# Estimate examinees' ability (proficiency) parameters

This function estimates examinees' latent ability parameters. Available
scoring methods include maximum likelihood estimation (ML), maximum
likelihood estimation with fences (MLF; Han, 2016), weighted likelihood
estimation (WL; Warm, 1989), maximum a posteriori estimation (MAP;
Hambleton et al., 1991), expected a posteriori estimation (EAP; Bock &
Mislevy, 1982), EAP summed scoring (Thissen et al., 1995; Thissen &
Orlando, 2001), and inverse test characteristic curve (TCC) scoring
(e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996).

## Usage

``` r
est_score(x, ...)

# Default S3 method
est_score(
  x,
  data,
  D = 1,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  fence.a = 3,
  fence.b = NULL,
  tol = 1e-04,
  max.iter = 100,
  se = TRUE,
  stval.opt = 1,
  intpol = TRUE,
  range.tcc = c(-7, 7),
  missing = NA,
  ncore = 1,
  ...
)

# S3 method for class 'est_irt'
est_score(
  x,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  fence.a = 3,
  fence.b = NULL,
  tol = 1e-04,
  max.iter = 100,
  se = TRUE,
  stval.opt = 1,
  intpol = TRUE,
  range.tcc = c(-7, 7),
  missing = NA,
  ncore = 1,
  ...
)

# S3 method for class 'est_item'
est_score(
  x,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  fence.a = 3,
  fence.b = NULL,
  tol = 1e-04,
  max.iter = 100,
  se = TRUE,
  stval.opt = 1,
  intpol = TRUE,
  range.tcc = c(-7, 7),
  missing = NA,
  ncore = 1,
  ...
)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.), or an object of class `est_irt`
  or `est_item` obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

  See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
  or [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- ...:

  Additional arguments passed to
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  when `ncore > 1`.

- data:

  A matrix of examinees' item responses, with rows for examinees and
  columns for items in the same order as the items in `x`. Responses to
  an item with K score categories must be integers from 0 to K - 1.
  Missing responses are coded as `NA` or, in functions that have a
  `missing` argument, as the value given in `missing`.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- method:

  A character string indicating the scoring method to use. Available
  options are:

  - `"ML"`: Maximum likelihood estimation

  - `"MLF"`: Maximum likelihood estimation with fences (Han, 2016)

  - `"WL"`: Weighted likelihood estimation (Warm, 1989)

  - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)

  - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)

  - `"EAP.SUM"`: Expected a posteriori summed scoring (Thissen et al.,
    1995; Thissen & Orlando, 2001)

  - `"INV.TCC"`: Inverse test characteristic curve scoring (e.g., Kolen
    & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996)

  Default is `"ML"`.

- range:

  A numeric vector of length two specifying the lower and upper bounds
  of the ability scale for the `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`
  methods. Estimates outside these bounds are set to the nearest bound.
  The bounds also define the search grid for `stval.opt = 1` and the
  default fence locations for `"MLF"`. Default is `c(-5, 5)`.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. For `"MAP"`, it defines
  the prior. For `"EAP"` and `"EAP.SUM"`, it is used to generate the
  Gaussian quadrature points and weights when `weights = NULL`. Ignored
  for `"ML"`, `"MLF"`, `"WL"`, and `"INV.TCC"`. Default is `c(0, 1)`.

- nquad:

  An integer indicating the number of Gaussian quadrature points to be
  generated from the normal prior distribution. Used only when `method`
  is `"EAP"` or `"EAP.SUM"` and `weights = NULL`. Default is `41`.

- weights:

  A two-column matrix or data frame containing the quadrature points (in
  the first column) and their corresponding weights (in the second
  column) for the latent variable prior distribution. The weights and
  points can be conveniently generated using the function
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).

  If `NULL` and `method` is either `"EAP"` or `"EAP.SUM"`, default
  quadrature values are generated based on the `norm.prior` and `nquad`
  arguments. Ignored if `method` is `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, or
  `"INV.TCC"`.

- fence.a:

  A numeric value specifying the slope (*a*) parameter of the two fence
  items used in MLF. See **Details** below. Default is 3.0.

- fence.b:

  A numeric vector of length two specifying the difficulty (*b*)
  parameters of the lower and upper fence items used in MLF. If `NULL`,
  the values in `range` are used. Default is `NULL`.

- tol:

  A positive numeric value specifying the convergence tolerance. For
  `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`, the Fisher scoring iterations
  stop when the absolute update is less than `tol`, and `tol` is also
  the tolerance of the root search used when the iterations do not
  converge. For `"INV.TCC"`, the bisection search stops when the search
  interval is no wider than `tol`. Default is 1e-4.

- max.iter:

  A positive integer specifying the maximum number of Fisher scoring
  iterations for `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. Default is 100.
  See **Details** for the case where the iterations do not converge.

- se:

  Logical. If `TRUE`, standard errors of ability estimates are computed.
  If `FALSE`, the standard errors are returned as `NA`. For `"EAP.SUM"`
  and `"INV.TCC"`, standard errors are always computed. Default is
  `TRUE`.

- stval.opt:

  A positive integer specifying the starting value option for the ML,
  MLF, WL, and MAP scoring methods. Available options are:

  - 1: Grid search (default)

  - 2: Based on observed sum scores

  - 3: Fixed at 0

  See **Details** below for more information.

- intpol:

  Logical. If `TRUE` and `method = "INV.TCC"`, ability estimates are
  assigned to sum scores that cannot be mapped through the TCC (sum
  scores less than or equal to the sum of the guessing parameters, and
  the maximum possible sum score). Default is `TRUE`. See **Details**
  below.

- range.tcc:

  A numeric vector of length two giving the ability estimates assigned
  to the lowest and the maximum possible sum scores when
  `method = "INV.TCC"` and `intpol = TRUE`. Default is `c(-7, 7)`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.
  See **Details** below.

- ncore:

  An integer specifying the number of logical CPU cores used for
  parallel scoring with `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, and `"EAP"`.
  Default is 1. See **Details** below.

## Value

For `method = "ML"`, `"MLF"`, `"WL"`, `"MAP"`, or `"EAP"`, a data frame
with one row per examinee and two columns:

- `est.theta`: Ability estimates.

- `se.theta`: Standard errors of the ability estimates (`NA` when
  `se = FALSE`).

For `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`, the standard error is set to
99.9999 when the estimate equals a bound of `range`. Examinees with all
responses missing receive `NA` in both columns.

For `method = "EAP.SUM"` or `"INV.TCC"`, a list with two data frames:

- `est.par`: One row per examinee with the columns `sum.score`,
  `est.theta`, and `se.theta`.

- `score.table`: One row per possible sum score, from 0 to the maximum
  possible score, with the same three columns.

For `"INV.TCC"`, sum scores without an estimate (see **Details**) have
`NA` in `est.theta` and `se.theta`.

## Details

For `"MAP"`, the prior is the normal distribution given by `norm.prior`.

Each response must be an integer from 0 to the number of categories of
the item minus 1, and `data` must have one column per item in `x`.
Missing responses must be coded as `NA` or declared with the `missing`
argument. For `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, and `"EAP"`, missing
responses are excluded from the likelihood, and examinees with all
responses missing receive `NA` with a warning. For `"EAP.SUM"` and
`"INV.TCC"`, missing responses are recoded as 0 with a warning.

In maximum likelihood estimation with fences (MLF; Han, 2016), two
imaginary 2PL items are added to every response pattern: a lower fence
item answered correctly and an upper fence item answered incorrectly.
The lower fence should have a *b*-parameter below, and the upper fence a
*b*-parameter above, all item difficulties in the test. Both should have
steep slopes (large *a*-parameters). If `fence.b = NULL`, the fences are
placed at the bounds given in `range`. See Han (2016) for details.

For `"INV.TCC"`, the ability estimate for a sum score X is the root of
TCC(theta) = X, found by the bisection method (Howard, 2017). No root
exists for the maximum possible sum score or for sum scores less than or
equal to the sum of the guessing parameters, which is only the score 0
when there are no guessing parameters. With `intpol = TRUE`, these
scores receive estimates as follows. Let \\\theta\_{min}\\ and
\\\theta\_{max}\\ be the first and second values of `range.tcc`, and let
\\\theta\_{X}\\ be the estimate for the smallest sum score X that is
greater than the sum of the guessing parameters. Scores below X are
mapped onto the line through \\(x = \theta\_{min}, y = 0)\\ and \\(x =
\theta\_{X}, y = X)\\, and the maximum possible sum score receives
\\\theta\_{max}\\. If \\\theta\_{min}\\ is above \\\theta\_{X}\\ or
\\\theta\_{max}\\ is below the largest root, a warning is issued and no
values are assigned. With `intpol = FALSE`, these scores receive `NA`.

For `"INV.TCC"`, the standard error for a sum score is the standard
deviation of the inverse TCC estimates over the conditional sum score
distribution at the estimate for that score, computed with the
Lord-Wingersky recursion (Lim et al., 2021; see
[`lwrc()`](https://hwangQ.github.io/irtQ/reference/lwrc.md)). When some
sum scores have no estimate (`intpol = FALSE`, or the warning described
above), the distribution is restricted to the sum scores with an
estimate and rescaled to sum to one. The implementation is based on a
modified version of `SNSequate::irt.eq.tse()` from the SNSequate package
(Gonzalez, 2014).

For the ML, MLF, WL, and MAP scoring methods, different strategies can
be used to determine the starting value for ability estimation based on
the `stval.opt` argument:

- When `stval.opt = 1` (default), the log-likelihood (log-posterior for
  `"MAP"`) is evaluated on a grid from `range[1]` to `range[2]` in steps
  of 0.1. The grid point at the highest local maximum is the starting
  value. If there is no local maximum on the grid, the starting value is
  0.

- When `stval.opt = 2`, the starting value is the log-odds of the
  observed sum score `obs.score` relative to the maximum possible score
  `max.score`, both computed over the non-missing items:
  `log(obs.score / (max.score - obs.score))`.

  - If `obs.score = 0`, the starting value is `log(1 / max.score)`.

  - If `obs.score = max.score`, the starting value is `log(max.score)`.

- When `stval.opt = 3`, the starting value is fixed at 0.

For the ML, MLF, WL, and MAP scoring methods, the estimate is found by
Fisher scoring. This is a Newton-type iteration that divides the first
derivative of the log-likelihood (of the log-posterior for `"MAP"`, or
Warm's weighted score function for `"WL"`) by the expected Fisher
information instead of the observed second derivative, and a single
update is capped at 1 in absolute value. The iterations stop when the
absolute update is less than `tol`. If they have not converged after
`max.iter` iterations (for example, when the iterations alternate
between two values for a low-scoring examinee on a test with 3PLM
items), the estimate is found as follows. Among the values visited
during the iterations, the function takes the closest pair in which the
score function (the derivative above) is positive at the lower value and
negative at the upper value, so that the pair encloses a local maximum,
and finds the root of the score function between them with
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) using the
tolerance `tol`. If there is no such pair, the last value of the
iterations is used. This search is not carried out when the iterations
converge. In all cases, the estimate is truncated to `range`, and the
standard error is computed at the final estimate.

For `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, and `"EAP"`, examinees can be
scored in parallel by setting `ncore` greater than 1. A warning is
issued when `ncore > 1` and there are fewer than 5,000 examinees,
because the parallel overhead then exceeds the computation time.

For `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`, the standard error is the
inverse square root of the expected Fisher information at the estimate.
The information includes the two fence items for `"MLF"` and the prior
information \\1/\sigma^2\\ for `"MAP"`. For `"EAP"` and `"EAP.SUM"`, the
standard error is the posterior standard deviation. For `"EAP.SUM"`, the
posterior for each sum score is computed with the Lord-Wingersky
recursion (see
[`lwrc()`](https://hwangQ.github.io/irtQ/reference/lwrc.md)).

For the implementation of the WL method, the function references the
`catR::Pi()`, `catR::Ji()`, and `catR::Ii()` functions from the catR
package (Magis & Barrada, 2017).

## Methods (by class)

- `est_score(default)`: Default method to estimate examinees' latent
  ability parameters using a data frame `x` containing the item
  metadata.

- `est_score(est_irt)`: Method for an object of class `est_irt`. The
  item parameter estimates, response data, and scaling constant `D` are
  taken from `x`, so supplying `data` or `D` stops with an error.

- `est_score(est_item)`: Method for an object of class `est_item`. The
  item parameter estimates, response data, and scaling constant `D` are
  taken from `x`, so supplying `data` or `D` stops with an error.

## References

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
in a microcomputer environment. *Applied Psychological Measurement,
6*(4), 431-444.
[doi:10.1177/014662168200600405](https://doi.org/10.1177/014662168200600405)
.

Gonzalez, J. (2014). SNSequate: Standard and nonstandard statistical
models and methods for test equating. *Journal of Statistical Software,
59*(7), 1-30.
[doi:10.18637/jss.v059.i07](https://doi.org/10.18637/jss.v059.i07) .

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals
of item response theory*. Newbury Park, CA: Sage.

Han, K. T. (2016). Maximum likelihood score estimation method with
fences for short-length tests and computerized adaptive tests. *Applied
Psychological Measurement, 40*(4), 289-301.
[doi:10.1177/0146621616631317](https://doi.org/10.1177/0146621616631317)
.

Howard, J. P. (2017). *Computational methods for numerical analysis with
R*. New York: Chapman and Hall/CRC.

Kolen, M. J., & Brennan, R. L. (2004). *Test equating, scaling, and
linking* (2nd ed.). Springer.

Kolen, M. J., & Tong, Y. (2010). Psychometric properties of IRT
proficiency estimates. *Educational Measurement: Issues and Practice,
29*(3), 8-14.
[doi:10.1111/j.1745-3992.2010.00179.x](https://doi.org/10.1111/j.1745-3992.2010.00179.x)
.

Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
approach to evaluate the performance of MST. *Journal of Educational
Measurement, 58*(2), 154-178.
[doi:10.1111/jedm.12276](https://doi.org/10.1111/jedm.12276) .

Magis, D., & Barrada, J. R. (2017). Computerized adaptive testing with
R: Recent updates of the package catR. *Journal of Statistical Software,
76*(Code Snippet 1), 1-19.
[doi:10.18637/jss.v076.c01](https://doi.org/10.18637/jss.v076.c01) .

Stocking, M. L. (1996). An alternative method for scoring adaptive
tests. *Journal of Educational and Behavioral Statistics, 21*(4),
365-389.
[doi:10.3102/10769986021004365](https://doi.org/10.3102/10769986021004365)
.

Thissen, D., & Orlando, M. (2001). Item response theory for items scored
in two categories. In D. Thissen & H. Wainer (Eds.), *Test scoring* (pp.
73-140). Mahwah, NJ: Lawrence Erlbaum.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995).
Item response theory for scores on tests including polytomous items with
ordered responses. *Applied Psychological Measurement, 19*(1), 39-49.
[doi:10.1177/014662169501900105](https://doi.org/10.1177/014662169501900105)
.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627) .

## See also

[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Read item parameters and convert them into item metadata
x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df

# Generate examinee ability values
set.seed(12)
theta <- rnorm(10)

# Simulate item response data based on the item metadata and abilities
data <- simdat(x, theta, D = 1)

# \donttest{
# Estimate abilities using maximum likelihood (ML)
est_score(x, data, D = 1, method = "ML", range = c(-4, 4), se = TRUE)
#>      est.theta  se.theta
#> 1  -2.37597568 0.5432999
#> 2   1.54041176 0.3411440
#> 3  -0.39456473 0.2715142
#> 4  -0.80001934 0.2861069
#> 5  -2.61351340 0.6288578
#> 6   0.03097513 0.2686083
#> 7   0.07766689 0.2691152
#> 8  -0.81133014 0.2866940
#> 9   0.38832795 0.2764596
#> 10  0.42226686 0.2776366

# Estimate abilities using weighted likelihood (WL)
est_score(x, data, D = 1, method = "WL", range = c(-4, 4), se = TRUE)
#>      est.theta  se.theta
#> 1  -2.25452427 0.5062071
#> 2   1.50571058 0.3382932
#> 3  -0.39564529 0.2715364
#> 4  -0.79478333 0.2858408
#> 5  -2.42123985 0.5582300
#> 6   0.02247278 0.2685341
#> 7   0.06952197 0.2690142
#> 8  -0.80657712 0.2864409
#> 9   0.37324895 0.2759575
#> 10  0.40825670 0.2771420

# Estimate abilities using MLF with default fences
# based on the `range` argument
est_score(x, data,
  D = 1, method = "MLF",
  fence.a = 3.0, fence.b = NULL, se = TRUE
)
#>      est.theta  se.theta
#> 1  -2.37558947 0.5429012
#> 2   1.54040048 0.3411375
#> 3  -0.39456453 0.2715141
#> 4  -0.80001852 0.2861065
#> 5  -2.61218863 0.6274710
#> 6   0.03097511 0.2686082
#> 7   0.07766686 0.2691152
#> 8  -0.81132937 0.2866936
#> 9   0.38832773 0.2764595
#> 10  0.42226664 0.2776365

# Estimate abilities using MLF with user-specified fences
est_score(x, data,
  D = 1, method = "MLF", fence.a = 3.0,
  fence.b = c(-7, 7), se = TRUE
)
#>      est.theta  se.theta
#> 1  -2.37597472 0.5432989
#> 2   1.54041173 0.3411440
#> 3  -0.39456473 0.2715142
#> 4  -0.80001934 0.2861069
#> 5  -2.61351009 0.6288543
#> 6   0.03097513 0.2686083
#> 7   0.07766689 0.2691152
#> 8  -0.81133014 0.2866940
#> 9   0.38832795 0.2764596
#> 10  0.42226686 0.2776366

# Estimate abilities using maximum a posteriori (MAP)
est_score(x, data,
  D = 1, method = "MAP", norm.prior = c(0, 1),
  se = TRUE
)
#>      est.theta  se.theta
#> 1  -1.90935697 0.3880537
#> 2   1.37883556 0.3122734
#> 3  -0.36718683 0.2615070
#> 4  -0.73763142 0.2723683
#> 5  -2.02185605 0.4072138
#> 6   0.02869450 0.2593945
#> 7   0.07265826 0.2598135
#> 8  -0.75565388 0.2731196
#> 9   0.35887600 0.2655943
#> 10  0.39324460 0.2666123

# Estimate abilities using expected a posteriori (EAP)
est_score(x, data,
  D = 1, method = "EAP", norm.prior = c(0, 1),
  nquad = 30, se = TRUE
)
#>      est.theta  se.theta
#> 1  -1.98741065 0.4164359
#> 2   1.39979503 0.3163898
#> 3  -0.36445308 0.2470032
#> 4  -0.76884012 0.2757507
#> 5  -2.12589479 0.4317501
#> 6   0.03088081 0.2998967
#> 7   0.09087299 0.2819628
#> 8  -0.79743301 0.2484775
#> 9   0.35318919 0.2606670
#> 10  0.37559150 0.2512978

# Estimate abilities using EAP summed scoring
est_score(x, data,
  D = 1, method = "EAP.SUM", norm.prior = c(0, 1),
  nquad = 30
)
#> $est.par
#>    sum.score   est.theta  se.theta
#> 1         15 -1.84945881 0.4193692
#> 2         55  1.01201284 0.3327285
#> 3         34 -0.41587865 0.2986418
#> 4         31 -0.62087399 0.3218871
#> 5         13 -2.04576501 0.4430296
#> 6         40 -0.04315730 0.3144901
#> 7         41  0.02805944 0.3150827
#> 8         29 -0.76191091 0.3134347
#> 9         46  0.35095185 0.2972605
#> 10        50  0.63635923 0.3261820
#> 
#> $score.table
#>    sum.score   est.theta  se.theta
#> 1          0 -3.12665588 0.4947348
#> 2          1 -3.06696536 0.4959440
#> 3          2 -3.00363668 0.4968896
#> 4          3 -2.93631019 0.4975675
#> 5          4 -2.86462002 0.4978945
#> 6          5 -2.78826682 0.4976399
#> 7          6 -2.70712182 0.4963949
#> 8          7 -2.62133160 0.4936339
#> 9          8 -2.53136890 0.4888928
#> 10         9 -2.43797666 0.4820164
#> 11        10 -2.34200326 0.4733294
#> 12        11 -2.24421762 0.4635516
#> 13        12 -2.14525247 0.4533911
#> 14        13 -2.04576501 0.4430296
#> 15        14 -1.94671250 0.4319439
#> 16        15 -1.84945881 0.4193692
#> 17        16 -1.75547421 0.4051694
#> 18        17 -1.66571260 0.3904342
#> 19        18 -1.58007691 0.3772410
#> 20        19 -1.49739146 0.3675672
#> 21        20 -1.41598904 0.3619472
#> 22        21 -1.33464611 0.3588298
#> 23        22 -1.25335609 0.3552883
#> 24        23 -1.17343058 0.3487230
#> 25        24 -1.09676958 0.3384327
#> 26        25 -1.02470912 0.3261176
#> 27        26 -0.95715839 0.3151243
#> 28        27 -0.89249371 0.3088107
#> 29        28 -0.82819457 0.3087049
#> 30        29 -0.76191091 0.3134347
#> 31        30 -0.69254754 0.3192458
#> 32        31 -0.62087399 0.3218871
#> 33        32 -0.54922749 0.3187256
#> 34        33 -0.48033020 0.3099879
#> 35        34 -0.41587865 0.2986418
#> 36        35 -0.35572290 0.2891027
#> 37        36 -0.29796720 0.2852870
#> 38        37 -0.23977641 0.2887030
#> 39        38 -0.17850221 0.2975782
#> 40        39 -0.11277486 0.3077038
#> 41        40 -0.04315730 0.3144901
#> 42        41  0.02805944 0.3150827
#> 43        42  0.09798419 0.3096104
#> 44        43  0.16459743 0.3010856
#> 45        44  0.22770650 0.2940215
#> 46        45  0.28894048 0.2923086
#> 47        46  0.35095185 0.2972605
#> 48        47  0.41626072 0.3069569
#> 49        48  0.48617270 0.3173587
#> 50        49  0.56020071 0.3244417
#> 51        50  0.63635923 0.3261820
#> 52        51  0.71230767 0.3235277
#> 53        52  0.78672325 0.3198855
#> 54        53  0.86006860 0.3192397
#> 55        54  0.93434155 0.3238112
#> 56        55  1.01201284 0.3327285
#> 57        56  1.09473828 0.3427152
#> 58        57  1.18255707 0.3503152
#> 59        58  1.27411421 0.3541637
#> 60        59  1.36784810 0.3559114
#> 61        60  1.46342571 0.3591339
#> 62        61  1.56252098 0.3666446
#> 63        62  1.66819684 0.3779267
#> 64        63  1.78281178 0.3892442
#> 65        64  1.90635687 0.3975784
#> 66        65  2.03833819 0.4047557
#> 67        66  2.18202764 0.4156153
#> 68        67  2.34491493 0.4325989
#> 69        68  2.53427044 0.4563203
#> 70        69  2.75296022 0.4866583
#> 71        70  3.00707957 0.5256893
#> 

# Estimate abilities using inverse TCC scoring
est_score(x, data,
  D = 1, method = "INV.TCC", intpol = TRUE,
  range.tcc = c(-7, 7)
)
#> $est.par
#>    sum.score   est.theta  se.theta
#> 1         15 -2.18410156 0.8408969
#> 2         55  1.12566406 0.3712260
#> 3         34 -0.45707031 0.3196167
#> 4         31 -0.66792969 0.3287452
#> 5         13 -2.55285156 1.0468414
#> 6         40 -0.04410156 0.3146460
#> 7         41  0.02519531 0.3152754
#> 8         29 -0.81285156 0.3382760
#> 9         46  0.38222656 0.3243666
#> 10        50  0.69035156 0.3392010
#> 
#> $score.table
#>    sum.score   est.theta  se.theta
#> 1          0 -7.00000000 1.1103831
#> 2          1 -6.74910156 1.1119776
#> 3          2 -6.49820313 1.1139931
#> 4          3 -6.24730469 1.1165384
#> 5          4 -5.99640625 1.1197472
#> 6          5 -5.74550781 1.1237801
#> 7          6 -5.49460938 1.1288239
#> 8          7 -5.24371094 1.1350844
#> 9          8 -4.99281250 1.1427654
#> 10         9 -4.74191406 1.1520241
#> 11        10 -3.62566406 1.2024653
#> 12        11 -3.12660156 1.1913471
#> 13        12 -2.79917969 1.1339512
#> 14        13 -2.55285156 1.0468414
#> 15        14 -2.35332031 0.9450795
#> 16        15 -2.18410156 0.8408969
#> 17        16 -2.03613281 0.7431657
#> 18        17 -1.90363281 0.6571415
#> 19        18 -1.78292969 0.5851965
#> 20        19 -1.67152344 0.5273994
#> 21        20 -1.56753906 0.4822437
#> 22        21 -1.46964844 0.4475187
#> 23        22 -1.37675781 0.4208710
#> 24        23 -1.28800781 0.4002422
#> 25        24 -1.20269531 0.3840154
#> 26        25 -1.12042969 0.3710377
#> 27        26 -1.04066406 0.3604623
#> 28        27 -0.96300781 0.3517186
#> 29        28 -0.88714844 0.3444143
#> 30        29 -0.81285156 0.3382760
#> 31        30 -0.73980469 0.3330989
#> 32        31 -0.66792969 0.3287452
#> 33        32 -0.59691406 0.3250969
#> 34        33 -0.52667969 0.3220755
#> 35        34 -0.45707031 0.3196167
#> 36        35 -0.38785156 0.3176693
#> 37        36 -0.31902344 0.3162006
#> 38        37 -0.25035156 0.3151809
#> 39        38 -0.18175781 0.3145908
#> 40        39 -0.11300781 0.3144157
#> 41        40 -0.04410156 0.3146460
#> 42        41  0.02519531 0.3152754
#> 43        42  0.09503906 0.3163010
#> 44        43  0.16542969 0.3177190
#> 45        44  0.23667969 0.3195335
#> 46        45  0.30894531 0.3217490
#> 47        46  0.38222656 0.3243666
#> 48        47  0.45683594 0.3274011
#> 49        48  0.53292969 0.3308677
#> 50        49  0.61066406 0.3347878
#> 51        50  0.69035156 0.3392010
#> 52        51  0.77207031 0.3441483
#> 53        52  0.85613281 0.3497026
#> 54        53  0.94285156 0.3559647
#> 55        54  1.03261719 0.3630787
#> 56        55  1.12566406 0.3712260
#> 57        56  1.22261719 0.3806898
#> 58        57  1.32394531 0.3918367
#> 59        58  1.43042969 0.4052066
#> 60        59  1.54308594 0.4215904
#> 61        60  1.66324219 0.4421928
#> 62        61  1.79253906 0.4689061
#> 63        62  1.93371094 0.5049851
#> 64        63  2.09019531 0.5558304
#> 65        64  2.26738281 0.6304313
#> 66        65  2.47363281 0.7424616
#> 67        66  2.72269531 0.9093035
#> 68        67  3.04050781 1.1450429
#> 69        68  3.48441406 1.4331872
#> 70        69  4.23816406 1.6207472
#> 71        70  7.00000000 0.7318590
#> 
# }
```
