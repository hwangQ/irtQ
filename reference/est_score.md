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

  Additional arguments passed to
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
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
  of the ability scale. This is used for the following scoring methods:
  `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. Default is `c(-5, 5)`.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. These values are used to
  generate the Gaussian quadrature points and weights. Ignored if
  `method` is `"ML"`, `"MLF"`, `"WL"`, or `"INV.TCC"`. Default is
  `c(0, 1)`.

- nquad:

  An integer indicating the number of Gaussian quadrature points to be
  generated from the normal prior distribution. Used only when `method`
  is `"EAP"` or `"EAP.SUM"`. Ignored for `"ML"`, `"MLF"`, `"WL"`,
  `"MAP"`, and `"INV.TCC"`. Default is `41`.

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

  A numeric value specifying the item slope parameter (i.e.,
  *a*-parameter) for the two imaginary items used in MLF. See
  **Details** below. Default is 3.0.

- fence.b:

  A numeric vector of length two specifying the lower and upper bounds
  of the item difficulty parameters (i.e., *b*-parameters) for the two
  imaginary items in MLF. If `fence.b = NULL`, the values specified in
  the `range` argument are used instead. Default is NULL.

- tol:

  A numeric value specifying the convergence tolerance for the ML, MLF,
  WL, MAP, and inverse TCC scoring methods. Newton-Raphson optimization
  is used for ML, MLF, WL, and MAP, while the bisection method is used
  for inverse TCC. Default is 1e-4.

- max.iter:

  A positive integer specifying the maximum number of iterations allowed
  for the Newton-Raphson optimization. Default is 100.

- se:

  Logical. If `TRUE`, standard errors of ability estimates are computed.
  If `method` is "EAP.SUM" or "INV.TCC", standard errors are always
  returned regardless of this setting. Default is `TRUE`.

- stval.opt:

  A positive integer specifying the starting value option for the ML,
  MLF, WL, and MAP scoring methods. Available options are:

  - 1: Brute-force search (default)

  - 2: Based on observed sum scores

  - 3: Fixed at 0

  See **Details** below for more information.

- intpol:

  Logical. If `TRUE` and `method = "INV.TCC"`, linear interpolation is
  applied to approximate ability estimates for sum scores that cannot be
  directly mapped using the TCC (e.g., when the observed sum score is
  less than the total of item guessing parameters). Default is `TRUE`.
  See **Details** below.

- range.tcc:

  A numeric vector of length two specifying the lower and upper bounds
  of ability estimates when `method = "INV.TCC"`. Default is `c(-7, 7)`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.
  See **Details** below.

- ncore:

  An integer specifying the number of logical CPU cores to use for
  parallel processing. Default is 1. See **Details** below.

## Value

When `method` is one of `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, or `"EAP"`, a
two-column data frame is returned:

- Column 1: Ability estimates

- Column 2: Standard errors of the ability estimates

When `method` is either `"EAP.SUM"` or `"INV.TCC"`, a list with two
components is returned:

- Object 1: A three-column data frame including:

  - Column 1: Observed sum scores

  - Column 2: Ability estimates

  - Column 3: Standard errors of the ability estimates

- Object 2: A score table showing possible raw sum scores and the
  corresponding ability and standard error estimates

## Details

For the MAP scoring method, only a normal prior distribution is
supported for the population distribution.

When there are missing responses in the data set, the missing value must
be explicitly specified using the `missing` argument. Missing data are
properly handled when using the ML, MLF, WL, MAP, or EAP methods.
However, when using the "EAP.SUM" or "INV.TCC" methods, any missing
responses are automatically treated as incorrect (i.e., recoded as 0s).

In the maximum likelihood estimation with fences (MLF; Han, 2016), two
imaginary items based on the 2PL model are introduced. The first
imaginary item functions as the lower fence, and its difficulty
parameter (*b*) should be smaller than any of the difficulty parameters
in the test form. Similarly, the second imaginary item serves as the
upper fence, and its *b* parameter should be greater than any difficulty
value in the test form. Both imaginary items should also have very steep
slopes (i.e., high *a*-parameter values). See Han (2016) for more
details. If `fence.b = NULL`, the function will automatically assign the
lower and upper fences based on the values provided in the `range`
argument.

When the "INV.TCC" method is used with the 3PL model, ability estimates
cannot be obtained for observed sum scores that are less than the sum of
the items' guessing parameters. In such cases, linear interpolation can
be applied by setting `intpol = TRUE`.

Let \\\theta\_{min}\\ and \\\theta\_{max}\\ denote the minimum and
maximum ability estimates, respectively, and let \\\theta\_{X}\\ be the
ability estimate corresponding to the smallest observed sum score, X,
that is greater than or equal to the sum of the guessing parameters.When
linear interpolation is applied, the first value in the `range.tcc`
argument is treated as \\\theta\_{min}\\. A line is then constructed
between the points \\(x = \theta\_{min}, y = 0)\\ and \\(x =
\theta\_{X}, y = X)\\. The second value in `range.tcc` is interpreted as
\\\theta\_{max}\\, which corresponds to the ability estimate for the
maximum observed sum score.

For the "INV.TCC" method, standard errors of ability estimates are
computed using the approach proposed by Lim et al. (2021). The
implementation of inverse TCC scoring in this function is based on a
modified version of the `SNSequate::irt.eq.tse()` function from the
SNSequate package (González, 2014).

For the ML, MLF, WL, and MAP scoring methods, different strategies can
be used to determine the starting value for ability estimation based on
the `stval.opt` argument:

- When `stval.opt = 1` (default), a brute-force search is performed by
  evaluating the log-likelihood at discrete theta values within the
  range specified by `range`, using 0.1 increments. The theta value
  yielding the highest log-likelihood is chosen as the starting value.

- When `stval.opt = 2`, the starting value is derived from the observed
  sum score using a logistic transformation. For example, if the maximum
  possible score (`max.score`) is 30 and the examinee’s observed sum
  score (`obs.score`) is 20, the starting value is
  `log(obs.score / (max.score - obs.score))`.

  - If all responses are incorrect (i.e., `obs.score = 0`), the starting
    value is `log(1 / max.score)`.

  - If all responses are correct (`obs.score = max.score`), the starting
    value is `log(max.score / 1)`.

- When `stval.opt = 3`, the starting value is fixed at 0.

To accelerate ability estimation using the ML, MLF, WL, MAP, and EAP
methods, this function supports parallel processing across multiple
logical CPU cores. The number of cores can be specified via the `ncore`
argument (default is 1).

Note that the standard errors of ability estimates are computed based on
the Fisher expected information for the ML, MLF, WL, and MAP methods.

For the implementation of the WL method, the function references the
`catR::Pi()`, `catR::Ji()`, and `catR::Ii()` functions from the catR
package (Magis & Barrada, 2017).

## Methods (by class)

- `est_score(default)`: Default method to estimate examinees' latent
  ability parameters using a data frame `x` containing the item
  metadata.

- `est_score(est_irt)`: An object created by the function
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

## References

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
in a microcomputer environment. *Psychometrika, 35*, 179-198.

González, J. (2014). SNSequate: Standard and nonstandard statistical
models and methods for test equating. *Journal of Statistical Software,
59*, 1-30.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).*Fundamentals
of item response theory*. Newbury Park, CA: Sage.

Han, K. T. (2016). Maximum likelihood score estimation method with
fences for short-length tests and computerized adaptive tests. *Applied
psychological measurement, 40*(4), 289-301.

Howard, J. P. (2017). *Computational methods for numerical analysis with
R*. New York: Chapman and Hall/CRC.

Kolen, M. J. & Brennan, R. L. (2004). *Test Equating, Scaling, and
Linking* (2nd ed.). New York: Springer

Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT
proficiency estimates. *Educational Measurement: Issues and Practice,
29*(3), 8-14.

Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
approach to evaluate the performance of MST. *Journal of Educational
Measurement, 58*(2), 154-178.

Magis, D., & Barrada, J. R. (2017). Computerized adaptive testing with
R: Recent updates of the package catR. *Journal of Statistical Software,
76*, 1-19.

Stocking, M. L. (1996). An alternative method for scoring adaptive
tests. *Journal of Educational and Behavioral Statistics, 21*(4),
365-389.

Thissen, D. & Orlando, M. (2001). Item response theory for items scored
in two categories. In D. Thissen & H. Wainer (Eds.), *Test scoring*
(pp.73-140). Mahwah, NJ: Lawrence Erlbaum.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item
Response Theory for Scores on Tests Including Polytomous Items with
Ordered Responses. *Applied Psychological Measurement, 19*(1), 39-49.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.

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
  nquad = 30, se = TRUE
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
#> 3         34 -0.41587865 0.2986419
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
#> 35        34 -0.41587865 0.2986419
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
