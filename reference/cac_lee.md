# Classification Accuracy and Consistency Using Lee's (2010) Approach

This function computes classification accuracy and consistency indices
for complex assessments using the item response theory (IRT) approach of
Lee (2010). Tests may contain dichotomous items, polytomous items, or
both.

## Usage

``` r
cac_lee(x, cutscore, theta = NULL, weights = NULL, D = 1, cut.obs = TRUE)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.). See
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- cutscore:

  A numeric vector of finite cut scores in strictly ascending order. The
  K - 1 cut scores divide examinees into K performance levels. With
  `cut.obs = TRUE`, the cut scores are on the observed summed score
  metric. With `cut.obs = FALSE`, they are on the ability (theta)
  metric. A score equal to a cut score is assigned to the higher level.

- theta:

  A numeric vector of ability estimates, one per examinee, used by the P
  method. It is ignored when `weights` is supplied. Missing values
  (`NA`) are excluded with a warning. Default is `NULL`.

- weights:

  A two-column data frame or matrix used by the D method. The first
  column holds the quadrature points (nodes), and the second column
  holds the corresponding weights, which should sum to 1; weights that
  do not sum to 1 are rescaled to sum to 1 with a warning (for example,
  frequencies can be given as weights). The weights must be nonnegative
  and finite and have a positive sum.
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)
  creates such a data frame. When `weights` is supplied, `theta` is
  ignored. Default is `NULL`.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- cut.obs:

  Logical. If `TRUE`, `cutscore` is on the observed summed score metric.
  If `FALSE`, `cutscore` is on the ability (theta) metric and is
  converted to the summed score metric with the test characteristic
  curve (TCC). Default is `TRUE`.

## Value

A list with the following elements:

- `confusion`: A K x K matrix of expected proportions, with the true
  levels in rows (dimension name `True`) and the expected (observed)
  levels in columns (dimension name `Expected`). The entries are rounded
  to 7 decimal places and sum to 1 up to rounding.

- `marginal`: A data frame with the columns `level`, `accuracy`, and
  `consistency`. The rows for levels 1 to K give the contribution of the
  ability values whose true level is that level, that is, the weighted
  sums of their conditional indices. These rows add up to the last row,
  `marginal`, which holds the marginal classification accuracy and
  consistency.

- `conditional`: A data frame with one row per ability value and the
  columns `theta`, `weights`, `true.score` (expected summed score),
  `level` (true level, a factor), `accuracy`, and `consistency`.

- `prob.level`: A data frame with the columns `theta`, `weights`,
  `true.score`, and `level`, followed by `p.level.1` to `p.level.K`, the
  probabilities of being classified into each level.

- `cutscore`: The cut scores used, on the observed summed score metric.
  When `cut.obs = FALSE`, these are the converted cut scores.

## Details

This function first validates the input arguments. If both `theta` and
`weights` are `NULL`, the function will stop and return an error
message. Either `theta` (a vector of ability estimates) or `weights` (a
quadrature-based weight matrix) must be specified.

If `cut.obs = FALSE`, the provided cut scores are assumed to be on the
theta (ability) scale, and they are internally converted to the observed
summed score scale using the test characteristic curve (TCC). This
transformation allows classification to be carried out on the summed
score metric, even if theta-based cut points are provided.

For each ability value (a quadrature point in the D method, or an
ability estimate in the P method), the conditional distribution of the
observed summed score is computed with the Lord-Wingersky recursion
(Lord & Wingersky, 1984) and its extension to polytomous items (see
[`lwrc()`](https://hwangQ.github.io/irtQ/reference/lwrc.md)). Summing
this distribution over the score ranges defined by the cut scores gives
the probability of being classified into each level. The true level of
an ability value is found by comparing its expected summed score (the
TCC value) with the same cut scores, so the observed cut scores also
serve as the true cut scores. The conditional classification accuracy is
the probability of being classified into the true level, and the
conditional classification consistency is the sum of the squared level
probabilities, which is the probability of the same classification on
two independent administrations (Lee, 2010).

In the D method (`weights` supplied), the marginal indices are the
weighted sums of the conditional indices over the quadrature points. In
the P method (`theta` supplied), each examinee receives the weight 1/N,
where N is the number of ability estimates (after the missing values are
excluded), so the marginal indices are the averages of the conditional
indices.

## References

Lee, W.-C. (2010). Classification consistency and accuracy for complex
assessments using item response theory. *Journal of Educational
Measurement, 47*(1), 1-17.
[doi:10.1111/j.1745-3984.2009.00096.x](https://doi.org/10.1111/j.1745-3984.2009.00096.x)
.

Lord, F. M., & Wingersky, M. S. (1984). Comparison of IRT true-score and
equipercentile observed-score "equatings." *Applied Psychological
Measurement, 8*(4), 453-461.
[doi:10.1177/014662168400800409](https://doi.org/10.1177/014662168400800409)
.

## See also

[`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md),
[`cac_rud()`](https://hwangQ.github.io/irtQ/reference/cac_rud.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
## --------------------------------------------------------------------------
## 1. When the ability distribution is given by quadrature points and
##    weights (D method)
## --------------------------------------------------------------------------

# Import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Read item parameter data and convert it to item metadata
x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df

# Set the cut scores on the observed summed score scale
cutscore <- c(10, 20, 30, 50)

# Create a data frame containing the quadrature points and corresponding weights
node <- seq(-4, 4, 0.25)
weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = node)

# Calculate classification accuracy and consistency
cac_1 <- cac_lee(x = x, cutscore = cutscore, weights = weights, D = 1)
print(cac_1)
#> $confusion
#>     Expected
#> True         1         2         3         4         5
#>    1 0.0000591 0.0000624 0.0000002 0.0000000 0.0000000
#>    2 0.0027978 0.0371775 0.0114124 0.0001049 0.0000000
#>    3 0.0000518 0.0247711 0.1301180 0.0588848 0.0000006
#>    4 0.0000000 0.0000592 0.0223073 0.4044774 0.0422750
#>    5 0.0000000 0.0000000 0.0000000 0.0389876 0.2264530
#> 
#> $marginal
#>     level     accuracy  consistency
#>         1 5.905558e-05 6.076492e-05
#>         2 3.717753e-02 3.110889e-02
#>         3 1.301180e-01 1.136090e-01
#>         4 4.044774e-01 3.686834e-01
#>         5 2.264530e-01 2.112560e-01
#>  marginal 7.982849e-01 7.247180e-01
#> 
#> $conditional
#>    theta      weights true.score level  accuracy consistency
#> 1  -4.00 3.345874e-05   9.528758     1 0.5147271   0.4995292
#> 2  -3.75 8.815204e-05   9.822381     1 0.4745603   0.4997196
#> 3  -3.50 2.181784e-04  10.204752     2 0.5732252   0.5085288
#> 4  -3.25 5.072800e-04  10.702236     2 0.6325676   0.5316859
#> 5  -3.00 1.108001e-03  11.347749     2 0.6998812   0.5744329
#> 6  -2.75 2.273471e-03  12.180947     2 0.7681380   0.6354594
#> 7  -2.50 4.382230e-03  13.247004     2 0.8221242   0.6960113
#> 8  -2.25 7.935194e-03  14.593125     2 0.8362343   0.7128993
#> 9  -2.00 1.349822e-02  16.262694     2 0.7778819   0.6430520
#> 10 -1.75 2.157009e-02  18.288173     2 0.6245373   0.5220154
#> 11 -1.50 3.238054e-02  20.684127     3 0.5802337   0.4952651
#> 12 -1.25 4.566389e-02  23.440794     3 0.7381357   0.5838746
#> 13 -1.00 6.049482e-02  26.520161     3 0.7032259   0.5575624
#> 14 -0.75 7.528702e-02  29.859615     3 0.4659769   0.4938488
#> 15 -0.50 8.801945e-02  33.382519     4 0.8041392   0.6846800
#> 16 -0.25 9.667045e-02  37.005127     4 0.9507335   0.9060773
#> 17  0.00 9.973910e-02  40.636849     4 0.9717552   0.9448413
#> 18  0.25 9.667045e-02  44.185216     4 0.8882987   0.8014650
#> 19  0.50 8.801945e-02  47.569664     4 0.6702547   0.5579645
#> 20  0.75 7.528702e-02  50.732806     5 0.6256424   0.5315718
#> 21  1.00 6.049482e-02  53.640666     5 0.8548586   0.7518493
#> 22  1.25 4.566389e-02  56.275042     5 0.9613577   0.9257019
#> 23  1.50 3.238054e-02  58.624498     5 0.9927425   0.9855904
#> 24  1.75 2.157009e-02  60.679812     5 0.9989876   0.9979773
#> 25  2.00 1.349822e-02  62.437625     5 0.9998878   0.9997757
#> 26  2.25 7.935194e-03  63.907982     5 0.9999894   0.9999788
#> 27  2.50 4.382230e-03  65.116198     5 0.9999991   0.9999982
#> 28  2.75 2.273471e-03  66.097088     5 0.9999999   0.9999998
#> 29  3.00 1.108001e-03  66.887482     5 1.0000000   1.0000000
#> 30  3.25 5.072800e-04  67.521427     5 1.0000000   1.0000000
#> 31  3.50 2.181784e-04  68.028320     5 1.0000000   1.0000000
#> 32  3.75 8.815204e-05  68.432707     5 1.0000000   1.0000000
#> 33  4.00 3.345874e-05  68.754747     5 1.0000000   1.0000000
#> 
#> $prob.level
#>    theta      weights true.score level    p.level.1    p.level.2    p.level.3
#> 1  -4.00 3.345874e-05   9.528758     1 5.147271e-01 4.843391e-01 9.337864e-04
#> 2  -3.75 8.815204e-05   9.822381     1 4.745603e-01 5.239369e-01 1.502750e-03
#> 3  -3.50 2.181784e-04  10.204752     2 4.241874e-01 5.732252e-01 2.587311e-03
#> 4  -3.25 5.072800e-04  10.702236     2 3.626588e-01 6.325676e-01 4.773254e-03
#> 5  -3.00 1.108001e-03  11.347749     2 2.907072e-01 6.998812e-01 9.410338e-03
#> 6  -2.75 2.273471e-03  12.180947     2 2.122214e-01 7.681380e-01 1.963481e-02
#> 7  -2.50 4.382230e-03  13.247004     2 1.353361e-01 8.221242e-01 4.251153e-02
#> 8  -2.25 7.935194e-03  14.593125     2 7.117453e-02 8.362343e-01 9.244310e-02
#> 9  -2.00 1.349822e-02  16.262694     2 2.861317e-02 7.778819e-01 1.926976e-01
#> 10 -1.75 2.157009e-02  18.288173     2 8.002132e-03 6.245373e-01 3.631612e-01
#> 11 -1.50 3.238054e-02  20.684127     3 1.395918e-03 3.976994e-01 5.802337e-01
#> 12 -1.25 4.566389e-02  23.440794     3 1.352758e-04 1.797482e-01 7.381357e-01
#> 13 -1.00 6.049482e-02  26.520161     3 6.523936e-06 5.091473e-02 7.032259e-01
#> 14 -0.75 7.528702e-02  29.859615     3 1.429429e-07 8.040194e-03 4.659769e-01
#> 15 -0.50 8.801945e-02  33.382519     4 1.332844e-09 6.449138e-04 1.950377e-01
#> 16 -0.25 9.667045e-02  37.005127     4 5.077770e-12 2.469659e-05 4.665223e-02
#> 17  0.00 9.973910e-02  40.636849     4 7.760281e-15 4.385747e-07 5.929865e-03
#> 18  0.25 9.667045e-02  44.185216     4 4.781963e-18 3.621528e-09 3.895648e-04
#> 19  0.50 8.801945e-02  47.569664     4 1.223175e-21 1.430172e-11 1.343626e-05
#> 20  0.75 7.528702e-02  50.732806     5 1.363265e-25 2.818451e-14 2.545469e-07
#> 21  1.00 6.049482e-02  53.640666     5 7.029159e-30 2.906191e-17 2.808023e-09
#> 22  1.25 4.566389e-02  56.275042     5 1.789119e-34 1.645343e-20 1.914634e-11
#> 23  1.50 3.238054e-02  58.624498     5 2.406061e-39 5.379924e-24 8.529518e-14
#> 24  1.75 2.157009e-02  60.679812     5 1.839603e-44 1.076124e-27 2.615692e-16
#> 25  2.00 1.349822e-02  62.437625     5 8.667755e-50 1.408680e-31 5.823896e-19
#> 26  2.25 7.935194e-03  63.907982     5 2.740372e-55 1.300251e-35 9.961582e-22
#> 27  2.50 4.382230e-03  65.116198     5 6.311039e-61 9.113146e-40 1.385279e-24
#> 28  2.75 2.273471e-03  66.097088     5 1.138001e-66 5.181708e-44 1.649561e-27
#> 29  3.00 1.108001e-03  66.887482     5 1.704715e-72 2.524923e-48 1.757377e-30
#> 30  3.25 5.072800e-04  67.521427     5 2.222555e-78 1.101209e-52 1.734857e-33
#> 31  3.50 2.181784e-04  68.028320     5 2.614151e-84 4.446247e-57 1.630519e-36
#> 32  3.75 8.815204e-05  68.432707     5 2.851046e-90 1.705679e-61 1.489018e-39
#> 33  4.00 3.345874e-05  68.754747     5 2.944232e-96 6.341813e-66 1.341229e-42
#>       p.level.4    p.level.5
#> 1  1.213588e-08 1.930802e-23
#> 2  3.304127e-08 1.524485e-22
#> 3  1.005527e-07 1.398092e-21
#> 4  3.445686e-07 1.522124e-20
#> 5  1.334044e-06 2.016135e-19
#> 6  5.816963e-06 3.328189e-18
#> 7  2.820067e-05 6.976833e-17
#> 8  1.480938e-04 1.869960e-15
#> 9  8.073718e-04 6.326350e-14
#> 10 4.299374e-03 2.595055e-12
#> 11 2.067106e-02 1.199693e-10
#> 12 8.198086e-02 5.627131e-09
#> 13 2.458526e-01 2.359746e-07
#> 14 5.259750e-01 7.766216e-06
#> 15 8.041392e-01 1.782156e-04
#> 16 9.507335e-01 2.589605e-03
#> 17 9.717552e-01 2.231446e-02
#> 18 8.882987e-01 1.113118e-01
#> 19 6.702547e-01 3.297318e-01
#> 20 3.743574e-01 6.256424e-01
#> 21 1.451413e-01 8.548586e-01
#> 22 3.864230e-02 9.613577e-01
#> 23 7.257495e-03 9.927425e-01
#> 24 1.012394e-03 9.989876e-01
#> 25 1.121609e-04 9.998878e-01
#> 26 1.060317e-05 9.999894e-01
#> 27 9.133383e-07 9.999991e-01
#> 28 7.563072e-08 9.999999e-01
#> 29 6.270749e-09 1.000000e+00
#> 30 5.351050e-10 1.000000e+00
#> 31 4.762291e-11 1.000000e+00
#> 32 4.422443e-12 1.000000e+00
#> 33 4.254596e-13 1.000000e+00
#> 
#> $cutscore
#> [1] 10 20 30 50
#> 

## -------------------------------------------------------------
## 2. When individual ability estimates are available (P method)
## -------------------------------------------------------------

# Randomly draw true ability values from N(0, 1)
set.seed(12)
theta <- rnorm(n = 1000, mean = 0, sd = 1)

# Simulate item response data
data <- simdat(x = x, theta = theta, D = 1)

# Estimate ability parameters using maximum likelihood (ML)
est_th <- est_score(
  x = x, data = data, D = 1, method = "ML",
  range = c(-4, 4), se = FALSE
)$est.theta

# Calculate classification accuracy and consistency
cac_2 <- cac_lee(x = x, cutscore = cutscore, theta = est_th, D = 1)

# Marginal indices and confusion matrix
cac_2$marginal
#>     level    accuracy consistency
#>         1 0.003088363 0.002997175
#>         2 0.038649080 0.033467238
#>         3 0.110618689 0.091887174
#>         4 0.442948314 0.404959837
#>         5 0.201415323 0.188708914
#>  marginal 0.796719769 0.722020337
cac_2$confusion
#>     Expected
#> True         1         2         3         4         5
#>    1 0.0030884 0.0029060 0.0000056 0.0000000 0.0000000
#>    2 0.0043690 0.0386491 0.0147554 0.0002265 0.0000000
#>    3 0.0000365 0.0207421 0.1106187 0.0406025 0.0000002
#>    4 0.0000000 0.0002357 0.0374015 0.4429483 0.0514145
#>    5 0.0000000 0.0000000 0.0000000 0.0305847 0.2014153

## ---------------------------------------------------------
## 3. When individual ability estimates are available,
##    but cut scores are specified on the IRT theta scale
## ---------------------------------------------------------
# Set the cut scores on the theta scale
cutscore <- c(-2, -0.4, 0.2, 1.0)

# Calculate classification accuracy and consistency
cac_3 <- cac_lee(
  x = x, cutscore = cutscore, theta = est_th, D = 1,
  cut.obs = FALSE
)

# Marginal indices and confusion matrix
cac_3$marginal
#>     level   accuracy consistency
#>         1 0.02778981  0.02547645
#>         2 0.25700457  0.23089145
#>         3 0.14645531  0.11728756
#>         4 0.17430960  0.14342915
#>         5 0.12397378  0.11530806
#>  marginal 0.72953307  0.63239267
cac_3$confusion
#>     Expected
#> True         1         2         3         4         5
#>    1 0.0277898 0.0072101 0.0000001 0.0000000 0.0000000
#>    2 0.0145173 0.2570046 0.0471682 0.0013095 0.0000003
#>    3 0.0000002 0.0423232 0.1464553 0.0488048 0.0004164
#>    4 0.0000000 0.0009595 0.0470240 0.1743096 0.0407069
#>    5 0.0000000 0.0000000 0.0001071 0.0199191 0.1239738
# }
```
