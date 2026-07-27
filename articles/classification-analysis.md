# Classification Accuracy and Consistency

## Overview

In many testing contexts, test scores are used to assign examinees to
performance categories — for example, pass/fail decisions or placement
into proficiency levels. Two indices quantify the quality of such
classification decisions:

- **Classification Accuracy (CA)**: the probability that an examinee is
  correctly classified based on their *true* latent ability.
- **Classification Consistency (CC)**: the probability that the same
  classification would result from two independent administrations of
  the same test.

Both indices are essential for supporting valid score interpretations in
high-stakes settings. Because data from repeated test administrations
are rarely available, single-administration estimation procedures based
on IRT models are commonly used.

**irtQ** implements two IRT-based estimation methods:

| Function | Method | Cut-score scale |
|----|----|----|
| [`cac_lee()`](https://hwangQ.github.io/irtQ/reference/cac_lee.md) | Lee (2010) — conditional summed-score distribution | Observed score or theta |
| [`cac_rud()`](https://hwangQ.github.io/irtQ/reference/cac_rud.md) | Rudner (2001); Rudner (2005) — normal approximation via test information function (TIF) | Theta scale only |

Both functions support two estimation approaches:

- **D-method (distribution-based)**: integrates conditional indices over
  a population ability distribution specified by quadrature points and
  weights. This approach is appropriate when the empirical ability
  distribution (e.g., posterior quadrature weights from IRT calibration)
  is available.
- **P-method (person-based)**: averages conditional indices over
  individual ability estimates from a sample of examinees.

Both functions return a list with five elements:

- `$confusion`: a confusion matrix (rows = true level, columns =
  expected level),
- `$marginal`: a data frame of marginal CA and CC indices per
  performance level plus an overall marginal row,
- `$conditional`: a data frame of conditional CA and CC for each theta
  (or node),
- `$prob.level`: a data frame of the probability of being assigned to
  each level for each theta (or node),
- `$cutscore`: the cut scores used in the analysis (on the observed
  score scale).

``` r

library(irtQ)
set.seed(2026)
```

------------------------------------------------------------------------

## Setup: Common Item Metadata and Response Data

We use two test forms throughout this vignette.

### Binary test (20 items, 3PLM)

A 20-item dichotomous test calibrated under the three-parameter logistic
model (3PLM).

``` r

meta_bin <- shape_df(
  par.drm = list(
    a = c(1.0, 1.2, 0.8, 1.4, 0.9, 1.1, 1.3, 0.7, 1.0, 1.2,
          0.9, 1.1, 1.4, 0.85, 1.0, 1.2, 0.8, 1.1, 1.3, 0.9),
    b = c(-2.0, -1.5, -1.0, -0.6, -0.2, 0.0, 0.4, 0.8, 1.1, 1.5,
          -1.3, -0.4,  0.5,  1.0,  1.8, -0.8, 0.2, 0.7, -1.1, 1.3),
    g = rep(0.15, 20)
  ),
  item.id = paste0("ITEM", 1:20),
  cats    = 2,
  model   = "3PLM"
)

theta_bin <- rnorm(1000, mean = 0, sd = 1)
resp_bin  <- simdat(x = meta_bin, theta = theta_bin, D = 1.702)
```

``` r

mod_bin <- est_irt(
  data       = resp_bin,
  D          = 1.702,
  model      = "3PLM",
  cats       = 2,
  item.id    = paste0("ITEM", 1:20),
  use.gprior = TRUE,
  gprior     = list(dist = "beta", params = c(4, 16)),
  EmpHist    = FALSE,
  verbose    = FALSE
)
meta_cal_bin <- mod_bin$par.est

# ML ability estimates (needed for the P-method)
score_bin <- est_score(
  x      = meta_cal_bin,
  data   = resp_bin,
  D      = 1.702,
  method = "ML",
  range  = c(-5, 5),
  se     = TRUE
)
```

### Mixed-format test (15 × 3PLM + 5 × GRM, 4 categories)

A mixed-format test combining 15 dichotomous (3PLM) and 5 polytomous
items (Graded Response Model, GRM; 4 ordered categories). The maximum
possible observed summed score is $`15 \times 1 + 5 \times 3 = 30`$.

``` r

meta_mix <- shape_df(
  par.drm = list(
    a = c(1.0, 1.2, 0.9, 1.1, 0.8, 1.3, 1.0, 0.9, 1.1, 1.2,
          0.85, 1.1, 1.3, 0.9, 1.0),
    b = c(-1.5, -1.0, -0.5, -0.2, 0.2, 0.5, 0.8, 1.1, 1.4, -1.2,
          -0.4,  0.1,  0.6,  1.0, -0.8),
    g = rep(0.15, 15)
  ),
  par.prm = list(
    a = c(1.5, 1.2, 1.0, 1.3, 0.9),
    d = list(
      c(-1.5, -0.3,  0.9),
      c(-1.2,  0.0,  1.1),
      c(-0.9,  0.4,  1.4),
      c(-1.1, -0.1,  1.0),
      c(-1.3,  0.3,  1.2)
    )
  ),
  item.id = c(paste0("DRM", 1:15), paste0("GRM", 1:5)),
  cats    = c(rep(2, 15), rep(4, 5)),
  model   = c(rep("3PLM", 15), rep("GRM", 5))
)

theta_mix <- rnorm(1000, mean = 0, sd = 1)
resp_mix  <- simdat(x = meta_mix, theta = theta_mix, D = 1.702)
```

``` r

mod_mix <- est_irt(
  data       = resp_mix,
  D          = 1.702,
  model      = c(rep("3PLM", 15), rep("GRM", 5)),
  cats       = c(rep(2, 15), rep(4, 5)),
  item.id    = c(paste0("DRM", 1:15), paste0("GRM", 1:5)),
  use.gprior = TRUE,
  gprior     = list(dist = "beta", params = c(4, 16)),
  EmpHist    = FALSE,
  verbose    = FALSE
)
meta_cal_mix <- mod_mix$par.est

# ML ability estimates
score_mix <- est_score(
  x      = meta_cal_mix,
  data   = resp_mix,
  D      = 1.702,
  method = "ML",
  range  = c(-5, 5),
  se     = TRUE
)
```

### Shared quadrature grid

Both methods use the same quadrature grid for D-method examples.

``` r

quad_nodes   <- seq(-4, 4, by = 0.25)
quad_weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = quad_nodes)
```

------------------------------------------------------------------------

## Part 1: Lee’s Method — `cac_lee()`

### Method overview

Lee (2010) proposed a general IRT framework for estimating CA and CC
indices from a single test administration. The method can be applied to
tests consisting of dichotomous items, polytomous items, or a mixture of
both.

The key quantity is the **conditional summed-score distribution**,
$`\Pr(X = x \mid \theta)`$, computed using the Lord–Wingersky recursive
algorithm (Lord and Wingersky 1984; Kolen and Brennan 2004). This
distribution gives the probability that an examinee with ability
$`\theta`$ obtains each possible observed summed score $`x`$.

**Cut scores** partition the range of observed summed scores into $`K`$
performance levels. The probability that an examinee at ability
$`\theta`$ is assigned to level $`k`$ is:
``` math
p_\theta(k) = \sum_{x \in \text{level } k} \Pr(X = x \mid \theta).
```

The **conditional CA index** is defined as the probability that the
examinee is correctly classified, that is, assigned to the level that
matches their true performance level $`\eta(\theta)`$:
``` math
\gamma_\theta = p_\theta(\eta(\theta)).
```
The true performance level $`\eta(\theta)`$ is determined by comparing
the expected summed score at $`\theta`$ to the true cut scores.

The **conditional CC index** is defined as the probability that two
independent test administrations yield the same classification decision:
``` math
\phi_\theta = \sum_{k=1}^{K} p_\theta(k)^2.
```

**Marginal** CA and CC indices are obtained by integrating the
conditional indices over the ability distribution $`g(\theta)`$:
``` math
\gamma = \int_{-\infty}^{\infty} \gamma_\theta \, g(\theta) \, d\theta, \qquad
\phi   = \int_{-\infty}^{\infty} \phi_\theta   \, g(\theta) \, d\theta.
```

If cut scores are specified on the **theta scale** (i.e.,
`cut.obs = FALSE`), they are internally converted to the observed
summed-score scale using the Test Characteristic Curve (TCC):
$`E(X \mid \theta^*) = \sum_i \sum_j j \, \Pr(U_i = j \mid \theta^*)`$.

where:

- $`X`$ is the total observed summed score across the entire test.
- $`\theta^*`$ is a specific cut score specified on the latent ability
  scale.
- $`i`$ is the index for items ($`i = 1, 2, \dots, N`$, where $`N`$ is
  the total number of items).
- $`j`$ is the index for the item score categories
  ($`j = 0, 1, \dots, K_i - 1`$, where $`K_i`$ is the number of response
  categories for item $`i`$). For dichotomous items, $`j \in \{0, 1\}`$.
- $`U_i`$ is the item response random variable for item $`i`$.
- $`\Pr(U_i = j \mid \theta^*)`$ is the probability of obtaining a score
  category $`j`$ on item $`i`$ given the ability level $`\theta^*`$,
  which is modeled by a specific IRT item response function (e.g., 3PLM
  or GRM).

### Key arguments in `cac_lee()`

| Argument | Description |
|----|----|
| `x` | Item metadata data frame |
| `cutscore` | Numeric vector of cut scores; defines $`K + 1`$ performance levels |
| `weights` | Two-column matrix of quadrature nodes and weights for D-method (use [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)) |
| `theta` | Numeric vector of individual ability estimates for P-method |
| `D` | Scaling constant — must match the value used during calibration |
| `cut.obs` | `TRUE` (default): cut scores on the observed summed-score scale; `FALSE`: on the theta scale (converted internally via TCC) |

Either `theta` or `weights` must be provided (but not both. If both are
supplied, `weights` takes priority and the D-method is applied).

------------------------------------------------------------------------

### Example 1: Binary test with `cac_lee()`

#### (1) D-method: integrate over the population ability distribution

Provide a quadrature grid via `weights`. Use
[`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)
to construct weights from a parametric distribution.

``` r

# Two cut scores on the observed summed-score scale (range: 0–20)
# Defines three performance levels: [0,8), [8,14), [14,20]
cutscore_obs <- c(8, 14)

# D-method: integrate over the N(0,1) quadrature grid
cac_l_d <- cac_lee(
  x        = meta_cal_bin,
  cutscore = cutscore_obs,
  weights  = quad_weights,     # quadrature grid (D-method)
  D        = 1.702,
  cut.obs  = TRUE              # cut scores on observed scale (default)
)

cac_l_d
#> $confusion
#>     Expected
#> True         1         2         3
#>    1 0.1379361 0.0521427 0.0000748
#>    2 0.0468864 0.4078148 0.0897047
#>    3 0.0000015 0.0271701 0.2382690
#> 
#> $marginal
#>     level  accuracy consistency
#>         1 0.1379361   0.1241428
#>         2 0.4078148   0.3581910
#>         3 0.2382690   0.2227162
#>  marginal 0.7840199   0.7050500
#> 
#> $conditional
#>    theta      weights true.score level  accuracy consistency
#> 1  -4.00 3.345874e-05   3.341994     1 0.9889093   0.9780645
#> 2  -3.75 8.815204e-05   3.380638     1 0.9881519   0.9765846
#> 3  -3.50 2.181784e-04   3.436371     1 0.9870018   0.9743415
#> 4  -3.25 5.072800e-04   3.516191     1 0.9852294   0.9708951
#> 5  -3.00 1.108001e-03   3.629290     1 0.9824451   0.9655065
#> 6  -2.75 2.273471e-03   3.787249     1 0.9779598   0.9568910
#> 7  -2.50 4.382230e-03   4.004115     1 0.9704978   0.9427364
#> 8  -2.25 7.935194e-03   4.296454     1 0.9576106   0.9188147
#> 9  -2.00 1.349822e-02   4.683088     1 0.9345632   0.8776896
#> 10 -1.75 2.157009e-02   5.182933     1 0.8926648   0.8083675
#> 11 -1.50 3.238054e-02   5.808922     1 0.8182603   0.7025563
#> 12 -1.25 4.566389e-02   6.560066     1 0.6965335   0.5771006
#> 13 -1.00 6.049482e-02   7.418926     1 0.5250879   0.5003034
#> 14 -0.75 7.528702e-02   8.357367     2 0.6653491   0.5520526
#> 15 -0.50 8.801945e-02   9.345581     2 0.8221135   0.7028173
#> 16 -0.25 9.667045e-02  10.361267     2 0.8944358   0.8056876
#> 17  0.00 9.973910e-02  11.395210     2 0.8606205   0.7562848
#> 18  0.25 9.667045e-02  12.445364     2 0.7189256   0.5943970
#> 19  0.50 8.801945e-02  13.502597     2 0.4948773   0.4997655
#> 20  0.75 7.528702e-02  14.546417     3 0.7379663   0.6132460
#> 21  1.00 6.049482e-02  15.555088     3 0.8997831   0.8196528
#> 22  1.25 4.566389e-02  16.505175     3 0.9734804   0.9483673
#> 23  1.50 3.238054e-02  17.356735     3 0.9951853   0.9904170
#> 24  1.75 2.157009e-02  18.065625     3 0.9993754   0.9987516
#> 25  2.00 1.349822e-02  18.616683     3 0.9999375   0.9998750
#> 26  2.25 7.935194e-03  19.026191     3 0.9999948   0.9999895
#> 27  2.50 4.382230e-03  19.322349     3 0.9999996   0.9999992
#> 28  2.75 2.273471e-03  19.532537     3 1.0000000   0.9999999
#> 29  3.00 1.108001e-03  19.679548     3 1.0000000   1.0000000
#> 30  3.25 5.072800e-04  19.781228     3 1.0000000   1.0000000
#> 31  3.50 2.181784e-04  19.850993     3 1.0000000   1.0000000
#> 32  3.75 8.815204e-05  19.898610     3 1.0000000   1.0000000
#> 33  4.00 3.345874e-05  19.931010     3 1.0000000   1.0000000
#> 
#> $prob.level
#>    theta      weights true.score level    p.level.1    p.level.2    p.level.3
#> 1  -4.00 3.345874e-05   3.341994     1 9.889093e-01 1.109058e-02 1.582213e-07
#> 2  -3.75 8.815204e-05   3.380638     1 9.881519e-01 1.184792e-02 1.808832e-07
#> 3  -3.50 2.181784e-04   3.436371     1 9.870018e-01 1.299796e-02 2.178853e-07
#> 4  -3.25 5.072800e-04   3.516191     1 9.852294e-01 1.477031e-02 2.810713e-07
#> 5  -3.00 1.108001e-03   3.629290     1 9.824451e-01 1.755450e-02 3.957164e-07
#> 6  -2.75 2.273471e-03   3.787249     1 9.779598e-01 2.203963e-02 6.213288e-07
#> 7  -2.50 4.382230e-03   4.004115     1 9.704978e-01 2.950104e-02 1.114837e-06
#> 8  -2.25 7.935194e-03   4.296454     1 9.576106e-01 4.238708e-02 2.344988e-06
#> 9  -2.00 1.349822e-02   4.683088     1 9.345632e-01 6.543090e-02 5.901569e-06
#> 10 -1.75 2.157009e-02   5.182933     1 8.926648e-01 1.073174e-01 1.782989e-05
#> 11 -1.50 3.238054e-02   5.808922     1 8.182603e-01 1.816767e-01 6.307509e-05
#> 12 -1.25 4.566389e-02   6.560066     1 6.965335e-01 3.032188e-01 2.476987e-04
#> 13 -1.00 6.049482e-02   7.418926     1 5.250879e-01 4.739040e-01 1.008033e-03
#> 14 -0.75 7.528702e-02   8.357367     2 3.306772e-01 6.653491e-01 3.973727e-03
#> 15 -0.50 8.801945e-02   9.345581     2 1.635252e-01 8.221135e-01 1.436128e-02
#> 16 -0.25 9.667045e-02  10.361267     2 5.986400e-02 8.944358e-01 4.570019e-02
#> 17  0.00 9.973910e-02  11.395210     2 1.535824e-02 8.606205e-01 1.240213e-01
#> 18  0.25 9.667045e-02  12.445364     2 2.621303e-03 7.189256e-01 2.784531e-01
#> 19  0.50 8.801945e-02  13.502597     2 2.842548e-04 4.948773e-01 5.048384e-01
#> 20  0.75 7.528702e-02  14.546417     3 1.896328e-05 2.620147e-01 7.379663e-01
#> 21  1.00 6.049482e-02  15.555088     3 7.663720e-07 1.002162e-01 8.997831e-01
#> 22  1.25 4.566389e-02  16.505175     3 1.876512e-08 2.651963e-02 9.734804e-01
#> 23  1.50 3.238054e-02  17.356735     3 2.855207e-10 4.814701e-03 9.951853e-01
#> 24  1.75 2.157009e-02  18.065625     3 2.879399e-12 6.245660e-04 9.993754e-01
#> 25  2.00 1.349822e-02  18.616683     3 2.112469e-14 6.251190e-05 9.999375e-01
#> 26  2.25 7.935194e-03  19.026191     3 1.237687e-16 5.238462e-06 9.999948e-01
#> 27  2.50 4.382230e-03  19.322349     3 6.247592e-19 3.928072e-07 9.999996e-01
#> 28  2.75 2.273471e-03  19.532537     3 2.874294e-21 2.763536e-08 1.000000e+00
#> 29  3.00 1.108001e-03  19.679548     3 1.254862e-23 1.882065e-09 1.000000e+00
#> 30  3.25 5.072800e-04  19.781228     3 5.349517e-26 1.265012e-10 1.000000e+00
#> 31  3.50 2.181784e-04  19.850993     3 2.271798e-28 8.486486e-12 1.000000e+00
#> 32  3.75 8.815204e-05  19.898610     3 9.743935e-31 5.717070e-13 1.000000e+00
#> 33  4.00 3.345874e-05  19.931010     3 4.259757e-33 3.879142e-14 1.000000e+00
#> 
#> $cutscore
#> [1]  8 14
```

The output contains:

- `$confusion`: confusion matrix (rows = true performance level, columns
  = expected level under CA). The diagonal entries correspond to correct
  classifications.
- `$marginal`: marginal CA and CC per level, plus an overall row
  labelled `"marginal"`.
- `$cutscore`: the cut scores used (always on the observed summed-score
  scale).

#### (2) P-method: average over individual ability estimates

Provide individual ML ability estimates via the `theta` argument.

``` r

# P-method: average over individual ML ability estimates
cac_l_p <- cac_lee(
  x        = meta_cal_bin,
  cutscore = cutscore_obs,
  theta    = score_bin$est.theta,   # individual estimates
  D        = 1.702,
  cut.obs  = TRUE
)

# Confusion matrix (rows = true level, columns = expected level)
cac_l_p$confusion
#>     Expected
#> True         1         2         3
#>    1 0.1639447 0.0579551 0.0001002
#>    2 0.0368065 0.3680726 0.0951209
#>    3 0.0000019 0.0257336 0.2522645

# Marginal CA and CC indices
cac_l_p$marginal
#>     level  accuracy consistency
#>         1 0.1639447   0.1507562
#>         2 0.3680726   0.3237276
#>         3 0.2522645   0.2380921
#>  marginal 0.7842818   0.7125758
```

#### (3) Theta-scale cut scores

When cut scores are naturally expressed on the IRT theta scale (e.g.,
standard-setting results), set `cut.obs = FALSE`. The function converts
them internally to the observed summed-score scale using the TCC before
computing CA and CC.

``` r

# Cut scores on the theta scale
cutscore_theta <- c(-0.5, 0.8)

cac_l_theta <- cac_lee(
  x        = meta_cal_bin,
  cutscore = cutscore_theta,
  theta    = score_bin$est.theta,
  D        = 1.702,
  cut.obs  = FALSE         # cut scores on the theta scale → converted via TCC
)

# Converted cut scores (now on the observed summed-score scale)
cac_l_theta$cutscore
#> [1]  9.345581 14.751596

# Marginal CA and CC indices
cac_l_theta$marginal
#>     level  accuracy consistency
#>         1 0.2709956   0.2504939
#>         2 0.3194857   0.2738093
#>         3 0.2021494   0.1889719
#>  marginal 0.7926307   0.7132751
```

------------------------------------------------------------------------

### Example 2: Mixed-format test with `cac_lee()`

For a mixed-format test (dichotomous + polytomous items), the maximum
possible observed summed score is $`\sum_j (K_j - 1)`$, where $`K_j`$ is
the number of score categories for item $`j`$. In our test: 15 binary
items (max 15) + 5 four-category items (max $`3 \times 5 = 15`$) → max =
30.

``` r

# Cut scores on the observed summed-score scale (range: 0–30)
cutscore_mix <- c(10, 22)

# D-method
cac_l_mix_d <- cac_lee(
  x        = meta_cal_mix,
  cutscore = cutscore_mix,
  weights  = quad_weights,
  D        = 1.702,
  cut.obs  = TRUE
)
cac_l_mix_d
#> $confusion
#>     Expected
#> True         1         2         3
#>    1 0.1581905 0.0319630 0.0000000
#>    2 0.0364469 0.4728394 0.0351196
#>    3 0.0000000 0.0401439 0.2252967
#> 
#> $marginal
#>     level  accuracy consistency
#>         1 0.1581905   0.1459529
#>         2 0.4728394   0.4356591
#>         3 0.2252967   0.2114101
#>  marginal 0.8563266   0.7930221
#> 
#> $conditional
#>    theta      weights true.score level  accuracy consistency
#> 1  -4.00 3.345874e-05   2.253728     1 0.9999875   0.9999751
#> 2  -3.75 8.815204e-05   2.283350     1 0.9999837   0.9999675
#> 3  -3.50 2.181784e-04   2.329305     1 0.9999768   0.9999536
#> 4  -3.25 5.072800e-04   2.400715     1 0.9999632   0.9999264
#> 5  -3.00 1.108001e-03   2.511670     1 0.9999334   0.9998668
#> 6  -2.75 2.273471e-03   2.683558     1 0.9998593   0.9997186
#> 7  -2.50 4.382230e-03   2.947785     1 0.9996472   0.9992946
#> 8  -2.25 7.935194e-03   3.347748     1 0.9989487   0.9978996
#> 9  -2.00 1.349822e-02   3.937048     1 0.9963770   0.9927802
#> 10 -1.75 2.157009e-02   4.769229     1 0.9864709   0.9733078
#> 11 -1.50 3.238054e-02   5.877340     1 0.9506156   0.9061089
#> 12 -1.25 4.566389e-02   7.253388     1 0.8443815   0.7371973
#> 13 -1.00 6.049482e-02   8.847094     1 0.6213434   0.5294483
#> 14 -0.75 7.528702e-02  10.591733     2 0.6716724   0.5589388
#> 15 -0.50 8.801945e-02  12.435242     2 0.8913060   0.8062151
#> 16 -0.25 9.667045e-02  14.344436     2 0.9780924   0.9570785
#> 17  0.00 9.973910e-02  16.290396     2 0.9831067   0.9667242
#> 18  0.25 9.667045e-02  18.243911     2 0.9155823   0.8453995
#> 19  0.50 8.801945e-02  20.174358     2 0.7123696   0.5902000
#> 20  0.75 7.528702e-02  22.037047     3 0.6031743   0.5212898
#> 21  1.00 6.049482e-02  23.764872     3 0.8577494   0.7559693
#> 22  1.25 4.566389e-02  25.283987     3 0.9675098   0.9371309
#> 23  1.50 3.238054e-02  26.543794     3 0.9949096   0.9898710
#> 24  1.75 2.157009e-02  27.533916     3 0.9993852   0.9987712
#> 25  2.00 1.349822e-02  28.278053     3 0.9999355   0.9998710
#> 26  2.25 7.935194e-03  28.817653     3 0.9999935   0.9999870
#> 27  2.50 4.382230e-03  29.198129     3 0.9999993   0.9999987
#> 28  2.75 2.273471e-03  29.460788     3 0.9999999   0.9999999
#> 29  3.00 1.108001e-03  29.639380     3 1.0000000   1.0000000
#> 30  3.25 5.072800e-04  29.759569     3 1.0000000   1.0000000
#> 31  3.50 2.181784e-04  29.839930     3 1.0000000   1.0000000
#> 32  3.75 8.815204e-05  29.893462     3 1.0000000   1.0000000
#> 33  4.00 3.345874e-05  29.929057     3 1.0000000   1.0000000
#> 
#> $prob.level
#>    theta      weights true.score level    p.level.1    p.level.2    p.level.3
#> 1  -4.00 3.345874e-05   2.253728     1 9.999875e-01 1.246405e-05 6.654498e-21
#> 2  -3.75 8.815204e-05   2.283350     1 9.999837e-01 1.626801e-05 3.455029e-20
#> 3  -3.50 2.181784e-04   2.329305     1 9.999768e-01 2.317736e-05 2.029941e-19
#> 4  -3.25 5.072800e-04   2.400715     1 9.999632e-01 3.679529e-05 1.387226e-18
#> 5  -3.00 1.108001e-03   2.511670     1 9.999334e-01 6.660666e-05 1.130533e-17
#> 6  -2.75 2.273471e-03   2.683558     1 9.998593e-01 1.406955e-04 1.123964e-16
#> 7  -2.50 4.382230e-03   2.947785     1 9.996472e-01 3.528008e-04 1.389746e-15
#> 8  -2.25 7.935194e-03   3.347748     1 9.989487e-01 1.051313e-03 2.163149e-14
#> 9  -2.00 1.349822e-02   3.937048     1 9.963770e-01 3.623040e-03 4.234861e-13
#> 10 -1.75 2.157009e-02   4.769229     1 9.864709e-01 1.352915e-02 1.020054e-11
#> 11 -1.50 3.238054e-02   5.877340     1 9.506156e-01 4.938435e-02 2.866605e-10
#> 12 -1.25 4.566389e-02   7.253388     1 8.443815e-01 1.556185e-01 8.611461e-09
#> 13 -1.00 6.049482e-02   8.847094     1 6.213434e-01 3.786563e-01 2.479428e-07
#> 14 -0.75 7.528702e-02  10.591733     2 3.283215e-01 6.716724e-01 6.135126e-06
#> 15 -0.50 8.801945e-02  12.435242     2 1.085756e-01 8.913060e-01 1.183608e-04
#> 16 -0.25 9.667045e-02  14.344436     2 2.027586e-02 9.780924e-01 1.631751e-03
#> 17  0.00 9.973910e-02  16.290396     2 2.017724e-03 9.831067e-01 1.487555e-02
#> 18  0.25 9.667045e-02  18.243911     2 1.056265e-04 9.155823e-01 8.431207e-02
#> 19  0.50 8.801945e-02  20.174358     2 2.985656e-06 7.123696e-01 2.876274e-01
#> 20  0.75 7.528702e-02  22.037047     3 4.838632e-08 3.968257e-01 6.031743e-01
#> 21  1.00 6.049482e-02  23.764872     3 4.932609e-10 1.422506e-01 8.577494e-01
#> 22  1.25 4.566389e-02  25.283987     3 3.562643e-12 3.249018e-02 9.675098e-01
#> 23  1.50 3.238054e-02  26.543794     3 2.070954e-14 5.090413e-03 9.949096e-01
#> 24  1.75 2.157009e-02  27.533916     3 1.085382e-16 6.147842e-04 9.993852e-01
#> 25  2.00 1.349822e-02  28.278053     3 5.574894e-19 6.452184e-05 9.999355e-01
#> 26  2.25 7.935194e-03  28.817653     3 2.944402e-21 6.493681e-06 9.999935e-01
#> 27  2.50 4.382230e-03  29.198129     3 1.627753e-23 6.712883e-07 9.999993e-01
#> 28  2.75 2.273471e-03  29.460788     3 9.409431e-26 7.428909e-08 9.999999e-01
#> 29  3.00 1.108001e-03  29.639380     3 5.638811e-28 8.994980e-09 1.000000e+00
#> 30  3.25 5.072800e-04  29.759569     3 3.470598e-30 1.203876e-09 1.000000e+00
#> 31  3.50 2.181784e-04  29.839930     3 2.177506e-32 1.787926e-10 1.000000e+00
#> 32  3.75 8.815204e-05  29.893462     3 1.385290e-34 2.944311e-11 1.000000e+00
#> 33  4.00 3.345874e-05  29.929057     3 8.904048e-37 5.346627e-12 1.000000e+00
#> 
#> $cutscore
#> [1] 10 22

# P-method
cac_l_mix_p <- cac_lee(
  x        = meta_cal_mix,
  cutscore = cutscore_mix,
  theta    = score_mix$est.theta,
  D        = 1.702,
  cut.obs  = TRUE
)
cac_l_mix_p$marginal
#>     level  accuracy consistency
#>         1 0.1702082   0.1598568
#>         2 0.4707304   0.4351404
#>         3 0.2145511   0.2030033
#>  marginal 0.8554897   0.7980006
```

------------------------------------------------------------------------

## Part 2: Rudner’s Method — `cac_rud()`

### Method overview

Rudner (2001) and Rudner (2005) proposed a simpler approach based on the
assumption that IRT ability estimates are normally distributed around
the true ability:
``` math
\hat{\theta} \mid \theta \sim N\!\left(\theta,\; \mathrm{SE}^2(\theta)\right),
```
where $`\mathrm{SE}(\theta) = 1 / \sqrt{I(\theta)}`$ is the standard
error of estimation derived from the TIF.

Under this assumption, the probability that an examinee with true
ability $`\theta`$ obtains an ability estimate falling in the cut-score
interval $`[c_{k-1}, c_k)`$ is:
``` math
p_\theta(k) = \Phi\!\left(\frac{c_k - \theta}{\mathrm{SE}(\theta)}\right)
            - \Phi\!\left(\frac{c_{k-1} - \theta}{\mathrm{SE}(\theta)}\right),
```
where $`\Phi(\cdot)`$ is the standard normal cumulative distribution
function, and the boundary cut scores are $`c_0 = -\infty`$ and
$`c_K = +\infty`$.

Conditional CA and CC are then computed in the same way as in Lee’s
method:
``` math
\gamma_\theta = p_\theta(\eta(\theta)), \qquad
\phi_\theta = \sum_{k=1}^{K} p_\theta(k)^2,
```
and marginal indices are obtained by averaging over the ability
distribution.

> **Important**: Unlike
> [`cac_lee()`](https://hwangQ.github.io/irtQ/reference/cac_lee.md), cut
> scores for
> [`cac_rud()`](https://hwangQ.github.io/irtQ/reference/cac_rud.md) must
> **always** be specified on the **IRT theta scale** (not the observed
> summed-score scale).

Standard errors can be supplied in two ways:

1.  Pass item metadata via `x` — the function computes SE from the TIF
    internally.
2.  Pass a pre-computed SE vector via `se` (same length as `theta` for
    P-method, or same length as the number of quadrature nodes for
    D-method).

Either `x` or `se` must be provided.

### Key arguments in `cac_rud()`

| Argument | Description |
|----|----|
| `x` | Item metadata data frame (optional if `se` is provided; used to compute SE from TIF) |
| `cutscore` | Numeric vector of cut scores on the **theta scale** |
| `theta` | Numeric vector of individual ability estimates (P-method) |
| `se` | Numeric vector of standard errors. If `NULL` and `x` is provided, SE is computed from the TIF |
| `weights` | Two-column matrix of quadrature nodes and weights (D-method) |
| `D` | Scaling constant |

------------------------------------------------------------------------

### Example 3: Binary test with `cac_rud()`

#### (1) D-method: SE computed from item metadata

``` r

# Cut scores on the theta scale (matching the metric used in cac_lee theta examples)
cutscore_th <- c(-0.5, 0.8)

# D-method: SE computed internally from the TIF using item metadata
cac_r_d <- cac_rud(
  x        = meta_cal_bin,     # item metadata → SE via TIF
  cutscore = cutscore_th,
  weights  = quad_weights,
  D        = 1.702
)

cac_r_d
#> $confusion
#>     Expected
#> True         1         2         3
#>    1 0.2303756 0.0349129 0.0001521
#>    2 0.0850413 0.3943581 0.0650065
#>    3 0.0000121 0.0268354 0.1633060
#> 
#> $marginal
#>     level  accuracy consistency
#>         1 0.2303756   0.2098833
#>         2 0.3943581   0.3479965
#>         3 0.1633060   0.1494771
#>  marginal 0.7880397   0.7073569
#> 
#> $conditional
#>    theta      weights level  accuracy consistency
#> 1  -4.00 3.345874e-05     1 0.7073627   0.5562361
#> 2  -3.75 8.815204e-05     1 0.7611902   0.6112320
#> 3  -3.50 2.181784e-04     1 0.8163547   0.6833776
#> 4  -3.25 5.072800e-04     1 0.8663885   0.7600463
#> 5  -3.00 1.108001e-03     1 0.9061324   0.8266633
#> 6  -2.75 2.273471e-03     1 0.9341427   0.8759656
#> 7  -2.50 4.382230e-03     1 0.9524072   0.9090811
#> 8  -2.25 7.935194e-03     1 0.9640000   0.9305317
#> 9  -2.00 1.349822e-02     1 0.9707328   0.9431667
#> 10 -1.75 2.157009e-02     1 0.9720558   0.9456707
#> 11 -1.50 3.238054e-02     1 0.9643817   0.9312996
#> 12 -1.25 4.566389e-02     1 0.9373092   0.8824770
#> 13 -1.00 6.049482e-02     1 0.8665594   0.7687230
#> 14 -0.75 7.528702e-02     1 0.7198456   0.5965787
#> 15 -0.50 8.801945e-02     2 0.4990237   0.4990256
#> 16 -0.25 9.667045e-02     2 0.7208676   0.5944490
#> 17  0.00 9.973910e-02     2 0.8639771   0.7593659
#> 18  0.25 9.667045e-02     2 0.8838826   0.7894726
#> 19  0.50 8.801945e-02     2 0.7697276   0.6429011
#> 20  0.75 7.528702e-02     2 0.5496171   0.5042292
#> 21  1.00 6.049482e-02     3 0.6941339   0.5753329
#> 22  1.25 4.566389e-02     3 0.8709842   0.7752572
#> 23  1.50 3.238054e-02     3 0.9523281   0.9092013
#> 24  1.75 2.157009e-02     3 0.9787137   0.9583335
#> 25  2.00 1.349822e-02     3 0.9862584   0.9728945
#> 26  2.25 7.935194e-03     3 0.9877023   0.9757069
#> 27  2.50 4.382230e-03     3 0.9862131   0.9728050
#> 28  2.75 2.273471e-03     3 0.9821978   0.9650213
#> 29  3.00 1.108001e-03     3 0.9751750   0.9515397
#> 30  3.25 5.072800e-04     3 0.9643879   0.9311237
#> 31  3.50 2.181784e-04     3 0.9491804   0.9028654
#> 32  3.75 8.815204e-05     3 0.9292921   0.8667527
#> 33  4.00 3.345874e-05     3 0.9050123   0.8239989
#> 
#> $prob.level
#>    theta      weights level    p.level.1  p.level.2    p.level.3
#> 1  -4.00 3.345874e-05     1 7.073627e-01 0.06552293 2.271144e-01
#> 2  -3.75 8.815204e-05     1 7.611902e-01 0.07874482 1.600650e-01
#> 3  -3.50 2.181784e-04     1 8.163547e-01 0.08550761 9.813767e-02
#> 4  -3.25 5.072800e-04     1 8.663885e-01 0.08247809 5.113342e-02
#> 5  -3.00 1.108001e-03     1 9.061324e-01 0.07124185 2.262575e-02
#> 6  -2.75 2.273471e-03     1 9.341427e-01 0.05716104 8.696281e-03
#> 7  -2.50 4.382230e-03     1 9.524072e-01 0.04464277 2.950054e-03
#> 8  -2.25 7.935194e-03     1 9.640000e-01 0.03514250 8.574977e-04
#> 9  -2.00 1.349822e-02     1 9.707328e-01 0.02906033 2.068985e-04
#> 10 -1.75 2.157009e-02     1 9.720558e-01 0.02789613 4.803937e-05
#> 11 -1.50 3.238054e-02     1 9.643817e-01 0.03560159 1.668758e-05
#> 12 -1.25 4.566389e-02     1 9.373092e-01 0.06267677 1.400760e-05
#> 13 -1.00 6.049482e-02     1 8.665594e-01 0.13340851 3.207722e-05
#> 14 -0.75 7.528702e-02     1 7.198456e-01 0.28000177 1.526424e-04
#> 15 -0.50 8.801945e-02     2 5.000000e-01 0.49902373 9.762737e-04
#> 16 -0.25 9.667045e-02     2 2.734345e-01 0.72086763 5.697842e-03
#> 17  0.00 9.973910e-02     2 1.107800e-01 0.86397714 2.524289e-02
#> 18  0.25 9.667045e-02     2 3.083264e-02 0.88388261 8.528475e-02
#> 19  0.50 8.801945e-02     2 5.802320e-03 0.76972763 2.244701e-01
#> 20  0.75 7.528702e-02     2 7.723930e-04 0.54961709 4.496105e-01
#> 21  1.00 6.049482e-02     3 7.032651e-05 0.30579582 6.941339e-01
#> 22  1.25 4.566389e-02     3 5.448798e-06 0.12901031 8.709842e-01
#> 23  1.50 3.238054e-02     3 9.428157e-07 0.04767095 9.523281e-01
#> 24  1.75 2.157009e-02     3 7.820366e-07 0.02128554 9.787137e-01
#> 25  2.00 1.349822e-02     3 2.185770e-06 0.01373939 9.862584e-01
#> 26  2.25 7.935194e-03     3 1.009047e-05 0.01228757 9.877023e-01
#> 27  2.50 4.382230e-03     3 5.050080e-05 0.01373637 9.862131e-01
#> 28  2.75 2.273471e-03     3 2.305769e-04 0.01757162 9.821978e-01
#> 29  3.00 1.108001e-03     3 8.953873e-04 0.02392962 9.751750e-01
#> 30  3.25 5.072800e-04     3 2.878648e-03 0.03273349 9.643879e-01
#> 31  3.50 2.181784e-04     3 7.651606e-03 0.04316803 9.491804e-01
#> 32  3.75 8.815204e-05     3 1.706328e-02 0.05364464 9.292921e-01
#> 33  4.00 3.345874e-05     3 3.265689e-02 0.06233080 9.050123e-01
#> 
#> $cutscore
#> [1] -0.5  0.8
```

#### (2) P-method: individual ability estimates + SE from item metadata

When individual ability estimates are available, the P-method averages
conditional indices over the sample of examinees.

``` r

# P-method: SE computed internally from item metadata
cac_r_p <- cac_rud(
  x        = meta_cal_bin,          # SE computed from TIF
  cutscore = cutscore_th,
  theta    = score_bin$est.theta,   # individual ML estimates
  D        = 1.702
)

# Alternatively, supply ML-based standard errors directly:
cac_r_p2 <- cac_rud(
  cutscore = cutscore_th,
  theta    = score_bin$est.theta,
  se       = score_bin$se.theta     # individual SEs from ML scoring
)

# Confusion matrix
cac_r_p$confusion
#>     Expected
#> True         1         2         3
#>    1 0.2617201 0.0518740 0.0054059
#>    2 0.0504975 0.3460603 0.0534422
#>    3 0.0016839 0.0376158 0.1917003
cac_r_p2$confusion
#>     Expected
#> True         1         2         3
#>    1 0.2611347 0.0517079 0.0061574
#>    2 0.0504975 0.3460603 0.0534422
#>    3 0.0052723 0.0369660 0.1887617

# Marginal CA and CC indices
cac_r_p$marginal
#>     level  accuracy consistency
#>         1 0.2617201   0.2395853
#>         2 0.3460603   0.3015966
#>         3 0.1917003   0.1761846
#>  marginal 0.7994807   0.7173665
cac_r_p2$marginal
#>     level  accuracy consistency
#>         1 0.2611347   0.2396159
#>         2 0.3460603   0.3015966
#>         3 0.1887617   0.1745779
#>  marginal 0.7959568   0.7157905
```

Note that `cac_r_p` and `cac_r_p2` yield **nearly** identical results
because ML-based SEs are theoretically equivalent to
$`1/\sqrt{I(\hat\theta)}`$. **However, slight discrepancies may occur in
practice due to how extreme boundary values are handled.** When using
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
with ML estimation, the SEs for extreme ability estimates (e.g.,
artificially bounded at -5 or 5 for all-correct or all-incorrect
responses)are internally capped at an arbitrary large value (e.g.,
99.99999) to prevent computational errors. In contrast, computing the SE
internally from the test information function (TIF) using
[`cac_rud()`](https://hwangQ.github.io/irtQ/reference/cac_rud.md)
calculates the exact analytical value at those bounded thetas.
Therefore, if the sample includes examinees with extreme scores, minor
differences in the final marginal indices might be observed.

------------------------------------------------------------------------

### Example 4: Mixed-format test with `cac_rud()`

``` r

# Cut scores on the theta scale
cutscore_th_mix <- c(-0.5, 0.7)

# D-method
cac_r_mix_d <- cac_rud(
  x        = meta_cal_mix,
  cutscore = cutscore_th_mix,
  weights  = quad_weights,
  D        = 1.702
)
cac_r_mix_d
#> $confusion
#>     Expected
#> True         1         2         3
#>    1 0.2456919 0.0197309 0.0000177
#>    2 0.0686227 0.3715068 0.0289894
#>    3 0.0000060 0.0456285 0.2198060
#> 
#> $marginal
#>     level  accuracy consistency
#>         1 0.2456919   0.2327622
#>         2 0.3715068   0.3380082
#>         3 0.2198060   0.2059670
#>  marginal 0.8370048   0.7767374
#> 
#> $conditional
#>    theta      weights level  accuracy consistency
#> 1  -4.00 3.345874e-05     1 0.8566147   0.7441107
#> 2  -3.75 8.815204e-05     1 0.8932248   0.8037172
#> 3  -3.50 2.181784e-04     1 0.9266828   0.8619148
#> 4  -3.25 5.072800e-04     1 0.9545778   0.9127087
#> 5  -3.00 1.108001e-03     1 0.9751163   0.9513859
#> 6  -2.75 2.273471e-03     1 0.9880033   0.9762881
#> 7  -2.50 4.382230e-03     1 0.9947104   0.9894765
#> 8  -2.25 7.935194e-03     1 0.9975914   0.9951943
#> 9  -2.00 1.349822e-02     1 0.9985758   0.9971557
#> 10 -1.75 2.157009e-02     1 0.9985441   0.9970924
#> 11 -1.50 3.238054e-02     1 0.9969301   0.9938791
#> 12 -1.25 4.566389e-02     1 0.9881681   0.9766162
#> 13 -1.00 6.049482e-02     1 0.9445119   0.8951817
#> 14 -0.75 7.528702e-02     1 0.7934408   0.6722146
#> 15 -0.50 8.801945e-02     2 0.4999660   0.4999660
#> 16 -0.25 9.667045e-02     2 0.7980060   0.6773199
#> 17  0.00 9.973910e-02     2 0.9447570   0.8947715
#> 18  0.25 9.667045e-02     2 0.9308300   0.8705262
#> 19  0.50 8.801945e-02     2 0.7514652   0.6263028
#> 20  0.75 7.528702e-02     3 0.5668379   0.5089235
#> 21  1.00 6.049482e-02     3 0.8346404   0.7239683
#> 22  1.25 4.566389e-02     3 0.9513866   0.9074998
#> 23  1.50 3.238054e-02     3 0.9849528   0.9703585
#> 24  1.75 2.157009e-02     3 0.9935915   0.9872651
#> 25  2.00 1.349822e-02     3 0.9958126   0.9916603
#> 26  2.25 7.935194e-03     3 0.9959260   0.9918853
#> 27  2.50 4.382230e-03     3 0.9946788   0.9894142
#> 28  2.75 2.273471e-03     3 0.9917223   0.9835804
#> 29  3.00 1.108001e-03     3 0.9861388   0.9726510
#> 30  3.25 5.072800e-04     3 0.9767029   0.9544173
#> 31  3.50 2.181784e-04     3 0.9622598   0.9270109
#> 32  3.75 8.815204e-05     3 0.9421717   0.8897937
#> 33  4.00 3.345874e-05     3 0.9166076   0.8439349
#> 
#> $prob.level
#>    theta      weights level    p.level.1   p.level.2    p.level.3
#> 1  -4.00 3.345874e-05     1 8.566147e-01 0.067092130 7.629317e-02
#> 2  -3.75 8.815204e-05     1 8.932248e-01 0.062503699 4.427150e-02
#> 3  -3.50 2.181784e-04     1 9.266828e-01 0.052247016 2.107015e-02
#> 4  -3.25 5.072800e-04     1 9.545778e-01 0.037848456 7.573735e-03
#> 5  -3.00 1.108001e-03     1 9.751163e-01 0.023039815 1.843928e-03
#> 6  -2.75 2.273471e-03     1 9.880033e-01 0.011727704 2.689928e-04
#> 7  -2.50 4.382230e-03     1 9.947104e-01 0.005268059 2.156149e-05
#> 8  -2.25 7.935194e-03     1 9.975914e-01 0.002407630 1.006886e-06
#> 9  -2.00 1.349822e-02     1 9.985758e-01 0.001424132 3.924484e-08
#> 10 -1.75 2.157009e-02     1 9.985441e-01 0.001455902 2.694048e-09
#> 11 -1.50 3.238054e-02     1 9.969301e-01 0.003069877 8.274931e-10
#> 12 -1.25 4.566389e-02     1 9.881681e-01 0.011831890 2.019165e-09
#> 13 -1.00 6.049482e-02     1 9.445119e-01 0.055488066 2.996825e-08
#> 14 -0.75 7.528702e-02     1 7.934408e-01 0.206558164 1.033168e-06
#> 15 -0.50 8.801945e-02     2 5.000000e-01 0.499965964 3.403613e-05
#> 16 -0.25 9.667045e-02     2 2.012604e-01 0.798006021 7.336209e-04
#> 17  0.00 9.973910e-02     2 4.605801e-02 0.944757011 9.184975e-03
#> 18  0.25 9.667045e-02     2 5.520485e-03 0.930829966 6.364955e-02
#> 19  0.50 8.801945e-02     2 3.358926e-04 0.751465218 2.481989e-01
#> 20  0.75 7.528702e-02     3 1.286862e-05 0.433149194 5.668379e-01
#> 21  1.00 6.049482e-02     3 5.771300e-07 0.165358987 8.346404e-01
#> 22  1.25 4.566389e-02     3 6.569883e-08 0.048613285 9.513866e-01
#> 23  1.50 3.238054e-02     3 2.944861e-08 0.015047125 9.849528e-01
#> 24  1.75 2.157009e-02     3 4.825563e-08 0.006408450 9.935915e-01
#> 25  2.00 1.349822e-02     3 1.985586e-07 0.004187199 9.958126e-01
#> 26  2.25 7.935194e-03     3 1.337600e-06 0.004072613 9.959260e-01
#> 27  2.50 4.382230e-03     3 1.035725e-05 0.005310819 9.946788e-01
#> 28  2.75 2.273471e-03     3 7.257494e-05 0.008205172 9.917223e-01
#> 29  3.00 1.108001e-03     3 4.045822e-04 0.013456618 9.861388e-01
#> 30  3.25 5.072800e-04     3 1.714398e-03 0.021582693 9.767029e-01
#> 31  3.50 2.181784e-04     3 5.553166e-03 0.032186985 9.622598e-01
#> 32  3.75 8.815204e-05     3 1.418039e-02 0.043647933 9.421717e-01
#> 33  4.00 3.345874e-05     3 2.968957e-02 0.053702855 9.166076e-01
#> 
#> $cutscore
#> [1] -0.5  0.7

# P-method
cac_r_mix_p <- cac_rud(
  x        = meta_cal_mix,
  cutscore = cutscore_th_mix,
  theta    = score_mix$est.theta,
  D        = 1.702
)
cac_r_mix_p$marginal
#>     level  accuracy consistency
#>         1 0.2749796   0.2574348
#>         2 0.3500664   0.3156845
#>         3 0.2149411   0.2002216
#>  marginal 0.8399872   0.7733409
```

------------------------------------------------------------------------

## Comparing Lee’s and Rudner’s Methods

``` r

# Side-by-side comparison of marginal CA and CC (binary test, P-method)
cat("=== Lee's method (P-method) ===\n")
#> === Lee's method (P-method) ===
print(cac_l_p$marginal)
#>     level  accuracy consistency
#>         1 0.1639447   0.1507562
#>         2 0.3680726   0.3237276
#>         3 0.2522645   0.2380921
#>  marginal 0.7842818   0.7125758

cat("\n=== Rudner's method (P-method, SE from TIF) ===\n")
#> 
#> === Rudner's method (P-method, SE from TIF) ===
print(cac_r_p$marginal)
#>     level  accuracy consistency
#>         1 0.2617201   0.2395853
#>         2 0.3460603   0.3015966
#>         3 0.1917003   0.1761846
#>  marginal 0.7994807   0.7173665
```

The two methods share the same conceptual framework — both estimate CA
and CC by computing, for each ability level, the probabilities of being
assigned to each performance category — but differ in how they model the
conditional score distribution and what metric the cut scores operate
on:

| Aspect | Lee (2010) | Rudner (2001, 2005) |
|----|----|----|
| Cut-score metric | Observed summed score (or theta, converted via TCC) | Theta scale only |
| Conditional distribution | Exact conditional summed-score distribution via Lord–Wingersky recursion | Normal approximation: $`\hat\theta \mid \theta \sim N(\theta, \text{SE}^2)`$ |
| SE source | Implicit (via IRT-based score distribution) | Explicit: $`\text{SE}(\theta) = 1/\sqrt{I(\theta)}`$ |
| Typical CA/CC values | Generally similar to Rudner’s method when IRT fits well | Generally similar to Lee’s method when IRT fits well |

**Practical guidance:**

- Use
  **[`cac_lee()`](https://hwangQ.github.io/irtQ/reference/cac_lee.md)**
  when cut scores are defined on the observed summed-score scale (e.g.,
  raw scores such as 70 out of 100), or when the exact conditional score
  distribution is desired. This method is the more rigorous of the two
  and is applicable to mixed-format assessments with any combination of
  IRT models.
- Use
  **[`cac_rud()`](https://hwangQ.github.io/irtQ/reference/cac_rud.md)**
  when cut scores are expressed on the theta (ability) scale, or when
  standard errors from ability estimation are already available. This
  method is simpler to implement and produces results very similar to
  Lee’s method when the normality assumption for ability estimates is
  reasonable.

------------------------------------------------------------------------

## References

Kolen, Michael J., and Robert L. Brennan. 2004. *Test Equating, Scaling,
and Linking*. 2nd ed. Springer.

Lee, Won-Chan. 2010. “Classification Consistency and Accuracy for
Complex Assessments Using Item Response Theory.” *Journal of Educational
Measurement* 47 (1): 1–17.
<https://doi.org/10.1111/j.1745-3984.2009.00096.x>.

Lord, Frederic M., and Marilyn S. Wingersky. 1984. “Comparison of IRT
True-Score and Equipercentile Observed-Score Equatings.” *Applied
Psychological Measurement* 8 (4): 453–61.
<https://doi.org/10.1177/014662168400800409>.

Rudner, Lawrence M. 2001. “Computing the Expected Proportions of
Misclassified Examinees.” *Practical Assessment, Research & Evaluation*
7 (14): 1–5. <https://ericae.net/pare/getvn.asp?v=7&n=14>.

Rudner, Lawrence M. 2005. “Expected Classification Accuracy.” *Practical
Assessment, Research & Evaluation* 10 (13): 1–4.
<https://doi.org/10.7275/56a5-bs64>.
