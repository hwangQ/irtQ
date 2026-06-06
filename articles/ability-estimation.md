# Ability Estimation

## Overview

After item calibration,
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
estimates examinee latent ability ($`\theta`$) from item response data.
The **irtQ** package supports two broad categories of IRT scoring
method, which differ in the information they use from the response data.

### IRT Pattern-Based Scoring

Pattern-based scoring uses each examinee’s *full item-response pattern*
— the complete vector of correct/incorrect (or polytomous) responses
across all items — to estimate $`\theta`$. Because it conditions on the
individual responses rather than only their sum, pattern-based scoring
can distinguish examinees who have the same total score but answered
different items correctly. In principle, this makes better use of the
available measurement information (Kolen & Tong, 2010).

Under the assumption of local independence, the likelihood of an
observed response pattern $`\mathbf{U} = (u_1, u_2, \ldots, u_n)`$ given
ability $`\theta`$ is

``` math
L(\theta \mid \mathbf{U}) = \prod_{i=1}^{n} P(U_i = u_i \mid \theta),
```

where $`P(U_i = u_i \mid \theta)`$ is the category characteristic
function for item $`i`$. The different pattern-based methods use this
likelihood in different ways:

| Method | Key | Description |
|----|----|----|
| Maximum Likelihood (ML) | `"ML"` | Finds the $`\theta`$ maximising $`L(\theta \mid \mathbf{U})`$. Unbiased in large samples but undefined for all-correct or all-incorrect response patterns (Kolen & Tong, 2010). |
| ML with Fences (MLF) | `"MLF"` | Augments the likelihood with imaginary “fence” items (a lower fence with a fixed correct response and an upper fence with a fixed incorrect response) to resolve the boundary-score problem while avoiding the shrinkage of Bayesian methods (Han, 2016). |
| Weighted Likelihood (WL) | `"WL"` | Multiplies the likelihood by a weighting function derived from the square root of test information before maximizing, yielding estimates that are nearly unbiased to order $`O(n^{-1})`$ — substantially less biased than both ML and Bayesian modal estimation over the full $`\theta`$ scale (Warm, 1989). |
| Maximum A Posteriori (MAP) | `"MAP"` | Returns the mode of the posterior $`p(\theta \mid \mathbf{U}) \propto L(\theta \mid \mathbf{U})\,g(\theta)`$, where $`g(\theta)`$ is a normal prior. Handles boundary scores but shrinks estimates toward the prior mean. |
| Expected A Posteriori (EAP) | `"EAP"` | Returns the mean of the posterior distribution via Gaussian quadrature integration over $`\theta`$. Also handles boundary scores; has the smallest conditional error variance among pattern-based methods but introduces the most shrinkage (Kolen & Tong, 2010). |

A key practical limitation of plain ML is that it yields no finite
solution when all responses to the scored items are correct or all are
incorrect (perfect and zero scores). MLF, WL, MAP, and EAP all handle
these boundary patterns, each through a different mechanism (Han, 2016).

### IRT Summed-Score Scoring

Summed-score methods map each possible *raw total score*
$`X_s = \sum_{i=1}^{n} U_i`$ to a $`\theta`$ estimate rather than using
the full response pattern. Because all examinees with the same total
score receive the same $`\theta`$ estimate, these methods are
transparent and easy for test users to understand (Kolen & Tong, 2010).
In particular, Kolen and Tong (2010) note that the statistical
difference between summed-score and pattern-score estimators is
typically smaller in practice than the difference between Bayesian and
non-Bayesian estimators — so the choice between summed-score and
pattern-score scoring need not be driven primarily by accuracy concerns.

An important practical consideration when using these methods is the
treatment of missing data. While pattern-based scoring methods handle
missing responses (`NA`) flexibly by conditioning only on the observed
items, summed-score methods (`"EAP.SUM"` and `"INV.TCC"`) automatically
treat any missing responses as incorrect (recoded as 0s) to calculate
the total raw score.

| Method | Key | Description |
|----|----|----|
| EAP for Summed Scores | `"EAP.SUM"` | Computes the Bayesian EAP estimate $`\hat{\theta}_{sEAP} = E(\theta \mid X_s)`$ for each possible summed-score value using the Lord–Wingersky recursive algorithm (Thissen et al., 1995; Thissen & Orlando, 2001). Returns a score table mapping every feasible raw score to a $`\theta`$ estimate and its SE. |
| Inverse TCC | `"INV.TCC"` | Solves $`\hat{\theta}_{TCF}`$ from the test characteristic function (TCF) equation $`\tau_s(\theta) = X_s`$ numerically. A non-Bayesian estimator that is monotonically related to $`X_s`$ and does not depend on the prior distribution (Kolen & Tong, 2010; Stocking, 1996). Standard errors are computed using a recursion-based analytical approach (**lim2020?**). Linear interpolation (`intpol = TRUE`) handles scores outside the range where the TCC is invertible (e.g., scores below the sum of guessing parameters in 3PLM). |

``` r

library(irtQ)
set.seed(2026)
```

------------------------------------------------------------------------

## Key `est_score()` Arguments

Before working through the examples, it is helpful to understand the
most important arguments of
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md):

| Argument | Description |
|----|----|
| `x` | An `est_irt` object (from [`irtQ::est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)), or a data frame of item metadata. When an `est_irt` object is provided, the response data embedded in it are used automatically. |
| `data` | Response matrix (rows = examinees, columns = items). **Required** when `x` is a metadata data frame rather than an `est_irt` object. |
| `D` | IRT scaling constant. Use `1` for the logistic metric or `1.702` to approximate the normal-ogive metric. |
| `method` | Scoring method: `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, `"EAP"`, `"EAP.SUM"`, or `"INV.TCC"`. |
| `range` | `c(lower, upper)` bounding the ability scale for iterative methods (`"ML"`, `"MLF"`, `"WL"`, `"MAP"`). Default `c(-5, 5)`. |
| `norm.prior` | `c(mean, sd)` of the normal prior distribution used by `"MAP"`, `"EAP"`, and `"EAP.SUM"`. Default `c(0, 1)`. |
| `nquad` | Number of Gaussian quadrature points for numerical integration (used by `"EAP"` and `"EAP.SUM"`). Default `41`. |
| `fence.a` | Discrimination parameter of the virtual fence items added in `"MLF"`. Default `3`. |
| `fence.b` | Location of the fence items on the $`\theta`$ scale in `"MLF"`. Defaults to the `range` bounds when `NULL`. |
| `tol` | Convergence tolerance for iterative methods. Default `1e-4`. |
| `max.iter` | Maximum Newton–Raphson iterations. Default `100`. |
| `se` | Logical; compute standard errors? Always returned for `"EAP.SUM"` and `"INV.TCC"`. Default `TRUE`. |
| `intpol` | Logical; apply linear interpolation in `"INV.TCC"` for extreme scores? Default `TRUE`. |
| `range.tcc` | Ability range used for the `"INV.TCC"` interpolation grid. Default `c(-7, 7)`. |

**Return values** differ by method:

- `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, `"EAP"`: a two-column data frame
  with columns `est.theta` and `se.theta` (one row per examinee).
- `"EAP.SUM"`, `"INV.TCC"`: a list with two elements — `$est.par`
  (examinee-level estimates including observed sum scores) and
  `$score.table` (the complete raw-score-to-$`\theta`$ mapping table).

------------------------------------------------------------------------

## Setup: Shared Data for All Examples

### Mixed-Format Test (15 Dichotomous + 10 Polytomous Items)

We define a 25-item test (15 dichotomous 3PLM items + 10 polytomous GRM
items with 4 categories) and simulate responses from 800 examinees. This
dataset is used in the **Part 1** and **Part 2** examples for the
mixed-format test.

``` r

# --- Item metadata: 15 3PLM + 10 GRM items ---
meta_mixed <- shape_df(
  par.drm = list(
    a = c(1.0, 1.2, 0.8, 1.4, 0.9, 1.1, 1.3, 0.7, 1.0, 1.2,
          0.9, 1.1, 1.4, 0.85, 1.0),
    b = c(-2.0, -1.5, -1.0, -0.6, -0.2, 0.0, 0.4, 0.8, 1.1, 1.5,
          -1.3, -0.4,  0.5, 1.0,  1.8),
    g = rep(0.15, 15)
  ),
  par.prm = list(
    a = c(1.5, 1.2, 1.0, 1.3, 0.9, 1.1, 0.8, 1.4, 1.2, 1.0),
    d = list(
      c(-1.5, -0.3,  0.9),
      c(-1.2,  0.0,  1.1),
      c(-0.9,  0.4,  1.4),
      c(-1.1, -0.1,  1.0),
      c(-1.3,  0.3,  1.2),
      c(-0.8,  0.5,  1.5),
      c(-1.0,  0.1,  0.9),
      c(-1.4, -0.2,  1.1),
      c(-0.7,  0.6,  1.3),
      c(-1.2,  0.0,  1.0)
    )
  ),
  item.id = c(paste0("BI", 1:15), paste0("PI", 1:10)),
  cats    = c(rep(2, 15), rep(4, 10)),
  model   = c(rep("3PLM", 15), rep("GRM", 10))
)

# --- Simulate 800 examinees from N(0, 1) ---
theta_mixed <- rnorm(800, mean = 0, sd = 1)
resp_mixed  <- simdat(x = meta_mixed, theta = theta_mixed, D = 1.702)

dim(resp_mixed)   # 800 examinees × 25 items
#> [1] 800  25
```

### Dichotomous-Only Test (30 Items, 3PLM)

For examples where a simpler, purely dichotomous test is informative —
particularly for summed-score scoring methods where the score table is
easiest to interpret — we additionally prepare a 30-item test with 3PLM
items.

``` r

# --- Item metadata: 30 3PLM items ---
meta_dich <- shape_df(
  par.drm = list(
    a = c(1.0, 1.2, 0.8, 1.4, 0.9, 1.1, 1.3, 0.7, 1.0, 1.2,
          0.9, 1.1, 1.4, 0.85, 1.0, 1.2, 0.8, 1.3, 1.0, 0.9,
          1.1, 1.3, 0.8, 1.0, 1.2, 0.9, 1.4, 1.1, 0.7, 1.0),
    b = c(-2.0, -1.5, -1.0, -0.6, -0.2,  0.0,  0.4,  0.8,  1.1,  1.5,
          -1.3, -0.4,  0.5,  1.0,  1.8, -0.8,  0.2,  0.7, -1.1,  1.3,
          -1.8, -0.9,  0.3,  1.2, -0.5,  0.6, -1.4,  0.1,  0.9, -0.3),
    g = rep(0.15, 30)
  ),
  cats  = rep(2, 30),
  model = rep("3PLM", 30)
)

# --- Simulate 500 examinees from N(0, 1) ---
theta_dich <- rnorm(500, mean = 0, sd = 1)
resp_dich  <- simdat(x = meta_dich, theta = theta_dich, D = 1.702)

dim(resp_dich)   # 500 examinees × 30 items
#> [1] 500  30
```

------------------------------------------------------------------------

## Part 1: IRT Pattern-Based Scoring

This section demonstrates all five pattern-based methods. For each
method, examples are shown for both the **mixed-format test** and the
**dichotomous-only test**.

### ML — Maximum Likelihood

ML estimation finds the $`\theta`$ that maximizes the log-likelihood of
the observed response pattern. It does not use a prior distribution,
making it purely data-driven. A `range` bound is required to prevent
$`|\hat{\theta}| \to \infty`$ for all-correct or all-incorrect response
patterns (Kolen & Tong, 2010).

``` r

# --- Mixed-format test (3PLM + GRM) ---
score_ml_mixed <- est_score(
  x        = meta_mixed,
  data     = resp_mixed,
  D        = 1.702,
  method   = "ML",
  range    = c(-5, 5),
  tol      = 0.0001,
  max.iter = 100,
  se       = TRUE
)

head(score_ml_mixed)
#>    est.theta  se.theta
#> 1  0.3613411 0.2521047
#> 2 -2.0907154 0.4083068
#> 3  0.4091581 0.2519701
#> 4 -0.2368663 0.2544738
#> 5 -0.3956027 0.2566540
#> 6 -4.3933779 3.3555512
```

``` r

# --- Dichotomous-only test (30 3PLM items) ---
score_ml_dich <- est_score(
  x        = meta_dich,
  data     = resp_dich,
  D        = 1.702,
  method   = "ML",
  range    = c(-5, 5),
  tol      = 0.0001,
  max.iter = 100,
  se       = TRUE
)

head(score_ml_dich)
#>    est.theta  se.theta
#> 1  0.8622801 0.3349916
#> 2 -0.2471560 0.3222841
#> 3 -1.4031575 0.4002763
#> 4  1.0747667 0.3502389
#> 5 -0.9934755 0.3510865
#> 6  0.1899347 0.3190865
```

### MLF — ML with Fences (Han, 2016)

MLF adds imaginary “fence” items with fixed responses at both ends of
the $`\theta`$ scale. This makes the log-likelihood unimodal,
eliminating the boundary-score problem of plain ML while producing
estimates that are not shrunk toward a prior mean — unlike MAP or EAP
(Han, 2016).

``` r

# --- Mixed-format test ---
score_mlf_mixed <- est_score(
  x        = meta_mixed,
  data     = resp_mixed,
  D        = 1.702,
  method   = "MLF",
  range    = c(-5, 5),
  fence.a  = 3.0,    # discrimination of fence items
  fence.b  = NULL,   # fence locations default to range bounds
  se       = TRUE
)

head(score_mlf_mixed)
#>    est.theta  se.theta
#> 1  0.3613411 0.2521047
#> 2 -2.0907150 0.4083063
#> 3  0.4091581 0.2519701
#> 4 -0.2368663 0.2544738
#> 5 -0.3956027 0.2566540
#> 6 -3.9498209 1.7591919
```

``` r

# --- Dichotomous-only test ---
score_mlf_dich <- est_score(
  x        = meta_dich,
  data     = resp_dich,
  D        = 1.702,
  method   = "MLF",
  range    = c(-5, 5),
  fence.a  = 3.0,
  fence.b  = NULL,
  se       = TRUE
)

head(score_mlf_dich)
#>    est.theta  se.theta
#> 1  0.8622801 0.3349916
#> 2 -0.2471560 0.3222841
#> 3 -1.4031574 0.4002763
#> 4  1.0747666 0.3502389
#> 5 -0.9934755 0.3510865
#> 6  0.1899347 0.3190865
```

### WL — Weighted Likelihood (Warm, 1989)

WL multiplies the likelihood by a weighting function $`w(\theta)`$
derived from the square root of test information before maximizing. The
resulting estimates are nearly unbiased to order $`O(n^{-1})`$ —
substantially less biased than both plain ML ($`O(n^{-1})`$ bias with a
positive correlation with $`\theta`$) and Bayesian estimators (also
$`O(n^{-1})`$ but with negative correlation) — across the entire
$`\theta`$ scale, making WL generally preferable to plain ML when an
unbiased non-Bayesian estimate is desired (Warm, 1989).

``` r

# --- Mixed-format test ---
score_wl_mixed <- est_score(
  x        = meta_mixed,
  data     = resp_mixed,
  D        = 1.702,
  method   = "WL",
  range    = c(-5, 5),
  se       = TRUE
)

head(score_wl_mixed)
#>    est.theta  se.theta
#> 1  0.3590491 0.2521112
#> 2 -1.9786439 0.3788110
#> 3  0.4072527 0.2519754
#> 4 -0.2374576 0.2544800
#> 5 -0.3945439 0.2566360
#> 6 -2.7427809 0.6990180
```

``` r

# --- Dichotomous-only test ---
score_wl_dich <- est_score(
  x        = meta_dich,
  data     = resp_dich,
  D        = 1.702,
  method   = "WL",
  range    = c(-5, 5),
  se       = TRUE
)

head(score_wl_dich)
#>    est.theta  se.theta
#> 1  0.8338661 0.3333307
#> 2 -0.2585616 0.3224005
#> 3 -1.3786223 0.3963382
#> 4  1.0353358 0.3470332
#> 5 -0.9875636 0.3505997
#> 6  0.1768412 0.3191358
```

### MAP — Maximum A Posteriori

MAP incorporates a normal prior $`g(\theta)`$ and returns the **mode**
of the posterior distribution. Compared with EAP, MAP shrinks estimates
less strongly toward the prior mean, and it can still produce estimates
outside the bulk of the prior when the data are informative enough.

``` r

# --- Mixed-format test ---
score_map_mixed <- est_score(
  x          = meta_mixed,
  data       = resp_mixed,
  D          = 1.702,
  method     = "MAP",
  range      = c(-5, 5),
  norm.prior = c(0, 1),
  se         = TRUE
)

head(score_map_mixed)
#>    est.theta  se.theta
#> 1  0.3376385 0.2445168
#> 2 -1.7942871 0.3222272
#> 3  0.3858189 0.2443928
#> 4 -0.2232262 0.2464886
#> 5 -0.3723269 0.2482450
#> 6 -2.2728080 0.4234079
```

``` r

# --- Dichotomous-only test ---
score_map_dich <- est_score(
  x          = meta_dich,
  data       = resp_dich,
  D          = 1.702,
  method     = "MAP",
  range      = c(-5, 5),
  norm.prior = c(0, 1),
  se         = TRUE
)

head(score_map_dich)
#>    est.theta  se.theta
#> 1  0.7784438 0.3136958
#> 2 -0.2250451 0.3065592
#> 3 -1.2410181 0.3526268
#> 4  0.9546880 0.3227622
#> 5 -0.8904420 0.3246810
#> 6  0.1717731 0.3040467
```

### EAP — Expected A Posteriori

EAP returns the **mean** of the posterior distribution, integrating over
a Gaussian quadrature grid. It has the smallest conditional error
variance among pattern-based estimators, but it introduces the most
shrinkage toward the prior mean and is sensitive to the choice of prior
distribution, particularly for short or less reliable tests (Kolen &
Tong, 2010).

``` r

# --- Mixed-format test ---
score_eap_mixed <- est_score(
  x          = meta_mixed,
  data       = resp_mixed,
  D          = 1.702,
  method     = "EAP",
  norm.prior = c(0, 1),
  nquad      = 41,
  se         = TRUE
)

head(score_eap_mixed)
#>    est.theta  se.theta
#> 1  0.3401543 0.2614950
#> 2 -1.8686950 0.3760946
#> 3  0.3941764 0.2350285
#> 4 -0.2249379 0.2601181
#> 5 -0.3916921 0.2436289
#> 6 -2.3955969 0.4613010
```

``` r

# --- Dichotomous-only test ---
score_eap_dich <- est_score(
  x          = meta_dich,
  data       = resp_dich,
  D          = 1.702,
  method     = "EAP",
  norm.prior = c(0, 1),
  nquad      = 41,
  se         = TRUE
)

head(score_eap_dich)
#>    est.theta  se.theta
#> 1  0.7804815 0.3207942
#> 2 -0.2439516 0.3084545
#> 3 -1.2990216 0.3562331
#> 4  0.9665304 0.3376843
#> 5 -0.9267918 0.3328926
#> 6  0.1588514 0.3154466
```

### Comparison: ML vs WL vs MAP vs EAP

The following code computes correlations with the true $`\theta`$ values
and mean absolute errors for both datasets, summarising the practical
differences among the four pattern-based methods.

``` r

# --- Mixed-format test: recovery statistics ---
cors_mixed <- c(
  ML  = cor(score_ml_mixed$est.theta,  theta_mixed),
  WL  = cor(score_wl_mixed$est.theta,  theta_mixed),
  MAP = cor(score_map_mixed$est.theta, theta_mixed),
  EAP = cor(score_eap_mixed$est.theta, theta_mixed)
)
round(cors_mixed, 4)
#>     ML     WL    MAP    EAP 
#> 0.9504 0.9640 0.9645 0.9645

mae_mixed <- c(
  ML  = mean(abs(score_ml_mixed$est.theta  - theta_mixed)),
  WL  = mean(abs(score_wl_mixed$est.theta  - theta_mixed)),
  MAP = mean(abs(score_map_mixed$est.theta - theta_mixed)),
  EAP = mean(abs(score_eap_mixed$est.theta - theta_mixed))
)
round(mae_mixed, 4)
#>     ML     WL    MAP    EAP 
#> 0.2318 0.2121 0.2083 0.2074
```

``` r

# --- Dichotomous-only test: recovery statistics ---
cors_dich <- c(
  ML  = cor(score_ml_dich$est.theta,  theta_dich),
  WL  = cor(score_wl_dich$est.theta,  theta_dich),
  MAP = cor(score_map_dich$est.theta, theta_dich),
  EAP = cor(score_eap_dich$est.theta, theta_dich)
)
round(cors_dich, 4)
#>     ML     WL    MAP    EAP 
#> 0.9108 0.9475 0.9503 0.9501

mae_dich <- c(
  ML  = mean(abs(score_ml_dich$est.theta  - theta_dich)),
  WL  = mean(abs(score_wl_dich$est.theta  - theta_dich)),
  MAP = mean(abs(score_map_dich$est.theta - theta_dich)),
  EAP = mean(abs(score_eap_dich$est.theta - theta_dich))
)
round(mae_dich, 4)
#>     ML     WL    MAP    EAP 
#> 0.3346 0.2703 0.2463 0.2478
```

### Providing an `est_irt` Object vs. Separate Metadata

[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
accepts either an `est_irt` object (which embeds the response data and
item parameters) or item metadata and response data provided separately.
When an `est_irt` object is passed to the `x` argument,
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
automatically extracts both the item parameters and the embedded
response data. This streamlines the post-calibration scoring workflow by
eliminating the need to explicitly supply a separate `data` matrix.

Both workflows produce identical results:

``` r

# --- Calibrate the mixed-format test first ---
mod_cal <- est_irt(
  data       = resp_mixed,
  D          = 1.702,
  model      = c(rep("3PLM", 15), rep("GRM", 10)),
  cats       = c(rep(2, 15), rep(4, 10)),
  use.gprior = TRUE,
  gprior     = list(dist = "beta", params = c(4, 16)),
  EmpHist    = FALSE,
  Etol       = 0.01,
  MaxE       = 150,
  se         = FALSE,
  verbose    = FALSE
)

# Score directly from the est_irt object (response data are automatically embedded)
score_from_obj <- est_score(
  x          = mod_cal,     # est_irt object
  method     = "EAP",
  norm.prior = c(0, 1),
  nquad      = 41
)

head(score_from_obj)
#>    est.theta  se.theta
#> 1  0.2817412 0.2637351
#> 2 -1.9220269 0.3786857
#> 3  0.3497562 0.2414748
#> 4 -0.2758520 0.2581984
#> 5 -0.4236150 0.2376667
#> 6 -2.4311991 0.4608047
```

------------------------------------------------------------------------

## Part 2: IRT Summed-Score Scoring

Summed-score methods assign the same $`\theta`$ estimate to all
examinees with the same total raw score. They are computationally
efficient, straightforward to communicate to test users, and — per Kolen
and Tong (2010) — typically give results that are statistically
comparable to pattern-based methods when the key choice is between
Bayesian and non-Bayesian estimation rather than between summed-score
and pattern-score approaches.

### EAP.SUM — EAP Based on Summed Scores

`"EAP.SUM"` computes the Bayesian EAP estimate $`\hat{\theta}_{sEAP}`$
for each possible summed-score value using the Lord–Wingersky recursive
algorithm (Thissen et al., 1995; Thissen & Orlando, 2001), then maps
each examinee’s observed sum score to the corresponding table entry.

#### Mixed-format test

``` r

score_eapsum_mixed <- est_score(
  x          = meta_mixed,
  data       = resp_mixed,
  D          = 1.702,
  method     = "EAP.SUM",
  norm.prior = c(0, 1),
  nquad      = 41
)

# Individual-level estimates (observed sum score + ability estimate)
head(score_eapsum_mixed$est.par)
#>   sum.score  est.theta  se.theta
#> 1        27  0.3392105 0.2636885
#> 2         6 -1.6988246 0.3685688
#> 3        28  0.4205068 0.2517488
#> 4        20 -0.2482483 0.2752892
#> 5        17 -0.5106020 0.2607398
#> 6         4 -2.0194923 0.4254111

# Score table: every feasible raw score mapped to a theta estimate
score_eapsum_mixed$score.table
#>    sum.score     est.theta  se.theta
#> 1          0 -2.6400959233 0.4729966
#> 2          1 -2.5026131153 0.4716512
#> 3          2 -2.3508016353 0.4632335
#> 4          3 -2.1879143687 0.4474294
#> 5          4 -2.0194923464 0.4254111
#> 6          5 -1.8536527623 0.3981374
#> 7          6 -1.6988246243 0.3685688
#> 8          7 -1.5578941210 0.3439696
#> 9          8 -1.4269102187 0.3289357
#> 10         9 -1.3021645129 0.3182185
#> 11        10 -1.1853933902 0.3045486
#> 12        11 -1.0791619459 0.2900754
#> 13        12 -0.9804716266 0.2830070
#> 14        13 -0.8826788980 0.2847774
#> 15        14 -0.7828603080 0.2865220
#> 16        15 -0.6849568417 0.2801191
#> 17        16 -0.5943130061 0.2681263
#> 18        17 -0.5106020161 0.2607398
#> 19        18 -0.4277526168 0.2640573
#> 20        19 -0.3400918980 0.2726108
#> 21        20 -0.2482482596 0.2752892
#> 22        21 -0.1586251209 0.2673844
#> 23        22 -0.0763635644 0.2553535
#> 24        23  0.0001459706 0.2503508
#> 25        24  0.0776097390 0.2568576
#> 26        25  0.1615809247 0.2677019
#> 27        26  0.2510607580 0.2714061
#> 28        27  0.3392104928 0.2636885
#> 29        28  0.4205067543 0.2517488
#> 30        29  0.4965738725 0.2475675
#> 31        30  0.5746977613 0.2558327
#> 32        31  0.6611728067 0.2685620
#> 33        32  0.7550654130 0.2733434
#> 34        33  0.8488276691 0.2664501
#> 35        34  0.9369105589 0.2573931
#> 36        35  1.0225444470 0.2591007
#> 37        36  1.1147070351 0.2728124
#> 38        37  1.2188917558 0.2866253
#> 39        38  1.3314507219 0.2909964
#> 40        39  1.4460344680 0.2917939
#> 41        40  1.5650918452 0.3022934
#> 42        41  1.6991744434 0.3233859
#> 43        42  1.8555106609 0.3455258
#> 44        43  2.0364523442 0.3686564
#> 45        44  2.2603373203 0.4061526
#> 46        45  2.5870193620 0.4804288
```

#### Dichotomous-only test

The score table is easiest to interpret for a purely dichotomous test,
because the possible raw scores are simple number-correct values from 0
to $`n`$.

``` r

score_eapsum_dich <- est_score(
  x          = meta_dich,
  data       = resp_dich,
  D          = 1.702,
  method     = "EAP.SUM",
  norm.prior = c(0, 1),
  nquad      = 41
)

head(score_eapsum_dich$est.par)
#>   sum.score  est.theta  se.theta
#> 1        22  0.6017808 0.3288213
#> 2        15 -0.3808020 0.3409553
#> 3         8 -1.4302687 0.4415852
#> 4        24  0.9185685 0.3393217
#> 5        11 -0.9554443 0.3832546
#> 6        19  0.1698180 0.3273205
score_eapsum_dich$score.table
#>    sum.score   est.theta  se.theta
#> 1          0 -2.56690400 0.4857530
#> 2          1 -2.46166527 0.4918490
#> 3          2 -2.34475143 0.4964471
#> 4          3 -2.21538122 0.4986044
#> 5          4 -2.07377453 0.4968411
#> 6          5 -1.92153943 0.4899382
#> 7          6 -1.76129481 0.4779021
#> 8          7 -1.59623330 0.4615404
#> 9          8 -1.43026866 0.4415852
#> 10         9 -1.26733209 0.4198054
#> 11        10 -1.10923160 0.3997283
#> 12        11 -0.95544429 0.3832546
#> 13        12 -0.80585928 0.3687075
#> 14        13 -0.66125667 0.3555731
#> 15        14 -0.52045066 0.3464583
#> 16        15 -0.38080205 0.3409553
#> 17        16 -0.24214263 0.3349628
#> 18        17 -0.10569416 0.3288854
#> 19        18  0.03054969 0.3269117
#> 20        19  0.16981801 0.3273205
#> 21        20  0.31147110 0.3255322
#> 22        21  0.45434422 0.3246167
#> 23        22  0.60178080 0.3288213
#> 24        23  0.75674342 0.3342290
#> 25        24  0.91856851 0.3393217
#> 26        25  1.08957949 0.3492458
#> 27        26  1.27581010 0.3642885
#> 28        27  1.48289586 0.3840249
#> 29        28  1.72057723 0.4133105
#> 30        29  2.00660976 0.4569571
#> 31        30  2.37111990 0.5231270
```

### INV.TCC — Inverse Test Characteristic Curve

`"INV.TCC"` solves numerically for the $`\hat{\theta}_{TCF}`$ that
satisfies $`\tau_s(\hat{\theta}_{TCF}) = X_s`$, where $`\tau_s(\theta)`$
is the test characteristic function (Kolen & Tong, 2010; Stocking,
1996). Because this is a non-Bayesian estimator, it does not depend on a
prior distribution and the resulting scores are monotonically related to
raw scores. Linear interpolation (`intpol = TRUE`) handles raw scores
that fall outside the invertible range of the TCC (for example, when
3PLM items have nonzero guessing parameters, a raw score equal to the
sum of all guessing parameters $`g`$ cannot be mapped without
interpolation).

#### Mixed-format test

``` r

score_invtcc_mixed <- est_score(
  x         = meta_mixed,
  data      = resp_mixed,
  D         = 1.702,
  method    = "INV.TCC",
  intpol    = TRUE,
  range.tcc = c(-7, 5)
)

head(score_invtcc_mixed$est.par)
#>   sum.score  est.theta  se.theta
#> 1        27  0.3585547 0.2687055
#> 2         6 -1.8713672 0.7200154
#> 3        28  0.4470703 0.2688782
#> 4        20 -0.2635547 0.2752142
#> 5        17 -0.5412891 0.2831045
#> 6         4 -2.3921484 1.1814220
score_invtcc_mixed$score.table
#>    sum.score    est.theta  se.theta
#> 1          0 -7.000000000 1.4806955
#> 2          1 -5.644335938 1.4806554
#> 3          2 -4.288671875 1.4791257
#> 4          3 -2.933007813 1.4044008
#> 5          4 -2.392148438 1.1814220
#> 6          5 -2.089960937 0.9333995
#> 7          6 -1.871367188 0.7200154
#> 8          7 -1.694492188 0.5610683
#> 9          8 -1.542148437 0.4551439
#> 10         9 -1.405429687 0.3902050
#> 11        10 -1.279414063 0.3519944
#> 12        11 -1.161132812 0.3291356
#> 13        12 -1.048632813 0.3145724
#> 14        13 -0.940742188 0.3045528
#> 15        14 -0.836601562 0.2971722
#> 16        15 -0.735664062 0.2914564
#> 17        16 -0.637304688 0.2868665
#> 18        17 -0.541289062 0.2831045
#> 19        18 -0.447226562 0.2799822
#> 20        19 -0.354726562 0.2773786
#> 21        20 -0.263554687 0.2752142
#> 22        21 -0.173398437 0.2734275
#> 23        22 -0.083945312 0.2719697
#> 24        23  0.004960938 0.2708039
#> 25        24  0.093554688 0.2699015
#> 26        25  0.181914063 0.2692478
#> 27        26  0.270195312 0.2688427
#> 28        27  0.358554687 0.2687055
#> 29        28  0.447070312 0.2688782
#> 30        29  0.535976562 0.2694295
#> 31        30  0.625507812 0.2704598
#> 32        31  0.715820312 0.2721045
#> 33        32  0.807460937 0.2745521
#> 34        33  0.900820313 0.2780529
#> 35        34  0.996601563 0.2829665
#> 36        35  1.095585938 0.2898366
#> 37        36  1.198945312 0.2996019
#> 38        37  1.308164062 0.3140597
#> 39        38  1.425273438 0.3369204
#> 40        39  1.553398437 0.3758617
#> 41        40  1.697070312 0.4450146
#> 42        41  1.864023438 0.5654411
#> 43        42  2.068398437 0.7578475
#> 44        43  2.342382813 1.0188835
#> 45        44  2.789179687 1.2306823
#> 46        45  5.000000000 0.3668820
```

#### Dichotomous-only test

``` r

score_invtcc_dich <- est_score(
  x         = meta_dich,
  data      = resp_dich,
  D         = 1.702,
  method    = "INV.TCC",
  intpol    = TRUE,
  range.tcc = c(-7, 5)
)

head(score_invtcc_dich$est.par)
#>   sum.score  est.theta  se.theta
#> 1        22  0.6709766 0.3527933
#> 2        15 -0.4023828 0.3625169
#> 3         8 -1.6483203 0.8035430
#> 4        24  1.0197266 0.3864184
#> 5        11 -1.0432422 0.4651407
#> 6        19  0.2016797 0.3426603
score_invtcc_dich$score.table
#>    sum.score   est.theta  se.theta
#> 1          0 -7.00000000 1.3094176
#> 2          1 -6.20410156 1.3093310
#> 3          2 -5.40820312 1.3090261
#> 4          3 -4.61230469 1.3078806
#> 5          4 -3.81640625 1.3029720
#> 6          5 -3.02050781 1.2766072
#> 7          6 -2.29722656 1.1497209
#> 8          7 -1.92207031 0.9790166
#> 9          8 -1.64832031 0.8035430
#> 10         9 -1.42253906 0.6520852
#> 11        10 -1.22425781 0.5393596
#> 12        11 -1.04324219 0.4651407
#> 13        12 -0.87378906 0.4198091
#> 14        13 -0.71214844 0.3922267
#> 15        14 -0.55566406 0.3745530
#> 16        15 -0.40238281 0.3625169
#> 17        16 -0.25089844 0.3540413
#> 18        17 -0.10019531 0.3481412
#> 19        18  0.05042969 0.3443825
#> 20        19  0.20167969 0.3426603
#> 21        20  0.35449219 0.3431309
#> 22        21  0.51019531 0.3462226
#> 23        22  0.67097656 0.3527933
#> 24        23  0.83957031 0.3646740
#> 25        24  1.01972656 0.3864184
#> 26        25  1.21660156 0.4299260
#> 27        26  1.43824219 0.5223865
#> 28        27  1.69941406 0.7067602
#> 29        28  2.03347656 1.0082549
#> 30        29  2.55042969 1.3141757
#> 31        30  5.00000000 0.4099288
```

------------------------------------------------------------------------

## References

Han, K. T. (2016). Maximum likelihood score estimation method with
fences for short-length tests and computerized adaptive tests. *Applied
Psychological Measurement*, *40*(4), 289–301.
<https://doi.org/10.1177/0146621616631317>

Kolen, M. J., & Tong, Y. (2010). Psychometric properties of IRT
proficiency estimates. *Educational Measurement: Issues and Practice*,
*29*(3), 8–14. <https://doi.org/10.1111/j.1745-3992.2010.00185.x>

Stocking, M. L. (1996). An alternative method for scoring adaptive
tests. *Journal of Educational and Behavioral Statistics*, *21*(4),
365–389. <https://doi.org/10.3102/10769986021004365>

Thissen, D., & Orlando, M. (2001). Item response theory for items scored
in two categories. In D. Thissen & H. Wainer (Eds.), *Test scoring* (pp.
73–140). Lawrence Erlbaum.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995).
Item response theory for scores on tests including polytomous items with
ordered responses. *Applied Psychological Measurement*, *19*(1), 39–49.
<https://doi.org/10.1177/014662169501900105>

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika*, *54*(3), 427–450.
<https://doi.org/10.1007/BF02294627>
