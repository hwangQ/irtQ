# DIF Detection

## Overview

Differential item functioning (DIF) occurs when examinees from different
groups who share the same underlying ability have systematically
different probabilities of responding correctly to an item. Detecting
DIF is an essential step in ensuring the fairness and validity of a
test.

**irtQ** provides three functions for DIF detection:

| Function | Method | Item Types | Groups |
|----|----|----|----|
| [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) | Residual-based DIF (RDIF) | Dichotomous (+ polytomous$`^*`$) | 2 |
| [`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md) | Generalized RDIF (GRDIF) | Dichotomous (+ polytomous$`^*`$) | ≥ 2 |
| [`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md) | CATSIB (modified SIBTEST for CAT) | Dichotomous only | 2 |

$`^*`$[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) and
[`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md) can accept
polytomous item response data, but their statistical performance for
polytomous items has not yet been formally validated. Research extending
and evaluating the RDIF and GRDIF frameworks to polytomous items is
currently ongoing.

All three functions require pooled item response data and group
membership labels. Pooled item parameter estimates and ability estimates
are also typically needed; however, for
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md), item
parameters (`x`) are not strictly required as long as ability estimates
and their standard errors are supplied externally. The `x` argument in
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md) is only
necessary when ability estimates need to be computed internally (i.e.,
`score = NULL`) or when purification is applied, since ability
re-estimation at each purification iteration requires item parameters.

``` r

library(irtQ)
set.seed(2026)
```

------------------------------------------------------------------------

## Setup: Simulating DIF Data

We simulate a 40-item dichotomous test (3PLM) with 1,000 examinees per
group. DIF is introduced into six items:

- **Items 1-2**: Uniform DIF, in which the focal group has a higher
  difficulty parameter ($`b + 0.5`$) while the discrimination parameter
  is unchanged.
- **Items 3-4**: Nonuniform DIF, in which the focal group has a lower
  discrimination parameter ($`a \times 0.5`$) while the difficulty
  parameter is unchanged.
- **Items 5-6**: Mixed DIF, in which the focal group has both a lower
  discrimination parameter ($`a \times 0.5`$) and a higher difficulty
  parameter ($`b + 0.5`$), making these items both less discriminating
  and harder for the focal group.
- **Items 7-40**: No DIF.

The focal group has a slightly lower mean ability ($`\mu = -0.3`$) than
the reference group ($`\mu = 0`$), representing a realistic impact
condition. The item difficulties are drawn from $`U(-2.5, 2.5)`$. Items
that are extremely easy or hard for the examinees can make the pooled
calibration unstable, so the range is kept moderate.

``` r

J   <- 40     # number of items
N_R <- 1000   # reference group size
N_F <- 1000   # focal group size

# Base item parameters (no DIF) shared by both groups
set.seed(123)
a_base <- round(runif(J, 0.7, 1.7), 2)
b_base <- round(runif(J, -2.5, 2.5), 2)
g_base <- rep(0.15, J)

# Reference group item metadata (no DIF)
meta_ref <- shape_df(
  par.drm = list(a = a_base, b = b_base, g = g_base),
  cats    = 2,
  model   = "3PLM"
)

# Focal group item metadata: DIF injected into items 1-6
a_focal      <- a_base
b_focal      <- b_base
b_focal[1:2] <- b_base[1:2] + 0.5          # uniform DIF: harder for focal
a_focal[3:4] <- a_base[3:4] * 0.5          # nonuniform DIF: less discriminating
a_focal[5:6] <- a_base[5:6] * 0.5          # mixed DIF: less discriminating ...
b_focal[5:6] <- b_base[5:6] + 0.5          # ... and harder for focal

meta_foc <- shape_df(
  par.drm = list(a = a_focal, b = b_focal, g = g_base),
  cats    = 2,
  model   = "3PLM"
)

# Simulate ability parameters: reference N(0,1), focal N(-0.3,1)
theta_R <- rnorm(N_R, mean =  0.0, sd = 1)
theta_F <- rnorm(N_F, mean = -0.3, sd = 1)

# Simulate item responses separately for each group
resp_R <- simdat(x = meta_ref, theta = theta_R, D = 1.702)
resp_F <- simdat(x = meta_foc, theta = theta_F, D = 1.702)

# Pool responses and define group membership vector (0 = reference, 1 = focal)
resp_pool <- rbind(resp_R, resp_F)
group_vec <- c(rep(0, N_R), rep(1, N_F))

cat("Pooled data:", nrow(resp_pool), "examinees x", ncol(resp_pool), "items\n")
#> Pooled data: 2000 examinees x 40 items
```

Next, we calibrate item parameters using the pooled data and estimate
ability scores. Because the RDIF framework requires ability estimates
based on pooled (aggregate) item parameters, both calibration and
scoring must be performed on the combined data regardless of group
membership (Lim et al., 2022). The element `test.1` of the calibration
result reports whether the EM algorithm converged; it should be checked
before the estimates are used.

``` r

# Calibrate item parameters from the pooled data
mod_pool <- est_irt(
  data       = resp_pool,
  D          = 1.702,
  model      = "3PLM",
  cats       = 2,
  use.gprior = TRUE,
  gprior     = list(dist = "beta", params = c(4, 16)),
  EmpHist    = TRUE,
  MaxE       = 1000,    # allow extra EM cycles
  verbose    = FALSE
)

# Check that the EM algorithm converged
mod_pool$test.1
#> [1] "Convergence criteria are satisfied."

meta_pool <- mod_pool$par.est

# Estimate ability using ML method
# ML is recommended for RDIF analysis as it provides point estimates
# unaffected by prior distribution assumptions (Lim et al., 2022)
score_pool <- est_score(
  x      = meta_pool,
  data   = resp_pool,
  D      = 1.702,
  method = "ML",
  range  = c(-5, 5)
)$est.theta
```

------------------------------------------------------------------------

## Part 1: Residual-Based DIF Detection (RDIF, `rdif()`)

### Statistical Framework

The RDIF framework (Lim et al., 2022) detects DIF by comparing
item-level residuals, defined as the difference between observed and
model-expected item scores, between the reference and focal groups. For
a dichotomously scored item, the residual for examinee $`h`$ is:

``` math
r_h = x_h - P_h(\hat{\theta}_h)
```

where $`x_h \in \{0, 1\}`$ is the observed response and
$`P_h(\hat{\theta}_h)`$ is the IRT model-predicted probability of a
correct response. Three statistics are computed from these residuals:

**RDIF$`_R`$**: targets **uniform DIF** via differences in mean raw
residuals:

``` math
\text{RDIF}_R = \frac{1}{N_F}\sum_{j=1}^{N_F} r_{Fj} - \frac{1}{N_R}\sum_{i=1}^{N_R} r_{Ri}
```

Under the null hypothesis of no DIF, $`\text{RDIF}_R`$ asymptotically
follows a normal distribution
$`\mathcal{N}(\mu_{\text{RDIF}_R},\, \sigma^2_{\text{RDIF}_R})`$.

**RDIF$`_S`$**: targets **nonuniform DIF** via differences in mean
squared residuals:

``` math
\text{RDIF}_S = \frac{1}{N_F}\sum_{j=1}^{N_F} r_{Fj}^2 - \frac{1}{N_R}\sum_{i=1}^{N_R} r_{Ri}^2
```

Under the null hypothesis, $`\text{RDIF}_S`$ also asymptotically follows
a normal distribution
$`\mathcal{N}(\mu_{\text{RDIF}_S},\, \sigma^2_{\text{RDIF}_S})`$.

**RDIF$`_{RS}`$**: a joint Wald-type statistic that detects **both
uniform and nonuniform DIF** simultaneously. It is based on the
bivariate normality of $`(\text{RDIF}_R, \text{RDIF}_S)`$ and under the
null hypothesis asymptotically follows a $`\chi^2`$ distribution with 2
degrees of freedom.

The analytic expressions for the means and variances of
$`\text{RDIF}_R`$ and $`\text{RDIF}_S`$ are derived from the IRT
model-predicted probabilities (see Lim et al. (2022), Appendix A for
details). Because these moments are computed analytically rather than
empirically, the RDIF framework is highly computationally efficient.

> **Practical guidance:** $`\text{RDIF}_{RS}`$ is the recommended
> primary detection criterion because it is sensitive to both types of
> DIF simultaneously (Lim et al., 2022). Use $`\text{RDIF}_R`$ and
> $`\text{RDIF}_S`$ to characterize the *type* of DIF after a
> significant $`\text{RDIF}_{RS}`$ flag.

### Key arguments of `rdif()`

- `x`: Item metadata data frame containing pooled item parameter
  estimates.
- `data`: Matrix of pooled item response data (rows = examinees, columns
  = items).
- `score`: Numeric vector of pooled ability estimates. If `NULL`,
  [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) estimates
  them internally using the method specified in `method`.
- `group`: Vector of group membership labels (length = number of
  examinees).
- `focal.name`: The label identifying the focal group in `group`.
- `D`: Scaling constant (use the same value as in calibration).
- `alpha`: Significance level for hypothesis testing (default `0.05`).
- `purify`: Logical; whether to apply iterative purification (default
  `FALSE`).
- `purify.by`: Statistic used for purification: `"rdifrs"`, `"rdifr"`,
  or `"rdifs"`.
- `max.iter`: Maximum number of purification iterations (default `10`).
- `method`: Scoring method for ability re-estimation during
  purification. `"ML"` is recommended (default).

### Example 1: RDIF without purification

``` r

rdif_npur <- rdif(
  x          = meta_pool,
  data       = resp_pool,
  score      = score_pool,
  group      = group_vec,
  focal.name = 1,          # 1 = focal group
  D          = 1.702,
  alpha      = 0.05,
  purify     = FALSE,
  verbose    = FALSE
)

# Summary output
print(rdif_npur)
#> 
#> Call:
#> rdif.default(x = meta_pool, data = resp_pool, score = score_pool, 
#>     group = group_vec, focal.name = 1, D = 1.702, alpha = 0.05, 
#>     purify = FALSE, verbose = FALSE)
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 5, 6, 8, 9, 32 
#>   - DIF Items identified by RDIF(S): 
#>     1, 3, 4, 5, 6, 9, 19, 23, 32 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 3, 4, 5, 6, 8, 9, 23, 32 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs      rdifrs p.rdifrs    
#> 1   V1  1000  1000 -0.072   0.000 ***  0.054   0.001 ***  29.038    0.000 ***
#> 2   V2  1000  1000 -0.121   0.000 ***  0.026   0.428      47.073    0.000 ***
#> 3   V3  1000  1000  0.024   0.191      0.054   0.000 ***  28.057    0.000 ***
#> 4   V4  1000  1000 -0.033   0.054   .  0.073   0.000 ***  35.487    0.000 ***
#> 5   V5  1000  1000 -0.136   0.000 ***  0.103   0.000 *** 104.017    0.000 ***
#> 6   V6  1000  1000 -0.121   0.000 ***  0.106   0.000 ***  79.526    0.000 ***
#> 7   V7  1000  1000  0.012   0.382      0.026   0.914       1.065    0.587    
#> 8   V8  1000  1000  0.050   0.004  **  0.006   0.113      10.803    0.004  **
#> 9   V9  1000  1000  0.034   0.018   *  0.010   0.032   *   6.455    0.040   *
#> 10 V10  1000  1000 -0.002   0.925     -0.016   0.469       1.910    0.385    
#> 11 V11  1000  1000 -0.002   0.707      0.009   0.750       0.963    0.618    
#> 12 V12  1000  1000 -0.019   0.290      0.018   0.912       1.120    0.571    
#> 13 V13  1000  1000  0.010   0.567     -0.009   0.785       0.448    0.799    
#> 14 V14  1000  1000  0.015   0.161      0.016   0.544       2.560    0.278    
#> 15 V15  1000  1000  0.017   0.396      0.000   0.309       1.890    0.389    
#> 16 V16  1000  1000 -0.013   0.253      0.032   0.176       2.027    0.363    
#> 17 V17  1000  1000  0.005   0.664      0.025   0.647       2.338    0.311    
#> 18 V18  1000  1000  0.023   0.237     -0.003   0.384       1.401    0.496    
#> 19 V19  1000  1000  0.034   0.051   .  0.014   0.048   *   4.038    0.133    
#> 20 V20  1000  1000  0.013   0.398      0.025   0.907       0.714    0.700    
#> 21 V21  1000  1000 -0.007   0.695     -0.018   0.250       1.403    0.496    
#> 22 V22  1000  1000  0.002   0.834      0.010   0.253       2.647    0.266    
#> 23 V23  1000  1000  0.021   0.194     -0.004   0.002  **  10.236    0.006  **
#> 24 V24  1000  1000  0.016   0.203      0.025   0.981       2.020    0.364    
#> 25 V25  1000  1000  0.001   0.965     -0.008   0.917       0.085    0.958    
#> 26 V26  1000  1000  0.026   0.132      0.010   0.375       2.774    0.250    
#> 27 V27  1000  1000  0.018   0.302     -0.003   0.486       1.189    0.552    
#> 28 V28  1000  1000  0.003   0.854     -0.009   0.694       1.120    0.571    
#> 29 V29  1000  1000 -0.007   0.698     -0.011   0.974       0.440    0.802    
#> 30 V30  1000  1000  0.001   0.973      0.005   0.134       2.403    0.301    
#> 31 V31  1000  1000 -0.003   0.866      0.001   0.576       1.354    0.508    
#> 32 V32  1000  1000  0.040   0.025   *  0.022   0.033   *   6.345    0.042   *
#> 33 V33  1000  1000  0.016   0.366      0.007   0.201       1.644    0.440    
#> 34 V34  1000  1000 -0.003   0.591      0.016   0.206       2.287    0.319    
#> 35 V35  1000  1000  0.027   0.182      0.012   0.758       1.780    0.411    
#> 36 V36  1000  1000  0.008   0.551      0.017   0.191       1.802    0.406    
#> 37 V37  1000  1000  0.024   0.127      0.018   0.509       2.425    0.298    
#> 38 V38  1000  1000  0.014   0.486     -0.005   0.471       1.479    0.477    
#> 39 V39  1000  1000  0.029   0.091   .  0.012   0.150       3.374    0.185    
#> 40 V40  1000  1000  0.010   0.403      0.019   0.627       0.768    0.681    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented.
```

``` r

# Full table of RDIF statistics for all items
rdif_npur$no_purify$dif_stat
#>     id   rdifr z.rdifr   rdifs z.rdifs   rdifrs p.rdifr p.rdifs p.rdifrs n.ref
#> 1   V1 -0.0715 -5.2903  0.0536  3.4352  29.0376  0.0000  0.0006   0.0000  1000
#> 2   V2 -0.1207 -6.8244  0.0264  0.7925  47.0729  0.0000  0.4281   0.0000  1000
#> 3   V3  0.0245  1.3072  0.0537  4.1500  28.0571  0.1911  0.0000   0.0000  1000
#> 4   V4 -0.0333 -1.9292  0.0730  5.9229  35.4870  0.0537  0.0000   0.0000  1000
#> 5   V5 -0.1358 -9.8322  0.1027  9.3258 104.0171  0.0000  0.0000   0.0000  1000
#> 6   V6 -0.1208 -7.3043  0.1063  8.8836  79.5258  0.0000  0.0000   0.0000  1000
#> 7   V7  0.0120  0.8744  0.0263 -0.1077   1.0652  0.3819  0.9142   0.5871  1000
#> 8   V8  0.0495  2.8712  0.0056 -1.5867  10.8034  0.0041  0.1126   0.0045  1000
#> 9   V9  0.0338  2.3768  0.0095 -2.1403   6.4545  0.0175  0.0323   0.0397  1000
#> 10 V10 -0.0017 -0.0947 -0.0155 -0.7249   1.9099  0.9246  0.4685   0.3848  1000
#> 11 V11 -0.0024 -0.3760  0.0089 -0.3192   0.9626  0.7069  0.7495   0.6180  1000
#> 12 V12 -0.0187 -1.0585  0.0180  0.1104   1.1205  0.2898  0.9121   0.5711  1000
#> 13 V13  0.0102  0.5732 -0.0094  0.2733   0.4481  0.5665  0.7846   0.7993  1000
#> 14 V14  0.0146  1.4017  0.0163 -0.6066   2.5600  0.1610  0.5441   0.2780  1000
#> 15 V15  0.0171  0.8480 -0.0003 -1.0172   1.8900  0.3964  0.3091   0.3887  1000
#> 16 V16 -0.0129 -1.1431  0.0325  1.3539   2.0269  0.2530  0.1758   0.3630  1000
#> 17 V17  0.0050  0.4344  0.0251  0.4579   2.3383  0.6640  0.6470   0.3106  1000
#> 18 V18  0.0233  1.1834 -0.0031  0.8707   1.4005  0.2366  0.3839   0.4965  1000
#> 19 V19  0.0340  1.9515  0.0140  1.9766   4.0379  0.0510  0.0481   0.1328  1000
#> 20 V20  0.0130  0.8443  0.0255 -0.1163   0.7139  0.3985  0.9074   0.6998  1000
#> 21 V21 -0.0070 -0.3923 -0.0176 -1.1504   1.4031  0.6949  0.2500   0.4958  1000
#> 22 V22  0.0020  0.2095  0.0100 -1.1436   2.6468  0.8341  0.2528   0.2662  1000
#> 23 V23  0.0213  1.2999 -0.0037 -3.1603  10.2360  0.1936  0.0016   0.0060  1000
#> 24 V24  0.0159  1.2729  0.0255  0.0238   2.0200  0.2030  0.9810   0.3642  1000
#> 25 V25  0.0008  0.0433 -0.0078 -0.1043   0.0849  0.9654  0.9169   0.9584  1000
#> 26 V26  0.0262  1.5075  0.0096 -0.8870   2.7742  0.1317  0.3751   0.2498  1000
#> 27 V27  0.0185  1.0317 -0.0027  0.6964   1.1893  0.3022  0.4862   0.5517  1000
#> 28 V28  0.0032  0.1844 -0.0089 -0.3929   1.1197  0.8537  0.6944   0.5713  1000
#> 29 V29 -0.0072 -0.3879 -0.0105  0.0331   0.4404  0.6981  0.9736   0.8023  1000
#> 30 V30  0.0006  0.0342  0.0053 -1.4989   2.4034  0.9727  0.1339   0.3007  1000
#> 31 V31 -0.0029 -0.1688  0.0011  0.5587   1.3541  0.8660  0.5764   0.5081  1000
#> 32 V32  0.0405  2.2397  0.0216  2.1288   6.3449  0.0251  0.0333   0.0419  1000
#> 33 V33  0.0163  0.9035  0.0067  1.2791   1.6438  0.3663  0.2008   0.4396  1000
#> 34 V34 -0.0034 -0.5368  0.0157  1.2634   2.2870  0.5914  0.2064   0.3187  1000
#> 35 V35  0.0266  1.3341  0.0117 -0.3086   1.7799  0.1822  0.7576   0.4107  1000
#> 36 V36  0.0080  0.5959  0.0168 -1.3079   1.8019  0.5512  0.1909   0.4062  1000
#> 37 V37  0.0241  1.5278  0.0183 -0.6608   2.4248  0.1266  0.5087   0.2975  1000
#> 38 V38  0.0138  0.6966 -0.0055 -0.7209   1.4787  0.4860  0.4710   0.4774  1000
#> 39 V39  0.0289  1.6905  0.0123 -1.4386   3.3742  0.0909  0.1503   0.1851  1000
#> 40 V40  0.0099  0.8360  0.0188 -0.4853   0.7679  0.4032  0.6275   0.6812  1000
#>    n.foc n.total
#> 1   1000    2000
#> 2   1000    2000
#> 3   1000    2000
#> 4   1000    2000
#> 5   1000    2000
#> 6   1000    2000
#> 7   1000    2000
#> 8   1000    2000
#> 9   1000    2000
#> 10  1000    2000
#> 11  1000    2000
#> 12  1000    2000
#> 13  1000    2000
#> 14  1000    2000
#> 15  1000    2000
#> 16  1000    2000
#> 17  1000    2000
#> 18  1000    2000
#> 19  1000    2000
#> 20  1000    2000
#> 21  1000    2000
#> 22  1000    2000
#> 23  1000    2000
#> 24  1000    2000
#> 25  1000    2000
#> 26  1000    2000
#> 27  1000    2000
#> 28  1000    2000
#> 29  1000    2000
#> 30  1000    2000
#> 31  1000    2000
#> 32  1000    2000
#> 33  1000    2000
#> 34  1000    2000
#> 35  1000    2000
#> 36  1000    2000
#> 37  1000    2000
#> 38  1000    2000
#> 39  1000    2000
#> 40  1000    2000

# Items flagged by each statistic
rdif_npur$no_purify$dif_item
#> $rdifr
#> [1]  1  2  5  6  8  9 32
#> 
#> $rdifs
#> [1]  1  3  4  5  6  9 19 23 32
#> 
#> $rdifrs
#>  [1]  1  2  3  4  5  6  8  9 23 32
```

The output contains a `dif_stat` data frame with the following columns
for each item: `rdifr` (RDIF$`_R`$ statistic), `z.rdifr` (standardized
RDIF$`_R`$), `rdifs` (RDIF$`_S`$ statistic), `z.rdifs` (standardized
RDIF$`_S`$), `rdifrs` (RDIF$`_{RS}`$ statistic), `p.rdifr`, `p.rdifs`,
`p.rdifrs` (corresponding p-values), and group sample sizes. Items with
p-values below `alpha` are flagged as DIF.

In this simulation, RDIF$`_{RS}`$ flags all six items with injected DIF
(items 1-6) and four DIF-free items (8, 9, 23, and 32). Among the items
with injected DIF, RDIF$`_R`$ flags the items with a difficulty shift
(items 1, 2, 5, and 6), and RDIF$`_S`$ flags the items with a lower
focal discrimination (items 3, 4, 5, and 6) together with item 1. Each
of these two statistics also flags some DIF-free items (RDIF$`_R`$: 8,
9, and 32; RDIF$`_S`$: 9, 19, 23, and 32). The DIF-free items that are
flagged are false positives, which the purification in the next example
is designed to remove.

### Example 2: RDIF with purification

When DIF items are present in the test, they can contaminate the ability
estimates used to compute the RDIF statistics, inflating Type I error
rates for DIF-free items. Iterative purification addresses this by
progressively removing flagged DIF items from ability re-estimation. At
each iteration, the most significant flagged item (the item with the
largest absolute standardized statistic for `"rdifr"` and `"rdifs"`, the
largest statistic for `"rdifrs"` and for the GRDIF statistics, and the
largest absolute `z.beta` for
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md)) is
removed, the abilities are re-estimated from the remaining items, and
the DIF tests are repeated. The procedure stops when no further item is
flagged or when `max.iter` is reached (Lim et al., 2022).

``` r

rdif_pur <- rdif(
  x          = meta_pool,
  data       = resp_pool,
  score      = score_pool,   # initial ability estimates (pre-purification)
  group      = group_vec,
  focal.name = 1,
  D          = 1.702,
  alpha      = 0.05,
  purify     = TRUE,
  purify.by  = "rdifrs",     # use RDIF_RS to drive purification
  max.iter   = 20,
  method     = "ML",         # re-estimate abilities with ML at each iteration
  range      = c(-5, 5),
  verbose    = FALSE
)

# Summary output
print(rdif_pur)
#> 
#> Call:
#> rdif.default(x = meta_pool, data = resp_pool, score = score_pool, 
#>     group = group_vec, focal.name = 1, D = 1.702, alpha = 0.05, 
#>     purify = TRUE, purify.by = "rdifrs", max.iter = 20, method = "ML", 
#>     range = c(-5, 5), verbose = FALSE)
#> 
#> DIF analysis using three RDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by RDIF(R): 
#>     1, 2, 5, 6, 8, 9, 32 
#>   - DIF Items identified by RDIF(S): 
#>     1, 3, 4, 5, 6, 9, 19, 23, 32 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 3, 4, 5, 6, 8, 9, 23, 32 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs      rdifrs p.rdifrs    
#> 1   V1  1000  1000 -0.072   0.000 ***  0.054   0.001 ***  29.038    0.000 ***
#> 2   V2  1000  1000 -0.121   0.000 ***  0.026   0.428      47.073    0.000 ***
#> 3   V3  1000  1000  0.024   0.191      0.054   0.000 ***  28.057    0.000 ***
#> 4   V4  1000  1000 -0.033   0.054   .  0.073   0.000 ***  35.487    0.000 ***
#> 5   V5  1000  1000 -0.136   0.000 ***  0.103   0.000 *** 104.017    0.000 ***
#> 6   V6  1000  1000 -0.121   0.000 ***  0.106   0.000 ***  79.526    0.000 ***
#> 7   V7  1000  1000  0.012   0.382      0.026   0.914       1.065    0.587    
#> 8   V8  1000  1000  0.050   0.004  **  0.006   0.113      10.803    0.004  **
#> 9   V9  1000  1000  0.034   0.018   *  0.010   0.032   *   6.455    0.040   *
#> 10 V10  1000  1000 -0.002   0.925     -0.016   0.469       1.910    0.385    
#> 11 V11  1000  1000 -0.002   0.707      0.009   0.750       0.963    0.618    
#> 12 V12  1000  1000 -0.019   0.290      0.018   0.912       1.120    0.571    
#> 13 V13  1000  1000  0.010   0.567     -0.009   0.785       0.448    0.799    
#> 14 V14  1000  1000  0.015   0.161      0.016   0.544       2.560    0.278    
#> 15 V15  1000  1000  0.017   0.396      0.000   0.309       1.890    0.389    
#> 16 V16  1000  1000 -0.013   0.253      0.032   0.176       2.027    0.363    
#> 17 V17  1000  1000  0.005   0.664      0.025   0.647       2.338    0.311    
#> 18 V18  1000  1000  0.023   0.237     -0.003   0.384       1.401    0.496    
#> 19 V19  1000  1000  0.034   0.051   .  0.014   0.048   *   4.038    0.133    
#> 20 V20  1000  1000  0.013   0.398      0.025   0.907       0.714    0.700    
#> 21 V21  1000  1000 -0.007   0.695     -0.018   0.250       1.403    0.496    
#> 22 V22  1000  1000  0.002   0.834      0.010   0.253       2.647    0.266    
#> 23 V23  1000  1000  0.021   0.194     -0.004   0.002  **  10.236    0.006  **
#> 24 V24  1000  1000  0.016   0.203      0.025   0.981       2.020    0.364    
#> 25 V25  1000  1000  0.001   0.965     -0.008   0.917       0.085    0.958    
#> 26 V26  1000  1000  0.026   0.132      0.010   0.375       2.774    0.250    
#> 27 V27  1000  1000  0.018   0.302     -0.003   0.486       1.189    0.552    
#> 28 V28  1000  1000  0.003   0.854     -0.009   0.694       1.120    0.571    
#> 29 V29  1000  1000 -0.007   0.698     -0.011   0.974       0.440    0.802    
#> 30 V30  1000  1000  0.001   0.973      0.005   0.134       2.403    0.301    
#> 31 V31  1000  1000 -0.003   0.866      0.001   0.576       1.354    0.508    
#> 32 V32  1000  1000  0.040   0.025   *  0.022   0.033   *   6.345    0.042   *
#> 33 V33  1000  1000  0.016   0.366      0.007   0.201       1.644    0.440    
#> 34 V34  1000  1000 -0.003   0.591      0.016   0.206       2.287    0.319    
#> 35 V35  1000  1000  0.027   0.182      0.012   0.758       1.780    0.411    
#> 36 V36  1000  1000  0.008   0.551      0.017   0.191       1.802    0.406    
#> 37 V37  1000  1000  0.024   0.127      0.018   0.509       2.425    0.298    
#> 38 V38  1000  1000  0.014   0.486     -0.005   0.471       1.479    0.477    
#> 39 V39  1000  1000  0.029   0.091   .  0.012   0.150       3.374    0.185    
#> 40 V40  1000  1000  0.010   0.403      0.019   0.627       0.768    0.681    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 6
#>   - RDIF statistic used for purification: RDIF(RS)
#>   - DIF Items identified by RDIF(RS): 
#>     1, 2, 3, 4, 5, 6 
#>   - RDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs      rdifrs
#> 1   V1      4  1000  1000 -0.079   0.000 ***  0.053   0.000 ***  34.429
#> 2   V2      2  1000  1000 -0.128   0.000 ***  0.027   0.325      53.515
#> 3   V3      5  1000  1000  0.009   0.651      0.052   0.000 ***  26.278
#> 4   V4      3  1000  1000 -0.046   0.008  **  0.071   0.000 ***  38.174
#> 5   V5      0  1000  1000 -0.136   0.000 ***  0.103   0.000 *** 104.017
#> 6   V6      1  1000  1000 -0.123   0.000 ***  0.107   0.000 ***  82.411
#> 7   V7      6  1000  1000  0.001   0.964      0.026   0.510       0.759
#> 8   V8      6  1000  1000  0.026   0.132      0.003   0.224       3.772
#> 9   V9      6  1000  1000  0.021   0.132      0.008   0.135       2.829
#> 10 V10      6  1000  1000 -0.006   0.736     -0.015   0.393       1.527
#> 11 V11      6  1000  1000 -0.004   0.493      0.009   0.973       0.985
#> 12 V12      6  1000  1000 -0.040   0.022   *  0.018   0.492       5.466
#> 13 V13      6  1000  1000  0.003   0.886     -0.009   0.990       0.054
#> 14 V14      6  1000  1000  0.009   0.388      0.017   0.944       2.074
#> 15 V15      6  1000  1000  0.001   0.945      0.000   0.356       0.869
#> 16 V16      6  1000  1000 -0.021   0.057   .  0.033   0.033   *   5.221
#> 17 V17      6  1000  1000 -0.001   0.923      0.023   0.424       1.619
#> 18 V18      6  1000  1000  0.014   0.484     -0.003   0.598       0.490
#> 19 V19      6  1000  1000  0.031   0.075   .  0.014   0.065   .   3.458
#> 20 V20      6  1000  1000 -0.007   0.649      0.027   0.260       1.341
#> 21 V21      6  1000  1000 -0.024   0.182     -0.016   0.181       2.318
#> 22 V22      6  1000  1000 -0.002   0.797      0.010   0.502       2.077
#> 23 V23      6  1000  1000  0.001   0.947     -0.003   0.022   *   5.527
#> 24 V24      6  1000  1000  0.004   0.758      0.024   0.528       0.810
#> 25 V25      6  1000  1000 -0.004   0.825     -0.008   0.794       0.068
#> 26 V26      6  1000  1000  0.004   0.816      0.012   0.999       0.055
#> 27 V27      6  1000  1000  0.013   0.487     -0.002   0.600       0.500
#> 28 V28      6  1000  1000 -0.002   0.897     -0.009   0.551       0.854
#> 29 V29      6  1000  1000 -0.016   0.389     -0.011   0.735       1.022
#> 30 V30      6  1000  1000 -0.017   0.365      0.006   0.394       2.133
#> 31 V31      6  1000  1000 -0.011   0.527      0.000   0.793       2.010
#> 32 V32      6  1000  1000  0.022   0.227      0.020   0.079   .   3.237
#> 33 V33      6  1000  1000  0.003   0.869      0.007   0.340       1.271
#> 34 V34      6  1000  1000 -0.005   0.413      0.015   0.179       2.005
#> 35 V35      6  1000  1000  0.011   0.567      0.012   0.896       0.402
#> 36 V36      6  1000  1000 -0.003   0.834      0.017   0.603       0.749
#> 37 V37      6  1000  1000  0.004   0.797      0.014   0.772       0.122
#> 38 V38      6  1000  1000 -0.002   0.920     -0.006   0.386       0.783
#> 39 V39      6  1000  1000  0.013   0.460      0.012   0.504       0.676
#> 40 V40      6  1000  1000  0.003   0.815      0.018   0.936       0.222
#>    p.rdifrs    
#> 1     0.000 ***
#> 2     0.000 ***
#> 3     0.000 ***
#> 4     0.000 ***
#> 5     0.000 ***
#> 6     0.000 ***
#> 7     0.684    
#> 8     0.152    
#> 9     0.243    
#> 10    0.466    
#> 11    0.611    
#> 12    0.065   .
#> 13    0.973    
#> 14    0.354    
#> 15    0.648    
#> 16    0.074   .
#> 17    0.445    
#> 18    0.783    
#> 19    0.178    
#> 20    0.511    
#> 21    0.314    
#> 22    0.354    
#> 23    0.063   .
#> 24    0.667    
#> 25    0.967    
#> 26    0.973    
#> 27    0.779    
#> 28    0.653    
#> 29    0.600    
#> 30    0.344    
#> 31    0.366    
#> 32    0.198    
#> 33    0.530    
#> 34    0.367    
#> 35    0.818    
#> 36    0.688    
#> 37    0.941    
#> 38    0.676    
#> 39    0.713    
#> 40    0.895    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05

# DIF items identified after purification
rdif_pur$with_purify$dif_item
#> [1] 1 2 3 4 5 6
```

The `with_purify` component contains the final results after
purification. Its `n.iter` field reports the total number of iterations
that were performed, and `complete` indicates whether purification
finished before reaching `max.iter`. Here the purification finishes
after 6 iterations and flags exactly the six items with injected DIF
(items 1-6); the false positives of the previous example are gone.

Do not confuse the `n.iter` field with the `n.iter` column of
`with_purify$dif_stat`, which gives, for each item, the iteration whose
result is reported: 0 for an item removed at the first iteration (its
statistics come from the initial analysis), and the final iteration for
the items that were never removed.

------------------------------------------------------------------------

## Part 2: Multiple-Group DIF Detection (GRDIF, `grdif()`)

### Statistical Framework

The GRDIF framework (Lim et al., 2024) generalizes the two-group RDIF
approach to simultaneously detect DIF across $`G \geq 2`$ groups. Three
statistics, $`\text{GRDIF}_R`$, $`\text{GRDIF}_S`$, and
$`\text{GRDIF}_{RS}`$, are constructed by applying a contrast matrix
$`\mathbf{C}`$ to the stacked vector of per-group mean raw residuals
(MRR) and mean squared residuals (MSR), exploiting their asymptotic
multivariate normality.

Under the null hypothesis that no DIF exists between any pair of groups,
the three statistics follow asymptotic $`\chi^2`$ distributions:

``` math
\text{GRDIF}_R \xrightarrow{d} \chi^2_{G-1}, \quad
  \text{GRDIF}_S \xrightarrow{d} \chi^2_{G-1}, \quad
  \text{GRDIF}_{RS} \xrightarrow{d} \chi^2_{2(G-1)}
```

where $`G`$ is the total number of groups. For the two-group special
case ($`G = 2`$), the GRDIF statistics reduce exactly to the original
RDIF statistics.

When a significant omnibus result is found, optional post-hoc pairwise
RDIF analyses (controlled via `post.hoc = TRUE`) can identify which
specific pairs of groups exhibit DIF.

> **Note on polytomous items:** Like
> [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
> [`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md) also
> accepts polytomous item response data. However, formal validation of
> the GRDIF framework for polytomous items has not yet been published;
> research on this extension is currently ongoing.

### Key arguments of `grdif()`

- `focal.name`: A *vector* of labels for all focal groups (e.g.,
  `focal.name = c("G2", "G3")` when there are two focal groups).
- `post.hoc`: Logical; whether to perform post-hoc pairwise RDIF tests
  on flagged items (default `TRUE`). Useful for identifying which group
  pairs drive the omnibus DIF signal.
- Other arguments (`x`, `data`, `score`, `group`, `D`, `alpha`,
  `purify`, `purify.by`, `max.iter`, `method`) follow the same
  conventions as
  [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md).

### Example 1: GRDIF without purification

We extend the simulation to a three-group scenario by adding a second
focal group. Items 1-2 exhibit uniform DIF (b-parameters shifted upward)
in both focal groups, with a larger shift for Focal Group 3 ($`+0.8`$
versus $`+0.5`$ for Focal Group 2).

``` r

N_g <- 700   # examinees per group

# Ability parameters: G1 = reference, G2 and G3 = focal groups (with impact)
theta_G1 <- rnorm(N_g,  0.0, 1)
theta_G2 <- rnorm(N_g, -0.2, 1)
theta_G3 <- rnorm(N_g, -0.4, 1)

# Item parameters per group
# G1: no DIF (base parameters from the Setup section)
# G2: items 1-2 harder for this focal group (uniform DIF, b + 0.5)
# G3: items 1-2 even harder for this focal group (uniform DIF, b + 0.8)
meta_G1 <- meta_ref
meta_G2 <- meta_ref
meta_G3 <- meta_ref
meta_G2$par.2[1:2] <- meta_ref$par.2[1:2] + 0.5
meta_G3$par.2[1:2] <- meta_ref$par.2[1:2] + 0.8

# Simulate responses for each group
resp_G1 <- simdat(x = meta_G1, theta = theta_G1, D = 1.702)
resp_G2 <- simdat(x = meta_G2, theta = theta_G2, D = 1.702)
resp_G3 <- simdat(x = meta_G3, theta = theta_G3, D = 1.702)

# Pool data and define group membership
resp_3g <- rbind(resp_G1, resp_G2, resp_G3)
grp_3g  <- c(rep("G1", N_g), rep("G2", N_g), rep("G3", N_g))

cat("Pooled data:", nrow(resp_3g), "examinees x", ncol(resp_3g), "items\n")
#> Pooled data: 2100 examinees x 40 items
```

``` r

# Calibrate using pooled data
mod_3g <- est_irt(
  data       = resp_3g,
  D          = 1.702,
  model      = "3PLM",
  cats       = 2,
  use.gprior = TRUE,
  gprior     = list(dist = "beta", params = c(4, 16)),
  EmpHist    = TRUE,
  MaxE       = 1000,    # allow extra EM cycles
  verbose    = FALSE
)

# Check that the EM algorithm converged
mod_3g$test.1
#> [1] "Convergence criteria are satisfied."

meta_3g  <- mod_3g$par.est

# Estimate pooled ability scores using ML
score_3g <- est_score(
  x      = meta_3g,
  data   = resp_3g,
  D      = 1.702,
  method = "ML",
  range  = c(-5, 5)
)$est.theta
```

``` r

grdif_npur <- grdif(
  x          = meta_3g,
  data       = resp_3g,
  score      = score_3g,
  group      = grp_3g,
  focal.name = c("G2", "G3"),   # both G2 and G3 are focal groups
  D          = 1.702,
  alpha      = 0.05,
  purify     = FALSE,
  post.hoc   = TRUE,            # run pairwise follow-up for flagged items
  verbose    = FALSE
)

# Summary output (omnibus test results)
print(grdif_npur)
#> 
#> Call:
#> grdif.default(x = meta_3g, data = resp_3g, score = score_3g, 
#>     group = grp_3g, focal.name = c("G2", "G3"), D = 1.702, alpha = 0.05, 
#>     purify = FALSE, post.hoc = TRUE, verbose = FALSE)
#> 
#> DIF analysis using three GRDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by GRDIF(R): 
#>     1, 2 
#>   - DIF Items identified by GRDIF(S): 
#>     1, 2 
#>   - DIF Items identified by GRDIF(RS): 
#>     1, 2 
#>   - GRDIF Statistics: 
#> 
#>     id n.ref n.foc1 n.foc2 grdifr p.grdifr     grdifs p.grdifs     grdifrs
#> 1   V1   700    700    700 44.868    0.000 *** 32.542    0.000 ***  46.683
#> 2   V2   700    700    700 67.755    0.000 ***  9.828    0.007  **  82.604
#> 3   V3   700    700    700  0.312    0.856      0.992    0.609       1.214
#> 4   V4   700    700    700  3.779    0.151      1.519    0.468       5.142
#> 5   V5   700    700    700  0.486    0.784      0.585    0.747       0.856
#> 6   V6   700    700    700  0.864    0.649      0.442    0.802       1.873
#> 7   V7   700    700    700  1.394    0.498      1.390    0.499       5.932
#> 8   V8   700    700    700  5.127    0.077   .  0.651    0.722       5.602
#> 9   V9   700    700    700  0.764    0.683      0.839    0.657       2.062
#> 10 V10   700    700    700  0.132    0.936      0.602    0.740       1.780
#> 11 V11   700    700    700  2.998    0.223      0.983    0.612       3.522
#> 12 V12   700    700    700  1.180    0.554      1.807    0.405       2.799
#> 13 V13   700    700    700  0.389    0.823      0.100    0.951       0.536
#> 14 V14   700    700    700  0.426    0.808      1.152    0.562       1.280
#> 15 V15   700    700    700  2.404    0.301      0.201    0.905       2.490
#> 16 V16   700    700    700  2.627    0.269      1.520    0.468       3.136
#> 17 V17   700    700    700  0.910    0.634      2.532    0.282       3.079
#> 18 V18   700    700    700  0.203    0.904      0.588    0.746       3.003
#> 19 V19   700    700    700  0.894    0.639      0.617    0.735       4.711
#> 20 V20   700    700    700  0.124    0.940      0.400    0.819       0.559
#> 21 V21   700    700    700  0.649    0.723      1.411    0.494       1.460
#> 22 V22   700    700    700  0.297    0.862      0.317    0.853       2.652
#> 23 V23   700    700    700  1.537    0.464      1.202    0.548       2.903
#> 24 V24   700    700    700  0.622    0.733      4.768    0.092   .   5.282
#> 25 V25   700    700    700  2.413    0.299      0.856    0.652       3.055
#> 26 V26   700    700    700  1.168    0.558      0.071    0.965       1.205
#> 27 V27   700    700    700  2.722    0.256      1.246    0.536       3.368
#> 28 V28   700    700    700  0.804    0.669      0.943    0.624       2.012
#> 29 V29   700    700    700  5.237    0.073   .  3.084    0.214       5.666
#> 30 V30   700    700    700  3.259    0.196      0.345    0.842       3.680
#> 31 V31   700    700    700  1.849    0.397      2.611    0.271       3.125
#> 32 V32   700    700    700  0.740    0.691      0.789    0.674       0.954
#> 33 V33   700    700    700  0.419    0.811      2.090    0.352       2.380
#> 34 V34   700    700    700  3.110    0.211      5.182    0.075   .   5.436
#> 35 V35   700    700    700  0.313    0.855      4.591    0.101       4.608
#> 36 V36   700    700    700  2.016    0.365      3.152    0.207       3.582
#> 37 V37   700    700    700  5.327    0.070   .  4.923    0.085   .   8.555
#> 38 V38   700    700    700  2.733    0.255      4.922    0.085   .   5.430
#> 39 V39   700    700    700  3.422    0.181      2.458    0.292       5.612
#> 40 V40   700    700    700  3.603    0.165      1.429    0.490       4.573
#>    p.grdifrs    
#> 1      0.000 ***
#> 2      0.000 ***
#> 3      0.876    
#> 4      0.273    
#> 5      0.931    
#> 6      0.759    
#> 7      0.204    
#> 8      0.231    
#> 9      0.724    
#> 10     0.776    
#> 11     0.475    
#> 12     0.592    
#> 13     0.970    
#> 14     0.865    
#> 15     0.646    
#> 16     0.535    
#> 17     0.545    
#> 18     0.557    
#> 19     0.318    
#> 20     0.968    
#> 21     0.834    
#> 22     0.618    
#> 23     0.574    
#> 24     0.260    
#> 25     0.549    
#> 26     0.877    
#> 27     0.498    
#> 28     0.734    
#> 29     0.226    
#> 30     0.451    
#> 31     0.537    
#> 32     0.917    
#> 33     0.666    
#> 34     0.245    
#> 35     0.330    
#> 36     0.466    
#> 37     0.073   .
#> 38     0.246    
#> 39     0.230    
#> 40     0.334    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented.

# Items flagged by each GRDIF statistic
grdif_npur$no_purify$dif_item
#> $grdifr
#> [1] 1 2
#> 
#> $grdifs
#> [1] 1 2
#> 
#> $grdifrs
#> [1] 1 2
```

``` r

# Post-hoc pairwise RDIF results for DIF-flagged items
# (identifies which specific group pairs drive the DIF signal)
grdif_npur$no_purify$post.hoc
#> $by.grdifr
#>   id group.pair   rdifr z.rdifr  rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs
#> 1 V1    G1 & G2 -0.0595 -3.3890 0.0568  3.3734 13.0361  0.0007  0.0007   0.0015
#> 2 V1    G1 & G3 -0.1194 -6.6780 0.0885  5.6615 46.2794  0.0000  0.0000   0.0000
#> 3 V1    G2 & G3 -0.0600 -3.1916 0.0317  2.2964 10.1998  0.0014  0.0217   0.0061
#> 4 V2    G1 & G2 -0.1116 -5.1701 0.0267  1.5508 30.2416  0.0000  0.1210   0.0000
#> 5 V2    G1 & G3 -0.1766 -8.0969 0.0495  3.1349 80.9867  0.0000  0.0017   0.0000
#> 6 V2    G2 & G3 -0.0650 -2.9333 0.0228  1.5540 13.4864  0.0034  0.1202   0.0012
#>   n.ref n.foc n.total
#> 1   700   700    1400
#> 2   700   700    1400
#> 3   700   700    1400
#> 4   700   700    1400
#> 5   700   700    1400
#> 6   700   700    1400
#> 
#> $by.grdifs
#>   id group.pair   rdifr z.rdifr  rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs
#> 1 V1    G1 & G2 -0.0595 -3.3890 0.0568  3.3734 13.0361  0.0007  0.0007   0.0015
#> 2 V1    G1 & G3 -0.1194 -6.6780 0.0885  5.6615 46.2794  0.0000  0.0000   0.0000
#> 3 V1    G2 & G3 -0.0600 -3.1916 0.0317  2.2964 10.1998  0.0014  0.0217   0.0061
#> 4 V2    G1 & G2 -0.1116 -5.1701 0.0267  1.5508 30.2416  0.0000  0.1210   0.0000
#> 5 V2    G1 & G3 -0.1766 -8.0969 0.0495  3.1349 80.9867  0.0000  0.0017   0.0000
#> 6 V2    G2 & G3 -0.0650 -2.9333 0.0228  1.5540 13.4864  0.0034  0.1202   0.0012
#>   n.ref n.foc n.total
#> 1   700   700    1400
#> 2   700   700    1400
#> 3   700   700    1400
#> 4   700   700    1400
#> 5   700   700    1400
#> 6   700   700    1400
#> 
#> $by.grdifrs
#>   id group.pair   rdifr z.rdifr  rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs
#> 1 V1    G1 & G2 -0.0595 -3.3890 0.0568  3.3734 13.0361  0.0007  0.0007   0.0015
#> 2 V1    G1 & G3 -0.1194 -6.6780 0.0885  5.6615 46.2794  0.0000  0.0000   0.0000
#> 3 V1    G2 & G3 -0.0600 -3.1916 0.0317  2.2964 10.1998  0.0014  0.0217   0.0061
#> 4 V2    G1 & G2 -0.1116 -5.1701 0.0267  1.5508 30.2416  0.0000  0.1210   0.0000
#> 5 V2    G1 & G3 -0.1766 -8.0969 0.0495  3.1349 80.9867  0.0000  0.0017   0.0000
#> 6 V2    G2 & G3 -0.0650 -2.9333 0.0228  1.5540 13.4864  0.0034  0.1202   0.0012
#>   n.ref n.foc n.total
#> 1   700   700    1400
#> 2   700   700    1400
#> 3   700   700    1400
#> 4   700   700    1400
#> 5   700   700    1400
#> 6   700   700    1400
```

The `post.hoc` component reports pairwise RDIF$`_R`$, RDIF$`_S`$, and
RDIF$`_{RS}`$ statistics and p-values for each pair of groups among the
items that were flagged by the omnibus test. This helps pinpoint whether
DIF occurs between the reference and one focal group, between the
reference and all focal groups, or between focal groups only.

In this simulation, all three GRDIF statistics flag items 1 and 2 and no
other item. The post-hoc results show significant RDIF$`_R`$ for every
pair of groups on both items, and the largest differences involve Focal
Group 3, whose difficulty shift is the largest.

### Example 2: GRDIF with purification

``` r

grdif_pur <- grdif(
  x          = meta_3g,
  data       = resp_3g,
  score      = score_3g,       # initial ability estimates (pre-purification)
  group      = grp_3g,
  focal.name = c("G2", "G3"),
  D          = 1.702,
  alpha      = 0.05,
  purify     = TRUE,
  purify.by  = "grdifr",       # use GRDIF_R to drive purification
  max.iter   = 20,
  method     = "ML",           # re-estimate abilities with ML at each iteration
  range      = c(-5, 5),
  post.hoc   = TRUE,
  verbose    = FALSE
)

# Summary output
print(grdif_pur)
#> 
#> Call:
#> grdif.default(x = meta_3g, data = resp_3g, score = score_3g, 
#>     group = grp_3g, focal.name = c("G2", "G3"), D = 1.702, alpha = 0.05, 
#>     purify = TRUE, purify.by = "grdifr", max.iter = 20, post.hoc = TRUE, 
#>     method = "ML", range = c(-5, 5), verbose = FALSE)
#> 
#> DIF analysis using three GRDIF statistics 
#> 
#>  1. Without purification 
#> 
#>   - DIF Items identified by GRDIF(R): 
#>     1, 2 
#>   - DIF Items identified by GRDIF(S): 
#>     1, 2 
#>   - DIF Items identified by GRDIF(RS): 
#>     1, 2 
#>   - GRDIF Statistics: 
#> 
#>     id n.ref n.foc1 n.foc2 grdifr p.grdifr     grdifs p.grdifs     grdifrs
#> 1   V1   700    700    700 44.868    0.000 *** 32.542    0.000 ***  46.683
#> 2   V2   700    700    700 67.755    0.000 ***  9.828    0.007  **  82.604
#> 3   V3   700    700    700  0.312    0.856      0.992    0.609       1.214
#> 4   V4   700    700    700  3.779    0.151      1.519    0.468       5.142
#> 5   V5   700    700    700  0.486    0.784      0.585    0.747       0.856
#> 6   V6   700    700    700  0.864    0.649      0.442    0.802       1.873
#> 7   V7   700    700    700  1.394    0.498      1.390    0.499       5.932
#> 8   V8   700    700    700  5.127    0.077   .  0.651    0.722       5.602
#> 9   V9   700    700    700  0.764    0.683      0.839    0.657       2.062
#> 10 V10   700    700    700  0.132    0.936      0.602    0.740       1.780
#> 11 V11   700    700    700  2.998    0.223      0.983    0.612       3.522
#> 12 V12   700    700    700  1.180    0.554      1.807    0.405       2.799
#> 13 V13   700    700    700  0.389    0.823      0.100    0.951       0.536
#> 14 V14   700    700    700  0.426    0.808      1.152    0.562       1.280
#> 15 V15   700    700    700  2.404    0.301      0.201    0.905       2.490
#> 16 V16   700    700    700  2.627    0.269      1.520    0.468       3.136
#> 17 V17   700    700    700  0.910    0.634      2.532    0.282       3.079
#> 18 V18   700    700    700  0.203    0.904      0.588    0.746       3.003
#> 19 V19   700    700    700  0.894    0.639      0.617    0.735       4.711
#> 20 V20   700    700    700  0.124    0.940      0.400    0.819       0.559
#> 21 V21   700    700    700  0.649    0.723      1.411    0.494       1.460
#> 22 V22   700    700    700  0.297    0.862      0.317    0.853       2.652
#> 23 V23   700    700    700  1.537    0.464      1.202    0.548       2.903
#> 24 V24   700    700    700  0.622    0.733      4.768    0.092   .   5.282
#> 25 V25   700    700    700  2.413    0.299      0.856    0.652       3.055
#> 26 V26   700    700    700  1.168    0.558      0.071    0.965       1.205
#> 27 V27   700    700    700  2.722    0.256      1.246    0.536       3.368
#> 28 V28   700    700    700  0.804    0.669      0.943    0.624       2.012
#> 29 V29   700    700    700  5.237    0.073   .  3.084    0.214       5.666
#> 30 V30   700    700    700  3.259    0.196      0.345    0.842       3.680
#> 31 V31   700    700    700  1.849    0.397      2.611    0.271       3.125
#> 32 V32   700    700    700  0.740    0.691      0.789    0.674       0.954
#> 33 V33   700    700    700  0.419    0.811      2.090    0.352       2.380
#> 34 V34   700    700    700  3.110    0.211      5.182    0.075   .   5.436
#> 35 V35   700    700    700  0.313    0.855      4.591    0.101       4.608
#> 36 V36   700    700    700  2.016    0.365      3.152    0.207       3.582
#> 37 V37   700    700    700  5.327    0.070   .  4.923    0.085   .   8.555
#> 38 V38   700    700    700  2.733    0.255      4.922    0.085   .   5.430
#> 39 V39   700    700    700  3.422    0.181      2.458    0.292       5.612
#> 40 V40   700    700    700  3.603    0.165      1.429    0.490       4.573
#>    p.grdifrs    
#> 1      0.000 ***
#> 2      0.000 ***
#> 3      0.876    
#> 4      0.273    
#> 5      0.931    
#> 6      0.759    
#> 7      0.204    
#> 8      0.231    
#> 9      0.724    
#> 10     0.776    
#> 11     0.475    
#> 12     0.592    
#> 13     0.970    
#> 14     0.865    
#> 15     0.646    
#> 16     0.535    
#> 17     0.545    
#> 18     0.557    
#> 19     0.318    
#> 20     0.968    
#> 21     0.834    
#> 22     0.618    
#> 23     0.574    
#> 24     0.260    
#> 25     0.549    
#> 26     0.877    
#> 27     0.498    
#> 28     0.734    
#> 29     0.226    
#> 30     0.451    
#> 31     0.537    
#> 32     0.917    
#> 33     0.666    
#> 34     0.245    
#> 35     0.330    
#> 36     0.466    
#> 37     0.073   .
#> 38     0.246    
#> 39     0.230    
#> 40     0.334    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 2
#>   - GRDIF statistic used for purification: GRDIF(R)
#>   - DIF Items identified by GRDIF(R): 
#>     1, 2 
#>   - GRDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc1 n.foc2 grdifr p.grdifr     grdifs p.grdifs    
#> 1   V1      1   700    700    700 47.438    0.000 *** 35.792    0.000 ***
#> 2   V2      0   700    700    700 67.755    0.000 ***  9.828    0.007  **
#> 3   V3      2   700    700    700  0.071    0.965      0.512    0.774    
#> 4   V4      2   700    700    700  5.967    0.051   .  1.492    0.474    
#> 5   V5      2   700    700    700  1.053    0.591      1.089    0.580    
#> 6   V6      2   700    700    700  0.982    0.612      0.875    0.646    
#> 7   V7      2   700    700    700  0.565    0.754      2.211    0.331    
#> 8   V8      2   700    700    700  3.007    0.222      0.626    0.731    
#> 9   V9      2   700    700    700  1.890    0.389      1.381    0.501    
#> 10 V10      2   700    700    700  0.068    0.966      0.354    0.838    
#> 11 V11      2   700    700    700  2.476    0.290      0.602    0.740    
#> 12 V12      2   700    700    700  0.442    0.802      0.838    0.658    
#> 13 V13      2   700    700    700  0.616    0.735      0.194    0.908    
#> 14 V14      2   700    700    700  0.509    0.775      1.280    0.527    
#> 15 V15      2   700    700    700  1.362    0.506      0.280    0.870    
#> 16 V16      2   700    700    700  1.706    0.426      0.728    0.695    
#> 17 V17      2   700    700    700  0.492    0.782      1.598    0.450    
#> 18 V18      2   700    700    700  0.254    0.881      0.753    0.686    
#> 19 V19      2   700    700    700  0.818    0.664      0.729    0.695    
#> 20 V20      2   700    700    700  0.183    0.912      0.736    0.692    
#> 21 V21      2   700    700    700  1.323    0.516      1.853    0.396    
#> 22 V22      2   700    700    700  0.522    0.770      0.479    0.787    
#> 23 V23      2   700    700    700  0.389    0.823      1.518    0.468    
#> 24 V24      2   700    700    700  0.635    0.728      3.609    0.165    
#> 25 V25      2   700    700    700  1.810    0.404      0.501    0.778    
#> 26 V26      2   700    700    700  0.899    0.638      0.161    0.923    
#> 27 V27      2   700    700    700  3.639    0.162      1.570    0.456    
#> 28 V28      2   700    700    700  1.116    0.572      1.314    0.518    
#> 29 V29      2   700    700    700  5.839    0.054   .  3.780    0.151    
#> 30 V30      2   700    700    700  2.234    0.327      0.276    0.871    
#> 31 V31      2   700    700    700  2.154    0.341      3.130    0.209    
#> 32 V32      2   700    700    700  0.228    0.892      0.417    0.812    
#> 33 V33      2   700    700    700  1.010    0.604      2.913    0.233    
#> 34 V34      2   700    700    700  2.730    0.255      4.404    0.111    
#> 35 V35      2   700    700    700  0.050    0.975      3.983    0.136    
#> 36 V36      2   700    700    700  2.476    0.290      3.534    0.171    
#> 37 V37      2   700    700    700  3.452    0.178      2.577    0.276    
#> 38 V38      2   700    700    700  2.202    0.332      4.404    0.111    
#> 39 V39      2   700    700    700  4.020    0.134      1.500    0.472    
#> 40 V40      2   700    700    700  3.722    0.156      1.209    0.546    
#>    grdifrs p.grdifrs    
#> 1   49.796     0.000 ***
#> 2   82.604     0.000 ***
#> 3    0.718     0.949    
#> 4    6.940     0.139    
#> 5    1.510     0.825    
#> 6    2.229     0.694    
#> 7    6.448     0.168    
#> 8    3.577     0.466    
#> 9    3.267     0.514    
#> 10   1.546     0.818    
#> 11   2.931     0.570    
#> 12   1.258     0.869    
#> 13   0.805     0.938    
#> 14   1.436     0.838    
#> 15   1.492     0.828    
#> 16   1.892     0.756    
#> 17   2.145     0.709    
#> 18   3.628     0.459    
#> 19   4.324     0.364    
#> 20   0.892     0.926    
#> 21   2.075     0.722    
#> 22   2.885     0.577    
#> 23   2.120     0.714    
#> 24   4.736     0.316    
#> 25   2.568     0.632    
#> 26   0.995     0.910    
#> 27   4.427     0.351    
#> 28   2.463     0.651    
#> 29   6.052     0.195    
#> 30   2.686     0.612    
#> 31   3.566     0.468    
#> 32   0.510     0.973    
#> 33   2.961     0.564    
#> 34   4.832     0.305    
#> 35   4.103     0.392    
#> 36   4.184     0.382    
#> 37   5.275     0.260    
#> 38   5.098     0.277    
#> 39   5.849     0.211    
#> 40   4.549     0.337    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05

# DIF items identified after purification
grdif_pur$with_purify$dif_item
#> [1] 1 2

# Post-hoc pairwise RDIF results after purification
grdif_pur$with_purify$post.hoc
#>   id group.pair   rdifr z.rdifr  rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs
#> 1 V2    G1 & G2 -0.1116 -5.1701 0.0267  1.5508 30.2416  0.0000  0.1210   0.0000
#> 2 V2    G1 & G3 -0.1766 -8.0969 0.0495  3.1349 80.9867  0.0000  0.0017   0.0000
#> 3 V2    G2 & G3 -0.0650 -2.9333 0.0228  1.5540 13.4864  0.0034  0.1202   0.0012
#> 4 V1    G1 & G2 -0.0617 -3.5088 0.0566  3.5106 14.0353  0.0005  0.0004   0.0009
#> 5 V1    G1 & G3 -0.1229 -6.8674 0.0887  5.9448 49.4401  0.0000  0.0000   0.0000
#> 6 V1    G2 & G3 -0.0612 -3.2639 0.0320  2.4372 10.7151  0.0011  0.0148   0.0047
#>   n.ref n.foc n.total n_iter
#> 1   700   700    1400      0
#> 2   700   700    1400      0
#> 3   700   700    1400      0
#> 4   700   700    1400      1
#> 5   700   700    1400      1
#> 6   700   700    1400      1
```

The purification, driven by GRDIF$`_R`$, finishes after 2 iterations
(one item is removed per iteration) and again flags items 1 and 2.

------------------------------------------------------------------------

## Part 3: CATSIB (`catsib()`)

### Statistical Framework

CATSIB (Nandakumar & Roussos, 2004) is a modified version of SIBTEST
(Shealy & Stout, 1993) adapted for computerized adaptive testing (CAT)
environments. The procedure estimates the DIF effect size
$`\hat{\beta}`$ by comparing the observed proportions of correct
responses between the reference and focal groups across matched ability
bins.

A key feature of CATSIB is its **regression correction**: because impact
(group mean ability differences) induces stochastic ordering of the two
groups on $`\hat{\theta}`$, naively matching examinees on
$`\hat{\theta}`$ inflates Type I error. To address this, CATSIB first
transforms each examinee’s ability estimate into a regression-corrected
score $`\hat{\theta}^*_G`$, which estimates the conditional expectation
$`E_G[\theta \mid \hat{\theta}]`$ separately for each group $`G`$(Shealy
& Stout, 1993):

``` math
\hat{\theta}^*_G = \bar{\hat{\theta}}_G + \hat{\rho}^2_G\left(\hat{\theta} - \bar{\hat{\theta}}_G\right)
```

where $`\bar{\hat{\theta}}_G`$ is the mean of the ability estimates in
group $`G`$ and
$`\hat{\rho}^2_G = 1 - \hat{\sigma}^2_{e,G} / \hat{\sigma}^2_{\hat{\theta},G}`$
is the estimated reliability in group $`G`$, with
$`\hat{\sigma}^2_{e,G}`$ being the mean squared standard error (SE) of
ability estimates and $`\hat{\sigma}^2_{\hat{\theta},G}`$ being the
observed variance of $`\hat{\theta}`$ in group $`G`$. Estimates at a
limit of `range` are excluded when these three quantities are computed.
Examinees are then matched on $`\hat{\theta}^*`$ rather than
$`\hat{\theta}`$.

The corrected scores $`\hat{\theta}^*`$ are divided into $`K`$
equal-width ability bins. Within each bin $`k`$, the observed
proportions of correct responses, $`\hat{P}_{R,k}`$ and
$`\hat{P}_{F,k}`$, are computed for the reference and focal groups,
respectively. The DIF effect size is estimated as:

``` math
\hat{\beta} = \sum_{k=1}^{K} \left[\hat{P}_{R,k} - \hat{P}_{F,k}\right] \hat{p}_k
```

where $`\hat{p}_k`$ is the weight for bin $`k`$. By default
(`weight.group = "comb"`), $`\hat{p}_k`$ is the observed proportion of
all examinees (both groups combined) classified into bin $`k`$,
following the recommendation of Nandakumar & Roussos (2004).
Alternatively, $`\hat{p}_k`$ can be defined using only the focal group
distribution (`weight.group = "foc"`) or only the reference group
distribution (`weight.group = "ref"`).

A positive $`\hat{\beta}`$ indicates that the item is easier for the
reference group than for the focal group (i.e., DIF favoring the
reference group). Under the null hypothesis of no DIF
($`H_0: \beta = 0`$), the standardized
$`\hat{\beta}/\widehat{SE}(\hat{\beta})`$ asymptotically follows a
standard normal distribution.

> **Note:** CATSIB requires standard errors (SE) of ability estimates in
> addition to point estimates, because $`\hat{\rho}^2_G`$ in the
> regression correction uses the measurement error variance. Supply
> these via the `se` argument. When `score` and `se` are provided
> externally, `x` (item metadata) is not required (`x = NULL`), unless
> purification is applied, in which case `x` must be provided for
> internal ability re-estimation.

### Key arguments of `catsib()`

- `x`: Item metadata data frame. Required when `score = NULL` (ability
  estimated internally) or when `purify = TRUE`. Can be set to `NULL` if
  `score` and `se` are supplied and `purify = FALSE`.
- `se`: Numeric vector of standard errors for the ability estimates
  (required for the regression correction). Obtain via
  `est_score(..., se = TRUE)$se.theta`.
- `n.bin`: A two-element vector `c(max_bins, min_bins)` controlling the
  range for the number of ability-scale intervals (default `c(80, 10)`).
- `min.binsize`: Minimum number of examinees required in each bin for
  both groups; bins failing this criterion are excluded from the
  computation (default `3`).
- `max.del`: Maximum allowable proportion of examinees excluded during
  the binning process (default `0.075`).
- `weight.group`: Target ability distribution used to compute
  $`\hat{p}_k`$: `"comb"` (combined reference + focal; default), `"foc"`
  (focal only), or `"ref"` (reference only).

### Example 1: CATSIB without purification

We use the same 40-item dichotomous data simulated in the Setup section.
Since ability estimates and their SEs are supplied externally and no
purification is applied, `x = NULL` is used.

``` r

# Estimate ability with SE (required for the regression correction)
score_se <- est_score(
  x      = meta_pool,
  data   = resp_pool,
  D      = 1.702,
  method = "ML",
  range  = c(-5, 5),
  se     = TRUE
)

# x = NULL: item metadata not required when score/SE are provided externally
# and purify = FALSE
catsib_npur <- catsib(
  x            = NULL,
  data         = resp_pool,
  score        = score_se$est.theta,
  se           = score_se$se.theta,   # SE required for regression correction
  group        = group_vec,
  focal.name   = 1,
  weight.group = "comb",
  D            = 1.702,
  alpha        = 0.05,
  purify       = FALSE,
  verbose      = FALSE
)

# Summary output
print(catsib_npur)
#> 
#> Call:
#> catsib(x = NULL, data = resp_pool, score = score_se$est.theta, 
#>     se = score_se$se.theta, group = group_vec, focal.name = 1, 
#>     D = 1.702, weight.group = "comb", alpha = 0.05, purify = FALSE, 
#>     verbose = FALSE)
#> 
#> DIF analysis using CATSIB method 
#> 
#>  1. Without purification 
#> 
#>   - Potential DIF Items: 
#>     1, 2, 5, 6, 8, 9, 19, 26, 32, 33, 35, 37, 39 
#>   - Test Statistic: 
#> 
#>         id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   item.1   939   973    1912  0.078 0.014  5.776 0.000 ***
#> 2   item.2   939   973    1912  0.112 0.018  6.070 0.000 ***
#> 3   item.3   939   973    1912 -0.033 0.020 -1.696 0.090   .
#> 4   item.4   939   973    1912  0.030 0.018  1.654 0.098   .
#> 5   item.5   939   973    1912  0.148 0.013 11.353 0.000 ***
#> 6   item.6   939   973    1912  0.126 0.017  7.380 0.000 ***
#> 7   item.7   939   973    1912 -0.008 0.014 -0.555 0.579    
#> 8   item.8   939   973    1912 -0.069 0.018 -3.935 0.000 ***
#> 9   item.9   939   973    1912 -0.029 0.015 -1.968 0.049   *
#> 10 item.10   939   973    1912 -0.014 0.019 -0.738 0.461    
#> 11 item.11   939   973    1912  0.007 0.006  1.196 0.232    
#> 12 item.12   939   973    1912  0.005 0.019  0.280 0.780    
#> 13 item.13   939   973    1912 -0.026 0.019 -1.406 0.160    
#> 14 item.14   939   973    1912 -0.008 0.011 -0.731 0.465    
#> 15 item.15   939   973    1912 -0.039 0.021 -1.839 0.066   .
#> 16 item.16   939   973    1912  0.021 0.011  1.885 0.060   .
#> 17 item.17   939   973    1912 -0.001 0.012 -0.086 0.932    
#> 18 item.18   939   973    1912 -0.039 0.021 -1.868 0.062   .
#> 19 item.19   939   973    1912 -0.046 0.018 -2.556 0.011   *
#> 20 item.20   939   973    1912 -0.019 0.016 -1.199 0.231    
#> 21 item.21   939   973    1912 -0.024 0.018 -1.306 0.192    
#> 22 item.22   939   973    1912  0.003 0.010  0.280 0.779    
#> 23 item.23   939   973    1912 -0.029 0.017 -1.683 0.092   .
#> 24 item.24   939   973    1912 -0.013 0.013 -0.961 0.337    
#> 25 item.25   939   973    1912 -0.020 0.019 -1.046 0.295    
#> 26 item.26   939   973    1912 -0.040 0.018 -2.200 0.028   *
#> 27 item.27   939   973    1912 -0.031 0.019 -1.650 0.099   .
#> 28 item.28   939   973    1912 -0.021 0.018 -1.142 0.254    
#> 29 item.29   939   973    1912 -0.011 0.020 -0.565 0.572    
#> 30 item.30   939   973    1912 -0.012 0.020 -0.596 0.551    
#> 31 item.31   939   973    1912 -0.024 0.018 -1.300 0.194    
#> 32 item.32   939   973    1912 -0.071 0.018 -3.871 0.000 ***
#> 33 item.33   939   973    1912 -0.046 0.019 -2.446 0.014   *
#> 34 item.34   939   973    1912  0.009 0.005  1.739 0.082   .
#> 35 item.35   939   973    1912 -0.042 0.021 -2.012 0.044   *
#> 36 item.36   939   973    1912 -0.006 0.014 -0.416 0.678    
#> 37 item.37   939   973    1912 -0.034 0.016 -2.087 0.037   *
#> 38 item.38   939   973    1912 -0.037 0.021 -1.778 0.075   .
#> 39 item.39   939   973    1912 -0.038 0.018 -2.116 0.034   *
#> 40 item.40   939   973    1912  0.000 0.012 -0.018 0.986    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Purification was not implemented.

# Items flagged as DIF
catsib_npur$no_purify$dif_item
#>  [1]  1  2  5  6  8  9 19 26 32 33 35 37 39
```

The `dif_stat` data frame contains the estimated $`\hat{\beta}`$
statistic, its standard error (`se`), the standardized $`\hat{\beta}`$
(`z.beta`), and the corresponding p-value for each item. Items with
$`|\hat{\beta}|`$ significantly different from zero (at the `alpha`
level) are flagged as DIF.

In this simulation, $`\hat{\beta}`$ is positive and significant for the
items with a difficulty shift against the focal group (items 1, 2, 5,
and 6), meaning that the items are easier for the reference group. The
items with only a lower focal discrimination (items 3 and 4) are not
flagged, which is consistent with the limited sensitivity of CATSIB to
nonuniform DIF discussed below. CATSIB also flags several DIF-free items
(for example, items 8, 9, 19, 26, and 32), all of them with small
negative $`\hat{\beta}`$ values.

### Example 2: CATSIB with purification

When purification is applied,
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md)
re-estimates ability parameters internally at each iteration using the
purified item set. This requires `x` to be provided. Although Nandakumar
and Roussos (2004) did not originally propose a purification procedure
for CATSIB,
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md)
implements one following the iterative scheme of Lim et al. (2022).

``` r

# x must be provided when purify = TRUE (needed for internal re-scoring)
catsib_pur <- catsib(
  x            = meta_pool,
  data         = resp_pool,
  score        = score_se$est.theta,
  se           = score_se$se.theta,
  group        = group_vec,
  focal.name   = 1,
  weight.group = "comb",
  D            = 1.702,
  alpha        = 0.05,
  purify       = TRUE,
  max.iter     = 20,
  method       = "ML",         # re-estimate abilities with ML at each iteration
  range        = c(-5, 5),
  verbose      = FALSE
)

# Summary output
print(catsib_pur)
#> 
#> Call:
#> catsib(x = meta_pool, data = resp_pool, score = score_se$est.theta, 
#>     se = score_se$se.theta, group = group_vec, focal.name = 1, 
#>     D = 1.702, weight.group = "comb", alpha = 0.05, purify = TRUE, 
#>     max.iter = 20, method = "ML", range = c(-5, 5), verbose = FALSE)
#> 
#> DIF analysis using CATSIB method 
#> 
#>  1. Without purification 
#> 
#>   - Potential DIF Items: 
#>     1, 2, 5, 6, 8, 9, 19, 26, 32, 33, 35, 37, 39 
#>   - Test Statistic: 
#> 
#>     id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   V1   939   973    1912  0.078 0.014  5.776 0.000 ***
#> 2   V2   939   973    1912  0.112 0.018  6.070 0.000 ***
#> 3   V3   939   973    1912 -0.033 0.020 -1.696 0.090   .
#> 4   V4   939   973    1912  0.030 0.018  1.654 0.098   .
#> 5   V5   939   973    1912  0.148 0.013 11.353 0.000 ***
#> 6   V6   939   973    1912  0.126 0.017  7.380 0.000 ***
#> 7   V7   939   973    1912 -0.008 0.014 -0.555 0.579    
#> 8   V8   939   973    1912 -0.069 0.018 -3.935 0.000 ***
#> 9   V9   939   973    1912 -0.029 0.015 -1.968 0.049   *
#> 10 V10   939   973    1912 -0.014 0.019 -0.738 0.461    
#> 11 V11   939   973    1912  0.007 0.006  1.196 0.232    
#> 12 V12   939   973    1912  0.005 0.019  0.280 0.780    
#> 13 V13   939   973    1912 -0.026 0.019 -1.406 0.160    
#> 14 V14   939   973    1912 -0.008 0.011 -0.731 0.465    
#> 15 V15   939   973    1912 -0.039 0.021 -1.839 0.066   .
#> 16 V16   939   973    1912  0.021 0.011  1.885 0.060   .
#> 17 V17   939   973    1912 -0.001 0.012 -0.086 0.932    
#> 18 V18   939   973    1912 -0.039 0.021 -1.868 0.062   .
#> 19 V19   939   973    1912 -0.046 0.018 -2.556 0.011   *
#> 20 V20   939   973    1912 -0.019 0.016 -1.199 0.231    
#> 21 V21   939   973    1912 -0.024 0.018 -1.306 0.192    
#> 22 V22   939   973    1912  0.003 0.010  0.280 0.779    
#> 23 V23   939   973    1912 -0.029 0.017 -1.683 0.092   .
#> 24 V24   939   973    1912 -0.013 0.013 -0.961 0.337    
#> 25 V25   939   973    1912 -0.020 0.019 -1.046 0.295    
#> 26 V26   939   973    1912 -0.040 0.018 -2.200 0.028   *
#> 27 V27   939   973    1912 -0.031 0.019 -1.650 0.099   .
#> 28 V28   939   973    1912 -0.021 0.018 -1.142 0.254    
#> 29 V29   939   973    1912 -0.011 0.020 -0.565 0.572    
#> 30 V30   939   973    1912 -0.012 0.020 -0.596 0.551    
#> 31 V31   939   973    1912 -0.024 0.018 -1.300 0.194    
#> 32 V32   939   973    1912 -0.071 0.018 -3.871 0.000 ***
#> 33 V33   939   973    1912 -0.046 0.019 -2.446 0.014   *
#> 34 V34   939   973    1912  0.009 0.005  1.739 0.082   .
#> 35 V35   939   973    1912 -0.042 0.021 -2.012 0.044   *
#> 36 V36   939   973    1912 -0.006 0.014 -0.416 0.678    
#> 37 V37   939   973    1912 -0.034 0.016 -2.087 0.037   *
#> 38 V38   939   973    1912 -0.037 0.021 -1.778 0.075   .
#> 39 V39   939   973    1912 -0.038 0.018 -2.116 0.034   *
#> 40 V40   939   973    1912  0.000 0.012 -0.018 0.986    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 12
#>   - Potential DIF Items: 
#>     1, 2, 4, 5, 6, 8, 16, 19, 21, 32, 33, 34 
#>   - Test Statistic: 
#> 
#>     id n.iter n.ref n.foc n.total   beta    se  z.beta     p    
#> 1   V1      3   950   968    1918  0.083 0.014   6.127 0.000 ***
#> 2   V2      2   950   969    1919  0.125 0.018   6.824 0.000 ***
#> 3   V3     12   939   973    1912 -0.010 0.020  -0.527 0.598    
#> 4   V4      4   950   975    1925  0.051 0.018   2.844 0.004  **
#> 5   V5      0   939   973    1912  0.148 0.013  11.353 0.000 ***
#> 6   V6      1   925   988    1913  0.134 0.017   7.992 0.000 ***
#> 7   V7     12   939   973    1912  0.003 0.014   0.244 0.807    
#> 8   V8     11   940   971    1911 -0.045 0.018  -2.555 0.011   *
#> 9   V9     12   939   973    1912 -0.021 0.015  -1.430 0.153    
#> 10 V10     12   939   973    1912 -0.006 0.018  -0.340 0.734    
#> 11 V11     12   939   973    1912  0.007 0.006   1.187 0.235    
#> 12 V12     12   939   973    1912  0.027 0.018   1.512 0.130    
#> 13 V13     12   939   973    1912 -0.013 0.018  -0.674 0.501    
#> 14 V14     12   939   973    1912 -0.004 0.011  -0.347 0.729    
#> 15 V15     12   939   973    1912 -0.018 0.021  -0.876 0.381    
#> 16 V16      9   946   975    1921  0.038 0.011   3.468 0.000 ***
#> 17 V17     12   939   973    1912  0.009 0.012   0.803 0.422    
#> 18 V18     12   939   973    1912 -0.029 0.021  -1.398 0.162    
#> 19 V19     10   959   968    1927 -0.044 0.018  -2.473 0.013   *
#> 20 V20     12   939   973    1912  0.002 0.016   0.109 0.913    
#> 21 V21      6   656  1000    1656 -0.169 0.021  -8.061 0.000 ***
#> 22 V22     12   939   973    1912  0.007 0.009   0.785 0.433    
#> 23 V23     12   939   973    1912 -0.010 0.017  -0.627 0.531    
#> 24 V24     12   939   973    1912 -0.002 0.013  -0.156 0.876    
#> 25 V25     12   939   973    1912 -0.012 0.018  -0.661 0.509    
#> 26 V26     12   939   973    1912 -0.013 0.018  -0.733 0.464    
#> 27 V27     12   939   973    1912 -0.025 0.018  -1.354 0.176    
#> 28 V28     12   939   973    1912 -0.021 0.018  -1.202 0.229    
#> 29 V29     12   939   973    1912  0.003 0.019   0.146 0.884    
#> 30 V30     12   939   973    1912  0.004 0.020   0.199 0.842    
#> 31 V31     12   939   973    1912 -0.014 0.018  -0.795 0.427    
#> 32 V32      5   662  1000    1662 -0.216 0.021 -10.378 0.000 ***
#> 33 V33      7   648  1000    1648 -0.238 0.020 -11.796 0.000 ***
#> 34 V34      8   747  1000    1747  0.080 0.006  14.342 0.000 ***
#> 35 V35     12   939   973    1912 -0.025 0.021  -1.197 0.231    
#> 36 V36     12   939   973    1912  0.007 0.014   0.541 0.589    
#> 37 V37     12   939   973    1912 -0.012 0.016  -0.746 0.455    
#> 38 V38     12   939   973    1912 -0.017 0.020  -0.824 0.410    
#> 39 V39     12   939   973    1912 -0.015 0.018  -0.848 0.396    
#> 40 V40     12   939   973    1912  0.005 0.012   0.395 0.693    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05

# DIF items identified after purification
catsib_pur$with_purify$dif_item
#>  [1]  1  2  4  5  6  8 16 19 21 32 33 34
```

In this run, the purification finishes after 12 iterations. Of the six
items with injected DIF, it flags items 1, 2, 4, 5, 6 (without
purification: items 1, 2, 5, 6) and misses item 3. The number of flagged
DIF-free items changes from 9 (items 8, 9, 19, 26, 32, 33, 35, 37, 39)
to 7 (items 8, 16, 19, 21, 32, 33, 34); items 16, 21, 34 are new false
positives that did not appear without purification. Purification
therefore changes the CATSIB results but does not remove all false
positives, so they should be interpreted with caution. Depending on the
platform,
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md) may also
warn during purification that the estimated reliability of a group fell
below 0.05 and was floored (see the considerations below); the warning
is displayed above when it occurs.

### Important Considerations for CATSIB

**Sensitivity to DIF type.** CATSIB, like its predecessor SIBTEST
(Shealy & Stout, 1993), was originally designed and validated for
detecting **uniform DIF**, the condition in which the direction of group
differences in item performance is consistent across all ability levels.
The $`\hat{\beta}`$ statistic accumulates the signed difference
$`\hat{P}_{R,k} - \hat{P}_{F,k}`$ across all ability bins, which means
that positive and negative bin-level differences cancel each other out.
As a result, CATSIB has limited sensitivity to **nonuniform DIF** (where
item discrimination differs between groups) or **mixed DIF** (where both
difficulty and discrimination differ), because in these cases the
bin-level differences may partially cancel and yield
$`\hat{\beta} \approx 0`$ even when DIF is present. If nonuniform or
mixed DIF is suspected,
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) with
RDIF$`_S`$ or RDIF$`_{RS}`$ is recommended instead.

**Numerical instability of the regression correction when the
reliability is low.** The regression correction (Equation 7 of
Nandakumar & Roussos (2004)) computes, for each group $`G`$, the
reliability estimate:

``` math
\hat{\rho}^2_G = 1 - \frac{\hat{\sigma}^2_{e,G}}{\hat{\sigma}^2_{\hat{\theta},G}}
```

where $`\hat{\sigma}^2_{e,G}`$ is the mean squared SE of ability
estimates and $`\hat{\sigma}^2_{\hat{\theta},G}`$ is the observed
variance of $`\hat{\theta}`$ in group $`G`$. The reliability estimate
falls toward or below zero when the mean squared SE approaches the
observed variance of $`\hat{\theta}`$, for example when few items remain
after purification. A negative $`\hat{\rho}^2_G`$ would reflect the
scores around the group mean and reverse their order, which severely
distorts the bin-level matching and leads to inflated Type I error
rates.

To prevent this,
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md) clamps
$`\hat{\rho}^2_G`$ to the interval $`[0.05,\,
1]`$:

``` math
\hat{\rho}^2_G \leftarrow \max\!\left(0.05,\; \min\!\left(1,\; 1 -
\frac{\hat{\sigma}^2_{e,G}}{\hat{\sigma}^2_{\hat{\theta},G}}\right)\right)
```

When the unclamped value falls below 0.05, a warning is issued. This
situation is most likely to occur during purification, when the removal
of DIF-flagged items reduces the number of items used for ability
re-estimation, causing standard errors to grow. If a warning is
triggered, results should be interpreted with caution, or
`purify = FALSE` should be considered.

------------------------------------------------------------------------

## Summary of Function Inputs

All three functions share a common calling pattern:

``` r
func(
  x          = <item metadata>,     # pooled item parameters (data frame)
  data       = <response matrix>,   # pooled responses (examinees x items)
  score      = <ability vector>,    # pooled ability estimates
  se         = <SE vector>,         # catsib() only: SEs of the ability estimates (required when score is given)
  group      = <group vector>,      # group membership labels
  focal.name = <focal label(s)>,    # which group(s) are focal
  D          = 1.702,               # scaling constant (match calibration)
  alpha      = 0.05,                # significance level
  purify     = TRUE / FALSE,        # apply iterative purification?
  verbose    = FALSE
)
```

Output always contains `$no_purify$dif_stat` (table of statistics for
all items) and `$no_purify$dif_item` (indices of flagged items). For
[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) and
[`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md),
`no_purify$dif_item` is a list of three vectors, one per statistic; for
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md) it is a
single vector. When `purify = TRUE`, corresponding results are also
available under `$with_purify`, where `dif_item` is a single vector: the
items flagged by the `purify.by` statistic across all purification
iterations (for
[`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md), by the
CATSIB statistic).

------------------------------------------------------------------------

## References

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement*, *59*(1), 80–104.
<https://doi.org/10.1111/jedm.12313>

Lim, H., Zhu, D., Choe, E. M., & Han, K. T. (2024). Detecting
differential item functioning among multiple groups using IRT residual
DIF framework. *Journal of Educational Measurement*, *61*(4), 656–681.
<https://doi.org/10.1111/jedm.12415>

Nandakumar, R., & Roussos, L. (2004). Evaluation of the CATSIB DIF
procedure in a pretest setting. *Journal of Educational and Behavioral
Statistics*, *29*(2), 177–199.
<https://doi.org/10.3102/10769986029002177>

Shealy, R. T., & Stout, W. F. (1993). A model-based standardization
approach that separates true bias/DIF from group ability differences and
detects test bias/DTF as well as item bias/DIF. *Psychometrika*,
*58*(2), 159–194. <https://doi.org/10.1007/BF02294572>
