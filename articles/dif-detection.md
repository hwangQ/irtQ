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

- **Items 1–2**: Uniform DIF — the focal group has a higher difficulty
  parameter ($`b + 0.5`$) while the discrimination parameter is
  unchanged.
- **Items 3–4**: Nonuniform DIF — the focal group has a lower
  discrimination parameter ($`a \times 0.5`$) while the difficulty
  parameter is unchanged.
- **Items 5–6**: Mixed DIF — the focal group has both a lower
  discrimination parameter ($`a \times 0.5`$) and a higher difficulty
  parameter ($`b + 0.5`$), making these items both less discriminating
  and harder for the focal group.
- **Items 7–40**: No DIF.

The focal group has a slightly lower mean ability ($`\mu = -0.3`$) than
the reference group ($`\mu = 0`$), representing a realistic impact
condition.

``` r

J   <- 40     # number of items
N_R <- 1000   # reference group size
N_F <- 1000   # focal group size

# Base item parameters (no DIF) shared by both groups
set.seed(100)
a_base <- round(runif(J, 0.7, 1.7), 2)
b_base <- round(runif(J, -4.0, 4.0), 2)
g_base <- rep(0.15, J)

# Reference group item metadata (no DIF)
meta_ref <- shape_df(
  par.drm = list(a = a_base, b = b_base, g = g_base),
  cats    = 2,
  model   = "3PLM"
)

# Focal group item metadata: DIF injected into items 1–6
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
membership (Lim et al. 2022).

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
  verbose    = FALSE
)

meta_pool <- mod_pool$par.est

# Estimate ability using ML method
# ML is recommended for RDIF analysis as it provides point estimates
# unaffected by prior distribution assumptions [@lim_etal2022]
score_pool <- est_score(
  x      = meta_pool,
  data   = resp_pool,
  D      = 1.702,
  method = "ML",
  range  = c(-5, 5)
)$est.theta
```

------------------------------------------------------------------------

## Part 1: RDIF — Residual-Based DIF (`rdif()`)

### Statistical Framework

The RDIF framework (Lim et al. 2022) detects DIF by comparing item-level
residuals—defined as the difference between observed and model-expected
item scores—between the reference and focal groups. For a dichotomously
scored item, the residual for examinee $`h`$ is:

``` math
r_h = x_h - P_h(\hat{\theta}_h)
```

where $`x_h \in \{0, 1\}`$ is the observed response and
$`P_h(\hat{\theta}_h)`$ is the IRT model-predicted probability of a
correct response. Three statistics are computed from these residuals:

**RDIF$`_R`$** — targets **uniform DIF** via differences in mean raw
residuals:

``` math
\text{RDIF}_R = \frac{1}{N_F}\sum_{j=1}^{N_F} r_{Fj} - \frac{1}{N_R}\sum_{i=1}^{N_R} r_{Ri}
```

Under the null hypothesis of no DIF, $`\text{RDIF}_R`$ asymptotically
follows a normal distribution
$`\mathcal{N}(\mu_{\text{RDIF}_R},\, \sigma^2_{\text{RDIF}_R})`$.

**RDIF$`_S`$** — targets **nonuniform DIF** via differences in mean
squared residuals:

``` math
\text{RDIF}_S = \frac{1}{N_F}\sum_{j=1}^{N_F} r_{Fj}^2 - \frac{1}{N_R}\sum_{i=1}^{N_R} r_{Ri}^2
```

Under the null hypothesis, $`\text{RDIF}_S`$ also asymptotically follows
a normal distribution
$`\mathcal{N}(\mu_{\text{RDIF}_S},\, \sigma^2_{\text{RDIF}_S})`$.

**RDIF$`_{RS}`$** — a joint Wald-type statistic that detects **both
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
> DIF simultaneously (Lim et al. 2022). Use $`\text{RDIF}_R`$ and
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
#>     1, 3, 4, 6, 25 
#>   - DIF Items identified by RDIF(S): 
#>     1, 3, 4, 5, 6, 18 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 3, 4, 5, 6 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs    
#> 1   V1  1000  1000 -0.106   0.000 ***  0.060   0.000 *** 44.231    0.000 ***
#> 2   V2  1000  1000  0.010   0.555      0.001   0.485      1.296    0.523    
#> 3   V3  1000  1000  0.051   0.004  **  0.021   0.004  **  8.536    0.014   *
#> 4   V4  1000  1000  0.141   0.000 ***  0.064   0.000 *** 56.025    0.000 ***
#> 5   V5  1000  1000  0.010   0.598      0.021   0.007  **  9.329    0.009  **
#> 6   V6  1000  1000 -0.073   0.000 ***  0.058   0.000 *** 47.996    0.000 ***
#> 7   V7  1000  1000  0.001   0.968      0.005   0.747      1.000    0.607    
#> 8   V8  1000  1000  0.026   0.118      0.013   0.101      2.702    0.259    
#> 9   V9  1000  1000 -0.008   0.296      0.016   0.412      1.149    0.563    
#> 10 V10  1000  1000 -0.012   0.418      0.033   0.386      0.808    0.668    
#> 11 V11  1000  1000  0.012   0.366      0.022   0.682      0.870    0.647    
#> 12 V12  1000  1000 -0.001   0.899      0.010   0.926      0.016    0.992    
#> 13 V13  1000  1000  0.021   0.059   .  0.009   0.135      3.625    0.163    
#> 14 V14  1000  1000  0.000   0.986      0.027   0.503      1.384    0.501    
#> 15 V15  1000  1000 -0.010   0.586     -0.006   0.639      0.338    0.845    
#> 16 V16  1000  1000  0.014   0.139      0.013   0.406      2.301    0.316    
#> 17 V17  1000  1000  0.000   0.971      0.009   0.818      0.537    0.765    
#> 18 V18  1000  1000  0.012   0.242      0.006   0.049   *  4.255    0.119    
#> 19 V19  1000  1000  0.009   0.624     -0.008   0.686      0.766    0.682    
#> 20 V20  1000  1000  0.002   0.786      0.009   0.712      0.144    0.930    
#> 21 V21  1000  1000 -0.020   0.261      0.018   0.917      1.305    0.521    
#> 22 V22  1000  1000  0.023   0.191      0.004   0.323      1.710    0.425    
#> 23 V23  1000  1000 -0.003   0.861     -0.002   0.915      0.360    0.835    
#> 24 V24  1000  1000  0.030   0.092   .  0.008   0.164      2.855    0.240    
#> 25 V25  1000  1000  0.038   0.031   *  0.020   0.831      4.738    0.094   .
#> 26 V26  1000  1000  0.004   0.825      0.019   0.390      1.042    0.594    
#> 27 V27  1000  1000  0.027   0.104      0.011   0.258      3.531    0.171    
#> 28 V28  1000  1000  0.016   0.340      0.031   0.269      2.377    0.305    
#> 29 V29  1000  1000 -0.013   0.127      0.017   0.737      4.118    0.128    
#> 30 V30  1000  1000 -0.001   0.939     -0.006   0.544      1.623    0.444    
#> 31 V31  1000  1000 -0.024   0.153      0.036   0.197      2.586    0.274    
#> 32 V32  1000  1000  0.005   0.679      0.027   0.757      0.569    0.752    
#> 33 V33  1000  1000  0.029   0.123      0.009   0.342      2.527    0.283    
#> 34 V34  1000  1000  0.005   0.750      0.004   0.750      0.101    0.951    
#> 35 V35  1000  1000  0.020   0.255     -0.003   0.370      1.296    0.523    
#> 36 V36  1000  1000  0.001   0.968     -0.007   0.667      0.367    0.832    
#> 37 V37  1000  1000  0.022   0.204      0.003   0.441      5.755    0.056   .
#> 38 V38  1000  1000  0.012   0.495      0.000   0.534      0.473    0.789    
#> 39 V39  1000  1000  0.000   0.998     -0.009   0.874      0.420    0.811    
#> 40 V40  1000  1000 -0.004   0.535      0.009   0.610      0.549    0.760    
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
#>     id   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs n.ref
#> 1   V1 -0.1057 -6.6338  0.0602  4.1145 44.2310  0.0000  0.0000   0.0000  1000
#> 2   V2  0.0098  0.5896  0.0012  0.6985  1.2964  0.5554  0.4848   0.5230  1000
#> 3   V3  0.0511  2.8434  0.0211  2.8962  8.5364  0.0045  0.0038   0.0140  1000
#> 4   V4  0.1413  7.4738  0.0645  7.3466 56.0249  0.0000  0.0000   0.0000  1000
#> 5   V5  0.0101  0.5268  0.0214  2.6853  9.3288  0.5983  0.0072   0.0094  1000
#> 6   V6 -0.0726 -3.6548  0.0579  6.0024 47.9960  0.0003  0.0000   0.0000  1000
#> 7   V7  0.0006  0.0397  0.0049  0.3232  1.0002  0.9683  0.7465   0.6065  1000
#> 8   V8  0.0257  1.5616  0.0129  1.6419  2.7024  0.1184  0.1006   0.2589  1000
#> 9   V9 -0.0083 -1.0458  0.0163  0.8195  1.1488  0.2957  0.4125   0.5630  1000
#> 10 V10 -0.0119 -0.8093  0.0333  0.8675  0.8076  0.4183  0.3856   0.6678  1000
#> 11 V11  0.0119  0.9037  0.0225 -0.4097  0.8702  0.3661  0.6821   0.6472  1000
#> 12 V12 -0.0008 -0.1269  0.0098  0.0923  0.0163  0.8990  0.9265   0.9919  1000
#> 13 V13  0.0214  1.8884  0.0088 -1.4948  3.6250  0.0590  0.1350   0.1632  1000
#> 14 V14 -0.0002 -0.0179  0.0269  0.6692  1.3839  0.9857  0.5033   0.5006  1000
#> 15 V15 -0.0098 -0.5447 -0.0063 -0.4695  0.3381  0.5859  0.6387   0.8445  1000
#> 16 V16  0.0139  1.4789  0.0132 -0.8300  2.3012  0.1392  0.4065   0.3164  1000
#> 17 V17 -0.0003 -0.0362  0.0093  0.2300  0.5368  0.9711  0.8181   0.7646  1000
#> 18 V18  0.0118  1.1704  0.0061 -1.9681  4.2549  0.2418  0.0491   0.1191  1000
#> 19 V19  0.0094  0.4901 -0.0076 -0.4039  0.7664  0.6241  0.6863   0.6817  1000
#> 20 V20  0.0020  0.2710  0.0089 -0.3696  0.1443  0.7864  0.7117   0.9304  1000
#> 21 V21 -0.0197 -1.1232  0.0179 -0.1046  1.3052  0.2614  0.9167   0.5207  1000
#> 22 V22  0.0234  1.3070  0.0035  0.9891  1.7098  0.1912  0.3226   0.4253  1000
#> 23 V23 -0.0028 -0.1747 -0.0024 -0.1062  0.3601  0.8613  0.9154   0.8352  1000
#> 24 V24  0.0301  1.6860  0.0079  1.3907  2.8554  0.0918  0.1643   0.2399  1000
#> 25 V25  0.0380  2.1627  0.0198 -0.2130  4.7376  0.0306  0.8313   0.0936  1000
#> 26 V26  0.0035  0.2212  0.0191 -0.8594  1.0418  0.8249  0.3901   0.5940  1000
#> 27 V27  0.0272  1.6236  0.0105 -1.1300  3.5308  0.1045  0.2585   0.1711  1000
#> 28 V28  0.0159  0.9532  0.0314  1.1062  2.3766  0.3405  0.2687   0.3047  1000
#> 29 V29 -0.0133 -1.5276  0.0174  0.3354  4.1182  0.1266  0.7373   0.1276  1000
#> 30 V30 -0.0014 -0.0762 -0.0058  0.6068  1.6234  0.9393  0.5440   0.4441  1000
#> 31 V31 -0.0241 -1.4284  0.0361  1.2901  2.5857  0.1532  0.1970   0.2745  1000
#> 32 V32  0.0051  0.4134  0.0267  0.3092  0.5686  0.6793  0.7571   0.7525  1000
#> 33 V33  0.0295  1.5423  0.0085  0.9498  2.5272  0.1230  0.3422   0.2826  1000
#> 34 V34  0.0050  0.3187  0.0036  0.3187  0.1007  0.7499  0.7499   0.9509  1000
#> 35 V35  0.0199  1.1381 -0.0032  0.8975  1.2958  0.2551  0.3695   0.5231  1000
#> 36 V36  0.0007  0.0398 -0.0075 -0.4307  0.3669  0.9683  0.6667   0.8324  1000
#> 37 V37  0.0216  1.2710  0.0030  0.7698  5.7549  0.2037  0.4414   0.0563  1000
#> 38 V38  0.0117  0.6827  0.0005  0.6222  0.4730  0.4948  0.5338   0.7894  1000
#> 39 V39  0.0000  0.0022 -0.0089 -0.1587  0.4197  0.9982  0.8739   0.8107  1000
#> 40 V40 -0.0039 -0.6200  0.0085  0.5100  0.5490  0.5353  0.6101   0.7600  1000
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
#> [1]  1  3  4  6 25
#> 
#> $rdifs
#> [1]  1  3  4  5  6 18
#> 
#> $rdifrs
#> [1] 1 3 4 5 6
```

The output contains a `dif_stat` data frame with the following columns
for each item: `rdifr` (RDIF$`_R`$ statistic), `z.rdifr` (standardized
RDIF$`_R`$), `rdifs` (RDIF$`_S`$ statistic), `z.rdifs` (standardized
RDIF$`_S`$), `rdifrs` (RDIF$`_{RS}`$ statistic), `p.rdifr`, `p.rdifs`,
`p.rdifrs` (corresponding p-values), and group sample sizes. Items with
p-values below `alpha` are flagged as DIF.

### Example 2: RDIF with purification

When DIF items are present in the test, they can contaminate the ability
estimates used to compute the RDIF statistics, inflating Type I error
rates for DIF-free items. Iterative purification addresses this by
progressively removing flagged DIF items from ability re-estimation
until the set of flagged items stabilizes (Lim et al. 2022).

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
#>     1, 3, 4, 6, 25 
#>   - DIF Items identified by RDIF(S): 
#>     1, 3, 4, 5, 6, 18 
#>   - DIF Items identified by RDIF(RS): 
#>     1, 3, 4, 5, 6 
#>   - RDIF Statistics: 
#> 
#>     id n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs    
#> 1   V1  1000  1000 -0.106   0.000 ***  0.060   0.000 *** 44.231    0.000 ***
#> 2   V2  1000  1000  0.010   0.555      0.001   0.485      1.296    0.523    
#> 3   V3  1000  1000  0.051   0.004  **  0.021   0.004  **  8.536    0.014   *
#> 4   V4  1000  1000  0.141   0.000 ***  0.064   0.000 *** 56.025    0.000 ***
#> 5   V5  1000  1000  0.010   0.598      0.021   0.007  **  9.329    0.009  **
#> 6   V6  1000  1000 -0.073   0.000 ***  0.058   0.000 *** 47.996    0.000 ***
#> 7   V7  1000  1000  0.001   0.968      0.005   0.747      1.000    0.607    
#> 8   V8  1000  1000  0.026   0.118      0.013   0.101      2.702    0.259    
#> 9   V9  1000  1000 -0.008   0.296      0.016   0.412      1.149    0.563    
#> 10 V10  1000  1000 -0.012   0.418      0.033   0.386      0.808    0.668    
#> 11 V11  1000  1000  0.012   0.366      0.022   0.682      0.870    0.647    
#> 12 V12  1000  1000 -0.001   0.899      0.010   0.926      0.016    0.992    
#> 13 V13  1000  1000  0.021   0.059   .  0.009   0.135      3.625    0.163    
#> 14 V14  1000  1000  0.000   0.986      0.027   0.503      1.384    0.501    
#> 15 V15  1000  1000 -0.010   0.586     -0.006   0.639      0.338    0.845    
#> 16 V16  1000  1000  0.014   0.139      0.013   0.406      2.301    0.316    
#> 17 V17  1000  1000  0.000   0.971      0.009   0.818      0.537    0.765    
#> 18 V18  1000  1000  0.012   0.242      0.006   0.049   *  4.255    0.119    
#> 19 V19  1000  1000  0.009   0.624     -0.008   0.686      0.766    0.682    
#> 20 V20  1000  1000  0.002   0.786      0.009   0.712      0.144    0.930    
#> 21 V21  1000  1000 -0.020   0.261      0.018   0.917      1.305    0.521    
#> 22 V22  1000  1000  0.023   0.191      0.004   0.323      1.710    0.425    
#> 23 V23  1000  1000 -0.003   0.861     -0.002   0.915      0.360    0.835    
#> 24 V24  1000  1000  0.030   0.092   .  0.008   0.164      2.855    0.240    
#> 25 V25  1000  1000  0.038   0.031   *  0.020   0.831      4.738    0.094   .
#> 26 V26  1000  1000  0.004   0.825      0.019   0.390      1.042    0.594    
#> 27 V27  1000  1000  0.027   0.104      0.011   0.258      3.531    0.171    
#> 28 V28  1000  1000  0.016   0.340      0.031   0.269      2.377    0.305    
#> 29 V29  1000  1000 -0.013   0.127      0.017   0.737      4.118    0.128    
#> 30 V30  1000  1000 -0.001   0.939     -0.006   0.544      1.623    0.444    
#> 31 V31  1000  1000 -0.024   0.153      0.036   0.197      2.586    0.274    
#> 32 V32  1000  1000  0.005   0.679      0.027   0.757      0.569    0.752    
#> 33 V33  1000  1000  0.029   0.123      0.009   0.342      2.527    0.283    
#> 34 V34  1000  1000  0.005   0.750      0.004   0.750      0.101    0.951    
#> 35 V35  1000  1000  0.020   0.255     -0.003   0.370      1.296    0.523    
#> 36 V36  1000  1000  0.001   0.968     -0.007   0.667      0.367    0.832    
#> 37 V37  1000  1000  0.022   0.204      0.003   0.441      5.755    0.056   .
#> 38 V38  1000  1000  0.012   0.495      0.000   0.534      0.473    0.789    
#> 39 V39  1000  1000  0.000   0.998     -0.009   0.874      0.420    0.811    
#> 40 V40  1000  1000 -0.004   0.535      0.009   0.610      0.549    0.760    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 5
#>   - RDIF statistic used for purification: RDIF(RS)
#>   - DIF Items identified by RDIF(RS): 
#>     1, 3, 4, 5, 6 
#>   - RDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc  rdifr p.rdifr      rdifs p.rdifs     rdifrs p.rdifrs
#> 1   V1      2  1000  1000 -0.106   0.000 ***  0.059   0.000 *** 44.880    0.000
#> 2   V2      5  1000  1000  0.009   0.606      0.001   0.512      1.714    0.424
#> 3   V3      4  1000  1000  0.048   0.007  **  0.021   0.005  **  7.782    0.020
#> 4   V4      0  1000  1000  0.141   0.000 ***  0.064   0.000 *** 56.025    0.000
#> 5   V5      3  1000  1000  0.004   0.827      0.022   0.008  ** 10.367    0.006
#> 6   V6      1  1000  1000 -0.072   0.000 ***  0.058   0.000 *** 47.423    0.000
#> 7   V7      5  1000  1000 -0.001   0.955      0.003   0.904      0.362    0.834
#> 8   V8      5  1000  1000  0.024   0.153      0.014   0.094   .  3.336    0.189
#> 9   V9      5  1000  1000 -0.010   0.236      0.016   0.309      1.411    0.494
#> 10 V10      5  1000  1000 -0.016   0.262      0.033   0.257      1.451    0.484
#> 11 V11      5  1000  1000  0.007   0.587      0.022   0.896      0.378    0.828
#> 12 V12      5  1000  1000 -0.002   0.801      0.010   0.800      0.071    0.965
#> 13 V13      5  1000  1000  0.019   0.096   .  0.008   0.190      2.822    0.244
#> 14 V14      5  1000  1000 -0.003   0.790      0.027   0.351      1.701    0.427
#> 15 V15      5  1000  1000 -0.019   0.293     -0.011   0.317      1.367    0.505
#> 16 V16      5  1000  1000  0.012   0.215      0.013   0.540      1.698    0.428
#> 17 V17      5  1000  1000 -0.001   0.889      0.009   0.730      0.643    0.725
#> 18 V18      5  1000  1000  0.010   0.349      0.007   0.108      2.871    0.238
#> 19 V19      5  1000  1000  0.001   0.941     -0.007   0.755      0.160    0.923
#> 20 V20      5  1000  1000  0.001   0.894      0.009   0.860      0.033    0.984
#> 21 V21      5  1000  1000 -0.029   0.099   .  0.018   0.871      2.718    0.257
#> 22 V22      5  1000  1000  0.016   0.359      0.004   0.336      1.024    0.599
#> 23 V23      5  1000  1000 -0.004   0.791     -0.003   0.809      0.127    0.939
#> 24 V24      5  1000  1000  0.024   0.182      0.010   0.134      2.319    0.314
#> 25 V25      5  1000  1000  0.029   0.093   .  0.021   0.859      3.115    0.211
#> 26 V26      5  1000  1000 -0.002   0.914      0.020   0.628      0.641    0.726
#> 27 V27      5  1000  1000  0.018   0.283      0.009   0.255      2.196    0.334
#> 28 V28      5  1000  1000  0.007   0.677      0.032   0.133      2.580    0.275
#> 29 V29      5  1000  1000 -0.015   0.085   .  0.018   0.563      4.394    0.111
#> 30 V30      5  1000  1000 -0.006   0.742     -0.005   0.545      2.853    0.240
#> 31 V31      5  1000  1000 -0.031   0.062   .  0.035   0.151      3.959    0.138
#> 32 V32      5  1000  1000  0.000   0.995      0.026   0.494      0.659    0.719
#> 33 V33      5  1000  1000  0.021   0.274      0.009   0.278      1.721    0.423
#> 34 V34      5  1000  1000  0.005   0.750      0.004   0.750      0.102    0.950
#> 35 V35      5  1000  1000  0.014   0.432     -0.001   0.311      1.025    0.599
#> 36 V36      5  1000  1000 -0.008   0.637     -0.007   0.680      0.243    0.886
#> 37 V37      5  1000  1000  0.020   0.245      0.003   0.454      4.469    0.107
#> 38 V38      5  1000  1000  0.009   0.596      0.000   0.568      0.327    0.849
#> 39 V39      5  1000  1000 -0.002   0.900     -0.009   0.835      0.146    0.929
#> 40 V40      5  1000  1000 -0.004   0.481      0.008   0.556      0.672    0.715
#>       
#> 1  ***
#> 2     
#> 3    *
#> 4  ***
#> 5   **
#> 6  ***
#> 7     
#> 8     
#> 9     
#> 10    
#> 11    
#> 12    
#> 13    
#> 14    
#> 15    
#> 16    
#> 17    
#> 18    
#> 19    
#> 20    
#> 21    
#> 22    
#> 23    
#> 24    
#> 25    
#> 26    
#> 27    
#> 28    
#> 29    
#> 30    
#> 31    
#> 32    
#> 33    
#> 34    
#> 35    
#> 36    
#> 37    
#> 38    
#> 39    
#> 40    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05

# DIF items identified after purification
rdif_pur$with_purify$dif_item
#> [1] 1 3 4 5 6
```

The `with_purify` component contains the final results after
purification converged. The `n.iter` field reports how many iterations
were performed, and `complete` indicates whether convergence was
achieved before reaching `max.iter`.

------------------------------------------------------------------------

## Part 2: GRDIF — Multiple-Group DIF Detection (`grdif()`)

### Statistical Framework

The GRDIF framework (Lim et al. 2024) generalizes the two-group RDIF
approach to simultaneously detect DIF across $`G \geq 2`$ groups. Three
statistics— $`\text{GRDIF}_R`$, $`\text{GRDIF}_S`$, and
$`\text{GRDIF}_{RS}`$—are constructed by applying a contrast matrix
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
focal group. Items 1–2 exhibit uniform DIF (b-parameters shifted upward)
in both focal groups, with a larger shift for Focal Group 2.

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
  verbose    = FALSE
)

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
#>     1, 13 
#>   - DIF Items identified by GRDIF(S): 
#>     1 
#>   - DIF Items identified by GRDIF(RS): 
#>     1, 23 
#>   - GRDIF Statistics: 
#> 
#>     id n.ref n.foc1 n.foc2 grdifr p.grdifr     grdifs p.grdifs     grdifrs
#> 1   V1   700    700    700 62.034    0.000 *** 22.045    0.000 ***  62.925
#> 2   V2   700    700    700  0.145    0.930      0.121    0.941       0.231
#> 3   V3   700    700    700  1.804    0.406      2.041    0.360       3.216
#> 4   V4   700    700    700  0.334    0.846      0.581    0.748       1.679
#> 5   V5   700    700    700  0.870    0.647      0.172    0.918       1.062
#> 6   V6   700    700    700  1.606    0.448      2.121    0.346       3.260
#> 7   V7   700    700    700  0.896    0.639      0.963    0.618       1.330
#> 8   V8   700    700    700  1.273    0.529      1.148    0.563       4.741
#> 9   V9   700    700    700  0.279    0.870      0.442    0.802       0.763
#> 10 V10   700    700    700  0.155    0.925      0.248    0.883       1.355
#> 11 V11   700    700    700  0.855    0.652      1.540    0.463       2.033
#> 12 V12   700    700    700  0.554    0.758      1.302    0.521       1.305
#> 13 V13   700    700    700  6.969    0.031   *  3.796    0.150       7.363
#> 14 V14   700    700    700  1.776    0.412      0.555    0.758       1.986
#> 15 V15   700    700    700  2.773    0.250      4.983    0.083   .   9.300
#> 16 V16   700    700    700  0.321    0.852      0.698    0.705       2.899
#> 17 V17   700    700    700  1.372    0.504      1.317    0.518       3.323
#> 18 V18   700    700    700  5.349    0.069   .  2.811    0.245       6.106
#> 19 V19   700    700    700  3.051    0.218      1.717    0.424       3.900
#> 20 V20   700    700    700  0.410    0.815      0.206    0.902       0.456
#> 21 V21   700    700    700  3.410    0.182      2.712    0.258       5.286
#> 22 V22   700    700    700  3.489    0.175      1.971    0.373       5.564
#> 23 V23   700    700    700  3.126    0.210      3.126    0.210      26.016
#> 24 V24   700    700    700  1.639    0.441      2.002    0.368       2.198
#> 25 V25   700    700    700  2.298    0.317      0.870    0.647       2.677
#> 26 V26   700    700    700  1.769    0.413      0.235    0.889       2.721
#> 27 V27   700    700    700  2.731    0.255      3.002    0.223       5.336
#> 28 V28   700    700    700  0.286    0.867      0.941    0.625       1.344
#> 29 V29   700    700    700  1.176    0.555      0.542    0.763       1.202
#> 30 V30   700    700    700  2.674    0.263      2.811    0.245       6.652
#> 31 V31   700    700    700  1.200    0.549      0.005    0.997       1.401
#> 32 V32   700    700    700  1.118    0.572      2.196    0.334       2.363
#> 33 V33   700    700    700  3.952    0.139      0.043    0.979       4.593
#> 34 V34   700    700    700  0.453    0.797      0.453    0.797       0.307
#> 35 V35   700    700    700  1.978    0.372      0.600    0.741       6.069
#> 36 V36   700    700    700  2.111    0.348      1.043    0.594       2.119
#> 37 V37   700    700    700  1.405    0.495      0.809    0.667       4.549
#> 38 V38   700    700    700  1.829    0.401      2.701    0.259       3.280
#> 39 V39   700    700    700  1.135    0.567      1.116    0.572       2.546
#> 40 V40   700    700    700  1.368    0.505      1.070    0.586       2.165
#>    p.grdifrs    
#> 1      0.000 ***
#> 2      0.994    
#> 3      0.522    
#> 4      0.794    
#> 5      0.900    
#> 6      0.515    
#> 7      0.856    
#> 8      0.315    
#> 9      0.943    
#> 10     0.852    
#> 11     0.730    
#> 12     0.861    
#> 13     0.118    
#> 14     0.738    
#> 15     0.054   .
#> 16     0.575    
#> 17     0.505    
#> 18     0.191    
#> 19     0.420    
#> 20     0.978    
#> 21     0.259    
#> 22     0.234    
#> 23     0.000 ***
#> 24     0.699    
#> 25     0.613    
#> 26     0.606    
#> 27     0.254    
#> 28     0.854    
#> 29     0.878    
#> 30     0.156    
#> 31     0.844    
#> 32     0.669    
#> 33     0.332    
#> 34     0.989    
#> 35     0.194    
#> 36     0.714    
#> 37     0.337    
#> 38     0.512    
#> 39     0.636    
#> 40     0.706    
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
#> [1]  1 13
#> 
#> $grdifs
#> [1] 1
#> 
#> $grdifrs
#> [1]  1 23
```

``` r

# Post-hoc pairwise RDIF results for DIF-flagged items
# (identifies which specific group pairs drive the DIF signal)
grdif_npur$no_purify$post.hoc
#> $by.grdifr
#>    id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1  V1    G1 & G2 -0.0898 -4.5465  0.0426  3.3077 21.4911  0.0000  0.0009
#> 2  V1    G1 & G3 -0.1578 -7.8153  0.0658  4.5342 61.7211  0.0000  0.0000
#> 3  V1    G2 & G3 -0.0680 -3.3093  0.0232  1.1267 11.2517  0.0009  0.2598
#> 4 V13    G1 & G2  0.0283  2.1288 -0.0088 -1.6613  4.6558  0.0333  0.0966
#> 5 V13    G1 & G3  0.0331  2.3526  0.0053 -1.6850  5.8566  0.0186  0.0920
#> 6 V13    G2 & G3  0.0048  0.3327  0.0141 -0.0651  0.2840  0.7394  0.9481
#>   p.rdifrs n.ref n.foc n.total
#> 1   0.0000   700   700    1400
#> 2   0.0000   700   700    1400
#> 3   0.0036   700   700    1400
#> 4   0.0975   700   700    1400
#> 5   0.0535   700   700    1400
#> 6   0.8676   700   700    1400
#> 
#> $by.grdifs
#>   id group.pair   rdifr z.rdifr  rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs
#> 1 V1    G1 & G2 -0.0898 -4.5465 0.0426  3.3077 21.4911   0e+00  0.0009   0.0000
#> 2 V1    G1 & G3 -0.1578 -7.8153 0.0658  4.5342 61.7211   0e+00  0.0000   0.0000
#> 3 V1    G2 & G3 -0.0680 -3.3093 0.0232  1.1267 11.2517   9e-04  0.2598   0.0036
#>   n.ref n.foc n.total
#> 1   700   700    1400
#> 2   700   700    1400
#> 3   700   700    1400
#> 
#> $by.grdifrs
#>    id group.pair   rdifr z.rdifr   rdifs z.rdifs  rdifrs p.rdifr p.rdifs
#> 1  V1    G1 & G2 -0.0898 -4.5465  0.0426  3.3077 21.4911  0.0000  0.0009
#> 2  V1    G1 & G3 -0.1578 -7.8153  0.0658  4.5342 61.7211  0.0000  0.0000
#> 3  V1    G2 & G3 -0.0680 -3.3093  0.0232  1.1267 11.2517  0.0009  0.2598
#> 4 V23    G1 & G2  0.0329  1.7629  0.0235  1.7629  3.1079  0.0779  0.0779
#> 5 V23    G1 & G3  0.0143  0.7665  0.0102  0.7665  0.5840  0.4434  0.4434
#> 6 V23    G2 & G3 -0.0186 -0.9964 -0.0133 -0.9964  0.9908  0.3190  0.3190
#>   p.rdifrs n.ref n.foc n.total
#> 1   0.0000   700   700    1400
#> 2   0.0000   700   700    1400
#> 3   0.0036   700   700    1400
#> 4   0.2114   700   700    1400
#> 5   0.7468   700   700    1400
#> 6   0.6093   700   700    1400
```

The `post.hoc` component reports pairwise RDIF$`_R`$, RDIF$`_S`$, and
RDIF$`_{RS}`$ statistics and p-values for each pair of groups among the
items that were flagged by the omnibus test. This helps pinpoint whether
DIF occurs between the reference and one focal group, between the
reference and all focal groups, or between focal groups only.

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
#>     1, 13 
#>   - DIF Items identified by GRDIF(S): 
#>     1 
#>   - DIF Items identified by GRDIF(RS): 
#>     1, 23 
#>   - GRDIF Statistics: 
#> 
#>     id n.ref n.foc1 n.foc2 grdifr p.grdifr     grdifs p.grdifs     grdifrs
#> 1   V1   700    700    700 62.034    0.000 *** 22.045    0.000 ***  62.925
#> 2   V2   700    700    700  0.145    0.930      0.121    0.941       0.231
#> 3   V3   700    700    700  1.804    0.406      2.041    0.360       3.216
#> 4   V4   700    700    700  0.334    0.846      0.581    0.748       1.679
#> 5   V5   700    700    700  0.870    0.647      0.172    0.918       1.062
#> 6   V6   700    700    700  1.606    0.448      2.121    0.346       3.260
#> 7   V7   700    700    700  0.896    0.639      0.963    0.618       1.330
#> 8   V8   700    700    700  1.273    0.529      1.148    0.563       4.741
#> 9   V9   700    700    700  0.279    0.870      0.442    0.802       0.763
#> 10 V10   700    700    700  0.155    0.925      0.248    0.883       1.355
#> 11 V11   700    700    700  0.855    0.652      1.540    0.463       2.033
#> 12 V12   700    700    700  0.554    0.758      1.302    0.521       1.305
#> 13 V13   700    700    700  6.969    0.031   *  3.796    0.150       7.363
#> 14 V14   700    700    700  1.776    0.412      0.555    0.758       1.986
#> 15 V15   700    700    700  2.773    0.250      4.983    0.083   .   9.300
#> 16 V16   700    700    700  0.321    0.852      0.698    0.705       2.899
#> 17 V17   700    700    700  1.372    0.504      1.317    0.518       3.323
#> 18 V18   700    700    700  5.349    0.069   .  2.811    0.245       6.106
#> 19 V19   700    700    700  3.051    0.218      1.717    0.424       3.900
#> 20 V20   700    700    700  0.410    0.815      0.206    0.902       0.456
#> 21 V21   700    700    700  3.410    0.182      2.712    0.258       5.286
#> 22 V22   700    700    700  3.489    0.175      1.971    0.373       5.564
#> 23 V23   700    700    700  3.126    0.210      3.126    0.210      26.016
#> 24 V24   700    700    700  1.639    0.441      2.002    0.368       2.198
#> 25 V25   700    700    700  2.298    0.317      0.870    0.647       2.677
#> 26 V26   700    700    700  1.769    0.413      0.235    0.889       2.721
#> 27 V27   700    700    700  2.731    0.255      3.002    0.223       5.336
#> 28 V28   700    700    700  0.286    0.867      0.941    0.625       1.344
#> 29 V29   700    700    700  1.176    0.555      0.542    0.763       1.202
#> 30 V30   700    700    700  2.674    0.263      2.811    0.245       6.652
#> 31 V31   700    700    700  1.200    0.549      0.005    0.997       1.401
#> 32 V32   700    700    700  1.118    0.572      2.196    0.334       2.363
#> 33 V33   700    700    700  3.952    0.139      0.043    0.979       4.593
#> 34 V34   700    700    700  0.453    0.797      0.453    0.797       0.307
#> 35 V35   700    700    700  1.978    0.372      0.600    0.741       6.069
#> 36 V36   700    700    700  2.111    0.348      1.043    0.594       2.119
#> 37 V37   700    700    700  1.405    0.495      0.809    0.667       4.549
#> 38 V38   700    700    700  1.829    0.401      2.701    0.259       3.280
#> 39 V39   700    700    700  1.135    0.567      1.116    0.572       2.546
#> 40 V40   700    700    700  1.368    0.505      1.070    0.586       2.165
#>    p.grdifrs    
#> 1      0.000 ***
#> 2      0.994    
#> 3      0.522    
#> 4      0.794    
#> 5      0.900    
#> 6      0.515    
#> 7      0.856    
#> 8      0.315    
#> 9      0.943    
#> 10     0.852    
#> 11     0.730    
#> 12     0.861    
#> 13     0.118    
#> 14     0.738    
#> 15     0.054   .
#> 16     0.575    
#> 17     0.505    
#> 18     0.191    
#> 19     0.420    
#> 20     0.978    
#> 21     0.259    
#> 22     0.234    
#> 23     0.000 ***
#> 24     0.699    
#> 25     0.613    
#> 26     0.606    
#> 27     0.254    
#> 28     0.854    
#> 29     0.878    
#> 30     0.156    
#> 31     0.844    
#> 32     0.669    
#> 33     0.332    
#> 34     0.989    
#> 35     0.194    
#> 36     0.714    
#> 37     0.337    
#> 38     0.512    
#> 39     0.636    
#> 40     0.706    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 1
#>   - GRDIF statistic used for purification: GRDIF(R)
#>   - DIF Items identified by GRDIF(R): 
#>     1 
#>   - GRDIF Statistics: 
#> 
#>     id n.iter n.ref n.foc1 n.foc2 grdifr p.grdifr     grdifs p.grdifs    
#> 1   V1      0   700    700    700 62.034    0.000 *** 22.045    0.000 ***
#> 2   V2      1   700    700    700  0.154    0.926      0.129    0.938    
#> 3   V3      1   700    700    700  1.855    0.396      2.009    0.366    
#> 4   V4      1   700    700    700  0.390    0.823      0.638    0.727    
#> 5   V5      1   700    700    700  0.759    0.684      0.168    0.919    
#> 6   V6      1   700    700    700  1.143    0.565      1.478    0.478    
#> 7   V7      1   700    700    700  0.880    0.644      0.939    0.625    
#> 8   V8      1   700    700    700  1.275    0.529      1.151    0.562    
#> 9   V9      1   700    700    700  0.092    0.955      0.315    0.855    
#> 10 V10      1   700    700    700  0.447    0.800      0.130    0.937    
#> 11 V11      1   700    700    700  1.405    0.495      1.682    0.431    
#> 12 V12      1   700    700    700  0.260    0.878      0.734    0.693    
#> 13 V13      1   700    700    700  5.671    0.059   .  2.873    0.238    
#> 14 V14      1   700    700    700  2.157    0.340      0.594    0.743    
#> 15 V15      1   700    700    700  2.196    0.334      4.644    0.098   .
#> 16 V16      1   700    700    700  0.156    0.925      0.470    0.791    
#> 17 V17      1   700    700    700  1.146    0.564      1.001    0.606    
#> 18 V18      1   700    700    700  4.822    0.090   .  2.231    0.328    
#> 19 V19      1   700    700    700  3.142    0.208      1.315    0.518    
#> 20 V20      1   700    700    700  0.471    0.790      0.421    0.810    
#> 21 V21      1   700    700    700  4.618    0.099   .  2.965    0.227    
#> 22 V22      1   700    700    700  3.484    0.175      2.199    0.333    
#> 23 V23      1   700    700    700  3.126    0.210      3.126    0.210    
#> 24 V24      1   700    700    700  1.871    0.392      2.192    0.334    
#> 25 V25      1   700    700    700  3.064    0.216      0.611    0.737    
#> 26 V26      1   700    700    700  1.248    0.536      0.085    0.958    
#> 27 V27      1   700    700    700  2.708    0.258      2.533    0.282    
#> 28 V28      1   700    700    700  1.111    0.574      0.770    0.680    
#> 29 V29      1   700    700    700  1.695    0.428      0.736    0.692    
#> 30 V30      1   700    700    700  2.388    0.303      2.360    0.307    
#> 31 V31      1   700    700    700  1.303    0.521      0.059    0.971    
#> 32 V32      1   700    700    700  0.991    0.609      1.286    0.526    
#> 33 V33      1   700    700    700  3.336    0.189      0.004    0.998    
#> 34 V34      1   700    700    700  0.453    0.797      0.453    0.797    
#> 35 V35      1   700    700    700  1.874    0.392      0.573    0.751    
#> 36 V36      1   700    700    700  1.815    0.404      0.947    0.623    
#> 37 V37      1   700    700    700  1.458    0.482      0.817    0.665    
#> 38 V38      1   700    700    700  1.846    0.397      2.718    0.257    
#> 39 V39      1   700    700    700  1.132    0.568      1.114    0.573    
#> 40 V40      1   700    700    700  1.214    0.545      0.908    0.635    
#>    grdifrs p.grdifrs    
#> 1   62.925     0.000 ***
#> 2    0.250     0.993    
#> 3    3.238     0.519    
#> 4    1.653     0.799    
#> 5    0.909     0.923    
#> 6    2.325     0.676    
#> 7    1.228     0.873    
#> 8    4.701     0.319    
#> 9    0.672     0.955    
#> 10   1.732     0.785    
#> 11   2.741     0.602    
#> 12   0.741     0.946    
#> 13   6.178     0.186    
#> 14   2.332     0.675    
#> 15   8.560     0.073   .
#> 16   2.543     0.637    
#> 17   3.078     0.545    
#> 18   5.732     0.220    
#> 19   3.806     0.433    
#> 20   0.597     0.963    
#> 21   6.634     0.156    
#> 22   5.463     0.243    
#> 23  91.654     0.000 ***
#> 24   2.470     0.650    
#> 25   3.264     0.515    
#> 26   2.090     0.719    
#> 27   4.977     0.290    
#> 28   2.057     0.725    
#> 29   1.787     0.775    
#> 30   5.994     0.200    
#> 31   1.543     0.819    
#> 32   1.658     0.798    
#> 33   3.861     0.425    
#> 34   0.476     0.976    
#> 35   5.881     0.208    
#> 36   1.819     0.769    
#> 37   4.744     0.315    
#> 38   3.270     0.514    
#> 39   2.531     0.639    
#> 40   2.008     0.734    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05

# DIF items identified after purification
grdif_pur$with_purify$dif_item
#> [1] 1

# Post-hoc pairwise RDIF results after purification
grdif_pur$with_purify$post.hoc
#>   id group.pair   rdifr z.rdifr  rdifs z.rdifs  rdifrs p.rdifr p.rdifs p.rdifrs
#> 1 V1    G1 & G2 -0.0898 -4.5465 0.0426  3.3077 21.4911   0e+00  0.0009   0.0000
#> 2 V1    G1 & G3 -0.1578 -7.8153 0.0658  4.5342 61.7211   0e+00  0.0000   0.0000
#> 3 V1    G2 & G3 -0.0680 -3.3093 0.0232  1.1267 11.2517   9e-04  0.2598   0.0036
#>   n.ref n.foc n.total n_iter
#> 1   700   700    1400      0
#> 2   700   700    1400      0
#> 3   700   700    1400      0
```

------------------------------------------------------------------------

## Part 3: CATSIB (`catsib()`)

### Statistical Framework

CATSIB (Nandakumar and Roussos 2004) is a modified version of SIBTEST
(Shealy and Stout 1993) adapted for computerized adaptive testing (CAT)
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
and Stout 1993):

``` math
\hat{\theta}^*_G = \bar{\theta}_G + \hat{\rho}^2_G\left(\hat{\theta} - \bar{\hat{\theta}}_G\right)
```

where $`\bar{\theta}_G`$ is the group mean of $`\hat{\theta}`$,
$`\bar{\hat{\theta}}_G`$ is also estimated by $`\bar{\theta}_G`$, and
$`\hat{\rho}^2_G = 1 - \hat{\sigma}^2_{e,G} / \hat{\sigma}^2_{\hat{\theta},G}`$
is the estimated reliability in group $`G`$, with
$`\hat{\sigma}^2_{e,G}`$ being the mean squared standard error (SE) of
ability estimates and $`\hat{\sigma}^2_{\hat{\theta},G}`$ being the
observed variance of $`\hat{\theta}`$ in group $`G`$. Examinees are then
matched on $`\hat{\theta}^*`$ rather than $`\hat{\theta}`$.

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
following the recommendation of Nandakumar and Roussos (2004).
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
> purification is applied — in which case `x` must be provided for
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
#>     1, 3, 4, 6, 13, 18, 24, 25, 27, 33 
#>   - Test Statistic: 
#> 
#>         id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   item.1   950   973    1923  0.110 0.016  6.708 0.000 ***
#> 2   item.2   950   973    1923 -0.017 0.017 -0.973 0.330    
#> 3   item.3   950   973    1923 -0.063 0.019 -3.368 0.001 ***
#> 4   item.4   950   973    1923 -0.151 0.019 -7.779 0.000 ***
#> 5   item.5   950   973    1923 -0.017 0.020 -0.856 0.392    
#> 6   item.6   950   973    1923  0.070 0.021  3.353 0.001 ***
#> 7   item.7   950   973    1923 -0.006 0.016 -0.388 0.698    
#> 8   item.8   950   973    1923 -0.028 0.017 -1.643 0.100   .
#> 9   item.9   950   973    1923  0.009 0.007  1.157 0.247    
#> 10 item.10   950   973    1923  0.011 0.016  0.720 0.471    
#> 11 item.11   950   973    1923 -0.017 0.014 -1.254 0.210    
#> 12 item.12   950   973    1923  0.002 0.006  0.428 0.668    
#> 13 item.13   950   973    1923 -0.023 0.012 -1.980 0.048   *
#> 14 item.14   950   973    1923 -0.001 0.012 -0.112 0.911    
#> 15 item.15   950   973    1923  0.004 0.019  0.223 0.823    
#> 16 item.16   950   973    1923 -0.017 0.010 -1.806 0.071   .
#> 17 item.17   950   973    1923  0.002 0.007  0.286 0.775    
#> 18 item.18   950   973    1923 -0.020 0.010 -1.995 0.046   *
#> 19 item.19   950   973    1923 -0.020 0.020 -0.999 0.318    
#> 20 item.20   950   973    1923 -0.003 0.007 -0.396 0.692    
#> 21 item.21   950   973    1923  0.004 0.018  0.224 0.823    
#> 22 item.22   950   973    1923 -0.032 0.018 -1.765 0.078   .
#> 23 item.23   950   973    1923  0.005 0.017  0.328 0.743    
#> 24 item.24   950   973    1923 -0.041 0.018 -2.229 0.026   *
#> 25 item.25   950   973    1923 -0.055 0.018 -3.056 0.002  **
#> 26 item.26   950   973    1923 -0.011 0.017 -0.666 0.505    
#> 27 item.27   950   973    1923 -0.049 0.017 -2.884 0.004  **
#> 28 item.28   950   973    1923 -0.032 0.017 -1.906 0.057   .
#> 29 item.29   950   973    1923  0.009 0.008  1.142 0.254    
#> 30 item.30   950   973    1923  0.002 0.020  0.076 0.940    
#> 31 item.31   950   973    1923  0.016 0.017  0.940 0.347    
#> 32 item.32   950   973    1923 -0.006 0.012 -0.522 0.601    
#> 33 item.33   950   973    1923 -0.043 0.020 -2.192 0.028   *
#> 34 item.34   950   973    1923 -0.008 0.017 -0.482 0.630    
#> 35 item.35   950   973    1923 -0.027 0.018 -1.472 0.141    
#> 36 item.36   950   973    1923 -0.013 0.018 -0.751 0.453    
#> 37 item.37   950   973    1923 -0.017 0.018 -0.939 0.348    
#> 38 item.38   950   973    1923 -0.020 0.018 -1.090 0.276    
#> 39 item.39   950   973    1923  0.005 0.017  0.288 0.773    
#> 40 item.40   950   973    1923  0.003 0.006  0.479 0.632    
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
#>  [1]  1  3  4  6 13 18 24 25 27 33
```

The `dif_stat` data frame contains the estimated $`\hat{\beta}`$
statistic, its standard error (`se`), the standardized $`\hat{\beta}`$
(`z.beta`), and the corresponding p-value for each item. Items with
$`|\hat{\beta}|`$ significantly different from zero (at the `alpha`
level) are flagged as DIF.

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
#>     1, 3, 4, 6, 13, 18, 24, 25, 27, 33 
#>   - Test Statistic: 
#> 
#>     id n.ref n.foc n.total   beta    se z.beta     p    
#> 1   V1   950   973    1923  0.110 0.016  6.708 0.000 ***
#> 2   V2   950   973    1923 -0.017 0.017 -0.973 0.330    
#> 3   V3   950   973    1923 -0.063 0.019 -3.368 0.001 ***
#> 4   V4   950   973    1923 -0.151 0.019 -7.779 0.000 ***
#> 5   V5   950   973    1923 -0.017 0.020 -0.856 0.392    
#> 6   V6   950   973    1923  0.070 0.021  3.353 0.001 ***
#> 7   V7   950   973    1923 -0.006 0.016 -0.388 0.698    
#> 8   V8   950   973    1923 -0.028 0.017 -1.643 0.100   .
#> 9   V9   950   973    1923  0.009 0.007  1.157 0.247    
#> 10 V10   950   973    1923  0.011 0.016  0.720 0.471    
#> 11 V11   950   973    1923 -0.017 0.014 -1.254 0.210    
#> 12 V12   950   973    1923  0.002 0.006  0.428 0.668    
#> 13 V13   950   973    1923 -0.023 0.012 -1.980 0.048   *
#> 14 V14   950   973    1923 -0.001 0.012 -0.112 0.911    
#> 15 V15   950   973    1923  0.004 0.019  0.223 0.823    
#> 16 V16   950   973    1923 -0.017 0.010 -1.806 0.071   .
#> 17 V17   950   973    1923  0.002 0.007  0.286 0.775    
#> 18 V18   950   973    1923 -0.020 0.010 -1.995 0.046   *
#> 19 V19   950   973    1923 -0.020 0.020 -0.999 0.318    
#> 20 V20   950   973    1923 -0.003 0.007 -0.396 0.692    
#> 21 V21   950   973    1923  0.004 0.018  0.224 0.823    
#> 22 V22   950   973    1923 -0.032 0.018 -1.765 0.078   .
#> 23 V23   950   973    1923  0.005 0.017  0.328 0.743    
#> 24 V24   950   973    1923 -0.041 0.018 -2.229 0.026   *
#> 25 V25   950   973    1923 -0.055 0.018 -3.056 0.002  **
#> 26 V26   950   973    1923 -0.011 0.017 -0.666 0.505    
#> 27 V27   950   973    1923 -0.049 0.017 -2.884 0.004  **
#> 28 V28   950   973    1923 -0.032 0.017 -1.906 0.057   .
#> 29 V29   950   973    1923  0.009 0.008  1.142 0.254    
#> 30 V30   950   973    1923  0.002 0.020  0.076 0.940    
#> 31 V31   950   973    1923  0.016 0.017  0.940 0.347    
#> 32 V32   950   973    1923 -0.006 0.012 -0.522 0.601    
#> 33 V33   950   973    1923 -0.043 0.020 -2.192 0.028   *
#> 34 V34   950   973    1923 -0.008 0.017 -0.482 0.630    
#> 35 V35   950   973    1923 -0.027 0.018 -1.472 0.141    
#> 36 V36   950   973    1923 -0.013 0.018 -0.751 0.453    
#> 37 V37   950   973    1923 -0.017 0.018 -0.939 0.348    
#> 38 V38   950   973    1923 -0.020 0.018 -1.090 0.276    
#> 39 V39   950   973    1923  0.005 0.017  0.288 0.773    
#> 40 V40   950   973    1923  0.003 0.006  0.479 0.632    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05 
#> 
#> 
#>  2. With purification 
#> 
#>   - Completion of purification: TRUE
#>   - Number of iterations: 9
#>   - Potential DIF Items: 
#>     1, 3, 4, 6, 24, 25, 27, 28, 33 
#>   - Test Statistic: 
#> 
#>     id n.iter n.ref n.foc n.total   beta    se z.beta     p    
#> 1   V1      1   949   974    1923  0.105 0.016  6.376 0.000 ***
#> 2   V2      9   941   968    1909 -0.016 0.018 -0.916 0.360    
#> 3   V3      3   945   970    1915 -0.059 0.019 -3.181 0.002  **
#> 4   V4      0   950   973    1923 -0.151 0.019 -7.779 0.000 ***
#> 5   V5      9   941   968    1909 -0.015 0.020 -0.748 0.455    
#> 6   V6      2   948   974    1922  0.069 0.021  3.331 0.001 ***
#> 7   V7      9   941   968    1909 -0.006 0.016 -0.401 0.688    
#> 8   V8      9   941   968    1909 -0.027 0.017 -1.565 0.118    
#> 9   V9      9   941   968    1909  0.005 0.007  0.705 0.481    
#> 10 V10      9   941   968    1909  0.006 0.015  0.386 0.700    
#> 11 V11      9   941   968    1909 -0.021 0.013 -1.568 0.117    
#> 12 V12      9   941   968    1909  0.002 0.006  0.274 0.784    
#> 13 V13      9   941   968    1909 -0.022 0.012 -1.946 0.052   .
#> 14 V14      9   941   968    1909 -0.003 0.012 -0.271 0.787    
#> 15 V15      9   941   968    1909 -0.013 0.018 -0.699 0.485    
#> 16 V16      9   941   968    1909 -0.016 0.009 -1.761 0.078   .
#> 17 V17      9   941   968    1909  0.001 0.007  0.137 0.891    
#> 18 V18      9   941   968    1909 -0.015 0.010 -1.470 0.142    
#> 19 V19      9   941   968    1909 -0.018 0.020 -0.879 0.380    
#> 20 V20      9   941   968    1909 -0.009 0.008 -1.123 0.261    
#> 21 V21      9   941   968    1909 -0.004 0.018 -0.206 0.837    
#> 22 V22      9   941   968    1909 -0.027 0.018 -1.489 0.136    
#> 23 V23      9   941   968    1909  0.010 0.017  0.562 0.574    
#> 24 V24      6   930   964    1894 -0.041 0.019 -2.172 0.030   *
#> 25 V25      4   940   968    1908 -0.048 0.018 -2.678 0.007  **
#> 26 V26      9   941   968    1909 -0.011 0.017 -0.671 0.502    
#> 27 V27      5   938   966    1904 -0.037 0.017 -2.185 0.029   *
#> 28 V28      8   938   966    1904 -0.037 0.016 -2.251 0.024   *
#> 29 V29      9   941   968    1909  0.010 0.008  1.220 0.222    
#> 30 V30      9   941   968    1909 -0.004 0.020 -0.180 0.857    
#> 31 V31      9   941   968    1909  0.016 0.018  0.928 0.353    
#> 32 V32      9   941   968    1909 -0.016 0.012 -1.344 0.179    
#> 33 V33      7   938   966    1904 -0.047 0.020 -2.420 0.016   *
#> 34 V34      9   941   968    1909  0.001 0.017  0.061 0.951    
#> 35 V35      9   941   968    1909 -0.032 0.018 -1.749 0.080   .
#> 36 V36      9   941   968    1909 -0.020 0.018 -1.145 0.252    
#> 37 V37      9   941   968    1909 -0.025 0.018 -1.367 0.172    
#> 38 V38      9   941   968    1909 -0.020 0.018 -1.105 0.269    
#> 39 V39      9   941   968    1909 -0.004 0.017 -0.262 0.793    
#> 40 V40      9   941   968    1909  0.003 0.006  0.411 0.681    
#> 
#> '***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1  
#> Significance level: 0.05

# DIF items identified after purification
catsib_pur$with_purify$dif_item
#> [1]  1  3  4  6 24 25 27 28 33
```

### Important Considerations for CATSIB

**Sensitivity to DIF type.** CATSIB, like its predecessor SIBTEST
(Shealy and Stout 1993), was originally designed and validated for
detecting **uniform DIF** — the condition in which the direction of
group differences in item performance is consistent across all ability
levels. The $`\hat{\beta}`$ statistic accumulates the signed difference
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

**Numerical instability of the regression correction under nonuniform or
mixed DIF.** The regression correction (Equation 7 of Nandakumar and
Roussos (2004)) computes, for each group $`G`$, the reliability
estimate:

``` math
\hat{\rho}^2_G = 1 - \frac{\hat{\sigma}^2_{e,G}}{\hat{\sigma}^2_{\hat{\theta},G}}
```

where $`\hat{\sigma}^2_{e,G}`$ is the mean squared SE of ability
estimates and $`\hat{\sigma}^2_{\hat{\theta},G}`$ is the observed
variance of $`\hat{\theta}`$ in group $`G`$. When DIF items with
substantially reduced discrimination are present, their large IRT
standard errors may inflate $`\hat{\sigma}^2_{e,G}`$ for the focal
group, causing $`\hat{\rho}^2_G`$ to approach or fall below zero. A
negative $`\hat{\rho}^2_G`$ would reverse the direction of the
regression correction — compressing scores away from the group mean
rather than toward it — which severely distorts the bin-level matching
and leads to inflated Type I error rates.

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
triggered, results should be interpreted with caution — or
`purify = FALSE` should be considered.

------------------------------------------------------------------------

## Summary of Function Inputs

All three functions share a common calling pattern:

``` r
func(
  x          = <item metadata>,     # pooled item parameters (data frame)
  data       = <response matrix>,   # pooled responses (examinees x items)
  score      = <ability vector>,    # pooled ability estimates
  group      = <group vector>,      # group membership labels
  focal.name = <focal label(s)>,    # which group(s) are focal
  D          = 1.702,               # scaling constant (match calibration)
  alpha      = 0.05,                # significance level
  purify     = TRUE / FALSE,        # apply iterative purification?
  verbose    = FALSE
)
```

Output always contains `$no_purify$dif_stat` (table of statistics for
all items) and `$no_purify$dif_item` (indices of flagged items). When
`purify = TRUE`, corresponding results are also available under
`$with_purify`.

------------------------------------------------------------------------

## References

Lim, Hwanggyu, Edison M. Choe, and Kyung T. Han. 2022. “A Residual-Based
Differential Item Functioning Detection Framework in Item Response
Theory.” *Journal of Educational Measurement* 59 (1): 80–104.
<https://doi.org/10.1111/jedm.12313>.

Lim, Hwanggyu, Danqi Zhu, Edison M. Choe, and Kyung T. Han. 2024.
“Detecting Differential Item Functioning Among Multiple Groups Using IRT
Residual DIF Framework.” *Journal of Educational Measurement* 61 (4):
656–81. <https://doi.org/10.1111/jedm.12415>.

Nandakumar, Ratna, and Louis Roussos. 2004. “Evaluation of the CATSIB
DIF Procedure in a Pretest Setting.” *Journal of Educational and
Behavioral Statistics* 29 (2): 177–99.
<https://doi.org/10.3102/10769986029002177>.

Shealy, Robin T., and William F. Stout. 1993. “A Model-Based
Standardization Approach That Separates True Bias/DIF from Group Ability
Differences and Detects Test Bias/DIF as Well as Item Bias/DIF.”
*Psychometrika* 58 (2): 159–94. <https://doi.org/10.1007/BF02294572>.
