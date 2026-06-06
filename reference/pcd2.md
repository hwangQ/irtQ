# Pseudo-count D2 method

This function calculates the Pseudo-count \\D^{2}\\ statistic to
evaluate item parameter drift, as described by Cappaert et al. (2018)
and Stone (2000). The Pseudo-count \\D^{2}\\ statistic is designed to
detect item parameter drift efficiently without requiring item
recalibration, making it especially valuable in computerized adaptive
testing (CAT) environments. This method compares observed and expected
response frequencies across quadrature points, which represent latent
ability levels. The expected frequencies are computed using the
posterior distribution of each examinee's ability (Stone, 2000),
providing a robust and sensitive measure of item parameter drift,
ensuring the stability and accuracy of the test over time.

## Usage

``` r
pcd2(
  x,
  data,
  D = 1,
  item.skip = NULL,
  missing = NA,
  Quadrature = c(49, 6),
  weights = NULL,
  group.mean = 0,
  group.var = 1,
  crit.val = NULL,
  min.resp = NULL,
  purify = FALSE,
  max.iter = 10,
  verbose = TRUE
)
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

- data:

  A matrix of examinees' item responses corresponding to the items
  specified in the `x` argument. Rows represent examinees and columns
  represent items.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- item.skip:

  A numeric vector of item indices to exclude from IPD analysis. If
  `NULL`, all items are included. Useful for omitting specific items
  based on prior insights.

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
  the first column) and their corresponding weights (in the second
  column) for the latent variable prior distribution. If not `NULL`, the
  scale of the latent ability distribution is fixed to match the scale
  of the provided quadrature points and weights. The weights and points
  can be conveniently generated using the function
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).

  If `NULL`, a normal prior density is used instead, based on the
  information provided in the `Quadrature`, `group.mean`, and
  `group.var` arguments. Default is `NULL`.

- group.mean:

  A numeric value specifying the mean of the latent variable prior
  distribution when `weights = NULL`. Default is 0. This value is fixed
  to resolve the indeterminacy of the item parameter scale during
  calibration.

- group.var:

  A positive numeric value specifying the variance of the latent
  variable prior distribution when `weights = NULL`. Default is 1. This
  value is fixed to resolve the indeterminacy of the item parameter
  scale during calibration.

- crit.val:

  A critical value applied in hypothesis testing using the Pseudo-count
  \\D^{2}\\ statistic. Default is `NULL`.

- min.resp:

  A positive integer specifying the minimum required number of responses
  for each evaluated item. Defaults to `NULL`.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- max.iter:

  A positive integer specifying the maximum number of iterations allowed
  for the purification process. Default is `10`.

- verbose:

  Logical. If `TRUE`, progress messages from the purification procedure
  will be displayed; if `FALSE`, the messages will be suppressed.
  Default is `TRUE`.

## Value

This function returns a list containing four main components:

- no_purify:

  A list containing the results of Pseudo-count \\D^{2}\\ analysis
  without applying the purification procedure. It includes:

  ipd_stat

  :   A data frame summarizing the Pseudo-count \\D^{2}\\ statistics for
      all items. The columns include: `id` (item ID), `pcd2` (the
      computed \\D^{2}\\ value), and `N` (the number of valid examinee
      responses per item).

  ipd_item

  :   A numeric vector of item indices that were flagged as exhibiting
      item parameter drift (IPD), based on the specified critical value
      `crit.val`. If no items are flagged or `crit.val = NULL`, this is
      `NULL`.

- purify:

  A logical value indicating whether the iterative purification
  procedure was applied (`TRUE`) or not (`FALSE`).

- with_purify:

  A list containing the results of Pseudo-count \\D^{2}\\ analysis after
  applying the purification procedure. This list is populated only when
  both `purify = TRUE` and `crit.val` is not `NULL`. It includes:

  ipd_stat

  :   A data frame reporting the final Pseudo-count \\D^{2}\\ statistics
      after purification. Columns include: `id` (item ID), `pcd2` (the
      computed \\D^{2}\\ value), `N` (the number of valid responses),
      and `n.iter` (the iteration number in which each item was
      evaluated).

  ipd_item

  :   A numeric vector of item indices flagged as IPD items during
      purification. Items are ordered by the iteration in which they
      were flagged.

  n.iter

  :   An integer indicating the number of purification iterations
      completed.

  complete

  :   A logical value indicating whether the purification procedure
      converged before reaching the maximum number of iterations
      (`max.iter`). If `FALSE`, the iteration limit was reached before
      convergence.

- crit.val:

  A numeric value indicating the critical threshold used to flag items
  for parameter drift. If not specified by the user, this will be
  `NULL`.

## Details

The Pseudo-count \\D^{2}\\ statistic quantifies item parameter drift
(IPD) by computing the weighted squared differences between the observed
and expected response frequencies for each score category across ability
levels. The expected frequencies are determined using the posterior
distribution of each examinee's ability (Stone, 2000).

The Pseudo-count \\D^{2}\\ statistic is calculated as: \$\$ Pseudo-count
D^{2} = \sum\_{k=1}^{Q} \left( \frac{r\_{0k} + r\_{1k}}{N}\right) \left(
\frac{r\_{1k}}{r\_{0k} + r\_{1k}} - E\_{1k} \right)^2 \$\$

where \\r\_{0k}\\ and \\r\_{1k}\\ are the pseudo-counts for the
incorrect and correct responses at each ability level \\k\\, \\E\_{1k}\\
is the expected proportion of correct responses at each ability level
\\k\\, calculated using item parameters from the item bank, and \\N\\ is
the total count of examinees who received each item.

**Critical Value (`crit.val`)**: The `crit.val` argument specifies the
threshold used to flag an item as exhibiting potential parameter drift.
If an item's Pseudo-count \\D^{2}\\ value exceeds this threshold, it is
identified as a drifted item. If `crit.val = NULL`, the function reports
the raw statistic without flagging.

**Minimum Response Count (`min.resp`)**: The `min.resp` argument sets a
lower bound on the number of responses required for an item to be
included in the analysis. Items with fewer responses than `min.resp` are
automatically excluded by replacing all their responses with `NA`. This
avoids unreliable estimates based on small sample sizes.

**Purification Procedure**: Although Cappaert et al. (2018) did not
incorporate purification into their method, `pcd2()` implements an
optional iterative purification process similar to Lim et al. (2022).
When `purify = TRUE` and a `crit.val` is provided:

- The procedure begins by identifying items flagged for drift using the
  initial Pseudo-count \\D^{2}\\ statistics.

- In each subsequent iteration, the item with the highest flagged
  Pseudo-count \\D^{2}\\ value is removed from the item set, and the
  statistics are recalculated using only the remaining items.

- The process continues until no additional items are flagged or the
  number of iterations reaches `max.iter`.

- All flagged items and statistics are saved, and convergence status is
  reported.

This process ensures that drift detection is not distorted by
already-flagged items, improving the robustness of the results.

## References

Cappaert, K. J., Wen, Y., & Chang, Y. F. (2018). Evaluating CAT-adjusted
approaches for suspected item parameter drift detection. *Measurement:
Interdisciplinary Research and Perspectives, 16*(4), 226-238.

Stone, C. A. (2000). Monte Carlo based null distribution for an
alternative goodness-of-fit test statistic in IRT models. *Journal of
educational measurement, 37*(1), 58-75.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1: No critical value specified
## Compute the Pseudo-count D² statistics for dichotomous items
## Import the "-prm.txt" output file generated by flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Extract metadata for the first 30 3PLM items
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[1:30, 1:6]

# Generate abilities for 500 examinees from N(0, 1)
set.seed(25)
score <- rnorm(500, mean = 0, sd = 1)

# Simulate response data using the item metadata and ability values
data <- simdat(x = x, theta = score, D = 1)

# Compute the Pseudo-count D² statistics (no purification applied)
ps_d2 <- pcd2(x = x, data = data)
print(ps_d2)
#> $no_purify
#> $no_purify$ipd_stat
#>       id         pcd2   N
#> 1   CMC1 0.0004106214 500
#> 2   CMC2 0.0002856148 500
#> 3   CMC3 0.0002623176 500
#> 4   CMC4 0.0003447557 500
#> 5   CMC5 0.0017572178 500
#> 6   CMC6 0.0009171795 500
#> 7   CMC7 0.0009759268 500
#> 8   CMC8 0.0014454403 500
#> 9   CMC9 0.0008666188 500
#> 10 CMC10 0.0010232559 500
#> 11 CMC11 0.0012532975 500
#> 12 CMC12 0.0011485885 500
#> 13 CMC13 0.0003172974 500
#> 14 CMC14 0.0008779335 500
#> 15 CMC15 0.0018646941 500
#> 16 CMC16 0.0005040721 500
#> 17 CMC17 0.0016632488 500
#> 18 CMC18 0.0014097661 500
#> 19 CMC19 0.0001394952 500
#> 20 CMC20 0.0001029313 500
#> 21 CMC21 0.0008133451 500
#> 22 CMC22 0.0010517606 500
#> 23 CMC23 0.0009793800 500
#> 24 CMC24 0.0027062555 500
#> 25 CMC25 0.0008871344 500
#> 26 CMC26 0.0027792672 500
#> 27 CMC27 0.0009563499 500
#> 28 CMC28 0.0008374350 500
#> 29 CMC29 0.0008160481 500
#> 30 CMC30 0.0018719723 500
#> 
#> $no_purify$ipd_item
#> NULL
#> 
#> 
#> $purify
#> [1] FALSE
#> 
#> $with_purify
#> $with_purify$ipd_stat
#> NULL
#> 
#> $with_purify$ipd_item
#> NULL
#> 
#> $with_purify$n.iter
#> NULL
#> 
#> $with_purify$complete
#> NULL
#> 
#> 
#> $crit.val
#> NULL
#> 
#> $call
#> pcd2(x = x, data = data)
#> 

## Example 2: Applying a critical value with purification
# Compute the Pseudo-count D² statistics with purification enabled
ps_d2_puri <- pcd2(x = x, data = data, crit.val = 0.002, purify = TRUE)
#> Purification started... 
#>  Iteration: 1 Iteration: 2 
#> Purification is finished. 
print(ps_d2_puri)
#> $no_purify
#> $no_purify$ipd_stat
#>       id         pcd2   N
#> 1   CMC1 0.0004106214 500
#> 2   CMC2 0.0002856148 500
#> 3   CMC3 0.0002623176 500
#> 4   CMC4 0.0003447557 500
#> 5   CMC5 0.0017572178 500
#> 6   CMC6 0.0009171795 500
#> 7   CMC7 0.0009759268 500
#> 8   CMC8 0.0014454403 500
#> 9   CMC9 0.0008666188 500
#> 10 CMC10 0.0010232559 500
#> 11 CMC11 0.0012532975 500
#> 12 CMC12 0.0011485885 500
#> 13 CMC13 0.0003172974 500
#> 14 CMC14 0.0008779335 500
#> 15 CMC15 0.0018646941 500
#> 16 CMC16 0.0005040721 500
#> 17 CMC17 0.0016632488 500
#> 18 CMC18 0.0014097661 500
#> 19 CMC19 0.0001394952 500
#> 20 CMC20 0.0001029313 500
#> 21 CMC21 0.0008133451 500
#> 22 CMC22 0.0010517606 500
#> 23 CMC23 0.0009793800 500
#> 24 CMC24 0.0027062555 500
#> 25 CMC25 0.0008871344 500
#> 26 CMC26 0.0027792672 500
#> 27 CMC27 0.0009563499 500
#> 28 CMC28 0.0008374350 500
#> 29 CMC29 0.0008160481 500
#> 30 CMC30 0.0018719723 500
#> 
#> $no_purify$ipd_item
#> [1] 24 26
#> 
#> 
#> $purify
#> [1] TRUE
#> 
#> $with_purify
#> $with_purify$ipd_stat
#>       id         pcd2   N n.iter
#> 1   CMC1 3.632177e-04 500      2
#> 2   CMC2 2.530213e-04 500      2
#> 3   CMC3 2.504885e-04 500      2
#> 4   CMC4 2.797221e-04 500      2
#> 5   CMC5 1.818836e-03 500      2
#> 6   CMC6 9.373936e-04 500      2
#> 7   CMC7 9.042253e-04 500      2
#> 8   CMC8 1.386178e-03 500      2
#> 9   CMC9 6.783198e-04 500      2
#> 10 CMC10 9.282039e-04 500      2
#> 11 CMC11 9.253406e-04 500      2
#> 12 CMC12 1.303119e-03 500      2
#> 13 CMC13 2.924213e-04 500      2
#> 14 CMC14 1.025110e-03 500      2
#> 15 CMC15 1.725936e-03 500      2
#> 16 CMC16 4.113741e-04 500      2
#> 17 CMC17 1.622536e-03 500      2
#> 18 CMC18 1.517624e-03 500      2
#> 19 CMC19 1.097570e-04 500      2
#> 20 CMC20 5.962552e-05 500      2
#> 21 CMC21 6.547380e-04 500      2
#> 22 CMC22 7.675185e-04 500      2
#> 23 CMC23 9.423791e-04 500      2
#> 24 CMC24 2.844727e-03 500      1
#> 25 CMC25 9.089136e-04 500      2
#> 26 CMC26 2.779267e-03 500      0
#> 27 CMC27 8.972233e-04 500      2
#> 28 CMC28 8.151431e-04 500      2
#> 29 CMC29 8.989447e-04 500      2
#> 30 CMC30 1.836625e-03 500      2
#> 
#> $with_purify$ipd_item
#> [1] 24 26
#> 
#> $with_purify$n.iter
#> [1] 2
#> 
#> $with_purify$complete
#> [1] TRUE
#> 
#> 
#> $crit.val
#> [1] 0.002
#> 
#> $call
#> pcd2(x = x, data = data, crit.val = 0.002, purify = TRUE)
#> 
```
