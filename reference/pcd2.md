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

# \donttest{
## ── Example 3: CAT-based IPD detection using simIPD ──────────────────────────
##
## The Pseudo-count D² statistic has no closed-form null distribution.
## Following Lim & Han (in press), the critical value is estimated empirically
## via bootstrap:
##   (1) Select drift-free (anchor) items to form the null D² distribution.
##   (2) Repeatedly resample from those values and take the 95th percentile.
##   (3) Average across bootstrap iterations to obtain the critical value.
##
## The simIPD dataset provides:
##   item_par  : original 360-item pool (3PLM, irtQ format)
##   foc_resp  : focal group CAT responses (N = 3,000; TL = 30; ~92% NA)
##   key_item  : indices of 90 highly-exposed key items (IPD analysis targets)
##   item.skip : indices of 270 non-key items (excluded from analysis)
##   ipd_item  : 18 truly drifted items (ground truth; a & b each −0.5)
## ─────────────────────────────────────────────────────────────────────────────

data(simIPD)

## ── Step 1. Select anchor items for bootstrap ─────────────────────────────
## From the full item pool, select items that (a) have sufficient observed
## responses (>= boot_size) and (b) are not known IPD items. These items
## serve as the empirical null distribution of D².
## boot_size = 300: chosen to match the minimum response count used in the
##   bootstrap procedure of the simulation study (Lim & Han, in press).
## In practice, exclude items you know or suspect have drifted; here the
## ground truth (simIPD$ipd_item) is used for illustration.
boot_size    <- 300
n_resp       <- colSums(!is.na(simIPD$foc_resp))   # observed responses per item
anchor_items <- setdiff(which(n_resp >= boot_size), simIPD$ipd_item)

## ── Step 2. Compute null D² values for the anchor items ──────────────────
pcd2_null <- pcd2(
  x        = simIPD$item_par[anchor_items, ],
  data     = simIPD$foc_resp[, anchor_items],
  D        = 1.7,
  crit.val = NULL,
  purify   = FALSE
)$no_purify$ipd_stat$pcd2

## ── Step 3. Bootstrap critical value (Lim & Han, in press) ──────────────
## For each bootstrap iteration: resample boot_size D² values from the null
## distribution and take the 95th percentile. The critical value is the mean
## of these percentiles across all iterations.
## n_boots = 500: reduced for demo speed; use 10,000 in practice for a
##   stable critical value estimate.
set.seed(2024)
n_boots <- 500

crit_val <- purrr::map_dbl(
  .x = 1:n_boots,
  .f = ~ sample(x = pcd2_null, size = boot_size, replace = TRUE) |>
           stats::quantile(probs = 0.95)
) |> mean()

cat("Bootstrap critical value (alpha = 0.05):", round(crit_val, 6), "\n")
#> Bootstrap critical value (alpha = 0.05): 0.003013 

## ── Step 4. Run PCD2 with the bootstrap critical value + purification ──
pcd2_result <- pcd2(
  x         = simIPD$item_par,
  data      = simIPD$foc_resp,
  D         = 1.7,
  item.skip = simIPD$item.skip,
  crit.val  = crit_val,
  purify    = TRUE,
  max.iter  = 30
)
#> Warning: Item(s) 17, Item(s) 19, Item(s) 22, Item(s) 24, Item(s) 27, Item(s) 30, Item(s) 32, Item(s) 34, Item(s) 44, Item(s) 47, Item(s) 50, Item(s) 54, Item(s) 56, Item(s) 65, Item(s) 73, Item(s) 76, Item(s) 77, Item(s) 82, Item(s) 86, Item(s) 93, Item(s) 102, Item(s) 105, Item(s) 108, Item(s) 111, Item(s) 112, Item(s) 120, Item(s) 122, Item(s) 131, Item(s) 137, Item(s) 140, Item(s) 144, Item(s) 155, Item(s) 156, Item(s) 164, Item(s) 165, Item(s) 174, Item(s) 176, Item(s) 182, Item(s) 193, Item(s) 196, Item(s) 201, Item(s) 206, Item(s) 219, Item(s) 222, Item(s) 223, Item(s) 225, Item(s) 232, Item(s) 241, Item(s) 255, Item(s) 265, Item(s) 269, Item(s) 271, Item(s) 278, Item(s) 279, Item(s) 284, Item(s) 286, Item(s) 289, Item(s) 291, Item(s) 292, Item(s) 293, Item(s) 295, Item(s) 305, Item(s) 313, Item(s) 315, Item(s) 316, Item(s) 339, Item(s) 342, Item(s) 343, Item(s) 344, Item(s) 345, Item(s) 352, Item(s) 353, Item(s) 359 has(have) no item response data. 
#> Purification started... 
#>  Iteration: 1 Iteration: 2 Iteration: 3 Iteration: 4 Iteration: 5 Iteration: 6 Iteration: 7 Iteration: 8 Iteration: 9 Iteration: 10 Iteration: 11 Iteration: 12 Iteration: 13 Iteration: 14 Iteration: 15 Iteration: 16 Iteration: 17 Iteration: 18 Iteration: 19 Iteration: 20 
#> Purification is finished. 

## ── Step 5. Compare detected items to ground truth ────────────────────
detected <- pcd2_result$with_purify$ipd_item
cat("Truly drifted items (ground truth):", simIPD$ipd_item, "\n")
#> Truly drifted items (ground truth): 26 42 55 58 97 129 138 145 149 181 184 203 242 246 252 310 334 336 
cat("PCD2-detected items:               ", detected, "\n")
#> PCD2-detected items:                26 42 55 58 97 129 138 145 149 181 184 203 236 242 246 252 261 310 334 336 
cat("True positives:", sum(detected %in% simIPD$ipd_item), "of",
    length(simIPD$ipd_item), "\n")
#> True positives: 18 of 18 
# }

```
