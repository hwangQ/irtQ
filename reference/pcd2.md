# Pseudo-Count D2 Method for Item Parameter Drift Detection

This function calculates the pseudo-count \\D^{2}\\ statistic to
evaluate item parameter drift, as described by Cappaert et al. (2018)
and Stone (2000). The pseudo-count \\D^{2}\\ statistic is designed to
detect item parameter drift efficiently without requiring item
recalibration, making it especially valuable in computerized adaptive
testing (CAT) environments. This method compares observed and expected
response frequencies across quadrature points, which represent latent
ability levels. The observed pseudo-frequencies are computed using the
posterior distribution of each examinee's ability, whereas the expected
frequencies are the response probabilities implied by the item
parameters in the item bank (Stone, 2000). The current version supports
dichotomous items only.

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
  function. Only dichotomous items are supported; the function stops if
  `x` contains GRM or GPCM items.

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

- item.skip:

  A numeric vector of item positions (row numbers of `x`), or a logical
  vector with one value per item, identifying the items to exclude from
  the IPD analysis. Skipped items still contribute to the posterior
  ability distributions; only their `pcd2` values are set to `NA`, and
  they are never flagged. If `NULL`, all items are analyzed. Default is
  `NULL`.

- missing:

  A value indicating missing responses in the data set. Default is `NA`.

- Quadrature:

  A numeric vector of length two giving the number of quadrature points
  and the symmetric bound (absolute value) of the points. For example,
  `c(49, 6)` specifies 49 evenly spaced points from -6 to 6. These
  points are used to compute each examinee's posterior ability
  distribution when `weights = NULL`. Default is `c(49, 6)`.

- weights:

  A two-column matrix or data frame of quadrature points (first column)
  and weights (second column) defining the prior ability distribution
  used to compute the posterior distributions. It can be created with
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).
  If `NULL`, a normal prior defined by `Quadrature`, `group.mean`, and
  `group.var` is used. Default is `NULL`.

- group.mean:

  A numeric value specifying the mean of the normal prior ability
  distribution when `weights = NULL`. Default is 0.

- group.var:

  A positive numeric value specifying the variance of the normal prior
  ability distribution when `weights = NULL`. Default is 1.

- crit.val:

  A numeric critical value for the pseudo-count \\D^{2}\\ statistic.
  Items whose statistic is greater than or equal to `crit.val` are
  flagged. Because the statistic has no known null distribution, the
  critical value is usually obtained by bootstrapping (see Examples). If
  `NULL` (default), the statistic is reported but no item is flagged and
  purification is not performed.

- min.resp:

  A single nonnegative number, usually a positive integer, specifying
  the minimum number of responses an item must have to be analyzed. All
  responses of items with fewer than `min.resp` (but at least one)
  responses are set to `NA`, so their `pcd2` value is `NaN` and `N` is
  0, and a warning lists those items. If `NULL`, no minimum is applied.
  Default is `NULL`.

- purify:

  Logical. Whether to apply the iterative purification procedure.
  Purification is performed only when `crit.val` is also supplied.
  Default is `FALSE`.

- max.iter:

  A positive integer specifying the maximum number of purification
  iterations. Default is 10. If the limit is reached while flagged items
  remain, a warning is issued, the items flagged at the last iteration
  are added to the flagged set, and the `complete` element of
  `with_purify` is `FALSE`. Because one item is removed per iteration,
  `max.iter` should be at least the number of items expected to drift.

- verbose:

  Logical. If `TRUE`, the progress of the purification procedure is
  printed to the console. Default is `TRUE`.

## Value

This function returns a list with the following elements:

- no_purify:

  A list containing the results of pseudo-count \\D^{2}\\ analysis
  without applying the purification procedure. It includes:

  ipd_stat

  :   A data frame summarizing the pseudo-count \\D^{2}\\ statistics for
      all items. The columns include: `id` (item ID), `pcd2` (the
      computed \\D^{2}\\ value), and `N` (the number of valid examinee
      responses per item).

  ipd_item

  :   A numeric vector of item indices that were flagged as exhibiting
      item parameter drift (IPD), based on the specified critical value
      `crit.val`. If no items are flagged or `crit.val = NULL`, this is
      `NULL`.

- purify:

  The value of the `purify` argument. Purification is actually performed
  only when `crit.val` is also supplied.

- with_purify:

  A list containing the results of pseudo-count \\D^{2}\\ analysis after
  applying the purification procedure. This list is populated only when
  both `purify = TRUE` and `crit.val` is not `NULL`; otherwise, all of
  its elements are `NULL`. It includes:

  ipd_stat

  :   A data frame of the pseudo-count \\D^{2}\\ statistics for all
      items, taken for each item from the analysis given in `n.iter`.
      Columns include: `id` (item ID), `pcd2` (the computed \\D^{2}\\
      value), `N` (the number of valid responses), and `n.iter` (for a
      flagged item, the iteration at which it was flagged, where 0 is
      the initial analysis; for the other items, the final iteration).

  ipd_item

  :   A numeric vector of item indices (row positions in `x`) flagged as
      IPD items during purification, sorted in ascending order. The
      iteration at which each item was flagged is given in the `n.iter`
      column of `ipd_stat`. It contains all items if every item is
      flagged. It is `NULL` when no item is flagged in the initial
      analysis.

  n.iter

  :   An integer indicating the number of purification iterations
      completed. The iteration that removes the last item is not
      counted, because the statistics are not recomputed.

  complete

  :   A logical value indicating whether the purification procedure
      converged before reaching the maximum number of iterations
      (`max.iter`). `FALSE` if the iteration limit was reached before
      convergence or if every item was flagged.

- crit.val:

  A numeric value indicating the critical threshold used to flag items
  for parameter drift. If not specified by the user, this will be
  `NULL`.

- call:

  The matched function call.

## Details

The pseudo-count \\D^{2}\\ statistic quantifies item parameter drift
(IPD) by computing, across quadrature points, the weighted squared
differences between the observed pseudo-proportion of correct responses
and the proportion expected under the item bank parameters, weighted by
the pseudo-counts. The observed pseudo-counts are obtained from the
posterior distribution of each examinee's ability, and the expected
proportions are the probabilities of a correct response computed from
the item bank parameters (Stone, 2000).

The pseudo-count \\D^{2}\\ statistic is calculated as: \$\$ D^{2}\_{PC}
= \sum\_{k=1}^{Q} \left( \frac{r\_{0k} + r\_{1k}}{N}\right) \left(
\frac{r\_{1k}}{r\_{0k} + r\_{1k}} - E\_{1k} \right)^2 \$\$

where \\Q\\ is the number of quadrature points, \\r\_{0k}\\ and
\\r\_{1k}\\ are the pseudo-counts of incorrect and correct responses at
quadrature point \\k\\ (the sums of the examinees' posterior
probabilities at that point over the examinees with each response),
\\E\_{1k}\\ is the probability of a correct response at point \\k\\
under the item bank parameters, and \\N\\ is the number of examinees who
responded to the item.

**Critical Value (`crit.val`)**: The `crit.val` argument specifies the
threshold used to flag an item as exhibiting potential parameter drift.
If an item's pseudo-count \\D^{2}\\ value is greater than or equal to
this threshold, the item is flagged. If `crit.val = NULL`, the function
reports the raw statistic without flagging.

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
  initial pseudo-count \\D^{2}\\ statistics.

- In each subsequent iteration, the item with the highest flagged
  pseudo-count \\D^{2}\\ value is removed from the item set, and the
  statistics are recalculated using only the remaining items.

- The process continues until no additional items are flagged or the
  number of iterations reaches `max.iter`.

- If `max.iter` is reached before convergence, a warning is issued, the
  items flagged at the last iteration are added to the flagged set, and
  `complete` is set to `FALSE`.

- If every item is flagged, the procedure stops with a warning, because
  no item is left to analyze, all items are returned as flagged, and
  `complete` is `FALSE`.

- All flagged items and statistics are saved, and convergence status is
  reported.

This process reduces the influence of already-flagged items on the
posterior distributions used for the remaining items.

## References

Cappaert, K. J., Wen, Y., & Chang, Y. F. (2018). Evaluating CAT-adjusted
approaches for suspected item parameter drift detection. *Measurement:
Interdisciplinary Research and Perspectives, 16*(4), 226-238.
[doi:10.1080/15366367.2018.1511199](https://doi.org/10.1080/15366367.2018.1511199)
.

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Lim, H., & Han, K. T. (2026). IRT residual-based approach to detecting
item parameter drift in CAT. *Journal of Educational and Behavioral
Statistics*.
[doi:10.3102/10769986261460852](https://doi.org/10.3102/10769986261460852)
.

Stone, C. A. (2000). Monte Carlo based null distribution for an
alternative goodness-of-fit test statistic in IRT models. *Journal of
Educational Measurement, 37*(1), 58-75.
[doi:10.1111/j.1745-3984.2000.tb01076.x](https://doi.org/10.1111/j.1745-3984.2000.tb01076.x)
.

## See also

[`ripd()`](https://hwangQ.github.io/irtQ/reference/ripd.md),
[simIPD](https://hwangQ.github.io/irtQ/reference/simIPD.md),
[`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1: No critical value specified
## Compute the pseudo-count D2 statistics for dichotomous items
## Import the "-prm.txt" output file generated by flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Extract metadata for the first 30 3PLM items
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[1:30, 1:6]

# Generate abilities for 500 examinees from N(0, 1)
set.seed(25)
score <- rnorm(500, mean = 0, sd = 1)

# Simulate response data using the item metadata and ability values
data <- simdat(x = x, theta = score, D = 1)

# Compute the pseudo-count D2 statistics (no purification applied)
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
# crit.val = 0.002 is an arbitrary value for illustration; no drift was
# simulated, so any flagged item is a false positive.
# Compute the pseudo-count D2 statistics with purification enabled
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
## -- Example 3: CAT-based IPD detection using simIPD --------------------------
##
## The pseudo-count D2 statistic has no closed-form null distribution.
## Following Cappaert et al. (2018), the critical value is estimated
## empirically via bootstrap:
##   (1) Select items that are not suspected of IPD (anchor items).
##   (2) Compute D2 for the whole item pool from all responses, and take
##       the D2 values of the anchor items as the null distribution.
##   (3) Resample the null distribution with replacement, take the 95th
##       percentile of each bootstrap sample, and average the percentiles
##       to obtain the critical value.
##   (4) Run pcd2() with this critical value.
##
## The simIPD dataset provides:
##   item_par  : original 360-item pool (3PLM, irtQ format)
##   foc_resp  : focal group CAT responses (N = 3,000; TL = 30; ~92% NA)
##   key_item  : indices of 90 highly-exposed key items (IPD analysis targets)
##   item.skip : indices of 270 non-key items (excluded from analysis)
##   ipd_item  : 18 truly drifted items (ground truth; a & b each -0.5)
## -----------------------------------------------------------------------------

data(simIPD)

## -- Step 1. Select anchor items for bootstrap -----------------------------
## Select items that (a) have at least min_resp_anchor observed responses
## and (b) are not suspected of IPD. In practice, exclude the items you
## suspect to have drifted; here the ground truth (simIPD$ipd_item) is
## excluded for illustration.
## The minimum response count follows Lim & Han (2026); it is a practical
## addition to the procedure of Cappaert et al. (2018), who did not set such
## a criterion. In a CAT pool with uneven exposure, items with few responses
## have large D2 values that inflate the tail of the null distribution. For
## the simIPD data, the critical value computed below is about 0.034 when
## all items with at least one response are used, about 0.008 with at least
## 50 responses, about 0.006 with at least 100, and about 0.004 with at
## least 300.
min_resp_anchor <- 300
n_resp          <- colSums(!is.na(simIPD$foc_resp))   # observed responses per item
anchor_items    <- setdiff(which(n_resp >= min_resp_anchor), simIPD$ipd_item)

## -- Step 2. Compute D2 for the whole pool and keep the anchor items -------
## D2 is computed once from all items and responses, so that the posterior
## ability distributions are the same as in the final analysis (Step 4).
## A warning lists the pool items that were never administered.
pcd2_all  <- pcd2(
  x    = simIPD$item_par,
  data = simIPD$foc_resp,
  D    = 1.7
)
#> Warning: Item(s) 17, 19, 22, 24, 27, 30, 32, 34, 44, 47, 50, 54, 56, 65, 73, 76, 77, 82, 86, 93, 102, 105, 108, 111, 112, 120, 122, 131, 137, 140, 144, 155, 156, 164, 165, 174, 176, 182, 193, 196, 201, 206, 219, 222, 223, 225, 232, 241, 255, 265, 269, 271, 278, 279, 284, 286, 289, 291, 292, 293, 295, 305, 313, 315, 316, 339, 342, 343, 344, 345, 352, 353, 359 have no response data.
pcd2_null <- pcd2_all$no_purify$ipd_stat$pcd2[anchor_items]

## -- Step 3. Bootstrap critical value (Cappaert et al., 2018) --------------
## For each bootstrap iteration: resample the null D2 values with replacement
## (the sample size equals the number of anchor items) and take the 95th
## percentile. The critical value is the mean of these percentiles.
## n_boots = 500: reduced for demo speed; use 10,000 in practice for a
##   stable critical value estimate.
set.seed(2024)
n_boots <- 500

crit_val <- mean(vapply(
  X   = seq_len(n_boots),
  FUN = function(b) {
    stats::quantile(sample(x = pcd2_null, size = length(pcd2_null), replace = TRUE),
                    probs = 0.95)
  },
  FUN.VALUE = numeric(1)
))

cat("Bootstrap critical value (alpha = 0.05):", round(crit_val, 6), "\n")
#> Bootstrap critical value (alpha = 0.05): 0.004059 

## -- Step 4. Run PCD2 with the bootstrap critical value + purification --
## The statistics of the first analysis (pcd2_result$no_purify) are the same
## as those of pcd2_all for the analyzed items, since the same responses are
## used. The same warning about the items that were never administered
## appears.
pcd2_result <- pcd2(
  x         = simIPD$item_par,
  data      = simIPD$foc_resp,
  D         = 1.7,
  item.skip = simIPD$item.skip,
  crit.val  = crit_val,
  purify    = TRUE,
  max.iter  = 30
)
#> Warning: Item(s) 17, 19, 22, 24, 27, 30, 32, 34, 44, 47, 50, 54, 56, 65, 73, 76, 77, 82, 86, 93, 102, 105, 108, 111, 112, 120, 122, 131, 137, 140, 144, 155, 156, 164, 165, 174, 176, 182, 193, 196, 201, 206, 219, 222, 223, 225, 232, 241, 255, 265, 269, 271, 278, 279, 284, 286, 289, 291, 292, 293, 295, 305, 313, 315, 316, 339, 342, 343, 344, 345, 352, 353, 359 have no response data.
#> Purification started... 
#>  Iteration: 1 Iteration: 2 Iteration: 3 Iteration: 4 Iteration: 5 Iteration: 6 Iteration: 7 Iteration: 8 Iteration: 9 Iteration: 10 Iteration: 11 Iteration: 12 Iteration: 13 Iteration: 14 Iteration: 15 Iteration: 16 Iteration: 17 Iteration: 18 Iteration: 19 
#> Purification is finished. 

## -- Step 5. Compare detected items to ground truth --------------------
detected <- pcd2_result$with_purify$ipd_item
cat("Truly drifted items (ground truth):", simIPD$ipd_item, "\n")
#> Truly drifted items (ground truth): 26 42 55 58 97 129 138 145 149 181 184 203 242 246 252 310 334 336 
cat("PCD2-detected items:               ", detected, "\n")
#> PCD2-detected items:                26 42 55 58 97 129 138 145 149 181 184 203 236 242 246 252 310 334 336 
cat("True positives:", sum(detected %in% simIPD$ipd_item), "of",
    length(simIPD$ipd_item), "\n")
#> True positives: 18 of 18 
# }

```
