# Simulated CAT Data for Item Parameter Drift (IPD) Detection

A simulated dataset designed to illustrate the use of the
[`ripd`](https://hwangQ.github.io/irtQ/reference/ripd.md) function for
detecting item parameter drift (IPD) in computerized adaptive testing
(CAT). The dataset represents one replication of a CAT-based IPD
simulation study in which 5% of a 360-item pool (18 items) had both
their discrimination (*a*) and difficulty (*b*) parameters drift
downward by 0.5, under a 30-item adaptive test where all focal group
examinees were exposed to the drifted items (100% exposure rate).

The data reflect a workflow of IPD detection using the residual-based
IPD (RIPD) framework (Lim & Han, in press): a focal group of examinees
takes a CAT using a drifted item pool, and a synthetic reference group
is created by re-running the CAT with the focal group's ability
estimates as true abilities but with the original (non-drifted) item
parameters. RIPD statistics are then used to detect which items have
drifted between the two groups.

## Usage

``` r
simIPD
```

## Format

A named list with eight elements:

- item_par:

  A data frame with 360 rows and 6 columns containing the **original**
  (pre-drift) item parameters in `irtQ` format (columns: `id`, `cats`,
  `model`, `par.1`, `par.2`, `par.3`). All items follow the 3PLM.

- key_item:

  An integer vector of length 90 giving the row indices (in `item_par`)
  of the *key items* - items selected by a preliminary CAT simulation as
  highly exposed and therefore most relevant for IPD analysis.

- item.skip:

  An integer vector of length 270 giving the row indices of non-key
  items that should be excluded from RIPD analysis (i.e., the complement
  of `key_item` in `1:360`). Pass this vector to the `item.skip`
  argument of
  [`ripd()`](https://hwangQ.github.io/irtQ/reference/ripd.md).

- ipd_item:

  An integer vector of length 18 giving the row indices of items that
  were subjected to IPD manipulation. These 18 items (5% of the 360-item
  pool) had both their discrimination (*a*) and difficulty (*b*)
  parameters decreased by 0.5. All 18 items are members of `key_item`.

- foc_resp:

  An integer matrix of dimensions 3000 x 360 containing the **focal
  group** CAT response data. Each row is one examinee; each column
  corresponds to an item in `item_par`. Because CAT administers only 30
  items per examinee, approximately 92% of entries are `NA`. Responses
  were generated using the *drifted* item parameters (all focal
  examinees were exposed to IPD items; exposure rate = 100%).

- foc_score:

  A numeric vector of length 3000 containing the **focal group** final
  maximum likelihood (ML) theta estimates obtained from the CAT.

- ref_resp:

  An integer matrix of dimensions 3000 x 360 containing the **synthetic
  reference group** CAT response data (same sparsity structure as
  `foc_resp`). The reference group was constructed by: (1) using
  `foc_score` as true ability values (1F scaling, i.e., the reference
  group has the same size as the focal group); (2) generating item
  responses from the *original* (non-drifted) item parameters; and (3)
  running an independent CAT simulation. This synthetic reference group
  mirrors the construction described in Lim & Han (in press).

- ref_score:

  A numeric vector of length 3000 containing the **synthetic reference
  group** final ML theta estimates.

## Details

**Simulation conditions:**

- Item pool: 360 three-parameter logistic model (3PLM) items

- Test length: 30 items per examinee

- Item selection: Maximum Fisher Information (MFI) with target exposure
  rate 0.30

- Focal group: \\n = 3{,}000\\; true abilities drawn from \\N(0, 1)\\

- Reference group: \\n = 3{,}000\\ (1F); true abilities = `foc_score`

- IPD items: 18 (5% of 360), randomly drawn from `key_item`

- IPD manipulation: both *a* and *b* decreased by 0.50

- Interim scoring: EAP; final scoring: ML

- Scaling constant: \\D = 1.7\\

- Random seed: 2024

**Note on reference group size:** A 1F reference group (same size as the
focal group) is used here for compactness. In practice, larger synthetic
reference groups (e.g., 3F - 8F) are recommended to improve RIPD
detection power (Lim & Han, in press). A larger reference group can be
created by replicating the focal theta estimates: e.g.,
`rep(foc_score, times = 3)` for a 3F group, then re-running the CAT
simulation with the original item parameters.

## References

Lim, H., & Han, K. T. (in press). IRT residual-based approach to
detecting item parameter drift in CAT. *Journal of Educational and
Behavioral Statistics*.

## See also

[`ripd`](https://hwangQ.github.io/irtQ/reference/ripd.md),
[`pcd2`](https://hwangQ.github.io/irtQ/reference/pcd2.md),
[`simCAT_DC`](https://hwangQ.github.io/irtQ/reference/simCAT_DC.md),
[`simCAT_MX`](https://hwangQ.github.io/irtQ/reference/simCAT_MX.md)

## Examples

``` r
data(simIPD)
str(simIPD, max.level = 1)
#> List of 8
#>  $ item_par :'data.frame':   360 obs. of  6 variables:
#>  $ key_item : num [1:90] 2 3 10 13 14 15 18 20 23 25 ...
#>  $ item.skip: int [1:270] 1 4 5 6 7 8 9 11 12 16 ...
#>  $ ipd_item : num [1:18] 26 42 55 58 97 129 138 145 149 181 ...
#>  $ foc_resp : num [1:3000, 1:360] NA NA NA NA NA NA NA NA NA NA ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ foc_score: num [1:3000] -0.667 0.771 -1.857 -1.003 1.934 ...
#>  $ ref_resp : num [1:3000, 1:360] NA NA NA NA NA NA NA NA NA NA ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ ref_score: num [1:3000] -0.452 0.413 -1.885 -1.33 1.79 ...

# Item parameter data frame (first 6 rows)
head(simIPD$item_par)
#>   id cats model par.1  par.2 par.3
#> 1 I1    2  3PLM 0.928 -0.606 0.218
#> 2 I2    2  3PLM 1.256 -0.562 0.173
#> 3 I3    2  3PLM 1.815  0.384 0.115
#> 4 I4    2  3PLM 1.128 -0.382 0.253
#> 5 I5    2  3PLM 0.873 -1.458 0.198
#> 6 I6    2  3PLM 0.790  2.210 0.121

# Focal group response matrix (sparse)
dim(simIPD$foc_resp)
#> [1] 3000  360
mean(is.na(simIPD$foc_resp))  # ~0.92 (92% NA due to CAT)
#> [1] 0.9166667
```
