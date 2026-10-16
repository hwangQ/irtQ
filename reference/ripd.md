# Residual-based Item Parameter Drift (RIPD) Detection Framework

This function computes three RIPD statistics, \\RIPD\_{R}\\,
\\RIPD\_{S}\\, and \\RIPD\_{RS}\\, for each item. \\RIPD\_{R}\\ captures
differences in mean raw residuals between groups, which is typically
indicative of uniform item parameter drift (IPD). \\RIPD\_{S}\\ captures
differences in mean squared residuals between groups, reflecting
nonuniform IPD. \\RIPD\_{RS}\\, a combined chi-square-based index, is
sensitive to both uniform and nonuniform IPD.

## Usage

``` r
ripd(x, ...)

# Default S3 method
ripd(
  x,
  data,
  score = NULL,
  group,
  focal.name,
  item.skip = NULL,
  D = 1,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("ripdrs", "ripdr", "ripds"),
  max.iter = 10,
  min.resp = NULL,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'est_irt'
ripd(
  x,
  score = NULL,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("ripdrs", "ripdr", "ripds"),
  max.iter = 10,
  min.resp = NULL,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'est_item'
ripd(
  x,
  group,
  focal.name,
  item.skip = NULL,
  alpha = 0.05,
  missing = NA,
  purify = FALSE,
  purify.by = c("ripdrs", "ripdr", "ripds"),
  max.iter = 10,
  min.resp = NULL,
  method = "ML",
  range = c(-5, 5),
  norm.prior = c(0, 1),
  nquad = 41,
  weights = NULL,
  ncore = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'ripd'
print(
  x,
  digits = max(2L, getOption("digits") - 5L),
  what = c("all", "no_purify", "with_purify"),
  ...
)
```

## Arguments

- x:

  A data frame containing item metadata (e.g., item parameters, number
  of categories, IRT model types, etc.), or an object of class `est_irt`
  or `est_item` obtained from
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) or
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).

  See [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)
  or [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for
  more details about the item metadata. This data frame can be easily
  created using the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function.

- ...:

  Additional arguments passed to
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  when abilities are estimated internally, for example `tol` or `se`.
  The
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  arguments `max.iter` and `missing` cannot be passed this way, because
  the function uses arguments with these names for its own purposes.

- data:

  A matrix of examinees' item responses, with rows for examinees and
  columns for items in the same order as the items in `x`. Responses to
  an item with K score categories must be integers from 0 to K - 1.
  Missing responses are coded as `NA` or, in functions that have a
  `missing` argument, as the value given in `missing`.

- score:

  A numeric vector of examinees' ability estimates (theta values) in the
  same row order as `data`. If `NULL`, the ability estimates are
  computed internally with
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  using `method`, `range`, `norm.prior`, `nquad`, `weights`, and
  `ncore`. When `purify = TRUE`, the ability estimates of all examinees
  in both groups are re-estimated internally at every purification
  iteration with these settings, so a supplied `score` is used only for
  the initial (non-purified) analysis. A missing value in `score`
  excludes the examinee from the analysis, also during purification,
  with a warning. Default is `NULL`.

- group:

  A numeric or character vector indicating the group membership of
  examinees. Its length must equal the number of rows in `data`.
  Examinees whose value equals `focal.name` form the focal group; all
  other examinees are pooled into the reference group (for example, the
  synthetic reference group described in Details). Examinees with `NA`
  in `group` are excluded from both groups.

- focal.name:

  A single numeric or character value specifying the focal group. For
  instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating the
  focal group, set `focal.name = 1`.

- item.skip:

  A numeric vector of item positions (row numbers of `x`), or a logical
  vector with one value per item, identifying the items to exclude from
  the IPD analysis. If `NULL`, all items are analyzed. Skipped items
  still contribute to the ability estimates; only their RIPD statistics
  and moments are set to `NA`, and they are never flagged. In CAT, this
  is typically used to skip low-exposure items (see Details). Default is
  `NULL`.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.702 is
  commonly used for this purpose. Default is 1.

- alpha:

  A numeric value specifying the significance level (\\\alpha\\) for
  hypothesis testing using the RIPD statistics. Default is `0.05`.

- missing:

  A value indicating missing values in the response data set. Default is
  NA.

- purify:

  Logical. Indicates whether to apply a purification procedure. Default
  is `FALSE`.

- purify.by:

  A character string specifying which RIPD statistic is used to perform
  the purification. Available options are "ripdrs" for \\RIPD\_{RS}\\,
  "ripdr" for \\RIPD\_{R}\\, and "ripds" for \\RIPD\_{S}\\. Used only
  when `purify = TRUE`. Default is `"ripdrs"`.

- max.iter:

  A positive integer specifying the maximum number of purification
  iterations. Default is 10. If the limit is reached while flagged items
  remain, a warning is issued, the items flagged at the last iteration
  are added to the flagged set, and the `complete` element of
  `with_purify` is `FALSE`. Because one item is removed per iteration,
  `max.iter` should be at least the number of items expected to drift
  (Lim & Han, 2026, used 80).

- min.resp:

  A single nonnegative number, usually a positive integer, specifying
  the minimum number of valid item responses required from a focal group
  examinee for the examinee's responses to be used. All responses of
  focal group examinees with fewer than `min.resp` (but at least one)
  responses are set to `NA` before the ability estimation, in the
  initial analysis and at every purification iteration, and these
  examinees are excluded from the analysis. A warning reports the number
  of examinees excluded in the initial analysis, and one warning at the
  end of the purification reports the examinees excluded during it. A
  large `min.resp` can leave the focal group without any ability
  estimate in a later purification iteration, which stops the
  purification with a warning. Reference group examinees are not
  affected. If `NULL`, no minimum is applied. Default is `NULL`.

- method:

  A character string indicating the scoring method used to estimate the
  abilities when `score = NULL` and in every purification iteration.
  Available options are:

  - `"ML"`: Maximum likelihood estimation

  - `"MLF"`: Maximum likelihood estimation with fences

  - `"WL"`: Weighted likelihood estimation

  - `"MAP"`: Maximum a posteriori estimation

  - `"EAP"`: Expected a posteriori estimation

  The methods `"EAP.SUM"` and `"INV.TCC"` of
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  are not supported, because they do not give an ability estimate for
  each response pattern. Default is `"ML"`.

- range:

  A numeric vector of length two specifying the lower and upper bounds
  of the ability scale. This is used for the following scoring methods:
  `"ML"`, `"WL"`, and `"MAP"`. Default is `c(-5, 5)`.

- norm.prior:

  A numeric vector of length two specifying the mean and standard
  deviation of the normal prior distribution. These values are used to
  generate the Gaussian quadrature points and weights. Ignored if
  `method` is `"ML"` or `"WL"`. Default is `c(0, 1)`.

- nquad:

  An integer indicating the number of Gaussian quadrature points to be
  generated from the normal prior distribution. Used only when `method`
  is `"EAP"`. Ignored for `"ML"`, `"WL"`, and `"MAP"`. Default is 41.

- weights:

  A two-column matrix or data frame containing the quadrature points (in
  the first column) and their corresponding weights (in the second
  column) for the latent variable prior distribution. The weights and
  points can be conveniently generated using the function
  [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md).

  If `NULL` and `method = "EAP"`, default quadrature values are
  generated based on the `norm.prior` and `nquad` arguments. Ignored if
  `method` is `"ML"`, `"WL"`, or `"MAP"`.

- ncore:

  An integer specifying the number of logical CPU cores used by
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  when abilities are estimated internally. The test statistics
  themselves are not computed in parallel. Default is `1`.

- verbose:

  Logical. If `TRUE`, the progress of the purification procedure is
  printed to the console. Default is `TRUE`.

- digits:

  Not used; kept for consistency with the other print methods.

- what:

  A character string selecting the results to print: `"all"` (default),
  `"no_purify"`, or `"with_purify"`.

## Value

This function returns an object of class `"ripd"`, a list with the
following elements:

- no_purify:

  A list of sub-objects containing the results of IPD analysis without
  applying a purification procedure. The sub-objects include:

  ipd_stat

  :   A data frame of RIPD results for all items with the columns `id`,
      `ripdr`, `z.ripdr` (standardized \\RIPD\_{R}\\), `ripds`,
      `z.ripds` (standardized \\RIPD\_{S}\\), `ripdrs`, `p.ripdr`,
      `p.ripds`, `p.ripdrs`, `n.ref`, `n.foc`, and `n.total`. `n.ref`
      and `n.foc` are the numbers of examinees in each group who
      responded to the item and have an ability estimate. Statistics and
      p-values are rounded to four decimal places, and items are flagged
      using the unrounded p-values; they are `NA` for items in
      `item.skip`. \\RIPD\_{RS}\\ has no standardized value because it
      is a \\\chi^{2}\\-based statistic.

  moments

  :   A data frame of the null means and standard deviations of the RIPD
      statistics with the columns `id`, `mu.ripdr`, `sigma.ripdr`,
      `mu.ripds`, `sigma.ripds`, and `covariance` (the covariance
      between \\RIPD\_{R}\\ and \\RIPD\_{S}\\). Values are `NA` for
      items in `item.skip`.

  ipd_item

  :   A list with the elements `ripdr`, `ripds`, and `ripdrs`, each
      giving the indices (row positions in `x`) of the items flagged by
      the corresponding statistic at level `alpha`. An element is `NULL`
      when no item is flagged.

  score

  :   A numeric vector of ability estimates used to compute the RIPD
      statistics.

- purify:

  A logical value indicating whether the purification procedure was
  applied.

- with_purify:

  A list of sub-objects containing the results of IPD analysis with a
  purification procedure. All elements are `NULL` when `purify = FALSE`.
  The sub-objects include:

  purify.by

  :   A character string indicating the RIPD statistic used for
      purification. Possible values are "ripdr", "ripds", and "ripdrs",
      corresponding to \\RIPD\_{R}\\, \\RIPD\_{S}\\, and \\RIPD\_{RS}\\,
      respectively.

  ipd_stat

  :   A data frame of the RIPD results for all items with the same
      columns as in `no_purify` plus `n.iter`, the analysis from which
      each row is taken. For a flagged item, `n.iter` is the iteration
      at which it was flagged (0 is the initial analysis); for the other
      items, it is the final iteration.

  moments

  :   A data frame of the null moments for all items with the same
      columns as in `no_purify` plus `n.iter`, taken from the same
      analysis as the corresponding row of `ipd_stat`.

  ipd_item

  :   A numeric vector of item indices (row positions in `x`) identified
      as IPD items across all purification iterations, sorted in
      ascending order. If `max.iter` is reached before convergence or
      the purification stops because a group has no ability estimate
      left, the items flagged at the last iteration (iteration `n.iter`)
      are also included. It contains all items if every item is flagged.
      `NULL` if no item is flagged in the initial analysis.

  n.iter

  :   An integer indicating the total number of iterations performed
      during the purification process. The iteration that removes the
      last item, or after which a group has no ability estimate left, is
      not counted, because the statistics are not recomputed.

  score

  :   A numeric vector of purified ability estimates used to compute the
      final RIPD statistics, from the last iteration counted in
      `n.iter`. `NULL` when no item is flagged in the first analysis,
      because purification is not carried out.

  complete

  :   A logical value indicating whether the purification process
      converged. `FALSE` if the procedure stopped because the maximum
      number of iterations was reached, because every item was flagged,
      or because a group had no examinee with an ability estimate left.

- alpha:

  A numeric value indicating the significance level (\\\alpha\\) used in
  hypothesis testing for RIPD statistics.

- call:

  The matched function call.

## Details

The current version supports dichotomous items only. The function stops
if `x` contains GRM or GPCM items or if `data` contain response values
other than 0 and 1.

**Theoretical Background: From RDIF to RIPD**

The RIPD framework directly adapts the residual-based differential item
functioning (RDIF) detection framework (Lim et al., 2022; Lim & Choe,
2023) to detect item parameter drift (IPD) in computerized adaptive
testing (CAT). Each RIPD statistic (\\RIPD_R\\, \\RIPD_S\\,
\\RIPD\_{RS}\\) mirrors its RDIF counterpart (\\RDIF_R\\, \\RDIF_S\\,
\\RDIF\_{RS}\\) in both computation and asymptotic theory, with two key
adaptations. First, whereas the original RDIF framework estimates item
parameters once from the pooled data of both groups, the RIPD framework
uses the *original, pre-calibrated* item parameters from the CAT pool
directly, eliminating the need for item recalibration. Second, the
reference group is a synthetic, drift-free group simulated from the
focal group's ability estimates (see CAT-Specific Workflow). This makes
RIPD practically scalable for operational CAT programs, where the
response data are sparse and recalibration is often impractical.

IPD in CAT can be viewed as a special case of DIF across time (Lord,
1980; Veerkamp & Glas, 2000): the focal group represents current
test-takers whose responses may reflect drift, while the reference group
represents a baseline under invariant item parameters.

**Residual Definition**

For examinee \\i\\ and item \\j\\, the raw residual is defined as
\$\$e\_{ij} = u\_{ij} - P_j(\hat{\theta}\_i),\$\$ where \\u\_{ij}\\ is
the observed binary response (0 or 1) and \\P_j(\hat{\theta}\_i)\\ is
the model-predicted probability of a correct response under the
*original* item parameters evaluated at the examinee's ability estimate
\\\hat{\theta}\_i\\. Because item responses are independent but not
identically distributed (due to varying predicted probabilities across
examinees with different ability levels), the asymptotic distributions
of the RIPD statistics are established via Lyapunov's central limit
theorem. Below, \\P\_{ij} = P_j(\hat{\theta}\_i)\\, and \\N_F\\ and
\\N_R\\ denote the numbers of focal and reference group examinees who
responded to item \\j\\; all sums run over these examinees only.

**Three RIPD Statistics and Their Asymptotic Distributions**

*RIPD\\\_R\\* (uniform drift): Measures the difference in mean raw
residuals between the focal (\\F\\) and reference (\\R\\) groups,
\$\$RIPD\_{R,j} = \frac{\sum\_{i \in F} e\_{ij}}{N_F} - \frac{\sum\_{i
\in R} e\_{ij}}{N_R}.\$\$ \\RIPD\_{R,j}\\ is most effective when the
item response functions (IRFs) of the two groups differ primarily in
location (e.g., drift in the difficulty \\b\\ or guessing \\c\\
parameter). Under the null hypothesis \\H_0\\ of no IPD, \\RIPD\_{R,j}\\
asymptotically follows a normal distribution with mean zero and
analytically derived variance \\\sigma^2\_{RIPD_R} = \frac{\sum\_{i \in
F} P\_{ij}(1 - P\_{ij})}{N_F^2} + \frac{\sum\_{i \in R} P\_{ij}(1 -
P\_{ij})}{N_R^2}\\. A standard Z-test is then applied as \\Z_R =
RIPD\_{R,j} / \sigma\_{RIPD_R}\\.

*RIPD\\\_S\\* (nonuniform drift): Measures the difference in mean
squared residuals between groups, \$\$RIPD\_{S,j} = \frac{\sum\_{i \in
F} e\_{ij}^2}{N_F} - \frac{\sum\_{i \in R} e\_{ij}^2}{N_R}.\$\$
\\RIPD\_{S,j}\\ captures drift that varies across the ability continuum,
such as a change in the discrimination parameter \\a\\, which
\\RIPD\_{R,j}\\ alone may not detect. Importantly, under \\H_0\\, the
null mean of \\RIPD\_{S,j}\\ is **not zero** in general; it equals
\$\$\mu\_{RIPD_S} = \frac{\sum\_{i \in F} P\_{ij}(1 - P\_{ij})}{N_F} -
\frac{\sum\_{i \in R} P\_{ij}(1 - P\_{ij})}{N_R},\$\$ which depends on
the predicted probabilities of the examinees in each group who responded
to item \\j\\. Even with a synthetic reference group, the ability
estimates and the administered items of the reference group differ from
those of the focal group, so this null mean is generally nonzero and is
computed item by item. The asymptotic variance is also analytically
derived. A standard Z-test is applied as \\Z_S = (RIPD\_{S,j} -
\mu\_{RIPD_S}) / \sigma\_{RIPD_S}\\.

*RIPD\\\_{RS}\\* (combined): A Wald-type statistic that jointly tests
for uniform and nonuniform drift within a single framework,
\$\$RIPD\_{RS,j} = (\boldsymbol{\nu}\_j - \boldsymbol{\mu}\_j)^\top
\hat{\boldsymbol{\Sigma}}\_j^{-1} (\boldsymbol{\nu}\_j -
\boldsymbol{\mu}\_j),\$\$ where \\\boldsymbol{\nu}\_j = (RIPD\_{R,j},\\
RIPD\_{S,j})^\top\\, \\\boldsymbol{\mu}\_j = (\mu\_{RIPD_R},\\
\mu\_{RIPD_S})^\top\\ is the vector of null means, and
\\\hat{\boldsymbol{\Sigma}}\_j\\ is the analytically derived \\2 \times
2\\ covariance matrix of \\(RIPD\_{R,j}, RIPD\_{S,j})\\. Under \\H_0\\,
\\RIPD\_{RS,j}\\ asymptotically follows a \\\chi^2\\ distribution with 2
degrees of freedom. \\RIPD\_{RS}\\ is the **recommended** statistic in
practice because it is sensitive to both uniform and nonuniform drift
and addresses potential inflation of the family-wise Type I error that
arises from applying \\RIPD_R\\ and \\RIPD_S\\ separately to the same
item.

In rare cases, the covariance matrix of \\RIPD_R\\ and \\RIPD_S\\ is
singular or nearly singular (its reciprocal condition number is 1e-10 or
less), for example when all examinees who responded to an item have the
same probability of a correct response, so that \\RIPD_S\\ is a linear
function of \\RIPD_R\\. Then \\RIPD\_{RS}\\ is computed with the
Moore-Penrose generalized inverse of the covariance matrix and compared
with a chi-square distribution whose degrees of freedom equal the rank
of the matrix (Moore, 1977), and a warning names the item. This handling
is not part of the original method (Lim et al., 2022; Lim & Han, 2026)
and is added in irtQ. In most data the matrix is not singular, and
\\RIPD\_{RS}\\ has two degrees of freedom.

**Diagnosing the Nature of Drift**

The pattern of flagging across the three statistics can help diagnose
the type of drift. If an item is flagged by \\RIPD_R\\ (and possibly
\\RIPD\_{RS}\\) but not \\RIPD_S\\, uniform drift is indicated (e.g.,
difficulty shift). If flagged by \\RIPD_S\\ (and possibly
\\RIPD\_{RS}\\) but not \\RIPD_R\\, nonuniform drift is suggested (e.g.,
discrimination change). Items flagged by all statistics or by
\\RIPD\_{RS}\\ alone likely exhibit mixed drift.

**CAT-Specific Workflow**

The RIPD procedure for CAT consists of three steps (Lim & Han, 2026):

1.  **Focal group CAT**: The current cohort of examinees takes the CAT
    using the operational (potentially drifted) item pool. Each examinee
    receives only a subset of items; the resulting response matrix is
    sparse. The final ability estimates \\\hat{\theta}\\ (e.g., ML
    estimates) from this step are retained.

2.  **Synthetic reference group construction**: The focal group's
    \\\hat{\theta}\\ values are treated as true abilities. A CAT
    simulation with the same item pool and algorithm settings as in step
    1 is run, generating item responses from the *original*
    (pre-calibrated, drift-free) item parameters. This creates a
    synthetic reference group that matches the focal group's ability
    distribution, so that any systematic differences in residuals can be
    attributed to IPD rather than ability confounds. The reference group
    size is controlled by a replication factor: e.g., 1F reuses the
    focal \\\hat{\theta}\\ values once (N_ref = N_foc), while kF
    replicates them k times. Larger reference groups reduce sampling
    variability in the reference residuals and improve detection power.
    In the simulation study of Lim and Han (2026), power gains
    diminished beyond 5F, and a reference group of at least 8 to 10
    times the focal group size (8F to 10F) was suggested as a practical
    guideline. Note that `ripd()` does not simulate the reference group;
    this step must be carried out with a CAT simulation program outside
    irtQ, and the simulated responses and ability estimates are then
    combined with the focal group data.

3.  **RIPD computation with purification**: Residuals are computed for
    both groups using the original item parameters, and the three RIPD
    statistics are evaluated item-by-item. Items whose Z-test (for
    \\RIPD_R\\ or \\RIPD_S\\) or chi-square test (for \\RIPD\_{RS}\\)
    exceeds the critical value at level `alpha` are flagged as drifting.

**Purification Procedure**

When `purify = TRUE`, an iterative purification procedure adapted from
Lim et al. (2022) is applied to mitigate the bias in ability estimates
caused by drifted items. This is analogous to the contaminating effect
of DIF items on matching variables. At each iteration:

1.  Among the remaining items flagged by the `purify.by` statistic, the
    item with the smallest p-value (compared on the log scale) is
    removed and added to the flagged set. For \\RIPD_R\\ and \\RIPD_S\\,
    this is the item with the largest absolute z value.

2.  Ability estimates of all examinees in both groups are recomputed
    with
    [`est_score`](https://hwangQ.github.io/irtQ/reference/est_score.md),
    excluding all currently flagged items and using the `method`,
    `range`, `norm.prior`, `nquad`, `weights`, and `ncore` arguments.

3.  RIPD statistics are recalculated using the updated ability estimates
    for the remaining items.

The process continues until no additional items are flagged
(convergence) or `max.iter` is reached. In the latter case, a warning is
issued, the items flagged at the last iteration are added to the flagged
set, and `complete` is set to `FALSE`. If every item is flagged, the
procedure stops with a warning, because no item is left to estimate the
abilities, all items are returned as flagged, and `complete` is `FALSE`.
The procedure also stops with a warning and `complete = FALSE` when,
after an item is removed, no examinee of the reference group or of the
focal group has an ability estimate left (see `min.resp`); the remaining
items then keep the results of the previous iteration, and the items
flagged in that iteration are added to the flagged items. Examinees
without any response have no ability estimate and are excluded. During
purification, an examinee whose remaining responses are all removed with
the flagged items, or who has fewer than `min.resp` responses to the
remaining items, is also excluded, and one warning at the end reports
the number of such examinees. The statistic used to drive purification
is specified by `purify.by`; `"ripdrs"` (\\RIPD\_{RS}\\) is recommended
as it is sensitive to both types of drift.

**Key Items and item.skip**

In CAT contexts, many items in the pool receive few responses and thus
lack sufficient data for reliable IPD analysis. The `item.skip` argument
excludes such items (e.g., low-exposure non-key items) from the IPD
tests. Skipped items are still used to compute the ability estimates. In
practice, key items for RIPD evaluation are typically identified through
a preliminary CAT simulation as those with the highest average exposure
frequencies under the intended operational settings.

**Printing**

The [`print()`](https://rdrr.io/r/base/print.html) method of a `"ripd"`
object shows the flagged items and the RIPD statistics of the analyses
without and with purification. Its argument `what` selects `"no_purify"`
or `"with_purify"` (default `"all"`). Items without statistics (skipped
items and items without responses in a group) are not shown in the
tables.

## Methods (by class)

- `ripd(default)`: Default method for a data frame `x` containing item
  metadata.

- `ripd(est_irt)`: Method for an object of class `est_irt` created by
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md). The
  response data and the scaling factor `D` are taken from the object.

- `ripd(est_item)`: Method for an object of class `est_item` created by
  [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md).
  The response data, ability estimates, and the scaling factor `D` are
  taken from the object.

## Methods (by generic)

- `print(ripd)`: Print method for an object of class `ripd`. It prints
  the flagged items and the RIPD statistics of the analyses without and
  with purification (see **Printing** in Details).

## References

Lim, H., & Choe, E. M. (2023). Detecting differential item functioning
in CAT using IRT residual DIF approach. *Journal of Educational
Measurement, 60*(4), 626-650.
[doi:10.1111/jedm.12366](https://doi.org/10.1111/jedm.12366) .

Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
item functioning detection framework in item response theory. *Journal
of Educational Measurement, 59*(1), 80-104.
[doi:10.1111/jedm.12313](https://doi.org/10.1111/jedm.12313) .

Lim, H., & Han, K. T. (2026). IRT residual-based approach to detecting
item parameter drift in CAT. *Journal of Educational and Behavioral
Statistics*.
[doi:10.3102/10769986261460852](https://doi.org/10.3102/10769986261460852)
.

Lord, F. M. (1980). Applications of item response theory to practical
testing problems. Lawrence Erlbaum Associates.

Moore, D. S. (1977). Generalized inverses, Wald's method, and the
construction of chi-squared tests of fit. *Journal of the American
Statistical Association, 72*(357), 131-137.
[doi:10.1080/01621459.1977.10479921](https://doi.org/10.1080/01621459.1977.10479921)
.

Veerkamp, W. J. J., & Glas, C. A. W. (2000). Detection of known items in
adaptive testing with a statistical quality control method. *Journal of
Educational and Behavioral Statistics, 25*(4), 373-389.
[doi:10.3102/10769986025004373](https://doi.org/10.3102/10769986025004373)
.

## See also

[`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md),
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md),
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md),
[`pcd2()`](https://hwangQ.github.io/irtQ/reference/pcd2.md),
[simIPD](https://hwangQ.github.io/irtQ/reference/simIPD.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## --- RIPD Example: Detecting IPD in CAT ---------------------------------
##
## Background (Lim & Han, 2026):
##   In CAT-based IPD detection using RIPD, the reference group is
##   "synthetic": it is created by re-administering a CAT to examinees whose
##   true abilities are set equal to the focal group's ML theta estimates,
##   using the ORIGINAL (pre-drift) item parameters. This eliminates the
##   need for recalibration and makes RIPD directly applicable to
##   operational CAT settings.
##
## The simIPD dataset contains:
##   - foc_resp / foc_score : focal group CAT responses + ML theta estimates
##                            (IPD items: a and b each drifted by -0.5)
##   - ref_resp / ref_score : synthetic reference group CAT responses + ML
##                            theta estimates (original parameters, 1F size)
##   - item_par  : original (non-drifted) 360-item pool in irtQ format
##   - key_item  : indices of 90 key items (highly exposed in CAT)
##   - item.skip : indices of 270 non-key items (excluded from RIPD)
##   - ipd_item  : indices of 18 truly drifted items (ground truth)
## -------------------------------------------------------------------------

data(simIPD)

## Step 1. Combine focal and synthetic reference group data
##         (the row order does not matter as long as group matches the rows)
data  <- rbind(simIPD$ref_resp,  simIPD$foc_resp)
score <- c(simIPD$ref_score,     simIPD$foc_score)
group <- c(rep(0, nrow(simIPD$ref_resp)),   # 0 = reference
           rep(1, nrow(simIPD$foc_resp)))   # 1 = focal

## Step 2. Run RIPD without purification
##         item.skip excludes the 270 non-key items from analysis
ripd_result <- ripd(
  x          = simIPD$item_par,
  data       = data,
  score      = score,
  group      = group,
  focal.name = 1,           # focal group is coded as 1
  item.skip  = simIPD$item.skip,
  D          = 1.7,
  alpha      = 0.05,
  purify     = FALSE,
  method     = "ML",
  range      = c(-5, 5)
)

## Step 3. Compare the items flagged by RIPD_RS to the ground truth
detected <- ripd_result$no_purify$ipd_item$ripdrs
cat("Truly drifted items (ground truth):", simIPD$ipd_item, "\n")
#> Truly drifted items (ground truth): 26 42 55 58 97 129 138 145 149 181 184 203 242 246 252 310 334 336 
cat("RIPD_RS-flagged items:             ", detected, "\n")
#> RIPD_RS-flagged items:              26 55 58 97 129 130 138 145 149 181 184 203 215 231 236 242 246 252 277 310 312 334 336 351 354 
cat("True positives:", sum(detected %in% simIPD$ipd_item), "of",
    length(simIPD$ipd_item), "\n")
#> True positives: 17 of 18 

## Step 4. Purification (not run: it takes several minutes because the
##         abilities of the examinees are re-estimated at every iteration)
##
##   ripd_purified <- ripd(
##     x = simIPD$item_par, data = data, score = score, group = group,
##     focal.name = 1, item.skip = simIPD$item.skip, D = 1.7, alpha = 0.05,
##     purify = TRUE, purify.by = "ripdrs", max.iter = 30,
##     method = "ML", range = c(-5, 5)
##   )
##   print(ripd_purified, what = "with_purify")

## -- Note on reference group size -----------------------------------------
## This example uses a 1F reference group (n_ref = n_foc = 3,000).
## In practice, a larger reference group improves detection power;
## Lim and Han (2026) suggest at least 8F to 10F. To create, for example,
## an 8F reference group, replicate the focal ability estimates and run
## the same CAT simulation with the original item parameters using a CAT
## simulation program (irtQ does not provide one):
##
##   theta_8F <- rep(simIPD$foc_score, times = 8)   # 24,000 examinees
##   # simulate a CAT for theta_8F with simIPD$item_par, then call ripd()
##   # with the simulated responses and final ability estimates
## -------------------------------------------------------------------------
```
