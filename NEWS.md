
<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# irtQ 1.4.0

## Bug Fixes

- `ctt()` returned a standardized alpha (`alpha_std`) that was too
  large, and could exceed 1, when a test included a constant item. The
  constant item is now kept in the item count as in the raw alpha, and
  the value is 0 when only one item varies, as is the raw alpha. An
  undefined standardized alpha or mean discrimination is `NA` instead of
  `NaN`.

- `ctt()`, `ctt_distr()`, and `freq_score()` silently left out of their
  counts a score within rounding error of a whole number, such as
  `(0.1 + 0.2) * 10`. Such a score is now counted as that whole number.

- `freq_score()` tabulated a factor by its level codes (1, 2, 3, ...)
  instead of its labels. A factor is now converted through its labels.

- `ctt_distr()` listed an option twice in selected-response mode when
  the responses wrote it in different ways, such as `"1"` and `"01"`,
  which inflated the frequencies and percentages. The option is now
  listed once.

- `score_resp()` and `ctt_distr()` did not treat a response with
  surrounding spaces, such as `" 9"`, as missing when `missing = "9"`.
  Responses are now compared with `missing` after trimming spaces.

- `score_resp()` and `ctt_distr()` silently scored every response to an
  item as wrong when its key value contained a comma, such as `"1,5"`.
  They now stop with an error, because every item must have exactly one
  correct option.

- `ctt()` did not flag an item whose discrimination is undefined, such
  as an item with a constant score, so the item could pass the
  `crit.dis` check. Such an item is now flagged with
  `"discrimination undefined"`.

- `shape_df_fipc()` misplaced the parameter columns of the new items
  when the fixed items included a polytomous item and every new item was
  dichotomous, which changed a threshold of a fixed polytomous item. The
  new items now follow the columns of the fixed items. FIPC results
  built from such metadata should be rerun.

- The GPCM probabilities of `prm()`, `traceline()`, and `info()` were
  distorted when an exponent exceeded 700, and the values at the other
  thetas and items of the same call changed with them. Each row is now
  rescaled by its own maximum. Results in the usual range change by less
  than 1e-12, and `prm()` no longer returns a GPCM probability of
  exactly 0.

- `info()` returned `NaN` instead of 0 for a dichotomous item at a theta
  far from the item location.

- `traceline()` returned a vector instead of a one-row matrix in
  `prob.cats` for a polytomous item when a single theta was given, which
  made `plot()` fail.

- `drm()` checked each probability against the guessing parameter of
  another item, and the expected item score of `traceline()` included
  categories that the item does not have. Probabilities and expected
  scores could differ by about 1e-10, and item parameter estimates by
  about 1e-9.

- `traceline()` and `info()` stopped with an error for a two-category
  GRM or GPCM item, such as one made by `bring.mirt()`, and for a 3PLM
  or DRM item with an `NA` guessing parameter, for which `drm()`
  returned `NA`. An `NA` guessing parameter is now treated as 0, so a
  two-category GRM or GPCM item is computed as a 2PLM item, also when
  the test has polytomous items of the same model.

- `shape_df()` used a single `g` value for the first dichotomous item
  only and set the guessing parameters of the other items to 0. A single
  `g` value is now used for every dichotomous item, as in `simdat()`,
  and a `g` vector of any other wrong length stops with an error.

- `est_score()` with `ncore > 1` failed when a chunk of examinees had a
  single row, for example with 3 examinees and `ncore = 2`, and it did
  not close the cluster when a worker failed. Each chunk is now kept as
  a matrix, also when it has one row, and the cluster is always closed.

- `est_score()` with `method = "ML"`, `"WL"`, `"MLF"`, or `"MAP"` could
  stop at a point that was not a solution when the iterations alternated
  between two values, which happened mostly for low-scoring examinees on
  tests with 3PLM items. The estimate then depended on `max.iter`. When
  the iterations do not converge, the root of the score function is now
  found within the range of the visited values.

- `est_score()` with `method = "INV.TCC"` returned wrong standard errors
  when `intpol = FALSE` or when the interpolation was skipped with a
  warning (a `range.tcc` that is too narrow), because the estimates were
  paired with the wrong sum scores of the score distribution. The
  standard deviation is now taken over the sum scores that have an
  estimate. Results with `intpol = TRUE` and no warning are unchanged.

- `est_score()` with `method = "MLF"` gave an examinee with all
  responses missing an estimate near 0 and a standard error of 316.2
  without a warning. The examinee now gets `NA` and the warning, as with
  the other methods.

- `est_score()` with `method = "EAP"` returned `NaN`, with a warning
  about missing responses, for tests of about 1,000 items or more,
  because the product of the item probabilities underflowed. The
  posterior is now computed from log-likelihoods, and results for usual
  tests change by less than 1e-14.

- `llike_score()` with `method = "MLF"` stopped with an error under the
  default `fence.b = NULL`. The fences are now placed at -5 and 5, the
  default `range` of `est_score()`. The function also merged examinees
  whose rows had the same row name, and an examinee without any
  response, who stopped the function (`"ML"`, `"MAP"`) or got values
  from the fence items alone (`"MLF"`), now gets `NA`.

- `lwrc()` stopped with an error for a single `theta` when the test had
  polytomous items.

- `sx2_fit()` counted some summed score groups twice in the tables of
  polytomous items when the test was so short that the groups merged at
  the two ends overlapped (for example, two items with five categories).
  The observed and expected frequencies and the degrees of freedom were
  too large. Each summed score group is now counted once.

- `sx2_fit()` and `irtfit()` reported a p-value of 0 (or `NaN` with a
  warning) for an item without degrees of freedom left after collapsing,
  which read as a strong misfit. The critical value and p-value of such
  an item are now `NA`, with a warning that names the item.

- `sx2_fit()` named the columns of the observed frequency and proportion
  tables (`obs_freq` and `obs_prop`) with an `exp_freq.` prefix
  (`exp_freq.score.0`, ...) when the test had a single polytomous item.
  The columns are now named `score.0`, `score.1`, ... as for tests with
  several polytomous items. The values are unchanged.

- `plot.irtfit()` ignored `xlab.text` when `type = "both"` and always
  used theta as the x-axis title. The given title is now used, as for
  the other types.

- `plot.irtfit()` drew the standardized residuals as blue crosses,
  instead of red circles, when every plotted residual exceeded `overSR`.
  Each kind of point now keeps its color and shape.

- `irtfit()` stopped with an error for a logical response matrix, and
  `irtfit()` and `sx2_fit()` stopped with an error when the groups of an
  item were merged into a single row (for example, with a large
  `min.collapse`). These cases now work; an item left without degrees of
  freedom gets `NA` for the critical value and p-value of the chi-square
  statistic, with a warning.

- `sx2_fit()` for an `est_irt` object computed the expected frequencies
  with a standard normal distribution even when a different latent
  distribution was used or estimated in the calibration (for example,
  with `EmpHist = TRUE` or FIPC). When `weights`, `norm.prior`, and
  `nquad` are all omitted, the latent distribution stored in the object
  is now used. The statistics change slightly for a calibration with a
  standard normal distribution and can change enough to alter the
  decision at `alpha` when the distribution was estimated. Set
  `norm.prior = c(0, 1)` to reproduce the earlier results.

- `run_mst()` stopped with an obscure error for a panel whose routing
  module is not module 1 when `ini_mod` was `NULL`, because the starting
  module was drawn with `sample()` from a single module number. The
  examinees now start in the routing module.

- `reval_mst()` took the score distribution of the first stage from
  module 1 even when the routing module had another number, which gave
  wrong results without an error. The routing module is now taken from
  `route_map`, and a first stage with more than one module stops with an
  error. `run_mst()` still runs such panels.

- `cac_lee()` and `cac_rud()` used weights that did not sum to 1 as
  given, so the marginal indices could exceed 1. Such weights are now
  rescaled to sum to 1, with a warning. Weights that are missing,
  infinite, or negative, or that have a non-positive sum, now stop with
  an error.

- `plot.find_cut()` with `theta_range` still stretched the x-axis to the
  crossing lines outside the range. The x-axis now covers the requested
  range.

- `ripd()` stopped with an obscure error when an examinee had no ability
  estimate: a missing value in `score`, an examinee without any
  response, a focal group examinee left out by `min.resp` (so `min.resp`
  could not be used), or an examinee whose responses were all removed
  during purification. Such examinees are now excluded, with a warning
  for the first and the third case, and the function stops when no
  examinee has an ability estimate. `min.resp` now also applies to the
  first analysis when `score` is given.

- `ripd()` and `pcd2()` stopped with an obscure error when purification
  left a single item. They now work, and stop with a clear message when
  every item is flagged.

- `ripd()` and `pcd2()` silently accepted responses outside the score
  categories (for example, an unrecorded missing code such as -9, which
  gave chi-square values in the hundred thousands), a `group` of the
  wrong length, and an `item.skip` that was not a set of item positions
  (an ID given to `pcd2()` added a spurious row and skipped nothing).
  They now stop with a clear message, a logical `item.skip` is read as
  item positions, also during purification, and character responses are
  read as numbers.

- `pcd2()` warned about excluded items even when `min.resp` excluded
  none, and worded the warning poorly. It now warns only when items are
  excluded.

- `ripd()` compared the p-values rounded to four decimals with `alpha`,
  so an item with a p-value just above `alpha` (for example, 0.05004 for
  `alpha = 0.05`) was flagged. The decision now uses the unrounded
  p-values, and the tables still show four decimals.

- `ripd()` computed `RIPD_RS` from arbitrary values when the covariance
  matrix of `RIPD_R` and `RIPD_S` was singular, which happens when all
  examinees who responded to an item have the same probability of a
  correct response (a tiny constant was added to every entry of the
  matrix). The statistic is now computed with a generalized inverse and
  the rank of the matrix as degrees of freedom, and a warning names the
  items. Results for regular data are unchanged.

- `rdif()` and `crdif()` stopped with an obscure error when an examinee
  had no ability estimate: a missing value in `score`, an examinee
  without any response, an examinee left out by `min.resp` (so
  `min.resp` could not be used when `score = NULL`), or an examinee
  whose responses were all removed during purification. Such examinees
  are now excluded, with a warning for the first and the third case, and
  the functions stop when no examinee has an ability estimate.

- `rdif()` and `crdif()` stopped with an obscure error when purification
  left a single item, and `rdif()` also stopped with "non-conformable
  arguments" when purification removed the item with the most score
  categories from a test with polytomous items. They now work, and stop
  with a clear message when every item is flagged.

- `rdif()` and `crdif()` silently accepted responses outside the score
  categories (for example, an unrecorded missing code such as -9, which
  flagged every item), a `group` of the wrong length, a `focal.name`
  that is not in `group`, an `item.skip` that was not a set of item
  positions (a logical `item.skip` skipped only the first item during
  purification), an `alpha` outside (0, 1), a `max.iter` that was not a
  positive whole number, and a `group` with more than two values (the
  values other than `focal.name` were pooled into the reference group).
  They now stop with a clear message; use `grdif()` to compare more than
  two groups. Character responses are read as numbers.

## Minor Improvements

- `ctt()` records the `correct` setting in its `crit` element, and the
  `summary()` report names the item-total correlation used for flagging.

- `ctt()` and `ctt_distr()` stop with a clear message when `crit.p`,
  `crit.dis`, `crit.distractor`, or `total` is not valid, when fewer
  than two examinees remain, or when `ctt_distr()` receives numeric and
  non-numeric keys together.

- Inputs that were silently misread now stop with a clear message.
  Functions that take item metadata reject parameter columns stored as
  text, a missing or too small `cats`, a dichotomous model with `cats`
  other than 2, and a polytomous item with more thresholds than
  `cats - 1`, including any value in `par.3` of a two-category GRM or
  GPCM item. `shape_df()` and `simdat()` reject parameter vectors that
  do not match the items, `shape_df_fipc()` checks `fix.loc`,
  `gen.weight()` checks `dist` and `theta`, and `prm()` rejects an
  unknown `pr.model`.

- `prm()` uses the GRM when `pr.model` is omitted, `simdat()` repeats a
  single `pr.model` for all polytomous items, `shape_df()` accepts
  `par.drm` without `g`, and item metadata without a `par.3` column is
  accepted when every item is 1PLM or 2PLM.

- `plot.info()` and `plot.traceline()` show item IDs as given
  (`"item-2"` was shown as `item.2`), and the panel titles of
  `plot.traceline()` read `Score: 0` instead of `Score: resp.0`.

- The warning for `"DRM"` items reads "All 'DRM' items are treated as
  '3PLM' items." and no longer mentions item parameter estimation,
  because it also appears in functions such as `info()` and `simdat()`.

- `est_score()` stops with a clear message when `data` has a number of
  columns other than the number of items (it scored the wrong items
  before), when a response is not an integer within the categories of
  its item, and when `method`, `stval.opt`, `tol`, or `max.iter` is not
  valid (a `max.iter` below 1 or not a whole number removed the
  iteration limit). For an `est_irt` object, supplying `data` or `D` is
  now an error instead of being ignored.

- `est_score()` now has a method for `est_item` objects, as its
  documentation already stated.

- `irtfit()` and `sx2_fit()` stop with a clear message when `data` has a
  number of columns other than the number of items or a response is not
  an integer within the categories of its item (such a response was
  counted in the wrong cells before), and they read character and factor
  responses as numbers instead of level codes. `irtfit()` also stops for
  a `score` of the wrong length, scores without variation, an invalid
  `range.score`, and a `loc.theta` other than `"average"` or `"middle"`
  (ignoring case; any other value was treated as `"average"`), leaves
  out examinees with a missing ability estimate instead of failing, and
  works when only one item remains after items with fewer than two
  responses are excluded. `plot.irtfit()` stops for a `type` other than
  `"both"`, `"icc"`, or `"sr"` and for an `item.loc` that is not an
  evaluated item. `sx2_fit()` stops for a test with a single item, which
  has no S-X2 statistic. The warnings for excluded items and for missing
  responses in `sx2_fit()` name the items, and the error for a response
  outside the categories names the item and its column.

- `bisection()` warns when the function values at `lb` and `ub` have the
  same sign, instead of returning a bound without notice, and it stops
  after exactly `max.it` iterations instead of `max.it + 1`. Its help
  page names the last returned element `delta`, as in the code.

- `est_score()` with `method = "INV.TCC"` stops with a clear message
  when the test characteristic curve does not reach a sum score within
  theta from -20 to 20, for example for items with very small slopes or
  guessing parameters that sum to nearly the maximum score.

- `ctt()` computes alpha with each item removed faster for large tests.

- `cac_lee()` and `cac_rud()` leave out examinees with a missing ability
  estimate (or standard error) with a warning. Before, `cac_lee()` and
  `cac_rud()` stopped with an obscure error for a missing ability
  estimate, and `cac_rud()` returned wrong indices for a missing
  standard error. Both stop for cut scores that are not in strictly
  ascending order, and `cac_rud()` stops for a standard error that is
  not positive. `cac_lee()` no longer stops when a performance level
  contains no integer summed score, and it is much faster with ability
  estimates.

- `panel_info()` stops with a clear message for a route map that does
  not define a panel (a module without transitions, a cycle, or a module
  outside the stages) instead of running forever, and it finds every
  stage when the modules are not numbered in the order of the stages.

- `reval_mst()` accepts item metadata with a missing guessing parameter,
  such as 2PLM items from `shape_df()`. It checks and completes the item
  metadata as `run_mst()` does, so `item.by.mod` and `item.by.path` hold
  plain data frames with upper-case model names and without parameter
  columns that are all `NA`. `reval_mst()` and `run_mst()` stop with a
  clear message when `cut_score` does not match the panel or `module`
  has a number of columns other than the number of modules, they read a
  `NULL` element of `cut_score` for a stage with a single module as an
  empty vector, and they assign the cut scores to the modules of a stage
  in the order of the module indices. `reval_mst()` also stops with a
  clear message when a sum score has no inverse TCC estimate (with
  `intpol = FALSE` or a `range.tcc` that does not cover the estimates);
  before, it stopped with an obscure error or returned `NA` for every
  row of `eval.tb`.

- `run_mst()` stops with a clear message when `tol` or `max.iter` in
  `route_score` or `final_score` is not valid (a `max.iter` of 0 or 2.5
  never ended), and it warns once about examinees without any observed
  response to the administered items.

- `find_cut()` no longer warns about two proper crossings when a
  crossing lies exactly on a grid point.

- `ripd()` results have a `print()` method that prints the flagged items
  and the RIPD statistics of the analyses without and with purification;
  `what` selects `"no_purify"` or `"with_purify"`. Before, `print()`
  listed the whole object, including the ability estimates of every
  examinee.

- The `est_irt` and `est_item` methods of `rdif()` and `crdif()` stop
  with a clear message when `data` or `D` is passed, because both are
  taken from the object. Before, they were silently ignored or caused an
  obscure error.

# irtQ 1.3.1

## Bug Fixes

- Fixed the routing in `run_mst()`. The ability estimate used to choose
  the next module was computed from the responses to the current module
  only. It is now computed from the responses to all modules
  administered so far, as intended. This applies to every routing method
  (`"bmat"`, `"mfi"`, and cut scores) and every `route_score` method, so
  simulated paths and results can differ from earlier versions, in
  particular from stage 3 onward.

- Fixed cut-score routing in `run_mst()` for route maps in which a
  module reaches only some modules of the next stage. The examinee could
  be sent to a module that did not match the cut scores. Only the cut
  scores that separate the reachable modules are now used, which is the
  rule applied by `reval_mst()`. Results are unchanged when every module
  of a stage can be reached.

- Fixed the final estimate of `run_mst()` for responses with missing
  values (`response` argument). The ML, WL, MLF, MAP, and EAP estimates
  used the parameters of other items instead of those of the observed
  items.

- `score_resp()` and `ctt_distr()` now match a data frame `key` to the
  item columns by item number when its `item` column is character or
  factor. Before, keys for tests with 10 or more items could be assigned
  to the wrong items.

- `shape_df()` with `default.par = TRUE` now repeats a single value of
  `cats` or `model` for all items. Before, a single `cats` value gave
  item IDs of "V1" for every item and a guessing parameter of 0 instead
  of 0.2 for 3PLM items.

- `irtfit()` now computes the observed category proportions directly
  from the frequencies. Before, they came from
  `janitor::adorn_percentages()`, which can include the total column in
  the denominator in some environments and halve the proportions, the
  residuals, and `overSR.prop`. Results are unchanged when the
  proportions were computed correctly.

- `catsib()` no longer overwrites `score` with `se` when `se` is given
  as a matrix or data frame.

- `cac_rud()` no longer stops with a dimnames error when a performance
  level has no examinees, and the label of the total row of `marginal`
  is now "marginal", as in `cac_lee()`.

- `simdat()` treats NA values in `g.drm` as zeros, as the item metadata
  input does. Before, the responses to those items were all NA.

- The Wald confidence intervals in `plot.irtfit()` now use the two-sided
  critical value, `qnorm(1 - alpha / 2)`. Before, they used
  `qnorm(1 - alpha)` and were 90% intervals at the default
  `alpha = 0.05`, while the Wilson intervals were 95% intervals.

## Minor Improvements

- `reval_mst()` now stops with an informative message when the modules
  in a stage differ in maximum sum score.

- In `run_mst()` and `reval_mst()`, an ability estimate equal to a cut
  score is now routed to the higher module, as in the classification
  rule of `cac_lee()` and `cac_rud()`. Results change only when an
  estimate equals a cut score exactly.

- `ctt()`, `ctt_distr()` (scored-category mode), and the CTT helper
  functions now stop with an informative error when an item score is not
  a whole number between 0 and `cats - 1`. When `cats` is inferred,
  every item has at least two categories, so an item that every examinee
  scores 0 on gets a difficulty of 0 and is flagged.

- `ctt()` reports the listwise deletion of incomplete rows once, and the
  item table has default row names.

- `plot.find_cut()` no longer passes an unused `inherit.aes` argument to
  `geom_vline()`, which caused warnings with some ggplot2 versions.

- `catsib()` now stops with an informative message when `score` is
  supplied without `se`. Before, it failed with a "missing value where
  TRUE/FALSE needed" error.

## Documentation

- Stated in `?reval_mst` that the recursion is based on inverse TCC
  estimates, that only inverse TCC scoring with cut-score routing is
  supported, and that modules in a stage must have the same maximum sum
  score; corrected the documented default of `theta` and the description
  of the returned list.

- In `?run_mst`, described the cumulative routing estimate and the
  relation to `reval_mst()`, and stated that `"EAP.SUM"` and `"INV.TCC"`
  count a missing response as 0 in the final sum score.

- Corrected the interpretation of alpha with the item removed in `?ctt`.
  The `ctt()` examples and the CTT article now use simulated IRT data.

- Corrected the description of the expected frequencies in `pcd2()` and
  stated that `crit.val = NULL` flags no items.

- Reworded the motivation of `ripd()`.

- In `?simMST`, noted that the item parameters are on the D = 1.702
  scale. In `?find_cut`, clarified that the cut scores in `simMST` were
  obtained with `find_cut()` and softened the statement about path
  reversals.

- Corrected the probability matrix of the first example in `?lwrc`, the
  GPCM formula note in `?irtQ`, and the class of `prob.cats` in
  `?traceline`. Corrected statements in the README, the vignette
  overview, and the articles (shrinkage of MAP and EAP, `range.score` in
  `irtfit()`, `fix.id` in `est_mg()`, the effect of `EmpHist` in FIPC,
  the fixed-slope 1PLM, fixed guessing, the CATSIB regression
  correction, the purification procedure, and the usage notes of the
  utility functions). The DIF article now simulates item difficulties in
  a narrower range so that the pooled calibrations converge.

# irtQ 1.3.0

## New Features

- New classical test theory (CTT) functions: `ctt()` for item- and
  test-level statistics, with `print()` and `summary()` methods;
  `freq_score()` for a total-score frequency table; `ctt_distr()` for
  option or score-category distributions with point-biserial
  correlations; and `score_resp()` for scoring raw selected-response
  data with an answer key. Item-total and option-total correlations are
  reported both with and without the item.

- Added a vignette, "Introduction to irtQ" (`vignette("irtQ")`), that
  walks through calibration, scoring, and model-fit evaluation.

## Bug Fixes

- Fixed the EM convergence check in `est_irt()` and `est_mg()`. It used
  the largest signed change instead of the largest absolute change, so
  the EM could stop before the `Etol` criterion was met, most often when
  some parameters were held constant (for example, the guessing
  parameter of 1PLM and 2PLM items, or slopes fixed by `fix.a.1pl` or
  `fix.a.gpcm`). Refitting a model can change the estimates and take
  more EM cycles. A fit from an earlier version that reported
  `maxpar.diff` equal to 0 most likely stopped early.

- `est_irt()` and `est_mg()` no longer warn that the convergence
  criteria are not satisfied when `fipc.method = "OEM"` or when the EM
  converges at the `MaxE`-th cycle. In these cases, `test.2` can now
  report a possible local maximum.

- `est_item()` now uses Beta(5, 16) as the default prior for the
  guessing parameter, as in `est_irt()` and `est_mg()`.

- `getirt()` now returns `posterior.dist` and `scale.D` for objects from
  `est_irt()` and `est_mg()`, as documented.

## Minor Improvements

- `find_cut()` now uses `D = 1` by default, like the other functions in
  the package (previously 1.702). Supply `D = 1.702` to reproduce
  earlier results.

- `shape_df_fipc()` now accepts a single value of `cats` or `model` for
  all new items.

- The mirt package moved from Imports to Suggests. Only `bring.mirt()`
  uses it, and it now stops with an informative message when mirt is not
  installed.

## Documentation

- Added examples to `print.ctt()`, `summary.ctt()`,
  `print.summary.ctt()`, and the data sets.

- Corrected the descriptions of `Etol`, `test.1`, `test.2`, and
  `maxpar.diff`, the return class of `est_mg()`, `dif_item` in the
  purified results of `rdif()`, `crdif()`, and `grdif()`, and the
  component names of `simCAT_DC` and `simCAT_MX`.

- Corrected bibliographic details, added DOIs, cited the polytomous RDIF
  extensions in `rdif()` and `crdif()`, and listed every reference cited
  in the help pages on the package help page (`?irtQ`).

# irtQ 1.2.0

## New Features

- Added a new function, `find_cut()`, which identifies TIF-crossing
  routing cut scores for MST panels. For each adjacent module pair
  within a stage, the function locates the theta where the two modules'
  test information functions (TIFs) intersect. A warning is issued when
  the mean difficulty order of modules within a stage differs from their
  input index order.

- Added a new S3 method, `plot.find_cut()`, which visualizes TIF curves
  and routing cut scores stage by stage using **ggplot2** facets.
  Proper, anomalous, and unselected cut scores are distinguished by line
  type. A `layout` argument (`"vertical"` / `"horizontal"`) controls the
  facet orientation.

- Added a new function, `run_mst()`, which simulates MST administrations
  for a given panel structure and returns response data along with
  ability and routing information for each simulated examinee.

- Added a new exported dataset, `simIPD`, which contains simulated CAT
  response data for illustrating IPD detection with `ripd()` and
  `pcd2()`. The dataset represents one replication of a CAT simulation
  (N = 3,000; test length = 30; 360-item 3PLM pool) in which 5% of items
  (18 items) had both $a$ and $b$ parameters decreased by 0.5.

## Minor Improvements

- `est_irt()`, `est_item()`, and `est_mg()` now accept a partial
  `control` list. Users can specify only the arguments they wish to
  override (e.g., `control = list(iter.max = 500)`); unspecified
  arguments fall back to their defaults via `modifyList()`.

- Reorganized the **pkgdown** reference page: `ripd()` and `pcd2()` are
  now grouped under a new *Item Parameter Drift (IPD)* section, and
  `reval_mst()`, `panel_info()`, `find_cut()`, `plot.find_cut()`, and
  `run_mst()` are grouped under a new *Multistage-Adaptive Test (MST)*
  section.

## Documentation

- Expanded the documentation for `ripd()` by adding a `@details` section
  that covers the theoretical background of the RIPD framework,
  asymptotic distributions of the three RIPD statistics, a drift-type
  diagnostic guide, the CAT-specific three-step workflow, and the
  purification procedure. Also added a `\donttest{}` example
  demonstrating a complete CAT-based IPD detection workflow using the
  `simIPD` dataset. Updated `@references` to include Lim & Choe (2023)
  and replaced the previous conference paper citation with the journal
  reference (Lim & Han, 2026).

- Added a `\donttest{}` example to `pcd2()` demonstrating CAT-based IPD
  detection using the `simIPD` dataset, including the bootstrap critical
  value procedure described in Lim & Han (2026).

- Updated the *MST Panel Evaluation and Simulation* article
  (`vignettes/articles/mst-panel-evaluation.Rmd`) to introduce
  `run_mst()` and extend the existing `reval_mst()` content with routing
  and scoring examples, a `find_cut()`-based principled cut score
  derivation, a side-by-side routing method comparison, and a Monte
  Carlo-vs-analytical validation example (Example 6).

## Bug Fixes

- Rebuilt the `simMST` dataset: the previous version had 9 items
  duplicated across non-adjacent modules because the original assembly
  only enforced no-overlap within a single routing pathway. The new
  version enforces a global no-overlap constraint across all 7 modules
  (56 unique items total) and adds a mean(*b*) per-module band
  constraint; cut scores were regenerated via `find_cut()` on the
  rebuilt modules.

# irtQ 1.1.0

## Major Improvements

- Improved the speed and reduced memory usage of item parameter
  estimation and standard error computation in `est_irt()`,
  `est_item()`, and `est_mg()`.
- Improved the computational speed of `est_score()` by up to 52% for
  dichotomous items (N = 10,000) and up to 27% for mixed-format tests,
  through a series of optimizations.
- Improved the computational speed of `sx2_fit()` substantially by
  replacing the O(J^2) Lord-Wingersky recursion with a forward-backward
  pass (up to 11x faster for mixed-format tests with J = 55 items) and
  vectorizing internal helper functions `expFreq()`, `obsFreq()`, and
  the PRM category-collapsing routine.

## New Features

- Added a unit test suite using the **testthat** 3rd edition
  (`testthat >= 3.0.0`). Tests cover core functions including `drm()`,
  `prm()`, `est_irt()`, `est_score()`, `est_mg()`, `rdif()`, `crdif()`,
  and `catsib()`, with the relevant tests across dichotomous,
  polytomous, and mixed-format item scenarios (355 tests total).
- Added a new function, `ripd()`, which implements the Residual-based
  Item Parameter Drift (RIPD) detection framework. The function computes
  three RIPD statistics: $RIPD_R$, $RIPD_S$, and $RIPD_{RS}$ (one for
  each item). $RIPD_R$ captures uniform item parameter drift (IPD) via
  differences in mean raw residuals between groups, $RIPD_S$ captures
  nonuniform IPD via differences in mean squared residuals, and
  $RIPD_{RS}$ is a combined chi-square-based statistic sensitive to both
  types of drift. An optional purification procedure is also supported.

## New Articles

- Launched the irtQ documentation website at
  <https://hwangQ.github.io/irtQ/>, built with **pkgdown**. The site
  includes a full function reference index and the following vignettes
  covering the complete irtQ workflow:
  - *Getting Started with irtQ*: an end-to-end overview of the package
    workflow.
  - *Item Parameter Estimation*: detailed guidance on `est_irt()`,
    `est_item()`, and `est_mg()`.
  - *Ability Estimation*: scoring methods available in `est_score()`.
  - *Model-Data Fit Evaluation*: using `irtfit()` and `sx2_fit()` to
    assess model fit.
  - *DIF Detection*: applying `rdif()`, `grdif()`, and `catsib()` to
    detect item bias.
  - *Classification Accuracy and Consistency*: computing indices via
    `cac_lee()` and `cac_rud()`.
  - *Utility Functions*: usage of `info()`, `traceline()`, `lwrc()`,
    `simdat()`, and related helpers.
  - *Evaluating MST Panels with `reval_mst()`*: measurement precision
    and bias evaluation for multistage adaptive tests.

## Bug Fixes

- Fixed a minor bug in `sx2_fit()` that caused incorrect cell collapsing
  between two adjacent score categories for polytomous items when
  computing the S-$X^2$ item fit statistic.
- Fixed minor bugs in `est_score()` and `info()`.
- Resolved an issue in `catsib()` where the function failed when all
  responses were missing (NA) in either the reference or focal group.
- Updated `cac_rud()` to include the `x` argument, allowing users to
  pass item metadata data frames directly.
- Revised default `control` parameters in `est_irt()`, `est_item()`, and
  `est_mg()`, and updated the documentation accordingly.
- Fixed a minor bug in `est_score()` function in terms of Newton-Raphson
  method.
- Fixed multiple stability issues in `catsib()`:
  - The final bin exclusion step in `catsib_item()` used a hardcoded
    threshold of 3 instead of the user-supplied `min.binsize` argument,
    causing inconsistent bin filtering behavior.
  - The reliability estimate `rho2` in `catsib_one()` was not clamped to
    $[0, 1]$, so when `errvar > sigma2` (e.g., very few items or
    purification cascade), a negative `rho2` reversed the regression
    correction direction, inflating the Type I error rate.
  - When `errvar >= sigma2` during purification, `rho2` collapsed to 0,
    causing all corrected scores to converge to the group mean. With
    group mean differences (impact), this produced empty bin data frames
    and an invalid purification result. A minimum floor of 0.05 is now
    enforced for `rho2` to preserve score spread.
- Fixed a critical bug in `covirt()` where the guessing parameter
  (`par[,3]`) was incorrectly passed as the difficulty parameter
  (`par[,2]`) to the `integrand()` function for DRM items. This caused
  the gradient computation to receive `c = b`, which zeroed out the
  $\partial P/\partial a$ and $\partial P/\partial b$ gradient
  components and produced a singular Fisher information matrix.
  Additionally, `NA` values in `par.3` for 1PLM and 2PLM items are now
  substituted with 0 prior to gradient evaluation to prevent `NA`
  propagation.

# irtQ 1.0.0

- The documentation for the `irtQ` package has been revised to reflect
  updates to function behavior, fix typos, and provide more relevant and
  detailed information for existing functions.

- A new function, `crdif()`, has been added. This function computes
  three statistics from the residual-based DIF detection framework using
  categorical residuals (RDIF-CR). It allows for the detection of global
  DIF, particularly in polytomously scored items.

- A new function, `shape_df_fipc()`, has been introduced. This function
  merges fixed-item metadata with automatically generated metadata for
  new items and produces a single data frame ordered by test position.
  It is designed to support fixed item parameter calibration (FIPC) via
  the `est_irt()` function.

- The `plot()` method has been enhanced to support the display of all
  item characteristic curves for a given item in a single panel.

- The `pcd2()` function has been updated to include a purification
  procedure.

- The `rdif()` and `catsib()` functions now include an `item.skip`
  argument. This allows users to specify a numeric vector of item
  indices to exclude from the DIF analysis.

# irtQ 0.2.1

- Enhanced functionality of the `bind.fill()` function by adding a new
  argument `fill`. The value in the argument is used to fill in missing
  data when aligning datasets.

- Fixed a bug within the `est_irt()` function that was previously unable
  to implement the fixed item parameter calibration (FIPC) when only
  freely estimating a single item given that all other items are fixed.

- Added a new function, `reval_mst()`, which evaluates the measurement
  precision and bias in Multistage-Adaptive Test (MST) panels using a
  recursion-based evaluation method introduced by Lim et al. (2020).

- Added a new function, `pcd2()`, which computes the Pseudo-count
  $D^{2}$ statistics (Cappaert et al., 2018; Stone, 2000) to detect item
  parameter drift.

# irtQ 0.2.0

- Introduced Warm's (1989) Weighted Likelihood (WL) estimation method to
  the `est_score()` function. This WL scoring method can now be utilized
  by setting `method = "WL"`.

- Enhanced the speed of ability parameter estimation in the
  `est_score()` function when using the ML, MLF, or MAP methods for the
  `method` argument. The updated version performs approximately three
  times faster than its predecessor.

- Addressed a bug within the `est_score()` function that was previously
  unable to accurately compute scores when only a single item data was
  provided. This issue was occurring with the EAP.SUM and INV.TCC
  estimation methods.

- Added two new functions for computing classification accuracy and
  consistency: `cac_rud()` and `cac_lee()`.

  - `cac_rud`: This function implements Rudner's (2001, 2005) method for
    computing classification accuracy and consistency. It takes cut
    scores, ability estimates, standard errors, and optional weights as
    inputs and returns a list containing a confusion matrix, marginal
    and conditional classification accuracy and consistency indices, the
    probability of being assigned to each level category, and the cut
    scores used in the analysis.
  - `cac_lee`: This function implements Lee's (2010) method for
    computing classification accuracy and consistency. It takes a data
    frame containing item metadata, cut scores, optional ability
    estimates, optional weights, a scaling factor, and a logical value
    indicating the cut score metric as inputs. It returns a list similar
    to `cac_rud`.

- Added a new function, `llike_score()`, which computes the
  loglikelihood of ability parameters given the item parameters and
  response data.

- Enhanced functionality of the `rdif()` and `grdif()` functions: Both
  now support the graded response model (GRM) and generalized partial
  credit model (GPCM).

- Fixed an issue in the `grdif()` function that inaccurately calculated
  the GRDIF statistics when group membership was specified in a
  non-standard way. Specifically, the problem arose when 0 wasn't used
  as the reference group and consecutive numbers (e.g., 1, 2, 3) weren't
  used to represent focal groups in the `group` argument.

# irtQ 0.1.1

- Resolved the misalignment issue of standard errors in the output of
  the `est_irt()` function when `fix.a.1pl = TRUE` is specified and the
  items are calibrated using the 1PLM.

- Added a new function, `grdif()`, to perform differential item
  functioning (DIF) analysis across multiple groups. This function
  calculates three generalized IRT residual DIF (GRDIF) statistics. For
  more information about the function and its usage, please refer to the
  accompanying documentation.

- Fixed several typos in the manual documentation

# irtQ 0.1.0

- Initial release on CRAN

- The `irtQ` package is a successor of the `irtplay` package which was
  retracted from R CRAN due to the intellectual property (IP) violation.
  All issues of the IP violation have been clearly resolved in the
  `irtQ` package.

- Most of the functions the `irtQ` package are identical in appearance
  and functionality to those of `irtplay` package except a few functions
  (e.g., `shape_df()`, `est_score()`). However, the computing speed of
  several functions (e.g., `est_irt()`, `est_score()`, `lwrc()`) in the
  `irtQ` package are faster than the previous ones in the `irtplay`
  package. Read the documentation carefully prior to using the
  functions.
