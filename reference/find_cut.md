# Find TIF-Crossing Cut Scores for MST Routing

Computes routing cut scores for a multistage test (MST) by identifying
the theta values at which adjacent modules' Test Information Functions
(TIFs) cross in the correct direction. The resulting cut scores can be
passed directly to
[`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md) as its
`cut_score` argument (with `route_method = NULL`), providing a
principled and psychometrically grounded alternative to heuristic
cut-score selection.

## Usage

``` r
find_cut(
  x,
  module,
  route_map,
  D = 1.702,
  theta_range = c(-6, 6),
  n_grid = 2001L,
  ref_theta = 0
)
```

## Arguments

- x:

  A data frame of item metadata in the standard irtQ format (columns:
  `id`, `cats`, `model`, `par.1`, `par.2`, ...). See
  [`shape_df`](https://hwangQ.github.io/irtQ/reference/shape_df.md) for
  details.

- module:

  A binary integer matrix of dimensions *J* \\\times\\ *M*, where *J* is
  the number of items and *M* is the number of modules. A value of 1 at
  row *j*, column *m* indicates that item *j* belongs to module *m*.
  This is the same `module` argument used in
  [`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md).

- route_map:

  A binary square matrix of dimensions *M* \\\times\\ *M* defining the
  MST transition structure. A value of 1 at row *i*, column *j*
  indicates that examinees can be routed from module *i* to module *j*.
  This is the same `route_map` argument used in
  [`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md).

- D:

  A numeric scaling constant for the IRT model. Default is `1.702`.

- theta_range:

  A numeric vector of length 2 specifying the theta range over which to
  search for TIF crossings. Default is `c(-6, 6)`.

- n_grid:

  An integer specifying the number of equally spaced theta points used
  for the initial sign-change scan. A finer grid reduces the chance of
  missing a crossing but increases computation time. Default is `2001`.

- ref_theta:

  A single numeric value used to break ties when multiple proper
  crossings are found for an adjacent module pair. The crossing closest
  to `ref_theta` is selected. Default is `0`.

## Value

An object of class `"find_cut"`, which is a named list with three
elements:

- `cut_score`:

  A list of length *S* - 1 (where *S* is the number of stages),
  compatible with the `cut_score` argument of
  [`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md). Each
  element is a numeric vector of cut scores for the transition into that
  stage, sorted in ascending order. Stages with only one module have an
  empty numeric vector (`numeric(0)`).

- `details`:

  A named list of per-stage diagnostic information, including module
  indices, mean item locations, difficulty-sorted module order, and
  per-pair crossing details (proper crossings, anomalous crossings, and
  the selected cut score).

- `tif_data`:

  A `tibble` with columns `theta`, `module`, `stage`, and `tif`,
  providing TIF curves for all modules at all stages (including stage
  1). Used by
  [`plot.find_cut`](https://hwangQ.github.io/irtQ/reference/plot.find_cut.md)
  to visualise TIF curves and crossing points.

## Details

### Motivation: the MFI routing problem

When [`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
uses Maximum Fisher Information (MFI) routing (`route_method = "mfi"`),
it selects the next-stage module that maximizes TIF at the examinee's
current ability estimate. In principle, a harder module should always
have a higher TIF for high-ability examinees, and a lower TIF for
low-ability examinees, making MFI routing work as intended. However, in
practice, a harder module's TIF curve can peak at a theta value that is
unexpectedly low, so that the harder module *also* has the highest TIF
near the low end of the ability scale. In this case, a low-ability
examinee is incorrectly routed into the harder module – the opposite of
what the MST is designed to do, which is called an *anomalous routing*
or *path reversal*.

### Solution: TIF-crossing cut scores

`find_cut()` resolves this problem by computing cut scores from the TIF
curves themselves, rather than relying on pure MFI at the time of
routing. Specifically, for each pair of adjacent modules (ordered by
mean item difficulty), the function finds the theta value at which the
harder module's TIF first overtakes the easier module's TIF *from below*
as theta increases. This crossing point is the natural boundary between
the two modules: below it, the easier module is more informative; above
it, the harder module is more informative. Using this point as a fixed
cut score pre-empts the path reversal problem entirely.

### Proper vs. anomalous crossings

Two TIF curves may cross more than once. `find_cut()` distinguishes two
types of crossing:

- **Proper crossing**: the difference
  \\\text{TIF}\_{\text{harder}}(\theta) -
  \text{TIF}\_{\text{easier}}(\theta)\\ changes sign from negative to
  positive (positive slope at the crossing). This is the correct,
  expected crossing in the middle of the theta scale.

- **Anomalous crossing**: the difference changes sign from positive to
  negative (negative slope at the crossing). This corresponds to the
  pathological situation in which the harder module temporarily
  dominates at the low end of the theta scale – exactly the pattern that
  causes path reversals under MFI routing.

Only proper crossings are used as cut scores. Anomalous crossings are
reported in `$details` for diagnostic purposes but excluded from the
returned `cut_score`. If no proper crossing is found for a module pair,
an informative error is thrown. If multiple proper crossings are found
(unusual for well-designed modules), a warning is issued and the one
closest to `ref_theta` is selected.

### Algorithm

The function proceeds stage by stage (from stage 2 onward). At each
stage:

1.  For each adjacent pair of modules in the order supplied by the user
    (ascending module index), the difference
    \\\text{TIF}\_{\text{right}}(\theta) -
    \text{TIF}\_{\text{left}}(\theta)\\ is evaluated on a fine theta
    grid of `n_grid` points, where "left" and "right" refer to the
    lower- and higher-index module in the pair. Sign changes on the grid
    are detected and then refined with
    [`uniroot`](https://rdrr.io/r/stats/uniroot.html).

2.  The slope at each refined root is estimated by a finite difference
    (step size 0.01) to classify it as proper or anomalous.

3.  The selected proper crossing is stored as the cut score for that
    pair.

### Module index vs. difficulty order

[`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
internally orders candidate modules by ascending module index and maps
routing rank 1 to the lowest-index module. For cut-score routing to
assign low-theta examinees to the easiest module, the modules at each
stage must be numbered in ascending difficulty order (easiest module =
lowest index). `find_cut()` checks this assumption and issues a warning
if the difficulty order does not match the index order. Regardless of
this warning, cut scores are always computed between consecutive module
pairs in the supplied index order (e.g., modules 2 & 3 and 3 & 4, not 2
& 4 and 4 & 3). The returned cut scores within each stage are sorted in
ascending order, as required by
[`cut`](https://rdrr.io/r/base/cut.html).

## See also

[`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md),
[`panel_info`](https://hwangQ.github.io/irtQ/reference/panel_info.md),
[`reval_mst`](https://hwangQ.github.io/irtQ/reference/reval_mst.md),
[`info`](https://hwangQ.github.io/irtQ/reference/info.md)

## Examples

``` r
## -- Setup: use the built-in simMST 1-3-3 panel --------------------------
## simMST is a 7-module, 3-stage MST panel (8 dichotomous 3PLM items each).
## Modules 1 (routing), 2-4 (stage 2: easy/medium/hard),
##         5-7 (stage 3: easy/medium/hard).
x         <- simMST$item_bank
module    <- simMST$module
route_map <- simMST$route_map

## -- Find TIF-crossing cut scores -----------------------------------------
## For each adjacent module pair at stages 2 and 3, find_cut() identifies
## the theta at which the harder module's TIF first exceeds the easier
## module's TIF (proper crossing), and returns it as a cut score.
cut_result <- find_cut(x = x, module = module, route_map = route_map)

## Print a summary: crossing points found, anomalous crossings excluded,
## and the final selected cut scores per stage transition.
print(cut_result)
#> MST TIF-Crossing Cut Score Results
#> ==================================================
#> 
#> stage.2
#> -----------------------------------
#>   Modules (index order)    : 2, 3, 4
#>   Modules (difficulty order): 2, 3, 4
#>   Mean item locations      : -0.87, -0.0574, 0.7858
#>   Selected cut score(s)    : -0.4451, 0.4589
#>   Pair details:
#>     [mod2_vs_mod3]
#>       Proper crossing(s)    : 1  [theta = -0.4451]
#>       Anomalous crossing(s) : 1  [theta = 3.3893]  (excluded)
#>       Selected cut score    : -0.4451
#>     [mod3_vs_mod4]
#>       Proper crossing(s)    : 1  [theta = 0.4589]
#>       Selected cut score    : 0.4589
#> 
#> stage.3
#> -----------------------------------
#>   Modules (index order)    : 5, 6, 7
#>   Modules (difficulty order): 5, 6, 7
#>   Mean item locations      : -0.5586, -0.042, 0.5066
#>   Selected cut score(s)    : -0.4171, 0.4532
#>   Pair details:
#>     [mod5_vs_mod6]
#>       Proper crossing(s)    : 1  [theta = -0.4171]
#>       Anomalous crossing(s) : 1  [theta = 1.617]  (excluded)
#>       Selected cut score    : -0.4171
#>     [mod6_vs_mod7]
#>       Proper crossing(s)    : 1  [theta = 0.4532]
#>       Anomalous crossing(s) : 1  [theta = -1.9874]  (excluded)
#>       Selected cut score    : 0.4532
#> ==================================================
#> Usage: run_mst(..., route_method = NULL, cut_score = <result>$cut_score)

## -- Visualise TIF curves and cut scores ---------------------------------
## plot() shows TIF curves for every module faceted by stage.
## Selected cut scores appear as solid black vertical lines with theta
## labels at the crossing point. Anomalous crossings (if any) appear
## as dashed red lines.
plot(cut_result)                         # stages stacked vertically (default)

plot(cut_result, layout = "horizontal")  # stages side by side

plot(cut_result, show_anomalous = FALSE) # hide anomalous crossing markers


## Inspect the cut_score element -- a list directly compatible with run_mst()
## cut_score[[1]]: two cut scores for the stage-1 -> stage-2 transition
## cut_score[[2]]: cut scores for the stage-2 -> stage-3 transition
cut_result$cut_score
#> $stage.2
#> [1] -0.4450901  0.4588774
#> 
#> $stage.3
#> [1] -0.4170988  0.4531584
#> 

## Compare with the manually specified cut scores stored in simMST
simMST$cut_score
#> $stage.2
#> [1] -0.4450901  0.4588774
#> 
#> $stage.3
#> [1] -0.4170988  0.4531584
#> 

## -- Use the cut scores in run_mst() --------------------------------------
## Pass cut_result$cut_score directly to run_mst() with route_method = NULL.
## This replaces pure MFI routing with TIF-crossing-based fixed cut scores,
## preventing anomalous path reversals at the extremes of the theta scale.
if (FALSE) { # \dontrun{
set.seed(1)
theta_true <- rnorm(500)
result <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true,
  route_method = NULL,
  cut_score    = cut_result$cut_score
)
} # }
```
