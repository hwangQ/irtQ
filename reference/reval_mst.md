# Recursion-based MST evaluation method

This function evaluates the measurement precision and bias in
Multistage-Adaptive Test (MST) panels using a recursion-based evaluation
method introduced by Lim et al. (2021). This function computes
conditional biases and standard errors of measurement (CSEMs) across a
range of IRT ability levels, facilitating efficient and accurate MST
panel assessments without extensive simulations.

## Usage

``` r
reval_mst(
  x,
  D = 1,
  route_map,
  module,
  cut_score,
  theta = seq(-5, 5, 1),
  intpol = TRUE,
  range.tcc = c(-7, 7),
  tol = 1e-04
)
```

## Arguments

- x:

  A data frame containing the metadata for the item bank, which includes
  important information for each item such as the number of score
  categories and the IRT model applied. This metadata is essential for
  evaluating the MST panel, with items selected based on the
  specifications in the `module` argument. To construct this item
  metadata efficiently, the
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)
  function is recommended. Further details on utilizing item bank
  metadata along with `module` for MST panel evaluation are provided
  below.

- D:

  A scaling constant used in IRT models to make the logistic function
  closely approximate the normal ogive function. A value of 1.7 is
  commonly used for this purpose. Default is 1.

- route_map:

  A binary square matrix that defines the MST structure, illustrating
  transitions between modules and stages. This concept and structure are
  inspired by the `transMatrix` argument in the `randomMST()` function
  from the mstR package (Magis et al., 2017), which provides a framework
  for representing MST pathways. For comprehensive understanding and
  examples of constructing `route_map`, refer to the mstR package (Magis
  et al., 2017) documentation. Also see below for details.

- module:

  A binary matrix that maps items from the item bank specified in `x` to
  modules within the MST framework. This parameter's structure is
  analogous to the `modules` argument in the `randomMST()` function of
  the mstR package, enabling precise item-to-module assignments for MST
  configurations. For detailed instructions on creating and utilizing
  the `module` matrix effectively, consult the documentation of the mstR
  package (Magis et al., 2017). Also see below for details.

- cut_score:

  A list defining cut scores for routing test takers through MST stages.
  Each list element is a vector of cut scores for advancing participants
  to subsequent stage modules. In a 1-3-3 MST configuration, for
  example, `cut_score` might be defined as
  `cut_score = list(c(-0.5, 0.5), c(-0.6, 0.6))`, where `c(-0.5, 0.5)`
  are thresholds for routing from the first to the second stage, and
  `c(-0.6, 0.6)` for routing from the second to the third stage. This
  setup facilitates tailored test progression based on performance.
  Further examples and explanations are available below.

- theta:

  A vector of ability levels (theta) at which the MST panel's
  performance is assessed. This allows for the evaluation of measurement
  precision and bias across a continuum of ability levels. The default
  range is `theta = seq(-5, 5, 0.1)`.

- intpol:

  A logical value to enable linear interpolation in the inverse test
  characteristic curve (TCC) scoring, facilitating ability estimate
  approximation for observed sum scores not directly obtainable from the
  TCC, such as those below the sum of item guessing parameters. Default
  is TRUE, applying interpolation to bridge gaps in the TCC. Refer to
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  for more details and consult Lim et al. (2021) for insights into the
  interpolation technique within inverse TCC scoring.

- range.tcc:

  A vector to define the range of ability estimates for inverse TCC
  scoring, expressed as the two numeric values for lower and upper
  bounds. Default is to c(-7, 7).

- tol:

  A numeric value of the convergent tolerance for the inverse TCC
  scoring. For the inverse TCC scoring, the bisection method is used for
  optimization. Default is 1e-4.

## Value

This function returns a list of seven internal objects. The four objects
are:

- panel.info:

  A list of several sub-objects containing detailed information about
  the MST panel configuration, including:

  config

  :   A nested list indicating the arrangement of modules across stages,
      showing which modules are included in each stage. For example, the
      first stage includes module 1, the second stage includes modules 2
      to 4, and so forth.

  pathway

  :   A matrix detailing all possible pathways through the MST panel.
      Each row represents a unique path a test taker might take, based
      on their performance and the cut scores defined.

  n.module

  :   A named vector indicating the number of modules available at each
      stage.

  n.stage

  :   A single numeric value representing the total number of stages in
      the MST panel.

- item.by.mod:

  A list where each entry represents a module in the MST panel,
  detailing the item metadata within that module. Each module's metadata
  includes item IDs, the number of categories, the IRT model used
  (model), and the item parameters (e.g., par.1, par.2, par.3).

- item.by.path:

  A list containing item metadata arranged according to the paths
  through the MST structure. This detailed breakdown allows for an
  analysis of item characteristics along specific MST paths. Each list
  entry corresponds to a testing stage and path, providing item
  metadata. This structure facilitates the examination of how items
  function within the context of each unique path through the MST.

- eq.theta:

  Estimated ability levels (\\\theta\\) corresponding to the observed
  scores, derived from the inverse TCC scoring method. This provides the
  estimated \\\theta\\ values for each potential pathway through the MST
  stages. For each stage, \\\theta\\ values are calculated for each
  path, indicating the range of ability levels across the test takers.
  For instance, in a three-stage MST, the `eq.theta` list may contain
  \\\theta\\ estimates for multiple paths within each stage, reflecting
  the progression of ability estimates as participants move through the
  test. The example below illustrates the structure of `eq.theta` output
  for a 1-3-3 MST panel with varying paths:

  stage.1

  :   `path.1` shows \\\theta\\ estimates ranging from -7 to +7,
      demonstrating the initial spread of abilities.

  stage.2

  :   Multiple paths (`path.1`, `path.2`, ...) each with their own
      \\\theta\\ estimates, indicating divergence in ability levels
      based on test performance.

  stage.3

  :   Further refinement of \\\theta\\ estimates across paths,
      illustrating the final estimation of abilities after the last
      stage.

- cdist.by.mod:

  A list where each entry contains the conditional distributions of the
  observed scores for each module given the ability levels.

- jdist.by.path:

  Joint distributions of observed scores for different paths at each
  stage in a MST panel. The example below outlines the organization of
  `jdist.by.path` data in a hypothetical 1-3-3 MST panel:

  stage.1

  :   Represents the distribution at the initial stage, indicating the
      broad spread of test-taker abilities at the outset.

  stage.2

  :   Represents the conditional joint distributions of the observed
      scores as test-takers move through different paths at the stage 2,
      based on their performance in earlier stages.

  stage.3

  :   Represents a further refinement of joint distribution of observed
      scores as test-takers move through different paths at the final
      stage 3, based on their performance in earlier stages.

- eval.tb:

  A table summarizing the measurement precision of the MST panel. It
  contains the true ability levels (`theta`) with the mean ability
  estimates (`mu`), variance (`sigma2`), bias, and conditional standard
  error of measurement (CSEM) given the true ability levels. This table
  highlights the MST panel's accuracy and precision across different
  ability levels, providing insights into its effectiveness in
  estimating test-taker abilities.

## Details

The `reval_mst()` function evaluates an MST panel by implementing a
recursion-based method to assess measurement precision across IRT
ability levels. This approach, detailed in Lim et al. (2021), enables
the computation of conditional biases and CSEMs efficiently, bypassing
the need for extensive simulations traditionally required for MST
evaluation.

The `module` argument, used in conjunction with the item bank metadata
`x`, systematically organizes items into modules for MST panel
evaluation. Each row of `x` corresponds to an item, detailing its
characteristics like score categories and IRT model. The `module`
matrix, structured with the same number of rows as `x` and columns
representing modules, indicates item assignments with 1s. This precise
mapping enables the `reval_mst()` function to evaluate the MST panel's
performance by analyzing how items within each module contribute to
measurement precision and bias, reflecting the tailored progression
logic inherent in MST designs.

The `route_map` argument is essential for defining the MST's structure
by indicating possible transitions between modules. Similar to the
`transMatrix()` in the mstR package (Magis et al., 2017), `route_map` is
a binary matrix that outlines which module transitions are possible
within an MST design. Each "1" in the matrix represents a feasible
transition from one module (row) to another (column), effectively
mapping the flow of test takers through the MST based on their
performance. For instance, a "1" at the intersection of row *i* and
column *j* indicates the possibility for test takers to progress from
the module corresponding to row *i* directly to the module denoted by
column *j*. This structure allows `reval_mst()` to simulate and evaluate
the dynamic routing of test takers through various stages and modules of
the MST panel.

To further detail the `cut_score` argument with an illustration: In a
1-3-3 MST configuration, the list
`cut_score = list(c(-0.5, 0.5), c(-0.6, 0.6))` operates as a decision
guide at each stage. Initially, all test takers start in the first
module. Upon completion, their scores determine their next stage module:
scores below -0.5 route to the first module of the next stage, between
-0.5 and 0.5 to the second, and above 0.5 to the third. This pattern
allows for dynamic adaptation, tailoring the test path to individual
performance levels.

## References

Magis, D., Yan, D., & Von Davier, A. A. (2017). *Computerized adaptive
and multistage testing with R: Using packages catR and mstR*. Springer.

Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
approach to evaluate the performance of MST. *Journal of Educational
Measurement, 58*(2), 154-178.

## See also

[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# \donttest{
## ------------------------------------------------------------------------------
# Evaluation of a 1-3-3 MST panel using simMST data.
# This panel was assembled using module-level target TIFs and design
# constraints similar to those used in Lim et al.'s (2021) simulation study
# (it is not the identical dataset from that study).
# Details:
#    (a) Panel configuration: 1-3-3 MST panel
#    (b) Test length: 24 items (each module contains 8 items across all stages)
#    (c) IRT model: 3-parameter logistic model (3PLM)
## ------------------------------------------------------------------------------
# Load the necessary library
library(dplyr)
library(tidyr)
library(ggplot2)

# Import item bank metadata
x <- simMST$item_bank

# Import module information
module <- simMST$module

# Import routing map
route_map <- simMST$route_map

# Import cut scores for routing to subsequent modules
cut_score <- simMST$cut_score

# Import ability levels (theta) for evaluating measurement precision
theta <- simMST$theta

# Evaluate MST panel using the reval_mst() function
eval <-
  reval_mst(x,
    D = 1.702, route_map = route_map, module = module,
    cut_score = cut_score, theta = theta, range.tcc = c(-7, 7)
  )

# Review evaluation results
# The evaluation result table below includes conditional biases and
# standard errors of measurement (CSEMs) across ability levels
print(eval$eval.tb)
#>    theta           mu     sigma2         bias      csem
#> 1   -4.0 -4.808121816 3.42224454 -0.808121816 1.8499310
#> 2   -3.9 -4.789018117 3.43809047 -0.889018117 1.8542089
#> 3   -3.8 -4.766227061 3.45633941 -0.966227061 1.8591233
#> 4   -3.7 -4.739003828 3.47720522 -1.039003828 1.8647266
#> 5   -3.6 -4.706457625 3.50082146 -1.106457625 1.8710482
#> 6   -3.5 -4.667529927 3.52716739 -1.167529927 1.8780754
#> 7   -3.4 -4.620974273 3.55595499 -1.220974273 1.8857240
#> 8   -3.3 -4.565341408 3.58646261 -1.265341408 1.8937958
#> 9   -3.2 -4.498975738 3.61729989 -1.298975738 1.9019201
#> 10  -3.1 -4.420032009 3.64609197 -1.320032009 1.9094743
#> 11  -3.0 -4.326524655 3.66908363 -1.326524655 1.9154852
#> 12  -2.9 -4.216425518 3.68069437 -1.316425518 1.9185136
#> 13  -2.8 -4.087827617 3.67311048 -1.287827617 1.9165361
#> 14  -2.7 -3.939190341 3.63608699 -1.239190341 1.9068526
#> 15  -2.6 -3.769671645 3.55723963 -1.169671645 1.8860646
#> 16  -2.5 -3.579531039 3.42319361 -1.079531039 1.8501875
#> 17  -2.4 -3.370550870 3.22192852 -0.970550870 1.7949731
#> 18  -2.3 -3.146374631 2.94638073 -0.846374631 1.7165025
#> 19  -2.2 -2.912612919 2.59873030 -0.712612919 1.6120578
#> 20  -2.1 -2.676548271 2.19389663 -0.576548271 1.4811808
#> 21  -2.0 -2.446322106 1.76007533 -0.446322106 1.3266783
#> 22  -1.9 -2.229645390 1.33450055 -0.329645390 1.1552058
#> 23  -1.8 -2.032321625 0.95460517 -0.232321625 0.9770390
#> 24  -1.7 -1.857088754 0.64771766 -0.157088754 0.8048091
#> 25  -1.6 -1.703286390 0.42419405 -0.103286390 0.6513018
#> 26  -1.5 -1.567529547 0.27740271 -0.067529547 0.5266903
#> 27  -1.4 -1.445069807 0.18987969 -0.045069807 0.4357519
#> 28  -1.3 -1.331207524 0.14154025 -0.031207524 0.3762184
#> 29  -1.2 -1.222225352 0.11571234 -0.022225352 0.3401652
#> 30  -1.1 -1.115711117 0.10132976 -0.015711117 0.3183234
#> 31  -1.0 -1.010467372 0.09223884 -0.010467372 0.3037085
#> 32  -0.9 -0.906250276 0.08548769 -0.006250276 0.2923828
#> 33  -0.8 -0.803427274 0.07989210 -0.003427274 0.2826519
#> 34  -0.7 -0.702526078 0.07520746 -0.002526078 0.2742398
#> 35  -0.6 -0.603713902 0.07166265 -0.003713902 0.2676988
#> 36  -0.5 -0.506419973 0.06957502 -0.006419973 0.2637708
#> 37  -0.4 -0.409361950 0.06900157 -0.009361950 0.2626815
#> 38  -0.3 -0.311031582 0.06954823 -0.011031582 0.2637200
#> 39  -0.2 -0.210392034 0.07043473 -0.010392034 0.2653954
#> 40  -0.1 -0.107394541 0.07079722 -0.007394541 0.2660775
#> 41   0.0 -0.003020654 0.07009583 -0.003020654 0.2647562
#> 42   0.1  0.101195029 0.06841021  0.001195029 0.2615535
#> 43   0.2  0.203935849 0.06642380  0.003935849 0.2577282
#> 44   0.3  0.304760735 0.06507124  0.004760735 0.2550906
#> 45   0.4  0.404268069 0.06505829  0.004268069 0.2550653
#> 46   0.5  0.503674818 0.06655934  0.003674818 0.2579910
#> 47   0.6  0.604113025 0.06926578  0.004113025 0.2631839
#> 48   0.7  0.706116775 0.07271246  0.006116775 0.2696525
#> 49   0.8  0.809566254 0.07664014  0.009566254 0.2768396
#> 50   0.9  0.913989404 0.08119260  0.013989404 0.2849431
#> 51   1.0  1.018936817 0.08697665  0.018936817 0.2949180
#> 52   1.1  1.124237650 0.09523922  0.024237650 0.3086085
#> 53   1.2  1.230138950 0.10846069  0.030138950 0.3293337
#> 54   1.3  1.337437576 0.13152988  0.037437576 0.3626705
#> 55   1.4  1.447691648 0.17341863  0.047691648 0.4164356
#> 56   1.5  1.563507131 0.24891023  0.063507131 0.4989090
#> 57   1.6  1.688793130 0.37944575  0.088793130 0.6159917
#> 58   1.7  1.828800578 0.59172170  0.128800578 0.7692345
#> 59   1.8  1.989742189 0.91275843  0.189742189 0.9553839
#> 60   1.9  2.177883977 1.36127926  0.277883977 1.1667387
#> 61   2.0  2.398208325 1.93741933  0.398208325 1.3919121
#> 62   2.1  2.652995666 2.61499495  0.552995666 1.6170946
#> 63   2.2  2.940808726 3.34098096  0.740808726 1.8278350
#> 64   2.3  3.256275361 4.04428294  0.956275361 2.0110403
#> 65   2.4  3.590775050 4.65138229  1.190775050 2.1567064
#> 66   2.5  3.933796614 5.10289537  1.433796614 2.2589589
#> 67   2.6  4.274529547 5.36492913  1.674529547 2.3162317
#> 68   2.7  4.603255939 5.43215059  1.903255939 2.3306974
#> 69   2.8  4.912272678 5.32338528  2.112272678 2.3072463
#> 70   2.9  5.196279204 5.07302346  2.296279204 2.2523373
#> 71   3.0  5.452320839 4.72181852  2.452320839 2.1729746
#> 72   3.1  5.679446946 4.30951268  2.579446946 2.0759366
#> 73   3.2  5.878241351 3.87022633  2.678241351 1.9672891
#> 74   3.3  6.050342862 3.43043317  2.750342862 1.8521429
#> 75   3.4  6.198025435 3.00881210  2.798025435 1.7345928
#> 76   3.5  6.323867445 2.61718809  2.823867445 1.6177726
#> 77   3.6  6.430512839 2.26193109  2.830512839 1.5039718
#> 78   3.7  6.520512752 1.94540173  2.820512752 1.3947766
#> 79   3.8  6.596230710 1.66722536  2.796230710 1.2912108
#> 80   3.9  6.659794430 1.42531046  2.759794430 1.1938637
#> 81   4.0  6.713079533 1.21660661  2.713079533 1.1029989

# Generate plots for biases and CSEMs
p_eval <-
  eval$eval.tb %>%
  dplyr::select(theta, bias, csem) %>%
  tidyr::pivot_longer(
    cols = c(bias, csem),
    names_to = "criterion", values_to = "value"
  ) %>%
  ggplot2::ggplot(mapping = ggplot2::aes(x = theta, y = value)) +
  ggplot2::geom_point(mapping = ggplot2::aes(shape = criterion), size = 3) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      color = criterion,
      linetype = criterion
    ),
    linewidth = 1.5
  ) +
  ggplot2::labs(x = expression(theta), y = NULL) +
  ggplot2::theme_classic() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.key.width = unit(1.5, "cm"))
print(p_eval)

# }
```
