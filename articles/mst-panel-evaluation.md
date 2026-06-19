# MST Panel Evaluation and Simulation

## What is a Multistage-Adaptive Test (MST)?

A **Multistage-Adaptive Test (MST)** is a computer-based adaptive
testing design that sits between fully adaptive Computerized Adaptive
Testing (CAT) and traditional linear fixed-form testing. In an MST, the
test is divided into *stages*, each containing one or more pre-assembled
groups of items called *modules*. Routing rules — based on performance
on earlier stages — determine which module a test taker receives at each
subsequent stage.

### The Basic Structure

A typical MST panel is described by its *stage-module configuration*.
For example, a **1-3-3 panel** has:

- **Stage 1**: 1 routing module — everyone starts with the same items
- **Stage 2**: 3 modules of varying difficulty (e.g., easy, medium,
  hard)
- **Stage 3**: 3 modules of varying difficulty (e.g., easy, medium,
  hard)

The diagram below illustrates the flow of test takers through a 1-3-3
panel (E = Easy, M = Medium, H = Hard):

![1-3-3 MST panel design. From Stage 1, all examinees begin at module
1M. After Stage 1, they are routed to one of three Stage 2 modules (2E,
2M, 2H) based on their estimated ability. After Stage 2, they are routed
to one of three Stage 3 modules (3E, 3M, 3H).](img/mst_133_design.png)

1-3-3 MST panel design. From Stage 1, all examinees begin at module 1M.
After Stage 1, they are routed to one of three Stage 2 modules (2E, 2M,
2H) based on their estimated ability. After Stage 2, they are routed to
one of three Stage 3 modules (3E, 3M, 3H).

A test taker’s *pathway* through the MST is the sequence of modules they
take, such as 1M → 2E → 3E (for a low-ability examinee) or 1M → 2H → 3H
(for a high-ability examinee). In a 1-3-3 panel above, there are 7 such
pathways.

### Advantages of MST Over Linear Tests and CAT

| Feature | Linear | CAT | MST |
|----|----|----|----|
| Adaptivity | None | Item-by-item | Stage-by-stage |
| Content review | Full review allowed | (Generally) Not allowed | Allowed within module |
| Item exposure control | Easy | Difficult | Moderate |
| Test assembly | Pre-assembled | On-the-fly | (Generally) Pre-assembled modules |
| Operational cost | Low | High | Moderate |

MST is especially popular in high-stakes certification and licensure
exams because it allows examinees to review and change answers within a
module (like a paper test) while still adapting the difficulty level to
each examinee.

------------------------------------------------------------------------

## Why Simulate MST Panel Performance?

Designing an operational MST program is rarely a one-shot exercise. Test
developers typically draft several **candidate panel designs** that
differ in the number of stages and modules, the difficulty targeting of
each module, the routing rule (fixed cut scores vs. an
IRT-information-based router such as b-matching or maximum Fisher
information), and the scoring method used at each stage. Before
committing to one design operationally, they need a way to judge which
candidate best serves the purpose of the testing program: a licensure
exam may need to maximize precision near the pass/fail cut score,
whereas a broad-range achievement test may need uniform precision across
the entire ability continuum. Comparing candidate panels against goals
like these has long relied on simulating how each design would behave
with real examinees (Luecht & Nungester, 1998; Zenisky et al., 2009):
generate examinees with known true abilities, route each one through a
candidate panel exactly as a real test taker would be routed, and
examine how closely the resulting ability estimates recover the true
values.

**irtQ** provides two functions for exactly this purpose.
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
simulates a full MST administration — examinee by examinee, stage by
stage — and reports the resulting ability estimates and routing
pathways.
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
conducts the MST evaluation analytically, via a recursion-based
method(Lim et al., 2021), without simulating any individual examinee.
This vignette covers both, starting with
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md). They
describe an MST panel using the same structural inputs, introduced next.

### Input Data: The Building Blocks of an MST Panel

[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md) and
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
describe an MST panel using the same core structural inputs, in addition
to the item metadata:

#### 1. Item Bank (`x`)

A standard **irtQ** item metadata data frame (see
[`?shape_df`](https://hwangQ.github.io/irtQ/reference/shape_df.md) or
the *Getting Started* vignette). Each row describes one item in the item
bank. The columns are `id`, `cats`, `model`, `par.1`, `par.2`, etc.

#### 2. Module Assignment Matrix (`module`)

A binary matrix with **rows = items** (same order as `x`) and **columns
= modules**. An entry of 1 in row $`i`$, column $`m`$ means item $`i`$
belongs to module $`m`$.

            M1  M2  M3  M4  M5  M6  M7
    item 1 [  1   0   0   0   0   0   0 ]  ← item 1 is in module 1
    item 2 [  1   0   0   0   0   0   0 ]
    ...
    item 9 [  0   1   0   0   0   0   0 ]  ← item 9 is in module 2
    ...

Each item belongs to exactly one module, so each row has exactly one 1
and all other entries are 0.

#### 3. Route Map (`route_map`)

A binary **square matrix** of dimension (total modules × total modules).
An entry of 1 in row $`i`$, column $`j`$ means test takers can be routed
from module $`i`$ directly to module $`j`$.

For a 1-3-3 panel with 7 modules (M1–M7):

             M1 M2 M3 M4 M5 M6 M7
    M1  ──→ [  0  1  1  1  0  0  0 ]   M1 routes to M2, M3, or M4
    M2  ──→ [  0  0  0  0  1  1  0 ]   M2 routes to M5 or M6
    M3  ──→ [  0  0  0  0  0  1  1 ]   M3 routes to M6 or M7
    M4  ──→ [  0  0  0  0  0  1  1 ]   M4 routes to M6 or M7
    M5  ──→ [  0  0  0  0  0  0  0 ]   terminal (Stage 3)
    M6  ──→ [  0  0  0  0  0  0  0 ]   terminal
    M7  ──→ [  0  0  0  0  0  0  0 ]   terminal

Stage 1 modules are identified automatically as those with all-zero
*columns* (no module routes *to* them). Terminal modules have all-zero
*rows*.

#### 4. Cut Scores (`cut_score`)

A list of numeric vectors — one element per routing stage transition.
Each vector contains the IRT $`\theta`$ cut points used to determine
which next-stage module a test taker receives.

For a 1-3-3 panel with 2 routing stages:

``` r

cut_score = list(
  c(-0.5, 0.5),   # Stage 1 → Stage 2: θ̂ < -0.5 → M2 (easy)
                  #                    -0.5 ≤ θ̂ < 0.5 → M3 (medium)
                  #                    θ̂ ≥ 0.5 → M4 (hard)
  c(-0.6, 0.6)    # Stage 2 → Stage 3: θ̂ < -0.6 → easy module
                  #                    -0.6 ≤ θ̂ < 0.6 → medium module
                  #                    θ̂ ≥ 0.6 → hard module
)
```

`cut_score` is required by
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md),
which always uses inverse-TCC scoring internally, and by
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
whenever `route_method = NULL`. When `route_method` is instead `"bmat"`
or `"mfi"`,
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
ignores `cut_score` entirely and routes from item difficulty or test
information directly, as shown below.

#### 5. True Ability or Observed Responses (`theta` / `response`)

[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
always treats `theta` as a fixed ability *grid* at which to evaluate
bias and CSEM analytically.
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
instead uses `theta` (or `response`) to simulate individual examinees:

- Supplying `theta`: a numeric vector of *true* ability values, one per
  simulated examinee.
  [`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
  generates each examinee’s item responses internally from these true
  values, `x`, and `D`.
- Supplying `response` instead: an already-generated $`N \times J`$ item
  response matrix (e.g., from
  [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md), or
  real observed data). `theta` can still be supplied alongside
  `response`, purely for RMSE comparison after the fact — it plays no
  role in generating the responses themselves.

How `theta` is constructed changes what the simulation tells you.
Drawing it randomly from an assumed population ability distribution —
`rnorm(1000, 0, 1)` in the examples below — evaluates the panel’s
**overall** performance across a single, population-representative batch
of examinees, much as a testing program would experience it
operationally. Fixing `theta` at one value and simulating many
replications at that single point instead evaluates the panel’s
**conditional** bias and precision at that specific ability level — the
same quantity
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
computes analytically, without any simulation. Repeating that
single-point simulation across a grid of ability values is exactly the
traditional Monte Carlo approach described later in this article, and
Example 6 demonstrates it directly with
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md).

#### 6. Scoring Specification (`route_score` / `final_score`)

[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
separates the scoring method used for routing decisions at intermediate
stages (`route_score`) from the method used to report each examinee’s
final score (`final_score`). Both arguments take a named list whose
`method` element selects one of **irtQ**’s ability estimators — `"ML"`,
`"WL"`, `"MAP"`, `"EAP"`, `"EAP.SUM"`, or `"INV.TCC"` — and whose
remaining elements supply that estimator’s own arguments, e.g. `range`
for `"ML"`/`"WL"`, or `norm.prior` and `nquad` for `"EAP"`:

``` r

route_score = list(method = "EAP", norm.prior = c(0, 1), nquad = 41)
final_score = list(method = "ML", range = c(-4, 4))
```

Routing and final scoring need not use the same method, and the more
principled choice runs this way: `"EAP"` is usually preferred for
routing, because the intermediate ability estimate must always be
well-defined, even from a short module. The Stage 1 routing module in
particular may have so few items that an examinee answers every one
correctly or every one incorrectly — a pattern for which `"ML"` has no
finite solution. `"EAP"`’s prior keeps the routing estimate finite and
stable no matter the response pattern. For the final reported score,
`"ML"` is often preferred instead: by the last stage examinees have
answered enough items that such extreme patterns are rare, and `"ML"`
avoids the shrinkage toward the population mean that `"EAP"`’s prior
introduces. `"INV.TCC"` recovers ability from the observed sum score via
the inverse Test Characteristic Curve — the same method
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
uses internally — and is useful when a testing program reports summed
scores directly.

### The `simMST` Dataset

**irtQ** includes `simMST`, a built-in dataset that packages all four
inputs described above. This dataset was used in the simulation study of
Lim et al. (2021) and represents a **1-3-3 MST panel** with the
following characteristics:

- 7 modules across 3 stages
- 8 items per module (24 items total)
- All items follow the 3-parameter logistic model (3PLM)
- Item parameters calibrated to span a broad ability range

``` r

library(irtQ)

# Inspect the simMST dataset
str(simMST, max.level = 1)
#> List of 5
#>  $ item_bank:'data.frame':   300 obs. of  6 variables:
#>  $ module   : num [1:300, 1:7] 0 0 0 1 0 0 1 0 0 0 ...
#>  $ route_map:'data.frame':   7 obs. of  7 variables:
#>  $ cut_score:List of 2
#>  $ theta    : num [1:81] -4 -3.9 -3.8 -3.7 -3.6 -3.5 -3.4 -3.3 -3.2 -3.1 ...
```

``` r

# Item bank: standard irtQ metadata format
head(simMST$item_bank, 10)
#>    id cats model     par.1       par.2      par.3
#> 1   1    2  3PLM 0.8750089 -0.01501279 0.03244385
#> 2   2    2  3PLM 0.8695412  2.11748844 0.13982330
#> 3   3    2  3PLM 0.9691471  1.56938982 0.04946834
#> 4   4    2  3PLM 1.5257721 -1.03235839 0.10875097
#> 5   5    2  3PLM 0.6737973 -0.65712968 0.05493992
#> 6   6    2  3PLM 1.0763360  1.23969146 0.05258577
#> 7   7    2  3PLM 1.4353704  0.91128180 0.06818777
#> 8   8    2  3PLM 0.8761569  0.68568074 0.08907224
#> 9   9    2  3PLM 1.1352444  0.60960604 0.16278267
#> 10 10    2  3PLM 1.3357048 -1.22224352 0.14586104
```

``` r

# Module matrix: 56 items × 7 modules
dim(simMST$module)
#> [1] 300   7

# First 16 rows (items 1-16, first 2 modules)
simMST$module[1:16, ]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#>  [1,]    0    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0    0
#>  [4,]    1    0    0    0    0    0    0
#>  [5,]    0    0    0    0    0    0    0
#>  [6,]    0    0    0    0    0    0    0
#>  [7,]    1    0    0    0    0    0    0
#>  [8,]    0    0    0    0    0    0    0
#>  [9,]    0    0    0    0    0    0    0
#> [10,]    0    0    0    0    0    0    0
#> [11,]    0    0    0    0    0    1    0
#> [12,]    0    0    0    0    0    0    0
#> [13,]    0    0    0    0    0    0    0
#> [14,]    0    0    0    0    0    0    0
#> [15,]    0    0    0    0    0    0    0
#> [16,]    0    0    0    0    0    0    0
```

``` r

# Route map: 7 × 7 transition matrix
simMST$route_map
#>   V1 V2 V3 V4 V5 V6 V7
#> 1  0  1  1  1  0  0  0
#> 2  0  0  0  0  1  1  0
#> 3  0  0  0  0  1  1  1
#> 4  0  0  0  0  0  1  1
#> 5  0  0  0  0  0  0  0
#> 6  0  0  0  0  0  0  0
#> 7  0  0  0  0  0  0  0
```

The `route_map` shows:

- Row 1 (Module 1, Stage 1): routes to columns 2, 3, 4 — the three Stage
  2 modules
- Rows 2–4 (Modules 2–4, Stage 2): each routes to two or three Stage 3
  modules
- Rows 5–7 (Modules 5–7, Stage 3): all zeros — terminal modules

``` r

# Cut scores: 2 routing transitions for 3 stages
simMST$cut_score
#> [[1]]
#> [1] -0.3954891  0.4422893
#> 
#> [[2]]
#> [1] -0.6611704  0.5130905
```

``` r

# Ability grid for evaluation
length(simMST$theta)
#> [1] 81
range(simMST$theta)
#> [1] -4  4
head(simMST$theta, 10)
#>  [1] -4.0 -3.9 -3.8 -3.7 -3.6 -3.5 -3.4 -3.3 -3.2 -3.1
```

### Simulating an MST Administration with `run_mst()`

**irtQ** provides
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md) to
carry out exactly this kind of simulation. Given an item bank, a
`route_map`, a `module` matrix, and a vector of true abilities (or a
pre-generated response matrix),
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
simulates the full MST administration for every examinee —
stage-by-stage routing and final scoring — and returns the resulting
ability estimates together with the module pathway each examinee
traveled.

The examples below use the same `simMST` 1-3-3 panel introduced above: a
routing module at Stage 1, followed by three modules of varying
difficulty at each of Stages 2 and 3.

[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
supports three routing strategies. Two require no fixed cut scores at
all: the default, `route_method = "bmat"`, routes each examinee to the
reachable module whose mean item difficulty is closest to their current
ability estimate, and `route_method = "mfi"` routes to the reachable
module with the highest test information function (TIF) value at their
current estimate. The third, cut-score routing (`route_method = NULL`),
compares the current estimate against a fixed set of cut scores. Start
with the default, `bmat`:

``` r

# Item bank and panel structure from simMST
x         <- simMST$item_bank
module    <- simMST$module
route_map <- simMST$route_map

# 1,000 simulated examinees with true abilities from N(0, 1)
set.seed(42)
theta_true <- rnorm(1000, mean = 0, sd = 1)

# bmat routing: no cut scores needed
sim_bmat <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true,
  D            = 1.702,
  route_method = "bmat",
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = TRUE
)
print(sim_bmat)
#> 
#> Call:
#> run_mst(x = x, route_map = route_map, module = module, theta = theta_true, 
#>     D = 1.702, route_method = "bmat", route_score = list(method = "EAP", 
#>         norm.prior = c(0, 1), nquad = 41), final_score = list(method = "ML", 
#>         range = c(-4, 4)), se = TRUE) 
#> 
#> MST Simulation Results
#> ======================================== 
#> 
#> Panel structure:
#>   Stages              : 3
#>   Modules per stage   : 1 - 3 - 3
#>   Valid pathways      : 7
#>   Routing method      : bmat (b-matching)
#> 
#> Number of examinees : 1000
#> 
#> Ability estimation:
#>   Routing method      : EAP
#>   Final method        : ML
#> 
#> Final ability estimates (est.theta):
#>   Mean : 0.002
#>   SD   : 1.121
#>   Min  : -4.000
#>   Max  : 4.000
#> 
#> Estimation accuracy (est.theta - true.theta):
#>   Bias : 0.028
#>   RMSE : 0.345
#> 
#> Module frequency by stage:
#>   Stage 1: Module 1: 1000 (100.0%)
#>   Stage 2: Module 2: 397 (39.7%),  Module 3: 436 (43.6%),  Module 4: 167 (16.7%)
#>   Stage 3: Module 5: 264 (26.4%),  Module 6: 341 (34.1%),  Module 7: 395 (39.5%)
```

A closely related alternative is **maximum Fisher information (MFI)
routing** (`route_method = "mfi"`): instead of comparing mean item
difficulty, it routes each examinee to the reachable module that
provides the highest TIF value at their current ability estimate.

``` r

sim_mfi <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true,
  D            = 1.702,
  route_method = "mfi",
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = TRUE
)
print(sim_mfi)
#> 
#> Call:
#> run_mst(x = x, route_map = route_map, module = module, theta = theta_true, 
#>     D = 1.702, route_method = "mfi", route_score = list(method = "EAP", 
#>         norm.prior = c(0, 1), nquad = 41), final_score = list(method = "ML", 
#>         range = c(-4, 4)), se = TRUE) 
#> 
#> MST Simulation Results
#> ======================================== 
#> 
#> Panel structure:
#>   Stages              : 3
#>   Modules per stage   : 1 - 3 - 3
#>   Valid pathways      : 7
#>   Routing method      : mfi (Maximum Fisher Information)
#> 
#> Number of examinees : 1000
#> 
#> Ability estimation:
#>   Routing method      : EAP
#>   Final method        : ML
#> 
#> Final ability estimates (est.theta):
#>   Mean : -0.020
#>   SD   : 1.048
#>   Min  : -4.000
#>   Max  : 4.000
#> 
#> Estimation accuracy (est.theta - true.theta):
#>   Bias : 0.006
#>   RMSE : 0.292
#> 
#> Module frequency by stage:
#>   Stage 1: Module 1: 1000 (100.0%)
#>   Stage 2: Module 2: 377 (37.7%),  Module 3: 348 (34.8%),  Module 4: 275 (27.5%)
#>   Stage 3: Module 5: 247 (24.7%),  Module 6: 456 (45.6%),  Module 7: 297 (29.7%)
```

MFI routing assumes each module’s TIF curve behaves as intended —
peaking higher for harder modules as theta increases, and higher for
easier modules as theta decreases. In practice this assumption can fail:
a harder module’s TIF curve may unexpectedly peak at a *low* theta
value, so that low-ability examinees who should be routed to an easy
module instead get routed into the harder one. This is called an
**anomalous routing**, or **path reversal**, and it is the opposite of
what an MST is designed to do.

Alternatively, setting `route_method = NULL` switches to traditional
**cut-score routing** — the design most commonly used in operational MST
programs — where each examinee’s intermediate ability estimate is
compared against a fixed set of cut scores to decide the next module:

``` r

# Fixed routing cut scores bundled with simMST
cut_score <- simMST$cut_score   # list(c(-0.5, 0.5), c(-0.6, 0.6))

sim_cut <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true,
  D            = 1.702,
  route_method = NULL,
  cut_score    = cut_score,
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = TRUE
)
print(sim_cut)
#> 
#> Call:
#> run_mst(x = x, route_map = route_map, module = module, theta = theta_true, 
#>     D = 1.702, route_method = NULL, cut_score = cut_score, route_score = list(method = "EAP", 
#>         norm.prior = c(0, 1), nquad = 41), final_score = list(method = "ML", 
#>         range = c(-4, 4)), se = TRUE) 
#> 
#> MST Simulation Results
#> ======================================== 
#> 
#> Panel structure:
#>   Stages              : 3
#>   Modules per stage   : 1 - 3 - 3
#>   Valid pathways      : 7
#>   Routing method      : cut-score based
#> 
#> Number of examinees : 1000
#> 
#> Ability estimation:
#>   Routing method      : EAP
#>   Final method        : ML
#> 
#> Final ability estimates (est.theta):
#>   Mean : 0.012
#>   SD   : 1.063
#>   Min  : -4.000
#>   Max  : 4.000
#> 
#> Estimation accuracy (est.theta - true.theta):
#>   Bias : 0.038
#>   RMSE : 0.314
#> 
#> Module frequency by stage:
#>   Stage 1: Module 1: 1000 (100.0%)
#>   Stage 2: Module 2: 351 (35.1%),  Module 3: 348 (34.8%),  Module 4: 301 (30.1%)
#>   Stage 3: Module 5: 248 (24.8%),  Module 6: 368 (36.8%),  Module 7: 384 (38.4%)
```

#### Deriving Principled Cut Scores with `find_cut()`

The cut scores used above (`c(-0.5, 0.5)` and `c(-0.6, 0.6)`) were
chosen heuristically when the `simMST` panel was built — a common
practice, but one that offers no guarantee that the cut points actually
fall where they should psychometrically. A more principled choice is to
derive the cut score directly from the modules’ TIF curves themselves:
the natural boundary between an easy and a hard module is the theta
value at which the hard module’s TIF first overtakes the easy module’s
TIF. Below that point the easier module is more informative; above it,
the harder module is. This is also exactly the fix for the
anomalous-routing problem described above for MFI routing — by fixing
the cut point at this *proper* crossing, low-ability examinees can never
be routed into the harder module no matter how its TIF curve behaves at
the extremes.

[`find_cut()`](https://hwangQ.github.io/irtQ/reference/find_cut.md)
automates this search. For every pair of adjacent modules at each stage
transition, it scans the theta scale for crossings between the two
modules’ TIF curves, keeps only the “proper” crossing (the TIF
difference changing sign from negative to positive), and returns the
result as a `cut_score` list that can be passed directly to
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md):

``` r

cut_result <- find_cut(
  x         = x,
  module    = module,
  route_map = route_map
)
print(cut_result)
#> MST TIF-Crossing Cut Score Results
#> ==================================================
#> 
#> stage.2
#> -----------------------------------
#>   Modules (index order)    : 2, 3, 4
#>   Modules (difficulty order): 2, 4, 3
#>   Mean item locations      : -0.7313, 0.1342, 0.1279
#>   Selected cut score(s)    : -0.3955, 0.4423
#>   Pair details:
#>     [mod2_vs_mod3]
#>       Proper crossing(s)    : 2  [theta = -0.3955, 3.0936]
#>       Anomalous crossing(s) : 1  [theta = 2.1328]  (excluded)
#>       Selected cut score    : -0.3955
#>     [mod3_vs_mod4]
#>       Proper crossing(s)    : 1  [theta = 0.4423]
#>       Anomalous crossing(s) : 2  [theta = -1.3914, 5.7208]  (excluded)
#>       Selected cut score    : 0.4423
#> 
#> stage.3
#> -----------------------------------
#>   Modules (index order)    : 5, 6, 7
#>   Modules (difficulty order): 5, 6, 7
#>   Mean item locations      : -0.8563, -0.1402, 0.4899
#>   Selected cut score(s)    : -0.6612, 0.5131
#>   Pair details:
#>     [mod5_vs_mod6]
#>       Proper crossing(s)    : 1  [theta = -0.6612]
#>       Selected cut score    : -0.6612
#>     [mod6_vs_mod7]
#>       Proper crossing(s)    : 1  [theta = 0.5131]
#>       Anomalous crossing(s) : 1  [theta = -1.7322]  (excluded)
#>       Selected cut score    : 0.5131
#> ==================================================
#> Usage: run_mst(..., route_method = NULL, cut_score = <result>$cut_score)

# Compare with the heuristic cut scores bundled in simMST
cut_result$cut_score
#> $stage.2
#> [1] -0.3954894  0.4422831
#> 
#> $stage.3
#> [1] -0.6611766  0.5130904
simMST$cut_score
#> [[1]]
#> [1] -0.3954891  0.4422893
#> 
#> [[2]]
#> [1] -0.6611704  0.5130905
```

Calling [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
`find_cut` object visualizes each stage’s module TIF curves, with the
selected cut score drawn as a solid vertical line and any excluded
anomalous crossings drawn as dashed red lines:

``` r

plot(cut_result)
```

![TIF curves for every module in the simMST panel, faceted by stage,
with TIF-crossing cut scores from
find_cut().](mst-panel-evaluation_files/figure-html/find-cut-plot-1.png)

TIF curves for every module in the simMST panel, faceted by stage, with
TIF-crossing cut scores from find_cut().

The derived cut scores plug into
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
exactly like any other fixed cut scores, by setting
`route_method = NULL`:

``` r

sim_findcut <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true,
  D            = 1.702,
  route_method = NULL,
  cut_score    = cut_result$cut_score,
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = TRUE
)
print(sim_findcut)
#> 
#> Call:
#> run_mst(x = x, route_map = route_map, module = module, theta = theta_true, 
#>     D = 1.702, route_method = NULL, cut_score = cut_result$cut_score, 
#>     route_score = list(method = "EAP", norm.prior = c(0, 1), 
#>         nquad = 41), final_score = list(method = "ML", range = c(-4, 
#>         4)), se = TRUE) 
#> 
#> MST Simulation Results
#> ======================================== 
#> 
#> Panel structure:
#>   Stages              : 3
#>   Modules per stage   : 1 - 3 - 3
#>   Valid pathways      : 7
#>   Routing method      : cut-score based
#> 
#> Number of examinees : 1000
#> 
#> Ability estimation:
#>   Routing method      : EAP
#>   Final method        : ML
#> 
#> Final ability estimates (est.theta):
#>   Mean : -0.013
#>   SD   : 1.107
#>   Min  : -4.000
#>   Max  : 4.000
#> 
#> Estimation accuracy (est.theta - true.theta):
#>   Bias : 0.012
#>   RMSE : 0.327
#> 
#> Module frequency by stage:
#>   Stage 1: Module 1: 1000 (100.0%)
#>   Stage 2: Module 2: 376 (37.6%),  Module 3: 335 (33.5%),  Module 4: 289 (28.9%)
#>   Stage 3: Module 5: 253 (25.3%),  Module 6: 379 (37.9%),  Module 7: 368 (36.8%)
```

#### Comparing Routing Methods on Identical Responses with `response`

All of the calls above let
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
simulate item responses internally from `theta`. When the goal is to
compare two routing or scoring strategies fairly, simulation noise in
the responses themselves should not be allowed to confound the
comparison. Generating the response matrix once with
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) and
feeding it through the `response` argument removes that source of noise
— and is also how
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md) would
be used with real, already-observed response data:

``` r

# Generate one response matrix from a fresh set of true abilities
set.seed(123)
theta_true2 <- rnorm(500, mean = 0, sd = 1)
resp_matrix <- simdat(x = x, theta = theta_true2, D = 1.702)
# resp_matrix is an N x J matrix (N examinees x J total items in the bank)

# bmat routing on the pre-generated responses
result_A <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true2,    # kept only for RMSE evaluation below
  response     = resp_matrix,    # pre-generated N x J response matrix
  D            = 1.702,
  route_method = "bmat",
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = FALSE
)

# find_cut()-based cut-score routing on the SAME responses
result_B <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = theta_true2,
  response     = resp_matrix,    # identical responses -> fair comparison
  D            = 1.702,
  route_method = NULL,
  cut_score    = cut_result$cut_score,
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = FALSE
)

# Any RMSE difference now reflects the routing rule alone, not response noise
rmse_A <- sqrt(mean((result_A$est.theta - theta_true2)^2))
rmse_B <- sqrt(mean((result_B$est.theta - theta_true2)^2))
cat(sprintf("RMSE (bmat):                 %.4f\n", rmse_A))
#> RMSE (bmat):                 0.3706
cat(sprintf("RMSE (find_cut cut scores):  %.4f\n", rmse_B))
#> RMSE (find_cut cut scores):  0.3432
```

All five calls above return the same kind of result regardless of the
routing rule or input format used: a final ability estimate
(`est.theta`), its standard error (`se.theta`), and the module pathway
taken by each examinee (`path`).

#### Inspecting the `run_mst()` Output

[`print()`](https://rdrr.io/r/base/print.html) reports only aggregated
information — the panel structure, the routing and final scoring
methods, summary statistics (mean, SD, min, max) of the final estimates,
and, when true abilities are supplied, overall bias and RMSE. Every
examinee’s individual results live in named components of the returned
object itself:

``` r

# Final ability estimate and its standard error, one value per examinee
head(sim_bmat$est.theta)
#> [1]  0.91354165 -0.45702570  0.42244620  0.76022539  0.45680206  0.07429823
head(sim_bmat$se.theta)
#> [1] 0.2717126 0.2503814 0.2453777 0.2596967 0.2480902 0.2542344

# Routing and final estimates at every stage: an N x n.stage matrix
head(sim_bmat$theta.route)
#>         stage.1     stage.2     stage.3
#> [1,]  1.7362572  0.58220759  0.91354165
#> [2,] -0.5054051 -0.19097482 -0.45702570
#> [3,]  1.0914435 -0.02701402  0.42244620
#> [4,]  1.7362572  0.23973623  0.76022539
#> [5,]  0.3594823  0.58220759  0.45680206
#> [6,]  0.3594823  0.28133615  0.07429823

# Module index administered at every stage: an N x n.stage matrix
head(sim_bmat$path)
#>      stage.1 stage.2 stage.3
#> [1,]       1       3       7
#> [2,]       1       2       6
#> [3,]       1       3       6
#> [4,]       1       3       7
#> [5,]       1       3       7
#> [6,]       1       3       7
```

`se.theta` in particular never appears in
[`print()`](https://rdrr.io/r/base/print.html) output, so accessing it
directly is the only way to see how estimation precision varies across
examinees — for example, by the Stage 3 (final) module each examinee
landed in:

``` r

# Mean standard error by terminal (Stage 3) module
tapply(sim_bmat$se.theta, sim_bmat$path[, ncol(sim_bmat$path)], mean)
#>         5         6         7 
#> 1.8259796 0.2648599 3.5765081
```

The full list of return components — including `panel` (reused from
[`panel_info()`](https://hwangQ.github.io/irtQ/reference/panel_info.md)),
`true.theta`, and `full.resp` (the complete item-level response matrix,
populated only when `return_full_resp = TRUE`) — is summarized in the
Function Reference section near the end of this article.

All five calls above drew `theta_true` (or `theta_true2`) from a single
population distribution, $`N(0, 1)`$ — appropriate for an **overall**,
population-representative look at how the panel performs across a
realistic mix of examinees, the way a testing program would experience
it operationally. The complementary question — how the panel performs
**conditionally**, at one specific ability level — requires fixing
$`\theta`$ instead of sampling it, and repeating the simulation many
times at that single point. Doing this across a whole grid of ability
levels is exactly the traditional approach described next.

------------------------------------------------------------------------

## The Challenge of Evaluating MST Panels

Before deploying an MST panel operationally, test developers need to
evaluate its **measurement quality**: how accurately and precisely does
the panel estimate examinee ability at each level of the latent trait
$`\theta`$?

The two key evaluation metrics are:

- **Conditional Bias**: The average difference between the ability
  estimate $`\hat{\theta}`$ and the true ability $`\theta`$, evaluated
  at each $`\theta`$ level. Ideally, bias should be near zero across the
  full ability range.
- **Conditional Standard Error of Measurement (CSEM)**: The standard
  deviation of ability estimates at each $`\theta`$ level. Lower CSEM
  means more precise measurement.

### The Traditional Approach: Monte Carlo Simulation

Conditional evaluation via simulation means repeating the
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
procedure shown above systematically across a fixed grid of ability
levels, with many replications at each grid point:

1.  Fix a set of true ability levels $`\theta_1, \theta_2, \ldots`$
2.  For each $`\theta_k`$, generate thousands of simulated response
    patterns
3.  Route each simulated examinee through the panel — exactly as
    [`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
    does above — using the chosen routing rule (cut scores, b-matching,
    or maximum information)
4.  Estimate ability from each simulated response and compute the mean
    and variance of the estimates
5.  Bias = mean($`\hat\theta`$) − $`\theta_k`$, CSEM =
    $`\sqrt{\text{Var}(\hat\theta)}`$

[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
already automates Steps 2–4 for a single $`\theta_k`$ and a single batch
of replications; a full evaluation just wraps it in a loop over the
ability grid with a sufficiently large number of replications per grid
point (Example 6 below demonstrates exactly this). This approach is
conceptually simple, but looping
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md) this
way inherits several practical drawbacks:

- **Computationally inefficient**: Thousands of response patterns must
  be generated and scored
- **Stochastic**: Results fluctuate across simulation replications
- **Time-consuming**: Evaluating many panel designs or cut score
  configurations requires many separate simulation runs

### The Recursion-Based Analytical Approach

Lim et al. (2021) proposed a fundamentally different approach: instead
of simulating individual examinees, the method **directly computes the
exact probability distribution of every possible observed score** at
every stage and pathway using a recursive algorithm.

The key insight is that the conditional distribution of the observed sum
score along any pathway can be built up stage by stage using the
**Lord–Wingersky recursion** (Lord & Wingersky, 1984). Given the
conditional score distribution of the modules visited so far, the joint
distribution at the next stage can be computed exactly — without any
random sampling.

Here is the core logic of the recursion:

1.  **Stage 1**: Compute $`P(X_1 = x \mid \theta)`$ for the routing
    module using the Lord–Wingersky recursion, where $`X_1`$ is the
    total score on Stage 1.

2.  **Routing**: For each possible score $`x`$ on Stage 1, determine
    which Stage 2 module a test taker would be routed to (using the cut
    scores). This converts the score distribution into a *pathway
    probability*.

3.  **Stage 2**: For each Stage 2 pathway, compute the **joint
    distribution** $`P(X_1 + X_2 = s \mid \theta, \text{path})`$ by
    convolving the Stage 1 score distribution with the conditional
    distribution of the assigned Stage 2 module.

4.  **Continue recursively** through all stages.

5.  **Ability estimation**: Convert each possible final sum score to a
    $`\hat\theta`$ estimate using the **inverse Test Characteristic
    Curve (TCC) method** — the method implied by IRT-based summed
    scoring.

    > **Note on linear interpolation (`intpol`)**: When items have
    > non-zero guessing parameters (e.g., 3PLM), the minimum
    > *achievable* expected sum score exceeds zero — meaning no valid
    > $`\hat\theta`$ exists for very low observed scores (those below
    > the sum of guessing parameters). With `intpol = TRUE` (the
    > default), the inverse TCC method applies **linear interpolation**
    > between the point $`(\theta_{\min},\, 0)`$ and the lowest
    > TCC-reachable score point, so that every possible sum score
    > receives a valid ability estimate rather than `NA`.

6.  **Evaluate**: Compute the conditional mean and variance of
    $`\hat\theta`$ given the true $`\theta`$, and derive conditional
    bias and CSEM.

This method is:

- **Exact**: The probability distributions are computed analytically,
  not approximated by random sampling
- **Fast**: Computation takes seconds, not minutes or hours
- **Deterministic**: Results are reproducible without any simulation
  noise
- **Comprehensive**: Any panel design, cut score configuration, or
  ability grid can be evaluated in a single function call

The
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
function in **irtQ** implements this recursion-based method, computing
the same conditional bias and CSEM that the Monte Carlo loop above
estimates by simulation — but analytically, in a single function call
and without any simulation noise.

------------------------------------------------------------------------

## Example 1: Evaluating the simMST Panel

With all inputs in place, running
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md) is
straightforward:

``` r

# Extract components from simMST
x         <- simMST$item_bank
module    <- simMST$module
route_map <- simMST$route_map
cut_score <- simMST$cut_score
theta     <- simMST$theta

# Evaluate the 1-3-3 MST panel
eval_result <- reval_mst(
  x          = x,
  D          = 1.702,
  route_map  = route_map,
  module     = module,
  cut_score  = cut_score,
  theta      = theta,
  range.tcc  = c(-5, 5)
)
```

The function returns a named list with 7 components. The most important
is `eval.tb`, the evaluation summary table:

``` r

# Evaluation table: theta, mu, sigma2, bias, csem
print(eval_result$eval.tb)
#>    theta           mu     sigma2          bias      csem
#> 1   -4.0 -3.315151711 0.81243639  0.6848482889 0.9013525
#> 2   -3.9 -3.308122145 0.81249024  0.5918778555 0.9013824
#> 3   -3.8 -3.299466478 0.81250415  0.5005335222 0.9013901
#> 4   -3.7 -3.288817058 0.81244163  0.4111829417 0.9013554
#> 5   -3.6 -3.275730371 0.81224431  0.3242696286 0.9012460
#> 6   -3.5 -3.259675887 0.81182079  0.2403241126 0.9010110
#> 7   -3.4 -3.240025903 0.81103090  0.1599740975 0.9005725
#> 8   -3.3 -3.216048229 0.80966409  0.0839517714 0.8998134
#> 9   -3.2 -3.186904502 0.80741090  0.0130954981 0.8985605
#> 10  -3.1 -3.151657917 0.80382781 -0.0516579167 0.8965644
#> 11  -3.0 -3.109295249 0.79829771 -0.1092952491 0.8934751
#> 12  -2.9 -3.058768592 0.78999264 -0.1587685917 0.8888153
#> 13  -2.8 -2.999061853 0.77785149 -0.1990618528 0.8819589
#> 14  -2.7 -2.929284858 0.76059315 -0.2292848582 0.8721199
#> 15  -2.6 -2.848793209 0.73679246 -0.2487932093 0.8583662
#> 16  -2.5 -2.757324677 0.70504694 -0.2573246765 0.8396707
#> 17  -2.4 -2.655133580 0.66425060 -0.2551335795 0.8150157
#> 18  -2.3 -2.543095304 0.61395916 -0.2430953040 0.7835555
#> 19  -2.2 -2.422746810 0.55477998 -0.2227468100 0.7448355
#> 20  -2.1 -2.296229108 0.48866012 -0.1962291085 0.6990423
#> 21  -2.0 -2.166108358 0.41890796 -0.1661083585 0.6472310
#> 22  -1.9 -2.035077758 0.34981064 -0.1350777582 0.5914479
#> 23  -1.8 -1.905584122 0.28583952 -0.1055841216 0.5346396
#> 24  -1.7 -1.779470568 0.23064795 -0.0794705684 0.4802582
#> 25  -1.6 -1.657751636 0.18624743 -0.0577516359 0.4315639
#> 26  -1.5 -1.540606085 0.15273247 -0.0406060853 0.3908100
#> 27  -1.4 -1.427583578 0.12865527 -0.0275835775 0.3586855
#> 28  -1.3 -1.317925176 0.11180081 -0.0179251758 0.3343663
#> 29  -1.2 -1.210862682 0.09995113 -0.0108626816 0.3161505
#> 30  -1.1 -1.105805063 0.09135406 -0.0058050633 0.3022483
#> 31  -1.0 -1.002391488 0.08486148 -0.0023914881 0.2913099
#> 32  -0.9 -0.900435000 0.07985646 -0.0004349999 0.2825889
#> 33  -0.8 -0.799797656 0.07608715  0.0002023439 0.2758390
#> 34  -0.7 -0.700256814 0.07346730 -0.0002568144 0.2710485
#> 35  -0.6 -0.601437729 0.07189144 -0.0014377292 0.2681258
#> 36  -0.5 -0.502856126 0.07113257 -0.0028561262 0.2667069
#> 37  -0.4 -0.404036888 0.07086214 -0.0040368885 0.2661994
#> 38  -0.3 -0.304627592 0.07074799 -0.0046275922 0.2659849
#> 39  -0.2 -0.204456986 0.07053962 -0.0044569862 0.2655930
#> 40  -0.1 -0.103545413 0.07009517 -0.0035454132 0.2647549
#> 41   0.0 -0.002086355 0.06937647 -0.0020863553 0.2633941
#> 42   0.1  0.099605094 0.06845181 -0.0003949059 0.2616330
#> 43   0.2  0.201186926 0.06750554  0.0011869262 0.2598183
#> 44   0.3  0.302412683 0.06682583  0.0024126833 0.2585069
#> 45   0.4  0.403238687 0.06675136  0.0032386872 0.2583629
#> 46   0.5  0.503865394 0.06758217  0.0038653939 0.2599657
#> 47   0.6  0.604676874 0.06948120  0.0046768744 0.2635929
#> 48   0.7  0.706084755 0.07241436  0.0060847548 0.2690992
#> 49   0.8  0.808347795 0.07618349  0.0083477949 0.2760136
#> 50   0.9  0.911476680 0.08057140  0.0114766796 0.2838510
#> 51   1.0  1.015292639 0.08556051  0.0152926392 0.2925073
#> 52   1.1  1.119613635 0.09157890  0.0196136351 0.3026201
#> 53   1.2  1.224482001 0.09979670  0.0244820009 0.3159062
#> 54   1.3  1.330368035 0.11256418  0.0303680346 0.3355059
#> 55   1.4  1.438341074 0.13403659  0.0383410744 0.3661101
#> 56   1.5  1.550220555 0.17081840  0.0502205547 0.4133018
#> 57   1.6  1.668675144 0.23213521  0.0686751437 0.4818041
#> 58   1.7  1.797161158 0.32880161  0.0971611579 0.5734122
#> 59   1.8  1.939556253 0.47041363  0.1395562528 0.6858671
#> 60   1.9  2.099420004 0.66100489  0.1994200038 0.8130221
#> 61   2.0  2.278999976 0.89466646  0.2789999764 0.9458681
#> 62   2.1  2.478299894 1.15350998  0.3782998942 1.0740158
#> 63   2.2  2.694589735 1.40985468  0.4945897346 1.1873730
#> 64   2.3  2.922593610 1.63250973  0.6225936097 1.2776970
#> 65   2.4  3.155313241 1.79472536  0.7553132408 1.3396736
#> 66   2.5  3.385198212 1.88048226  0.8851982121 1.3713068
#> 67   2.6  3.605288505 1.88687626  1.0052885051 1.3736361
#> 68   2.7  3.810040895 1.82247511  1.1100408947 1.3499908
#> 69   2.8  3.995722980 1.70319022  1.1957229797 1.3050633
#> 70   2.9  4.160415693 1.54765003  1.2604156927 1.2440458
#> 71   3.0  4.303751347 1.37351113  1.3037513472 1.1719689
#> 72   3.1  4.426528159 1.19527061  1.3265281586 1.0932843
#> 73   3.2  4.530311684 1.02346632  1.3303116839 1.0116651
#> 74   3.3  4.617089274 0.86483582  1.3170892740 0.9299655
#> 75   3.4  4.689004875 0.72297630  1.2890048747 0.8502801
#> 76   3.5  4.748175629 0.59915902  1.2481756286 0.7740536
#> 77   3.6  4.796578286 0.49309273  1.1965782860 0.7022056
#> 78   3.7  4.835988848 0.40354394  1.1359888485 0.6352511
#> 79   3.8  4.867959380 0.32879397  1.0679593803 0.5734056
#> 80   3.9  4.893818703 0.26694965  0.9938187025 0.5166717
#> 81   4.0  4.914687017 0.21613853  0.9146870173 0.4649070
```

Each row corresponds to one true ability level $`\theta`$. The columns
are:

| Column | Meaning |
|----|----|
| `theta` | True ability level |
| `mu` | Conditional mean of ability estimates $`E[\hat\theta \mid \theta]`$ |
| `sigma2` | Conditional variance of ability estimates $`\text{Var}[\hat\theta \mid \theta]`$ |
| `bias` | Conditional bias = $`\mu - \theta`$ |
| `csem` | Conditional SEM = $`\sqrt{\sigma^2}`$ |

A well-designed panel will show:

- `bias` values close to zero across the full ability range
- `csem` values that are relatively small and stable (or slightly
  U-shaped, higher at the extremes where fewer items provide
  information)

------------------------------------------------------------------------

## Example 2: Visualizing Bias and CSEM

Plotting the evaluation results provides an intuitive picture of the
panel’s measurement quality.

``` r

eval_tb <- eval_result$eval.tb

# Side-by-side plots: bias (left) and CSEM (right)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Bias plot
plot(
  eval_tb$theta, eval_tb$bias,
  type = "b", pch = 16, col = "steelblue", lwd = 2,
  xlab = expression(theta),
  ylab = "Conditional Bias",
  main = "Conditional Bias",
  ylim = c(-0.3, 0.3)
)
abline(h = 0, col = "red", lty = 2, lwd = 1.5)
grid()

# CSEM plot
plot(
  eval_tb$theta, eval_tb$csem,
  type = "b", pch = 16, col = "darkorange", lwd = 2,
  xlab = expression(theta),
  ylab = "CSEM",
  main = "Conditional SEM (CSEM)",
  ylim = c(0, max(eval_tb$csem) * 1.2)
)
abline(h = mean(eval_tb$csem), col = "red", lty = 2, lwd = 1.5)
legend("topright", legend = "Mean CSEM", lty = 2, col = "red", lwd = 1.5)
grid()
```

![Conditional bias (left) and CSEM (right) for the simMST 1-3-3 panel
across the ability
scale.](mst-panel-evaluation_files/figure-html/plot-bias-csem-1.png)

Conditional bias (left) and CSEM (right) for the simMST 1-3-3 panel
across the ability scale.

``` r

par(mfrow = c(1, 1))
```

You can also create a combined `ggplot2` figure, as shown in the
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
help page:

``` r

library(ggplot2)
library(tidyr)
library(dplyr)

eval_tb %>%
  dplyr::select(theta, bias, csem) %>%
  tidyr::pivot_longer(
    cols      = c(bias, csem),
    names_to  = "criterion",
    values_to = "value"
  ) %>%
  ggplot2::ggplot(aes(x = theta, y = value)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_line(aes(colour = criterion, linetype = criterion), linewidth = 1.2) +
  ggplot2::geom_point(aes(shape = criterion), size = 2.5) +
  ggplot2::scale_colour_manual(values = c(bias = "steelblue", csem = "darkorange")) +
  ggplot2::scale_linetype_manual(values = c(bias = "solid", csem = "dashed")) +
  ggplot2::labs(
    x       = expression(theta),
    y       = NULL,
    title   = "1-3-3 MST Panel: Conditional Bias and CSEM",
    colour  = NULL,
    linetype = NULL,
    shape   = NULL
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.key.width = unit(1.5, "cm"))
```

![Bias and CSEM plotted together using
ggplot2.](mst-panel-evaluation_files/figure-html/plot-ggplot2-1.png)

Bias and CSEM plotted together using ggplot2.

------------------------------------------------------------------------

## Example 3: Exploring Other Output Components

Beyond `eval.tb`,
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
returns intermediate objects that can help diagnose how the panel
operates.

### Panel Structure (`panel.info`)

``` r

# Panel configuration: which modules belong to each stage
eval_result$panel.info$config
#> $stage.1
#> [1] 1
#> 
#> $stage.2
#> [1] 2 3 4
#> 
#> $stage.3
#> [1] 5 6 7

# All valid pathways through the panel
eval_result$panel.info$pathway
#>        stage.1 stage.2 stage.3
#> path.1       1       2       5
#> path.2       1       2       6
#> path.3       1       3       5
#> path.4       1       3       6
#> path.5       1       3       7
#> path.6       1       4       6
#> path.7       1       4       7

# Number of modules per stage
eval_result$panel.info$n.module
#> stage.1 stage.2 stage.3 
#>       1       3       3

# Total number of stages
eval_result$panel.info$n.stage
#> [1] 3
```

### Items per Module (`item.by.mod`)

``` r

# Item metadata for Module 1 (the routing module)
eval_result$item.by.mod$m.1
#>    id cats model    par.1       par.2      par.3
#> 1   4    2  3PLM 1.525772 -1.03235839 0.10875097
#> 2   7    2  3PLM 1.435370  0.91128180 0.06818777
#> 3  41    2  3PLM 1.307673  0.84073482 0.04531530
#> 4  79    2  3PLM 1.627854 -0.04267089 0.08726219
#> 5  85    2  3PLM 1.593441  0.12433962 0.04843426
#> 6  97    2  3PLM 1.425347  1.62529602 0.08793450
#> 7 120    2  3PLM 1.629243 -1.55679267 0.08627866
#> 8 178    2  3PLM 1.292352 -1.37542472 0.03876443
```

``` r

# Item metadata for Module 5 (Stage 3, first terminal module)
eval_result$item.by.mod$m.5
#>    id cats model    par.1       par.2      par.3
#> 1  20    2  3PLM 1.115701 -0.84652958 0.06801272
#> 2  31    2  3PLM 1.414366 -1.28095919 0.11286449
#> 3  42    2  3PLM 1.457161 -0.75722758 0.06916838
#> 4 109    2  3PLM 1.622737  0.06477017 0.05118737
#> 5 123    2  3PLM 1.398542 -1.18015541 0.13934292
#> 6 163    2  3PLM 1.033430 -0.88532141 0.08320088
#> 7 194    2  3PLM 1.241976 -1.93284746 0.13121892
#> 8 287    2  3PLM 1.559166 -0.03225149 0.08974932
```

### Inverse-TCC Ability Estimates (`eq.theta`)

The `eq.theta` component contains the IRT $`\theta`$ estimates
corresponding to each possible observed sum score, computed via the
inverse TCC method. These are the score-to-$`\theta`$ mappings used for
routing and final ability reporting.

``` r

# eq.theta[[stage]][[path]] gives a vector of theta estimates,
# one per possible observed sum score on that partial path

# Stage 1, Path 1 (routing module only — 8 items, scores 0-8)
cat("Theta estimates for Stage 1 (8 items, scores 0-8):\n")
#> Theta estimates for Stage 1 (8 items, scores 0-8):
round(eval_result$eq.theta$stage.1[, 1], 3)
#> [1] -5.000 -2.052 -1.335 -0.753 -0.173  0.323  0.832  1.441  5.000
```

``` r

# Stage 3 has multiple columns — one per complete pathway
cat("Dimensions of eq.theta at Stage 3 (rows = possible scores, cols = pathways):\n")
#> Dimensions of eq.theta at Stage 3 (rows = possible scores, cols = pathways):
dim(eval_result$eq.theta$stage.3)
#> [1] 25  7
```

The number of rows equals the number of possible sum scores (0 through
maximum score) for items along that pathway. Each column corresponds to
one complete pathway through the MST.

### Inspecting Test Information by Pathway

Comparing item information across pathways reveals why bias and CSEM
vary with ability level:

``` r

# Retrieve item metadata for the easiest and hardest complete pathways
# (pathway 1 = low ability; last pathway = high ability)
n_paths <- nrow(eval_result$panel.info$pathway)

meta_low  <- eval_result$item.by.path$stage.3$path.1          # low-ability path
meta_high <- eval_result$item.by.path[[3]][[n_paths]]          # high-ability path

theta_grid <- seq(-4, 4, 0.1)

tif_low  <- info(x = meta_low,  theta = theta_grid, D = 1.702)$tif
tif_high <- info(x = meta_high, theta = theta_grid, D = 1.702)$tif

par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
plot(
  theta_grid, tif_low,
  type = "l", col = "steelblue", lwd = 2,
  xlab = expression(theta), ylab = "Test Information",
  main = "Test Information: Low-Ability vs. High-Ability Pathway",
  ylim = c(0, max(c(tif_low, tif_high)) * 1.1)
)
lines(theta_grid, tif_high, col = "darkorange", lwd = 2, lty = 2)
legend(
  "topright",
  legend  = c("Low-ability pathway (path 1)", "High-ability pathway (last path)"),
  col     = c("steelblue", "darkorange"),
  lty     = c(1, 2), lwd = 2
)
grid()
```

![Test information functions for two contrasting pathways in the 1-3-3
panel.](mst-panel-evaluation_files/figure-html/pathway-info-1.png)

Test information functions for two contrasting pathways in the 1-3-3
panel.

The low-ability pathway peaks at negative $`\theta`$ values, while the
high-ability pathway peaks at positive values — exactly what good MST
design achieves.

------------------------------------------------------------------------

## Example 4: Building a Custom MST Panel

To evaluate a panel design from scratch, you need to construct `x`,
`module`, `route_map`, and `cut_score` yourself. This example shows how
to do this for a **1-2-2 MST panel** with 5 modules and 6 items per
module (30 items total).

### Design Overview

![1-2-2 MST panel design. From Stage 1, all examinees begin at module
1M. After Stage 1, they are routed to either 2E (easy) or 2H (hard).
After Stage 2, they can be routed to either 3E or 3H regardless of which
Stage 2 module they took.](img/mst_122_design.png)

1-2-2 MST panel design. From Stage 1, all examinees begin at module 1M.
After Stage 1, they are routed to either 2E (easy) or 2H (hard). After
Stage 2, they can be routed to either 3E or 3H regardless of which Stage
2 module they took.

Routing:

- Stage 1 → Stage 2: $`\hat\theta < 0`$ → M2 (Easy),
  $`\hat\theta \geq 0`$ → M3 (Hard)
- Stage 2 → Stage 3: $`\hat\theta < 0`$ → M4 (Easy),
  $`\hat\theta \geq 0`$ → M5 (Hard)

This gives **4 pathways**: M1-M2-M4, M1-M2-M5, M1-M3-M4, M1-M3-M5.

### Step 1: Build the Item Bank

Items in easy modules have lower difficulty; items in hard modules have
higher difficulty. We use the 3PLM for all items.

``` r

# 6 items per module × 5 modules = 30 items total
n_per_mod <- 6

# Helper: create item metadata for one module
make_mod_items <- function(a_mean, b_vec, g_val, id_prefix) {
  shape_df(
    par.drm  = list(
      a = rep(a_mean, n_per_mod),
      b = b_vec,
      g = rep(g_val, n_per_mod)
    ),
    item.id  = paste0(id_prefix, 1:n_per_mod),
    cats     = 2,
    model    = "3PLM"
  )
}

# Module 1 (Routing, Stage 1): moderate difficulty
items_m1 <- make_mod_items(
  a_mean   = 1.2,
  b_vec    = c(-0.5, -0.2,  0.0,  0.2,  0.4,  0.6),
  g_val    = 0.15,
  id_prefix = "M1_I"
)

# Module 2 (Easy, Stage 2): lower difficulty
items_m2 <- make_mod_items(
  a_mean   = 1.1,
  b_vec    = c(-2.0, -1.6, -1.3, -1.0, -0.7, -0.4),
  g_val    = 0.15,
  id_prefix = "M2_I"
)

# Module 3 (Hard, Stage 2): higher difficulty
items_m3 <- make_mod_items(
  a_mean   = 1.3,
  b_vec    = c( 0.4,  0.6,  0.9,  1.2,  1.5,  1.8),
  g_val    = 0.15,
  id_prefix = "M3_I"
)

# Module 4 (Easy, Stage 3): lowest difficulty
items_m4 <- make_mod_items(
  a_mean   = 1.0,
  b_vec    = c(-2.5, -2.2, -1.9, -1.6, -1.3, -1.0),
  g_val    = 0.15,
  id_prefix = "M4_I"
)

# Module 5 (Hard, Stage 3): highest difficulty
items_m5 <- make_mod_items(
  a_mean   = 1.4,
  b_vec    = c( 0.8,  1.1,  1.4,  1.7,  2.0,  2.3),
  g_val    = 0.15,
  id_prefix = "M5_I"
)

# Combine into single item bank
item_bank_122 <- dplyr::bind_rows(items_m1, items_m2, items_m3, items_m4, items_m5)
cat("Item bank dimensions:", nrow(item_bank_122), "items ×", ncol(item_bank_122), "columns\n")
#> Item bank dimensions: 30 items × 6 columns
print(item_bank_122)
#>       id cats model par.1 par.2 par.3
#> 1  M1_I1    2  3PLM   1.2  -0.5  0.15
#> 2  M1_I2    2  3PLM   1.2  -0.2  0.15
#> 3  M1_I3    2  3PLM   1.2   0.0  0.15
#> 4  M1_I4    2  3PLM   1.2   0.2  0.15
#> 5  M1_I5    2  3PLM   1.2   0.4  0.15
#> 6  M1_I6    2  3PLM   1.2   0.6  0.15
#> 7  M2_I1    2  3PLM   1.1  -2.0  0.15
#> 8  M2_I2    2  3PLM   1.1  -1.6  0.15
#> 9  M2_I3    2  3PLM   1.1  -1.3  0.15
#> 10 M2_I4    2  3PLM   1.1  -1.0  0.15
#> 11 M2_I5    2  3PLM   1.1  -0.7  0.15
#> 12 M2_I6    2  3PLM   1.1  -0.4  0.15
#> 13 M3_I1    2  3PLM   1.3   0.4  0.15
#> 14 M3_I2    2  3PLM   1.3   0.6  0.15
#> 15 M3_I3    2  3PLM   1.3   0.9  0.15
#> 16 M3_I4    2  3PLM   1.3   1.2  0.15
#> 17 M3_I5    2  3PLM   1.3   1.5  0.15
#> 18 M3_I6    2  3PLM   1.3   1.8  0.15
#> 19 M4_I1    2  3PLM   1.0  -2.5  0.15
#> 20 M4_I2    2  3PLM   1.0  -2.2  0.15
#> 21 M4_I3    2  3PLM   1.0  -1.9  0.15
#> 22 M4_I4    2  3PLM   1.0  -1.6  0.15
#> 23 M4_I5    2  3PLM   1.0  -1.3  0.15
#> 24 M4_I6    2  3PLM   1.0  -1.0  0.15
#> 25 M5_I1    2  3PLM   1.4   0.8  0.15
#> 26 M5_I2    2  3PLM   1.4   1.1  0.15
#> 27 M5_I3    2  3PLM   1.4   1.4  0.15
#> 28 M5_I4    2  3PLM   1.4   1.7  0.15
#> 29 M5_I5    2  3PLM   1.4   2.0  0.15
#> 30 M5_I6    2  3PLM   1.4   2.3  0.15
```

### Step 2: Build the Module Matrix

The `module` matrix has the same number of rows as `item_bank_122` and
one column per module. Each row has exactly one 1 — in the column for
the module that item belongs to.

``` r

n_mods  <- 5
n_items <- nrow(item_bank_122)   # 30

module_122 <- matrix(0L, nrow = n_items, ncol = n_mods)
colnames(module_122) <- paste0("M", 1:n_mods)

# Items 1-6 → M1, items 7-12 → M2, ..., items 25-30 → M5
for (m in 1:n_mods) {
  idx <- ((m - 1) * n_per_mod + 1):(m * n_per_mod)
  module_122[idx, m] <- 1L
}

# Verify: each row sums to 1, each column sums to n_per_mod
cat("Row sums (all should be 1):\n"); print(table(rowSums(module_122)))
#> Row sums (all should be 1):
#> 
#>  1 
#> 30
cat("Column sums (all should be", n_per_mod, "):\n"); print(colSums(module_122))
#> Column sums (all should be 6 ):
#> M1 M2 M3 M4 M5 
#>  6  6  6  6  6
```

### Step 3: Build the Route Map

``` r

route_map_122 <- matrix(0L, nrow = n_mods, ncol = n_mods)
rownames(route_map_122) <- colnames(route_map_122) <- paste0("M", 1:n_mods)

# Stage 1 → Stage 2
route_map_122[1, 2] <- 1L   # M1 → M2 (low ability)
route_map_122[1, 3] <- 1L   # M1 → M3 (high ability)

# Stage 2 → Stage 3 (both M2 and M3 can route to either M4 or M5)
route_map_122[2, 4] <- 1L   # M2 → M4
route_map_122[2, 5] <- 1L   # M2 → M5
route_map_122[3, 4] <- 1L   # M3 → M4
route_map_122[3, 5] <- 1L   # M3 → M5

print(route_map_122)
#>    M1 M2 M3 M4 M5
#> M1  0  1  1  0  0
#> M2  0  0  0  1  1
#> M3  0  0  0  1  1
#> M4  0  0  0  0  0
#> M5  0  0  0  0  0
```

### Step 4: Define Cut Scores

``` r

# 1 cut point for each of the 2 routing transitions (2-way branching)
cut_score_122 <- list(
  c(0),    # Stage 1 → Stage 2: θ̂ < 0 → M2, θ̂ ≥ 0 → M3
  c(0)     # Stage 2 → Stage 3: θ̂ < 0 → M4, θ̂ ≥ 0 → M5
)
```

### Step 5: Evaluate the Custom Panel

``` r

# Evaluation grid: -3 to 3 in steps of 0.5 (coarser for speed)
theta_grid_122 <- seq(-3, 3, 0.5)

eval_122 <- reval_mst(
  x         = item_bank_122,
  D         = 1.702,
  route_map = route_map_122,
  module    = module_122,
  cut_score = cut_score_122,
  theta     = theta_grid_122,
  range.tcc = c(-5, 5)
)

# Evaluation table
print(eval_122$eval.tb)
#>    theta          mu    sigma2         bias      csem
#> 1   -3.0 -3.19623125 0.8592354 -0.196231249 0.9269495
#> 2   -2.5 -2.72861584 0.7871656 -0.228615843 0.8872235
#> 3   -2.0 -2.11902994 0.5155952 -0.119029937 0.7180496
#> 4   -1.5 -1.52846765 0.2785520 -0.028467654 0.5277803
#> 5   -1.0 -0.99289282 0.2024505  0.007107183 0.4499450
#> 6   -0.5 -0.47947559 0.1729555  0.020524412 0.4158792
#> 7    0.0 -0.01380778 0.1486202 -0.013807784 0.3855129
#> 8    0.5  0.46107384 0.1534703 -0.038926159 0.3917529
#> 9    1.0  0.99556805 0.1267071 -0.004431952 0.3559594
#> 10   1.5  1.52787128 0.1395430  0.027871283 0.3735546
#> 11   2.0  2.15212965 0.5041920  0.152129651 0.7100648
#> 12   2.5  3.12555584 1.5165325  0.625555842 1.2314757
#> 13   3.0  4.12958178 1.4429566  1.129581775 1.2012313
```

``` r

eval_tb_122 <- eval_122$eval.tb

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

plot(
  eval_tb_122$theta, eval_tb_122$bias,
  type = "b", pch = 16, col = "steelblue", lwd = 2,
  xlab = expression(theta), ylab = "Conditional Bias",
  main = "1-2-2 Panel: Conditional Bias",
  ylim = range(c(-0.4, 0.4, eval_tb_122$bias))
)
abline(h = 0, col = "red", lty = 2, lwd = 1.5)
grid()

plot(
  eval_tb_122$theta, eval_tb_122$csem,
  type = "b", pch = 16, col = "darkorange", lwd = 2,
  xlab = expression(theta), ylab = "CSEM",
  main = "1-2-2 Panel: CSEM",
  ylim = c(0, max(eval_tb_122$csem) * 1.2)
)
grid()
```

![Conditional bias and CSEM for the custom 1-2-2 MST
panel.](mst-panel-evaluation_files/figure-html/plot-custom-1.png)

Conditional bias and CSEM for the custom 1-2-2 MST panel.

``` r


par(mfrow = c(1, 1))
```

### Step 6: Inspect the Panel Pathways

``` r

# Confirmed pathways through the 1-2-2 panel
eval_122$panel.info$pathway
#>        stage.1 stage.2 stage.3
#> path.1       1       2       4
#> path.2       1       2       5
#> path.3       1       3       4
#> path.4       1       3       5
```

The 1-2-2 panel has 4 pathways — considerably fewer than the 1-3-3
panel’s 7 pathways. The simpler branching structure is appropriate for
smaller-scale tests or when fewer stage-3 difficulty levels are needed.

------------------------------------------------------------------------

## Example 5: Comparing Two Cut Score Configurations

Cut score placement directly affects measurement bias. Placing cuts too
close to the centre can cause mid-ability examinees to be routed
sub-optimally. Here we compare the original `c(-0.5, 0.5)` cut scores
with a wider configuration `c(-1.0, 1.0)` for the 1-3-3 panel.

``` r

# Alternative cut scores: wider routing bands
cut_score_wide <- list(c(-1.0, 1.0), c(-1.0, 1.0))

eval_wide <- reval_mst(
  x         = x,
  D         = 1.702,
  route_map = route_map,
  module    = module,
  cut_score = cut_score_wide,
  theta     = theta,
  range.tcc = c(-5, 5)
)
```

``` r

tb_orig <- eval_result$eval.tb
tb_wide <- eval_wide$eval.tb

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Bias comparison
ylim_bias <- range(c(tb_orig$bias, tb_wide$bias, -0.4, 0.4))
plot(
  tb_orig$theta, tb_orig$bias,
  type = "b", pch = 16, col = "steelblue", lwd = 2,
  xlab = expression(theta), ylab = "Conditional Bias",
  main = "Bias Comparison", ylim = ylim_bias
)
lines(tb_wide$theta, tb_wide$bias,
      type = "b", pch = 17, col = "darkorange", lwd = 2, lty = 2)
abline(h = 0, col = "grey50", lty = 3)
legend("topright",
       legend = c("Original cuts (±0.5)", "Wider cuts (±1.0)"),
       col    = c("steelblue", "darkorange"),
       pch    = c(16, 17), lty = c(1, 2), lwd = 2, cex = 0.85)
grid()

# CSEM comparison
ylim_csem <- c(0, max(c(tb_orig$csem, tb_wide$csem)) * 1.2)
plot(
  tb_orig$theta, tb_orig$csem,
  type = "b", pch = 16, col = "steelblue", lwd = 2,
  xlab = expression(theta), ylab = "CSEM",
  main = "CSEM Comparison", ylim = ylim_csem
)
lines(tb_wide$theta, tb_wide$csem,
      type = "b", pch = 17, col = "darkorange", lwd = 2, lty = 2)
legend("topright",
       legend = c("Original cuts (±0.5)", "Wider cuts (±1.0)"),
       col    = c("steelblue", "darkorange"),
       pch    = c(16, 17), lty = c(1, 2), lwd = 2, cex = 0.85)
grid()
```

![CSEM comparison: original vs. wider cut scores in the 1-3-3
panel.](mst-panel-evaluation_files/figure-html/compare-cuts-plot-1.png)

CSEM comparison: original vs. wider cut scores in the 1-3-3 panel.

``` r


par(mfrow = c(1, 1))
```

Wide cuts concentrate test takers in the middle module of Stage 2 and 3
for the most ability levels, which tends to reduce bias near the centre
but may increase CSEM at the ability extremes.

------------------------------------------------------------------------

## Example 6: Validating `reval_mst()` Against a `run_mst()` Monte Carlo Simulation

The “Why Simulate MST Panel Performance?” section introduced
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md) and
noted that turning a single simulation run into a full panel evaluation
means looping it over a $`\theta`$ grid with many replications per point
— the traditional Monte Carlo approach. Lim et al. (2021) validated
their recursion-based analytical method against exactly this kind of
simulation, and the closeness of the agreement they found depended on
the scoring method used. When the simulation used the same scoring the
recursion assumes throughout — an equated-number-correct (ENC) score, a
deterministic score-to-$`\theta`$ transformation conceptually analogous
to
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)’s
internal inverse-TCC scoring — the simulated and analytical results
matched closely at every $`\theta`$, and converged further as
replications increased. When the simulation instead used a combination
more common in practice — EAP for routing, maximum likelihood (MLE) for
the final score — the two methods still agreed closely in the central
ability range, but diverged toward the extremes, where the item pool
carries little information and MLE behaves quite differently from the
inverse-TCC-based recursion. That region of close central agreement
widens as test length increases.

We reproduce that second condition directly here for the `simMST` panel,
using the same scoring combination Lim et al. (2021) used in their own
validation: `route_score = "EAP"` for routing at each interim stage, and
`final_score = "ML"` for the reported final score. Neither method
matches
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)’s
internal inverse-TCC logic, so the two results are not expected to
coincide everywhere; the comparison below should reproduce the same
central-agreement-with-tail-divergence pattern Lim et al. (2021)
reported.

Because the exact values stored in `simMST$theta` are not needed for
this comparison, we define our own explicit evaluation grid so that the
Monte Carlo simulation and the analytical computation are guaranteed to
line up point-for-point:

``` r

# An explicit, independent evaluation grid
theta_grid <- seq(-5, 5, 0.5)
n_reps     <- 300   # simulated examinees per grid point

# Repeat each grid point n_reps times into one long vector of true abilities
set.seed(2026)
true_theta <- rep(theta_grid, each = n_reps)
```

**Step 1 — Monte Carlo simulation with
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md).** All
`n_reps` replications at every grid point are simulated in a single call
by passing the full `true_theta` vector. Routing uses the panel’s fixed
cut scores (matching the routing rule
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
assumes), with `EAP` routing scores and `ML` final scores — the scoring
combination Lim et al. (2021) used in their own Monte Carlo validation:

``` r

sim_mc <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = true_theta,
  D            = 1.702,
  route_method = NULL,
  cut_score    = cut_score,
  route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
  final_score  = list(method = "ML", range = c(-4, 4)),
  se           = FALSE,
  verbose      = FALSE
)

# Empirical bias and RMSE at each grid point
mc_tb <- data.frame(theta = true_theta, est.theta = sim_mc$est.theta)
diff_mc  <- mc_tb$est.theta - mc_tb$theta
bias_mc  <- tapply(diff_mc, mc_tb$theta, mean)
rmse_mc  <- tapply(diff_mc, mc_tb$theta, function(d) sqrt(mean(d^2)))

mc_result <- data.frame(
  theta = as.numeric(names(bias_mc)),
  bias  = as.numeric(bias_mc),
  rmse  = as.numeric(rmse_mc)
)
mc_result <- mc_result[order(mc_result$theta), ]
```

**Step 2 — Analytical computation with
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)**
on the *same* `theta_grid`, so the comparison is point-for-point:

``` r

eval_grid <- reval_mst(
  x         = x,
  D         = 1.702,
  route_map = route_map,
  module    = module,
  cut_score = cut_score,
  theta     = theta_grid,
  range.tcc = c(-5, 5)
)

tb_grid      <- eval_grid$eval.tb
tb_grid$rmse <- sqrt(tb_grid$sigma2 + tb_grid$bias^2)
```

**Step 3 — Overlay the two results:**

``` r

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Bias: Monte Carlo (points) vs. analytical (line)
plot(
  tb_grid$theta, tb_grid$bias,
  type = "l", lwd = 2, col = "steelblue",
  xlab = expression(theta), ylab = "Bias", main = "Bias: Simulation vs. Analytical",
  ylim = range(c(tb_grid$bias, mc_result$bias, -0.3, 0.3))
)
points(mc_result$theta, mc_result$bias, pch = 16, col = "darkorange")
abline(h = 0, col = "grey50", lty = 3)
legend("topright",
       legend = c("reval_mst() (analytical)", "run_mst() (Monte Carlo)"),
       col = c("steelblue", "darkorange"), lty = c(1, NA), pch = c(NA, 16),
       lwd = c(2, NA), cex = 0.8)
grid()

# RMSE: Monte Carlo (points) vs. analytical (line)
plot(
  tb_grid$theta, tb_grid$rmse,
  type = "l", lwd = 2, col = "steelblue",
  xlab = expression(theta), ylab = "RMSE", main = "RMSE: Simulation vs. Analytical",
  ylim = c(0, max(c(tb_grid$rmse, mc_result$rmse)) * 1.2)
)
points(mc_result$theta, mc_result$rmse, pch = 16, col = "darkorange")
legend("topright",
       legend = c("reval_mst() (analytical)", "run_mst() (Monte Carlo)"),
       col = c("steelblue", "darkorange"), lty = c(1, NA), pch = c(NA, 16),
       lwd = c(2, NA), cex = 0.8)
grid()
```

![Empirical bias and RMSE from a run_mst() Monte Carlo simulation
(points) overlaid on the analytical reval_mst() curve
(line).](mst-panel-evaluation_files/figure-html/mc-vs-analytical-plot-1.png)

Empirical bias and RMSE from a run_mst() Monte Carlo simulation (points)
overlaid on the analytical reval_mst() curve (line).

``` r


par(mfrow = c(1, 1))
```

In the central ability range — which Lim et al. (2021) reported as
roughly $`-1.5 \le \theta \le 1.5`$ for this scoring combination,
widening for longer tests — the Monte Carlo points should track the
analytical curve closely, with the small remaining deviations
attributable to sampling error at `n_reps = 300` replications per grid
point (these shrink further as `n_reps` increases). Toward the extremes,
the two are expected to diverge: both `EAP` routing and `ML` final
scoring depart from
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)’s
internal inverse-TCC logic, and that departure grows where the item pool
carries little information and EAP/MLE estimates behave quite
differently from the inverse-TCC-based recursion. This mirrors the
pattern Lim et al. (2021) reported for this identical
EAP-routing/MLE-final condition. Matching the recursion’s own scoring
exactly — `route_score = "INV.TCC"` and `final_score = "INV.TCC"`
instead of `"EAP"`/`"ML"` — would close most of this gap.
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md) is
therefore best read as predicting the performance of a panel routed and
scored via inverse-TCC; whenever a program uses a different scoring
approach operationally, a
[`run_mst()`](https://hwangQ.github.io/irtQ/reference/run_mst.md)
simulation like this one remains the more direct check, particularly for
ability levels far from the panel’s region of peak information.

**Step 4 — Confirm the gap closes under matched scoring.** As a final
check, re-run the same simulation with `route_score = "INV.TCC"` and
`final_score = "INV.TCC"`, matching
[`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)’s
internal logic exactly:

``` r

sim_mc_matched <- run_mst(
  x            = x,
  route_map    = route_map,
  module       = module,
  theta        = true_theta,
  D            = 1.702,
  route_method = NULL,
  cut_score    = cut_score,
  route_score  = list(method = "INV.TCC", range.tcc = c(-5, 5)),
  final_score  = list(method = "INV.TCC", range.tcc = c(-5, 5)),
  se           = FALSE,
  verbose      = FALSE
)

# Empirical bias and RMSE at each grid point, matched-scoring condition
diff_matched <- sim_mc_matched$est.theta - true_theta
bias_matched <- tapply(diff_matched, true_theta, mean)
rmse_matched <- tapply(diff_matched, true_theta, function(d) sqrt(mean(d^2)))

mc_result_matched <- data.frame(
  theta = as.numeric(names(bias_matched)),
  bias  = as.numeric(bias_matched),
  rmse  = as.numeric(rmse_matched)
)
mc_result_matched <- mc_result_matched[order(mc_result_matched$theta), ]
```

Comparing all three at the tail grid points, where the `EAP`/`ML`
condition diverged most from the analytical curve:

``` r

compare_tb <- data.frame(
  theta              = theta_grid,
  analytical_bias    = round(tb_grid$bias, 3),
  EAP_ML_bias        = round(mc_result$bias, 3),
  INVTCC_INVTCC_bias = round(mc_result_matched$bias, 3),
  analytical_rmse    = round(tb_grid$rmse, 3),
  EAP_ML_rmse        = round(mc_result$rmse, 3),
  INVTCC_INVTCC_rmse = round(mc_result_matched$rmse, 3)
)
compare_tb[compare_tb$theta <= -3 | compare_tb$theta >= 3, ]
#>    theta analytical_bias EAP_ML_bias INVTCC_INVTCC_bias analytical_rmse
#> 1   -5.0           1.658       1.549              1.670           1.887
#> 2   -4.5           1.165       1.057              1.260           1.473
#> 3   -4.0           0.685       0.479              0.732           1.132
#> 4   -3.5           0.240       0.130              0.229           0.933
#> 5   -3.0          -0.109      -0.179             -0.109           0.900
#> 17   3.0           1.304       0.585              1.225           1.753
#> 18   3.5           1.248       0.372              1.265           1.469
#> 19   4.0           0.915      -0.054              0.887           1.026
#> 20   4.5           0.472      -0.527              0.491           0.544
#> 21   5.0          -0.009      -1.005             -0.009           0.157
#>    EAP_ML_rmse INVTCC_INVTCC_rmse
#> 1        1.706              1.900
#> 2        1.255              1.565
#> 3        0.822              1.151
#> 4        0.747              0.919
#> 5        0.787              0.890
#> 17       0.924              1.726
#> 18       0.571              1.469
#> 19       0.298              1.034
#> 20       0.566              0.514
#> 21       1.010              0.151
```

The same comparison, overlaid across the full ability range:

``` r

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Bias: matched-scoring Monte Carlo (points) vs. analytical (line)
plot(
  tb_grid$theta, tb_grid$bias,
  type = "l", lwd = 2, col = "steelblue",
  xlab = expression(theta), ylab = "Bias",
  main = "Bias: Simulation vs. Analytical",
  ylim = range(c(tb_grid$bias, mc_result_matched$bias, -0.3, 0.3))
)
points(mc_result_matched$theta, mc_result_matched$bias, pch = 16, col = "forestgreen")
abline(h = 0, col = "grey50", lty = 3)
legend("topright",
       legend = c("reval_mst() (analytical)", "run_mst() (INV.TCC/INV.TCC)"),
       col = c("steelblue", "forestgreen"), lty = c(1, NA), pch = c(NA, 16),
       lwd = c(2, NA), cex = 0.8)
grid()

# RMSE: matched-scoring Monte Carlo (points) vs. analytical (line)
plot(
  tb_grid$theta, tb_grid$rmse,
  type = "l", lwd = 2, col = "steelblue",
  xlab = expression(theta), ylab = "RMSE",
  main = "RMSE: Simulation vs. Analytical",
  ylim = c(0, max(c(tb_grid$rmse, mc_result_matched$rmse)) * 1.2)
)
points(mc_result_matched$theta, mc_result_matched$rmse, pch = 16, col = "forestgreen")
legend("topright",
       legend = c("reval_mst() (analytical)", "run_mst() (INV.TCC/INV.TCC)"),
       col = c("steelblue", "forestgreen"), lty = c(1, NA), pch = c(NA, 16),
       lwd = c(2, NA), cex = 0.8)
grid()
```

![Empirical bias and RMSE from the INV.TCC/INV.TCC run_mst() Monte Carlo
simulation (points) overlaid on the analytical reval_mst() curve
(line).](mst-panel-evaluation_files/figure-html/mc-vs-analytical-plot-matched-1.png)

Empirical bias and RMSE from the INV.TCC/INV.TCC run_mst() Monte Carlo
simulation (points) overlaid on the analytical reval_mst() curve (line).

``` r


par(mfrow = c(1, 1))
```

Both the table and the plot confirm the same conclusion: the
`INV.TCC`/`INV.TCC` condition tracks the analytical bias and RMSE
closely across the *entire* ability range, including the extremes where
the `EAP`/`ML` condition diverged.

------------------------------------------------------------------------

## Function Reference

### `reval_mst()` Arguments

| Argument | Type | Default | Description |
|----|----|----|----|
| `x` | data.frame | — | Item bank metadata (irtQ format) |
| `D` | numeric | 1 | Scaling constant (use 1.702 for normal-ogive approximation) |
| `route_map` | matrix | — | Binary square matrix of module transitions |
| `module` | matrix | — | Binary matrix mapping items to modules |
| `cut_score` | list | — | List of routing cut score vectors (one per stage transition) |
| `theta` | numeric | `seq(-5, 5, 1)` | Ability grid for evaluation |
| `intpol` | logical | TRUE | Linear interpolation for out-of-range TCC scores |
| `range.tcc` | numeric(2) | `c(-7, 7)` | Ability range for inverse TCC scoring |
| `tol` | numeric | 1e-4 | Convergence tolerance for bisection in inverse TCC |

### `reval_mst()` Return Value

| Component | Description |
|----|----|
| `panel.info` | Panel structure: `$config`, `$pathway`, `$n.module`, `$n.stage` |
| `item.by.mod` | List of item metadata data frames, one per module |
| `item.by.path` | List of cumulative item metadata per stage and pathway |
| `eq.theta` | Inverse-TCC $`\theta`$ estimates for each possible score, by stage and pathway |
| `cdist.by.mod` | Conditional score distributions per module, indexed by $`\theta`$ |
| `jdist.by.path` | Joint conditional score distributions, indexed by $`\theta`$ and stage |
| `eval.tb` | Evaluation table: `theta`, `mu`, `sigma2`, `bias`, `csem` |

### `run_mst()` Arguments

| Argument | Type | Default | Description |
|----|----|----|----|
| `x` | data.frame | — | Item bank metadata (irtQ format) |
| `route_map` | matrix | — | Binary square matrix of module transitions |
| `module` | matrix | — | Binary matrix mapping items to modules |
| `theta` | numeric | — | True ability for each simulated examinee (or supply `response` instead) |
| `response` | matrix | `NULL` | Pre-generated response matrix, used instead of simulating from `theta` |
| `D` | numeric | 1 | Scaling constant (use 1.702 for normal-ogive approximation) |
| `ini_mod` | integer | `NULL` | Fixed Stage-1 module for all examinees; `NULL` assigns each independently at random |
| `route_method` | character | `"bmat"` | `"bmat"`, `"mfi"`, or `NULL` for fixed cut-score routing |
| `cut_score` | list | `NULL` | List of routing cut score vectors; used only when `route_method = NULL` |
| `route_score` | list | — | Scoring method used at intermediate stages (`ML`, `WL`, `MAP`, `EAP`, `EAP.SUM`, `INV.TCC`, …) |
| `final_score` | list | — | Scoring method used for the final reported score |
| `se` | logical | `TRUE` | Whether to compute standard errors of the final score |
| `missing` | scalar | `NA` | Value representing a not-administered response in `response` |
| `verbose` | logical | `TRUE` | Whether to print simulation progress messages to the console |
| `return_full_resp` | logical | `FALSE` | If `TRUE`, also returns the full *N* x *J* item-level response matrix as `full.resp` |

### `run_mst()` Return Value

| Component | Description |
|----|----|
| `est.theta` | Final ability estimate for each simulated examinee |
| `se.theta` | Standard error of the final ability estimate |
| `theta.route` | Intermediate ability estimate used for routing at each stage |
| `path` | Module pathway taken by each examinee |
| `true.theta` | True ability supplied as input (`NA` when `response` was supplied instead) |
| `panel` | Panel structure information (from [`panel_info()`](https://hwangQ.github.io/irtQ/reference/panel_info.md)) |
| `full.resp` | Full item-level response matrix, if `return_full_resp = TRUE` |

------------------------------------------------------------------------

## References

Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
approach to evaluate the performance of MST. *Journal of Educational
Measurement*, *58*(2), 154–178. <https://doi.org/10.1111/jedm.12276>

Lord, F. M., & Wingersky, M. S. (1984). Comparison of IRT true-score and
equipercentile observed-score equatings. *Applied Psychological
Measurement*, *8*(4), 453–461.
<https://doi.org/10.1177/014662168400800409>

Luecht, R. M., & Nungester, R. J. (1998). Some practical examples of
computer-adaptive sequential testing. *Journal of Educational
Measurement*, *35*(3), 229–249.
<https://doi.org/10.1111/j.1745-3984.1998.tb00537.x>

Zenisky, A., Hambleton, R. K., & Luecht, R. M. (2009). Multistage
testing: Issues, designs, and research. In W. J. van der Linden & C. A.
W. Glas (Eds.), *Elements of adaptive testing* (pp. 355–372). Springer.
<https://doi.org/10.1007/978-0-387-85461-8_18>
