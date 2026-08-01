# Option/Category Response Distribution and Distractor Analysis

Computes, for each item, the percentage of examinees choosing each
response option (or, for already-scored polytomous data, each score
category), together with each option's/category's point-biserial
correlation with the total score. When an answer key is supplied, the
function operates in "selected-response" mode - treating `data` as raw,
unscored option responses (as consumed by
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)) -
and additionally reports omission/double-marking rates and flags
distractors (incorrect options) whose correlation with the total score
is positive, a common indicator of a miskeyed or poorly written item.
When no key is supplied, the function operates in "scored-category"
mode - treating `data` as an already-scored response matrix (0, 1, 2,
..., cats - 1) - and simply reports each score category's response
proportion and category-total correlation, with no notion of a "correct"
option to flag distractors against.

## Usage

``` r
ctt_distr(
  data,
  item.id = NULL,
  key = NULL,
  opt = NULL,
  cats = NULL,
  total = NULL,
  correct = FALSE,
  missing = NA,
  crit.distractor = 0
)
```

## Arguments

- data:

  In selected-response mode (`key` supplied), a data frame or matrix of
  raw item responses in the same format as the `data` argument of
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
  (a selected option number, a missing-response indicator, or a
  comma-separated string such as `"1,5"` for a double-marked response).
  In scored-category mode (`key = NULL`), a data frame or matrix of
  already-scored item responses with scores 0 to `cats[j] - 1` for each
  item j. Either way, examinees are rows and items are columns.

- item.id:

  A character vector of item identifiers, in the same order as the
  columns of `data`. If `NULL` (default), item IDs are generated
  automatically as `paste0("V", 1:ncol(data))`, following the convention
  used elsewhere in irtQ (e.g.,
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md),
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md)).
  Note that this generated ID does not fall back to `colnames(data)`;
  pass `item.id` explicitly to label items using the column names of
  `data` or any other identifier scheme.

- key:

  `NULL` (default) for scored-category mode, or - to use
  selected-response mode - either a vector of correct options in item
  order, or a data frame with columns `item` and `key`, exactly as
  accepted by
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
  As in
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md),
  each item's option-coding scheme (numeric, Latin-letter, or general
  label such as a Korean syllable label) is inferred independently from
  its own key value; see **Details**.

- opt:

  A vector giving the full set of possible response option values in
  selected-response mode (e.g., `1:5` for a five-option item, or
  `c("A", "B", "C", "D", "E")` for a letter-coded item), applied to
  every item. Only used when `key` is supplied; it is silently ignored
  in scored-category mode (`key = NULL`), which instead derives its
  category set from `cats`. If `NULL` (default), the option set is
  inferred as the sorted set of all single (non-blank,
  non-double-marked) response values observed anywhere in `data`, using
  each item's own option-coding scheme to decide what counts as a single
  response (see **Details**). Because `opt` is shared across every item,
  `ctt_distr()` assumes a test form uses one consistent option-coding
  scheme throughout; see **Details** for what happens if items actually
  differ.

- cats:

  A numeric vector giving the number of score categories per item,
  following the same convention as
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md). Only used
  in scored-category mode (`key = NULL`); silently ignored when `key` is
  supplied. If `NULL` (default), inferred per item as the observed
  maximum score plus one.

- total:

  An optional numeric vector of total scores, one per examinee, aligned
  with the rows of `data`. If `NULL` (default), the total score is
  computed internally: via
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
  in selected-response mode, or as `rowSums(data)` in scored-category
  mode (see **Details** for how `total` interacts with
  `missing`/listwise deletion in scored-category mode).

- correct:

  Logical. Both the raw (uncorrected) item-total correlation - where an
  item (or, in `ctt_distr()`, a response option/category) is correlated
  with the total score that includes its own contribution - and the
  corrected item-total correlation - excluding its own contribution -
  are always computed and reported as separate columns (see
  **Details**/**Value**). This argument only controls which of the two
  is used for flagging: if `TRUE`, flagging (in
  [`ctt_item()`](https://hwangQ.github.io/irtQ/reference/ctt_item.md)'s
  `flag` column,
  [`ctt_alpha()`](https://hwangQ.github.io/irtQ/reference/ctt_alpha.md)'s
  summary discrimination figure, and, in `ctt_distr()`, the
  `crit.distractor` check) is based on the corrected value; if `FALSE`
  (default), on the raw value.

- missing:

  A value indicating missing responses, analogous to the `missing`
  argument in
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) and
  [`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
  Its effect differs by mode (see **Details**): in selected-response
  mode, a cell equal to `missing` is recoded to `NA` and then treated as
  a blank (omitted) response, contributing to that item's omission rate,
  without removing the examinee from the analysis. In scored-category
  mode, a cell equal to `missing` is recoded to `NA` and the examinee's
  row is then excluded listwise, as in
  [`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md). Default is
  `NA`.

- crit.distractor:

  A single numeric value; in selected-response mode, a distractor
  (non-key option) is flagged when its point-biserial correlation with
  the total score exceeds this value. Default is `0`.

## Value

A list with these elements:

- distr:

  A data frame in long format with one row per item x option
  (selected-response mode) or item x category (scored-category mode),
  containing `item`, `option` (the option number, option label, or score
  category; numeric for a numeric-coded or scored-category test form,
  character for a letter- or general-label-coded test form), `is_key`
  (logical; whether this option is the item's correct answer, or `NA` in
  scored-category mode), `freq` (the raw number of examinees choosing
  this option/category), `pct` (the same count expressed as a percentage
  of examinees), `pb_raw` (raw, uncorrected point-biserial correlation
  with the total score), `pb_corrected` (corrected, item-excluded
  point-biserial correlation with the total score), and `flag` (a
  description of `crit.distractor` violations, based on whichever of
  `pb_raw`/`pb_corrected` is selected by `correct`, in selected-response
  mode, or `""` otherwise).

- omit:

  In selected-response mode, a data frame with one row per item giving
  the blank/omission and double-marked response percentages
  (`pct_blank`, `pct_double`). `NULL` in scored-category mode.

## Details

In selected-response mode, each item's option-coding scheme is inferred
from its own `key` value, using exactly the same three-way rule as
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md):
a key value that parses as a number is compared numerically; a key value
consisting of one or more Latin letters is compared as letters,
case-insensitively, and a response must itself consist only of Latin
letters to count as a single valid option; any other key value (e.g., a
Korean syllable label, a Roman numeral written in a non-Latin numeral
script, or a circled-number symbol) is treated as a general option
label, under which any non-blank, non-double-marked response is accepted
as a single valid option and compared to the key (and to `opt`) with a
case-insensitive exact string match. This ensures `ctt_distr()`
interprets a given `data`/`key` pair identically to
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
Because `opt` (the option universe used for every item's frequency
table) is a single vector shared across all items, `ctt_distr()` is
intended for a test form that uses one consistent option-coding scheme
throughout (e.g., every item numbered 1-5, or every item lettered A-E);
mixing schemes across items within a single call is not well-supported,
since the shared `opt` cannot represent two different option universes
at once.

The two modes handle missing responses differently because they
represent different kinds of data. In selected-response mode, an omitted
response is itself a meaningful, common outcome for a selected-response
item (the examinee saw the item and chose not to answer it), so
omissions are tabulated as their own category rather than causing the
examinee to be dropped, exactly as in
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md).
In scored-category mode, `data` is assumed to already be a scored
response matrix (as used by
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md)), where a
missing cell instead means the item was not administered to, or not
observed for, that examinee, so the examinee's row is excluded listwise
for that analysis. When `total` is supplied together with
scored-category mode's listwise deletion, `total` is assumed to align
with the *original* rows of `data`; it is subset using the same
completeness mask before use, so its length need only match the original
`nrow(data)`.

Point-biserial correlations are computed as the Pearson correlation
between a 0/1 indicator (whether the examinee chose that particular
option, or scored in that particular category) and the total score. As
with the item-level discrimination statistics computed by
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md), both the raw
(uncorrected, item/option-included) and corrected (item-excluded)
versions are always computed and returned as separate columns
(`pb_raw`/`pb_corrected`); `correct` only selects which of the two feeds
the `crit.distractor` flagging criterion in selected-response mode. A
point-biserial correlation is undefined, and reported as `NA`, whenever
an option/category was chosen by zero or all examinees (a constant
indicator).

## Note

Selected-response mode (`key` supplied) is restricted to **dichotomous,
single-key items only**, exactly like
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md),
which it calls internally to obtain the 0/1 scored responses and total
score: a selected-response item has exactly one correct option and
therefore cannot yield a partial-credit score by construction, so a
mixed dichotomous/polytomous test form cannot be analyzed in this mode.
Scored-category mode (`key = NULL`) has no such restriction: because
`cats` is always a per-item vector, a single call fully supports a
genuinely mixed dichotomous/polytomous test form (see **Examples**). If
raw selected-response data and already-scored polytomous data need to be
analyzed together, score the selected-response items first (via
[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md)
or this function's own selected-response mode), then combine the result
with the polytomous scores and call `ctt_distr()` again in
scored-category mode on the combined matrix.

## See also

[`score_resp()`](https://hwangQ.github.io/irtQ/reference/score_resp.md),
[`ctt()`](https://hwangQ.github.io/irtQ/reference/ctt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# Selected-response mode: raw option responses + an answer key
raw <- data.frame(
  V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
  V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
)
key <- c(1, 4)
out <- ctt_distr(data = raw, key = key)
out$distr
#>   item option is_key freq  pct pb_raw pb_corrected flag
#> 1   V1      1   TRUE    4 50.0  0.626       -0.258     
#> 2   V1      2  FALSE    2 25.0 -0.602       -0.149     
#> 3   V1      3  FALSE    1 12.5 -0.079        0.293     
#> 4   V1      4  FALSE    1 12.5 -0.079        0.293     
#> 5   V2      1  FALSE    1 12.5 -0.709       -0.378     
#> 6   V2      2  FALSE    1 12.5 -0.079        0.378     
#> 7   V2      3  FALSE    0  0.0     NA           NA     
#> 8   V2      4   TRUE    5 62.5  0.592       -0.258     
out$omit
#>   item pct_blank pct_double
#> 1   V1         0        0.0
#> 2   V2         0       12.5

# Scored-category mode: an already-scored polytomous matrix, no key
set.seed(1)
scored <- data.frame(matrix(sample(0:2, 8 * 3, replace = TRUE), nrow = 8))
out2 <- ctt_distr(data = scored)
out2$distr
#>   item option is_key freq  pct pb_raw pb_corrected flag
#> 1   V1      0     NA    3 37.5  0.026        0.613     
#> 2   V1      1     NA    2 25.0 -0.533       -0.493     
#> 3   V1      2     NA    3 37.5  0.450       -0.172     
#> 4   V2      0     NA    3 37.5 -0.609       -0.098     
#> 5   V2      1     NA    3 37.5 -0.185       -0.488     
#> 6   V2      2     NA    2 25.0  0.889        0.655     
#> 7   V3      0     NA    4 50.0 -0.513        0.107     
#> 8   V3      1     NA    2 25.0  0.415        0.309     
#> 9   V3      2     NA    2 25.0  0.178       -0.433     

# Scored-category mode also handles a genuinely mixed dichotomous/
# polytomous test form in a single call, since `cats` is always a
# per-item vector (see Note): item 3 below contributes 3 category rows
# while items 1-2 contribute 2 each
mixed <- data.frame(
  I1 = c(1, 0, 1, 1, 0, 1, 0, 1),   # dichotomous item (cats = 2)
  I2 = c(0, 1, 0, 1, 1, 0, 0, 1),   # dichotomous item (cats = 2)
  I3 = c(2, 0, 1, 2, 1, 0, 2, 1)    # polytomous item (cats = 3)
)
out2b <- ctt_distr(data = mixed, cats = c(2, 2, 3))
out2b$distr
#>   item option is_key freq  pct pb_raw pb_corrected flag
#> 1   V1      0     NA    3 37.5 -0.467        0.038     
#> 2   V1      1     NA    5 62.5  0.467       -0.038     
#> 3   V2      0     NA    4 50.0 -0.258        0.258     
#> 4   V2      1     NA    4 50.0  0.258       -0.258     
#> 5   V3      0     NA    2 25.0 -0.745       -0.120     
#> 6   V3      1     NA    3 37.5  0.067        0.269     
#> 7   V3      2     NA    3 37.5  0.600       -0.162     

# Selected-response mode with letter-coded (A-E) options: each item's
# coding scheme is inferred from its own key value, exactly as in
# irtQ::score_resp()
raw3 <- data.frame(
  V1 = c("A", "B", "A", "A", "C", "B", "A", "D"),
  V2 = c("D", "D", "B", "B,D", "D", "A", "D", "D")
)
key3 <- c("A", "D")
out3 <- ctt_distr(data = raw3, key = key3)
out3$distr
#>   item option is_key freq  pct pb_raw pb_corrected flag
#> 1   V1      A   TRUE    4 50.0  0.626       -0.258     
#> 2   V1      B  FALSE    2 25.0 -0.602       -0.149     
#> 3   V1      C  FALSE    1 12.5 -0.079        0.293     
#> 4   V1      D  FALSE    1 12.5 -0.079        0.293     
#> 5   V2      A  FALSE    1 12.5 -0.709       -0.378     
#> 6   V2      B  FALSE    1 12.5 -0.079        0.378     
#> 7   V2      C  FALSE    0  0.0     NA           NA     
#> 8   V2      D   TRUE    5 62.5  0.592       -0.258     
```
