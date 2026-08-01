# Score Selected-Response Item Data Against an Answer Key

This function converts raw selected-response (e.g., multiple-choice)
item data - coded as the selected option (an option number, e.g., 1-5
for a five-option item; a single-letter option label, e.g., `"A"`-`"E"`;
or any other option label in a different script or notation, such as
Korean syllable labels or Roman numerals), a missing-response indicator
for an omitted response, or a comma-separated string (e.g., `"1,5"` or
`"A,C"`) for a double-marked (multiple-option) response - into a
dichotomously scored item-response matrix (0 = incorrect, 1 = correct)
using a supplied answer key. The option coding scheme is detected
independently for each item from its own `key` value, so a single test
form may freely mix numerically coded, letter-coded, and other
label-coded items. Letter (and other case-bearing) option labels are
matched case-insensitively (e.g., a response of `"a"` matches a key of
`"A"`). Omitted and double-marked responses are both scored as incorrect
(0), and their frequencies are tabulated separately for reporting.

## Usage

``` r
score_resp(data, key, missing = NA)
```

## Arguments

- data:

  A data frame or matrix of raw item responses, with examinees in rows
  and items in columns (in the same left-to-right order as `key`). Do
  not include non-item columns (e.g., group or examinee identifiers);
  subset those out before calling `score_resp()`. Each cell should
  contain either a single selected option - an option number (e.g.,
  `2`), a single-letter option label (e.g., `"C"`), or any other option
  label (e.g., a Korean syllable label, or a Roman numeral such as
  `"II"`) - matching whichever coding scheme is used for that item's
  `key` value, a missing-response indicator (see `missing`), or a
  character string of comma-separated options (e.g., `"1,5"` or `"A,C"`)
  for a double-marked response. To keep double-marked responses as
  character strings rather than having them coerced to `NA` on import,
  read the source file with all item columns imported as text (e.g.,
  `readxl::read_excel(..., col_types = "text")`) before passing `data`
  to this function.

- key:

  Either (a) a vector of correct options in item order (length must
  equal `ncol(data)`), where each element is a numeric/integer option
  number (e.g., `4`), a single-letter option label (e.g., `"D"`), or any
  other non-blank option label (e.g., a Korean syllable label, or a
  Roman numeral such as `"II"`), or (b) a data frame with columns `item`
  and `key` giving the item number and its correct option, respectively
  (the `key` column may likewise mix option formats across rows). Each
  item's coding scheme is inferred independently from its own key value,
  so items using different option formats may be freely mixed within the
  same vector or data frame; letter (and other case-bearing) values are
  matched case-insensitively. When a data frame is supplied, it is
  internally sorted by `item` before use, so its row order does not need
  to match the column order of `data`.

- missing:

  A value indicating missing (omitted) responses in `data`, analogous to
  the `missing` argument in
  [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) and
  [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md).
  Any cell equal to `missing` is recoded to `NA` before scoring. Default
  is `NA`, meaning `data` is assumed to already use `NA` (or an empty
  string) for omitted responses, and no recoding is performed. Set this
  to whatever sentinel value a particular data set uses for a missing
  response (e.g., `-9`, `"9"`) when it differs from `NA`.

## Value

A list with two elements:

- scored:

  A data frame with the same dimensions as `data`, containing the
  dichotomously scored (0/1) item responses. Column names follow
  `colnames(data)` (or `V1, V2, ...` when `data` has no column names).

- resp_summary:

  A data frame with one row per item, reporting the number and
  percentage of blank and double-marked responses, the number of
  unrecognized/invalid response tokens (`n_invalid`), and the number of
  examinees scored correct/incorrect, for each item.

## Details

For each item, a response is scored as correct (1) only when it is a
single, non-missing option that matches the corresponding value in
`key`. All other cases - an omitted response, a double-marked response,
or a single but incorrect option - are scored as incorrect (0). Missing
and double-marked responses are tallied separately in `resp_summary` so
that omission and double-marking rates can be reported independently,
even though both are scored as 0.

Each item's option coding scheme is determined independently from its
own `key` value, using one of three rules:

- If the key value parses as a number (e.g., `4`), responses to that
  item are compared numerically, exactly as in earlier versions of this
  function.

- If the key value consists of one or more Latin letters (e.g., `"D"`),
  responses are compared as letters, case-insensitively (`"d"` and `"D"`
  are treated as the same option), and a response must itself consist
  only of Latin letters to count as a single valid option (see below for
  what happens otherwise).

- Otherwise (e.g., a Korean syllable label, a Roman numeral written with
  non-Latin numeral characters, a circled-number symbol, or any other
  non-numeric, non-Latin-letter label), the key value is treated as a
  general option label: any non-blank, non-double-marked response to
  that item is accepted as a single valid option, and is compared to the
  key with a case-insensitive
  ([`toupper()`](https://rdrr.io/r/base/chartr.html)-based, which only
  affects any embedded Latin characters) exact string match. There is no
  universal, script-independent rule for "a well-formed option label"
  analogous to the numeric or Latin-letter checks above, so this scheme
  cannot distinguish a genuinely incorrect option from a garbled
  response token; both are simply scored 0, and `n_invalid` is always 0
  for items using this scheme (see below).

This per-item detection means a single test form may freely mix items
using any of the three schemes - only the value supplied in `key` for a
given item determines how that item's column in `data` is interpreted. A
`key` value that is blank or `NA` causes `score_resp()` to stop with an
error, since no correct option was actually supplied for that item.

A response cell is classified as double-marked when, after coercion to
character, it contains a comma (e.g., `"1,5"`, `"B,D"`). A response cell
is classified as blank/omitted when it equals `missing` (recoded to `NA`
beforehand), or is otherwise `NA` or an empty/whitespace-only string.
For numeric- and Latin-letter-coded items, any response that is none of
blank, double-marked, or a single valid option in the coding scheme used
for that item (e.g., a numeric token where a letter was expected for a
letter-coded item, or stray text from a data-entry artifact) is
classified as invalid; it is still scored 0, but tallied separately in
`n_invalid` (with a [`warning()`](https://rdrr.io/r/base/warning.html)
raised) so it is visible rather than silently merged into the "wrong"
count. As noted above, general label-coded items never populate
`n_invalid`, since any non-blank, non-double-marked token is accepted as
a single valid (though possibly incorrect) option for that scheme. Note
that `n_correct + n_wrong` always equals `n`, but `n_blank`, `n_double`,
and `n_invalid` are reported as separate diagnostic tallies, not as a
four-way disjoint partition of `n_wrong`.

This function performs dichotomous (single-key, selected-response)
scoring only, since a selected-response item has exactly one correct
option by design and therefore only two possible score categories
(correct/incorrect). It is not generalized to items with more than two
score categories (e.g., partial-credit polytomous items). If such data
become available, they should be supplied directly, already scored, to
downstream item-analysis functions that support polytomous data,
bypassing `score_resp()`.

## Note

This function applies to **dichotomous, single-key, selected-response
items only** (i.e., items with exactly one correct option, scored into
exactly two categories: 0 = incorrect and 1 = correct). It does not
generalize to items with more than two score categories (e.g.,
partial-credit polytomous items scored 0, 1, 2, ...), since a
selected-response item cannot yield a partial-credit score by
construction. If polytomous item response data are available, supply an
already-scored response matrix (with categories 0, 1, 2, ..., cats - 1)
directly to downstream item-analysis functions instead of using
`score_resp()`. See **Details** for further discussion.

## See also

[`ctt_distr()`](https://hwangQ.github.io/irtQ/reference/ctt_distr.md),
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md),
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# A small 5-examinee, 3-item toy example with a blank and a double-mark
raw <- data.frame(
  V1 = c("1", "2", NA, "1", "3"),
  V2 = c("4", "4", "4", "2,4", "4"),
  V3 = c("2", "1", "3", "3", NA)
)
key <- c(1, 4, 3)
out <- score_resp(data = raw, key = key)
out$scored
#>   V1 V2 V3
#> 1  1  1  0
#> 2  0  1  0
#> 3  0  1  1
#> 4  1  0  1
#> 5  0  1  0
out$resp_summary
#>   item key n n_correct n_wrong n_blank n_double n_invalid pct_blank pct_double
#> 1   V1   1 5         2       3       1        0         0        20          0
#> 2   V2   4 5         4       1       0        1         0         0         20
#> 3   V3   3 5         2       3       1        0         0        20          0

# Same data, but omitted responses are coded as "9" instead of NA
raw2 <- data.frame(
  V1 = c("1", "2", "9", "1", "3"),
  V2 = c("4", "4", "4", "2,4", "4"),
  V3 = c("2", "1", "3", "3", "9")
)
out2 <- score_resp(data = raw2, key = key, missing = "9")
identical(out$scored, out2$scored)
#> [1] TRUE

# A letter-coded (A-E) example, including a double-marked response and a
# lowercase response ("a"), which is matched case-insensitively
raw3 <- data.frame(
  V1 = c("A", "B", NA, "a", "C"),
  V2 = c("D", "D", "D", "B,D", "D"),
  V3 = c("B", "A", "C", "C", NA)
)
key3 <- c("A", "D", "C")
out3 <- score_resp(data = raw3, key = key3)
out3$scored
#>   V1 V2 V3
#> 1  1  1  0
#> 2  0  1  0
#> 3  0  1  1
#> 4  1  0  1
#> 5  0  1  0

# A mixed-format test: item 1 is numerically coded, items 2-3 are
# letter-coded; the coding scheme is inferred separately for each item
# from its own key value
raw4 <- data.frame(
  V1 = c("1", "2", "1", "1", "3"),
  V2 = c("D", "D", "D", "B,D", "D"),
  V3 = c("B", "A", "C", "C", NA)
)
key4 <- c(1, "D", "C")
out4 <- score_resp(data = raw4, key = key4)
out4$scored
#>   V1 V2 V3
#> 1  1  1  0
#> 2  0  1  0
#> 3  1  1  1
#> 4  1  0  1
#> 5  0  1  0

# A full five-option (A-E) multiple-choice item, scored against key "C"
raw5 <- data.frame(
  V1 = c("A", "B", "C", "D", "E", "c")
)
key5 <- c("C")
out5 <- score_resp(data = raw5, key = key5)
out5$scored
#>   V1
#> 1  0
#> 2  0
#> 3  1
#> 4  0
#> 5  0
#> 6  1

# A general-label-coded example using non-Latin option labels: the Korean
# syllables romanized "ga"/"na"/"da"/"ra"/"ma", written here with \uXXXX
# escapes for ASCII portability. Because these labels are neither numbers
# nor Latin letters, they are handled by the general label-coding scheme
# described in Details
raw6 <- data.frame(
  V1 = c("\uAC00", "\uB098", "\uB2E4", "\uB2E4", NA) # ga, na, da, da, NA
)
key6 <- c("\uB2E4") # da
out6 <- score_resp(data = raw6, key = key6)
out6$scored
#>   V1
#> 1  0
#> 2  0
#> 3  1
#> 4  1
#> 5  0
```
