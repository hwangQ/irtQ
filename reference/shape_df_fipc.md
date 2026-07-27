# Combine fixed and new item metadata for fixed-item parameter calibration (FIPC)

This function merges existing fixed-item metadata with automatically
generated metadata for new items, producing a single data frame ordered
by specified test positions, to facilitate fixed item parameter
calibration using
[`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md).

## Usage

``` r
shape_df_fipc(x, fix.loc = NULL, item.id = NULL, cats, model)
```

## Arguments

- x:

  A data.frame of metadata for items whose parameters remain fixed
  (e.g., output from
  [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)).

- fix.loc:

  An integer vector specifying the row positions in the final output
  where fixed items should be placed.

- item.id:

  A character vector of IDs for new items whose parameters will be
  estimated.If `NULL`, default IDs (e.g., "V1", "V2", ...) are assigned
  automatically.

- cats:

  An integer vector indicating the number of response categories for
  each new item; order must match `item.id`.

- model:

  A character vector of IRT model names for each new item. Valid options
  for dichotomous items: "1PLM", "2PLM", "3PLM", "DRM"; for polytomous
  items: "GRM", "GPCM".

## Value

A data.frame containing combined metadata for all items (fixed and new),
ordered by test position.

## Details

To use this function, first prepare a metadata frame `x` containing only
fixed items - either created by
[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md) or
imported from external software (e.g., via
[`bring.flexmirt()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md)),
which must include columns `id`, `cats`, `model`, and all relevant
parameter columns (`par.1`, `par.2`, etc.). The `fix.loc` argument
should then specify the exact row positions in the final test form where
these fixed items should remain. The length of `fix.loc` must match the
number of rows in `x`, and the order of positions in `fix.loc`
determines where each fixed-item row is placed.

Next, provide information for the new items whose parameters will be
estimated. Supply vectors for `item.id`, `cats`, and `model` matching
the number of new items (equal to total form length minus length of
`fix.loc`). If `item.id` is `NULL`, unique IDs are generated
automatically.

## See also

[`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Import the flexMIRT parameter output file
prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_fixed <- bring.flexmirt(file = prm_file, "par")$Group1$full_df

## Define positions of fixed items in the test form
fixed_pos <- c(1:40, 43:57)

## Specify IDs, models, and category counts for new items
new_ids <- paste0("NI", 1:6)
new_models <- c("3PLM", "1PLM", "2PLM", "GRM", "GRM", "GPCM")
new_cats <- c(2, 2, 2, 4, 5, 6)

## Generate combined metadata for FIPC
shape_df_fipc(x = x_fixed, fix.loc = fixed_pos, item.id = new_ids,
  cats = new_cats, model = new_models)
#>       id cats model     par.1       par.2        par.3      par.4      par.5
#> 1   CMC1    2  3PLM 0.7627842  1.46195307  0.261505353         NA         NA
#> 2   CMC2    2  3PLM 1.9212842 -1.04997761  0.175294361         NA         NA
#> 3   CMC3    2  3PLM 0.9266599  0.39475508  0.099746565         NA         NA
#> 4   CMC4    2  3PLM 1.0526419 -0.40704222  0.201193373         NA         NA
#> 5   CMC5    2  3PLM 0.8663350 -0.12481696  0.160420403         NA         NA
#> 6   CMC6    2  3PLM 1.6956204  0.62610659  0.072168572         NA         NA
#> 7   CMC7    2  3PLM 0.9061455  1.01573864  0.119458135         NA         NA
#> 8   CMC8    2  3PLM 0.8442812  0.80047702  0.114479251         NA         NA
#> 9   CMC9    2  3PLM 0.8541021  0.85236906  0.255248657         NA         NA
#> 10 CMC10    2  3PLM 1.5320152  0.09036640  0.135385451         NA         NA
#> 11 CMC11    2  3PLM 1.0018535 -0.46145320  0.132401405         NA         NA
#> 12 CMC12    2  3PLM 0.8784359  1.18144944  0.087687982         NA         NA
#> 13 CMC13    2  3PLM 1.4566944  1.40553715  0.181275480         NA         NA
#> 14 CMC14    2  3PLM 1.5142828  0.18435004  0.248064790         NA         NA
#> 15 CMC15    2  3PLM 1.2965707 -0.22914423  0.111920040         NA         NA
#> 16 CMC16    2  3PLM 2.0469491 -0.08536470  0.050890225         NA         NA
#> 17 CMC17    2  3PLM 1.3978190 -0.13421645  0.178526460         NA         NA
#> 18 CMC18    2  3PLM 1.6957073  1.24709400  0.273106846         NA         NA
#> 19 CMC19    2  3PLM 2.3117795 -1.01102523  0.179638721         NA         NA
#> 20 CMC20    2  3PLM 1.4489340 -1.65436976  0.190693574         NA         NA
#> 21 CMC21    2  3PLM 1.6281382 -1.19437963  0.121929251         NA         NA
#> 22 CMC22    2  3PLM 0.8320715 -0.67723050  0.197432042         NA         NA
#> 23 CMC23    2  3PLM 0.9768588 -0.25632988  0.130014303         NA         NA
#> 24 CMC24    2  3PLM 1.1418671  1.67963811  0.254011205         NA         NA
#> 25 CMC25    2  3PLM 0.7944843 -1.39426607  0.264555718         NA         NA
#> 26 CMC26    2  3PLM 1.0923706 -1.84806805  0.166491117         NA         NA
#> 27 CMC27    2  3PLM 1.1730201  0.07444365  0.131922668         NA         NA
#> 28 CMC28    2  3PLM 2.1546851 -0.08787145  0.208148649         NA         NA
#> 29 CMC29    2  3PLM 1.2810740 -1.38040199  0.201295671         NA         NA
#> 30 CMC30    2  3PLM 1.3505976  0.82351324  0.323145357         NA         NA
#> 31 CMC31    2  3PLM 0.8203985  0.70908638  0.082856907         NA         NA
#> 32 CMC32    2  3PLM 1.5235060 -0.88999000  0.257420730         NA         NA
#> 33 CMC33    2  3PLM 1.2711786 -1.31320312  0.186172388         NA         NA
#> 34 CMC34    2  3PLM 1.3107828  0.18690381  0.158109796         NA         NA
#> 35 CMC35    2  3PLM 1.4665995 -0.13801287  0.226328899         NA         NA
#> 36 CMC36    2  3PLM 0.8918806  1.09548554  0.126283949         NA         NA
#> 37 CMC37    2  3PLM 1.7287653 -0.40939097  0.104544594         NA         NA
#> 38 CMC38    2  3PLM 0.7327085 -0.37360246  0.218883035         NA         NA
#> 39  CFR1    5   GRM 1.9134815 -1.86983653 -1.238979525 -0.7140496 -0.2301269
#> 40  CFR2    5   GRM 1.2784658 -0.72406919 -0.068767894  0.5689640  1.0726913
#> 41   NI1    2  3PLM 1.0000000  0.00000000  0.200000000         NA         NA
#> 42   NI2    2  1PLM 1.0000000  0.00000000           NA         NA         NA
#> 43  AMC1    2  3PLM 1.4655489  0.64088029  0.225388960         NA         NA
#> 44  AMC2    2  3PLM 1.7609502 -1.52832130  0.158898368         NA         NA
#> 45  AMC3    2  3PLM 1.4430127  0.54078083  0.138825943         NA         NA
#> 46  AMC4    2  3PLM 0.9784240 -0.36582964  0.126958991         NA         NA
#> 47  AMC5    2  3PLM 0.9915383  2.37079859  0.164822203         NA         NA
#> 48  AMC6    2  3PLM 2.2732723  1.62427761  0.178409401         NA         NA
#> 49  AMC7    2  3PLM 1.2256893 -0.07172021  0.131317231         NA         NA
#> 50  AMC8    2  3PLM 1.6399644  0.16526401  0.182665802         NA         NA
#> 51  AMC9    2  3PLM 1.2090633  0.23974948  0.081366545         NA         NA
#> 52 AMC10    2  3PLM 1.3201851  1.33993665  0.079323850         NA         NA
#> 53 AMC11    2  3PLM 1.7372043 -0.99805101  0.246773029         NA         NA
#> 54 AMC12    2  3PLM 0.9746876 -0.72824708  0.223540392         NA         NA
#> 55  AFR1    5   GRM 1.1378595 -0.37401164  0.215487852  0.8485578  1.3826066
#> 56  AFR2    5   GRM 1.2333756 -2.07872217 -1.347637330 -0.7054699 -0.1163185
#> 57  AFR3    5   GRM 0.8762246 -0.75520683 -0.005929073  0.6650100  1.2474287
#> 58   NI3    2  2PLM 1.0000000  0.00000000           NA         NA         NA
#> 59   NI4    4   GRM 1.0000000  0.00000000  0.000000000  0.0000000         NA
#> 60   NI5    5   GRM 1.0000000  0.00000000  0.000000000  0.0000000  0.0000000
#> 61   NI6    6  GPCM 1.0000000  0.00000000  0.000000000  0.0000000  0.0000000
#>    par.6
#> 1     NA
#> 2     NA
#> 3     NA
#> 4     NA
#> 5     NA
#> 6     NA
#> 7     NA
#> 8     NA
#> 9     NA
#> 10    NA
#> 11    NA
#> 12    NA
#> 13    NA
#> 14    NA
#> 15    NA
#> 16    NA
#> 17    NA
#> 18    NA
#> 19    NA
#> 20    NA
#> 21    NA
#> 22    NA
#> 23    NA
#> 24    NA
#> 25    NA
#> 26    NA
#> 27    NA
#> 28    NA
#> 29    NA
#> 30    NA
#> 31    NA
#> 32    NA
#> 33    NA
#> 34    NA
#> 35    NA
#> 36    NA
#> 37    NA
#> 38    NA
#> 39    NA
#> 40    NA
#> 41    NA
#> 42    NA
#> 43    NA
#> 44    NA
#> 45    NA
#> 46    NA
#> 47    NA
#> 48    NA
#> 49    NA
#> 50    NA
#> 51    NA
#> 52    NA
#> 53    NA
#> 54    NA
#> 55    NA
#> 56    NA
#> 57    NA
#> 58    NA
#> 59    NA
#> 60    NA
#> 61     0
```
