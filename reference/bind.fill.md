# Bind Fill

This function creates a matrix using either row-wise (`rbind`) or
column-wise (`cbind`) binding of a list of numeric vectors with varying
lengths. Shorter vectors are padded with a specified fill value.

## Usage

``` r
bind.fill(List, type = c("rbind", "cbind"), fill = NA)
```

## Arguments

- List:

  A list containing numeric vectors of possibly different lengths.

- type:

  A character string indicating the type of binding to perform. Options
  are `"rbind"` or `"cbind"`.

- fill:

  A value used to fill missing elements when aligning the vectors. For
  `type = "cbind"`, this fills missing rows in shorter columns; for
  `type = "rbind"`, this fills missing columns in shorter rows. Accepts
  any R object (e.g., numeric, character, logical). Default is `NA`.

## Value

A matrix formed by binding the elements of the list either row-wise or
column-wise, with shorter vectors padded by the specified `fill` value.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# Sample list
score_list <- list(item1 = 0:3, item2 = 0:2, item3 = 0:5, item4 = 0:4)

# 1) Create a row-bound matrix (rbind)
bind.fill(score_list, type = "rbind")
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    1    2    3   NA   NA
#> [2,]    0    1    2   NA   NA   NA
#> [3,]    0    1    2    3    4    5
#> [4,]    0    1    2    3    4   NA

# 2) Create a column-bound matrix (cbind)
bind.fill(score_list, type = "cbind")
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    1    1    1    1
#> [3,]    2    2    2    2
#> [4,]    3   NA    3    3
#> [5,]   NA   NA    4    4
#> [6,]   NA   NA    5   NA

# 3) Create a column-bound matrix and fill missing values with 0
bind.fill(score_list, type = "cbind", fill = 0L)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    1    1    1    1
#> [3,]    2    2    2    2
#> [4,]    3    0    3    3
#> [5,]    0    0    4    4
#> [6,]    0    0    5    0
```
