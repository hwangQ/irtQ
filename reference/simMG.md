# Simulated Multiple-Group Data

A simulated data set for three groups that take three test forms linked
by common items, for illustrating multiple-group calibration.

## Usage

``` r
simMG
```

## Format

A list with three elements:

- item.prm:

  A list of three data frames of item metadata, named `Group1`,
  `Group2`, and `Group3`. The forms of Groups 1 and 2 have 50 items (47
  3PLM items and three GRM items with five score categories), and the
  form of Group 3 has 38 items (37 3PLM items and one GRM item). Forms 1
  and 2 share 12 common items (`C1I1` to `C1I12`), forms 2 and 3 share
  10 common items (`C2I1` to `C2I10`), and forms 1 and 3 share no item.
  The item parameters are the true values used to simulate the
  responses, on the logistic metric (D = 1).

- res.dat:

  A list of three item response matrices, named `Group1` to `Group3`,
  with 2,000 rows (examinees) and 50, 50, and 38 columns (items in the
  order of `item.prm`). The abilities were drawn from \\N(0, 1)\\ (Group
  1), \\N(0.5, 0.8^{2})\\ (Group 2), and \\N(-0.3, 1.3^{2})\\ (Group 3).

- group.name:

  A character vector of the group names, `"Group1"`, `"Group2"`, and
  `"Group3"`.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# structure of the data
str(simMG, max.level = 1)
#> List of 3
#>  $ item.prm  :List of 3
#>  $ res.dat   :List of 3
#>  $ group.name: chr [1:3] "Group1" "Group2" "Group3"

# \donttest{
# multiple-group calibration of the three groups (see est_mg() for details)
x <- simMG$item.prm
fit <- est_mg(
  data = simMG$res.dat, group.name = simMG$group.name,
  model = list(x$Group1$model, x$Group2$model, x$Group3$model),
  cats = list(x$Group1$cats, x$Group2$cats, x$Group3$cats),
  item.id = list(x$Group1$id, x$Group2$id, x$Group3$id),
  D = 1, free.group = c(2, 3), use.gprior = TRUE,
  gprior = list(dist = "beta", params = c(5, 16)),
  group.mean = 0, group.var = 1, EmpHist = TRUE, Etol = 0.001, MaxE = 500,
  verbose = FALSE
)

# group means and variances
fit$group.par
#> $Group1
#>           mu sigma2 sigma
#> estimates  0      1     1
#> se        NA     NA    NA
#> 
#> $Group2
#>                   mu    sigma2      sigma
#> estimates 0.48545465 0.5575252 0.74667607
#> se        0.01669618 0.0176349 0.01180894
#> 
#> $Group3
#>                    mu     sigma2      sigma
#> estimates -0.36960845 1.95177274 1.39705860
#> se         0.03123918 0.06173591 0.02209496
#> 
# }
```
