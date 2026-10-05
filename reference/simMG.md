# Simulated multiple-group data

This data set has a list consisting of item metadata, item response
data, and group names of three simulated groups.

## Usage

``` r
simMG
```

## Format

This data set includes a list of three internal objects: (1) a list of
item metadata (item.prm) for three groups, (2) a list of item response
data (res.dat) for the three groups, and (3) a vector of group names
(group.name) for the three groups.

The first internal object (item.prm) contains a list of item metadata of
three test forms for the three groups. In terms of test forms, the test
forms for the first and second groups have fifty items consisting of
forty seven 3PLM items and three GRM items. The test form for the third
group has thirty eight items consisting of thirty seven 3PLM items and
one GRM item. Among the three forms, the first and second test forms
share twelve common items (C1I1 through C1I12) and the second and third
test forms share ten common items (C2I1 through c2I10). There is no
common item between the first and third forms. The item parameters in
the item metadata were used to simulate the item response data sets for
the three groups (see the second object of the list).

Regrading the second internal object, all three response data sets were
simulated with 2,000 latent abilities randomly sampled from \\N(0, 1)\\
(Group 1), \\N(0.5, 0.8^{2})\\ (Group 2), and \\N(-0.3, 1.3^{2})\\
(Group 3), respectively, using the true item parameters provided in the
item metadata.

The third internal object is a vector of three group names which are
"Group1", "Group2", and "Group3".

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
