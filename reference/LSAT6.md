# LSAT6 Data

A well-known dichotomous response dataset from the Law School Admission
Test (LSAT), Section 6, as used in Thissen (1982).

## Usage

``` r
LSAT6
```

## Format

An integer matrix with 1,000 rows (examinees) and 5 columns (items
`item.1` to `item.5`). Each entry is a dichotomously scored response (0
= incorrect, 1 = correct). The 1,000 rows contain 30 distinct response
patterns.

## References

Thissen, D. (1982). Marginal maximum likelihood estimation for the
one-parameter logistic model. *Psychometrika, 47*, 175-186.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# structure of the data
head(LSAT6)
#>   item.1 item.2 item.3 item.4 item.5
#> 1      0      0      0      0      0
#> 2      0      0      0      0      0
#> 3      0      0      0      0      0
#> 4      0      0      0      0      1
#> 5      0      0      0      0      1
#> 6      0      0      0      0      1

# fit the 2PL model to the LSAT6 data
est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
#> 
#> Call:
#> est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
#> 
#> Item parameter estimation using MMLE-EM. 
#> 50 E-step cycles were completed using 49 quadrature points.
#> First-order test: Convergence criteria are satisfied.
#> Second-order test: Solution is a possible local maximum.
#> Computation of variance-covariance matrix: 
#>   Variance-covariance matrix of item parameter estimates is obtainable.
#> 
#> Log-likelihood: -2466.653
#> 
```
