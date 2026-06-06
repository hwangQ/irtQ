# Simulated Mixed-Item Format CAT Data

A simulated dataset for computerized adaptive testing (CAT), containing
an item pool, sparse response data, and examinee ability estimates. The
item pool includes both dichotomous and polytomous items.

## Usage

``` r
simCAT_MX
```

## Format

A list of length three:

- item_pool:

  A data frame in item metadata format consisting of 200 dichotomous
  items and 30 polytomous items.

  - Dichotomous items: Calibrated using the IRT 3PL model.

  - Polytomous items: Calibrated using the Generalized Partial Credit
    Model (GPCM), with three score categories (0, 1, 2).

- response_data:

  A sparse matrix of item responses from 30,000 examinees.

- theta_estimates:

  A numeric vector of ability estimates for the 30,000 examinees.

## Author

Hwanggyu Lim <hglim83@gmail.com>
