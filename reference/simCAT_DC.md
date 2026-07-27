# Simulated Single-Item Format CAT Data

A simulated dataset containing an item pool, sparse response data, and
examinee ability estimates, designed for single-item computerized
adaptive testing (CAT).

## Usage

``` r
simCAT_DC
```

## Format

A list of length three:

- item_pool:

  A data frame in item metadata format containing 100 dichotomous items.

  - Items 1-90: Generated and calibrated under the IRT 2PL model.

  - Items 91-100: Generated under the IRT 3PL model but calibrated using
    the 2PL model.

- response_data:

  A sparse matrix of item responses from 10,000 examinees.

- theta_estimates:

  A numeric vector of ability estimates for the 10,000 examinees.

## Author

Hwanggyu Lim <hglim83@gmail.com>
