# Simulated 1-3-3 MST Panel Data

A simulated Multistage-Adaptive Test (MST) dataset representing a 1-3-3
panel structure, assembled using module-level target test information
functions (TIFs) and design constraints similar to those used in the
simulation study by Lim et al. (2021).

## Usage

``` r
simMST
```

## Format

A list containing five internal objects:

- item_bank:

  A data frame of item metadata for the 56 items assigned to the panel
  (one row per item; columns follow the irtQ item metadata format: `id`,
  `cats`, `model`, `par.1`, `par.2`, `par.3`).

- module:

  A binary matrix that maps items in the item bank to MST modules. This
  structure specifies the item-to-module assignment, similar to the
  `modules` argument in the `randomMST()` function from the mstR package
  (Magis et al., 2017).

- route_map:

  A binary square matrix that defines the MST transition structure,
  showing module pathways across stages. This corresponds to the
  `transMatrix` argument in the `randomMST()` function in the mstR
  package.

- cut_score:

  A list of numeric vectors specifying the routing cut scores between
  MST stages, computed via
  [`find_cut()`](https://hwangQ.github.io/irtQ/reference/find_cut.md)
  based on TIF-crossing points of the assembled modules.

- theta:

  A numeric vector of ability (theta) values used to evaluate the
  panel's measurement precision across the latent trait continuum.

This 1-3-3 MST panel includes 7 modules across 3 stages, drawn from 56
unique items (8 per module) with no item shared across modules. Each
module contains 8 dichotomously scored items calibrated under the IRT
3-parameter logistic (3PL) model.

## References

Magis, D., Yan, D., & von Davier, A. A. (2017). *Computerized adaptive
and multistage testing with R: Using packages catR and mstR*. Springer.

Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
approach to evaluate the performance of MST. *Journal of Educational
Measurement, 58*(2), 154-178.

## Author

Hwanggyu Lim <hglim83@gmail.com>
