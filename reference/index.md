# Package index

## Getting Started

New to irtQ? Start here with core data structures.

- [`irtQ-package`](https://hwangQ.github.io/irtQ/reference/irtQ-package.md)
  [`irtQ`](https://hwangQ.github.io/irtQ/reference/irtQ-package.md) :
  irtQ: Unidimensional Item Response Theory Modeling
- [`shape_df()`](https://hwangQ.github.io/irtQ/reference/shape_df.md) :
  Create a Data Frame of Item Metadata
- [`shape_df_fipc()`](https://hwangQ.github.io/irtQ/reference/shape_df_fipc.md)
  : Combine fixed and new item metadata for fixed-item parameter
  calibration (FIPC)

## Item Parameter Estimation

Calibrate item parameters using marginal maximum likelihood (MMLE-EM)
and fixed ability methods.

- [`est_irt()`](https://hwangQ.github.io/irtQ/reference/est_irt.md) :
  Item parameter estimation using MMLE-EM algorithm
- [`est_mg()`](https://hwangQ.github.io/irtQ/reference/est_mg.md) :
  Multiple-group item calibration using MMLE-EM algorithm
- [`est_item()`](https://hwangQ.github.io/irtQ/reference/est_item.md) :
  Fixed ability parameter calibration
- [`getirt()`](https://hwangQ.github.io/irtQ/reference/getirt.md) :
  Extract Components from 'est_irt', 'est_mg', or 'est_item' Objects

## Ability Estimation & Scoring

Estimate examinee ability scores from item response data.

- [`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md)
  : Estimate examinees' ability (proficiency) parameters
- [`llike_score()`](https://hwangQ.github.io/irtQ/reference/llike_score.md)
  : Log-Likelihood of Ability Parameters

## Model-Data Fit Evaluation

Assess how well the fitted IRT model describes the observed data.

- [`irtfit()`](https://hwangQ.github.io/irtQ/reference/irtfit.md) :
  Traditional IRT Item Fit Statistics
- [`sx2_fit()`](https://hwangQ.github.io/irtQ/reference/sx2_fit.md) :
  S-X2 Fit Statistic
- [`plot(`*`<irtfit>`*`)`](https://hwangQ.github.io/irtQ/reference/plot.irtfit.md)
  : Draw Raw and Standardized Residual Plots

## Differential Item Functioning (DIF)

Detect item bias across examinee groups using residual-based and CATSIB
methods.

- [`rdif()`](https://hwangQ.github.io/irtQ/reference/rdif.md) : IRT
  Residual-Based Differential Item Functioning (RDIF) Detection
  Framework
- [`crdif()`](https://hwangQ.github.io/irtQ/reference/crdif.md) :
  Residual-Based DIF Detection Framework Using Categorical Residuals
  (RDIF-CR)
- [`grdif()`](https://hwangQ.github.io/irtQ/reference/grdif.md) :
  Generalized IRT residual-based DIF detection framework for multiple
  groups (GRDIF)
- [`catsib()`](https://hwangQ.github.io/irtQ/reference/catsib.md) :
  CATSIB DIF Detection Procedure
- [`ripd()`](https://hwangQ.github.io/irtQ/reference/ripd.md) :
  Residual-based Item Parameter Drift (RIPD) Detection Framework
- [`pcd2()`](https://hwangQ.github.io/irtQ/reference/pcd2.md) :
  Pseudo-count D2 method

## Classification Accuracy & Consistency

Evaluate the reliability of cut-score-based pass/fail classifications.

- [`cac_lee()`](https://hwangQ.github.io/irtQ/reference/cac_lee.md) :
  Classification Accuracy and Consistency Using Lee's (2010) Approach
- [`cac_rud()`](https://hwangQ.github.io/irtQ/reference/cac_rud.md) :
  Classification Accuracy and Consistency Based on Rudner's (2001, 2005)
  Approach

## Information & Characteristic Functions

Compute and visualize item/test information functions and characteristic
curves.

- [`info()`](https://hwangQ.github.io/irtQ/reference/info.md) : Item and
  Test Information Function
- [`traceline()`](https://hwangQ.github.io/irtQ/reference/traceline.md)
  : Compute Item/Test Characteristic Functions
- [`plot(`*`<info>`*`)`](https://hwangQ.github.io/irtQ/reference/plot.info.md)
  : Plot Item and Test Information Functions
- [`plot(`*`<traceline>`*`)`](https://hwangQ.github.io/irtQ/reference/plot.traceline.md)
  : Plot Item and Test Characteristic Curves

## Data Simulation & Utilities

Simulate item response data and compute supporting quantities.

- [`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) :
  Simulated Response Data
- [`lwrc()`](https://hwangQ.github.io/irtQ/reference/lwrc.md) :
  Lord-Wingersky Recursion Formula
- [`gen.weight()`](https://hwangQ.github.io/irtQ/reference/gen.weight.md)
  : Generate Weights
- [`covirt()`](https://hwangQ.github.io/irtQ/reference/covirt.md) :
  Asymptotic Variance-Covariance Matrices of Item Parameter Estimates
- [`bind.fill()`](https://hwangQ.github.io/irtQ/reference/bind.fill.md)
  : Bind Fill
- [`bisection()`](https://hwangQ.github.io/irtQ/reference/bisection.md)
  : The Bisection Method to Find a Root
- [`reval_mst()`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
  : Recursion-based MST evaluation method
- [`summary()`](https://hwangQ.github.io/irtQ/reference/summary.md) :
  Summary of Item Calibration Results

## Importing External Software Output

Read item parameter estimates from IRT software output files.

- [`bring.flexmirt()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md)
  [`bring.bilog()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md)
  [`bring.parscale()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md)
  [`bring.mirt()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md)
  : Import Item and Ability Parameters from IRT Software
- [`run_flexmirt()`](https://hwangQ.github.io/irtQ/reference/run_flexmirt.md)
  : Run flexMIRT from Within R
- [`write.flexmirt()`](https://hwangQ.github.io/irtQ/reference/write.flexmirt.md)
  : Write a "-prm.txt" File for flexMIRT

## IRT Model Probability Functions

Core IRT model functions for computing item response probabilities.

- [`drm()`](https://hwangQ.github.io/irtQ/reference/drm.md) :
  Dichotomous Response Model (DRM) Probabilities
- [`prm()`](https://hwangQ.github.io/irtQ/reference/prm.md) : Polytomous
  Response Model (PRM) Probabilities (GRM and GPCM)

## Datasets

Example datasets bundled with the package.

- [`LSAT6`](https://hwangQ.github.io/irtQ/reference/LSAT6.md) : LSAT6
  Data
- [`simCAT_DC`](https://hwangQ.github.io/irtQ/reference/simCAT_DC.md) :
  Simulated Single-Item Format CAT Data
- [`simCAT_MX`](https://hwangQ.github.io/irtQ/reference/simCAT_MX.md) :
  Simulated Mixed-Item Format CAT Data
- [`simMG`](https://hwangQ.github.io/irtQ/reference/simMG.md) :
  Simulated multiple-group data
- [`simMST`](https://hwangQ.github.io/irtQ/reference/simMST.md) :
  Simulated 1-3-3 MST Panel Data
