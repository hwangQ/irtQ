# Import Item and Ability Parameters from IRT Software

These functions import item and/or ability parameters from BILOG-MG 3,
PARSCALE 4, flexMIRT, and the mirt R package.

## Usage

``` r
bring.flexmirt(
  file,
  type = c("par", "sco"),
  rePar = TRUE,
  rePar.gpc = TRUE,
  n.factor = 1
)

bring.bilog(file, type = c("par", "sco"))

bring.parscale(file, type = c("par", "sco"))

bring.mirt(x)
```

## Arguments

- file:

  A file name (including the full path) containing the item or ability
  parameter estimates.

- type:

  A character string indicating the type of output file. Available
  options are `"par"` for item parameter files and `"sco"` for ability
  parameter files.

- rePar:

  Logical. If `TRUE`, and when a dichotomous IRT model (e.g., 3PLM) or
  the graded response model (GRM) is fit, the item intercept and the
  logit of the guessing parameter are reparameterized into the item
  difficulty and guessing parameters, respectively. Default is `TRUE`.

- rePar.gpc:

  Logical. If `TRUE`, and when the generalized partial credit model
  (GPCM) is fit, the nominal model parameters in the flexMIRT output are
  reparameterized into GPCM slope and difficulty parameters. Default is
  `TRUE`.

- n.factor:

  A numeric value indicating the number of latent traits (factors)
  estimated. This argument must be specified when `type = "sco"`.
  Default is 1.

- x:

  An object returned by the function
  [`mirt::mirt()`](https://philchalmers.github.io/mirt/reference/mirt.html).

## Value

These functions return a list containing several components. For
flexMIRT output files, results from multiple-group analyses can be
handled; in such cases, each list element corresponds to the estimation
results for a separate group.

## Details

The `bring.flexmirt()` function was developed by modifying the
`read.flexmirt()` function (Pritikin & Falk, 2020). Similarly,
`bring.bilog()` and `bring.parscale()` were based on modified versions
of the `read.bilog()` and `read.parscale()` functions (Weeks, 2010),
respectively.

The supported file extensions for item and ability parameter files are:
".par" and ".sco" for BILOG-MG and PARSCALE, and "-prm.txt" and
"-sco.txt" for flexMIRT. For mirt, the user provides the object name
directly.

Although `bring.flexmirt()` can extract multidimensional item and
ability parameter estimates, the irtQ package is designed exclusively
for unidimensional IRT applications.

For polytomous items, both `bring.flexmirt()` and `bring.mirt()` can
import item parameters for the graded response model (GRM) and the
generalized partial credit model (GPCM).

## Note

For item parameter files from any IRT software, only the internal object
`"full_df"` in the returned list is required for various functions in
the irtQ package. This object is a data frame containing item metadata
(e.g., item parameters, number of categories, IRT model types). See
[`info()`](https://hwangQ.github.io/irtQ/reference/info.md) or
[`simdat()`](https://hwangQ.github.io/irtQ/reference/simdat.md) for more
details on item metadata.

In addition, when item parameters are estimated using the partial credit
model (PCM) or the generalized partial credit model (GPCM), item step
parameters are included in the `"full_df"` object. These step parameters
are calculated by subtracting the category threshold parameters from the
overall item difficulty (or location) parameter. See the **IRT Models**
section in
[irtQ-package](https://hwangQ.github.io/irtQ/reference/irtQ-package.md)
for further details on the parameterization of the GPCM.

## Sample Output Files of IRT software

To illustrate how to import item parameter estimate files from PARSCALE
4 and flexMIRT using `bring.parscale()` and `bring.flexmirt()`, two
example output files are included in this package.

One file is from PARSCALE 4 with a ".PAR" extension (i.e.,
"parscale_sample.PAR"), and the other is from flexMIRT with a "-prm.txt"
extension (i.e., "flexmirt_sample-prm.txt").

Both files contain item parameter estimates from a mixed-format test
with 55 items: fifty dichotomous items following the 3PL model and five
polytomous items with five response categories modeled using the graded
response model (GRM). The examples below demonstrate how to import these
output files.

## References

Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item
analysis and test scoring (Computer software). Chapel Hill, NC: Vector
Psychometric Group.

Chalmers, R. P. (2012). mirt: A multidimensional item response theory
package for the R environment. *Journal of Statistical Software, 48*(6),
1-29.

Weeks, J. P. (2010). plink: An R Package for Linking Mixed-Format Tests
Using IRT-Based Methods. *Journal of Statistical Software, 35*(12),
1-33. URL http://www.jstatsoft.org/v35/i12/.

Pritikin, J. (2018). *rpf: Response Probability Functions*. R package
version 0.59. https://CRAN.R-project.org/package=rpf.

Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research
environment for item response theory method development. *Applied
Psychological Measurement, 44*(7-8), 561-562.

Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test
scoring for rating scale data (Computer Software). Chicago, IL:
Scientific Software International. URL http://www.ssicentral.com

Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003).
BILOG-MG 3: Multiple-group IRT analysis and test maintenance for binary
items (Computer Software). Chicago, IL: Scientific Software
International. URL http://www.ssicentral.com

## See also

[irtQ-package](https://hwangQ.github.io/irtQ/reference/irtQ-package.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1
# Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Read item parameters and convert them to item metadata
bring.flexmirt(file = flex_sam, "par")$Group1$full_df
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
#> 41  AMC1    2  3PLM 1.4655489  0.64088029  0.225388960         NA         NA
#> 42  AMC2    2  3PLM 1.7609502 -1.52832130  0.158898368         NA         NA
#> 43  AMC3    2  3PLM 1.4430127  0.54078083  0.138825943         NA         NA
#> 44  AMC4    2  3PLM 0.9784240 -0.36582964  0.126958991         NA         NA
#> 45  AMC5    2  3PLM 0.9915383  2.37079859  0.164822203         NA         NA
#> 46  AMC6    2  3PLM 2.2732723  1.62427761  0.178409401         NA         NA
#> 47  AMC7    2  3PLM 1.2256893 -0.07172021  0.131317231         NA         NA
#> 48  AMC8    2  3PLM 1.6399644  0.16526401  0.182665802         NA         NA
#> 49  AMC9    2  3PLM 1.2090633  0.23974948  0.081366545         NA         NA
#> 50 AMC10    2  3PLM 1.3201851  1.33993665  0.079323850         NA         NA
#> 51 AMC11    2  3PLM 1.7372043 -0.99805101  0.246773029         NA         NA
#> 52 AMC12    2  3PLM 0.9746876 -0.72824708  0.223540392         NA         NA
#> 53  AFR1    5   GRM 1.1378595 -0.37401164  0.215487852  0.8485578  1.3826066
#> 54  AFR2    5   GRM 1.2333756 -2.07872217 -1.347637330 -0.7054699 -0.1163185
#> 55  AFR3    5   GRM 0.8762246 -0.75520683 -0.005929073  0.6650100  1.2474287

## Example 2
# Import the ".PAR" output file from PARSCALE
pscale_sam <- system.file("extdata", "parscale_sample.PAR", package = "irtQ")

# Read item parameters and convert them to item metadata
bring.parscale(file = pscale_sam, "par")$full_df
#>      id cats model   par.1    par.2    par.3    par.4    par.5
#> 1  0001    2   DRM 0.78668  1.45392  0.26424       NA       NA
#> 2  0002    2   DRM 1.95256 -1.03590  0.16790       NA       NA
#> 3  0003    2   DRM 0.92423  0.36744  0.08911       NA       NA
#> 4  0004    2   DRM 1.06017 -0.41665  0.19257       NA       NA
#> 5  0005    2   DRM 0.86593 -0.15680  0.14700       NA       NA
#> 6  0006    2   DRM 1.70781  0.62067  0.06983       NA       NA
#> 7  0007    2   DRM 0.91006  0.99884  0.11532       NA       NA
#> 8  0008    2   DRM 0.84174  0.77277  0.10594       NA       NA
#> 9  0009    2   DRM 0.86645  0.84402  0.25407       NA       NA
#> 10 0010    2   DRM 1.54383  0.08762  0.13141       NA       NA
#> 11 0011    2   DRM 1.00673 -0.47776  0.12047       NA       NA
#> 12 0012    2   DRM 0.87475  1.15744  0.08077       NA       NA
#> 13 0013    2   DRM 1.48901  1.39467  0.18195       NA       NA
#> 14 0014    2   DRM 1.53152  0.18416  0.24621       NA       NA
#> 15 0015    2   DRM 1.30079 -0.24027  0.10225       NA       NA
#> 16 0016    2   DRM 2.06064 -0.08604  0.04616       NA       NA
#> 17 0017    2   DRM 1.41035 -0.13591  0.17396       NA       NA
#> 18 0018    2   DRM 1.73702  1.23955  0.27421       NA       NA
#> 19 0019    2   DRM 2.35588 -0.99336  0.17489       NA       NA
#> 20 0020    2   DRM 1.47074 -1.64613  0.17122       NA       NA
#> 21 0021    2   DRM 1.65023 -1.18579  0.10906       NA       NA
#> 22 0022    2   DRM 0.83234 -0.71780  0.17873       NA       NA
#> 23 0023    2   DRM 0.97745 -0.28252  0.11614       NA       NA
#> 24 0024    2   DRM 1.19351  1.66222  0.25751       NA       NA
#> 25 0025    2   DRM 0.80015 -1.42399  0.24565       NA       NA
#> 26 0026    2   DRM 1.10754 -1.84767  0.14416       NA       NA
#> 27 0027    2   DRM 1.17525  0.06034  0.12368       NA       NA
#> 28 0028    2   DRM 2.17855 -0.08451  0.20623       NA       NA
#> 29 0029    2   DRM 1.29244 -1.39057  0.17961       NA       NA
#> 30 0030    2   DRM 1.37376  0.82118  0.32351       NA       NA
#> 31 0031    2   DRM 0.81518  0.67339  0.07125       NA       NA
#> 32 0032    2   DRM 1.54708 -0.87791  0.25249       NA       NA
#> 33 0033    2   DRM 1.28305 -1.32022  0.16690       NA       NA
#> 34 0034    2   DRM 1.31798  0.17924  0.15299       NA       NA
#> 35 0035    2   DRM 1.48301 -0.13656  0.22323       NA       NA
#> 36 0036    2   DRM 0.89426  1.07673  0.12151       NA       NA
#> 37 0037    2   DRM 1.74214 -0.40933  0.09754       NA       NA
#> 38 0038    2   DRM 0.73472 -0.40711  0.20611       NA       NA
#> 39 0039    5   GRM 1.94931 -1.82855 -1.21040 -0.69548 -0.22021
#> 40 0040    5   GRM 1.30020 -0.70553 -0.06094  0.56671  1.06264
#> 41 0041    2   DRM 1.48353  0.63767  0.22449       NA       NA
#> 42 0042    2   DRM 1.79481 -1.50864  0.14541       NA       NA
#> 43 0043    2   DRM 1.45567  0.53570  0.13648       NA       NA
#> 44 0044    2   DRM 0.98254 -0.38343  0.11552       NA       NA
#> 45 0045    2   DRM 1.05231  2.31660  0.16867       NA       NA
#> 46 0046    2   DRM 2.41144  1.60547  0.18062       NA       NA
#> 47 0047    2   DRM 1.23217 -0.08089  0.12396       NA       NA
#> 48 0048    2   DRM 1.65413  0.16364  0.17983       NA       NA
#> 49 0049    2   DRM 1.21099  0.22634  0.07378       NA       NA
#> 50 0050    2   DRM 1.33178  1.32871  0.07789       NA       NA
#> 51 0051    2   DRM 1.76740 -0.98230  0.24200       NA       NA
#> 52 0052    2   DRM 0.97768 -0.75322  0.20821       NA       NA
#> 53 0053    5   GRM 1.15715 -0.36099  0.21902  0.84214  1.36791
#> 54 0054    5   GRM 1.25620 -2.03425 -1.31648 -0.68590 -0.10716
#> 55 0055    5   GRM 0.89213 -0.73461  0.00167  0.66112  1.23362
```
