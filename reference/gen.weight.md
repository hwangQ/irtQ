# Generate Weights

This function generates a set of normalized weights based on theta
(ability) values to be used in functions such as
[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md),
[`sx2_fit()`](https://hwangQ.github.io/irtQ/reference/sx2_fit.md), and
[`covirt()`](https://hwangQ.github.io/irtQ/reference/covirt.md).

## Usage

``` r
gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1, l = -4, u = 4, theta)
```

## Arguments

- n:

  An integer specifying the number of theta (node) values for which
  weights are to be generated. Default is 41.

- dist:

  A character string indicating the distribution type used to generate
  weights. Available options are `"norm"` for a normal distribution,
  `"unif"` for a uniform distribution, and `"emp"` for an empirical
  distribution.

  - If `dist = "norm"`, either `n` or `theta` must be provided.

  - If `dist = "unif"`, only `n` is applicable.

  - If `dist = "emp"`, only `theta` must be specified.

- mu, sigma:

  Mean and standard deviation of the normal distribution (used when
  `dist = "norm"`).

- l, u:

  Lower and upper bounds of the uniform distribution (used when
  `dist = "unif"`).

- theta:

  A numeric vector of empirical theta (node) values for which weights
  are generated.

## Value

A data frame with two columns:

- `theta`: The theta (node) values.

- `weight`: The corresponding normalized weights.

## Details

If `theta` is not specified, *n* equally spaced quadrature points and
corresponding weights are generated from either the normal or uniform
distribution:

- When `dist = "norm"`, Gaussian quadrature points and weights are
  computed using `gauss.quad.prob()` from the statmod package.

- When `dist = "unif"`, equally spaced points are drawn from the
  specified interval \[`l`, `u`\], and weights are proportional to the
  uniform density.

If `theta` is specified:

- When `dist = "norm"`, the weights are proportional to the normal
  density evaluated at each theta value and normalized to sum to 1.

- When `dist = "emp"`, equal weights are assigned to each provided theta
  value.

## See also

[`est_score()`](https://hwangQ.github.io/irtQ/reference/est_score.md),
[`sx2_fit()`](https://hwangQ.github.io/irtQ/reference/sx2_fit.md),
[`covirt()`](https://hwangQ.github.io/irtQ/reference/covirt.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example 1:
## Generate 41 Gaussian quadrature points and weights from the normal distribution
gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1)
#>            theta       weight
#> 1  -1.161494e+01 2.257864e-30
#> 2  -1.064754e+01 8.308559e-26
#> 3  -9.843433e+00 2.746891e-22
#> 4  -9.123070e+00 2.326384e-19
#> 5  -8.456099e+00 7.655982e-17
#> 6  -7.826882e+00 1.220335e-14
#> 7  -7.226023e+00 1.077818e-12
#> 8  -6.647308e+00 5.769853e-11
#> 9  -6.086349e+00 1.994795e-09
#> 10 -5.539884e+00 4.667348e-08
#> 11 -5.005397e+00 7.658186e-07
#> 12 -4.480878e+00 9.058609e-06
#> 13 -3.964684e+00 7.894719e-05
#> 14 -3.455432e+00 5.158014e-04
#> 15 -2.951937e+00 2.561642e-03
#> 16 -2.453159e+00 9.777903e-03
#> 17 -1.958171e+00 2.893721e-02
#> 18 -1.466125e+00 6.684766e-02
#> 19 -9.762388e-01 1.211489e-01
#> 20 -4.877686e-01 1.728495e-01
#> 21 -2.512499e-16 1.945450e-01
#> 22  4.877686e-01 1.728495e-01
#> 23  9.762388e-01 1.211489e-01
#> 24  1.466125e+00 6.684766e-02
#> 25  1.958171e+00 2.893721e-02
#> 26  2.453159e+00 9.777903e-03
#> 27  2.951937e+00 2.561642e-03
#> 28  3.455432e+00 5.158014e-04
#> 29  3.964684e+00 7.894719e-05
#> 30  4.480878e+00 9.058609e-06
#> 31  5.005397e+00 7.658186e-07
#> 32  5.539884e+00 4.667348e-08
#> 33  6.086349e+00 1.994795e-09
#> 34  6.647308e+00 5.769853e-11
#> 35  7.226023e+00 1.077818e-12
#> 36  7.826882e+00 1.220335e-14
#> 37  8.456099e+00 7.655982e-17
#> 38  9.123070e+00 2.326384e-19
#> 39  9.843433e+00 2.746891e-22
#> 40  1.064754e+01 8.308559e-26
#> 41  1.161494e+01 2.257864e-30

## Example 2:
## Generate 41 theta values and weights from the uniform distribution,
## given a minimum value of -4 and a maximum value of 4
gen.weight(n = 41, dist = "unif", l = -4, u = 4)
#>    theta     weight
#> 1   -4.0 0.02439024
#> 2   -3.8 0.02439024
#> 3   -3.6 0.02439024
#> 4   -3.4 0.02439024
#> 5   -3.2 0.02439024
#> 6   -3.0 0.02439024
#> 7   -2.8 0.02439024
#> 8   -2.6 0.02439024
#> 9   -2.4 0.02439024
#> 10  -2.2 0.02439024
#> 11  -2.0 0.02439024
#> 12  -1.8 0.02439024
#> 13  -1.6 0.02439024
#> 14  -1.4 0.02439024
#> 15  -1.2 0.02439024
#> 16  -1.0 0.02439024
#> 17  -0.8 0.02439024
#> 18  -0.6 0.02439024
#> 19  -0.4 0.02439024
#> 20  -0.2 0.02439024
#> 21   0.0 0.02439024
#> 22   0.2 0.02439024
#> 23   0.4 0.02439024
#> 24   0.6 0.02439024
#> 25   0.8 0.02439024
#> 26   1.0 0.02439024
#> 27   1.2 0.02439024
#> 28   1.4 0.02439024
#> 29   1.6 0.02439024
#> 30   1.8 0.02439024
#> 31   2.0 0.02439024
#> 32   2.2 0.02439024
#> 33   2.4 0.02439024
#> 34   2.6 0.02439024
#> 35   2.8 0.02439024
#> 36   3.0 0.02439024
#> 37   3.2 0.02439024
#> 38   3.4 0.02439024
#> 39   3.6 0.02439024
#> 40   3.8 0.02439024
#> 41   4.0 0.02439024

## Example 3:
## Generate normalized weights from the standard normal distribution,
## given a user-defined set of theta values
theta <- seq(-4, 4, by = 0.1)
gen.weight(dist = "norm", mu = 0, sigma = 1, theta = theta)
#>    theta       weight
#> 1   -4.0 1.338370e-05
#> 2   -3.9 1.986656e-05
#> 3   -3.8 2.919618e-05
#> 4   -3.7 4.248019e-05
#> 5   -3.6 6.119330e-05
#> 6   -3.5 8.727271e-05
#> 7   -3.4 1.232282e-04
#> 8   -3.3 1.722657e-04
#> 9   -3.2 2.384209e-04
#> 10  -3.1 3.266985e-04
#> 11  -3.0 4.432074e-04
#> 12  -2.9 5.952835e-04
#> 13  -2.8 7.915854e-04
#> 14  -2.7 1.042146e-03
#> 15  -2.6 1.358366e-03
#> 16  -2.5 1.752919e-03
#> 17  -2.4 2.239567e-03
#> 18  -2.3 2.832848e-03
#> 19  -2.2 3.547640e-03
#> 20  -2.1 4.398583e-03
#> 21  -2.0 5.399371e-03
#> 22  -1.9 6.561915e-03
#> 23  -1.8 7.895417e-03
#> 24  -1.7 9.405386e-03
#> 25  -1.6 1.109265e-02
#> 26  -1.5 1.295242e-02
#> 27  -1.4 1.497351e-02
#> 28  -1.3 1.713773e-02
#> 29  -1.2 1.941959e-02
#> 30  -1.1 2.178633e-02
#> 31  -1.0 2.419830e-02
#> 32  -0.9 2.660988e-02
#> 33  -0.8 2.897063e-02
#> 34  -0.7 3.122698e-02
#> 35  -0.6 3.332415e-02
#> 36  -0.5 3.520832e-02
#> 37  -0.4 3.682889e-02
#> 38  -0.3 3.814072e-02
#> 39  -0.2 3.910626e-02
#> 40  -0.1 3.969727e-02
#> 41   0.0 3.989626e-02
#> 42   0.1 3.969727e-02
#> 43   0.2 3.910626e-02
#> 44   0.3 3.814072e-02
#> 45   0.4 3.682889e-02
#> 46   0.5 3.520832e-02
#> 47   0.6 3.332415e-02
#> 48   0.7 3.122698e-02
#> 49   0.8 2.897063e-02
#> 50   0.9 2.660988e-02
#> 51   1.0 2.419830e-02
#> 52   1.1 2.178633e-02
#> 53   1.2 1.941959e-02
#> 54   1.3 1.713773e-02
#> 55   1.4 1.497351e-02
#> 56   1.5 1.295242e-02
#> 57   1.6 1.109265e-02
#> 58   1.7 9.405386e-03
#> 59   1.8 7.895417e-03
#> 60   1.9 6.561915e-03
#> 61   2.0 5.399371e-03
#> 62   2.1 4.398583e-03
#> 63   2.2 3.547640e-03
#> 64   2.3 2.832848e-03
#> 65   2.4 2.239567e-03
#> 66   2.5 1.752919e-03
#> 67   2.6 1.358366e-03
#> 68   2.7 1.042146e-03
#> 69   2.8 7.915854e-04
#> 70   2.9 5.952835e-04
#> 71   3.0 4.432074e-04
#> 72   3.1 3.266985e-04
#> 73   3.2 2.384209e-04
#> 74   3.3 1.722657e-04
#> 75   3.4 1.232282e-04
#> 76   3.5 8.727271e-05
#> 77   3.6 6.119330e-05
#> 78   3.7 4.248019e-05
#> 79   3.8 2.919618e-05
#> 80   3.9 1.986656e-05
#> 81   4.0 1.338370e-05

## Example 4:
## Generate equal normalized weights for theta values
## randomly sampled from the standard normal distribution
theta <- rnorm(100)
gen.weight(dist = "emp", theta = theta)
#>             theta weight
#> 1   -2.283046e+00   0.01
#> 2    1.652605e+00   0.01
#> 3    1.159808e+00   0.01
#> 4   -8.819287e-01   0.01
#> 5    9.452881e-01   0.01
#> 6    9.765185e-02   0.01
#> 7   -2.937563e+00   0.01
#> 8   -2.027510e-02   0.01
#> 9    3.900717e-01   0.01
#> 10  -1.942910e+00   0.01
#> 11  -6.381566e-01   0.01
#> 12  -3.076665e-01   0.01
#> 13   7.297180e-01   0.01
#> 14  -9.544751e-01   0.01
#> 15   8.580355e-01   0.01
#> 16  -1.481392e+00   0.01
#> 17  -9.303329e-01   0.01
#> 18  -1.016615e+00   0.01
#> 19  -2.116883e-01   0.01
#> 20   1.652866e+00   0.01
#> 21   7.252966e-01   0.01
#> 22  -2.322362e-01   0.01
#> 23   5.845080e-01   0.01
#> 24   2.415454e+00   0.01
#> 25   2.069770e-02   0.01
#> 26   3.373461e-01   0.01
#> 27   1.589362e-01   0.01
#> 28   9.256376e-01   0.01
#> 29   8.395591e-01   0.01
#> 30  -1.093925e+00   0.01
#> 31  -2.008624e-01   0.01
#> 32  -1.076350e+00   0.01
#> 33   2.599281e-01   0.01
#> 34  -7.896609e-01   0.01
#> 35  -2.002830e+00   0.01
#> 36   5.922348e-01   0.01
#> 37  -4.093031e-01   0.01
#> 38  -1.459266e+00   0.01
#> 39   1.501300e-01   0.01
#> 40  -7.663523e-01   0.01
#> 41   6.369605e-01   0.01
#> 42   7.179695e-02   0.01
#> 43   5.206760e-01   0.01
#> 44   7.046772e-01   0.01
#> 45   1.560625e+00   0.01
#> 46   1.113369e+00   0.01
#> 47  -5.690960e-01   0.01
#> 48  -9.182966e-01   0.01
#> 49   3.676878e-05   0.01
#> 50  -1.076635e-01   0.01
#> 51  -1.954848e+00   0.01
#> 52  -3.137719e-01   0.01
#> 53  -3.648914e-01   0.01
#> 54  -2.691230e-01   0.01
#> 55   5.202546e-01   0.01
#> 56  -4.154702e-01   0.01
#> 57  -1.004545e+00   0.01
#> 58   1.667266e-01   0.01
#> 59   2.820770e-01   0.01
#> 60   1.243522e+00   0.01
#> 61  -2.021780e-01   0.01
#> 62  -1.566144e+00   0.01
#> 63   2.234027e-01   0.01
#> 64  -1.260001e+00   0.01
#> 65   1.569565e+00   0.01
#> 66  -3.048343e-01   0.01
#> 67  -1.491114e+00   0.01
#> 68   2.326835e-01   0.01
#> 69   5.657022e-01   0.01
#> 70  -3.264018e-01   0.01
#> 71   9.234372e-01   0.01
#> 72  -1.140258e+00   0.01
#> 73  -7.381785e-01   0.01
#> 74  -1.055464e+00   0.01
#> 75   4.982216e-01   0.01
#> 76  -8.638460e-01   0.01
#> 77  -7.854393e-01   0.01
#> 78   4.949613e-01   0.01
#> 79  -1.552748e-01   0.01
#> 80   6.260879e-01   0.01
#> 81  -9.382887e-01   0.01
#> 82  -3.088562e-02   0.01
#> 83   2.190185e-01   0.01
#> 84   1.345403e-01   0.01
#> 85   8.457396e-01   0.01
#> 86  -2.089090e-01   0.01
#> 87   3.340083e-02   0.01
#> 88   8.123233e-01   0.01
#> 89   2.401086e+00   0.01
#> 90   6.459132e-01   0.01
#> 91   5.415901e-01   0.01
#> 92   7.025061e-01   0.01
#> 93  -1.138278e+00   0.01
#> 94   6.196648e-03   0.01
#> 95  -1.260213e+00   0.01
#> 96  -1.251115e+00   0.01
#> 97  -4.388664e-01   0.01
#> 98   1.697829e-01   0.01
#> 99  -6.257547e-01   0.01
#> 100  1.630853e+00   0.01
```
