# Plot Item and Test Information Functions

This method plots item or test information functions for a specified set
of theta values. It can also display the conditional standard error of
estimation (CSEE) at the test level.

## Usage

``` r
# S3 method for class 'info'
plot(
  x,
  item.loc = NULL,
  overlap = FALSE,
  csee = FALSE,
  xlab.text,
  ylab.text,
  main.text,
  lab.size = 15,
  main.size = 15,
  axis.size = 15,
  line.color,
  line.size = 1,
  layout.col = 4,
  strip.size = 12,
  ...
)
```

## Arguments

- x:

  x An object of class `info` obtained from
  [`info()`](https://hwangQ.github.io/irtQ/reference/info.md).

- item.loc:

  A numeric vector indicating which item information functions to plot,
  specified by item position (e.g., 1 for the first item). If `NULL`
  (default), the test information function for the entire test form is
  plotted.

- overlap:

  Logical. If `TRUE`, multiple item information functions are plotted in
  a single panel. If `FALSE` (default), each item information function
  is displayed in a separate panel.

- csee:

  Logical. If `TRUE`, plots the conditional standard error of estimation
  (CSEE) at the test level. Note that the CSEE plot is only available at
  the test level, not for individual items. If `FALSE` (default), item
  or test information functions are plotted.

- xlab.text, ylab.text:

  Character strings specifying the labels for the x and y axes,
  respectively.

- main.text:

  Character string specifying the overall title of the plot.

- lab.size:

  Numeric value specifying the font size of axis titles. Default is 15.

- main.size:

  Numeric value specifying the font size of the plot title. Default is
  15.

- axis.size:

  Numeric value specifying the font size of axis tick labels. Default is
  15.

- line.color:

  A character string specifying the color of the plot lines. See
  <http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/> for available
  color names.

- line.size:

  Numeric value specifying the thickness of plot lines. Default is 1.

- layout.col:

  Integer. Number of columns to use when faceting multiple item
  information functions. Used only when `overlap = FALSE`. Default is 4.

- strip.size:

  Numeric value specifying the font size of facet labels when multiple
  items are displayed.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  from the ggplot2 package.

## Value

This method function displays the item or test information function
plot. When `csee = TRUE`, the CSEE is returned at the test level.

## Details

All of the plots are drawn using the ggplot2 package. The object of
class `info` can be obtained from the function
[`info()`](https://hwangQ.github.io/irtQ/reference/info.md).

## See also

[`info()`](https://hwangQ.github.io/irtQ/reference/info.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r

if (FALSE) { # \dontrun{
## Example using a "-prm.txt" file exported from flexMIRT

# Import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Read the item parameters and convert them to item metadata
test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df

# Define a sequence of theta values
theta <- seq(-4, 4, 0.1)

# Compute item and test information values for the given theta values
x <- info(x = test_flex, theta = theta, D = 1, tif = TRUE)

# Plot the test information function
plot(x)

# Plot the item information function for the second item
plot(x, item.loc = 2)

# Plot multiple item information functions, each in a separate panel
plot(x, item.loc = 1:8, overlap = FALSE)

# Plot multiple item information functions in a single panel
plot(x, item.loc = 1:8, overlap = TRUE)

# Plot the conditional standard error of estimation (CSEE) at the test level
plot(x, csee = TRUE)

} # }
```
