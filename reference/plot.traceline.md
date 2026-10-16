# Plot Item and Test Characteristic Curves

This method visualizes item characteristic curves (ICCs), item score
curves, or the test characteristic curve (TCC) using the ggplot2
package. ICCs or item score curves can be plotted for one or more
selected items, while the TCC is plotted for the entire test form.

## Usage

``` r
# S3 method for class 'traceline'
plot(
  x,
  item.loc = NULL,
  score.curve = FALSE,
  overlap = FALSE,
  layout.col = 2,
  xlab.text,
  ylab.text,
  main.text,
  lab.size = 15,
  main.size = 15,
  axis.size = 15,
  line.color,
  line.size = 1,
  strip.size = 12,
  ...
)
```

## Arguments

- x:

  An object of class `traceline` obtained from
  [`traceline()`](https://hwangQ.github.io/irtQ/reference/traceline.md).

- item.loc:

  A numeric vector specifying the position(s) of the item(s) to plot. If
  `NULL` (default), the test characteristic curve (TCC) for the entire
  test form is plotted.

- score.curve:

  Logical. If `TRUE`, plots the item score curve, which is the expected
  item score (the sum of category scores weighted by their
  probabilities) as a function of theta. If `FALSE`, plots the item
  characteristic curves (ICCs) of all score categories. For dichotomous
  items, the item score curve equals the ICC for score category 1.
  Ignored when `item.loc = NULL`. Default is `FALSE`.

- overlap:

  Logical. Determines how multiple curves are displayed when plotting
  ICCs or item score curves. If `TRUE`, curves are overlaid in different
  colors: the category ICCs of each item share one panel per item when
  `score.curve = FALSE`, and the score curves of all items share a
  single panel when `score.curve = TRUE`. If `FALSE`, each curve is
  drawn in its own panel: one panel per score category of a single item
  when `score.curve = FALSE`, or one panel per item when
  `score.curve = TRUE`. Default is `FALSE`.

- layout.col:

  An integer value indicating the number of panel columns when multiple
  panels are drawn. Not used when `overlap = TRUE` and
  `score.curve = TRUE`. Default is 2.

- xlab.text, ylab.text:

  Character strings specifying the labels for the x and y axes,
  respectively. If missing, default labels are used.

- main.text:

  Character string specifying the overall title of the plot. If missing,
  a default title is used.

- lab.size:

  Numeric value specifying the font size of axis titles. Default is 15.

- main.size:

  Numeric value specifying the font size of the plot title. Default is
  15.

- axis.size:

  Numeric value specifying the font size of axis tick labels. Default is
  15.

- line.color:

  A character string specifying the color of the plot lines. Default is
  `"#F8766D"`. Not used when curves are overlaid (`overlap = TRUE` with
  `item.loc` specified), where colors are assigned automatically. See
  <http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/> for available
  color names.

- line.size:

  Numeric value specifying the thickness of plot lines. Default is 1.

- strip.size:

  Numeric. Font size of the panel (facet) labels. Default is 12.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  from the ggplot2 package.

## Value

This method prints the plot of item characteristic curves (ICCs), item
score curves, or the test characteristic curve (TCC), and invisibly
returns the ggplot2 object.

## Details

All plots are drawn with the ggplot2 package. If `item.loc = NULL`, the
test characteristic curve (TCC) of the entire test form is plotted.
Otherwise, `item.loc` gives the positions of the items to plot; for
example, `item.loc = 1:3` plots the first three items.

## See also

[`traceline()`](https://hwangQ.github.io/irtQ/reference/traceline.md)

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
## Example using a "-prm.txt" file exported from flexMIRT

# Import the "-prm.txt" output file from flexMIRT
flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Read the item parameters and convert them to item metadata
test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df

# Define a sequence of theta values
theta <- seq(-3, 3, 0.1)

# Compute item category probabilities and item/test characteristic functions
x <- traceline(x = test_flex, theta, D = 1)

# Plot the test characteristic curve (TCC) for the full test form
plot(x, item.loc = NULL)


# \donttest{
# Plot ICCs for the first item (dichotomous),
# with a separate panel for each score category
plot(x, item.loc = 1, score.curve = FALSE, layout.col = 2)


# Plot ICCs for the first item in a single panel
# (all score categories overlaid)
plot(x, item.loc = 1, score.curve = FALSE, overlap = TRUE)


# Plot ICCs for multiple items (both dichotomous and polytomous),
# with each item's ICCs shown in a single panel
plot(x, item.loc = c(1:3, 53:55), score.curve = FALSE, overlap = TRUE)


# Plot the item score curve for the first item (dichotomous)
plot(x, item.loc = 1, score.curve = TRUE)


# Plot item score curves for the first six dichotomous items
# using multiple panels
plot(x, item.loc = 1:6, score.curve = TRUE, overlap = FALSE)


# Plot item score curves for the first six dichotomous items
# overlaid in a single panel
plot(x, item.loc = 1:6, score.curve = TRUE, overlap = TRUE)


# Plot ICCs for the last item (polytomous),
# with each score category in a separate panel
plot(x, item.loc = 55, score.curve = FALSE, layout.col = 2)


# Plot the item score curve for the last item (polytomous)
plot(x, item.loc = 55, score.curve = TRUE)


# Plot item score curves for the last three polytomous items
# using multiple panels
plot(x, item.loc = 53:55, score.curve = TRUE, overlap = FALSE)


# Plot item score curves for the last three polytomous items
# overlaid in a single panel
plot(x, item.loc = 53:55, score.curve = TRUE, overlap = TRUE)

# }
```
