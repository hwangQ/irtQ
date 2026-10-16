# Run flexMIRT from Within R

This function runs flexMIRT (Cai, 2017) from within R by executing a
model specified in a flexMIRT syntax file (i.e., \*.flexmirt). flexMIRT
runs on Windows only, so this function works only on Windows. To use
this function, the flexMIRT software must be installed on your system.
This interface is especially useful for conducting simulation studies or
automating batch analyses involving flexMIRT.

## Usage

``` r
run_flexmirt(file.syntax, dir.flex = NULL, show.output.on.console = FALSE, ...)
```

## Arguments

- file.syntax:

  A single string or character vector specifying the path(s) to one or
  more flexMIRT syntax files (with extension \*.flexmirt) to be
  executed. For example: `"C:/Users/Data/irtmodel.flexmirt"`.

- dir.flex:

  A character string giving the folder where flexMIRT is installed
  (e.g., `"C:/Program Files/flexMIRT3"`). If `NULL` (default), the
  function looks in `"C:/Program Files"` for a folder whose name
  contains "flexMIRT"; give the folder explicitly when there is no such
  folder or more than one.

- show.output.on.console:

  Logical. If `TRUE`, the output of the system command is printed to the
  R console. Default is `FALSE`. See
  [`base::system()`](https://rdrr.io/r/base/system.html).

- ...:

  Additional arguments passed to
  [`base::system()`](https://rdrr.io/r/base/system.html).

## Value

Called for its side effect. flexMIRT writes its output files (e.g.,
"-prm.txt", "-irt.txt") to the folder of each syntax file, and they can
be read with
[`bring.flexmirt()`](https://hwangQ.github.io/irtQ/reference/bring.flexmirt.md).
The value returned to R is not meant to be used.

## Details

When using a version of flexMIRT earlier than 3.6, the directory
specified in `dir.flex` must contain the following six files:

- `WinFlexMIRT.exe`

- `FlexMIRT_x64.exe`

- `FlexMIRT_x86.exe`

- `vpg.dll`

- `vpg.licensing.client.dll`

- `vpg.licensing.dll`

For flexMIRT version 3.6 or later, the directory must include the
following five files:

- `WinFlexMIRT.exe`

- `vpg.dll`

- `vpg.licensing.client.dll`

- `vpg.licensing.dll`

- `VPGLicenseClientNet.dll`

along with a subdirectory named `Resources` that contains the following
two files:

- `flexMIRT_x64_AVX.exe`

- `flexMIRT_x86_AVX.exe`

The syntax files "2PLM_example.flexmirt" and "3PLM_example.flexmirt" in
the extdata folder of the package fit the 2PL and 3PL models to the
responses of 2,844 examinees to 12 items in "g341-19.dat". The folder
also contains the resulting "-prm.txt" and "-irt.txt" files.

## References

Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item
analysis and test scoring (Computer software). Chapel Hill, NC: Vector
Psychometric Group.

## Author

Hwanggyu Lim <hglim83@gmail.com>

## Examples

``` r
# Examples below will run if the flexMIRT software is installed
# in the default directory "C:/Program Files/flexMIRT3".
# Otherwise, specify the directory where flexMIRT is installed
# using the 'dir.flex' argument.

if (FALSE) { # \dontrun{
# (1) Run a single syntax file
# Copy an example syntax file and its data file to a temporary folder
ext <- system.file("extdata", package = "irtQ")
file.copy(file.path(ext, c("2PLM_example.flexmirt", "g341-19.dat")),
          tempdir())
file.syntax <- file.path(tempdir(), "2PLM_example.flexmirt")

# Run flexMIRT to estimate item parameters for the 2PL model
run_flexmirt(file.syntax = file.syntax, dir.flex = NULL,
             show.output.on.console = TRUE)

# Read the item parameter file written by flexMIRT
out.file <- file.path(tempdir(), "2PLM_example-prm.txt")
bring.flexmirt(out.file, type = "par")

# (2) Run multiple syntax files
# Copy two example syntax files and the data file to a temporary folder
file.copy(file.path(ext, c("3PLM_example.flexmirt", "g341-19.dat")),
          tempdir())
file.syntax1 <- file.path(tempdir(), "2PLM_example.flexmirt")
file.syntax2 <- file.path(tempdir(), "3PLM_example.flexmirt")

# Run flexMIRT to estimate item parameters for both models
run_flexmirt(file.syntax = c(file.syntax1, file.syntax2), dir.flex = NULL,
             show.output.on.console = FALSE)

# Read the item parameter files written by flexMIRT
out.file1 <- file.path(tempdir(), "2PLM_example-prm.txt")
out.file2 <- file.path(tempdir(), "3PLM_example-prm.txt")
bring.flexmirt(out.file1, type = "par")
bring.flexmirt(out.file2, type = "par")
} # }
```
