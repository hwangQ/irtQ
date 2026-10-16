# Extract Panel Structure from an MST Route Map

Given a binary transition matrix that defines the Multistage-Adaptive
Test (MST) structure, this function extracts key panel information: the
module configuration by stage, all valid test pathways, the number of
modules per stage, and the total number of stages.

## Usage

``` r
panel_info(route_map)
```

## Arguments

- route_map:

  A binary square matrix defining the MST transition structure, where a
  value of 1 at row *i*, column *j* indicates that a test taker can be
  routed from module *i* to module *j*. This is equivalent to the
  `transMatrix` argument in the `randomMST()` function of the mstR
  package (Magis et al., 2017). See
  [`reval_mst`](https://hwangQ.github.io/irtQ/reference/reval_mst.md)
  for details on how to construct this matrix.

## Value

A named list with four elements:

- `config`:

  A named list of length equal to the number of stages. Each element is
  an integer vector of module indices belonging to that stage. Names are
  `"stage.1"`, `"stage.2"`, etc.

- `pathway`:

  An integer matrix of all valid test pathways, with one row per pathway
  and one column per stage. Row names are `"path.1"`, `"path.2"`, etc.,
  and column names are `"stage.1"`, `"stage.2"`, etc. Each entry is the
  index of the module administered at that stage.

- `n.module`:

  A named integer vector giving the number of modules in each stage.

- `n.stage`:

  An integer scalar giving the total number of stages.

## Details

The function identifies stage-1 modules as those with no incoming
transitions (column sum of `route_map` equals zero). It then traverses
the route map stage by stage, collecting modules reachable from the
previous stage. All logically possible pathways are enumerated via
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html), and pathways
that violate the transition constraints in `route_map` are removed. The
remaining pathways are sorted and returned as a matrix.

The route map must be a square matrix of 0s and 1s, and the modules do
not need to be numbered in the order of the stages. The function stops
with an error when the route map defines a single stage, contains a
cycle, has a module with no transition to or from another module, has a
module that cannot be reached from a stage-1 module or that belongs to
more than one stage, or has a module before the last stage without a
transition to a module of the next stage. A logical route map is read as
a 0/1 matrix.

## References

Magis, D., Yan, D., & von Davier, A. A. (2017). *Computerized adaptive
and multistage testing with R: Using packages catR and mstR*. Springer.
[doi:10.1007/978-3-319-69218-0](https://doi.org/10.1007/978-3-319-69218-0)
.

## See also

[`reval_mst`](https://hwangQ.github.io/irtQ/reference/reval_mst.md),
[`run_mst`](https://hwangQ.github.io/irtQ/reference/run_mst.md)

## Examples

``` r
# Use the 1-3-3 simMST panel route map
route_map <- simMST$route_map
pinfo <- panel_info(route_map)

# Stage configuration: which modules belong to which stage
pinfo$config
#> $stage.1
#> [1] 1
#> 
#> $stage.2
#> [1] 2 3 4
#> 
#> $stage.3
#> [1] 5 6 7
#> 

# All valid pathways (rows = pathways, columns = stages)
pinfo$pathway
#>        stage.1 stage.2 stage.3
#> path.1       1       2       5
#> path.2       1       2       6
#> path.3       1       3       5
#> path.4       1       3       6
#> path.5       1       3       7
#> path.6       1       4       6
#> path.7       1       4       7

# Number of modules per stage
pinfo$n.module
#> stage.1 stage.2 stage.3 
#>       1       3       3 

# Total number of stages
pinfo$n.stage
#> [1] 3
```
