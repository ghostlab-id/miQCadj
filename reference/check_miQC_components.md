# Visualize fitted miQC components and posterior classifications

Creates diagnostic scatter plots of mitochondrial percentage versus the
complexity metric used in
[`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md),
including keep/discard classifications based on a chosen posterior
cutoff, and optionally, fitted mixture component aggregation lines.

## Usage

``` r
check_miQC_components(
  flexmix_list,
  sample_name = NULL,
  use_median_columns = FALSE,
  x_values = "model_column",
  posterior_cutoff = 0.75,
  point_size = 1,
  point_alpha = 0.8,
  color_lines = c(discard = "black", keep = "blue"),
  color_points = c(discard = "dimgrey", keep = "cornflowerblue"),
  combine_plots = TRUE
)
```

## Arguments

- flexmix_list:

  The `list` object returned by
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md)
  with required `flexmix_data` and `model_info` slots.

- sample_name:

  Character scalar, optional. Sample name to be used for the plots. A
  prefix `paste0(sample_name, ": ")` will either be used for each plot,
  if `combine_plots=FALSE`, or for the global plot only, if
  `combine_plots=TRUE`. Defaults to `NULL`.

- use_median_columns:

  Logical Scalar. Whether component classification based on aggregated
  low-quality mitochondrial should be used. Only possible, if
  median-related component selection columns `"_median"` are available
  in the `flexmix_data` data frame. Defaults to `FALSE`.

- x_values:

  Character Scalar. Available options include `model_column` (default)
  or `complexity_column`. Will return identical plots for linear,
  spline, polynomial, and polynomial_raw models. For one-dimensional
  models `model_column` will return plots with decreasing mitochondrial
  percentage using ranked cell indices in the x-axis. If
  `complexity_column` had been provided to
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md)
  values from `complexity_column` can be used as `x_values`.

- posterior_cutoff:

  Numeric Scalar between 0 and 1. Cutoff used for keep/discard
  classification. Defaults to `0.75`.

- point_size:

  Non-negative Numeric Scalar. Point size used for scatter plots.
  Defaults to `1`.

- point_alpha:

  Non-negative Numeric Scalar. Transparency used for scatter plots.
  Defaults to `0.8`.

- color_lines:

  Named character vector, optional. Specifying the colors for `keep` and
  `discard` classifications, will plot lines for keep/discard
  classification curves based on aggregated medians. Defaults to
  `c(discard = "black", keep = "blue")`.

- color_points:

  Named character vector for keep/discard point color values with names
  `c("discard", "keep")`. Defaults to
  `c(discard = "dimgrey", keep = "cornflowerblue")`.

- combine_plots:

  Logical Scalar. Whether the list of plots should be combined into a
  2-column-grid using `patchwork`. Defaults to `TRUE`.

## Value

If `combine_plots=FALSE` named list of `ggplot` objects:

- classificationPlot with miQC classification based on posterior
  probability thresholds..

- probabilityPlot with miQC probability for each cell being related to
  low-quality classification.

If `combine_plots=TRUE` (default) the items will be combined into a
1x2-grid using patchwork with `sample_name` being used for the global
title, only.

## Examples

``` r
if (FALSE) { # \dontrun{
fm <- run_miQCadj(input_obj = obj)

plots <- check_miQC_components(
  flexmix_list = fm,
  sample_name = "Sample")

# show combined plot grid
print(plots)
} # }
```
