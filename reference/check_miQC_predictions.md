# Visualize miQC classifications after geometry validation

Creates diagnostic scatter plots for miQC predictions after
[`adjust_classification()`](https://ghostlab-id.github.io/miQCadj/reference/adjust_classification.md).

## Usage

``` r
check_miQC_predictions(
  flexmix_list,
  sample_name = NULL,
  add_unadjusted_prediction = FALSE,
  add_compared_prediction = FALSE,
  add_probability = FALSE,
  adjustment_method = c("original_check", "geometry_check"),
  x_values = "model_column",
  point_size = 1,
  point_alpha = 0.8,
  color_lines = c(discard = "black", keep = "blue"),
  color_points = c(discard = "dimgrey", keep = "cornflowerblue"),
  shape_geom = c(unchanged = 1, geom_pass = 15, geom_fail = 4),
  color_regions = list(geometry_check = c(valid = "#99D696", converge = "#F2DB68",
    invalid = "#D6A696"), original_check = c(unchanged = "#99D696", below_boundary =
    "#F2DB68", left_cutoff = "#D6A696")),
  color_compare = c(both_keep = "#6495ED", both_discard = "#A9A9A9", keep_orig =
    "#6932C2", keep_geom = "#123EA3"),
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

- add_unadjusted_prediction:

  Logical Scalar. Whether the unadjusted prediction plot for
  `"miQC.keep"`, should be re-created. Defaults to `FALSE`.

- add_compared_prediction:

  Logical Scalar. Whether the plot comparing the post-hoc prediction
  adjustments `"miQC.keep_adj_compared"` created by
  [`adjust_classification()`](https://ghostlab-id.github.io/miQCadj/reference/adjust_classification.md),
  if both adjustment methods were applied, should should be created.
  Defaults to `FALSE`.

- add_probability:

  Logical Scalar. Whether the plot with the probability for each cell
  belonging to the compromised compartment `"miQC.prob_invalid"`, should
  be re-created. Defaults to `FALSE`.

- adjustment_method:

  Character Scalar or Vector, optional. Name of method(s) used for
  geometric prediction adjustments or `NULL`. Available options include
  `original_check` with required columns `miQC.keep_orig_adj`,
  `miQC.orig_adj`, and `miQC.orig_region` and `geometry_check` with
  `miQC.keep_geom_adj`, `miQC.geom_adj`, and `miQC.geom_region` in
  `flexmix_data`. By default plots for both methods will be created.

- x_values:

  Character Scalar. Available options include `model_column` (default)
  or `complexity_column`. Will return identical plots for linear,
  spline, polynomial, and polynomial_raw models. For one-dimensional
  models `model_column` will return plots with decreasing mitochondrial
  percentage using ranked cell indices in the x-axis. If
  `complexity_column` had been provided to
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md)
  values from `complexity_column` can be used as `x_values`.

- point_size:

  Positive Numeric. Specification of the point size. Defaults to `1`.

- point_alpha:

  Positive Numeric. Specification of the transparency of the points.
  Defaults to `1`.

- color_lines:

  Named character vector, optional. Specifying the colors for `keep` and
  `discard` classifications, will plot lines for keep/discard
  classification curves based on aggregated medians. Defaults to
  `c(discard = "black", keep = "blue")`.

- color_points:

  Named character vector for keep/discard point color values with names
  `c("discard", "keep")`. Defaults to
  `c(discard = "dimgrey", keep = "cornflowerblue")`

- shape_geom:

  Named character vector specifying the shapes for geometry-based
  adjustments. Defaults to
  `c(unchanged = 1, geom_pass = 15, geom_fail = 4)`.

- color_regions:

  Named character vector specifying the colors for geometry-based
  regions. Defaults to
  `c(valid="#81F168", invalid="#E28861", converge="gold")`.

- color_compare:

  Named character vector specifying the colors for comparing predictions
  by `original_check` and `geometry_check`. Defaults to
  `c(both_keep = "#6495ED", both_discard = "#A9A9A9", keep_orig= "#6932C2", keep_geom="#123EA3")`

- combine_plots:

  Logical Scalar. Whether the list of plots should be combined into a
  2-column-grid using `patchwork`. Defaults to `TRUE`.

## Value

If `combine_plots=FALSE` named list of `ggplot` objects. Depending on
the requested adjustment_method, and the additional plots, it may
include:

- `unadjusted_classification`Original keep/discard classifications.

- `prediction_comparison`Comparing keep/discard classifications from
  `original_check` and `geometry_check`. Can only be created, if both
  adjustment methods were applied, though the plot can be created
  independently.

- `probability`Posterior invalid probabilities.

- `original_adjusted_classification`Adjusted classifications using
  original prediction adjustments.

- `original_region_adjustments`Region-related adjustments performed by
  original prediction adjustments.

- `geometry_adjusted_classification`Adjusted classifications using
  geometruc prediction adjustments.

- `geometry_region_adjustments`Region-related adjustments performed by
  geometric prediction adjustments.

If `combine_plots=TRUE` (default) the items will be combined into a
2-column grid using patchwork with `sample_name` being used for the
global title, only.

## Examples

``` r
if (FALSE) { # \dontrun{
plots <- check_miQC_predictions(flexmix_list, sample_name = "Sample")
# show combined plot grid
print(plots)
} # }
```
