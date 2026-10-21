# Post-hoc Adjustments of cell classifications

The function adds adjusted classifications and diagnostic measures to
the `data.frame` stored in `flexmix_data` and adds thresholds and the
number of reassigned cells in `flexmix_list[["adjustment_method"]]`.

## Usage

``` r
adjust_classification(
  flexmix_list,
  adjustment_method = c("original_check", "geometry_check"),
  min_index_type = "all",
  enforce_left_cutoff = TRUE,
  verbose = TRUE
)
```

## Arguments

- flexmix_list:

  The `list` object returned by
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md)
  with required `flexmix_data` and `model_info` slots.

- adjustment_method:

  Character Scalar or Vector. Available options include `original_check`
  performing the original miQC filtering adjustments with
  `keep_all_below_boundary=TRUE` and `geometry_check` performing the
  modified geometry check using the default parameters with
  `rescue_quantile = 0.5`, `discard_quantile = 0.95`,
  `converge_quantile = 0.05`,`safe_keep_dist = 1.25`. For further
  adjustments the specific functions should be used. By default both
  checks will be performed.

- min_index_type:

  Character Scalar. The minimum of mitochondrial percentage within
  discarded cells. Whether all cells should be checked for this minimum
  `min_index_type = "all"` or whether to restrict the check to discarded
  cells `min_index_type = "discard"`. Since the original wrapper used
  all cells, the default is `"all"`.

- enforce_left_cutoff:

  Logical Scalar. Whether cells with low library complexity and cells
  with high mitochondrial percentage should be discarded. Will only be
  used, if `adjustment_method = original_check`. Defaults to `TRUE`.

- verbose:

  Logical Scalar. Whether information regarding reassignments should be
  printed. Defaults to `TRUE`.

## Value

`flexmix_list` with updated slot information for:

- flexmix_data:

  Updated data frame with adjustment results.

- adjustment_method:

  Method-specific adjustment information.

## See also

[`original_adjustments()`](https://ghostlab-id.github.io/miQCadj/reference/original_adjustments.md),
[`geometry_adjustments()`](https://ghostlab-id.github.io/miQCadj/reference/geometry_adjustments.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# run both geometry checks
flexmix_list <- adjust_classification(flexmix_list)
} # }
```
