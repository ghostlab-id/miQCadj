# Adjusted geometric validation of selected miQC classifications

Performs a post hoc geometry-based plausibility check of a selected miQC
classification. Cells may be reassigned when their current keep/discard
status is inconsistent with the fitted model curves.

## Usage

``` r
geometry_adjustments(
  df,
  model_type,
  rescue_quantile = 0.5,
  discard_quantile = 0.95,
  converge_quantile = 0.05,
  safe_keep_dist = 1.25,
  verbose = TRUE
)
```

## Arguments

- df:

  The `data.frame` object for `flexmix_data`.

- model_type:

  Character Scalar describing the model type used in
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md).
  Available options include `"one_dimensional"`, `"linear"`, `"spline"`,
  `"polynomial"`, or `"polynomial_raw"`.

- rescue_quantile:

  Numeric between 0 and 1. Lower threshold used for rescuing discarded
  cells. Will be ignored for one-dimensional model. Defaults to `0.5`.

- discard_quantile:

  Numeric between 0 and 1. Upper threshold used for reassigning
  implausible kept cells. Will be ignored for one-dimensional model.
  Defaults to `0.95`.

- converge_quantile:

  Numeric between 0 and 1. Lower quantile of curve gaps used to define
  convergence regions. Will be ignored for one-dimensional model.
  Defaults to `0.05`.

- safe_keep_dist:

  Non-negative Numeric Scalar. Distance ratio favoring keeping cells.
  Will be ignored for one-dimensional model. Defaults to `1.25`).

- verbose:

  Logical Scalar. Whether information regarding reassignments should be
  printed. Defaults to `TRUE`.

## Value

Named list with slot information for:

- df:

  Updated `data.frame` with geometry-adjusted `miQC.keep_geom_adj`
  classification column and geometry-related columns `miQC.geom_adj`,
  `miQC.geom_region`, `miQC.dist_up`, and `miQC.dist_low`.

- geometry_check:

  Adjustment thresholds used for geometry checks, including the
  specified function arguments.

## Details

The function adds adjusted classifications and diagnostic measures to
the provided `data.frame` and records thresholds and the number of
reassigned cells in `geometry_check`.

## Examples

``` r
if (FALSE) { # \dontrun{
df_list <- geometry_adjustments(df=flexmix_data, model_type = "linear")
} # }
```
