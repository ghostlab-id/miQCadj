# Fit a two-component mixture model Fit a two-component mixture model for mitochondrial percentage versus a complexity metric using [flexmix](https://CRAN.R-project.org/package=flexmix).

Fit a two-component mixture model

Fit a two-component mixture model for mitochondrial percentage versus a
complexity metric using
[flexmix](https://CRAN.R-project.org/package=flexmix).

## Usage

``` r
run_miQCadj(
  input_obj,
  percent_mt = NULL,
  complexity_column = NULL,
  model_type = "linear",
  exclude_invalid_values = TRUE,
  use_seed = 123,
  add_median_selection = FALSE
)
```

## Arguments

- input_obj:

  Provided input object can either be a `Seurat` object or a
  `SummarizedExperiment`-derived object.

- percent_mt:

  Character Scalar, optional. Metadata column with mitochondrial
  percentage information. If `NULL` (default) it will be resolved to
  `percent.mt`, if the incoming object is a Seurat object, and to
  `subsets_mito_percent`, if the incoming object is a
  `SummarizedExperiment`- derived object. If the data uses distinct
  column names for the mitochondrial percentage information, it should
  be specified.

- complexity_column:

  Character Scalar, optional. Metadata column with feature information
  (default) used as predictor variable for fitting the `flexmix` model.
  If `NULL` (default) it will be resolved to `nFeature_RNA`, if the
  incoming object is a Seurat object, and to `detected`, if the provided
  object is a `SummarizedExperiment`-derived object. Should be
  specified, if a distinct metadata column should be used, i.e. UMI
  counts, or if distinct naming conventions apply to the default
  columns, i.e. the feature column being called `nfeatures`. In
  one-dimensional models ignored for model fitting, but will be required
  for original post-hoc adjustments.

- model_type:

  Character Scalar. `flexmix` model type used for fitting the model.
  Available options are `"linear"`, `"one_dimensional"`, `"spline"`,
  `"polynomial"`, or `"polynomial_raw"`. Defaults to `"linear"`

- exclude_invalid_values:

  Logical Scalar. Whether cells with missing or non-finite values in
  `percent_mt` or `complexity_column` should be excluded from model
  fitting. Defaults to `TRUE`.

- use_seed:

  Positive Numeric, optional. Seed control for improving
  reproducibility. Defaults to `123`. Set `NULL`, if no seed should be
  used.

- add_median_selection:

  Logical Scalar. Whether component classification should be
  additionally performed based on aggregated low-quality mitochondrial
  median values. Will only be performed, if model type is `spline`,
  `polynomial`, or `polynomial_raw` and the number of cells used for
  fitting the model is at least 50. Defaults to `FALSE`.

## Value

A named list containing:

- flexmix_modelThe model fitted using `flexmix`.

- flexmix_data`data.frame` with required meta data columns for fitting
  the model, scores for fitted components, and posterior probabilities
  for component classifications.

- model_infoInformation regarding the column names, settings, component
  selection, and excluded cells when fitting the `flexmix` model.

## Examples

``` r
if (FALSE) { # \dontrun{
flexmix_model <- run_miQCadj(input_obj = obj)
} # }
```
