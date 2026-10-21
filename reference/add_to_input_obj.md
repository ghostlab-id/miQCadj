# Transfer selected miQC predictions to the input object

Transfers selected miQC classifications and associated diagnostic
information from `flexmix_list` to the original input object. Cells
without model predictions receive missing values. No cells are removed
from the input object.

## Usage

``` r
add_to_input_obj(
  flexmix_list,
  input_obj,
  adjustment_method,
  transfer_columns = "final",
  include_info = TRUE,
  include_model = FALSE,
  model_name = "miQCadj",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- flexmix_list:

  List returned by
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md)
  and optionally updated by
  [`adjust_classification()`](https://ghostlab-id.github.io/miQCadj/reference/adjust_classification.md).

- input_obj:

  Original Seurat or SummarizedExperiment-related object.

- adjustment_method:

  Character Scalar. Available options include `"unadjusted"`
  `"original_check"` or `"geometry_check"`.

- transfer_columns:

  Character scalar. Which columns should be added from `flexmix_data` to
  the metadata for `input_obj`. Available options include `"all"` adding
  all columns from `flexmix_data`, `"method_specific"` adding the
  method-specific columns from the corresponding `adjustment_method`
  from `flexmix_data`, or `"final"` adding only the final prediction
  column as selected with `adjustment_method` to `input_obj` metadata.
  Defaults to `"final"`.

- include_info:

  Logical scalar. Whether model information from `flexmix_list` should
  be added to `input_obj`. If a `Seurat` object was provided, the model
  information will be available in `input_obj@misc`. For a
  `SummarizedExperiment` object, the information will be available as
  metadata. Defaults to `TRUE`.

- include_model:

  Logical scalar. Whether the flexmix model from `flexmix_list` should
  be added to `input_obj`. If a Seurat object was provided, the model
  information will be available in `input_obj@misc`. For a
  SummarizedExperiment object, the information will be available as
  metadata. Defaults to `FALSE`.

- model_name:

  Character scalar. Name of model information slot. Required for adding
  model information to `input_obj`. Defaults to `"miQCadj"`.

- overwrite:

  Logical scalar. Whether existing destination metadata columns may be
  replaced. Defaults to `FALSE`.

- verbose:

  Logical scalar. Whether transfer information should be printed.

## Value

The input object with additional cell-level metadata. If requested, the
flexmix model and / or information are added to the input object.
