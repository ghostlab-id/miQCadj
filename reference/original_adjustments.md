# Original geometric validation of selected miQC classifications

Performs a post hoc geometric plausibility check of an existing miQC
classification. Cells may be reassigned based on the fitted
keep-component prediction and, optionally, a complexity-dependent left
cutoff.

## Usage

``` r
original_adjustments(
  df,
  keep_all_below_boundary = TRUE,
  min_index_type = "all",
  enforce_left_cutoff = TRUE,
  verbose = TRUE
)
```

## Arguments

- df:

  The `data.frame` object for `flexmix_data`.

- keep_all_below_boundary:

  Logical Scalar. Whether cells with mitochondrial percentage below the
  fitted keep component prediction should be retained regardless of
  their posterior probabilities. Defaults to `TRUE`.

- min_index_type:

  Character Scalar. Whether all cells should be checked for the minimum
  of mitochondrial percentage within discarded cells
  `min_index_type = "all"` or whether to restrict the check to discarded
  cells `min_index_type = "discard"`. Since the original wrapper used
  all cells, the default is `"all"`.

- enforce_left_cutoff:

  Logical Scalar. Whether cells with low library complexity and cells
  with high mitochondrial percentage should be discarded. Requires
  `complexity_column` within `flexmix_data` and at least one cell being
  discarded after applying `keep_all_below_boundary`. Defaults to
  `TRUE`.

- verbose:

  Logical Scalar. Whether information regarding reassignments should be
  printed. Defaults to `TRUE`.

## Value

Named list with slot information for:

- df:

  Updated `data.frame` with geometry-adjusted `miQC.keep_orig_adj`
  classification column and geometry-related columns
  `miQC.below_boundary`, `miQC.left_cutoff`, `miQC.orig_adj`,
  `miQC.orig_region`.

- original_check:

  Adjustment thresholds used for geometry checks, including the
  specified function arguments.

## Details

The function adds adjusted classifications and diagnostic measures to
the provided `data.frame` and records thresholds and the number of
reassigned cells in `original_check`.

## Examples

``` r
if (FALSE) { # \dontrun{
df_list <- original_adjustments(df=flexmix_data)
} # }
```
