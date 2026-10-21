# Define components and posterior cutoff for cell classification

Define components and posterior cutoff for cell classification

## Usage

``` r
define_classification(
  flexmix_list,
  component_selection = "intercept",
  posterior_cutoff = 0.75
)
```

## Arguments

- flexmix_list:

  The `list` object returned by
  [`run_miQCadj()`](https://ghostlab-id.github.io/miQCadj/reference/run_miQCadj.md)
  with required `flexmix_data` and `model_info` slots.

- component_selection:

  Character scalar. Available options are `median` or `intercept`
  (default) resulting in the corresponding component selection.

- posterior_cutoff:

  Numeric between 0 and 1. Defaults to `0.75`.

## Value

`flexmix_list` with updated slot information for:

- `flexmix_data`Updated `data.frame`: If `comp_selection == intercept`
  median-related columns are discarded. If `comp_selection == median`
  the median-related columns replace the intercept-derived
  `miQC.prob_invalid`, `miQC.prob_valid`, `miQC.keep_comp`, and
  `miQC.discard_comp` columns. A new column `miQC.keep` is added,
  comprising cell classifications based on \# posterior thresholds and
  the updated `miQC.prob_invalid` column.

- `model_info`Updated list with `component_selection` and
  `posterior_cutoff` as specified in the function arguments.

## Examples

``` r
if (FALSE) { # \dontrun{
flexmix_list <- define_classification(
  flexmix_list,
  component_selection = "median",
  posterior_cutoff = 0.75
)
} # }
```
