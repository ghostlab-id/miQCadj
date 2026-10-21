
# miQCadj

`miQCadj` provides extensions and utilities for
[miQC](https://bioconductor.org/packages//release/bioc/html/miQC.html)
quality-control cell filtering in
[Seurat](https://satijalab.org/seurat/) and
[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)
workflows.

The original `miQC` package uses a two-component mixture model fitted
with [flexmix](https://CRAN.R-project.org/package=flexmix) to model
mitochondrial transcript percentage as a function of cellular complexity
and identify potentially low-quality cells in single-cell RNA-sequencing
datasets.

`miQCadj` extends this workflow with additional model specifications,
configurable component selection, posterior-probability-based
classification, and post-hoc geometric validation. It provides both an
implementation of the original `miQC` classification adjustments and an
alternative geometry-based approach, allowing users to inspect and
compare adjusted cell classifications before transferring the selected
results back to their input objects.

------------------------------------------------------------------------

## Installation

You can install `miQCadj` from
[GitHub](https://github.com/ghostlab-id/miQCadj/) using

``` r
# install.packages("remotes")
remotes::install_github("ghostlab-id/miQCadj")
```

Alternatively, `devtools` can be used for installation.

------------------------------------------------------------------------

## Exported Functions

Exported functions comprise (alphabetically ordered):

| Function                   | Purpose                                                        |
|:---------------------------|:---------------------------------------------------------------|
| `add_to_input_obj()`       | Add results to the Input object                                |
| `adjust_classification()`  | Post-hoc Adjustments of cell classifications                   |
| `check_miQC_components()`  | Visualize fitted miQC components and posterior classifications |
| `check_miQC_predictions()` | Visualize cell classifications after geometry validation       |
| `define_classification()`  | Define components and posterior cutoff for cell classification |
| `geometry_adjustments()`   | Adjusted geometric validation of selected miQC classifications |
| `original_adjustments()`   | Original geometric validation of selected miQC classifications |
| `run_miQCadj()`            | Fit a two-component mixture model                              |

------------------------------------------------------------------------

## General Workflow

``` mermaid
flowchart TD
  A("run_miQCadj()") -->B("check_miQC_components()")
  A -->C("define_classification()")
  B -->C("define_classification()")
  C -->D("adjust_classification()")
  C -->E("check_miQC_predictions()")
  D -->E("check_miQC_predictions()")
  E -->F("add_to_input_obj()")

style A color:#000000,fill:#ffffff,stroke:#808080
style B color:#000000,fill:#8eb6ed,stroke:#1e4985
style C color:#000000,fill:#ffffff,stroke:#808080
style D color:#000000,fill:#ffffff,stroke:#808080
style E color:#000000,fill:#8eb6ed,stroke:#1e4985
style F color:#000000,fill:#ffffff,stroke:#808080
```

------------------------------------------------------------------------

## Basic workflow

Load package and fit a flexmix mixture model:

``` r
library(miQCadj)
flexmix_list <- run_miQCadj(sc_data)
```

By default a linear model (model_type = “linear”) facilitating
intercept-based component selection will be used.

------------------------------------------------------------------------

Visualize component identities and posterior probabilities of cells
belonging to the low-quality component:

``` r
plots <- check_miQC_components(flexmix_list)
```

By default the `ggplot` objects are combined into a 1x2 plot grid. To
return a list of `ggplot` objects set `combine_plots = FALSE`.

------------------------------------------------------------------------

Define component selection mode and posterior probability threshold
filtering:

``` r
flexmix_list <- define_classification(flexmix_list)
```

By default the component selection using intercepts
(`component_selection = "intercept"`) and a posterior probability
threshold of `posterior_cutoff=0.75` will be used.

------------------------------------------------------------------------

Adjust predictions using geometric plausibility checks:

``` r
flexmix_list <- adjust_classification(flexmix_list)
```

By default post-hoc adjustments as implemented in the original `miQC`
package as well as adjusted geometric adjustments are performed
(`adjustment_method = c("original_check", "geometry_check")`). Both
adjustment methods can be used independently, but should be added by
`adjust_classification()`, if requested in downstream analysis.

------------------------------------------------------------------------

Inspect adjusted predictions:

``` r
plots <- check_miQC_predictions(flexmix_list)
```

By default plots for both post-hoc adjustment methods
(`adjustment_method = c("original_check", "geometry_check")`) are
created and combined into a 2x2 patchwork grid. The unmodified
prediction plots can be added, setting
`add_unadjusted_prediction = TRUE`.

Select final predictions and add data back to the input object:

``` r
sc_data <- add_to_input_obj(flexmix_list = flexmix_list, 
input_obj = sc_data, 
adjustment_method="geometry_check")
```

The final prediction based on the selected flexmix model, component
mode, posterior cutoff, and adjustment method will be available in the
input object meta data columns, enabling further validation or
filtering.

By default the related information will also be saved to the input
object, and will be available in the `misc` slot for Seurat objects or
in the `metadata` slot for SummarizedExperiment objects.

------------------------------------------------------------------------

## Vignettes

For further information please check the package vignette:
[miQC_pbmcsca](https://ghostlab-id.github.io/miQCadj/doc/miQC_pbmcsca/miQC_pbmcsca.html)

------------------------------------------------------------------------

## Notes

miQCadj is designed for exploratory and reproducible QC workflows.
Though, the predictions made by miQCadj should never be trusted without
being checked.

------------------------------------------------------------------------

## Development status

Initial public release: **v0.1.0**

------------------------------------------------------------------------

## Citation

If you use `miQCadj` in published work, please cite the package as well
as
[`miQC`](https://bioconductor.org/packages//release/bioc/html/miQC.html),
and [flexmix](https://CRAN.R-project.org/package=flexmix) packages.

Package citation information is available in R using:

``` r
citation("miQCadj")
citation("miQC")
citation("flexmix")
```

Alternatively, the citation for `miQCadj`, `miQC`, and `flexmix` can be
copied below:

**miQCadj**

    Spies E (2026). "miQCadj: Extensions and Utilities for miQC-based Cell Filtering." R package version 0.1.0.

**miQC**

    Hippen A, Falco M, Weber L, Erkan E, Zhang K, Doherty J, Vähärautio A, Greene C, Hicks S (2021). “miQC: An adaptive probabilistic framework for quality control of single-cell RNA-sequencing data.” PLoS Comput Biol. https://doi.org/10.1371/journal.pcbi.1009290.

**flexmix**

    Gruen B, Leisch F (2026). flexmix: Flexible Mixture Modeling. doi:10.32614/CRAN.package.flexmix. R package version 2.3-21, https://CRAN.R-project.org/package=flexmix.
