## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE)

## ----load-data----------------------------------------------------------------
library(miQCadj)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(ggplot2)})
data("pbmc_sub2", package = "miQCadj")

## ----table-sample-tech--------------------------------------------------------
stats::addmargins(table(pbmc_sub2$sample_tech))

## ----define-samples-----------------------------------------------------------
pbmc1 <- subset(pbmc_sub2, subset=sample_tech=="pbmc1_Smartseq2")
pbmc2 <- subset(pbmc_sub2, subset=sample_tech=="pbmc2_CELSeq2")
pbmc2 <- Seurat::as.SingleCellExperiment(pbmc2)

## ----run-pbmc1----------------------------------------------------------------
fm1 <- run_miQCadj(pbmc1)

## ----fm1-check-components, fig.height=5, fig.width=11, out.width="95%", fig.align="center", fig.cap="**pbmc1_Smartseq2: miQC Classification and Posterior Probability of a Cell being low Quality**"----
check_miQC_components(flexmix_list = fm1, sample_name = "pbmc1_Smartseq2")

## ----fm1-define-classification------------------------------------------------
fm1 <-  define_classification(flexmix_list = fm1)

## ----fm1-adjust-classification------------------------------------------------
fm1 <-  adjust_classification(flexmix_list = fm1)

## ----fm1-check-predictions, fig.height=5, fig.width=12, out.width="95%", fig.align="center", fig.cap="**pbmc1_Smartseq2: Unmodified and Adjusted Cell Predictions**"----
check_miQC_predictions(flexmix_list = fm1, sample_name = "pbmc1_Smartseq2",
                       add_unadjusted_prediction = TRUE,  
                       add_compared_prediction = TRUE,
                       adjustment_method = NULL)

## ----fm1-prediction-details, fig.height=8, fig.width=11, out.width="95%", fig.align="center", fig.cap="**pbmc1_Smartseq2: Original and Adjusted Cell Classification**"----
check_miQC_predictions(flexmix_list = fm1, sample_name = "pbmc1_Smartseq2")

## ----pbmc1-add-to-seurat------------------------------------------------------
pbmc1 <- add_to_input_obj(flexmix_list = fm1, input_obj = pbmc1, 
                          adjustment_method = "original_check")

## ----pbmc1-cell-counts--------------------------------------------------------
meta <- as.data.frame(pbmc1@meta.data)
stats::addmargins(table(meta$CellType, meta[["miQC.final_prediction"]]))

## ----pbmc1-cell-prop-counts---------------------------------------------------
round(prop.table(table(meta$CellType, meta[["miQC.final_prediction"]]), 
                 margin = 1), digits = 3)

## ----pbmc1-cell-type-plot, fig.height=5, fig.width=7, out.width="95%", fig.align="center", fig.cap="**pbmc1_Smartseq2: Cell Prediction status colored by Cell Type**"----
meta$CellType <- as.factor(meta$CellType)
meta[["miQC.final_prediction"]] <- as.factor(meta[["miQC.final_prediction"]])
ggplot2::ggplot(meta, ggplot2::aes(x = nFeature_RNA, 
                                   y = percent.mt, 
                                   colour = CellType, 
                                   shape = miQC.final_prediction)) +
  ggplot2::geom_point() + 
  ggplot2::scale_color_brewer(palette = "Paired") + 
  ggplot2::labs(
    title="pbmc1_Smartseq2: Cell Prediction Status colored by Cell Type") + 
  ggplot2::theme_bw()

## ----remaining-megakaryocytes-five-percent------------------------------------
# Calculate the 5th percentile of detected genes
feature_cutoff <- unname(quantile(meta$nFeature_RNA, probs = 0.05))
sum(meta$CellType == "Megakaryocyte" & meta$nFeature_RNA > feature_cutoff)

## ----fm2-polynomial, fig.height=5, fig.width=11, out.width="95%", fig.align="center", fig.cap="**Intercept component selection for pbmc2_CELSeq2: miQC Classification and Posterior Probabilities of a Cell being low Quality**"----
fm2 <- run_miQCadj(input_obj = pbmc2, percent_mt = "percent.mt", 
                   complexity_column = "nFeature_RNA", 
                   model_type = "polynomial", add_median_selection = TRUE)
check_miQC_components(flexmix_list = fm2, 
                      sample_name = "Intercept-based components for pbmc2_CELSeq2", 
                      combine_plots=TRUE)

## ----fm2-check-components-median, fig.height=5, fig.width=11, out.width="95%", fig.align="center", fig.cap="**Median component selection for pbmc2_CELSeq2: miQC Classification and Posterior Probabilities of a Cell being low Quality**"----
check_miQC_components(flexmix_list = fm2, 
                      sample_name = "Median-based components for pbmc2_CELSeq2", 
                      use_median_columns = TRUE)

## ----fm2-define-classification------------------------------------------------
fm2 <-  define_classification(flexmix_list = fm2, 
                              component_selection = "median")

## ----fm2-adjust-classification-check-prediction, fig.height=4.5, fig.width=12, out.width="95%", fig.align="center", fig.cap="**Unmodified and Adjusted Cell Predictions for pbmc2_CELSeq2**"----
fm2 <-  adjust_classification(flexmix_list = fm2)

# check predictions
check_miQC_predictions(flexmix_list = fm2,  sample_name = "pbmc2_CELSeq2",
                       add_unadjusted_prediction = TRUE, 
                       add_compared_prediction = TRUE,
                       adjustment_method = NULL)

## ----fm2-table-check-counts---------------------------------------------------
addmargins(table(fm2$flexmix_data[["miQC.keep_adj_compared"]]))

## ----fm2-adjust-prediction-no-left-cutoff, fig.height=4.5, fig.width=12, out.width="95%", fig.align="center", fig.cap="**Unmodified and Adjusted Cell Predictions for pbmc2_CELSeq2 (without left cutoff)**"----
# adjust classification
fm2 <-  adjust_classification(flexmix_list = fm2, enforce_left_cutoff = FALSE)

# check predictions
check_miQC_predictions(flexmix_list = fm2, sample_name = "pbmc2_CELSeq2",
                       add_unadjusted_prediction = TRUE, 
                       add_compared_prediction = TRUE,
                       adjustment_method = NULL)

## ----fm2-prediction-details,fig.height=8, fig.width=11, out.width="95%", fig.align="center", fig.cap="**Median-based Component Selection and Region-Adjustments for pbmc2_CELSeq2 (without left cutoff)**"----
check_miQC_predictions(fm2, sample_name = "pbmc2_CELSeq2")

## ----pbmc2-add-to-sce---------------------------------------------------------
pbmc2 <- add_to_input_obj(flexmix_list = fm2, input_obj = pbmc2,  
                          transfer_columns = "all", 
                          adjustment_method = "geometry_check")

## ----pbmc2-cell-counts--------------------------------------------------------
col_data <- as.data.frame(SummarizedExperiment::colData(pbmc2))
stats::addmargins(table(col_data$CellType, col_data[["miQC.final_prediction"]]))

## ----pbmc2-cell-prop-counts---------------------------------------------------
round(prop.table(table(col_data$CellType, col_data[["miQC.final_prediction"]]), 
                 margin = 1), digits = 3)

## ----pbmc2-cell-type-plot, fig.height=5, fig.width=7, out.width="95%", fig.align="center", fig.cap="**pbmc2_CELSeq2: Cell Prediction status colored by Cell Type**"----
col_data$CellType <- as.factor(col_data$CellType)
col_data[["miQC.final_prediction"]] <- as.factor(col_data[["miQC.final_prediction"]])
ggplot2::ggplot(col_data, ggplot2::aes(x = nFeature_RNA, 
                                       y = percent.mt, 
                                       colour = CellType, 
                                       shape = miQC.final_prediction)) +
  ggplot2::geom_point() + ggplot2::scale_color_brewer(palette = "Paired") + 
  ggplot2::labs(
    title="pbmc2_CELSeq2: Cell Prediction Status colored by Cell Type") +
  ggplot2::theme_bw()

## ----session-info-------------------------------------------------------------
session_info <- utils::sessionInfo()
session_info$loadedOnly <- NULL
print(session_info)

