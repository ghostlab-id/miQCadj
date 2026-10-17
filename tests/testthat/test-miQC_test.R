
# model slots
expected_slots <- c("flexmix_model", "flexmix_data", "model_info")
model_slots <- c("input_obj", "model_type", "column_names", "exclusion_info", "components_intercept", "components_median", "seed_used", "time_stamp", "flexmix_version")
column_slots <- c("percent_mt", "complexity_column", "model_column")
exclusion_slots <- c("exclude_invalid_values", "excluded_cells", "n_cells_total", "n_cells_valid" , "n_cells_invalid")
comp_intercept <- c("intercept_comp1", "intercept_comp2", "keep_comp", "discard_comp")
comp_median_general <- c("median_selection_requested", "median_selection_added")
comp_median_specific <- c("median_comp1", "median_comp2", "keep_median", "discard_median", "q_low_complexity", "n_low_complexity_cells", "low_complexity_cut", "low_complexity_flag")

flexmix_cols <- c("miQC.prob_invalid", "miQC.prob_valid", "miQC.keep_comp", "miQC.discard_comp")

# check components
plots_components <- c("classification","probability")

# define classification
# new model_info
comp_info <- c("component_selection", "posterior_cutoff")
# miQC.keep in flexmix_data

# adjustments
check_methods <- c("original_check", "geometry_check")

geom_cols <- list(original_check = c("miQC.keep_orig_adj", "miQC.orig_adj", "miQC.orig_region", "miQC.below_boundary", "miQC.left_cutoff"),
               geometry_check = c("miQC.keep_geom_adj", "miQC.geom_adj", "miQC.geom_region", "miQC.dist_up", "miQC.dist_low"))

# check predictions
plots_predictions_optional <- c("unadjusted_classification","prediction_comparison", "probability") 
plots_predictions_specific <- list(original_check= c("original_adjusted_classification", "original_region_adjustments"),
                               geometry_check= c("geometry_adjusted_classification", "geometry_region_adjustments"))

geom_add <- list(unadjusted = flexmix_cols,
                original_check = setdiff(geom_cols[["original_check"]], "miQC.keep_orig_adj"),
                geometry_check = setdiff(geom_cols[["geometry_check"]], "miQC.keep_geom_adj"))

final_slots <- c("flexmix_model", "model_info", "adjustment_method", "selected_classification")

# function to create a seurat object with a small test data set
create_test_data <- function(n_cells = 500, n_genes = 100, type = NULL) {
  
  type <-  match.arg(type, c("Seurat", "SingleCellExperiment"))
  set.seed(42)
  
  count_mtx <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes,
    ncol = n_cells)
  
  rownames(count_mtx) <- paste0("Gene", seq_len(n_genes))
  colnames(count_mtx) <- paste0("Cell", seq_len(n_cells))
  
  # create bimodal QC signal
  features <- c(
    rnorm(n_cells/2, mean = 100, sd = 5),
    rnorm(n_cells/2, mean = 50, sd = 5))
  
  counts <- round(features * runif(n_cells, 2, 5))
  
   pct_mt <- c(
    runif(n_cells/2, 0.01, 0.06),
    runif(n_cells/2, 0.12, 0.25))
  
  if(type == "Seurat"){
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(count_mtx))
  obj$nFeature_RNA <- features
  obj$nCount_RNA <- counts
  obj$percent.mt <- pct_mt
  }
  if (type == "SingleCellExperiment"){
  obj <- SingleCellExperiment::SingleCellExperiment(count_mtx)
  obj$detected <- features
  obj$counts <- counts
  obj$subsets_mito_percent <- pct_mt}
  
  obj}

######################################

# basic tests
test_that("main functions exist", {
  expect_true(is.function(run_miQCadj))
  expect_true(is.function(check_miQC_components))
  expect_true(is.function(define_classification))
  expect_true(is.function(adjust_classification))
  expect_true(is.function(original_adjustments))
  expect_true(is.function(geometry_adjustments))
  expect_true(is.function(check_miQC_predictions))
  expect_true(is.function(add_to_input_obj))})

######################################

# function output
test_that("Test 1: Default functions return expected objects", {
  
  for (obj_type in c("Seurat", "SingleCellExperiment")){
  # create data
  obj <- create_test_data(type = obj_type)
  
  # run miQCadj
  fm <- run_miQCadj(input_obj = obj)
  expect_type(fm, "list")
  # check list slots
  expect_true(all(expected_slots %in% names(fm)))
  # check model_info list
  fm_model <- fm[["flexmix_model"]]
  expect_s4_class(fm_model, "flexmix")
  # check flexmix_data columns
  fm_data <- fm[["flexmix_data"]]
  expect_true(all(flexmix_cols %in% colnames(fm_data)))
  # check model_info slots
  fm_info <- fm[["model_info"]]
  expect_type(fm_info, "list")
  expect_true(all(model_slots %in% names(fm_info)))
  expect_type(fm_info[["column_names"]], "list")
  expect_true(all(column_slots %in% names(fm_info[["column_names"]])))
  expect_type(fm_info[["exclusion_info"]], "list")
  expect_true(all(exclusion_slots %in% names(fm_info[["exclusion_info"]])))
  expect_type(fm_info[["components_intercept"]], "list")
  expect_true(all(comp_intercept %in% names(fm_info[["components_intercept"]])))
  expect_type(fm_info[["components_median"]], "list")
  expect_true(all(comp_median_general %in% names(fm_info[["components_median"]])))
  expect_false(any(comp_median_specific %in% names(fm_info[["components_median"]])))  
  expect_true(isFALSE(fm_info[["components_median"]][["median_selection_added"]]))
  
  # verify selected component identities
  intercept_info <- fm_info[["components_intercept"]]
  expect_true(intercept_info[["keep_comp"]] %in% 1:2)
  expect_true(intercept_info[["discard_comp"]] %in% 1:2)
  expect_false(identical(intercept_info[["keep_comp"]], intercept_info[["discard_comp"]]))
  
  # posterior probabilities should sum to 1
  prob_sum <- fm_data[["miQC.prob_invalid"]] + fm_data[["miQC.prob_valid"]]
  expect_equal(unname(prob_sum), rep(1, length(prob_sum)), tolerance = 1e-6)
  
  # check components
  plots <- check_miQC_components(fm, combine_plots = FALSE)
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(plots_components %in% names(plots)))
  
  # define classifications
  fm <- define_classification(fm)
  expect_type(fm, "list")
  expect_true(all(expected_slots %in% names(fm)))
  # check availability of new info
  fm_data <- fm[["flexmix_data"]]
  expect_true("miQC.keep" %in% colnames(fm_data))
  # check model_info list
  fm_info <- fm[["model_info"]]
  expect_type(fm_info, "list")
  expect_true(all(comp_info %in% names(fm_info)))
  expect_identical(fm_info[["component_selection"]], "intercept")
  
  # expected classifications
  expected_class <- ifelse(fm_data[["miQC.prob_invalid"]] <= 0.75, "keep", "discard")
  expect_equal(as.character(fm_data[["miQC.keep"]]), expected_class)
  
  # check geometry
  for (method in check_methods){
    fm <- adjust_classification(fm, adjustment_method = method)
    expect_type(fm, "list")
    expect_true(all(c(expected_slots, "adjustment_method") %in% names(fm)))
    expect_true(method %in% names(fm[["adjustment_method"]])) 
    # check availability of new info
    fm_data <- fm[["flexmix_data"]]
    expect_true(all(geom_cols[[method]] %in% colnames(fm_data)))
    
    # check predictions
    plots <- check_miQC_predictions(fm, adjustment_method = method, combine_plots = FALSE)
    expect_type(plots, "list")
    expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
    expect_length(plots, 2)
    expect_true(all(plots_predictions_specific[[method]] %in% names(plots)))
    expect_false(any(plots_predictions_optional %in% names(plots)))}

  # check availability of combined checks
  expect_true(all(c(geom_cols[["original_check"]], geom_cols[["geometry_check"]]) %in% colnames(fm_data)))
  
  # check predictions all
  plots <- check_miQC_predictions(fm, add_unadjusted_prediction = TRUE, add_compared_prediction = TRUE, add_probability = TRUE, combine_plots = FALSE)
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_length(plots, 7)
  expect_true(all(c(plots_predictions_optional,plots_predictions_specific[["original_check"]], plots_predictions_specific[["geometry_check"]]) %in% names(plots)))
  
  # re-adding to input object, final column
  for (method in c("unadjusted", check_methods)){
    updated_obj <- add_to_input_obj(flexmix_list = fm, input_obj = obj, adjustment_method = method, overwrite = TRUE, verbose=FALSE)
    expect_s4_class(updated_obj, obj_type)
    if (obj_type == "Seurat"){
    expect_true("miQC.final_prediction" %in% colnames(updated_obj@meta.data))
    expect_true("miQCadj" %in% names(updated_obj@misc))
    stored <- updated_obj@misc[["miQCadj"]]}
    if (obj_type =="SingleCellExperiment"){
      expect_true("miQC.final_prediction" %in% colnames(SummarizedExperiment::colData(updated_obj)))
      expect_true("miQCadj" %in% names(S4Vectors::metadata(updated_obj)))
      stored <- S4Vectors::metadata(updated_obj)[["miQCadj"]]}
    if (method == "unadjusted"){
      expect_true(all(setdiff(final_slots, c("adjustment_method", "flexmix_model")) %in% names(stored)))
    }else{
      expect_true(all(setdiff(final_slots, "flexmix_model") %in% names(stored)))}
    expect_false("flexmix_model" %in% names(stored))
    if (method != "unadjusted"){
     expect_true(method %in% names(stored[["adjustment_method"]]))
     expect_false(setdiff(c("original_check", "geometry_check"), method) %in% names(stored[["adjustment_method"]]))}}
  
  # only model, no info
  obj <- add_to_input_obj(flexmix_list = fm, input_obj = obj, adjustment_method = "unadjusted", include_model = TRUE, include_info = FALSE, verbose=FALSE)
  if (obj_type =="Seurat"){
  expect_true("miQCadj" %in% names(obj@misc))
  stored <- obj@misc[["miQCadj"]]}
  if (obj_type =="SingleCellExperiment"){
    expect_true("miQCadj" %in% names(S4Vectors::metadata(obj)))
    stored <- S4Vectors::metadata(obj)[["miQCadj"]]}
  expect_false(any(setdiff(final_slots, "flexmix_model") %in% names(stored)))
  expect_true("flexmix_model" %in% names(stored))
  # overwrite issue
  expect_error(add_to_input_obj(flexmix_list = fm, input_obj = obj, adjustment_method = "unadjusted", overwrite = FALSE, verbose=TRUE))}})

######################################


test_that("Test 2: Nonlinear median selection workflow", {
  
  obj <- create_test_data(type = "Seurat")
  
  median_cols <- paste0(flexmix_cols, "_median")
  
  for (model_type in c("spline", "polynomial", "polynomial_raw")) {
    
    # fit nonlinear mixture model
    fm <- run_miQCadj(input_obj = obj, model_type = model_type, 
                      add_median_selection = TRUE, use_seed = 123)
    expect_type(fm, "list")
    # check list slots
    expect_true(all(expected_slots %in% names(fm)))
    # check model_info list
    fm_model <- fm[["flexmix_model"]]
    expect_s4_class(fm_model, "flexmix")
    # check flexmix_data columns
    fm_data <- fm[["flexmix_data"]]
    expect_true(all(flexmix_cols %in% colnames(fm_data)) && all(median_cols %in% colnames(fm_data)))
    # check model_info list
    fm_info <- fm[["model_info"]]
    expect_type(fm_info, "list")
    expect_true(all(model_slots %in% names(fm_info)))
    expect_true(all(column_slots %in% names(fm_info[["column_names"]])))
    expect_true(all(exclusion_slots %in% names(fm_info[["exclusion_info"]])))
    expect_true(all(comp_intercept %in% names(fm_info[["components_intercept"]])))
    expect_true(all(comp_median_general %in% names(fm_info[["components_median"]])))
    expect_true(all(comp_median_specific %in% names(fm_info[["components_median"]]))) 
    
    # check median component information
    expect_type(fm_info[["components_median"]], "list")
    expect_true(isTRUE(fm_info[["components_median"]][["median_selection_added"]]))

    # verify selected component identities
    median_info <- fm_info[["components_median"]]
    expect_true(median_info[["keep_median"]] %in% 1:2)
    expect_true(median_info[["discard_median"]] %in% 1:2)
    expect_false(identical(median_info[["keep_median"]], median_info[["discard_median"]]))
    
    # posterior probabilities should sum to 1
    prob_sum <- fm_data[["miQC.prob_invalid_median"]] + fm_data[["miQC.prob_valid_median"]]
    expect_equal(unname(prob_sum), rep(1, length(prob_sum)), tolerance = 1e-6)
    
    # plot median-based components
    plots <- check_miQC_components(fm, use_median_columns = TRUE, combine_plots = FALSE)
    expect_true(all(plots_components %in% names(plots)))
    expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
    
    # classify using median-based components
    fm <- define_classification(fm, component_selection = "median", posterior_cutoff = 0.75)
    fm_data <- fm[["flexmix_data"]]
    expect_true("miQC.keep" %in% colnames(fm_data))
    expect_false(any(grepl("_median$", colnames(fm_data))))
    expect_identical(fm[["model_info"]][["component_selection"]], "median")
    
    # expected classifications
    expected_class <- ifelse(fm_data[["miQC.prob_invalid"]] <= 0.75, "keep", "discard")
    expect_equal(as.character(fm_data[["miQC.keep"]]), expected_class)
    
    # check geometry
    for (method in check_methods){
      fm <- adjust_classification(fm, adjustment_method = method)
      expect_type(fm, "list")
      expect_true(all(c(expected_slots, "adjustment_method") %in% names(fm)))
      expect_true(method %in% names(fm[["adjustment_method"]])) 
      # check availability of new info
      fm_data <- fm[["flexmix_data"]]
      expect_true(all(geom_cols[[method]] %in% colnames(fm_data)))
      
      # check predictions
      plots <- check_miQC_predictions(fm, adjustment_method = method, combine_plots = FALSE)
      expect_type(plots, "list")
      expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
      expect_length(plots, 2)
      expect_true(all(plots_predictions_specific[[method]] %in% names(plots)))
      expect_false(any(plots_predictions_optional %in% names(plots)))}
    
    # check availability of combined checks
    expect_true(all(c(geom_cols[["original_check"]], geom_cols[["geometry_check"]]) %in% colnames(fm_data)))

    # check predictions all
    plots <- check_miQC_predictions(fm,  add_unadjusted_prediction = TRUE, add_compared_prediction = TRUE, add_probability = TRUE, combine_plots = FALSE)
    expect_type(plots, "list")
    expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
    expect_length(plots, 7)
    expect_true(all(c(plots_predictions_optional, plots_predictions_specific[["original_check"]], plots_predictions_specific[["geometry_check"]]) %in% names(plots)))
    
    # re-adding to seurat, all columns
    for (method in c("unadjusted", check_methods)){
      updated_obj <- add_to_input_obj(flexmix_list = fm, input_obj = obj, transfer_columns = "all", adjustment_method = method, overwrite = TRUE, verbose=FALSE)
      expect_s4_class(updated_obj, "Seurat")
      expect_true(all(c("miQC.final_prediction", geom_add[["unadjusted"]], geom_add[["original_check"]], geom_add[["geometry_check"]]) %in% colnames(updated_obj@meta.data)))
      expect_true("miQCadj" %in% names(updated_obj@misc))}}})


######################################


test_that("Test 3: One-dimensional workflow", {
  
  # create data
  obj <- create_test_data(type = "Seurat")
  
  complex_col <- "nFeature_RNA"
  # run miQCadj
  fm <- run_miQCadj(input_obj = obj, model_type = "one_dimensional")
  expect_type(fm, "list")
  # check list slots
  expect_true(all(expected_slots %in% names(fm)))
  # check model_info list
  expect_s4_class(fm[["flexmix_model"]], "flexmix")
  # check flexmix_data columns
  expect_true(all(flexmix_cols %in% colnames(fm[["flexmix_data"]])))
  # check model_info list
  expect_type(fm[["model_info"]], "list")
  expect_equal(fm[["model_info"]][["column_names"]][["complexity_column"]], complex_col)
  
  # check components
  plots <- check_miQC_components(fm, combine_plots = FALSE)
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(plots_components %in% names(plots)))
  plots <- check_miQC_components(fm, combine_plots = FALSE, x_values = "complexity_column")
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(plots_components %in% names(plots)))
  
  # define classifications
  fm <- define_classification(fm)
  expect_type(fm, "list")
  expect_true(all(expected_slots %in% names(fm)))
  # check availability of new info
  expect_true("miQC.keep" %in% colnames(fm[["flexmix_data"]]))
  # check model_info list
  expect_type(fm[["model_info"]], "list")
  expect_true(all(comp_info %in% names(fm[["model_info"]])))
  
  # check geometry
  for (method in check_methods){
    fm <- adjust_classification(fm, adjustment_method = method)
    expect_type(fm, "list")
    expect_true(all(c(expected_slots, "adjustment_method") %in% names(fm)))
    expect_true(method %in% names(fm[["adjustment_method"]])) 
    # check availability of new info
    expect_true(all(geom_cols[[method]] %in% colnames(fm[["flexmix_data"]])))

    # check predictions
    plots <- check_miQC_predictions(fm, adjustment_method = method, combine_plots = FALSE)
    expect_type(plots, "list")
    expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
    expect_length(plots, 2)
    expect_true(all(plots_predictions_specific[[method]] %in% names(plots)))
    expect_false(any(plots_predictions_optional %in% names(plots)))}
  
  # check region
  expect_identical(unique(fm[["flexmix_data"]][["miQC.geom_region"]]), "valid")
  
  # check availability of combined checks
  expect_true(all(c(geom_cols[["original_check"]], geom_cols[["geometry_check"]]) %in% colnames(fm[["flexmix_data"]])))
  
  # check predictions all
  plots <- check_miQC_predictions(fm, combine_plots = FALSE)
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_length(plots, 4)
  expect_true(all(c(plots_predictions_specific[["original_check"]], plots_predictions_specific[["geometry_check"]]) %in% names(plots)))

  plots <- check_miQC_predictions(fm, combine_plots = FALSE, x_values = "complexity_column")
  expect_type(plots, "list")
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_length(plots, 4)
  expect_true(all(c(plots_predictions_specific[["original_check"]], plots_predictions_specific[["geometry_check"]]) %in% names(plots)))
  
  # re-adding to seurat
  for (method in c("unadjusted", check_methods)){
    updated_obj <- add_to_input_obj(flexmix_list = fm, input_obj = obj, adjustment_method = method, overwrite = TRUE, verbose=FALSE)
    expect_s4_class(updated_obj, "Seurat")
    expect_true("miQC.final_prediction" %in% colnames(updated_obj@meta.data))
    expect_true("miQCadj" %in% names(updated_obj@misc))}})


######################################

test_that("Predictions are matched by cell ID", {
  
  obj <- create_test_data(type = "Seurat")
  
  fm <- run_miQCadj(obj)
  fm <- define_classification(fm)
  
  # Reverse input cell order
  obj_reordered <- obj[, rev(colnames(obj))]
  
  updated_obj <- add_to_input_obj(flexmix_list = fm, input_obj = obj_reordered,
                                  adjustment_method = "unadjusted",
                                  transfer_columns = "final",
                                  include_info = FALSE, include_model = FALSE)
  
  expected <- as.character(fm[["flexmix_data"]][colnames(updated_obj), "miQC.keep"])
  
  expect_identical(as.character(updated_obj$miQC.final_prediction), expected)})

