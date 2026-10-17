#' Fit a two-component mixture model 
#'  
#' Fit a two-component mixture model for mitochondrial percentage versus a 
#' complexity metric using [flexmix](https://CRAN.R-project.org/package=flexmix). 
#'
#' @param input_obj Provided input object can either be a `Seurat` object or a
#'  `SummarizedExperiment`-derived object.
#' @param percent_mt  Character Scalar, optional. Metadata column with 
#'  mitochondrial percentage information. If `NULL` (default) it will be 
#'  resolved to `percent.mt`, if the incoming object is a Seurat object, and to 
#'  `subsets_mito_percent`, if the incoming object is a `SummarizedExperiment`-
#'  derived object. If the data uses distinct column names for the mitochondrial
#'  percentage information, it should be specified.
#' @param complexity_column Character Scalar, optional. Metadata column with 
#' feature information (default) used as predictor variable for fitting the 
#' `flexmix` model. If `NULL` (default) it will be resolved to `nFeature_RNA`, 
#'  if the incoming object is a Seurat object, and to `detected`, if the 
#'  provided object is a `SummarizedExperiment`-derived object. Should be 
#'  specified, if a distinct metadata column should be used, i.e. UMI counts, 
#'  or if distinct naming conventions apply to the default columns, i.e. the 
#'  feature column being called `nfeatures`. In one-dimensional models ignored
#'  for model fitting, but will be required for original post-hoc adjustments. 
#' @param model_type Character Scalar. `flexmix` model type used for fitting 
#'  the model. Available options are \code{"linear"}, \code{"one_dimensional"}, 
#'  \code{"spline"}, \code{"polynomial"}, or  \code{"polynomial_raw"}. Defaults
#'  to \code{"linear"}
#' @param exclude_invalid_values Logical Scalar. Whether cells with missing or
#'   non-finite values in `percent_mt` or `complexity_column` should be excluded
#'   from model fitting. Defaults to \code{TRUE}.
#' @param use_seed Positive Numeric, optional. Seed control for improving 
#'  reproducibility. Defaults to \code{123}. Set \code{NULL}, if no seed should 
#'  be used.
#' @param add_median_selection Logical Scalar. Whether component classification
#'  should be additionally performed based on aggregated low-quality 
#'  mitochondrial median values. Will only be performed, if model type is 
#'  `spline`, `polynomial`, or `polynomial_raw` and the number of cells used
#'  for fitting the model is at least 50. Defaults to \code{FALSE}.
#' 
#' @return A named list containing:
#' \itemize{
#'   \item{flexmix_model}{The model fitted using `flexmix`.}
#'   \item{flexmix_data}{`data.frame` with required meta data columns for 
#'   fitting the model, scores for fitted components, and posterior probabilities 
#'   for component classifications.}
#'   \item{model_info}{Information regarding the column names, settings, 
#'   component selection, and excluded cells when fitting the `flexmix` model.}
#' }
#'
#' @examples
#' \dontrun{
#' flexmix_model <- run_miQCadj(input_obj = obj)
#' }
#'
#' @export
run_miQCadj <- function (input_obj, percent_mt = NULL, complexity_column = NULL,
                         model_type = "linear", exclude_invalid_values=TRUE,
                         use_seed=123, add_median_selection = FALSE){
  
  #----------- check input ------------#
  
  n_min_model <- 20
  
  model_type <- match.arg(model_type, choices= c("one_dimensional", "linear", "spline", "polynomial", "polynomial_raw"))

  complexity_provided <- !is.null(complexity_column) 
  
  # check logical arguments
  logical_args <- list(exclude_invalid_values = exclude_invalid_values,
                       add_median_selection = add_median_selection)

  lapply(names(logical_args), function(x){
    check_x <- logical_args[[x]]
    if (!is.logical(check_x) ||
        length(check_x) != 1 ||
        is.na(check_x)) {
      stop("`",x,"` must be one non-missing logical value.")}})
  
  # set seed
  if (!is.null(use_seed)) {
    if (!is.numeric(use_seed) ||
        length(use_seed) != 1L ||
        is.na(use_seed) ||
        !is.finite(use_seed) ||
        use_seed < 0 ||
        use_seed != as.integer(use_seed)) {
      stop("`use_seed` must be a single non-negative integer or `NULL`.")}
    set.seed(use_seed)}

  #----------- prepare input data ------------#

  if (inherits(input_obj, "Seurat")) {
    
    if (!identical(colnames(input_obj), rownames(input_obj@meta.data))){
      stop("Column names don't match meta data rownames in the provided Seurat object.")}
    
    df <- input_obj@meta.data
    if (is.null(percent_mt)){
      percent_mt <- "percent.mt"}
    if (is.null(complexity_column)){
      complexity_column <- "nFeature_RNA"}
    
  }else if (inherits(input_obj, "SummarizedExperiment")) {
    df <- as.data.frame(SummarizedExperiment::colData(input_obj))
    rownames(df) <- colnames(input_obj)
    if (is.null(percent_mt)){
      percent_mt <- "subsets_mito_percent"}
    if (is.null(complexity_column)){
      complexity_column <- "detected"}
    
  }else{
    stop("`input_obj` must be a Seurat or SummarizedExperiment-derived object.")}

    percent_mt <- as.character(percent_mt)
    if (length(percent_mt) != 1L ||
        is.na(percent_mt) ||
        !nzchar(percent_mt)) {
      stop("`percent_mt` must be `NULL` or a non-empty character scalar.")}
    
    if (!is.null(complexity_column)){
    complexity_column <- as.character(complexity_column)
    if (length(complexity_column) != 1L ||
        is.na(complexity_column) ||
        !nzchar(complexity_column)) {
      stop("`complexity_column` must be `NULL` or a non-empty character scalar.")}}

    # Metadata columns to retain
    model_column <- if (model_type == "one_dimensional") NULL else complexity_column
    required_cols <- c(percent_mt = percent_mt)
    
    if (model_type == "one_dimensional"){
      if (!complexity_column %in% colnames(df)){
        if(complexity_provided){
        stop("Specified `complexity_column=",complexity_column,"` not found.")
          }else{
      message("`complexity_column=",complexity_column,"` not found. ",
        "Complexity-dependent post-hoc adjustments will not be available.")
      complexity_column <- NULL}}
      if (!is.null(complexity_column) && complexity_column %in% colnames(df)){
          message("`complexity_column=",complexity_column,"` will not be used as a model predictor, ",
                  "\nbut will be retained in `flexmix_data` and `model_info` for post-hoc adjustments.")}}
    
    if (!is.null(complexity_column)) {
      if (identical(percent_mt, complexity_column)) {
        stop("`percent_mt` and `complexity_column` must refer to different metadata columns.")}
      
      required_cols <- c(percent_mt = percent_mt, 
                         complexity_column = complexity_column)}
        
  # check columns
  lapply(names(required_cols), function(x){
    if (miQC_check_col(required_cols[[x]], colnames(df))){
      stop(required_cols[[x]], " column provided for ",x, " not found in metadata.")}})

  # get metadata
  n_cells_total <- nrow(df)
  md <- df[, unname(required_cols), drop = FALSE]
  colnames(md) <- names(required_cols)
  
  # ensure numeric columns
  for (col_id in names(required_cols)){
    md[[col_id]] <- suppressWarnings(as.numeric(as.character(md[[col_id]])))}
  
  md$valid_numeric <- rowSums(!is.finite(as.matrix(md[, names(required_cols), drop = FALSE]))) == 0L
  excl_invalid <- rownames(md)[!md$valid_numeric]
  
  if (length(excl_invalid) > 0L) {
    warning("`percent_mt` and `complexity_column` contain ", length(excl_invalid), " NA or non-finite values.")
    if (isFALSE(exclude_invalid_values)){
      stop("Modelling variables contain missing or non-finite values.")}}
  
  # exclude invalid cells
  if (length(excl_invalid) > 0L){
    md <- md[!rownames(md) %in% excl_invalid, , drop=FALSE]
    message("Excluding ", length(excl_invalid), " cells with missing or non-finite values from model fitting.")
  }else{
      md <- md}
  
  if (nrow(md) < n_min_model) {
    stop("Only ", nrow(md), " cells with valid numeric values for model fitting.
         \nFor reliable modelling `input_obj` should comprise at least ", n_min_model, " valid cells.")}
  
  #----------- fit model ------------#
  
  # run model
  miQC_model <- tryCatch({
    if (model_type == "linear") {
      miQC_model <- flexmix::flexmix(percent_mt ~ complexity_column, 
                                     data = md, k = 2)
    }else if (model_type == "one_dimensional") {
      miQC_model <- flexmix::flexmix(percent_mt ~ 1, data = md, k = 2)    
    }else if (model_type == "spline") {
      miQC_model <- flexmix::flexmix(percent_mt ~ splines::bs(complexity_column), 
                                     data = md, k = 2)
    }else if (model_type == "polynomial") {
      miQC_model<- flexmix::flexmix(percent_mt ~ poly(complexity_column, 
                                                      degree = 2), data = md, k = 2)
    }else if (model_type == "polynomial_raw") {
      md$feature_z <- as.numeric(scale(md$complexity_column))
      miQC_model<- flexmix::flexmix(percent_mt ~ feature_z + I(feature_z^2), data = md, k = 2)}},
    error = function(e) {
      stop("Error fitting ", model_type, " `flexmix` model: ", conditionMessage(e), 
           "\nConsider changing `model_type` or applying pre-filtering.", call. = FALSE)})
  
  #----------- select components ------------#
  
  # check model
  if (!inherits(miQC_model, "flexmix")) {
    stop("`flexmix` returned unexpected object of class ", class(miQC_model), ".")}

  if (ncol(flexmix::parameters(miQC_model)) != 2) {
    stop("Expected 2 columns for `flexmix` parameters data frame. 
         \nGot ", ncol(flexmix::parameters(miQC_model)), " column(s).")}
  
  # get components
  comp_df <- as.data.frame(flexmix::fitted(miQC_model))
  if (ncol(comp_df) != 2) {
    stop("Expected exactly 2 fitted mixture components.")}
  rownames(comp_df) <- rownames(md)
  
  # get posterior probabilities
  miQC_posterior <- as.data.frame(flexmix::posterior(miQC_model))
  if (ncol(miQC_posterior) != 2) {
    stop("Expected exactly 2 posterior probability columns.")}
  rownames(miQC_posterior) <- rownames(md)
  
  #----------- intercept ------------#
  
  # intercept selection
  intercept1 <- as.numeric(flexmix::parameters(miQC_model, component = 1)[1])
  intercept2 <- as.numeric(flexmix::parameters(miQC_model, component = 2)[1])
  
  if (intercept1 != intercept2){
    discard_comp <- if (intercept1 > intercept2) 1L else 2L
    keep_comp <- if (discard_comp == 1L) 2L else 1L
    df <- data.frame(miQC.keep_comp = comp_df[,keep_comp],
                     miQC.discard_comp = comp_df[,discard_comp],
                     miQC.prob_invalid = miQC_posterior[, discard_comp],
                     miQC.prob_valid = miQC_posterior[,keep_comp])
  }else{
    discard_comp <- NULL
    keep_comp <- NULL
      warning("Component identification not possible using intercepts, adding NA columns.")
    df <- data.frame(miQC.keep_comp = rep(NA_real_, nrow(md)),
                       miQC.discard_comp = rep(NA_real_, nrow(md)),
                       miQC.prob_invalid = rep(NA_real_, nrow(md)),
                       miQC.prob_valid = rep(NA_real_, nrow(md)))}
    
  rownames(df) <- rownames(md)
  md <- cbind(md, df)
    
  #----------- median ------------#
  
  # median selection
  median_components <- NULL
  median_selection_added <- FALSE
  if (isTRUE(add_median_selection)){
    if (model_type %in% c("linear", "one_dimensional") || nrow(md) < 50){
      stop("Median Selection should not be used for linear or one-dimensional models or datasets with less than 50 cells.")}

      n_min_low_complexity <- 10
      q_low_complexity <- 0.05
      
      # define low-complexity region
      low_complexity_flag <- FALSE
      low_complexity_cut <- as.numeric(stats::quantile(md$complexity_column, probs = q_low_complexity, na.rm = TRUE))
      idx_low_complexity <- md$complexity_column <= low_complexity_cut
      
      adj_low_complexity <- q_low_complexity
      while (sum(idx_low_complexity, na.rm = TRUE) < n_min_low_complexity && adj_low_complexity < 0.2) {
        adj_low_complexity <- min(adj_low_complexity * 2, 0.2)
        low_complexity_cut <- as.numeric(stats::quantile(md$complexity_column, probs = adj_low_complexity, na.rm = TRUE))
        idx_low_complexity <- md$complexity_column <= low_complexity_cut}
      
      if (sum(idx_low_complexity, na.rm = TRUE) < n_min_low_complexity){
        warning("Very low number of cells in the sample. Please check results.")
        low_complexity_flag <- TRUE}
        
      if (adj_low_complexity != q_low_complexity){
        message("For reliable median calculation at least ",n_min_low_complexity," cells should be available in the low-complexity region.
                \nAdjusted low-complexity quantile to ", adj_low_complexity, ".")}
    
      median_mt_comp1 <- stats::median(comp_df[,1][idx_low_complexity], na.rm = TRUE)
      median_mt_comp2 <- stats::median(comp_df[,2][idx_low_complexity], na.rm = TRUE)
      
      if (median_mt_comp1 == median_mt_comp2){
        discard_median <- NULL
        keep_median <- NULL
      }else{
        discard_median <- if (median_mt_comp1 > median_mt_comp2) 1L else 2L
        keep_median <- if (discard_median == 1L) 2L else 1L
        df_median <- data.frame(miQC.keep_comp = comp_df[,keep_median],
                         miQC.discard_comp = comp_df[,discard_median],
                         miQC.prob_invalid = miQC_posterior[, discard_median],
                         miQC.prob_valid = miQC_posterior[,keep_median])
        colnames(df_median) <- paste0(colnames(df_median), "_median")
        rownames(df_median) <- rownames(md)
        md <- cbind(md, df_median)
        median_selection_added <- TRUE}
    
      median_components <- list(median_comp1 = median_mt_comp1,
           median_comp2 = median_mt_comp2,
           keep_median = keep_median,
           discard_median = discard_median,
           q_low_complexity = adj_low_complexity,
           n_low_complexity_cells = sum(idx_low_complexity, na.rm=TRUE),
           low_complexity_cut = low_complexity_cut,
           low_complexity_flag = low_complexity_flag)}

  #----------- prepare output ------------#
  
  # save model info
  global_info <- list(input_obj = if (inherits(input_obj, "Seurat")) "Seurat" else "SummarizedExperiment",
                      model_type = model_type,
                      column_names = list(percent_mt = percent_mt,
                                          complexity_column = complexity_column,
                                          model_column = model_column),
                      exclusion_info = list(exclude_invalid_values = exclude_invalid_values,
                                            excluded_cells = excl_invalid,
                                            n_cells_total = n_cells_total,
                                            n_cells_valid = nrow(md),
                                            n_cells_invalid = length(excl_invalid)),
                      components_intercept = list(intercept_comp1 = intercept1,
                                                  intercept_comp2 = intercept2,
                                                  keep_comp = keep_comp,
                                                  discard_comp = discard_comp),
                      components_median = c(median_components,
                                            list(median_selection_requested = add_median_selection,
                                                 median_selection_added = median_selection_added)),
                      seed_used=use_seed,
                      time_stamp = Sys.time(),
                      flexmix_version = utils::packageVersion("flexmix"))
  
  # return model
  flexmix_model <- list("flexmix_model"= miQC_model,
                        "flexmix_data" = md,
                        "model_info"= global_info)
  
  return(flexmix_model)}

