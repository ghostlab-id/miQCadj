#' Transfer selected miQC predictions to the input object
#'
#' Transfers selected miQC classifications and associated diagnostic
#' information from `flexmix_list` to the original input object.
#' Cells without model predictions receive missing values.
#' No cells are removed from the input object.
#'
#' @param flexmix_list List returned by `run_miQCadj()` and optionally updated 
#'  by `adjust_classification()`.
#' @param input_obj Original Seurat or SummarizedExperiment-related object.
#' @param adjustment_method Character Scalar. Available options 
#'  include `"unadjusted"` `"original_check"` or `"geometry_check"`.
#' @param transfer_columns Character scalar. Which columns should be added
#'  from `flexmix_data` to the metadata for `input_obj`. Available options 
#'  include \code{"all"} adding all columns from `flexmix_data`,
#'  \code{"method_specific"} adding the method-specific columns from the 
#'  corresponding `adjustment_method` from `flexmix_data`, or \code{"final"} 
#'  adding only the final prediction column as selected with `adjustment_method`
#'  to `input_obj` metadata. Defaults to \code{"final"}.
#' @param include_info Logical scalar. Whether model information from 
#' `flexmix_list` should be added to `input_obj`. If a `Seurat` object was 
#'  provided, the model information will be available in `input_obj@misc`. For 
#'  a `SummarizedExperiment` object, the information will be available as 
#'  metadata. Defaults to `TRUE`.
#' @param include_model Logical scalar. Whether the flexmix model from 
#' `flexmix_list` should be added to `input_obj`. If a Seurat object was 
#'  provided, the model information will be available in `input_obj@misc`. For 
#'  a SummarizedExperiment object, the information will be available as 
#'  metadata. Defaults to `FALSE`.
#' @param model_name Character scalar. Name of model information slot.
#'  Required for adding model information to `input_obj`. Defaults to 
#'  \code{"miQCadj"}.
#' @param overwrite Logical scalar. Whether existing destination
#'   metadata columns may be replaced. Defaults to `FALSE`.
#' @param verbose Logical scalar. Whether transfer information
#'   should be printed.
#'
#' @return The input object with additional cell-level metadata. If requested, 
#'  the flexmix model and / or information are added to the input object.
#'
#' @export
add_to_input_obj <- function(flexmix_list, input_obj, adjustment_method,
                             transfer_columns="final", 
                             include_info = TRUE,
                             include_model = FALSE,
                             model_name="miQCadj",
                             overwrite = FALSE, verbose =TRUE) {
  
  # check input
  adjustment_method <- match.arg(adjustment_method, choices= c("unadjusted", "original_check", "geometry_check"))
  transfer_columns <- match.arg(transfer_columns, choices= c("all", "method_specific", "final"))
  
  # logical arguments
  logical_args <- list(include_info = include_info,
                       include_model = include_model,
                       overwrite = overwrite,
                       verbose=verbose)
  
  lapply(names(logical_args), function(x){
    check_x <- logical_args[[x]]
    if (!is.logical(check_x) ||
        length(check_x) != 1 ||
        is.na(check_x)) {
      stop("`",x,"` must be one non-missing logical value.")}})
  
  # input_obj
  if (inherits(input_obj, "Seurat")) {
    input_type <- "Seurat"
    if (!identical(colnames(input_obj), rownames(input_obj@meta.data))){
      stop("Column names don't match meta data rownames in the provided Seurat object.")}
  }else if (inherits(input_obj, "SummarizedExperiment")) {
    input_type <- "SummarizedExperiment"
  }else{
    stop("`input_obj` must be a Seurat or SummarizedExperiment object.")}
  
  # flexmix_list
  if (!inherits(flexmix_list, "list")) {
    stop("`flexmix_list` must be a list object.")}
  
  # model & data
  required_slots <- c("flexmix_data",
    if (isTRUE(include_info)) "model_info",
    if (isTRUE(include_model)) "flexmix_model",
    if (include_info && adjustment_method != "unadjusted") "adjustment_method")
  missing_slots <- setdiff(required_slots, names(flexmix_list))
  if (length(missing_slots) > 0L){
    stop("`flexmix_list` is missing required slot information: ",paste(missing_slots, collapse = ", "))}
  
  # flexmix data
  flexmix_data <- flexmix_list[["flexmix_data"]]
  if (!inherits(flexmix_data, "data.frame")) {
    stop("`flexmix_data` must be a data.frame object.")}
  
  # check model name availability
  if (isTRUE(include_info) || isTRUE(include_model)){
    model_name <- as.character(model_name)
    if (length(model_name) != 1L ||
        is.na(model_name) ||
        !nzchar(model_name)) {
      stop("If model or adjustment information should be added to `input_obj`, 
           \n`model_name` must a non-empty character scalar.")}}

  # info slots
  if (include_info){ 
  if (adjustment_method != "unadjusted"){
    if (!"adjustment_method" %in% names(flexmix_list)){
      stop("`flexmix_list` is missing `adjustment_method` slot.")}
    adj_info <- flexmix_list[["adjustment_method"]]
    if (!is.list(adj_info)) {
      stop("`adjustment_method` must contain a list.")}
    if (!adjustment_method %in% names(adj_info)) {
      stop("No adjustment information available for `", adjustment_method, "`.")}
      adj_info <- adj_info[adjustment_method]}}
  
  # ---------------- select predictions ----------------
  
  df <- flexmix_data
  
  classification_cols <- c(unadjusted = "miQC.keep",
                           original_check = "miQC.keep_orig_adj",
                           geometry_check = "miQC.keep_geom_adj")
  
  # check availability & proper values for selected classification column
  selected_col <- unname(classification_cols[[adjustment_method]])
  if (!selected_col %in% colnames(df)) {
    stop("Selected classification column `", selected_col, "` is missing from `flexmix_data`.")}
  
  selected_values <- as.character(df[[selected_col]])
  unexpected <- setdiff(stats::na.omit(unique(selected_values)), c("keep", "discard"))
  if (length(unexpected) > 0L) {
    stop("Unexpected classification values: ", paste(unexpected, collapse = ", "))}
  
  # add final prediction
  df[["miQC.final_prediction"]] <- as.character(df[[selected_col]])
  
  geom_cols <- list(unadjusted = c("miQC.prob_invalid", "miQC.prob_valid", "miQC.keep_comp", "miQC.discard_comp"),
                    original_check = c("miQC.orig_adj", "miQC.orig_region", "miQC.below_boundary", "miQC.left_cutoff"),
                    geometry_check = c("miQC.geom_adj", "miQC.geom_region", "miQC.dist_up", "miQC.dist_low"))
                     
  transfer_cols <- switch(transfer_columns,
                          all = colnames(df),
                          method_specific = c("miQC.final_prediction", geom_cols[[adjustment_method]]),
                          final = "miQC.final_prediction")
  missing_cols <- setdiff(transfer_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop("Missing required prediction columns: ", paste(missing_cols, collapse = ", "))}
  
  # ---------------- identify input cells ----------------
  
  cells_input <- colnames(input_obj)
  cells_model <- rownames(df)

  if (anyNA(cells_input) || anyDuplicated(cells_input) || is.null(cells_input) ||
      anyNA(cells_model) || anyDuplicated(cells_model) || is.null(cells_model)) {
    stop("Unique, non-missing cell identifiers are required for input_obj and flexmix_list[[\"flexmix_data\"]].")}
  
  unexpected_cells <- setdiff(cells_model, cells_input)
  if (length(unexpected_cells) > 0L) {
    stop(length(unexpected_cells), " modelled cells are missing from `input_obj`.")}
  
  unexpected_cells <- setdiff(cells_input, cells_model)
  if (length(unexpected_cells) > 0L) {
    if (isTRUE(verbose)){
    message(length(unexpected_cells), " cells were omitted from fitting the model. 
            \nRelated values with be NA.")}}
  
  # ---------------- prepare transfer ----------------
  
  # add final prediction column
  transfer_df <- df[, transfer_cols, drop = FALSE]
  
  # Match to input_obj
  transfer_df <- transfer_df[match(cells_input, cells_model), ,drop = FALSE]
  rownames(transfer_df) <- cells_input
  
  # ---------------- check existing metadata ----------------
  
  existing_cols <- if (input_type == "Seurat") colnames(input_obj@meta.data) else colnames(SummarizedExperiment::colData(input_obj))
  
  collisions <- intersect(colnames(transfer_df), existing_cols)
  
  if (length(collisions) > 0L){
    if (!isTRUE(overwrite)){
    stop("Destination metadata columns already exist: ",
      paste(collisions, collapse = ", "),". 
      \nSet `overwrite = TRUE` to replace them.")}
    if (isTRUE(verbose)){
    message("Overwriting metadata columns: ", paste(collisions, collapse = ", "))}}
  
  if (isTRUE(include_info) || isTRUE(include_model)){
    existing_names <- if (input_type == "Seurat") names(input_obj@misc) else names(S4Vectors::metadata(input_obj))
    
  if (model_name %in% existing_names){
    if (!isTRUE(overwrite)){
      stop("Destination model slot for ",model_name," already exists. 
      \nSet `overwrite = TRUE` to replace it.")}
    if (isTRUE(verbose)){
    message("Overwriting model slot for ",model_name,".")}}}
  
  # ---------------- transfer metadata ----------------
  
  for (col_id in colnames(transfer_df)){
    if (input_type == "Seurat"){
      input_obj@meta.data[[col_id]] <- transfer_df[[col_id]]
    }else if (input_type == "SummarizedExperiment"){
      SummarizedExperiment::colData(input_obj)[[col_id]] <- transfer_df[[col_id]]}}
  
  if (isTRUE(verbose)){
    message("Transferred ", length(transfer_cols)," column(s) to metadata: ",paste(transfer_cols, collapse = ", "),". 
            \nFinal prediction `miQC.final_prediction` was added for `adjustment_method=", adjustment_method, "`.")}
  
  # ---------------- transfer model info ----------------

  info_to_add <- list()
  
  if (isTRUE(include_model)) {
    info_to_add[["flexmix_model"]] <- flexmix_list[["flexmix_model"]]}
  if (isTRUE(include_info)) {
    info_to_add[["model_info"]] <- flexmix_list[["model_info"]]
    info_to_add[["adjustment_method"]] <- if (adjustment_method != "unadjusted") adj_info else NULL
    info_to_add[["selected_classification"]] <- list(adjustment_method = adjustment_method,
                                                     source_column = selected_col,
                                                     destination_column = "miQC.final_prediction")}
  
  if (length(info_to_add) > 0L){
    if (input_type == "Seurat"){
      input_obj@misc[[model_name]] <- info_to_add
    }else if (input_type == "SummarizedExperiment"){
      S4Vectors::metadata(input_obj)[[model_name]] <- info_to_add}
    if (isTRUE(verbose)){
      type_info <- if (input_type == "Seurat") " to `misc` slot" else " to metadata"
      message("Added ", paste0(paste(names(info_to_add), collapse = ", "), type_info),".")}}

  return(input_obj)}

