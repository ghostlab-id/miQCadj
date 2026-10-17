#' Post-hoc Adjustments of cell classifications
#'
#' The function adds adjusted classifications and diagnostic measures to the
#' `data.frame` stored in `flexmix_data` and adds thresholds and the number of 
#' reassigned cells in `flexmix_list[["adjustment_method"]]`.
#'
#' @param flexmix_list The `list` object returned by `run_miQCadj()` with 
#'  required `flexmix_data` and `model_info` slots.
#' @param adjustment_method Character Scalar or Vector. Available options include 
#'  `original_check` performing the original miQC filtering adjustments with 
#'  `keep_all_below_boundary=TRUE` and `geometry_check` performing the modified 
#'  geometry check using the default parameters with `rescue_quantile = 0.5`,
#'  `discard_quantile = 0.95`, `converge_quantile = 0.05`,`safe_keep_dist = 1.25`. 
#'  For further adjustments the specific functions should be used. By default both
#'  checks will be performed.
#' @param min_index_type Character Scalar. The minimum of mitochondrial 
#'  percentage within discarded cells. Whether all cells should be checked for 
#'  this minimum `min_index_type = "all"` or whether to restrict the check to 
#'  discarded cells `min_index_type = "discard"`. Since the original wrapper used
#'  all cells, the default is \code{"all"}.
#' @param enforce_left_cutoff Logical Scalar. Whether cells with low library 
#'  complexity and cells with high mitochondrial percentage should be 
#'  discarded. Will only be used, if `adjustment_method = original_check`.
#'  Defaults to \code{TRUE}.
#' @param verbose Logical Scalar. Whether information regarding reassignments
#'  should be printed. Defaults to \code{TRUE}.
#' 
#' @return `flexmix_list` with updated slot information for:
#' \describe{
#'   \item{flexmix_data}{Updated data frame with adjustment results.}
#'   \item{adjustment_method}{Method-specific adjustment information.}
#' }
#' 
#' @seealso [miQCadj::original_adjustments()], [miQCadj::geometry_adjustments()]
#'
#' @examples
#' \dontrun{
#' # run both geometry checks
#' flexmix_list <- adjust_classification(flexmix_list)
#' }
#'
#' @export
adjust_classification <- function(flexmix_list, 
                                  adjustment_method = c("original_check", "geometry_check"),
                                  min_index_type = "all",
                                  enforce_left_cutoff = TRUE,
                                  verbose = TRUE){

  # check arguments
  logical_args <- list(enforce_left_cutoff = enforce_left_cutoff,
                       verbose = verbose)
  
  lapply(names(logical_args), function(x){
    check_x <- logical_args[[x]]
    if (!is.logical(check_x) ||
        length(check_x) != 1 ||
        is.na(check_x)) {
      stop("`",x,"` must be one non-missing logical value.")}})
  
  adjustment_method <- match.arg(adjustment_method, choices= c("original_check", "geometry_check"), several.ok = TRUE)
  
  if (!inherits(flexmix_list, "list")) {
    stop("`flexmix_list` must be a list object.")}
  
  required_slots <- c("flexmix_data", "model_info")
  missing_slots <- setdiff(required_slots, names(flexmix_list))
  if (length(missing_slots) > 0L){
    stop("`flexmix_list` is missing required slot information: ",paste(missing_slots, collapse = ", "))}
  
  # flexmix data
  flexmix_data <- flexmix_list[["flexmix_data"]]
  if (!inherits(flexmix_data, "data.frame")) {
    stop("`flexmix_data` must be a data.frame object.")}
  
  # model info
  model_info <- flexmix_list[["model_info"]]
  if (!"model_type" %in% names(model_info)){
    stop("Missing `model_type` information in `model_info`.")}
  model_type <- model_info[["model_type"]]
  
  if (model_type == "one_dimensional" &&
      "original_check" %in% adjustment_method &&
      isTRUE(enforce_left_cutoff) &&
      !"complexity_column" %in% colnames(flexmix_data)) {
    stop("`enforce_left_cutoff = TRUE` requires `complexity_column` for one-dimensional models.
         \nEither provide a complexity column for `run_miQCadj()` or set `enforce_left_cutoff=FALSE`.")}
  
  # ----------------  run geometry check --------------------#
    
  prev_cols <- list(original_check = c("miQC.keep_orig_adj", "miQC.orig_adj", "miQC.orig_region", "miQC.below_boundary","miQC.left_cutoff"),
                    geometry_check = c("miQC.keep_geom_adj", "miQC.geom_adj", "miQC.geom_region", "miQC.dist_up", "miQC.dist_low"))
  
  if (is.null(flexmix_list[["adjustment_method"]])){
    flexmix_list[["adjustment_method"]] <- list()}
  
  for (method in adjustment_method){
    drop_cols <- intersect(prev_cols[[method]], colnames(flexmix_data))
    if (length(drop_cols) > 0L) {
      flexmix_data[drop_cols] <- NULL
      if (isTRUE(verbose)){
      message("Recalculating ", method, " column and  ",method," information.")}}
  flexmix_list[["adjustment_method"]][[method]] <- NULL
      orig_cols <- colnames(flexmix_data)
    if (method == "original_check"){
      md_list <- original_adjustments(df=flexmix_data, keep_all_below_boundary = TRUE, min_index_type = min_index_type,
                                 enforce_left_cutoff = enforce_left_cutoff, verbose = verbose)}
    if (method == "geometry_check"){
      md_list <- geometry_adjustments(df=flexmix_data, model_type = model_type, verbose = verbose)}
    
    if (!inherits(md_list, "list") || !all(c("df", method) %in% names(md_list))){
      stop("Expected list object with `df` and `",method,"` slots returned by `",gsub("check", "", method),"_adjustments()` function.")}
    flexmix_list[["adjustment_method"]][[method]] <- md_list[[method]]
    df <- md_list[["df"]]
    new_cols <- setdiff(colnames(df), orig_cols)
    # Verify expected columns
    missing_cols <- setdiff(prev_cols[[method]], new_cols)
    if (length(missing_cols) > 0L) {
      stop("Missing ",method," adjustment columns: ", paste(missing_cols, collapse = ", "))}
    # verify expected cell ids
    if (anyDuplicated(rownames(df)) || !identical(rownames(df), rownames(flexmix_data))){
      stop("Adjustment results must contain exactly the same cell IDs as `flexmix_data`.")}
    flexmix_data <- cbind(flexmix_data, df[ , prev_cols[[method]], drop=FALSE])}
  
  shared_vector <- rep(NA_character_, nrow(flexmix_data))
  if (all(c("miQC.keep_orig_adj", "miQC.keep_geom_adj") %in% colnames(flexmix_data))){
    orig <- as.character(flexmix_data[["miQC.keep_orig_adj"]])
    geom <- as.character(flexmix_data[["miQC.keep_geom_adj"]])
    non_na <- !is.na(orig) & !is.na(geom)
    
    shared_vector[non_na & orig == "keep" & geom == "keep"] <- "both_keep"
    shared_vector[non_na & orig == "discard" & geom == "discard"] <- "both_discard"
    shared_vector[non_na & orig == "keep" & geom == "discard"] <- "keep_orig"
    shared_vector[non_na & orig == "discard" & geom == "keep"] <- "keep_geom"}
  
  flexmix_data[["miQC.keep_adj_compared"]] <- shared_vector

  flexmix_list[["flexmix_data"]] <- flexmix_data
  return(flexmix_list)}

