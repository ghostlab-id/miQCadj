#' Original geometric validation of selected miQC classifications
#' 
#' Performs a post hoc geometric plausibility check of an existing miQC 
#' classification. Cells may be reassigned based on the fitted keep-component 
#' prediction and, optionally, a complexity-dependent left cutoff.
#'
#' The function adds adjusted classifications and diagnostic measures to the
#' provided `data.frame` and records thresholds and the number of reassigned 
#' cells in `original_check`.
#'
#' @param df The `data.frame` object for `flexmix_data`.
#' @param keep_all_below_boundary Logical Scalar. Whether cells with 
#'  mitochondrial percentage below the fitted keep component prediction should 
#'  be retained regardless of their posterior probabilities. Defaults to 
#'  \code{TRUE}.
#' @param min_index_type Character Scalar. Whether all cells should be checked
#'  for the minimum of mitochondrial percentage within discarded cells 
#'  `min_index_type = "all"` or whether to restrict the check to discarded cells 
#'  `min_index_type = "discard"`. Since the original wrapper used all cells, 
#'  the default is \code{"all"}.
#' @param enforce_left_cutoff Logical Scalar. Whether cells with low library 
#'  complexity and cells with high mitochondrial percentage should be 
#'  discarded. Requires `complexity_column` within `flexmix_data` and at least 
#'  one cell being discarded after applying `keep_all_below_boundary`. 
#'  Defaults to \code{TRUE}.
#' @param verbose Logical Scalar. Whether information regarding reassignments
#'  should be printed. Defaults to \code{TRUE}.
#' 
#' @return Named list with slot information for:
#' \describe{
#' \item{df}{Updated `data.frame` with geometry-adjusted `miQC.keep_orig_adj`
#'  classification column and geometry-related columns `miQC.below_boundary`, 
#'  `miQC.left_cutoff`, `miQC.orig_adj`, `miQC.orig_region`.}
#'  \item{original_check}{Adjustment thresholds used for geometry checks,
#'  including the specified function arguments.} 
#' }
#'
#' @examples
#' \dontrun{
#' df_list <- original_adjustments(df=flexmix_data)
#' }
#'
#' @export
original_adjustments <- function (df, keep_all_below_boundary = TRUE, min_index_type="all",
                                  enforce_left_cutoff = TRUE, verbose = TRUE){
  
  #----------- check input -------------#

  # check logical arguments
  logical_args <- list(keep_all_below_boundary = keep_all_below_boundary,
                       enforce_left_cutoff = enforce_left_cutoff,
                       verbose = verbose)
  
  lapply(names(logical_args), function(x){
    check_x <- logical_args[[x]]
    if (!is.logical(check_x) ||
        length(check_x) != 1 ||
        is.na(check_x)) {
      stop("`",x,"` must be one non-missing logical value.")}})
  
  min_index_type <- match.arg(min_index_type, choices= c("all", "discard"))
  
  # df
  if (!inherits(df, "data.frame")) {
    stop("`df` must be a data.frame object.")}
  
  orig_df <- df
  if (anyDuplicated(rownames(orig_df)) || anyNA(rownames(orig_df))){
    stop("Unique, non-NA row identifiers required.")}
  
  # required columns
  # Complexity is optional for diagnostics,
  # but required when enforcing the left cutoff
  complexity_column <- if ("complexity_column" %in% colnames(df)) "complexity_column" else  NULL
  if (isTRUE(enforce_left_cutoff) && is.null(complexity_column)) {
    stop("`enforce_left_cutoff=TRUE` requires `complexity_column` in `df`.")}
  
  keep_col <- "miQC.keep"
  keep_comp <- "miQC.keep_comp"
  
  required_cols <- c("percent_mt", complexity_column, keep_col, keep_comp)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))}
  
  # ensure numeric values
  numeric_cols <- c("percent_mt", complexity_column, keep_comp)
  for (col_id in numeric_cols){
    df[[col_id]] <- suppressWarnings(as.numeric(as.character(df[[col_id]])))}
  
  # check for invalid numeric values
  check_numeric <- df[, numeric_cols, drop = FALSE]
  if (any(!is.finite(as.matrix(check_numeric)))) {
    stop( "Unexpected missing or non-finite numeric values.")}
  
  # check keep_col labels
  keep_vals <- as.character(df[[keep_col]])
  unexpected <- setdiff(stats::na.omit(unique(keep_vals)), c("keep", "discard"))
  if (length(unexpected)> 0L) {
    stop("Got unexpected value for `miQC.keep`: ",unexpected, ". 
         \nExpected values are keep/discard.")}
  
  # subset to classified cells
  df <- df[keep_vals %in% c("keep", "discard"), , drop=FALSE]
  
  if (nrow(df) == 0L) {
    stop("No classified cells available for geometry adjustment.")}
 
  #----------- adjust classifications -------------#
  
  # keep below boundary
  df[["miQC.keep_orig_adj"]] <- as.character(df[["miQC.keep"]] )
  df[["miQC.below_boundary"]] <- df[["percent_mt"]] < df[["miQC.keep_comp"]]

  if (isTRUE(keep_all_below_boundary)) {
    df[["miQC.keep_orig_adj"]][df[["miQC.below_boundary"]]] <-  "keep"}

  n_keep_below_boundary <- sum(df[["miQC.keep_orig_adj"]]=="keep" & df[["miQC.below_boundary"]], na.rm = TRUE) - 
    sum(df[["miQC.keep"]]=="keep" & df[["miQC.below_boundary"]], na.rm = TRUE)
  
  # Initialize left-cutoff diagnostics
  df[["miQC.left_cutoff"]] <- rep(NA, nrow(df))
  min_discard <- NULL
  min_index_all <- NULL
  min_index_discard <- NULL
  min_index <- NULL
  lib_complexity <- NULL
    
  # discard by left cutoff    
  discard_idx <- df[["miQC.keep_orig_adj"]] =="discard"
  
  adj_info <- if (isTRUE(enforce_left_cutoff)) " Skipping left cutoff adjustments." else NULL
  if (!any(discard_idx)){
    if (isTRUE(verbose)){
      message("No discarded cells available for determining  the left-cutoff threshold.")}

  }else{
    min_discard <- min(df[["percent_mt"]][discard_idx])
    min_index_all <- which(df[["percent_mt"]] == min_discard)[1]
    min_index_discard <- which(discard_idx & df[["percent_mt"]] == min_discard)[1L]
    min_index <- if (min_index_type == "all") min_index_all else min_index_discard
    
    min_index_cellID <- if (!is.null(min_index)) rownames(df)[min_index] else NULL
    
    if (!is.null(complexity_column)){
      lib_complexity <- df[["complexity_column"]][min_index]
      
      low_complexity <- df[["complexity_column"]] <= lib_complexity
      high_mt <- df[["percent_mt"]] >= min_discard
      df[["miQC.left_cutoff"]] <- low_complexity & high_mt
      
      if (isTRUE(enforce_left_cutoff)) {  
        n_discard_left_cutoff <- sum(df[["miQC.left_cutoff"]] & df[["miQC.keep_orig_adj"]] == "keep", na.rm = TRUE)
        df[["miQC.keep_orig_adj"]][df[["miQC.left_cutoff"]] %in% TRUE] <-"discard"}}}
  
  # Identify reassignments
  override_discard <- (df[["miQC.keep"]] == "discard" & df[["miQC.keep_orig_adj"]] == "keep")
  override_keep <- (df[["miQC.keep"]] == "keep" & df[["miQC.keep_orig_adj"]] == "discard")
  
  df[["miQC.orig_adj"]] <- rep("unchanged", nrow(df))
  df[["miQC.orig_adj"]][override_keep] <- "geom_fail"
  df[["miQC.orig_adj"]][override_discard] <- "geom_pass"
  
  # Geometry regions
  df[["miQC.orig_region"]] <- rep("unchanged", nrow(df))
  df[["miQC.orig_region"]][df[["miQC.below_boundary"]] %in% TRUE] <- "below_boundary"
  df[["miQC.orig_region"]][df[["miQC.left_cutoff"]] %in% TRUE] <- "left_cutoff"
  
  if (isTRUE(verbose)) {
    message(paste("Original geometry-related Reassignments: ",
                  paste(capture.output(
                    table(df[,c("miQC.below_boundary","miQC.left_cutoff")], useNA = "ifany")), 
                    collapse = "\n"), sep = "\n"))}
  
  original_keep <- as.character(df[["miQC.keep"]])
  adjusted_keep <- as.character(df[["miQC.keep_orig_adj"]])
  
  n_discard_left_cutoff <- sum(df[["miQC.keep_orig_adj"]] == "discard" & df[["miQC.left_cutoff"]], na.rm = TRUE) - 
    sum(df[["miQC.keep"]] == "discard" & df[["miQC.left_cutoff"]], na.rm = TRUE)
  
  new_cols <- setdiff(colnames(df), colnames(orig_df))
  df <- df[rownames(orig_df), new_cols, drop = FALSE]
  orig_df <- cbind(orig_df, df)
  
  return(list(df= orig_df,
                original_check = list(adjustment_method = "original",
                                      keep_all_below_boundary = keep_all_below_boundary,
                                      enforce_left_cutoff = enforce_left_cutoff,
                                      min_discard = min_discard,
                                      min_index = list(min_index_type = min_index_type,
                                                       min_index_used = min_index,
                                                       min_index_cellID =  min_index_cellID,
                                                       min_index_all = min_index_all,
                                                       min_index_discard = min_index_discard),
                                      lib_complexity = lib_complexity,
                                      n_reassigned_total = sum(override_discard, na.rm = TRUE) + sum(override_keep, na.rm=TRUE),
                                      n_override_discard = sum(override_discard, na.rm = TRUE),
                                      n_below_boundary = sum(df[["miQC.below_boundary"]], na.rm = TRUE),
                                      n_keep_below_boundary = n_keep_below_boundary,
                                      n_override_keep = sum(override_keep, na.rm = TRUE),
                                      n_left_cutoff = sum(df[["miQC.left_cutoff"]], na.rm = TRUE),
                                      n_discard_left_cutoff = n_discard_left_cutoff)))}

