#' Define components and posterior cutoff for cell classification
#'
#' @param flexmix_list The `list` object returned by `run_miQCadj()` with 
#'  required `flexmix_data` and `model_info` slots.
#' @param component_selection Character scalar. Available options are `median` 
#'  or `intercept` (default) resulting in the corresponding component selection. 
#' @param posterior_cutoff Numeric between 0 and 1. Defaults to \code{0.75}.
#'
#' @return `flexmix_list` with updated slot information for:
#' \itemize{
#' \item{`flexmix_data`}{Updated `data.frame`: If `comp_selection == intercept`
#'  median-related columns are discarded. If `comp_selection == median` the 
#'  median-related columns replace the intercept-derived `miQC.prob_invalid`, 
#'  `miQC.prob_valid`, `miQC.keep_comp`, and `miQC.discard_comp` columns.
#'  A new column `miQC.keep` is added, comprising cell classifications based on #
#'  posterior thresholds and the updated `miQC.prob_invalid` column.}
#'  \item{`model_info`}{Updated list with `component_selection` and 
#'  `posterior_cutoff` as specified in the function arguments.} 
#' }
#'
#' @examples
#' \dontrun{
#' flexmix_list <- define_classification(
#'   flexmix_list,
#'   component_selection = "median",
#'   posterior_cutoff = 0.75
#' )
#' }
#'
#' @export
define_classification <- function(flexmix_list, 
                                  component_selection = "intercept", 
                                  posterior_cutoff = 0.75){
  
  # checks
  if (!inherits(flexmix_list, "list")) {
    stop("`flexmix_list` must be a list object.")}

  required_slots <- c("flexmix_data", "model_info")
  missing_slots <- setdiff(required_slots, names(flexmix_list))
  if (length(missing_slots) > 0L){
    stop("`flexmix_list` is missing required slot information: ",paste(missing_slots, collapse = ", "))}
  
  # model info
  model_info <- flexmix_list[["model_info"]]
  if (!inherits(model_info, "list")) {
    stop("`model_info` must be a list object.")}
  
  # flexmix data
  flexmix_data <- flexmix_list[["flexmix_data"]]
  if (!inherits(flexmix_data, "data.frame")) {
    stop("`flexmix_data` must be a data.frame object.")}

  component_selection <- match.arg(component_selection, choices = c("median", "intercept"))
  
  if (!is.numeric(posterior_cutoff) ||
      length(posterior_cutoff) != 1 ||
      is.na(posterior_cutoff) ||
      posterior_cutoff < 0 ||
      posterior_cutoff > 1) {
    stop("posterior_cutoff must be between 0 and 1.")}

  select_cols <- c("miQC.prob_invalid",
                "miQC.prob_valid",
                "miQC.keep_comp",
                "miQC.discard_comp")

  median_cols <- paste0(select_cols, "_median")
  
  md <- flexmix_data  
  missing_cols <- setdiff(select_cols, colnames(md))
 if (length(missing_cols)> 0L){
   stop("`flexmix_data` is missing required columns: ", paste(missing_cols, collapse = ", "))}
  
  if (component_selection == "median"){
    missing_cols <- setdiff(median_cols, colnames(md))
    if (length(missing_cols)> 0L){
      stop("`flexmix_data` is missing required columns: ", paste(missing_cols, collapse = ", "))}
    md[, select_cols] <- md[, median_cols]}
  
   if(any(median_cols %in% colnames(md))){
      md[, median_cols] <- NULL}

  # select component
  md[["miQC.keep"]] <- ifelse(md[["miQC.prob_invalid"]] <= posterior_cutoff, "keep", "discard")
  
  if (all(md[["miQC.keep"]] == "keep")){
    if ("seed_used" %in% names(model_info)){
      current_seed <- paste0(" (current seed ", model_info[["seed_used"]],")")
    }else{
      warning("No seed information available in `model_info` slot.")}
    warning("No cells were filtered using `posterior_cutoff=",posterior_cutoff,"`.
         \nCheck plots returned by `check_miQC_components()`. Adjust the `posterior_cutoff` or 
         \nrerun `run_miQCadj()` using a distinct `model_type`, `complexity_column` and/or seed",current_seed,".")}
  
  # update data
  flexmix_list[["flexmix_data"]] <- md
  
  model_info[["component_selection"]] <- component_selection
  model_info[["posterior_cutoff"]] <- posterior_cutoff
  flexmix_list[["model_info"]] <- model_info
  
  return(flexmix_list)}


