#' Visualize fitted miQC components and posterior classifications
#'
#' Creates diagnostic scatter plots of mitochondrial percentage versus the
#' complexity metric used in \code{run_miQCadj()}, including keep/discard 
#' classifications based on a chosen posterior cutoff, and optionally, fitted 
#' mixture component aggregation lines.
#'
#' @param flexmix_list The `list` object returned by `run_miQCadj()` with 
#'  required `flexmix_data` and `model_info` slots.
#' @param sample_name Character scalar, optional. Sample name to be used for 
#'  the plots. A prefix \code{paste0(sample_name, ": ")} will either be used for 
#'  each plot, if `combine_plots=FALSE`, or for the global plot only, if 
#'  `combine_plots=TRUE`. Defaults to \code{NULL}.
#' @param use_median_columns Logical Scalar. Whether component classification
#'  based on aggregated low-quality mitochondrial should be used. Only possible,
#'  if median-related component selection columns \code{"_median"} are available 
#'  in the `flexmix_data` data frame. Defaults to \code{FALSE}.
#' @param x_values Character Scalar. Available options include `model_column` 
#' (default) or `complexity_column`. Will return identical plots for linear, 
#'  spline, polynomial, and polynomial_raw models. For one-dimensional models 
#'  `model_column` will return plots with decreasing mitochondrial percentage 
#'  using ranked cell indices in the x-axis. If `complexity_column` had been 
#'  provided to `run_miQCadj()` values from `complexity_column` can be used as 
#'  `x_values`.
#' @param posterior_cutoff Numeric Scalar between 0 and 1. Cutoff used for 
#'  keep/discard classification. Defaults to \code{0.75}.
#' @param point_size Non-negative Numeric Scalar. Point size used for scatter
#'  plots. Defaults to \code{1}.
#' @param point_alpha Non-negative Numeric Scalar. Transparency used for scatter 
#' plots. Defaults to \code{0.8}.
#' @param color_lines Named character vector, optional. Specifying the colors 
#'  for `keep` and `discard` classifications, will plot lines for 
#'  keep/discard classification curves based on aggregated medians. Defaults to 
#'  \code{c(discard = "black", keep = "blue")}.
#' @param color_points Named character vector for keep/discard point 
#' color values with names \code{c("discard", "keep")}.
#' Defaults to \code{c(discard = "dimgrey", keep = "cornflowerblue")}.
#' @param combine_plots Logical Scalar. Whether the list of plots should be 
#'  combined into a 2-column-grid using `patchwork`. Defaults to \code{TRUE}.
#'
#' @return If `combine_plots=FALSE` named list of `ggplot` objects:
#' \itemize{
#'   \item{classification}{Plot with miQC classification based on posterior 
#'   probability thresholds.}.
#'   \item{probability}{Plot with miQC probability for each cell being related 
#'   to low-quality classification.}
#'   }
#' If `combine_plots=TRUE` (default) the items will be combined into a 1x2-grid 
#' using patchwork with `sample_name` being used for the global title, only.
#'
#' @examples
#' \dontrun{
#' fm <- run_miQCadj(input_obj = obj)
#' 
#' plots <- check_miQC_components(
#'   flexmix_list = fm,
#'   sample_name = "Sample")
#'
#' # show combined plot grid
#' print(plots)
#' }
#'
#' @export
check_miQC_components <- function(flexmix_list, sample_name = NULL, 
                                  use_median_columns = FALSE, x_values = "model_column",
                                  posterior_cutoff = 0.75, point_size = 1, point_alpha = 0.8,
                                  color_lines = c(discard = "black", keep = "blue"),
                                  color_points = c(discard = "dimgrey", keep = "cornflowerblue"),
                                  combine_plots=TRUE){
  # checks
  x_values <- match.arg(x_values, choices= c("model_column", "complexity_column"))
  
  if (!inherits(flexmix_list, "list")) {
    stop("`flexmix_list` must be a list object.")}
  
  required_slots <- c("flexmix_data", "model_info")
  missing_slots <- setdiff(required_slots, names(flexmix_list))
    if (length(missing_slots) > 0L){
      stop("`flexmix_list` is missing required slot information: ",paste(missing_slots, collapse = ", "))}
  
  # check logical arguments
  logical_args <- list(use_median_columns = use_median_columns,
                       combine_plots = combine_plots)
  
  lapply(names(logical_args), function(x){
    check_x <- logical_args[[x]]
    if (!is.logical(check_x) ||
        length(check_x) != 1 ||
        is.na(check_x)) {
      stop("`",x,"` must be one non-missing logical value.")}})
  
  # flexmix data
  flexmix_data <- flexmix_list[["flexmix_data"]]
  if (!inherits(flexmix_data, "data.frame")) {
    stop("`flexmix_data` must be a data.frame object.")}
  
  # model info
  model_info <- flexmix_list[["model_info"]]
  required_info <- c("model_type", "column_names", "components_median")
  missing_info <- setdiff(required_info, names(model_info))
  if (length(missing_info) > 0L){
    stop("Missing required slot information in `model_info`: ", paste(missing_info, collapse = ", "))}
  
  # model_type
  model_type <- model_info[["model_type"]]
  
  # column_names
  missing_info <- setdiff("complexity_column", names(model_info[["column_names"]]))
  if (length(missing_info) > 0L){
    stop("Missing required slot information in `model_info$column_names`: ", paste(missing_info, collapse = ", "))}
  
  # median selection
  missing_info <- setdiff("median_selection_added", names(model_info[["components_median"]]))
  if (length(missing_info) > 0L){
    stop("Missing required slot information in `model_info$components_median`: ", paste(missing_info, collapse = ", "))}
  
  if (isTRUE(use_median_columns) && !isTRUE(model_info[["components_median"]][["median_selection_added"]])){
    stop("Median-based component selection wasn't added by `run_miQCadj()`.")}

  if (!is.numeric(posterior_cutoff) ||
      length(posterior_cutoff) != 1L ||
      is.na(posterior_cutoff) ||
      posterior_cutoff < 0 ||
      posterior_cutoff > 1) {
    stop("`posterior_cutoff` must be a numeric value between 0 and 1.")}
  
  if (!is.null(sample_name)){
    sample_name <- as.character(sample_name)
      if (length(sample_name) != 1L ||
      is.na(sample_name) ||
      !nzchar(sample_name)) {
    stop("`sample_name` must be `NULL` or a non-empty character scalar.")}}
  
  if (isTRUE(combine_plots)){
    prefix <- ""
    global_prefix <- if (!is.null(sample_name)) paste0(sample_name, ": ") else ""  
  }else{
    prefix <- if (!is.null(sample_name)) paste0(sample_name, ": ") else "" }
  
  # columns added by `run_miQCadj()`
  p_invalid <- "miQC.prob_invalid"
  p_invalid <- if (isTRUE(use_median_columns)) paste0(p_invalid, "_median") else p_invalid
  
  keep_comp <- "miQC.keep_comp"
  keep_comp  <- if (isTRUE(use_median_columns)) paste0(keep_comp, "_median") else keep_comp
  discard_comp <- "miQC.discard_comp"
  discard_comp  <- if (isTRUE(use_median_columns)) paste0(discard_comp, "_median") else discard_comp
  
  if (!is.null(color_lines)){
    model_required <- c(p_invalid, keep_comp, discard_comp)
  }else{
    model_required <- p_invalid}
  
  # required columns  
  if (x_values == "model_column") {
    complexity_column <- if (model_type == "one_dimensional") NULL else "complexity_column"
  } else {
    if (is.null(model_info[["column_names"]][["complexity_column"]])) {
      stop("No complexity column was provided to `run_miQCadj()`.")}
    complexity_column <- "complexity_column"}
  required_cols <- c("percent_mt", complexity_column, model_required)
  missing_cols <- setdiff(required_cols, colnames(flexmix_data))
  if (length(missing_cols) > 0L) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))}
  
  #-------------------------------------------------------#
  
  # get metdata
  md <- flexmix_data[ , required_cols, drop = FALSE]
  rownames(md) <- rownames(flexmix_data)

  # ensure numeric values
  for (col_id in required_cols){
      md[[col_id]] <- suppressWarnings(as.numeric(as.character(md[[col_id]])))}
  
  # check for invalid numeric values
  check_numeric <- md[, required_cols, drop = FALSE]
  if (any(!is.finite(as.matrix(check_numeric)))) {
    stop( "Unexpected missing or non-finite numeric values.")}
  
  # invalid probability
  md[["P(LowQuality)"]] <- md[[p_invalid]]
    
  # classification
  md$class <- ifelse(md[["P(LowQuality)"]] <= posterior_cutoff, "keep", "discard")
  md$class <- factor(md$class, levels = c("discard", "keep"))
  
  # sort indices
  if (is.null(complexity_column)) {
    md <- md[order(md[["percent_mt"]], decreasing = TRUE), , drop = FALSE]
    complexity_column <- "ranked_cell_index"
    md[[complexity_column]] <- seq_len(nrow(md))}
  
  # curve df
  curve_df <- NULL
  if (!is.null(color_lines)){
    curve_df <- data.frame(feature = md[[complexity_column]],
                           keep_y = md[[keep_comp]],
                           discard_y = md[[discard_comp]])
    curve_df <- curve_df[order(curve_df$feature), , drop = FALSE]
    curve_df <- stats::aggregate(cbind(keep_y, discard_y) ~ feature, data = curve_df, FUN = stats::median)}
  
  # create plots  
  plots <- list()
  plots[["classification"]] <- miQC_build_plot(df = md, col_var = "class",
        complexity_column = complexity_column, percent_mt = "percent_mt",
        title_txt = paste0(prefix, "miQC Classification"),
        color_title= "Classification",
        point_alpha = point_alpha, point_size = point_size,
        discrete = TRUE, curve_df = curve_df, color_lines = color_lines,
        color_points = color_points)
  
  plots[["probability"]] <- miQC_build_plot(df = md, col_var = "P(LowQuality)",
        complexity_column = complexity_column, percent_mt = "percent_mt",
        title_txt = paste0(prefix, "miQC Probability"),
        color_title= "P(LowQuality)",
        point_alpha = point_alpha, point_size = point_size,
        discrete = FALSE, curve_df = curve_df, color_lines = color_lines,
        color_points = color_points)
  
  if (isTRUE(combine_plots)){
    plots <- patchwork::wrap_plots(plots, ncol=2, guides = "collect") +
      patchwork::plot_annotation(title= paste0(global_prefix, "Cell Classification and posterior Probabilities"))}
  
  return(plots)}


