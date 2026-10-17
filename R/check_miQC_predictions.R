#' Visualize miQC classifications after geometry validation
#'
#' Creates diagnostic scatter plots for miQC predictions after 
#' `adjust_classification()`. 
#'
#' @param flexmix_list The `list` object returned by `run_miQCadj()` with 
#'  required `flexmix_data` and `model_info` slots.
#' @param sample_name Character scalar, optional. Sample name to be used for 
#'  the plots. A prefix \code{paste0(sample_name, ": ")} will either be used for 
#'  each plot, if `combine_plots=FALSE`, or for the global plot only, if 
#'  `combine_plots=TRUE`. Defaults to \code{NULL}.
#' @param point_size Positive Numeric. Specification of the point size. 
#'  Defaults to \code{1}.
#' @param add_unadjusted_prediction Logical Scalar. Whether the unadjusted 
#'  prediction plot for `"miQC.keep"`, should be re-created. Defaults to 
#'  \code{FALSE}.
#'  @param add_compared_prediction Logical Scalar. Whether the plot comparing 
#'  the post-hoc prediction adjustments `"miQC.keep_adj_compared"` created by 
#'  `adjust_classification()`, if both adjustment methods were applied, should 
#'  should be created. Defaults to \code{FALSE}.
#'  @param add_probability Logical Scalar. Whether the plot with the 
#'  probability for each cell belonging to the compromised compartment 
#'  `"miQC.prob_invalid"`, should be re-created. Defaults to \code{FALSE}.
#' @param adjustment_method Character Scalar or Vector, optional. Name of 
#' method(s) used for geometric prediction adjustments or \code{NULL}. Available 
#'  options include `original_check` with required columns `miQC.keep_orig_adj`, 
#'  `miQC.orig_adj`, and `miQC.orig_region` and `geometry_check` with 
#'  `miQC.keep_geom_adj`, `miQC.geom_adj`, and `miQC.geom_region` in 
#'  `flexmix_data`. By default plots for both methods will be created.
#' @param x_values Character Scalar. Available options include `model_column` 
#' (default) or `complexity_column`. Will return identical plots for linear, 
#'  spline, polynomial, and polynomial_raw models. For one-dimensional models 
#'  `model_column` will return plots with decreasing mitochondrial percentage 
#'  using ranked cell indices in the x-axis. If `complexity_column` had been 
#'  provided to `run_miQCadj()` values from `complexity_column` can be used as 
#'  `x_values`.
#' @param point_alpha Positive Numeric. Specification of the transparency of 
#'  the points. Defaults to \code{1}.
#' @param color_lines Named character vector, optional. Specifying the colors 
#' for `keep` and `discard` classifications, will plot lines for 
#' keep/discard classification curves based on aggregated medians. Defaults to 
#' \code{c(discard = "black", keep = "blue")}.
#' @param color_points Named character vector for keep/discard point 
#' color values with names \code{c("discard", "keep")}. Defaults to 
#' \code{c(discard = "dimgrey", keep = "cornflowerblue")}
#' @param shape_geom Named character vector specifying the shapes for 
#'  geometry-based adjustments. Defaults to \code{c(unchanged = 1, 
#'  geom_pass = 15, geom_fail = 4)}.
#' @param color_regions Named character vector specifying the colors for 
#'  geometry-based regions. Defaults to \code{c(valid="#81F168", 
#'  invalid="#E28861", converge="gold")}.
#' @param color_compare Named character vector specifying the colors for 
#'  comparing predictions by `original_check` and `geometry_check`.  
#'  Defaults to \code{c(both_keep = "#6495ED", both_discard = "#A9A9A9", 
#'  keep_orig= "#6932C2", keep_geom="#123EA3")}
#' @param combine_plots Logical Scalar. Whether the list of plots should be 
#'  combined into a 2-column-grid using `patchwork`. Defaults to \code{TRUE}.
#'
#' @return  If `combine_plots=FALSE` named list of `ggplot` objects. Depending 
#' on the requested adjustment_method, and the additional plots, it may include:
#' \itemize{
#'   \item{`unadjusted_classification`}{Original keep/discard classifications.}
#'   \item{`prediction_comparison`}{Comparing keep/discard classifications from
#'   `original_check` and `geometry_check`. Can only be created, if both 
#'   adjustment methods were applied, though the plot can be created
#'   independently.}
#'   \item{`probability`}{Posterior invalid probabilities.}
#'   \item{`original_adjusted_classification`}{Adjusted classifications using 
#'   original prediction adjustments.}
#'   \item{`original_region_adjustments`}{Region-related adjustments performed
#'   by original prediction adjustments.}
#'   \item{`geometry_adjusted_classification`}{Adjusted classifications using 
#'   geometruc prediction adjustments.}
#'   \item{`geometry_region_adjustments`}{Region-related adjustments performed
#'   by geometric prediction adjustments.}
#' }
#' If `combine_plots=TRUE` (default) the items will be combined into a 2-column
#' grid using patchwork with `sample_name` being used for the global title, only.
#' 
#' 
#' @examples
#' \dontrun{
#' plots <- check_miQC_predictions(flexmix_list, sample_name = "Sample")
#' # show combined plot grid
#' print(plots)
#' }
#'
#' @export
check_miQC_predictions <- function(flexmix_list, sample_name = NULL, 
                                   add_unadjusted_prediction=FALSE,
                                   add_compared_prediction = FALSE,
                                   add_probability=FALSE, 
                                   adjustment_method= c("original_check", "geometry_check"),
                                   x_values = "model_column",
                                   point_size = 1, point_alpha = 0.8,
                                   color_lines = c(discard = "black", keep = "blue"),
                                   color_points = c(discard = "dimgrey", keep = "cornflowerblue"),
                                   shape_geom = c(unchanged = 1, geom_pass = 15, geom_fail = 4),
                                   color_regions= list(geometry_check = c(valid="#99D696", converge="#F2DB68",invalid="#D6A696"),
                                                       original_check = c(unchanged="#99D696", below_boundary="#F2DB68",left_cutoff="#D6A696")),
                                   color_compare = c(both_keep = "#6495ED", both_discard = "#A9A9A9", keep_orig= "#6932C2", keep_geom="#123EA3"),
                                   combine_plots =TRUE){
  
  # ------------ sanity check ----------------
  # checks
  x_values <- match.arg(x_values, choices= c("model_column", "complexity_column"))
  if (!is.null(adjustment_method)){
  adjustment_method <- match.arg(adjustment_method, choices= c("original_check", "geometry_check"), several.ok = TRUE)}
  
  # check logical arguments
  logical_args <- list(add_unadjusted_prediction = add_unadjusted_prediction,
                       add_compared_prediction = add_compared_prediction,
                       add_probability = add_probability,
                       combine_plots = combine_plots)
  
  lapply(names(logical_args), function(x){
    check_x <- logical_args[[x]]
    if (!is.logical(check_x) ||
        length(check_x) != 1 ||
        is.na(check_x)) {
      stop("`",x,"` must be one non-missing logical value.")}})
  
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
  
  # required flexmix_data columns
  p_invalid <- if (isTRUE(add_probability)) "miQC.prob_invalid" else NULL
  keep_compare <- if (isTRUE(add_compared_prediction)) "miQC.keep_adj_compared" else NULL
  keep_col <- "miQC.keep"
  keep_comp <- "miQC.keep_comp"
  discard_comp <- "miQC.discard_comp"
  
  keep_adj <- list(geometry_check = "miQC.keep_geom_adj", 
                   original_check = "miQC.keep_orig_adj")
  
  geom_adj <- list(geometry_check = "miQC.geom_adj", 
                   original_check = "miQC.orig_adj")
  
  geom_region <- list(geometry_check = "miQC.geom_region", 
                      original_check = "miQC.orig_region")

  # required columns  
  comp_cols <- if (!is.null(color_lines)) c(keep_comp, discard_comp) else NULL
  general_cols <- c(p_invalid, keep_col, comp_cols, keep_compare)
  specific_cols <- NULL
  if (!is.null(adjustment_method)){
  specific_cols <- unlist(c(keep_adj[adjustment_method], 
                            geom_adj[adjustment_method], 
                            geom_region[adjustment_method]), use.names = FALSE)}
  model_required <- c(general_cols, specific_cols)
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
  
  # get metadata
  md <- flexmix_data[ , required_cols, drop = FALSE]
  rownames(md) <- rownames(flexmix_data)
  
  # check keep_col labels
  keep_vals <- as.character(md[[keep_col]])
  unexpected <- setdiff(stats::na.omit(unique(keep_vals)), c("keep", "discard"))
  if (length(unexpected)> 0L) {
    stop("Got an unexpected value for cell prediction: ",paste0(unexpected, collapse = ", "), ". 
         \nExpected values are keep/discard.")}
  
  md <- md[keep_vals %in% c("keep", "discard"), , drop=FALSE]
  
  # ensure numeric values for classified cells
  numeric_cols <- c("percent_mt", complexity_column, comp_cols)
  for (col_id in numeric_cols){
    md[[col_id]] <- suppressWarnings(as.numeric(as.character(md[[col_id]])))}
  
  # check for invalid numeric values
  check_numeric <- md[, numeric_cols, drop = FALSE]
  if (any(!is.finite(as.matrix(check_numeric)))) {
    stop( "Unexpected missing or non-finite numeric values.")}
  
  # check availability of miQC.keep_adj_compared values, when requested
  if (isTRUE(add_compared_prediction)) {
  if (all(is.na(as.character(md[["miQC.keep_adj_compared"]])))) {
    stop("No comparison results are available. 
           \nApply both `original_check` and `geometry_check` before requesting an adjustment comparison.")}}
  
  # ---------------- curve df ----------------
  
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
  
  # ---------------- create plots ----------------
  
  # general required column: percent_mt
  plots <- list()
  
  # ---------------- unadjusted prediction ----------------
  
  if(isTRUE(add_unadjusted_prediction)){
  # unadjusted prediction
  md[["miQC.keep"]] <- factor(md[["miQC.keep"]], levels = c("discard", "keep"))
  plots[["unadjusted_classification"]] <- miQC_build_plot(df = md, col_var = "miQC.keep",
                                               complexity_column = complexity_column, percent_mt = "percent_mt",
                                               title_txt = paste0(prefix, "Unadjusted Classification"),
                                               color_title= "Classification",
                                               point_alpha = point_alpha, point_size = point_size,
                                               discrete = TRUE, curve_df = curve_df, color_lines = color_lines,
                                               color_points = color_points)}
  
  # ---------------- compared prediction ----------------
  
  # required column miQC.keep_adj_compared
  if (add_compared_prediction){
     md[["miQC.keep_adj_compared"]] <- factor(md[["miQC.keep_adj_compared"]], levels = c("both_keep", "both_discard", "keep_orig", "keep_geom"))
     plots[["prediction_comparison"]] <- miQC_build_plot(df = md, col_var = "miQC.keep_adj_compared",
                                                             complexity_column = complexity_column, percent_mt = "percent_mt",
                                                             title_txt = paste0(prefix, "Adjustment Comparison"),
                                                             color_title= "Keep Status",
                                                             point_alpha = point_alpha, point_size = point_size,
                                                             discrete = TRUE, curve_df = curve_df, color_lines = color_lines,
                                                             color_points = color_compare)}

  # ---------------- invalid probability ----------------
  
  # required column: p_invalid
  if(add_probability){
  md[["P(LowQuality)"]] <- md[[p_invalid]]
  
  plots[["probability"]] <- miQC_build_plot(df = md, col_var = "P(LowQuality)",
                                            complexity_column = complexity_column, percent_mt = "percent_mt",
                                            title_txt = paste0(prefix, "Low Quality Probability"),
                                            color_title= "P(LowQuality)",
                                            point_alpha = point_alpha, point_size = point_size,
                                            discrete = FALSE, curve_df = curve_df, color_lines = color_lines,
                                            color_points = color_points)}
  
  
  # ---------------- adjustments ----------------
  
  for (method in adjustment_method){
    adj_keep <- keep_adj[[method]]
    adj_geom <- geom_adj[[method]]
    adj_region <- geom_region[[method]]
    region_color <- color_regions[[method]]
    method_prefix <- gsub("_check", "", method)
    prefix2 <- paste0(toupper(substr(method_prefix, 1, 1)), substr(method_prefix, 2, nchar(method_prefix)), " ")
    
  # adjusted prediction
    # required columns: keep_adj & geom_adj
    md[[adj_geom]] <- factor(md[[adj_geom]], levels = c("unchanged", "geom_pass", "geom_fail"))
    md[[adj_keep]] <- factor(md[[adj_keep]], levels = c("discard", "keep"))
    
    plots[[paste0(method_prefix, "_adjusted_classification")]]  <- miQC_build_plot(df = md, col_var = adj_keep,
                                                      complexity_column = complexity_column, percent_mt = "percent_mt",
                                                      title_txt = paste0(prefix, prefix2, "Adjusted Classification"),
                                                      shape_var= adj_geom, color_title = "Classification",
                                                      point_alpha = point_alpha, point_size = point_size,
                                                      discrete = TRUE, curve_df = curve_df, color_lines = color_lines,
                                                      color_points = color_points, shape_geom= shape_geom)
    

  # geometry actions
  # required columns: keep_adj & geom_region
    if (method == "original_check"){
      md[[adj_region]] <- factor(md[[adj_region]], levels = c("unchanged", "below_boundary", "left_cutoff"))}
    if (method == "geometry_check"){
    md[[adj_region]] <- factor(md[[adj_region]], levels = c("valid", "converge", "invalid"))}
    
    plots[[paste0(method_prefix, "_region_adjustments")]]  <- miQC_build_plot(df = md, col_var = adj_region,
                                                      complexity_column = complexity_column, percent_mt = "percent_mt",
                                                      title_txt = paste0(prefix, prefix2, "Region-related Adjustments"),
                                                      shape_var = adj_geom, color_title= paste0(prefix2, "Regions"),
                                                      point_alpha = point_alpha, point_size = point_size,
                                                      discrete = TRUE, curve_df = curve_df, color_lines = color_lines,
                                                      color_points = region_color, shape_geom = shape_geom)}
  
  if (length(plots) == 0L){
    stop("No plots were created.")}
  
    if (isTRUE(combine_plots)){
        prefix3 <- if (isTRUE(add_unadjusted_prediction)) "Unmodified and " else ""
      plots <- patchwork::wrap_plots(plots, ncol=2) + patchwork::plot_annotation(title= paste0(global_prefix, prefix3, "Adjusted Cell Classification"))}
    
    return(plots)}

