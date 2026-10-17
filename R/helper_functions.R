#' Internal helper functions
#' @keywords internal
miQC_label_cutoff <- function(x){ 
  sprintf("%.2f", round(x,2))}

#' @keywords internal
miQC_is_valid_color <- function(x){
  all(vapply(x, function(y){
    !inherits(try(grDevices::col2rgb(y), silent = TRUE), "try-error")}, logical(1)))}

#' @keywords internal
utils::globalVariables(c(
  ".data", "feature", "keep_y", "discard_y"))

#' @keywords internal
miQC_check_col <- function(col, col_names){
  length(col) != 1L || is.null(col) || is.na(col) || !col %in% col_names}

# helper plotting function
#' @keywords internal
miQC_build_plot <- function(df, col_var, title_txt,
                       complexity_column = "nFeature_RNA", percent_mt = "percent.mt", 
                       discrete = TRUE, curve_df = NULL, shape_var = NULL,
                       color_title=NULL, ylim_max = NULL,
                       point_alpha = 0.8, point_size = 1, 
                       color_lines = c(discard = "black", keep = "blue"),
                       color_points = c(discard = "dimgrey", keep = "cornflowerblue"),
                       shape_geom = NULL) {
  
  if (!inherits(df, "data.frame")) {
    stop("`df` must be a data.frame object.")}
  
  required_cols <- c(col_var, complexity_column, percent_mt, shape_var)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L){
    stop("`df` is missing required columns: ", paste(missing_cols, collapse = ", "))}
  
  if (!is.numeric(point_alpha) ||
      length(point_alpha) != 1L ||
      !is.finite(point_alpha) ||
      point_alpha < 0 ||
      point_alpha > 1) {
    stop("`point_alpha` must be a numeric value between 0 and 1.")}
  
  if (!is.numeric(point_size) ||
      length(point_size) != 1L ||
      !is.finite(point_size) ||
      point_size < 0) {
    stop("`point_size` must be a non-negative numeric value.")}
  
  if (is.null(ylim_max)){
      ylim_max <- as.numeric(max(df[[percent_mt]], na.rm = TRUE) + 0.01)}
  
  if (!is.numeric(ylim_max) ||
      length(ylim_max) != 1L ||
      !is.finite(ylim_max) ||
      ylim_max <= 0) {
    stop("`ylim_max` must be a positive finite numeric scalar.")}
  
  if (!is.null(color_title)){
    color_title <- as.character(color_title)
    if (length(color_title) != 1L ||
        is.na(color_title) ||
        !nzchar(color_title)) {
      stop("`color_title` must be `NULL` or a non-empty character scalar.")}}

    if (!is.logical(discrete) || length(discrete) != 1L || is.na(discrete)) {
      stop("`discrete` must be one non-missing logical value.")}
  
  # ---------------- create plots ----------------
  
  # validate colors
  if (isTRUE(discrete)){
    point_vals <- unique(as.character(stats::na.omit(df[[col_var]])))
    missing_cols <- setdiff(point_vals, names(color_points))
    if (length(missing_cols) > 0L) {
      stop("Missing discrete color specification for `color_points`: ",paste(missing_cols, collapse= ", "))}
    color_points <- color_points[names(color_points) %in% point_vals]
    valid_cols <- vapply(unname(color_points), miQC_is_valid_color, logical(1))
    if (!all(valid_cols)) {
      stop("Invalid discrete color specification for `color_points`.")}}
  
  # validate shapes
  if (!is.null(shape_var)){
    geom_vals <- unique(as.character(stats::na.omit(df[[shape_var]])))
    missing_shapes <- setdiff(geom_vals, names(shape_geom))
    if (length(missing_shapes) > 0L) {
      stop("Missing shape specification for `shape_geom`: ",paste(missing_shapes, collapse= ", "))}
    selected_shapes <- shape_geom[geom_vals]
    if (!is.numeric(selected_shapes) ||
        any(!is.finite(selected_shapes)) ||
        any(selected_shapes < 0 | selected_shapes > 25)) {
      stop("`shape_geom` must contain valid numeric ggplot2 shape codes.")}}
  
  if (!isTRUE(discrete)) {
    # check continuous color variable
    df[[col_var]] <- as.numeric(df[[col_var]])
    invalid_col <- !is.finite(df[[col_var]])
    if (any(invalid_col)) {
      stop("Found ", sum(invalid_col)," cells with invalid values for `", col_var, "`.")}}
  
  # create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[complexity_column]],
                                        y = .data[[percent_mt]],
                                        colour = .data[[col_var]]))

  if (!is.null(shape_var)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_var]]),
                                 size = point_size, alpha = point_alpha) +
      ggplot2::scale_shape_manual(values = selected_shapes)
  }else{
    p <- p + ggplot2::geom_point(size = point_size, alpha = point_alpha)}
  
  p <- p + ggplot2::labs(x = complexity_column,
                         y = percent_mt,
                         title = title_txt,
                         colour = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, ylim_max)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 11),
                   legend.text  = ggplot2::element_text(size = 10))
  
  if (!is.null(curve_df)) {
    if (!inherits(curve_df, "data.frame")) {
      stop("`curve_df` must be a data.frame object.")}
    
    required_cols <- c("feature", "keep_y", "discard_y")
    missing_cols <- setdiff(required_cols, colnames(curve_df))
    if (length(missing_cols) > 0L){
      stop("curve_df is missing required columns: ", paste(missing_cols, collapse = ", "))}
    
    # check line colors
    required_names <- c("keep", "discard")
    valid_lines <- !is.null(color_lines) &&
      all(required_names %in% names(color_lines)) &&
      all(vapply(unname(color_lines[required_names]), miQC_is_valid_color, logical(1)))
      
    if (!valid_lines) {
      message("Invalid line colours. Using defaults.")
      color_lines <- c(discard = "black", keep = "blue")}
    
    p <- p + ggplot2::geom_line(data = curve_df,
        ggplot2::aes(x = feature,
                     y = keep_y,
                     linetype = "keep"),
        inherit.aes = FALSE,
        linewidth = 1,
        colour = unname(color_lines["keep"])) +
      
      ggplot2::geom_line(data = curve_df,
        ggplot2::aes(x = feature,
                     y = discard_y,
                     linetype = "discard"),
        inherit.aes = FALSE,
        linewidth = 1,
        colour = unname(color_lines["discard"])) +
      
      ggplot2::scale_linetype_manual(
        name = "Fitted components",
        breaks = c("discard", "keep"),
        values = c(discard = "solid", keep = "solid"),
        labels = c(discard = "Discard", keep = "Keep")) +
      
      ggplot2::guides(
        linetype = ggplot2::guide_legend(
          override.aes = list(
            colour = unname(color_lines[c("discard", "keep")]),
            linewidth = 1)))}
  
  if (isTRUE(discrete)) {
    p <- p + ggplot2::scale_colour_manual(values = color_points, name=color_title)
  }else{
    p <- p + ggplot2::scale_colour_viridis_c(limits = c(0, 1), option = "D", name=color_title, na.value = "lightgrey")}
  
  return(p)}
