#' Adjusted geometric validation of selected miQC classifications
#'
#' Performs a post hoc geometry-based plausibility check of a selected
#' miQC classification. Cells may be reassigned when their current
#' keep/discard status is inconsistent with the fitted model curves.
#'
#' The function adds adjusted classifications and diagnostic measures to the
#' provided `data.frame` and records thresholds and the number of reassigned 
#' cells in `geometry_check`.
#'
#' @param df The `data.frame` object for `flexmix_data`.
#' @param model_type Character Scalar describing the model type used in 
#' `run_miQCadj()`. Available options include \code{"one_dimensional"},
#'  \code{"linear"}, \code{"spline"}, \code{"polynomial"}, or 
#'   \code{"polynomial_raw"}. 
#' @param rescue_quantile Numeric between 0 and 1. Lower threshold used for
#'  rescuing discarded cells. Will be ignored for one-dimensional model.
#'  Defaults to \code{0.5}.
#' @param discard_quantile Numeric between 0 and 1. Upper threshold used for 
#'  reassigning implausible kept cells. Will be ignored for one-dimensional 
#'  model. Defaults to \code{0.95}.
#' @param converge_quantile Numeric between 0 and 1. Lower quantile of curve 
#'  gaps used to define convergence regions. Will be ignored for one-dimensional 
#'  model. Defaults to \code{0.05}.
#' @param safe_keep_dist Non-negative Numeric Scalar. Distance ratio favoring 
#'  keeping cells. Will be ignored for one-dimensional model. Defaults to 
#'  \code{1.25}).
#' @param verbose Logical Scalar. Whether information regarding reassignments
#'  should be printed. Defaults to \code{TRUE}.
#' 
#' @return Named list with slot information for:
#' \describe{
#' \item{df}{Updated `data.frame` with geometry-adjusted `miQC.keep_geom_adj` 
#'  classification column and geometry-related columns `miQC.geom_adj`, 
#' `miQC.geom_region`, `miQC.dist_up`, and `miQC.dist_low`.}
#'  \item{geometry_check}{Adjustment thresholds used for geometry checks,
#'  including the specified function arguments.} 
#' }
#'
#' @examples
#' \dontrun{
#' df_list <- geometry_adjustments(df=flexmix_data, model_type = "linear")
#' }
#'
#' @export
geometry_adjustments <- function(df, model_type, 
                               rescue_quantile=0.5, discard_quantile=0.95, 
                               converge_quantile=0.05, safe_keep_dist=1.25,
                               verbose = TRUE){
  
  #-------------------------------------------------------#
  
  n_min_check <- 10
  
  model_type <- match.arg(model_type, choices = c("one_dimensional", "linear", "spline", "polynomial", "polynomial_raw"))
  
  # checks
  if (!is.logical(verbose) ||
      length(verbose) != 1 ||
      is.na(verbose)) {
    stop("`verbose` must be one non-missing logical value.")}
  
  # df
  if (!inherits(df, "data.frame")) {
    stop("`df` must be a data.frame object.")}
  
  orig_df <- df
  if (anyDuplicated(rownames(orig_df)) || anyNA(rownames(orig_df))){
    stop("Unique, non-NA row identifiers required.")}
  
  # required df columns
  keep_col <- "miQC.keep"
  keep_comp <- "miQC.keep_comp"
  discard_comp <- "miQC.discard_comp"
  
  # required columns
  complexity_column <- if (model_type == "one_dimensional") NULL else "complexity_column"
  required_cols <- c("percent_mt", complexity_column, keep_col, keep_comp, discard_comp)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))}
  
  # ensure numeric values
  numeric_cols <- c("percent_mt", complexity_column, keep_comp, discard_comp)
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
    stop("Got an unexpected value for `miQC.keep`: ",unexpected, ". 
         \nExpected values are keep/discard.")}
  
  # subset to classified cells
  df <- df[keep_vals %in% c("keep", "discard"), , drop=FALSE]
  
  # get keep & discard group
  keep_idx <- df[[keep_col]] == "keep"
  discard_idx <- df[[keep_col]] == "discard"
  
  if (nrow(df) == 0L) {
    stop("No classified cells available for geometry adjustment.")}
  
  #---------------geometry check -----------------#
  
  # initialize thresholds
  threshold_invalid_region_keep <- NULL
  threshold_invalid_region_discard <- NULL
  
  if (model_type != "one_dimensional"){
    
    if (nrow(df) < n_min_check) {
      stop("Found only ",nrow(df)," non-excluded keep/discard cells. 
           \nAt least ",n_min_check," cells are required for geometry check.")}
    
    # provided quantiles
    check_quantiles <- list(rescue_quantile = rescue_quantile,
                            discard_quantile = discard_quantile,
                            converge_quantile = converge_quantile)
    
    lapply(names(check_quantiles), function(x){
      check_x <- check_quantiles[[x]]
      if (!is.numeric(check_x) ||
          length(check_x) != 1L ||
          !is.finite(check_x) ||
          check_x < 0 ||
          check_x > 1) {
        stop("`",x,"` must be a numeric value between 0 and 1.")}})
    
    if (!is.numeric(safe_keep_dist) ||
        length(safe_keep_dist) != 1L ||
        !is.finite(safe_keep_dist) ||
        safe_keep_dist < 0) {
      stop("`safe_keep_dist` must be a non-negative numeric value.")}
  
    ###############################################
    
    # run geometry check
    group_keep <- if (sum(keep_idx) >= 10) keep_idx else rep(TRUE, nrow(df))
    group_discard <- if (sum(discard_idx) >= 10) discard_idx else rep(TRUE, nrow(df))
    
    # aggregate classification by median along features
    curve_df <- data.frame(
      feature = df[["complexity_column"]],
      keep_y = df[[keep_comp]],
      discard_y = df[[discard_comp]])
    curve_df <- curve_df[order(curve_df$feature), , drop = FALSE]
    curve_df <- stats::aggregate(
      cbind(keep_y, discard_y) ~ feature,
      data = curve_df,
      FUN = stats::median)
    
    if (nrow(curve_df) < 2L) {
      stop("At least two distinct complexity values are required for geometry validation.")}
    
    f_keep <- stats::approxfun(curve_df$feature, curve_df$keep_y, rule = 2)
    f_discard <- stats::approxfun(curve_df$feature, curve_df$discard_y, rule = 2)
    
    # get curves
    keep_curve <- f_keep(df[["complexity_column"]])
    discard_curve <- f_discard(df[["complexity_column"]])
    
    upper_curve <- pmax(keep_curve, discard_curve)
    lower_curve <- pmin(keep_curve, discard_curve)
    
    mid_curve <- (upper_curve+lower_curve)/2
    
    # define valid region (discard > keep)
    invalid_region <- discard_curve <= keep_curve
    valid_region <- discard_curve > keep_curve
    
    # define converge region & add to invalid region
    gap <- abs(upper_curve - lower_curve)
    valid_gap <- gap[valid_region]
    
    if (length(valid_gap) < 5) {
      message("Less than 5 cells in valid_region. Convergence threshold estimated from all gaps.")
      converge_val <- stats::quantile(gap, converge_quantile, na.rm = TRUE)
    } else {
      converge_val <- stats::quantile(valid_gap, converge_quantile, na.rm = TRUE)}
    
    converge_region <- gap <= converge_val
    valid_region <- (valid_region & !converge_region)
    invalid_region <- (invalid_region | converge_region)
    
    # define distances
    dist_up <- abs(df[["percent_mt"]] - upper_curve)
    dist_low <- abs(df[["percent_mt"]] - lower_curve)
    
    # require distances for redefining cell classification
    safe_keep <- (dist_low*safe_keep_dist) < dist_up  
    
    # cell localization
    below_keep <-  df[["percent_mt"]] < keep_curve
    below_low <-  df[["percent_mt"]] < lower_curve
    
    # define tolerance based on related cell pool
    vertical_up <- df[["percent_mt"]] - upper_curve
    tol_up_val <- if (sum(group_keep & is.finite(vertical_up)) >= 5) stats::quantile(vertical_up[group_keep & is.finite(vertical_up)], rescue_quantile, na.rm = TRUE) else stats::quantile(vertical_up[is.finite(vertical_up)], rescue_quantile, na.rm = TRUE)
    too_high_val <- if(sum(group_discard & is.finite(vertical_up)) >= 5) stats::quantile(vertical_up[group_discard & is.finite(vertical_up)], discard_quantile, na.rm = TRUE) else stats::quantile(vertical_up[is.finite(vertical_up)], discard_quantile, na.rm = TRUE)
    below_upper_tol <- vertical_up <= tol_up_val
    above_upper_tol <- vertical_up > too_high_val
    
    vertical_up_mid <- df[["percent_mt"]] - mid_curve
    tol_up_mid_val <- if (sum(group_keep & is.finite(vertical_up_mid)) >= 5) stats::quantile(vertical_up_mid[group_keep & is.finite(vertical_up_mid)], rescue_quantile, na.rm = TRUE) else stats::quantile(vertical_up_mid[is.finite(vertical_up_mid)], rescue_quantile, na.rm = TRUE) 
    too_high_mid_val <- if(sum(group_discard & is.finite(vertical_up_mid)) >= 5) stats::quantile(vertical_up_mid[group_discard & is.finite(vertical_up_mid)], discard_quantile, na.rm = TRUE) else stats::quantile(vertical_up_mid[is.finite(vertical_up_mid)], discard_quantile, na.rm = TRUE)
    below_upper_mid_tol <- vertical_up_mid <= tol_up_mid_val
    above_upper_mid_tol <- vertical_up_mid > too_high_mid_val
    
    if (model_type == "linear"){  
      # values above both lines (exclude very high mt values)  
      valid_cell <- ((valid_region & below_keep) | (invalid_region & below_upper_tol))
      invalid_cell <-  (invalid_region & above_upper_tol)
      threshold_invalid_region_keep <- tol_up_val
      threshold_invalid_region_discard <- too_high_val
    }else{    
      valid_cell <- ((valid_region & (below_low | safe_keep))  | (invalid_region & below_upper_mid_tol))
      invalid_cell <- (invalid_region & above_upper_mid_tol)
      threshold_invalid_region_keep <- tol_up_mid_val
      threshold_invalid_region_discard <- too_high_mid_val}
    
  }else if (model_type == "one_dimensional"){
    rescue_quantile <- NULL
    discard_quantile <- NULL
    converge_quantile <- NULL
    safe_keep_dist <- NULL
    
    # get lines
    keep_curve <- stats::median(df[[keep_comp]])
    discard_curve <- stats::median(df[[discard_comp]])
    if (keep_curve == discard_curve){
      stop("Identical components for keep and discard lines in one-dimensional model.")}
    if (keep_curve > discard_curve) {
      stop("Keep component exceeds discard component in one-dimensional model.")}
    
    lower_curve <- min(keep_curve, discard_curve)
    upper_curve <- max(keep_curve, discard_curve)
    
    # basic check
    valid_cell <- df[["percent_mt"]] <= lower_curve
    invalid_cell <- df[["percent_mt"]] > upper_curve
    
    # define distances
    dist_up <- abs(df[["percent_mt"]] - upper_curve)
    dist_low <- abs(df[["percent_mt"]] - lower_curve)
    
    # regions
    converge_region <- rep(FALSE, nrow(df))
    valid_region <- rep(TRUE, nrow(df))}
  
    
  ###############################################
  
  # set NA FALSE
  valid_cell[is.na(valid_cell)] <- FALSE
  invalid_cell[is.na(invalid_cell)] <- FALSE
  
  # override discard -> rescue discarded cells
  override_discard <- discard_idx & valid_cell
  
  # override keep -> reassign suspicious kept cells
  override_keep <- keep_idx & invalid_cell

  # add to df
  df[["miQC.keep_geom_adj"]] <- as.character(df[["miQC.keep"]])
  df[["miQC.keep_geom_adj"]][override_keep] <- "discard"
  df[["miQC.keep_geom_adj"]][override_discard] <- "keep"
  
  df[["miQC.geom_adj"]] <- rep("unchanged", nrow(df))
  df[["miQC.geom_adj"]][override_keep] <- "geom_fail"
  df[["miQC.geom_adj"]][override_discard] <- "geom_pass"
  
  df[["miQC.geom_region"]] <- rep(NA_character_, nrow(df))
  df[["miQC.geom_region"]][converge_region] <- "converge"
  df[["miQC.geom_region"]][!converge_region & valid_region] <- "valid"
  df[["miQC.geom_region"]][!converge_region & !valid_region] <- "invalid"
  
  df[["miQC.dist_up"]] <- dist_up
  df[["miQC.dist_low"]] <- dist_low
  
  if (isTRUE(verbose)){
    message(paste("Updated geometry-related Reassignments: ",
                  paste(utils::capture.output(
                    table(df[,c("miQC.geom_region","miQC.geom_adj")], useNA = "ifany")), 
                    collapse = "\n"), sep = "\n"))}
  
  
  new_cols <- setdiff(colnames(df), colnames(orig_df))
  df <- df[rownames(orig_df), new_cols, drop = FALSE]
  orig_df <- cbind(orig_df, df)
  
return(list(df = orig_df,
            geometry_check = list(adjustment_method= "adjusted",
                                 model_type = model_type,
                                 rescue_quantile = rescue_quantile,
                                 discard_quantile = discard_quantile,
                                 converge_quantile = converge_quantile,
                                 safe_keep_dist = safe_keep_dist,
                                 threshold_invalid_region_keep = threshold_invalid_region_keep,
                                 threshold_invalid_region_discard = threshold_invalid_region_discard,
                                 n_reassigned_total = sum(override_keep) + sum(override_discard),
                                 n_override_keep = sum(override_keep),
                                 n_override_discard = sum(override_discard),
                                 n_converge_region = sum(converge_region),
                                 n_valid_region = sum(df[["miQC.geom_region"]] == "valid"),
                                 n_invalid_region = sum(df[["miQC.geom_region"]] == "invalid"))))}

