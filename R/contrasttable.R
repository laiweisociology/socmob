#' Creating counterfactual mobility table based on a comparison case
#'
#' @param data Data frame containing the origin and destination variables
#' @param o Origin variable name, a string
#' @param d Destination variable name, a string
#' @param contrast A comparison mobility table with the same origin and destination dimensions
#'
#' @returns A counterfactual mobility table with odds ratios from the comparison case and marginal O and D distributions from data
#' @export
#'
#' @examples
#' toy_data <- data.frame(
#'   padeg = factor(sample(1:3, 100, replace = TRUE)),
#'   degree = factor(sample(1:3, 100, replace = TRUE))
#' )
#'
#' comparison_data <- data.frame(
#'   padeg = factor(sample(1:3, 100, replace = TRUE)),
#'   degree = factor(sample(1:3, 100, replace = TRUE))
#' )
#'
#'
#'comparison_table <- table(comparison_data$padeg,comparison_data$degree)
#'
#' contrasttable(toy_data, o = "padeg", d = "degree", contrast = comparison_table)
#'

contrasttable <- function(data, o, d, contrast) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  if (!all(c(o, d) %in% names(data))) {
    stop("Origin (`o`) and destination (`d`) variables must be present in the data.")
  }

  # Check if 'o' and 'd' are factors
  if (!is.factor(data[[o]])) {
    stop(paste("Origin variable `", o, "` must be a factor.", sep = ""))
  }

  if (!is.factor(data[[d]])) {
    stop(paste("Destination variable `", d, "` must be a factor.", sep = ""))
  }

  # Check for missing values in 'o' and 'd'
  missing_o <- sum(is.na(data[[o]]))
  missing_d <- sum(is.na(data[[d]]))

  if (missing_o > 0 || missing_d > 0) {
    stop(paste(
      "Missing values detected:",
      if (missing_o > 0) paste(missing_o, "in origin variable `", o, "`.", sep = "") else "",
      if (missing_d > 0) paste(missing_d, "in destination variable `", d, "`.", sep = "") else "",
      sep = " "
    ))
  }

  # Create observed OD table
  O <- as.matrix(table(data[[o]], data[[d]]))

  # Add 0.01 to zero cells to prevent issues with zero counts
  O[O == 0] <- 0.01

  # Validate contrast table
  if (!is.matrix(contrast) && !is.table(contrast)) {
    stop("`contrast` must be a matrix or table.")
  }

  # Check if contrast table has the same number of rows and columns as observed table
  if (nrow(contrast) != nrow(O)) {
    stop("Contrast table must have the same number of rows as the observed table.")
  }

  if (ncol(contrast) != ncol(O)) {
    stop("Contrast table must have the same number of columns as the observed table.")
  }

  # Check if contrast table has the same row and column names as observed table
  if (!all(rownames(contrast) == rownames(O))) {
    stop("Row names of `contrast` table must match row names of the observed table.")
  }

  if (!all(colnames(contrast) == colnames(O))) {
    stop("Column names of `contrast` table must match column names of the observed table.")
  }

  # Convert contrast to matrix if it's a table
  if (is.table(contrast)) {
    contrast <- as.matrix(contrast)
  }

  # Add 0.01 to zero cells in contrast table
  contrast[contrast == 0] <- 0.01

  # Compute proportions from contrast table
  contrast_prop <- contrast / sum(contrast)

  # Initialize counterfactual table with proportions scaled to observed grand total
  grand_total <- sum(O)
  C_initial <- contrast_prop * grand_total

  # Apply IPF to adjust C_initial to match observed marginals
  row_totals <- rowSums(O)
  col_totals <- colSums(O)

  C_cf <- C_initial
  max_iter <- 1000
  tol <- 1e-6
  for (iter in 1:max_iter) {
    # Adjust rows
    row_scaling <- row_totals / rowSums(C_cf)
    C_cf <- sweep(C_cf, 1, row_scaling, FUN = "*")

    # Adjust columns
    col_scaling <- col_totals / colSums(C_cf)
    C_cf <- sweep(C_cf, 2, col_scaling, FUN = "*")

    # Check convergence
    row_diff <- max(abs(rowSums(C_cf) - row_totals))
    col_diff <- max(abs(colSums(C_cf) - col_totals))
    if (row_diff < tol && col_diff < tol) {
      break
    }
    if (iter == max_iter) {
      warning("IPF did not converge within the maximum number of iterations.")
    }
  }

  # Return the counterfactual table
  return(C_cf)
}
