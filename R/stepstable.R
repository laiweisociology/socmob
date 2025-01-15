#' Creating counterfactual mobility table based on specified interventional strength
#'
#' @param data Data frame containing the origin and destination variables
#' @param o Origin variable name, a string
#' @param d Destination variable name, a string
#' @param strength Interventional strength as defined by proportional logged odds-ratio distance from the independence table, range from 0-1
#'
#' @returns A counterfactual mobility table as implied by the interventional strength
#' @export
#'
#' @examples
#' toy_data <- data.frame(
#'   padeg = factor(sample(1:3, 100, replace = TRUE)),
#'   degree = factor(sample(1:3, 100, replace = TRUE))
#' )
#' stepstable(toy_data, o = "padeg", d = "degree", strength = 1)
#'
stepstable <- function(data, o, d, strength) {
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


  if (!is.numeric(strength) || strength < 0 || strength > 1) {
    stop("`strength` must be a numeric value between 0 and 1.")
  }

  # Create observed OD table
  O <- as.matrix(table(data[[o]], data[[d]]))

  # Add 0.01 to zero cells to prevent issues with zero counts
  zero_cells <- which(O == 0, arr.ind = TRUE)
  if (nrow(zero_cells) > 0) {
    O[zero_cells] <- O[zero_cells] + 0.01
  }

  # Compute marginal totals
  row_totals <- rowSums(O)
  col_totals <- colSums(O)
  total <- sum(O)

  # Create independence table I
  I <- outer(row_totals, col_totals) / total

  # Handle zero cells in independence table to avoid division by zero or log issues
  if(any(I == 0)) {
    stop("Independence table has zero cells after adjustment, cannot proceed.")
  }

  # Compute initial counterfactual table using geometric mean interpolation
  # C_initial = I^(strength) * O^(1 - strength)
  C_initial <- I^strength * O^(1 - strength)

  # Replace any NaN or Inf with small positive number to avoid issues
  C_initial[!is.finite(C_initial)] <- 1e-10

  # Iterative Proportional Fitting (IPF) to adjust C_initial to match O's margins
  # Initialize
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








