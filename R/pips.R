#' Estimating post-intervention propensity score
#'
#' @param data A dataframe
#' @param d The destination or the treatment variable, a string that specifies the name of the destination variable
#' @param o The origin variable, a string that specifies the name of the origin variable
#' @param x The confounds between Y and D except for o, a vector of strings containing control variable names
#' @param custom.table Custom post-intervention mobility table, a "table" object where rows are origin and destination column
#' @param strength Strength of intervention captured by proportional log odds-ratio change, 0 indicates no intervention and 1 indicates an intervention that sets origin-destination to be independent
#' @param contrast A hypothetical contrast table, a "table" object where the odds-ratios in the table will be extracted to form a interventional table
#'
#' @returns a dataframe that lists the propensity of receiving treatment for all indiduals, where the first column is probability of receiving destination category 1, the second column the probability of category 2, and so on
#' @export
#'
#' @examples
#' # Example data
#' toy <- data.frame(
#'   padeg = factor(rep(c("A", "B", "C"), each = 3)),
#'   degree = factor(c("X", "Y", "Z", "X", "Y", "Z", "X", "Y", "Z")),
#'   sibs = rnorm(9),
#'   ethnic = sample(c("group1", "group2"), 9, replace = TRUE)
#' )
#'
#' # Example usage with strength
#' pips(data = toy, o = "padeg", d = "degree", x = c("sibs", "ethnic"), strength = 1)
#'
#' # Example usage with a custom mobility table
#' custom_table <- matrix(
#'   c(0.3, 0.4, 0.3,
#'     0.2, 0.5, 0.3,
#'     0.1, 0.6, 0.3),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(c("A", "B", "C"), c("X", "Y", "Z"))
#' )
#' pips(data = toy, o = "padeg", d = "degree", x = c("sibs", "ethnic"), custom.table = custom_table)
#'
#' # Example usage with a contrast table
#' contrast_table <- matrix(
#'   c(1, 2, 1,
#'     1, 1, 2,
#'     2, 1, 1),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(c("A", "B", "C"), c("X", "Y", "Z"))
#' )
#' pips(data = toy, o = "padeg", d = "degree", x = c("sibs", "ethnic"), contrast = contrast_table)
pips <- function(data, o, d, x, custom.table = NULL, strength = NULL, contrast = NULL){

  # Load required package
  if(!requireNamespace("nnet", quietly = TRUE)) {
    stop("Package 'nnet' is required but not installed.")
  }

  # Validate that only one of custom.table, strength, contrast is provided
  params_provided <- sum(!is.null(custom.table), !is.null(strength), !is.null(contrast))
  if (params_provided > 1) {
    stop("Only one of `custom.table`, `strength`, or `contrast` must be specified.")
  }

  # If none are provided, set strength = 1
  if (params_provided == 0) {
    strength <- 1
  }

  # Check for missing values in o, d, and x
  required_vars <- unique(c(o, d, x))
  missing_vars <- sapply(required_vars, function(var) any(is.na(data[[var]])))
  if (any(missing_vars)) {
    vars_with_na <- names(missing_vars)[missing_vars]
    stop(paste("The following variables contain missing values:", paste(vars_with_na, collapse = ", ")))
  }

  # Internal function to compute delta
  intervention.delta <- function(intervention, observed, reference_col){
    # Adding a small value to avoid division by zero
    intervention <- intervention + 0.01
    observed <- observed + 0.01

    # Compute proportions based on reference column
    prop_int <- sweep(intervention, 2, intervention[, reference_col], FUN = "/")
    prop_obs <- sweep(observed, 2, observed[, reference_col], FUN = "/")

    # Compute delta as ratio of proportions
    delta <- prop_int / prop_obs

    return(delta)
  }

  # Estimate observed propensity scores using multinomial logistic regression
  formula <- as.formula(paste(d, "~", paste(c(x, o), collapse = " + ")))
  pscore_model <- nnet::multinom(formula, data = data, trace = FALSE)
  pscores <- predict(pscore_model, newdata = data, type = "probs")

  # Generate observed mobility table
  observed_table <- prop.table(table(data[[o]], data[[d]]), 2)  # Column proportions

  # Generate counterfactual table
  if (!is.null(custom.table)) {
    if (!is.matrix(custom.table)) {
      stop("`custom.table` must be a matrix.")
    }
    if (!all(rownames(custom.table) == levels(data[[o]])) || !all(colnames(custom.table) == levels(data[[d]]))) {
      stop("Row and column names of `custom.table` must match those of the observed table.")
    }
    counterfactual_table <- custom.table
  } else if (!is.null(contrast)) {
    counterfactual_table <- contrasttable(data, o, d, contrast)
  } else if (!is.null(strength)) {
    # stepstable function is defined externally
    counterfactual_table <- stepstable(data, o, d, strength)
  }

  # Ensure counterfactual_table is a matrix with same dimensions as observed
  if (!is.matrix(counterfactual_table)) {
    counterfactual_table <- as.matrix(counterfactual_table)
  }

  # Define reference column (e.g., first category)
  reference_col <- 1

  # Compute delta
  delta <- intervention.delta(counterfactual_table, observed_table, reference_col)

  # Assign row and column names to delta to match tables
  rownames(delta) <- rownames(observed_table)
  colnames(delta) <- colnames(observed_table)

  # Adjust propensity scores
  origins <- data[[o]]
  unique_origins <- levels(origins)

  # Initialize pips matrix
  pips_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(pscores))
  colnames(pips_matrix) <- paste0("pips.", colnames(pscores))

  for (orig in unique_origins) {
    idx <- which(origins == orig)
    pips_matrix[idx, ] <- pscores[idx, ] * delta[orig, ]
    # Normalize to sum to 1
    row_sums <- rowSums(pips_matrix[idx, , drop = FALSE])
    pips_matrix[idx, ] <- sweep(pips_matrix[idx, , drop = FALSE], 1, row_sums, FUN = "/")
  }

  # Convert pips_matrix to dataframe
  pips_df <- as.data.frame(pips_matrix)

  return(pips_df)
}
