# File: R/conditionaleffect.R

#' Regression-based estimators for conditional average treatment effect
#'
#' @description
#' This function fits a Super Learner (through \code{sl3}) to estimate the outcome model,
#' then predicts outcomes for each treatment level to derive conditional average
#' treatment effects (CATEs). It can optionally bootstrap these estimates.
#'
#' @param data A dataframe
#' @param y The outcome variable (string), numeric or factor
#' @param d The destination/treatment variable (string), must be a factor with >=2 levels
#' @param o The origin variable (string), must be a factor
#' @param x A vector of strings for additional control variables (confounders)
#' @param estimator A learner name or vector of learner names from \code{sl3}. If multiple,
#'   a Super Learner is built. Default is "Lrnr_glm".
#' @param boot_reps The number of bootstrap repetitions for estimating standard errors. Default is 500.
#'
#' @returns A list containing:
#' \item{predicted_outcomes}{A dataframe of predicted outcomes, each column corresponding to a treatment level.}
#' \item{cate_summary}{A dataframe of conditional average treatment effects (CATEs) with bootstrapped standard errors.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   result <- conditionaleffect(
#'     data = my_data,
#'     y = "outcome",
#'     o = "origin",
#'     d = "treatment",
#'     x = c("control1", "control2"),
#'     estimator = c("Lrnr_glm", "Lrnr_randomForest"),
#'     boot_reps = 1000
#'   )
#' }
conditionaleffect <- function(
    data,
    y,
    o,
    d,
    x,
    estimator = "Lrnr_glm",
    boot_reps = 500
) {
  # Load required packages
  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' is required but not installed.")
  }
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required but not installed.")
  }

  # 1. Input Validation
  required_vars <- unique(c(y, d, o, x))
  missing_vars  <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are missing from data: ",
         paste(missing_vars, collapse = ", "))
  }

  # Check for missing values
  vars_with_na <- required_vars[sapply(required_vars, function(var) any(is.na(data[[var]])))]
  if (length(vars_with_na) > 0) {
    stop("The following variables contain missing values: ",
         paste(vars_with_na, collapse = ", "))
  }

  # 2. Ensure 'd' and 'o' are factors
  data[[d]] <- as.factor(data[[d]])
  data[[o]] <- as.factor(data[[o]])

  # 3. Check that 'd' has at least two levels
  if (nlevels(data[[d]]) < 2) {
    stop("The treatment variable 'd' must have at least two levels.")
  }


  # 5. Build the sl3 learner (single or multiple)
  if (length(estimator) == 1) {
    sl <- sl3::make_learner(estimator)
  } else {
    # multiple => super learner
    stack <- do.call("make_learner_stack", as.list(estimator))
    sl <- sl3::make_learner("Lrnr_sl", learners = stack,
                            metalearner = sl3::make_learner("Lrnr_glm"))
  }

  # 6. Determine outcome_type from the actual class of data[[y]]
  #    to ensure sl3 respects the user-specified numeric/factor outcome
  outcome_vec <- data[[y]]
  if (is.numeric(outcome_vec)) {
    # treat as continuous
    outcome_type <- "continuous"
  } else if (is.factor(outcome_vec)) {
    if (nlevels(outcome_vec) == 2) {
      outcome_type <- "binomial"
    } else {
      # multi-level factor => "categorical"
      outcome_type <- "categorical"
    }
  } else {
    stop("Unsupported outcome class: only numeric or factor are allowed.")
  }

  # 7. Create a Task and fit the SuperLearner to the data
  task <- sl3::make_sl3_Task(
    data          = data,
    covariates    = c(d, o, x),
    outcome       = y,
    outcome_type  = outcome_type
  )
  sl_fit <- sl$train(task)

  # 8. Generate predictions for each treatment level
  treatment_levels <- levels(data[[d]])
  predicted_outcomes <- data.frame(matrix(
    NA, nrow = nrow(data), ncol = length(treatment_levels)
  ))
  colnames(predicted_outcomes) <- treatment_levels

  for (level in treatment_levels) {
    new_data <- data
    new_data[[d]] <- factor(level, levels = treatment_levels)
    new_task <- sl3::make_sl3_Task(
      data         = new_data,
      covariates   = c(d, o, x),
      outcome      = y,
      outcome_type = outcome_type
    )
    pred <- sl_fit$predict(new_task)
    predicted_outcomes[[level]] <- as.numeric(pred)
  }

  # 9. Calculate Conditional Average Treatment Effects (CATEs)
  pairwise_combinations <- combn(treatment_levels, 2)
  cate_list <- list()

  for (i in seq_len(ncol(pairwise_combinations))) {
    treat1 <- pairwise_combinations[1, i]
    treat2 <- pairwise_combinations[2, i]
    cate_name <- paste0("CATE_", treat1, "_vs_", treat2)
    cate_list[[cate_name]] <- predicted_outcomes[[treat1]] - predicted_outcomes[[treat2]]
  }
  cate_df <- as.data.frame(cate_list)

  # Example aggregator
  cate_df[[o]] <- data[[o]]
  cate_df <- cate_df %>%
    dplyr::group_by(across(all_of(o))) %>%
    dplyr::summarise(across(everything(), mean)) %>%
    tidyr::pivot_longer(
      cols = starts_with("CATE"),
      names_to = "NAME",
      values_to = "Estimate"
    )

  # 10. Bootstrap for Standard Errors
  # Define a function to compute CATEs on bootstrap samples
  boot_cate <- function(dframe, indices) {
    boot_data <- dframe[indices, ]

    # re-create the outcome_type logic here
    boot_outcome_type <- outcome_type
    # or optionally re-check is.numeric(...) in boot_data[[y]]

    boot_task <- sl3::make_sl3_Task(
      data         = boot_data,
      covariates   = c(d, o, x),
      outcome      = y,
      outcome_type = boot_outcome_type
    )
    boot_sl_fit <- sl$train(boot_task)

    # Predict outcomes for each treatment level in the boot sample
    predicted_outcomes_boot <- data.frame(matrix(
      NA, nrow = nrow(dframe), ncol = length(treatment_levels)
    ))
    colnames(predicted_outcomes_boot) <- treatment_levels

    for (level in treatment_levels) {
      new_data_boot <- boot_data
      new_data_boot[[d]] <- factor(level, levels = treatment_levels)
      new_task_boot <- sl3::make_sl3_Task(
        data         = new_data_boot,
        covariates   = c(d, o, x),
        outcome      = y,
        outcome_type = boot_outcome_type
      )
      pred_boot <- boot_sl_fit$predict(new_task_boot)
      if (is.list(pred_boot)) {
        pred_boot <- unlist(pred_boot)
      }
      predicted_outcomes_boot[[level]] <- as.numeric(pred_boot)
    }

    # Calculate CATEs for the bootstrap sample
    cate_boot <- list()
    for (i in seq_len(ncol(pairwise_combinations))) {
      treat1 <- pairwise_combinations[1, i]
      treat2 <- pairwise_combinations[2, i]
      cname <- paste0("CATE_", treat1, "_vs_", treat2)
      cate_boot[[cname]] <- predicted_outcomes_boot[[treat1]] -
        predicted_outcomes_boot[[treat2]]
    }

    boot_df <- as.data.frame(cate_boot)
    boot_df[[o]] <- boot_data[[o]]
    boot_df <- boot_df %>%
      dplyr::group_by(across(all_of(o))) %>%
      dplyr::summarise(across(everything(), mean)) %>%
      tidyr::pivot_longer(
        cols = starts_with("CATE"),
        names_to = "NAME",
        values_to = "Estimate"
      )

    return(boot_df$Estimate)
  }

  # Perform the bootstrap
  set.seed(123)  # For reproducibility
  boot_results <- boot::boot(
    data     = data,
    statistic= boot_cate,
    R        = boot_reps
  )

  # Calculate bootstrap standard errors
  cate_se <- apply(boot_results$t, 2, sd, na.rm = TRUE)

  # Attach SE to the final data frame
  cate_df$SE <- cate_se

  # Return the results
  list(
    predicted_outcomes = predicted_outcomes,
    cate_summary       = cate_df
  )
}
