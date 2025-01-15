# File: test-stepstable.R

library(testthat)

# Assuming stepstable and calculate_odds_ratios functions are defined in stepstable.R
# You might need to source the function if it's not part of a package
# source("path_to_stepstable.R")

# Sample Test Data
set.seed(123)  # For reproducibility
toy_data <- data.frame(
  padeg = factor(sample(1:3, 100, replace = TRUE)),
  degree = factor(sample(1:3, 100, replace = TRUE))
)

# Observed OD Table
observed_O <- as.matrix(table(toy_data$padeg, toy_data$degree))

# Function to Compute Independence Table
compute_independence_table <- function(O) {
  row_totals <- rowSums(O)
  col_totals <- colSums(O)
  total <- sum(O)
  I <- outer(row_totals, col_totals) / total
  return(I)
}

# Function to Compute Log Odds Ratios
compute_log_odds_ratios <- function(tbl) {
  or <- calculate_odds_ratios(tbl)
  # Take log, handle NA by returning NA
  log_or <- log(or)
  return(log_or)
}

# Test Suite for stepstable Function


test_that("stepstable returns the independence table when strength = 1", {
  strength <- 1
  counterfactual_O <- stepstable(toy_data, o = "padeg", d = "degree", strength = strength)

  # Compute independence table
  independence_I <- compute_independence_table(observed_O)

  # Check that counterfactual_O is equal to independence_I
  expect_equal(counterfactual_O, independence_I, tolerance = 1e-6)
})

test_that("stepstable with strength = 0.5 produces log odds ratios approximately halfway between observed and independence tables", {
  strength <- 0.5
  counterfactual_O <- stepstable(toy_data, o = "padeg", d = "degree", strength = strength)

  # Compute log odds ratios
  log_or_observed <- compute_log_odds_ratios(observed_O)
  log_or_counterfactual <- compute_log_odds_ratios(counterfactual_O)
  log_or_independence <- compute_log_odds_ratios(compute_independence_table(observed_O))

  # Calculate expected log odds ratios
  # Since log(a^0.5 * b^0.5) = 0.5*log(a) + 0.5*log(b)
  expected_log_or <- 0.5 * log_or_observed + 0.5 * log_or_independence

  # Compare log_or_counterfactual to expected_log_or
  # Allowing a small tolerance due to iterative fitting
  # We'll use expect_true with all(abs difference) < some threshold
  difference <- abs(log_or_counterfactual - expected_log_or)

  # Define a reasonable tolerance, e.g., 0.1
  tolerance <- 0.1

  # Check that all non-NA differences are within tolerance
  expect_true(all(is.na(difference) | difference < tolerance))
})




test_that("test a value other than 0.5, like 0.3", {
  strength <- 0.3
  counterfactual_O <- stepstable(toy_data, o = "padeg", d = "degree", strength = strength)

  # Compute log odds ratios
  log_or_observed <- compute_log_odds_ratios(observed_O)
  log_or_counterfactual <- compute_log_odds_ratios(counterfactual_O)
  log_or_independence <- compute_log_odds_ratios(compute_independence_table(observed_O))

  # Calculate expected log odds ratios
  # Since log(a^0.5 * b^0.5) = 0.5*log(a) + 0.5*log(b)
  expected_log_or <- 0.7 * log_or_observed + 0.3 * log_or_independence

  # Compare log_or_counterfactual to expected_log_or
  # Allowing a small tolerance due to iterative fitting
  # We'll use expect_true with all(abs difference) < some threshold
  difference <- abs(log_or_counterfactual - expected_log_or)

  # Define a reasonable tolerance, e.g., 0.1
  tolerance <- 0.1

  # Check that all non-NA differences are within tolerance
  expect_true(all(is.na(difference) | difference < tolerance))
})
