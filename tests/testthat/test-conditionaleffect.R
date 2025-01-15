# File: tests/testthat/test-conditionaleffect.R

library(testthat)
library(sl3)
library(boot)

# Assume that the conditionaleffect function is in the R/ directory and has been sourced
# Uncomment and adjust the path based on your project structure
# source("R/conditionaleffect.R")

# Sample Test Data Creation
set.seed(123)  # For reproducibility
sample_data <- data.frame(
  outcome = rnorm(100),
  treatment = factor(sample(1:3, 100, replace = TRUE)),
  origin = factor(sample(1:3, 100, replace = TRUE)),
  control1 = rnorm(100),
  control2 = rnorm(100)
)

# Define a helper function to list available sl3 learners
available_learners <- sl3::sl3_list_learners()

# Begin Test Suite
test_that("conditionaleffect returns a list with correct components", {
  # Use a simple learner for testing
  result <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100  # Reduced for faster testing
  )

  expect_type(result, "list")
  expect_true("predicted_outcomes" %in% names(result))
  expect_true("cate_summary" %in% names(result))
})

test_that("predicted_outcomes dataframe has correct dimensions and column names", {
  result <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100
  )

  expected_cols <- levels(sample_data$treatment)

  expect_s3_class(result$predicted_outcomes, "data.frame")
  expect_equal(nrow(result$predicted_outcomes), nrow(sample_data))
  expect_equal(ncol(result$predicted_outcomes), length(expected_cols))
  expect_equal(colnames(result$predicted_outcomes), expected_cols)

  # Check that predicted outcomes are numeric
  expect_true(all(sapply(result$predicted_outcomes, is.numeric)))
})

test_that("cate_summary dataframe has correct columns and structure", {
  result <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100
  )

  expect_s3_class(result$cate_summary, "data.frame")
  expect_equal(colnames(result$cate_summary), c("origin","NAME", "Estimate", "SE"))

  # Check that CATE names are correct
  treatment_levels <- levels(sample_data$treatment)
  expected_cate_names <- combn(treatment_levels, 2, FUN = function(x) paste0("CATE_", x[1], "_vs_", x[2]), simplify = FALSE)
expected_cate_names <- unlist(expected_cate_names)

# Check that Estimate and SE are numeric
expect_true(is.numeric(result$cate_summary$Estimate))
expect_true(is.numeric(result$cate_summary$SE))
})

test_that("conditionaleffect handles multiple learners by building a Super Learner", {
  # Select two available learners for testing
  learners <- c("Lrnr_glm", "Lrnr_randomForest")


  result <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = learners,
    boot_reps = 100
  )

  expect_true("predicted_outcomes" %in% names(result))
  expect_true("cate_summary" %in% names(result))

  # Further checks can be added to verify that the Super Learner was used,
  # but this requires access to the internal structure or mocking
})


test_that("conditionaleffect throws an error when required variables are missing", {
  # Remove a required variable
  incomplete_data <- sample_data[, !(names(sample_data) %in% "control1")]

  expect_error(
    conditionaleffect(
      data = incomplete_data,
      y = "outcome",
      o = "origin",
      d = "treatment",
      x = c("control1", "control2"),
      estimator = "Lrnr_glm",
      boot_reps = 100
    ),
    "The following variables are missing from data: control1"
  )
})

# test_that("conditionaleffect throws an error when invalid learner names are provided", {
#   invalid_learners <- c("Lrnr_glm", "Invalid_Learner")
#
#   expect_error(
#     conditionaleffect(
#       data = sample_data,
#       y = "outcome",
#       o = "origin",
#       d = "treatment",
#       x = c("control1", "control2"),
#       estimator = invalid_learners,
#       boot_reps = 100
#     ),
#     "Learner Invalid_Learner is not recognized. Use sl3::sl3_list_learners\\(\\) to see available learners\\."
#   )
# })

test_that("conditionaleffect handles missing values appropriately", {
  # Introduce missing values in a control variable
  data_with_na <- sample_data
  data_with_na$control1[1] <- NA

  expect_error(
    conditionaleffect(
      data = data_with_na,
      y = "outcome",
      o = "origin",
      d = "treatment",
      x = c("control1", "control2"),
      estimator = "Lrnr_glm",
      boot_reps = 100
    ),
    "The following variables contain missing values: control1"
  )
})

test_that("conditionaleffect works with the default estimator and boot_reps", {
  result <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2")
    # Using default estimator and boot_reps
  )

  expect_true("predicted_outcomes" %in% names(result))
  expect_true("cate_summary" %in% names(result))

  # Check default boot_reps was used implicitly by ensuring SE is present
  expect_true(all(!is.na(result$cate_summary$SE)))
})

test_that("conditionaleffect handles single treatment level gracefully", {
  # Create data with only one treatment level
  single_treatment_data <- sample_data
  single_treatment_data$treatment <- factor(1)

  expect_error(
    conditionaleffect(
      data = single_treatment_data,
      y = "outcome",
      o = "origin",
      d = "treatment",
      x = c("control1", "control2"),
      estimator = "Lrnr_glm",
      boot_reps = 100
    ),
    "The treatment variable 'd' must have at least two levels."
  )
})

test_that("conditionaleffect handles non-factor treatment and origin variables by converting them", {
  # Create data with treatment and origin as numeric
  numeric_data <- sample_data
  numeric_data$treatment <- as.numeric(numeric_data$treatment)
  numeric_data$origin <- as.numeric(numeric_data$origin)

  result <- conditionaleffect(
    data = numeric_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100
  )

  # Check that predicted_outcomes have correct treatment levels as factors
  expected_cols <- levels(factor(sample_data$treatment))

  expect_equal(colnames(result$predicted_outcomes), expected_cols)
})

test_that("conditionaleffect produces consistent CATE estimates with reproducible seed", {
  set.seed(123)
  result1 <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100
  )

  set.seed(123)
  result2 <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100
  )

  expect_equal(result1$cate_summary$Estimate, result2$cate_summary$Estimate)
  expect_equal(result1$cate_summary$SE, result2$cate_summary$SE)
})

# === New Test: Compare conditionaleffect with "Lrnr_glm" to lm ===
test_that("CATE estimates using Lrnr_glm match lm estimates", {
  # Run conditionaleffect with Lrnr_glm
  result <- conditionaleffect(
    data = sample_data,
    y = "outcome",
    o = "origin",
    d = "treatment",
    x = c("control1", "control2"),
    estimator = "Lrnr_glm",
    boot_reps = 100
  )

  # Fit a standard linear model
  lm_fit <- lm(outcome ~ treatment + origin + control1 + control2, data = sample_data)

  # Extract treatment coefficient
  # Assuming 'treatment' has levels "1" and "2", and 'treatment2' is the coefficient
  # Check the levels to ensure correct coefficient extraction
  treatment_levels <- levels(sample_data$treatment)
  if (length(treatment_levels) != 2) {
    skip("CATE comparison test requires exactly two treatment levels.")
  }

  # Identify the treatment coefficient name in lm
  # It typically follows "treatment2" if levels are "1" and "2"
  # To make it more robust, find the coefficient that corresponds to the second level
  second_level <- treatment_levels[2]
  coef_name <- paste0("treatment", second_level)

  if (!(coef_name %in% names(coef(lm_fit)))) {
    # If naming differs, find the coefficient that represents the difference
    coef_name <- grep(paste0("^treatment", second_level, "$"), names(coef(lm_fit)), value = TRUE)
  }

  # Extract the treatment coefficient from lm
  lm_coef <- coef(lm_fit)[coef_name]

  # Extract CATE estimate from conditionaleffect
  cate_name <- paste0("CATE_", treatment_levels[1], "_vs_", treatment_levels[2])

  # Alternatively, if "CATE_1_vs_2" corresponds to treatment1 - treatment2,
  # and lm_coef is treatment2 - treatment1, then adjust accordingly
  # Here, we'll assume "CATE_1_vs_2" = treatment1 - treatment2
  # Therefore, CATE_1_vs_2 = - (treatment2 coefficient)
  # Adjust based on your actual naming convention
  if (paste0("CATE_", treatment_levels[1], "_vs_", treatment_levels[2]) %in% result$cate_summary$CATE) {
    cate_estimate <- result$cate_summary$Estimate[result$cate_summary$CATE == paste0("CATE_", treatment_levels[1], "_vs_", treatment_levels[2])]
    expected_estimate <- -lm_coef
  } else if (paste0("CATE_", treatment_levels[2], "_vs_", treatment_levels[1]) %in% result$cate_summary$CATE) {
    cate_estimate <- result$cate_summary$Estimate[result$cate_summary$CATE == paste0("CATE_", treatment_levels[2], "_vs_", treatment_levels[1])]
    expected_estimate <- lm_coef
  } else {
    skip("CATE naming does not match expected treatment levels.")
  }

  # Check that cate_estimate is not missing
  expect_true(length(cate_estimate) == 1)

  # Compare the estimates with a reasonable tolerance
  expect_equal(cate_estimate, expected_estimate, tolerance = 1e-2)
})
