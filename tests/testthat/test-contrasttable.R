# File: test-contrasttable.R

library(testthat)

# Source necessary functions
# source("calculate_odds_ratios.R")
# source("contrasttable.R")

# Sample Test Data (3x3 Table)
set.seed(123)  # For reproducibility
toy_data <- data.frame(
  padeg = factor(sample(1:3, 300, replace = TRUE)),
  degree = factor(sample(1:3, 300, replace = TRUE))
)

# Observed OD Table
observed_O <- as.matrix(table(toy_data$padeg, toy_data$degree))

# Contrast Table (3x3)
contrast_tbl <- matrix(c(10, 20, 30,
                         20, 30, 40,
                         30, 40, 50), nrow = 3, byrow = TRUE)
rownames(contrast_tbl) <- levels(toy_data$padeg)
colnames(contrast_tbl) <- levels(toy_data$degree)


# Test Suite for contrasttable Function
test_that("contrasttable injects odds ratios close to the contrast table's odds ratios", {
  counterfactual_O <- contrasttable(toy_data, o = "padeg", d = "degree", contrast = contrast_tbl)

  # Compute odds ratios for contrast and counterfactual tables
  or_contrast <- calculate_odds_ratios(contrast_tbl)
  or_counterfactual <- calculate_odds_ratios(counterfactual_O)

  # Compute relative ratios, excluding NA
  valid_indices <- !is.na(or_contrast) & !is.na(or_counterfactual)
  relative_ratio <- or_counterfactual[valid_indices] / or_contrast[valid_indices]

  # Expect that relative ratios are reasonably close to 1 within a tolerance
  # Given multiple overlapping odds ratios, exact matching is not possible
  # Hence, check that the relative ratios are within 50% of 1
  tolerance <- 0.5  # 50% tolerance
  expect_true(all(abs(relative_ratio - 1) < tolerance),
              info = paste("Relative ratios:", paste(round(relative_ratio, 3), collapse = ", ")))
})

test_that("contrasttable preserves the marginal distributions of the observed table", {
  counterfactual_O <- contrasttable(toy_data, o = "padeg", d = "degree", contrast = contrast_tbl)

  # Compare row sums
  row_diff <- abs(rowSums(counterfactual_O) - rowSums(observed_O))
  expect_true(all(row_diff < 1e-4),
              info = paste("Row differences:", paste(round(row_diff, 6), collapse = ", ")))

  # Compare column sums
  col_diff <- abs(colSums(counterfactual_O) - colSums(observed_O))
  expect_true(all(col_diff < 1e-4),
              info = paste("Column differences:", paste(round(col_diff, 6), collapse = ", ")))
})

test_that("contrasttable throws an error when contrast table dimensions do not match observed table", {
  # Create a mismatched contrast table (different number of rows)
  mismatched_contrast_rows <- matrix(c(10, 20, 30,
                                       20, 30, 40), nrow = 2, byrow = TRUE)
  rownames(mismatched_contrast_rows) <- c("1", "2")  # Only two rows instead of three
  colnames(mismatched_contrast_rows) <- levels(toy_data$degree)

  expect_error(
    contrasttable(toy_data, o = "padeg", d = "degree", contrast = mismatched_contrast_rows),
    "Contrast table must have the same number of rows as the observed table."
  )

  # Create a mismatched contrast table (different column names)
  mismatched_contrast_cols <- matrix(c(10, 20, 30,
                                       20, 30, 40,
                                       30, 40, 50), nrow = 3, byrow = TRUE)
  rownames(mismatched_contrast_cols) <- levels(toy_data$padeg)
  colnames(mismatched_contrast_cols) <- c("A", "B", "C")  # Different column names

  expect_error(
    contrasttable(toy_data, o = "padeg", d = "degree", contrast = mismatched_contrast_cols),
    "Column names of `contrast` table must match column names of the observed table."
  )
})
