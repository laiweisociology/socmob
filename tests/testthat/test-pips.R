# File: test-pips.R

library(testthat)

# Source necessary functions
# source("R/calculate_odds_ratios.R")
# source("R/contrasttable.R")
# source("R/stepstable.R")
# source("R/pips.R")


# Sample Test Data (3x3 Table)
set.seed(123)  # For reproducibility
toy_data <- data.frame(
  padeg = factor(sample(1:3, 300, replace = TRUE)),
  degree = factor(sample(1:3, 300, replace = TRUE)),
  controls = rnorm(300)
)

# Contrast Table (3x3)
contrast_tbl <- matrix(c(0.2, 0.3, 0.5,
                         0.1, 0.4, 0.5,
                         0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)
rownames(contrast_tbl) <- levels(toy_data$padeg)
colnames(contrast_tbl) <- levels(toy_data$degree)

# Custom Table (3x3)
custom_tbl <- matrix(c(0.25, 0.25, 0.5,
                       0.2, 0.3, 0.5,
                       0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)
rownames(custom_tbl) <- levels(toy_data$padeg)
colnames(custom_tbl) <- levels(toy_data$degree)

# Test Suite for pips Function
test_that("pips returns a dataframe with correct dimensions and probabilities", {
  pips_df <- pips(data = toy_data,   o = "padeg", d = "degree",
                  x = c("controls"), contrast = contrast_tbl)

  expect_s3_class(pips_df, "data.frame")
  expect_equal(nrow(pips_df), nrow(toy_data))
  expect_equal(ncol(pips_df), ncol(contrast_tbl))
  expect_true(all(abs(rowSums(pips_df) - 1) < 1e-6))
})


test_that("pips throws an error when multiple intervention parameters are specified", {
  expect_error(
    pips(data = toy_data,   o = "padeg", d = "degree",
         x = c("controls"), custom.table = custom_tbl, contrast = contrast_tbl),
    "Only one of `custom.table`, `strength`, or `contrast` must be specified."
  )

  expect_error(
    pips(data = toy_data,   o = "padeg", d = "degree",
         x = c("controls"), strength = 0.5, contrast = contrast_tbl),
    "Only one of `custom.table`, `strength`, or `contrast` must be specified."
  )
})

test_that("pips defaults to strength=1 when no intervention parameter is specified", {
  pips_df <- pips(data = toy_data,   o = "padeg", d = "degree",
                  x = c("controls"))

  expect_s3_class(pips_df, "data.frame")
  expect_equal(nrow(pips_df), nrow(toy_data))
  expect_equal(ncol(pips_df), ncol(table(toy_data$padeg, toy_data$degree)))
  expect_true(all(abs(rowSums(pips_df) - 1) < 1e-6))
})


test_that("pips works with custom.table", {
  pips_df <- pips(data = toy_data,   o = "padeg", d = "degree",
                  x = c("controls"), custom.table = custom_tbl)

  expect_s3_class(pips_df, "data.frame")
  expect_equal(nrow(pips_df), nrow(toy_data))
  expect_equal(ncol(pips_df), ncol(custom_tbl))
  expect_true(all(abs(rowSums(pips_df) - 1) < 1e-6))
})

test_that("pips works with strength parameter", {
  # Assuming strength parameter modifies towards uniform distribution
  pips_df <- pips(data = toy_data,   o = "padeg", d = "degree",
                  x = c("controls"), strength = 0.5)

  expect_s3_class(pips_df, "data.frame")
  expect_equal(nrow(pips_df), nrow(toy_data))
  expect_equal(ncol(pips_df), ncol(table(toy_data$padeg, toy_data$degree)))
  expect_true(all(abs(rowSums(pips_df) - 1) < 1e-6))
})

test_that("pips throws an error when y, d, o, or x contain missing values", {
  # Introduce missing values
  toy_data_na <- toy_data
  toy_data_na$degree[1] <- NA
  toy_data_na$controls[2] <- NA

  expect_error(
    pips(data = toy_data_na,   o = "padeg", d = "degree",
         x = c("controls"), contrast = contrast_tbl),
    "The following variables contain missing values: degree, controls"
  )

  toy_data_na2 <- toy_data
  toy_data_na2$padeg[3] <- NA

  expect_error(
    pips(data = toy_data_na2,   o = "padeg", d = "degree",
         x = c("controls"), custom.table = custom_tbl),
    "The following variables contain missing values: padeg"
  )

  toy_data_na3 <- toy_data
  toy_data_na3$controls[4] <- NA

  expect_error(
    pips(data = toy_data_na3,   o = "padeg", d = "degree",
         x = c("controls"), strength = 0.3),
    "The following variables contain missing values: controls"
  )
})
