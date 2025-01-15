# File: tests/testthat/test-socmob.R

library(testthat)
library(boot)

# source("R/conditionaleffect.R")
# source("R/pips.R")
# source("R/stepstable.R")
# source("R/contrasttable.R")
# source("R/calculate_odds_ratios.R")
# source("R/socmob.R")


set.seed(101)
toy_data <- data.frame(
  padeg   = factor(sample(1:3, 50, replace = TRUE)),
  degree  = factor(sample(1:3, 50, replace = TRUE)),
  sibs    = rnorm(50),
  ethnic  = sample(c("group1", "group2"), 50, replace = TRUE),
  outcome = rnorm(50)
)

test_that("socmob returns correct structure with default boot_reps=200", {
  result <- socmob(
    data       = toy_data,
    y          = "outcome",
    o          = "padeg",
    d          = "degree",
    x          = c("sibs","ethnic"),
    strength   = 0.5,
    boot_reps  = 50  # for a quicker test
  )

  expect_type(result, "list")
  for (nm in c("predicted_outcomes","pips_observed","pips_post",
               "observed_outcome","post_outcome","mobility_effect",
               "origin_specific","destination_specific")) {
    expect_true(nm %in% names(result))
  }

  # Check mobility_effect => columns measure, estimate, se
  expect_true(all(c("measure","estimate","se") %in% names(result$mobility_effect)))

  # Check origin_specific => columns pre_outcome, post_outcome, effect + their se
  for (cc in c("pre_outcome","post_outcome","effect",
               "pre_outcome_se","post_outcome_se","effect_se")) {
    expect_true(cc %in% names(result$origin_specific))
  }
  # Check dimension matches factor levels
  expect_equal(nrow(result$origin_specific), nlevels(toy_data$padeg))

  # Similar check for destination_specific
  for (cc in c("pre_outcome","post_outcome","effect",
               "pre_outcome_se","post_outcome_se","effect_se")) {
    expect_true(cc %in% names(result$destination_specific))
  }
  expect_equal(nrow(result$destination_specific), nlevels(toy_data$degree))
})

test_that("socmob with strength=0 yields near-zero difference between observed and post-outcome", {
  # With strength=0 => post distribution ~ observed distribution => difference ~ 0 on average
  result0 <- socmob(
    data      = toy_data,
    y         = "outcome",
    o         = "padeg",
    d         = "degree",
    x         = c("sibs","ethnic"),
    strength  = 0,
    boot_reps = 30
  )

  diffs <- result0$post_outcome - result0$observed_outcome
  expect_lt(abs(mean(diffs)), 1e-6) # or a small threshold
})

test_that("socmob handles custom pips properly", {
  # Create a custom pips with exactly 3 columns, each row sums to 1
  # e.g., everyone is put in the first category => post_outcome ~ predicted_outcomes col 1
  custom_pips_df <- data.frame(
    pips.1 = rep(1, nrow(toy_data)),
    pips.2 = rep(0, nrow(toy_data)),
    pips.3 = rep(0, nrow(toy_data))
  )

  result_custom <- socmob(
    data       = toy_data,
    y          = "outcome",
    o          = "padeg",
    d          = "degree",
    x          = c("sibs"),
    custom.pips= custom_pips_df,
    boot_reps  = 20
  )

  # Now post_outcome should match predicted_outcomes col 1
  pcol1 <- result_custom$predicted_outcomes[,1]
  expect_equal(result_custom$post_outcome, pcol1, tolerance=1e-6)
})

test_that("socmob runs with custom.table and no 'strength' or 'contrast'", {
  # minimal custom.table
  tab_custom <- matrix(
    c(0.3, 0.3, 0.4,
      0.2, 0.5, 0.3,
      0.1, 0.6, 0.3),
    nrow=3, byrow=TRUE,
    dimnames = list(c("1","2","3"), c("1","2","3"))
  )
  # As the data has factor levels '1','2','3' for both origin & destination
  result_tab <- socmob(
    data         = toy_data,
    y            = "outcome",
    o            = "padeg",
    d            = "degree",
    x            = c("sibs"),
    custom.table = tab_custom,
    boot_reps    = 15
  )
  # Just check we get a result with no error
  expect_type(result_tab, "list")
})

