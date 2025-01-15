# File: R/socmob.R

#' Estimating social mobility effect (full bootstrap over conditionaleffect & pips)
#'
#' @description
#' This function:
#'   1. Calls \code{\link{conditionaleffect}} on the full data (with \code{boot_reps = 1}) to get a baseline predicted outcome structure.
#'   2. Uses \code{pips(..., strength=0)} to compute the model-based "observed" distribution (and outcome).
#'   3. Uses user-supplied or newly computed \code{pips(...)} for the post-intervention distribution.
#'   4. Computes the difference in outcomes (post - observed) across individuals.
#'   5. Runs a bootstrap of size \code{boot_reps} that re-fits the outcome model (conditionaleffect) and re-computes pips for each replicate:
#'      - i.e., for each replicate we do:
#'         (a) sample data with replacement,
#'         (b) conditionaleffect(boot_data, ..., boot_reps=1),
#'         (c) pips(boot_data),
#'         (d) compute observed & post-outcomes,
#'         (e) compute overall/gainer/giver + by-origin/destination.
#'
#' The final standard errors for all mobility effects are derived from that external bootstrap.
#'
#' @param data A dataframe
#' @param y The outcome variable (string)
#' @param d The destination/treatment variable (string)
#' @param o The origin variable (string)
#' @param x A character vector of additional control variables (confounders)
#' @param estimator A learner name or vector of learner names from sl3 (default "Lrnr_glm")
#' @param custom.pips Optionally, a dataframe of post-intervention propensities (columns = treatment levels)
#' @param custom.table A matrix/table for the post-intervention scenario (rows = origin, columns = destination)
#' @param strength A numeric \eqn{\in [0,1]} specifying the log-odds ratio change in \code{\link{stepstable}}
#' @param contrast A table specifying odds ratios to form a new interventional scenario via \code{\link{contrasttable}}
#' @param boot_reps The number of external bootstrap repetitions for final mobility effect SE. Default 200.
#'
#' @return A named list with:
#' \item{predicted_outcomes}{Dataframe of predicted outcomes (from full data, boot_reps=1 in conditionaleffect)}
#' \item{pips_observed}{Observed (pre-intervention) propensity from pips(..., strength=0)}
#' \item{pips_post}{Post-intervention propensity (computed or user-supplied)}
#' \item{observed_outcome}{Model-based observed outcome = sum of predicted outcomes * pips_observed}
#' \item{post_outcome}{Post-intervention outcome = sum of predicted outcomes * pips_post}
#' \item{mobility_effect}{Overall effect (overall, gainer, giver) + boot SE}
#' \item{origin_specific}{By-origin difference (pre vs. post) + boot SE}
#' \item{destination_specific}{By-destination difference (pre vs. post) + boot SE}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- socmob(
#'   data     = toy_data,
#'   y        = "child_outcome",
#'   o        = "origin_var",
#'   d        = "treatment_var",
#'   x        = c("confound1","confound2"),
#'   strength = 1,
#'   boot_reps= 200
#' )
#' }
socmob <- function(
    data, y, o, d, x,
    estimator    = "Lrnr_glm",
    custom.pips  = NULL,
    custom.table = NULL,
    strength     = NULL,
    contrast     = NULL,
    boot_reps    = 200
) {

  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required but not installed.")
  }

  # 1) Fit conditionaleffect once on the full data (boot_reps=1 => no internal bootstrap).
  cond_full <- conditionaleffect(
    data      = data,
    y         = y,
    o         = o,
    d         = d,
    x         = x,
    estimator = estimator,
    boot_reps = 1
  )
  predicted_outcomes <- cond_full$predicted_outcomes

  # 2) pips(..., strength=0) => model-based "observed" distribution.
  #    Then sum across predicted_outcomes to get "observed_outcome".
  pips_observed <- pips(data, o = o, d = d, x = x, strength=0)
  treat_levels <- colnames(predicted_outcomes)
  n_obs <- nrow(data)
  observed_outcome <- numeric(n_obs)
  for (i in seq_len(n_obs)) {
    row_val <- 0
    for (lvl in seq_along(treat_levels)) {
      row_val <- row_val + predicted_outcomes[i, lvl] * pips_observed[i, lvl]
    }
    observed_outcome[i] <- row_val
  }

  # 3) Either accept custom pips or compute new pips(...).
  if (!is.null(custom.pips)) {
    pips_post <- custom.pips
  } else {
    pips_post <- pips(
      data         = data,
      o            = o,
      d            = d,
      x            = x,
      custom.table = custom.table,
      strength     = strength,
      contrast     = contrast
    )
  }
  # Post-outcome
  post_outcome <- numeric(n_obs)
  for (i in seq_len(n_obs)) {
    row_val <- 0
    for (lvl in seq_along(treat_levels)) {
      row_val <- row_val + predicted_outcomes[i, lvl] * pips_post[i, lvl]
    }
    post_outcome[i] <- row_val
  }

  # Observed vs. post difference
  diffs <- post_outcome - observed_outcome
  overall_est <- mean(diffs)
  gainer_est  <- mean(diffs[diffs > 0])
  giver_est   <- mean(diffs[diffs < 0])

  origin_var    <- data[[o]]
  origin_levels <- levels(origin_var)
  origin_effect <- data.frame(
    origin       = origin_levels,
    pre_outcome  = numeric(length(origin_levels)),
    post_outcome = numeric(length(origin_levels)),
    effect       = numeric(length(origin_levels)),
    stringsAsFactors = FALSE
  )
  for (jj in seq_along(origin_levels)) {
    idx <- which(origin_var == origin_levels[jj])
    origin_effect$pre_outcome[jj]  <- mean(observed_outcome[idx])
    origin_effect$post_outcome[jj] <- mean(post_outcome[idx])
    origin_effect$effect[jj]       <- origin_effect$post_outcome[jj] - origin_effect$pre_outcome[jj]
  }

  destination_var     <- data[[d]]
  destination_levels  <- levels(destination_var)
  destination_effect  <- data.frame(
    destination   = destination_levels,
    pre_outcome   = numeric(length(destination_levels)),
    post_outcome  = numeric(length(destination_levels)),
    effect        = numeric(length(destination_levels)),
    stringsAsFactors=FALSE
  )
  for (jj in seq_along(destination_levels)) {
    idx <- which(destination_var == destination_levels[jj])
    destination_effect$pre_outcome[jj]  <- mean(observed_outcome[idx])
    destination_effect$post_outcome[jj] <- mean(post_outcome[idx])
    destination_effect$effect[jj]       <- destination_effect$post_outcome[jj] - destination_effect$pre_outcome[jj]
  }

  # 4) A single pass "point estimate" is done. Next we do the external bootstrap (size=boot_reps).
  #    Each replicate:
  #      - sample with replacement => boot_data
  #      - re-run conditionaleffect(boot_data, boot_reps=1)
  #      - re-run pips(boot_data, strength=0) => pips_obs_boot
  #      - re-run pips(boot_data, strength|table|contrast) => pips_post_boot
  #      - compute observed vs. post difference => overall, gainer, giver, origin, destination

  boot_fun <- function(index) {
    # index is a vector of row indices for the bootstrap sample
    boot_data <- data[index, , drop=FALSE]

    # (a) conditionaleffect with no internal boot
    cond_boot <- conditionaleffect(
      data      = boot_data,
      y         = y,
      o         = o,
      d         = d,
      x         = x,
      estimator = estimator,
      boot_reps = 1
    )
    pred_boot <- cond_boot$predicted_outcomes

    # (b) pips(..., strength=0) => observed distribution
    pips_obs_boot <- pips(boot_data, o=o, d=d, x=x, strength=0)

    # observed outcome
    n_boot <- nrow(boot_data)
    obs_out_boot <- numeric(n_boot)
    treat_lvls   <- colnames(pred_boot)
    for (i in seq_len(n_boot)) {
      row_val <- 0
      for (lvl in seq_along(treat_lvls)) {
        row_val <- row_val + pred_boot[i, lvl] * pips_obs_boot[i, lvl]
      }
      obs_out_boot[i] <- row_val
    }

    # (c) pips post
    if (!is.null(custom.pips)) {
      # If custom pips was provided, we can't re-fit it for the bootstrap sample,
      # but we might re-sample the same custom pips row indices => not typical though
      # We'll assume custom pips is not re-fit
      pips_post_boot <- custom.pips[index, , drop=FALSE]
    } else {
      pips_post_boot <- pips(
        data         = boot_data,
        o            = o,
        d            = d,
        x            = x,
        custom.table = custom.table,
        strength     = strength,
        contrast     = contrast
      )
    }

    # post outcome
    post_out_boot <- numeric(n_boot)
    for (i in seq_len(n_boot)) {
      row_val <- 0
      for (lvl in seq_along(treat_lvls)) {
        row_val <- row_val + pred_boot[i, lvl] * pips_post_boot[i, lvl]
      }
      post_out_boot[i] <- row_val
    }

    # diffs
    diffs_bs   <- post_out_boot - obs_out_boot
    overall_bs <- mean(diffs_bs)
    gainer_bs  <- mean(diffs_bs[diffs_bs > 0])
    giver_bs   <- mean(diffs_bs[diffs_bs < 0])

    # origin
    origin_bs <- factor(boot_data[[o]], levels=levels(data[[o]]))
    n_orig <- length(levels(data[[o]]))
    preO  <- numeric(n_orig)
    postO <- numeric(n_orig)
    effO  <- numeric(n_orig)
    for (kk in seq_len(n_orig)) {
      lvl_name <- levels(data[[o]])[kk]
      idx2 <- which(origin_bs == lvl_name)
      if (length(idx2) == 0) {
        preO[kk]  <- NA
        postO[kk] <- NA
        effO[kk]  <- NA
      } else {
        preO[kk]  <- mean(obs_out_boot[idx2])
        postO[kk] <- mean(post_out_boot[idx2])
        effO[kk]  <- postO[kk] - preO[kk]
      }
    }

    # destination
    dest_bs <- factor(boot_data[[d]], levels=levels(data[[d]]))
    n_dest <- length(levels(data[[d]]))
    preD  <- numeric(n_dest)
    postD <- numeric(n_dest)
    effD  <- numeric(n_dest)
    for (kk in seq_len(n_dest)) {
      lvl_name <- levels(data[[d]])[kk]
      idx2 <- which(dest_bs == lvl_name)
      if (length(idx2) == 0) {
        preD[kk]  <- NA
        postD[kk] <- NA
        effD[kk]  <- NA
      } else {
        preD[kk]  <- mean(obs_out_boot[idx2])
        postD[kk] <- mean(post_out_boot[idx2])
        effD[kk]  <- postD[kk] - preD[kk]
      }
    }

    # Combine into vector:
    #  3 => overall, gainer, giver
    #  3*n_orig => origin pre, post, effect
    #  3*n_dest => destination pre, post, effect
    c(
      overall_bs,
      gainer_bs,
      giver_bs,
      preO,
      postO,
      effO,
      preD,
      postD,
      effD
    )
  }

  # We now run the external bootstrap
  set.seed(123)
  n <- nrow(data)
  boot_out <- boot::boot(
    data      = seq_len(n),  # pass row indices
    statistic = function(xx, idx) boot_fun(idx), # 'xx' is ignored
    R         = boot_reps
  )

  # Summarize SE
  # The dimension of each replicate vector:
  #   3 + 3*n_orig + 3*n_dest
  origin_levels <- levels(data[[o]])
  destination_levels <- levels(data[[d]])
  n_orig <- length(origin_levels)
  n_dest <- length(destination_levels)

  block_len <- 3 + 3*n_orig + 3*n_dest
  if (ncol(boot_out$t) != block_len) {
    warning("Mismatch in length of the bootstrap replicate vector. Check code.")
  }
  se_vec <- apply(boot_out$t, 2, sd, na.rm=TRUE)

  # parse
  overall_se <- se_vec[1]
  gainer_se  <- se_vec[2]
  giver_se   <- se_vec[3]

  # origin => 3 blocks => pre, post, eff
  idx4 <- 4
  idxO_pre  <- idx4:(3 + n_orig)
  idxO_post <- (3 + n_orig + 1):(3 + 2*n_orig)
  idxO_eff  <- (3 + 2*n_orig + 1):(3 + 3*n_orig)

  # dest => next 3 blocks => pre, post, eff
  idxD_pre  <- (3 + 3*n_orig + 1):(3 + 3*n_orig + n_dest)
  idxD_post <- (3 + 3*n_orig + n_dest + 1):(3 + 3*n_orig + 2*n_dest)
  idxD_eff  <- (3 + 3*n_orig + 2*n_dest + 1):(3 + 3*n_orig + 3*n_dest)

  mobility_df <- data.frame(
    measure  = c("overall", "gainer", "giver"),
    estimate = c(overall_est, gainer_est, giver_est),
    se       = c(overall_se, gainer_se, giver_se),
    stringsAsFactors=FALSE
  )

  origin_effect$pre_outcome_se  <- se_vec[idxO_pre]
  origin_effect$post_outcome_se <- se_vec[idxO_post]
  origin_effect$effect_se       <- se_vec[idxO_eff]

  destination_effect$pre_outcome_se  <- se_vec[idxD_pre]
  destination_effect$post_outcome_se <- se_vec[idxD_post]
  destination_effect$effect_se       <- se_vec[idxD_eff]

  # Return final
  list(
    predicted_outcomes   = predicted_outcomes,
    pips_observed        = pips_observed,
    pips_post            = pips_post,
    observed_outcome     = observed_outcome,
    post_outcome         = post_outcome,
    mobility_effect      = mobility_df,
    origin_specific      = origin_effect,
    destination_specific = destination_effect
  )
}
