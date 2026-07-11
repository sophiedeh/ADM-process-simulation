# =============================================================================
# Helpers.R
# Shared paths, loaders and utilities used by every analysis script.
# =============================================================================

library(dplyr)
library(ggplot2)
library(here)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
SIMULATION_DIR <- here("Simulation", "simulated")
RESULTS_DIR  <- here("Results")

# -----------------------------------------------------------------------------
# Load / save simulation output
# -----------------------------------------------------------------------------
# Load a previously saved simulation for one training programme.

load_sim <- function(training_name) {
  path <- file.path(SIMULATION_DIR, paste0(training_name, ".RData"))
  if (!file.exists(path)) {
    stop("No simulation file found at: ", path,
         "\nRun Simulation.R first.")
  }
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env$sim
}

# Save a simulation list to Simulation/simulated/<training_name>.RData
save_sim <- function(sim) {
  dir.create(SIMULATION_DIR, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(SIMULATION_DIR, paste0(sim$name, ".RData"))
  save(sim, file = path)
  message("Saved simulation to: ", path)
}

# -----------------------------------------------------------------------------
# Results helpers
# -----------------------------------------------------------------------------
# Save a ggplot to Results/<training_name>/<filename>

save_plot <- function(plot, filename, training_name, decision_name = NULL, ...) {
  dir <- if (!is.null(decision_name)) {
    file.path(RESULTS_DIR, training_name, decision_name)
  } else {
    file.path(RESULTS_DIR, training_name)
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  full_path <- file.path(dir, filename)
  ggsave(full_path, plot, ...)
  message("Saved plot to: ", full_path)
}

# Save a table
save_table <- function(df, filename, training_name, decision_name = NULL) {
  dir <- if (!is.null(decision_name)) {
    file.path(RESULTS_DIR, training_name, decision_name)
  } else {
    file.path(RESULTS_DIR, training_name)
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  full_path <- file.path(dir, filename)
  write.csv(df, full_path, row.names = FALSE)
  message("Saved table to: ", full_path)
}

# -----------------------------------------------------------------------------
# Shared binning helper (used by binned_random and distribution_comparison)
# -----------------------------------------------------------------------------
# Add a risk_bin column to a data frame using consistent breaks.

add_risk_bins <- function(df, bin_width = 0.1, breaks = NULL) {
  if (is.null(breaks)) {
    breaks <- seq(
      floor(min(df$risk_score_log,   na.rm = TRUE) / bin_width) * bin_width,
      ceiling(max(df$risk_score_log, na.rm = TRUE) / bin_width) * bin_width,
      by = bin_width
    )
  }
  df$risk_bin <- cut(df$risk_score_log, breaks = breaks,
                     include.lowest = TRUE, right = FALSE)
  df
}

# Build a short run tag (decision + capacity + randomization + threshold)
# for use inside saved filenames.
run_tag <- function(decision_name, capacity, pct_random, threshold) {
  sprintf(
    "%s_cap%s_ran%s_thres%s",
    decision_name,
    format(capacity,   nsmall = 2),
    format(pct_random, nsmall = 2),
    format(threshold,  nsmall = 2)
  )
}

# -----------------------------------------------------------------------------
# Significance label helper (used by binned_random and IV scripts)
# -----------------------------------------------------------------------------

sig_label <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# -----------------------------------------------------------------------------
# Representative control group construction
# -----------------------------------------------------------------------------
# Implementation of method to construct a valid representative control group.
build_representative_control <- function(sim, almp, seed = NULL) {
  threshold <- sim$params$score_threshold
  
  # ── Build the two halves of the control group ────────────────────────────────
  # not_elig : below-threshold units that were never reviewed
  # not_above: above-threshold units that lost the algorithmic lottery
  not_elig  <- sim$not_selected[sim$not_selected$risk_score_log < threshold, ]
  not_above <- sim$not_selected_above_threshold
  
  # Align columns: simulation_core adds human_decision/source only to not_above
  not_elig$human_decision <- 0L
  not_elig$source         <- NA_character_
  
  current_control <- rbind(not_elig, not_above)
  
  # ── K-S test: is the current control representative of the full population? ──
  # Tests whether the two risk-score distributions are drawn from the same
  # underlying distribution. Rejection (p < 0.05) triggers the adjustment below.
  ks_orig <- suppressWarnings(
    ks.test(current_control$risk_score_log, almp$risk_score_log)
  )
  
  ks_fails <- ks_orig$p.value < 0.05
  
  if (!ks_fails) {
    return(list(
      control        = current_control,
      adjusted       = FALSE,
      ks_p_original  = ks_orig$p.value,
      ks_p_adjusted  = NA_real_,
      pct_algo       = NA_real_,
      n_removed      = NA_integer_
    ))
  }
  
  cat(sprintf("  Adjustment triggered by: K-S p = %.4f < 0.05\n", ks_orig$p.value))
  
  # ── Adjustment ───────────────────────────────────────────────────────────────
  n_eligible_remaining <- nrow(
    sim$remaining[sim$remaining$risk_score_log >= threshold, ]
  )
  pct_algo <- sim$params$n_algo / n_eligible_remaining
  n_remove <- round(pct_algo * nrow(not_elig))
  
  if (!is.null(seed)) set.seed(seed)
  idx_remove       <- sample(nrow(not_elig), size = n_remove, replace = FALSE)
  trimmed_elig     <- not_elig[-idx_remove, ]
  adjusted_control <- rbind(trimmed_elig, not_above)
  
  # ── K-S test on the adjusted group ───────────────────────────────────────────
  ks_adj <- suppressWarnings(
    ks.test(adjusted_control$risk_score_log, almp$risk_score_log)
  )
  
  ks_adj_passes <- ks_adj$p.value >= 0.05
  
  cat(sprintf(
    "  Control group adjusted: removed %d ineligible units (%.1f%% of ineligible pool).\n",
    n_remove, 100 * n_remove / nrow(not_elig)
  ))
  cat(sprintf(
    "  K-S:  before p = %.4f (FAIL) | after p = %.4f (%s)\n",
    ks_orig$p.value,
    ks_adj$p.value,
    ifelse(ks_adj_passes, "pass", "FAIL")
  ))
  if (ks_adj_passes) {
    cat("  => Adjusted control group passes K-S test: representative.\n")
  } else {
    warning(
      "Adjusted control group does not pass the K-S representativeness check after adjustment.\n",
      "  Consider wider binning or a different adjustment strategy.",
      call. = FALSE
    )
  }
  
  list(
    control        = adjusted_control,
    adjusted       = TRUE,
    ks_p_original  = ks_orig$p.value,
    ks_p_adjusted  = ks_adj$p.value,
    pct_algo       = pct_algo,
    n_removed      = n_remove
  )
}