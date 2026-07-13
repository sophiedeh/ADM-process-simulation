# =============================================================================
# scenario_1.R  —  RDD evaluation (vocational_training)
# =============================================================================
# Estimands (with robust CIs and null-hypothesis tests for each):
#   - Sharp RDD  -> ITT  (effect of becoming eligible)
#   - Fuzzy RDD  -> LATE (effect of treatment for compliers)
# =============================================================================

library(here)
library(rdrobust)
library(ggplot2)
library(dplyr)

source(here("Simulation", "Helpers.R"))

TRAINING <- "vocational_training"

# -----------------------------------------------------------------------------
# Function that prints a null-hypothesis test
# -----------------------------------------------------------------------------
print_rd_test <- function(obj, label, h0 = "effect = 0") {
  cat(sprintf("\n=== %s — NULL HYPOTHESIS TEST (H0: %s) ===\n", label, h0))
  if (is.null(obj)) {
    cat("Not available.\n")
    return(invisible(NULL))
  }
  cat("Estimate (conventional):  ", round(obj$coef[1], 6), "\n")
  cat("Robust SE:                ", round(obj$se[3],   6), "\n")
  cat("z-statistic (robust):     ", round(obj$z[3],    4), "\n")
  cat("p-value (conventional):   ", round(obj$pv[1],   6), "\n")
  cat("p-value (robust, 2-sided):", round(obj$pv[3],   6), "\n")
  cat("95% CI (robust):          [",
      round(obj$ci[3, 1], 6), ",", round(obj$ci[3, 2], 6), "]\n")
  if (obj$pv[3] < 0.05) {
    cat("=> Reject H0 at 5%: effect significantly different from zero.\n")
  } else {
    cat("=> Fail to reject H0 at 5%: effect not significantly different from zero.\n")
  }
}

# =============================================================================
# Load simulation output and build the analysis sample according to scenario 4
# =============================================================================
sim             <- load_sim(TRAINING)
decision_name   <- sim$params$decision_fn_name
threshold       <- sim$params$score_threshold

# Tags to create output dirs per variable variation
cap_tag <- sprintf("cap_%s", format(sim$params$capacity_multiplier, nsmall = 2))
ran_tag <- sprintf("ran_%s", format(sim$params$pct_random, nsmall = 2))
threshold_tag <- sprintf("thres_%s", format(sim$params$score_threshold, nsmall = 2))
out_dir <- file.path(RESULTS_DIR, "scenario_1", decision_name, cap_tag, ran_tag, threshold_tag, TRAINING)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TAG <- run_tag(decision_name, sim$params$capacity_multiplier,
               sim$params$pct_random, sim$params$score_threshold)

vt_score        <- sim$selected[sim$selected$source == "score", ]
vt_not_selected <- sim$not_selected[sim$not_selected$risk_score_log < threshold, ]
vt_not_selected$source         <- NA
vt_not_selected$human_decision <- NA
vt_combined <- rbind(vt_score, vt_not_selected)

# Assignment: 1 if human approved, 0 otherwise (including unreviewed)
vt_combined$assignment <- ifelse(
  !is.na(vt_combined$human_decision) & vt_combined$human_decision == 1,
  1L, 0L
)

# Outcomes are assigned inside the simulation_score
if (anyNA(vt_combined$outcome)) {
  stop(sprintf("[%s] %d NA outcome(s) in the analysis sample — the simulation should assign every outcome.",
               TRAINING, sum(is.na(vt_combined$outcome))))
}

R <- vt_combined$risk_score_log

cat("\n\n##############################################################\n")
cat("##  RDD EVALUATION:", TRAINING, "\n")
cat("##############################################################\n")

# -----------------------------------------------------------------------------
# RDD plot
# -----------------------------------------------------------------------------
plot_path <- file.path(out_dir, sprintf("rdd_plot_%s.png", TAG))

png(plot_path, width = 800, height = 600)
rdp <- rdplot(y = vt_combined$outcome, x = R, c = threshold, p = 1,
              x.label = "R (-)", y.label = "Y (-)", hide = TRUE)
print(rdp$rdplot)   # rdplot returns a ggplot; inside a function it must be
# printed explicitly (top-level auto-print does not apply)
dev.off()
message("Saved RDD plot to: ", plot_path)

# -----------------------------------------------------------------------------
# Sharp RDD — ITT 
# -----------------------------------------------------------------------------
sharp <- tryCatch(
  rdrobust(y = vt_combined$outcome, x = R, c = threshold, p = 1),
  error = function(e) { message("Sharp RDD failed: ", e$message); NULL }
)

print_rd_test(sharp, "SHARP RDD (ITT)", h0 = "ITT = 0")
if (!is.null(sharp)) {
  cat("Effective N (below/above):", sharp$N_h[1], "/", sharp$N_h[2], "\n")
}

# -----------------------------------------------------------------------------
# First stage
# -----------------------------------------------------------------------------
fs <- tryCatch(
  rdrobust(y = vt_combined$assignment, x = R, c = threshold, p = 1),
  error = function(e) { message("First stage failed: ", e$message); NULL }
)
fs_jump <- if (!is.null(fs)) fs$coef[1] else NA_real_

print_rd_test(fs, "FIRST STAGE", h0 = "first-stage jump = 0")

# -----------------------------------------------------------------------------
# Fuzzy RDD — LATE
# -----------------------------------------------------------------------------
fuzzy_obj <- if (!is.null(fs) && !is.na(fs_jump) && abs(fs_jump) > 0.01) {
  tryCatch(
    rdrobust(y = vt_combined$outcome, x = R, c = threshold, p = 1,
             fuzzy = vt_combined$assignment),
    error = function(e) { message("Fuzzy RDD failed: ", e$message); NULL }
  )
} else {
  message("First stage too weak for fuzzy RDD (jump = ", round(fs_jump, 4), ")")
  NULL
}
fuzzy_est <- if (!is.null(fuzzy_obj)) fuzzy_obj$coef[1] else NA_real_

print_rd_test(fuzzy_obj, "FUZZY RDD (LATE)", h0 = "LATE = 0")

# =============================================================================
# One-row results summary
# =============================================================================
rdd_results <- data.frame(
  training  = TRAINING,
  decision  = decision_name,
  threshold = threshold,

  # --- Sharp RDD (ITT) ---
  sharp_estimate  = if (!is.null(sharp)) sharp$coef[1]   else NA_real_,
  sharp_se_robust = if (!is.null(sharp)) sharp$se[3]     else NA_real_,
  sharp_p_robust  = if (!is.null(sharp)) sharp$pv[3]     else NA_real_,
  sharp_ci_lo     = if (!is.null(sharp)) sharp$ci[3, 1]  else NA_real_,
  sharp_ci_hi     = if (!is.null(sharp)) sharp$ci[3, 2]  else NA_real_,
  n_below         = if (!is.null(sharp)) sharp$N_h[1]    else NA_real_,
  n_above         = if (!is.null(sharp)) sharp$N_h[2]    else NA_real_,
  bandwidth_mse   = if (!is.null(sharp)) sharp$bws[1, 1] else NA_real_,

  # --- First stage ---
  fs_estimate  = fs_jump,
  fs_se_robust = if (!is.null(fs)) fs$se[3]    else NA_real_,
  fs_p_robust  = if (!is.null(fs)) fs$pv[3]    else NA_real_,
  fs_ci_lo     = if (!is.null(fs)) fs$ci[3, 1] else NA_real_,
  fs_ci_hi     = if (!is.null(fs)) fs$ci[3, 2] else NA_real_,

  # --- Fuzzy RDD (LATE) ---
  fuzzy_estimate  = fuzzy_est,
  fuzzy_se_robust = if (!is.null(fuzzy_obj)) fuzzy_obj$se[3]    else NA_real_,
  fuzzy_p_robust  = if (!is.null(fuzzy_obj)) fuzzy_obj$pv[3]    else NA_real_,
  fuzzy_ci_lo     = if (!is.null(fuzzy_obj)) fuzzy_obj$ci[3, 1] else NA_real_,
  fuzzy_ci_hi     = if (!is.null(fuzzy_obj)) fuzzy_obj$ci[3, 2] else NA_real_
)

write.csv(rdd_results, file.path(out_dir, "scenario1_rdd_results.csv"),
          row.names = FALSE)
message("Saved results for ", TRAINING, " to: ", out_dir)

# ---------------------------------------------------------------------------
# Readable summary of the two estimands
# ---------------------------------------------------------------------------
cat("\n\n##############################################################\n")
cat("##  SCENARIO 1 — SUMMARY\n")
cat("##############################################################\n")

summary_tbl <- data.frame(
  training      = rdd_results$training,
  ITT_estimate  = round(rdd_results$sharp_estimate, 4),
  ITT_CI        = sprintf("[%.4f, %.4f]", rdd_results$sharp_ci_lo, rdd_results$sharp_ci_hi),
  ITT_p         = round(rdd_results$sharp_p_robust, 4),
  LATE_estimate = round(rdd_results$fuzzy_estimate, 4),
  LATE_CI       = sprintf("[%.4f, %.4f]", rdd_results$fuzzy_ci_lo, rdd_results$fuzzy_ci_hi),
  LATE_p        = round(rdd_results$fuzzy_p_robust, 4),
  stringsAsFactors = FALSE
)

print(summary_tbl, row.names = FALSE)
cat("\nFull results saved to:", out_dir, "\n")
