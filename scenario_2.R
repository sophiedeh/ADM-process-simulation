# =============================================================================
# scenario_2.R  —  Random vs. score selection: difference in approval rate
# =============================================================================
# From the vocational_training simulation, compare the two selected groups:
#   - source == "random"  (drawn at random from the full population)
#   - source == "score"   (drawn algorithmically from score-eligible units)
# =============================================================================
library(here)

source(here("Simulation", "Helpers.R"))

TRAINING <- "vocational_training"

# -----------------------------------------------------------------------------
# Load simulation output
# -----------------------------------------------------------------------------
sim           <- load_sim(TRAINING)
decision_name <- sim$params$decision_fn_name

cap_tag <- sprintf("cap_%s", format(sim$params$capacity_multiplier, nsmall = 2))
ran_tag <- sprintf("ran_%s", format(sim$params$pct_random, nsmall = 2))
threshold_tag <- sprintf("thres_%s", format(sim$params$score_threshold, nsmall = 2))
out_dir <- file.path(RESULTS_DIR, "scenario_2", decision_name, cap_tag, ran_tag, threshold_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Build the two groups from the selected units
# -----------------------------------------------------------------------------
selected <- sim$selected

# Outcome = human reviewer decision (1 = approved, 0 = not approved)
random_grp <- selected$human_decision[selected$source == "random"]
score_grp  <- selected$human_decision[selected$source == "score"]

# Check binary human decision
stopifnot(
  length(random_grp) > 0, length(score_grp) > 0,
  !anyNA(random_grp), !anyNA(score_grp),
  all(random_grp %in% c(0, 1)), all(score_grp %in% c(0, 1))
)

# -----------------------------------------------------------------------------
# Group means (approval rates) 
# -----------------------------------------------------------------------------
group_means <- data.frame(
  source        = c("score", "random"),
  n             = c(length(score_grp), length(random_grp)),
  n_approved    = c(sum(score_grp),    sum(random_grp)),
  approval_rate = c(mean(score_grp),   mean(random_grp))
)

cat("\n=== SCENARIO 2: APPROVAL RATE BY SELECTION SOURCE ===\n")
cat("Training:", TRAINING, "| Decision function:", decision_name, "\n\n")
print(group_means, row.names = FALSE)

diff_means <- mean(score_grp) - mean(random_grp)
cat(sprintf("\nDifference in mean outcome (score - random): %.6f\n", diff_means))

# -----------------------------------------------------------------------------
# Two-sample Kolmogorov-Smirnov test (H0: identical distributions)
# -----------------------------------------------------------------------------
cat("\n=== TWO-SAMPLE KOLMOGOROV-SMIRNOV TEST (H0: equal distributions) ===\n")

ks_res <- suppressWarnings(ks.test(score_grp, random_grp))

cat("D statistic:", round(unname(ks_res$statistic), 6), "\n")
cat("p-value:    ", round(ks_res$p.value, 6), "\n")
if (ks_res$p.value < 0.05) {
  cat("=> Reject H0 at 5%: approval-rate distributions differ between score and random groups.\n")
} else {
  cat("=> Fail to reject H0 at 5%: no significant difference in approval-rate distributions.\n")
}

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------
results <- rbind(
  data.frame(
    metric = "group_mean", source = "score",
    n = length(score_grp), n_approved = sum(score_grp),
    approval_rate = mean(score_grp),
    diff_score_minus_random = NA_real_,
    ks_D = NA_real_, ks_p = NA_real_, ks_reject_H0 = NA
  ),
  data.frame(
    metric = "group_mean", source = "random",
    n = length(random_grp), n_approved = sum(random_grp),
    approval_rate = mean(random_grp),
    diff_score_minus_random = NA_real_,
    ks_D = NA_real_, ks_p = NA_real_, ks_reject_H0 = NA
  ),
  data.frame(
    metric = "ks_test", source = NA_character_,
    n = NA_integer_, n_approved = NA_integer_, approval_rate = NA_real_,
    diff_score_minus_random = diff_means,
    ks_D = unname(ks_res$statistic), ks_p = ks_res$p.value,
    ks_reject_H0 = ks_res$p.value < 0.05
  )
)

write.csv(results, file.path(out_dir, "scenario_2_results.csv"),
          row.names = FALSE)

message("Saved scenario 2 results to: ", out_dir)
