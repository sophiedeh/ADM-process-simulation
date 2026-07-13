# =============================================================================
# 04_prioritarian_check.R  —  Prioritarian
# =============================================================================
# This script consists of the following parts:
#   1. Representativeness check: K-S test comparing the never-reviewed risk
#      score distribution against the full population.
#   2. Variation in Y(0) for never-reviewed individuals, expressed as outcome
#      proportions (with the implied Bernoulli variance p(1 - p)).
#   3. AUC of the risk score as a predictor of the no-program
#      outcome Y(0) among never-reviewed individuals.
#   4. Plot of risk score (R) against Y(0) for never-reviewed individuals.
#   5. Two "proportion treated by risk-rank" plots:
#        a. algorithmic route only        (remaining population)
#        b. full population               (algorithmic + random)
# =============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(scales)
library(pROC)

source(here("Simulation", "Helpers.R"))

TRAINING <- "vocational_training"

# -----------------------------------------------------------------------------
# Load simulation output
# -----------------------------------------------------------------------------
sim           <- load_sim(TRAINING)
DECISION_NAME <- sim$params$decision_fn_name

# Tags to create output dirs per variable variation
cap_tag       <- sprintf("cap_%s",   format(sim$params$capacity_multiplier, nsmall = 2))
ran_tag       <- sprintf("ran_%s",   format(sim$params$pct_random,          nsmall = 2))
threshold_tag <- sprintf("thres_%s", format(sim$params$score_threshold,     nsmall = 2))
out_dir <- file.path(RESULTS_DIR, "prioritarian", DECISION_NAME, cap_tag, ran_tag, threshold_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TAG <- run_tag(DECISION_NAME, sim$params$capacity_multiplier,
               sim$params$pct_random, sim$params$score_threshold)

# In simulation_core, only the selected individuals from algorithm get a source
vt_score   <- sim$selected[sim$selected$source == "score", ]   # algorithmic route
vt_not_sel <- sim$not_selected                                  # never reviewed

N_RANK_BINS <- 30   # shared binning for the two rank plots

# -----------------------------------------------------------------------------
# 1. Representativeness check and control group construction
# -----------------------------------------------------------------------------
almp <- read.csv(here("Datasets", "ALMP_riskscore_treatment_clean.csv"))
cat("=== 1. REPRESENTATIVENESS CHECK: control group vs full population ===\n")
cat(sprintf("  N (almp): %d\n", nrow(almp)))
ctrl <- build_representative_control(sim, almp, seed = sim$params$seed + 1L)
vt_not_sel <- ctrl$control   # used in all steps below
cat(sprintf("  N (control group): %d  (adjusted: %s)\n",
            nrow(vt_not_sel), ctrl$adjusted))

# ── K-S test results ─────────────────────────────────────────────────────────
# The K-S test checks for any distributional difference between the control
# group's risk scores and the full population. Rejection (p < 0.05) is the
# condition that triggers the adjustment in build_representative_control().
# Non-rejection means no significant distributional difference was detected 
# (which might still not be a full guarantee that the datasets are identical).

cat(sprintf("  K-S p-value (original): %.4f  => %s\n",
            ctrl$ks_p_original,
            ifelse(ctrl$ks_p_original < 0.05,
                   "reject H0 — distributional difference detected",
                   "fail to reject H0 — no significant distributional difference")))
if (ctrl$adjusted)
  cat(sprintf("  K-S p-value (adjusted): %.4f  => %s\n",
              ctrl$ks_p_adjusted,
              ifelse(ctrl$ks_p_adjusted < 0.05,
                     "still significantly different — interpret with caution",
                     "no significant difference after adjustment")))

# Save representativeness summary (one row per run, all test statistics included)
representativeness <- data.frame(
  training              = TRAINING,
  decision              = DECISION_NAME,
  n_almp                = nrow(almp),
  n_control             = nrow(vt_not_sel),
  control_adjusted      = ctrl$adjusted,
  pct_algo              = ctrl$pct_algo,
  n_removed             = ctrl$n_removed,
  # K-S test
  ks_p_original         = ctrl$ks_p_original,
  ks_p_adjusted         = ctrl$ks_p_adjusted,
  ks_reject_H0          = ctrl$ks_p_original < 0.05,
  ks_reject_H0_adjusted = ifelse(ctrl$adjusted, ctrl$ks_p_adjusted < 0.05, NA)
)
write.csv(representativeness, file.path(out_dir, "representativeness.csv"),
          row.names = FALSE)
# -----------------------------------------------------------------------------
# 2. Variation in Y(0) for never-reviewed
#    For a binary outcome, we determine variation by checking that a proportion
#    is not zero or one
# -----------------------------------------------------------------------------
p_hat <- mean(vt_not_sel$outcome) # because binary outcome

cat("=== 2. VARIATION IN Y(0) (never reviewed) ===\n")
cat(sprintf("  P(Y(0) = 1): %.4f  (n = %d)\n", p_hat, nrow(vt_not_sel)))
if (p_hat == 0 || p_hat == 1)
  warning("No variation in Y(0) — proportion is degenerate.")

variation_Y0 <- data.frame(
  training     = TRAINING,
  decision     = DECISION_NAME,
  n            = nrow(vt_not_sel),
  n_Y0_eq_1    = sum(vt_not_sel$outcome == 1),
  prop_Y0_eq_1 = p_hat
)
write.csv(variation_Y0, file.path(out_dir, "variation_Y0.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. AUC
# -----------------------------------------------------------------------------
roc_obj <- roc(vt_not_sel$outcome, vt_not_sel$risk_score_log, direction = "<")
auc_val <- auc(roc_obj)
auc_ci  <- ci.auc(roc_obj)

cat("=== 3. RISK SCORE vs Y(0): AUC (never reviewed) ===\n")
cat(sprintf("  AUC:    %.4f\n", as.numeric(auc_val)))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n\n", auc_ci[1], auc_ci[3]))

# Save the AUC point estimate and 95% CI
auc_results <- data.frame(
  training = TRAINING,
  decision = DECISION_NAME,
  n        = nrow(vt_not_sel),
  auc      = as.numeric(auc_val),
  ci_lo    = auc_ci[1],
  ci_hi    = auc_ci[3]
)
write.csv(auc_results, file.path(out_dir, "auc_risk_vs_Y0.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Risk score (R) vs no-program outcome Y(0) — never reviewed
# -----------------------------------------------------------------------------
p_scatter <- vt_not_sel %>%
  mutate(risk_bin = cut(risk_score_log, breaks = 30)) %>%
  group_by(risk_bin) %>%
  summarise(
    mean_outcome = mean(outcome),
    mid          = mean(risk_score_log),
    .groups      = "drop"
  ) %>%
  ggplot(aes(x = mid, y = mean_outcome)) +
  geom_point(color = "#E07B54") +
  geom_smooth(method = "lm", color = "#185FA5", se = TRUE) +
  labs(x = "R (-)", y = "Y(0) (-)") +
  theme_minimal(base_size = 13)

ggsave(file.path(out_dir, sprintf("risk_vs_outcome_%s.png", TAG)), p_scatter,
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Proportion treated (D = 1) by risk-score rank bin
# -----------------------------------------------------------------------------
# Helper: bin by rank of R and compute the treated proportion per bin.
rank_proportion <- function(df, n_bins = N_RANK_BINS) {
  df %>%
    mutate(
      rank_R   = rank(risk_score_log, ties.method = "average"),
      rank_bin = cut(rank_R, breaks = n_bins)
    ) %>%
    group_by(rank_bin) %>%
    summarise(
      n_bin        = n(),
      n_reviewed   = sum(reviewed),
      n_treated    = sum(D == 1),
      prop_treated = n_treated  / n_bin,
      review_rate  = n_reviewed / n_bin,
      mid_rank     = mean(rank_R),
      .groups      = "drop"
    )
}

# Helper: rank position of the threshold within a given population's R
# distribution.
threshold_rank <- function(df, threshold) {
  ranks <- rank(df$risk_score_log, ties.method = "average")
  below <- df$risk_score_log < threshold
  if (!any(below) || all(below)) return(NA_real_)
  # midpoint between the highest below-threshold rank and the lowest at/above
  (max(ranks[below]) + min(ranks[!below])) / 2
}

# Helper: bar plot of treated proportion against risk-rank, with optional
# threshold marker and optional review-rate alpha encoding.
rank_plot <- function(prop_df, threshold_rank_pos = NA_real_, threshold_val = NA_real_,
                      show_review_rate = TRUE) {
  if (show_review_rate) {
    p <- ggplot(prop_df, aes(x = mid_rank, y = prop_treated, alpha = review_rate)) +
      geom_col(fill = "#185FA5") +
      scale_alpha_continuous(
        name   = "Review rate\n(fraction seen by human)",
        range  = c(0.15, 0.9),
        labels = percent_format(accuracy = 1)
      )
  } else {
    p <- ggplot(prop_df, aes(x = mid_rank, y = prop_treated)) +
      geom_col(fill = "#185FA5", alpha = 0.8)
  }
  
  p <- p +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      x = "Risk score rank (left = lowest R, right = highest R)",
      y = "Proportion treated within rank-bin (D = 1) (%)"
    ) +
    theme_minimal(base_size = 13)
  
  if (!is.na(threshold_rank_pos)) {
    p <- p +
      geom_vline(xintercept = threshold_rank_pos,
                 linetype = "dashed", color = "gray40") +
      annotate("text", x = threshold_rank_pos, y = Inf,
               label = sprintf("threshold (R = %.2f)", threshold_val),
               hjust = -0.05, vjust = 1.5, size = 3.5, color = "gray30")
  }
  p
}

# 5a. Algorithmic route only --------------------------------------------------
population_algo <- merge(
  sim$remaining,
  vt_score[, c("ID", "human_decision")],
  by    = "ID",
  all.x = TRUE
) %>%
  mutate(
    D        = ifelse(!is.na(human_decision) & human_decision == 1, 1L, 0L),
    reviewed = !is.na(human_decision)
  )

prop_algo <- rank_proportion(population_algo)
trank_algo <- threshold_rank(population_algo, sim$params$score_threshold)
tval_algo  <- sim$params$score_threshold

ggsave(file.path(out_dir, sprintf("rank_proportion_treated_%s.png", TAG)),
       rank_plot(prop_algo, trank_algo, tval_algo, show_review_rate = TRUE),
       width = 10, height = 6, dpi = 300)

ggsave(file.path(out_dir, sprintf("rank_proportion_treated_no_alpha_%s.png", TAG)),
       rank_plot(prop_algo, trank_algo, tval_algo, show_review_rate = FALSE),
       width = 10, height = 6, dpi = 300)

# 5b. Full population (algorithmic + random) ----------------------------------
population_full <- bind_rows(
  sim$selected     %>% mutate(D = human_decision, reviewed = TRUE),
  sim$not_selected %>% mutate(D = 0L,             reviewed = FALSE)
)

prop_full  <- rank_proportion(population_full)
trank_full <- threshold_rank(population_full, sim$params$score_threshold)
tval_full  <- sim$params$score_threshold

ggsave(file.path(out_dir, sprintf("rank_proportion_treated_full_%s.png", TAG)),
       rank_plot(prop_full, trank_full, tval_full, show_review_rate = TRUE),
       width = 10, height = 6, dpi = 300)

ggsave(file.path(out_dir, sprintf("rank_proportion_treated_full_no_alpha_%s.png", TAG)),
       rank_plot(prop_full, trank_full, tval_full, show_review_rate = FALSE),
       width = 10, height = 6, dpi = 300)

message("Saved prioritarian-check results to: ", out_dir)