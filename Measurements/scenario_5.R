# =============================================================================
# scenario_5.R  —  Average Decision Effect (ADE)
# =============================================================================
# Compares the randomly selected units (source == "random") versus the
# not-selected population, over the entire risk range.
#
#   ADE = E[Y | random] - E[Y | not_selected], 
#
# and AGAS (algorithmic ADE - random ADE).
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here("Simulation", "Helpers.R"))

TRAINING <- "vocational_training" 

# -----------------------------------------------------------------------------
# Two-proportion z-test 
# -----------------------------------------------------------------------------
prop_test_diff <- function(succ, n, conf = 0.95) {
  diff <- succ[1] / n[1] - succ[2] / n[2]
  pt   <- suppressWarnings(prop.test(succ, n, correct = FALSE, conf.level = conf))
  list(estimate = diff, ci_lo = pt$conf.int[1], ci_hi = pt$conf.int[2],
       p_value = pt$p.value, method = "two_proportion_z")
}

# -----------------------------------------------------------------------------
# Load simulation output
# -----------------------------------------------------------------------------
sim           <- load_sim(TRAINING)
decision_name <- sim$params$decision_fn_name

cap_tag <- sprintf("cap_%s", format(sim$params$capacity_multiplier, nsmall = 2))
ran_tag <- sprintf("ran_%s", format(sim$params$pct_random, nsmall = 2))
threshold_tag <- sprintf("thres_%s", format(sim$params$score_threshold, nsmall = 2))
out_dir <- file.path(RESULTS_DIR, "scenario_5", decision_name, cap_tag, ran_tag, threshold_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TAG <- run_tag(decision_name, sim$params$capacity_multiplier,
               sim$params$pct_random, sim$params$score_threshold)

# -----------------------------------------------------------------------------
# Build the two groups
# -----------------------------------------------------------------------------
# Load population (needed for KS-test)
almp <- read.csv(here("Datasets", "ALMP_riskscore_treatment_clean.csv"))

# Build control group — adjusted if the K-S test rejects representativeness.
ctrl       <- build_representative_control(sim, almp, seed = sim$params$seed + 1L)
vt_not_sel <- ctrl$control

if (ctrl$adjusted) {
  message(sprintf(
    "[scenario_5] Control group adjusted for representativeness: %d units removed from ineligible pool.",
    ctrl$n_removed
  ))
}

vt_random <- sim$selected[sim$selected$source == "random", ]

combined <- bind_rows(
  vt_random  %>% mutate(group = "random"),
  vt_not_sel %>% mutate(group = "not_selected")
)

if (sum(combined$group == "random") == 0)
  stop("No source == 'random' units in ", TRAINING,
       " — scenario 5 needs the random sub-sample (use a dataset with pct_random > 0).")
stopifnot(
  sum(combined$group == "not_selected") > 0,
  all(combined$outcome %in% c(0, 1)), !anyNA(combined$outcome)
)

cat("\n=== SCENARIO 5: AVERAGE DECISION EFFECT (random vs not_selected) ===\n")
thr <- sim$params$score_threshold
cat("Training:", TRAINING, "| Decision:", decision_name, "| Threshold:", thr, "\n")
cat("N random:", sum(combined$group == "random"),
    "| N not_selected:", sum(combined$group == "not_selected"), "\n")
cat("Mean R  random:", round(mean(combined$risk_score_log[combined$group == "random"]), 4),
    "| not_selected:", round(mean(combined$risk_score_log[combined$group == "not_selected"]), 4),
    " (differ -> overall ADE is not R-adjusted)\n")

# -----------------------------------------------------------------------------
# Overall ADE 
# H0: E[Y | random] = E[Y | not_selected]
# -----------------------------------------------------------------------------
y_rand <- combined$outcome[combined$group == "random"]
y_nsel <- combined$outcome[combined$group == "not_selected"]
n_r <- length(y_rand); n_n <- length(y_nsel)
p_r <- mean(y_rand);   p_n <- mean(y_nsel)

ade_test <- prop_test_diff(succ = c(sum(y_rand), sum(y_nsel)), n = c(n_r, n_n))
ADE <- ade_test$estimate

cat("\n=== OVERALL ADE (two-proportion z-test) ===\n")
cat("E[Y|random]:", round(p_r, 4), "| E[Y|not_selected]:", round(p_n, 4), "\n")
cat("ADE:", round(ADE, 6), "| 95% CI: [", round(ade_test$ci_lo, 6), ",",
    round(ade_test$ci_hi, 6), "] | p =", round(ade_test$p_value, 6), "\n")

# -----------------------------------------------------------------------------
# AGAS = algorithmic ADE - random ADE
# -----------------------------------------------------------------------------
# ADE_algorithmic = the ITT above the threshold (scenario_4): the outcome gap
# between algorithmic-lottery winners (source == "score") and losers
# (not_selected_above_threshold). AGAS subtracts the random-selection ADE
# computed above.
# -----------------------------------------------------------------------------
y_score <- sim$selected$outcome[sim$selected$source == "score"]   # Z = 1
y_nsat  <- sim$not_selected_above_threshold$outcome               # Z = 0

if (length(y_score) > 0 && length(y_nsat) > 0) {
  n_s <- length(y_score); p_s <- mean(y_score)
  n_a <- length(y_nsat);  p_a <- mean(y_nsat)
  ADE_alg <- p_s - p_a

  AGAS <- ADE_alg - ADE   # ADE = random-selection ADE computed above

  var_AGAS <- p_s * (1 - p_s) / n_s + p_a * (1 - p_a) / n_a +
              p_r * (1 - p_r) / n_r + p_n * (1 - p_n) / n_n
  se_AGAS  <- sqrt(var_AGAS)
  agas_ci  <- AGAS + c(-1, 1) * qnorm(0.975) * se_AGAS
  agas_p   <- 2 * pnorm(-abs(AGAS / se_AGAS))

  cat("\n=== AGAS (algorithmic ADE - random ADE) ===\n")
  cat("ADE algorithmic (ITT above threshold):", round(ADE_alg, 6),
      sprintf("  [E[Y|score]=%.4f, E[Y|nsat]=%.4f]\n", p_s, p_a))
  cat("ADE random (this script):             ", round(ADE, 6), "\n")
  cat("AGAS = ADE_alg - ADE_random:          ", round(AGAS, 6), "\n")
  cat("95% CI:", sprintf("[%.6f, %.6f]", agas_ci[1], agas_ci[2]),
      " | p =", round(agas_p, 6), "\n")

  agas_out <- data.frame(
    training        = TRAINING,
    decision        = decision_name,
    ADE_algorithmic = ADE_alg,
    ADE_random      = ADE,
    AGAS            = AGAS,
    se_AGAS         = se_AGAS,
    AGAS_ci_lo      = agas_ci[1],
    AGAS_ci_hi      = agas_ci[2],
    AGAS_p          = agas_p
  )
  write.csv(agas_out, file.path(out_dir, "scenario_5_agas.csv"), row.names = FALSE)
} else {
  cat("\n=== AGAS ===\n")
  cat("Skipped: need both score-selected and not_selected_above_threshold units",
      "(algorithmic ITT undefined in this dataset).\n")
}

# -----------------------------------------------------------------------------
# Per-bin ADE 
# -----------------------------------------------------------------------------
breaks <- if (exists("GLOBAL_RISK_BREAKS")) {
  GLOBAL_RISK_BREAKS
} else {
  rng <- range(combined$risk_score_log, na.rm = TRUE)
  seq(floor(rng[1] / 0.1) * 0.1, ceiling(rng[2] / 0.1) * 0.1, by = 0.1)
}
combined <- add_risk_bins(combined, breaks = breaks)

# -----------------------------------------------------------------------------
# Overlap / positivity diagnostic 
# -----------------------------------------------------------------------------
vt_score_all <- sim$selected[sim$selected$source == "score", ]

overlap_src <- bind_rows(
  combined     %>% transmute(risk_score_log, arm = group,  human_decision = human_decision),
  vt_score_all %>% transmute(risk_score_log, arm = "score", human_decision = human_decision)
)
overlap_src <- add_risk_bins(overlap_src, breaks = breaks)

overlap_by_bin <- overlap_src %>%
  group_by(risk_bin, arm) %>%
  summarise(
    n          = n(),
    n_approved = sum(human_decision == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = arm,
    values_from = c(n, n_approved),
    values_fill = 0
  ) %>%
  filter(!is.na(risk_bin)) %>%
  arrange(risk_bin)

# Guarantee all columns exist even if an arm is absent
for (g in c("n_random", "n_not_selected", "n_score",
            "n_approved_random", "n_approved_score"))
  if (!g %in% names(overlap_by_bin)) overlap_by_bin[[g]] <- 0L

cat("\n=== OVERLAP / POSITIVITY BY BIN (raw counts) ===\n")
print(as.data.frame(
  overlap_by_bin[, c("risk_bin",
                     "n_random",       "n_approved_random",
                     "n_not_selected",
                     "n_score",        "n_approved_score")]
))

ade_by_bin <- combined %>%
  group_by(risk_bin, group) %>%
  summarise(p = mean(outcome), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = c(p, n)) %>%
  filter(!is.na(risk_bin), !is.na(p_random), !is.na(p_not_selected),
         n_random > 0, n_not_selected > 0) %>%
  mutate(
    ADE       = p_random - p_not_selected,
    se_ADE    = sqrt(p_random * (1 - p_random) / n_random +
                       p_not_selected * (1 - p_not_selected) / n_not_selected),
    ci_lo     = ADE - qnorm(0.975) * se_ADE,
    ci_hi     = ADE + qnorm(0.975) * se_ADE,
    p_val     = 2 * pnorm(-abs(ADE / se_ADE)),
    sig_label = sig_label(p_val),
    N         = n_random + n_not_selected
  )

cat("\n=== ADE BY RISK BIN ===\n")
print(ade_by_bin[, c("risk_bin", "n_random", "n_not_selected",
                     "ADE", "ci_lo", "ci_hi", "p_val")])

# -----------------------------------------------------------------------------
# Plot: ADE per bin
# -----------------------------------------------------------------------------
size_breaks <- pretty(ade_by_bin$N, n = 4)
size_breaks <- size_breaks[size_breaks >= min(ade_by_bin$N) &
                           size_breaks <= max(ade_by_bin$N)]

p_ade <- ggplot(ade_by_bin, aes(x = risk_bin, y = ADE)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.3,
                color = "#185FA5", alpha = 0.7) +
  geom_point(aes(size = N), color = "#185FA5") +
  geom_text(aes(y = ci_hi + 0.02, label = sig_label), size = 3.5, color = "black") +
  scale_size_binned(name = "N in bin", range = c(2, 6),
                    limits = range(ade_by_bin$N, na.rm = TRUE),
                    breaks = size_breaks) +
  labs(x = "Risk score bin", y = expression(ADE ~ "/" ~ italic(D)[i] ~ tau[i])) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, sprintf("scenario_5_ade_by_bin_%s.png", TAG)), p_ade,
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Save tables
# -----------------------------------------------------------------------------
overall <- data.frame(
  training         = TRAINING,
  decision         = decision_name,
  n_random         = n_r,
  n_not_selected   = n_n,
  mean_random      = p_r,
  mean_not_selected = p_n,
  ADE              = ADE,
  ADE_ci_lo        = ade_test$ci_lo,
  ADE_ci_hi        = ade_test$ci_hi,
  ADE_p            = ade_test$p_value
)

write.csv(overall,    file.path(out_dir, "scenario_5_overall.csv"),    row.names = FALSE)
write.csv(ade_by_bin, file.path(out_dir, "scenario_5_ade_by_bin.csv"), row.names = FALSE)
write.csv(overlap_by_bin, file.path(out_dir, "scenario_5_overlap_by_bin.csv"), row.names = FALSE)

message("Saved scenario 5 results to: ", out_dir)
