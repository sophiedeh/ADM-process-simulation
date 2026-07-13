# =============================================================================
# scenario_4.R  —  ITT and LATE above the threshold, from the algorithmic lottery
# =============================================================================
# Estimands, each overall and per risk-score bin:
#   ITT  = E[Y|Z=1] - E[Y|Z=0]
#   LATE = ITT / first_stage          (Wald), first_stage = E[D|Z=1] - E[D|Z=0]
# =============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(AER)        # for ivreg() used to cross-check the Wald LATE

source(here("Simulation", "Helpers.R"))

TRAINING <- "vocational_training"

WEAK_FS_THRESHOLD <- 0.1   # minimum |first stage| to report a LATE

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
threshold     <- sim$params$score_threshold

cap_tag <- sprintf("cap_%s", format(sim$params$capacity_multiplier, nsmall = 2))
ran_tag <- sprintf("ran_%s", format(sim$params$pct_random, nsmall = 2))
threshold_tag <- sprintf("thres_%s", format(sim$params$score_threshold, nsmall = 2))
out_dir <- file.path(RESULTS_DIR, "scenario_4", decision_name, cap_tag, ran_tag, threshold_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TAG <- run_tag(decision_name, sim$params$capacity_multiplier,
               sim$params$pct_random, sim$params$score_threshold)

# -----------------------------------------------------------------------------
# Build the above-threshold analysis frame from the algorithmic lottery
# -----------------------------------------------------------------------------
# Z = 1: algorithmic-lottery winners (eligible/above-threshold, drawn)
vt_score   <- sim$selected[sim$selected$source == "score", ]
vt_score$Z <- 1L

# Z = 0: eligible individuals who lost the algorithmic lottery.
# outcome, human_decision and source are already set in simulation_core.R.
vt_control   <- sim$not_selected_above_threshold
vt_control$Z <- 0L

dat   <- rbind(vt_score, vt_control)
dat$D <- dat$human_decision   

cat("\n=== SCENARIO 4: ABOVE-THRESHOLD ITT / LATE FROM ALGORITHMIC LOTTERY ===\n")
cat("Training:", TRAINING, "| Decision:", decision_name,
    "| Threshold:", threshold, "\n")
cat("N lottery winners (Z=1):", sum(dat$Z == 1),
    "| N lottery losers (Z=0):", sum(dat$Z == 0), "\n")

# Balance check (means should be similar if the lottery was random)
cat("Mean R  Z=1:", round(mean(dat$risk_score_log[dat$Z == 1]), 4),
    "| Z=0:", round(mean(dat$risk_score_log[dat$Z == 0]), 4), "\n")

# -----------------------------------------------------------------------------
# Overall ITT 
# H0: E[Y|Z=1] = E[Y|Z=0]
# -----------------------------------------------------------------------------
n_Z1 <- sum(dat$Z == 1);  n_Z0 <- sum(dat$Z == 0)
p_Z1 <- mean(dat$outcome[dat$Z == 1])
p_Z0 <- mean(dat$outcome[dat$Z == 0])
se_ITT <- sqrt(p_Z1 * (1 - p_Z1) / n_Z1 + p_Z0 * (1 - p_Z0) / n_Z0)

itt_test <- prop_test_diff(
  succ = c(sum(dat$outcome[dat$Z == 1]), sum(dat$outcome[dat$Z == 0])),
  n    = c(n_Z1, n_Z0)
)
ITT <- itt_test$estimate

cat("\n=== OVERALL ITT (two-proportion z-test) ===\n")
cat("E[Y|Z=1]:", round(p_Z1, 4), "| E[Y|Z=0]:", round(p_Z0, 4), "\n")
cat("ITT:", round(ITT, 6), "| 95% CI: [", round(itt_test$ci_lo, 6), ",",
    round(itt_test$ci_hi, 6), "] | p =", round(itt_test$p_value, 6), "\n")

# Cross-check
hand_ci <- ITT + c(-1, 1) * qnorm(0.975) * se_ITT
itt_ci_matches <- isTRUE(all.equal(hand_ci, c(itt_test$ci_lo, itt_test$ci_hi),
                                    tolerance = 1e-6))
cat("Hand-CI vs prop.test-CI match:", itt_ci_matches, "\n")

# -----------------------------------------------------------------------------
# First stage 
# -----------------------------------------------------------------------------
d_Z1 <- mean(dat$D[dat$Z == 1])
d_Z0 <- mean(dat$D[dat$Z == 0])
se_pi <- sqrt(d_Z1 * (1 - d_Z1) / n_Z1 + d_Z0 * (1 - d_Z0) / n_Z0)

fs_test <- prop_test_diff(
  succ = c(sum(dat$D[dat$Z == 1]), sum(dat$D[dat$Z == 0])),
  n    = c(n_Z1, n_Z0)
)
first_stage <- fs_test$estimate

cat("\n=== FIRST STAGE (two-proportion z-test) ===\n")
cat("E[D|Z=1]:", round(d_Z1, 4), "| E[D|Z=0]:", round(d_Z0, 4), "\n")
cat("First stage:", round(first_stage, 6), "| 95% CI: [", round(fs_test$ci_lo, 6), ",",
    round(fs_test$ci_hi, 6), "] | p =", round(fs_test$p_value, 6), "\n")

# -----------------------------------------------------------------------------
# Overall LATE
# -----------------------------------------------------------------------------
cat("\n=== OVERALL LATE (Wald + delta method) ===\n")
late_matches_ivreg <- NA
if (abs(first_stage) >= WEAK_FS_THRESHOLD) {
  LATE    <- ITT / first_stage
  se_LATE <- (1 / abs(first_stage)) * sqrt(se_ITT^2 + LATE^2 * se_pi^2)
  late_ci <- LATE + c(-1, 1) * qnorm(0.975) * se_LATE
  late_p  <- 2 * pnorm(-abs(LATE / se_LATE))

  cat("LATE:", round(LATE, 6), "| 95% CI: [", round(late_ci[1], 6), ",",
      round(late_ci[2], 6), "] | p =", round(late_p, 6), "\n")

  # Cross-check against 2SLS
  iv_model <- ivreg(outcome ~ D | Z, data = dat)
  iv_coef  <- coef(iv_model)[["D"]]
  late_matches_ivreg <- isTRUE(all.equal(LATE, iv_coef, tolerance = 1e-6))

  cat("LATE point matches ivreg (2SLS)?", late_matches_ivreg,
      sprintf(" [hand = %.6f, ivreg = %.6f]\n", LATE, iv_coef))
} else {
  LATE <- NA_real_; se_LATE <- NA_real_; late_ci <- c(NA_real_, NA_real_); late_p <- NA_real_
  cat("First stage too weak (", round(first_stage, 4), ") — LATE not reported.\n")
}

# -----------------------------------------------------------------------------
# Per-bin setup
# -----------------------------------------------------------------------------
breaks <- if (exists("GLOBAL_RISK_BREAKS")) {
  GLOBAL_RISK_BREAKS
} else {
  rng <- range(dat$risk_score_log, na.rm = TRUE)
  seq(floor(rng[1] / 0.1) * 0.1, ceiling(rng[2] / 0.1) * 0.1, by = 0.1)
}
dat <- add_risk_bins(dat, breaks = breaks)

# -----------------------------------------------------------------------------
# Overlap / positivity diagnostic
# -----------------------------------------------------------------------------
overlap_by_bin <- dat %>%
  group_by(risk_bin) %>%
  summarise(
    n_Z1        = sum(Z == 1),
    n_Z0        = sum(Z == 0),
    n_treated   = sum(D == 1),
    n_untreated = sum(D == 0),
    .groups = "drop"
  ) %>%
  filter(!is.na(risk_bin))

cat("\n=== OVERLAP / POSITIVITY BY BIN (raw counts) ===\n")
print(as.data.frame(overlap_by_bin))

# -----------------------------------------------------------------------------
# Per-bin ITT 
# -----------------------------------------------------------------------------
itt_by_bin <- dat %>%
  group_by(risk_bin) %>%
  summarise(
    n_Z1 = sum(Z == 1), n_Z0 = sum(Z == 0),
    p_Z1 = mean(outcome[Z == 1]), p_Z0 = mean(outcome[Z == 0]),
    .groups = "drop"
  ) %>%
  filter(!is.na(risk_bin), n_Z1 > 0, n_Z0 > 0) %>%
  mutate(
    ITT       = p_Z1 - p_Z0,
    se_ITT    = sqrt(p_Z1 * (1 - p_Z1) / n_Z1 + p_Z0 * (1 - p_Z0) / n_Z0),
    ci_lo     = ITT - qnorm(0.975) * se_ITT,
    ci_hi     = ITT + qnorm(0.975) * se_ITT,
    p_val     = 2 * pnorm(-abs(ITT / se_ITT)),
    sig_label = sig_label(p_val),
    N         = n_Z1 + n_Z0
  )

cat("\n=== ITT BY RISK BIN ===\n")
print(itt_by_bin[, c("risk_bin", "n_Z1", "n_Z0", "ITT", "ci_lo", "ci_hi", "p_val")])

# -----------------------------------------------------------------------------
# Per-bin LATE
# -----------------------------------------------------------------------------
late_by_bin <- dat %>%
  group_by(risk_bin) %>%
  summarise(
    n_Z1 = sum(Z == 1), n_Z0 = sum(Z == 0),
    p_Y_Z1 = mean(outcome[Z == 1]), p_Y_Z0 = mean(outcome[Z == 0]),
    d_Z1   = mean(D[Z == 1]),        d_Z0   = mean(D[Z == 0]),
    .groups = "drop"
  ) %>%
  filter(!is.na(risk_bin), n_Z1 > 0, n_Z0 > 0) %>%
  mutate(
    ITT_k     = p_Y_Z1 - p_Y_Z0,
    se_ITT_k  = sqrt(p_Y_Z1 * (1 - p_Y_Z1) / n_Z1 + p_Y_Z0 * (1 - p_Y_Z0) / n_Z0),
    pi_k      = d_Z1 - d_Z0,
    se_pi_k   = sqrt(d_Z1 * (1 - d_Z1) / n_Z1 + d_Z0 * (1 - d_Z0) / n_Z0),
    weak_fs   = abs(pi_k) < WEAK_FS_THRESHOLD,
    LATE_k    = ifelse(weak_fs, NA_real_, ITT_k / pi_k),
    se_LATE_k = ifelse(weak_fs, NA_real_,
                       (1 / abs(pi_k)) * sqrt(se_ITT_k^2 + LATE_k^2 * se_pi_k^2)),
    ci_lo     = LATE_k - qnorm(0.975) * se_LATE_k,
    ci_hi     = LATE_k + qnorm(0.975) * se_LATE_k,
    p_val     = 2 * pnorm(-abs(LATE_k / se_LATE_k)),
    sig_label = sig_label(p_val),
    N         = n_Z1 + n_Z0
  )

cat("\n=== LATE BY RISK BIN ===\n")
print(late_by_bin[, c("risk_bin", "n_Z1", "n_Z0", "pi_k",
                      "ITT_k", "LATE_k", "ci_lo", "ci_hi", "weak_fs")])

n_weak <- sum(late_by_bin$weak_fs, na.rm = TRUE)
if (n_weak > 0)
  message(n_weak, " bin(s) omitted from the LATE plot (|first stage| < ",
          WEAK_FS_THRESHOLD, ").")

# Shared point-size scale
n_counts   <- itt_by_bin$N
size_breaks <- pretty(n_counts, n = 4)
size_breaks <- size_breaks[size_breaks >= min(n_counts) & size_breaks <= max(n_counts)]
size_scale <- scale_size_binned(
  name   = "N in bin",
  range  = c(2, 6),
  limits = range(n_counts, na.rm = TRUE),
  breaks = size_breaks
)

# -----------------------------------------------------------------------------
# Plot 1: ITT per bin
# -----------------------------------------------------------------------------
p_itt <- ggplot(itt_by_bin, aes(x = risk_bin, y = ITT)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.3,
                color = "#185FA5", alpha = 0.7) +
  geom_point(aes(size = N), color = "#185FA5") +
  geom_text(aes(y = ci_hi + 0.02, label = sig_label), size = 3.5, color = "black") +
  size_scale +
  labs(x = "Risk score bin", y = expression(italic(D)[i] ~ tau[i])) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, sprintf("scenario_4_itt_by_bin_%s.png", TAG)), p_itt,
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Plot 2: LATE per bin (weak first-stage bins omitted)
# -----------------------------------------------------------------------------
p_late <- ggplot(late_by_bin %>% filter(!weak_fs), aes(x = risk_bin, y = LATE_k)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.3,
                color = "#185FA5", alpha = 0.7) +
  geom_point(aes(size = N), color = "#185FA5") +
  geom_text(aes(y = ci_hi + 0.02, label = sig_label), size = 3.5, color = "black") +
  size_scale +
  labs(x = "Risk score bin", y = expression(tau[complier])) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, sprintf("scenario_4_late_by_bin_%s.png", TAG)), p_late,
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Save tables
# -----------------------------------------------------------------------------
overall <- data.frame(
  training         = TRAINING,
  decision         = decision_name,
  threshold        = threshold,
  n_offered        = n_Z1,
  n_control        = n_Z0,
  ITT              = ITT,
  ITT_ci_lo        = itt_test$ci_lo,
  ITT_ci_hi        = itt_test$ci_hi,
  ITT_p            = itt_test$p_value,
  first_stage      = first_stage,
  first_stage_ci_lo = fs_test$ci_lo,
  first_stage_ci_hi = fs_test$ci_hi,
  first_stage_p    = fs_test$p_value,
  LATE             = LATE,
  LATE_ci_lo       = late_ci[1],
  LATE_ci_hi       = late_ci[2],
  LATE_p           = late_p,
  itt_ci_matches_hand_calc = itt_ci_matches,
  late_matches_ivreg       = late_matches_ivreg
)

write.csv(overall,     file.path(out_dir, "scenario_4_overall.csv"),     row.names = FALSE)
write.csv(itt_by_bin,  file.path(out_dir, "scenario_4_itt_by_bin.csv"),  row.names = FALSE)
write.csv(late_by_bin, file.path(out_dir, "scenario_4_late_by_bin.csv"), row.names = FALSE)
write.csv(overlap_by_bin, file.path(out_dir, "scenario_4_overlap_by_bin.csv"), row.names = FALSE)

message("Saved scenario 4 results to: ", out_dir)
