# =============================================================================
# decision_functions.R
# Human reviewer decision models.
#
# Every function has the same contract:
#   INPUT : df    — the data frame of units under review (must contain
#                   whatever column(s) the specific function needs, e.g.
#                   risk_score_log or y_exit12)
#   OUTPUT: numeric vector of approval probabilities in [0, 1], same length
#           and row order as df
# =============================================================================

# -----------------------------------------------------------------------------
# 1. No decision:
#    Approves if and only if score > threshold:"no human discretion".
# -----------------------------------------------------------------------------
decision_no_decision <- function(df, threshold = 0.4) {
  as.numeric(df$risk_score_log > threshold)
}

# -----------------------------------------------------------------------------
# 2. Logistic:
#    Transition around the threshold.
#    `steepness` controls how sharp the transition is.
# -----------------------------------------------------------------------------
decision_logistic <- function(df, threshold = 0.4, steepness = 8) {
  1 / (1 + exp(-steepness * (df$risk_score_log - threshold)))
}

# -----------------------------------------------------------------------------
# 3. Noisy:
#    Noise is added to the threshold, sometimes lead to no acceptance above
#    score and acceptance below score.
# -----------------------------------------------------------------------------
decision_noisy <- function(df, noise_sd = 0.3, threshold = 0.4) {
  base_prob <- as.numeric(df$risk_score_log > threshold)
  noisy     <- base_prob + rnorm(nrow(df), mean = 0, sd = noise_sd)
  pmin(pmax(noisy, 0), 1)   # clamp to [0, 1]
}

# -----------------------------------------------------------------------------
# 4. Outcome-oracle reviewer
#    Approves (with probability p) iff y_exit12 == 1; approves with
#    probability (1 - p) otherwise.
# -----------------------------------------------------------------------------
decision_oracle <- function(df, p = 1) {
  stopifnot(
    "y_exit12" %in% names(df),
    p >= 0, p <= 1,
    !anyNA(df$y_exit12),
    all(df$y_exit12 %in% c(0, 1))
  )
  if (p == 1) {
    return(as.numeric(df$y_exit12 == 1))
  }
  ifelse(df$y_exit12 == 1, p, 1 - p) 
}

# -----------------------------------------------------------------------------
# 5. Outcome-oracle reviewer with p = 0.75
#    Returns approval probability 0.75 when yexit12 == 1.
# -----------------------------------------------------------------------------
decision_oracle_p75 <- function(df) {
  decision_oracle(df, p = 0.75) 
}
