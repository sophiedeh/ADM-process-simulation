# =============================================================================
# simulation_core.R
# =============================================================================
# This script runs the ADM process based on the Sollicitatiescan.
# It is the script sourced by Simulation.R.

run_simulation <- function(
    data,
    training_name,
    iapo_col,
    treatment_label     = NULL,   
    base_capacity       = NULL,   
    capacity_multiplier = 1.1,
    pct_random          = 0.3,
    score_threshold     = 0.4,
    decision_fn         = decision_oracle,
    seed                = 123
) {

  # ── Validate inputs ─────────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(data),
    is.character(training_name), nchar(training_name) > 0,
    iapo_col %in% names(data),
    "iapo_no_program" %in% names(data),
    "risk_score_log"  %in% names(data),
    pct_random > 0, pct_random < 1,
    score_threshold > 0, score_threshold < 1
  )

  # Seeded once
  set.seed(seed)

  # ── Step 1: Capacity ─────────────────────────────────────────────────────────
  # Derive base_capacity from the treatment6 column indicating the observed treament

  if (!is.null(treatment_label)) {
    if (!"treatment6" %in% names(data)) {
      stop("Column 'treatment6' not found. Either add it to your data or ",
           "supply base_capacity directly.")
    }
    base_capacity <- sum(data$treatment6 == treatment_label, na.rm = TRUE)
    if (base_capacity == 0) {
      stop("No rows found in treatment6 matching label: '", treatment_label,
           "'. Check the label against your data.")
    }
    message(sprintf(
      "[%s] Observed count from treatment6 ('%s'): %d",
      training_name, treatment_label, base_capacity
    ))
  } else if (is.null(base_capacity)) {
    stop("Provide either treatment_label (to derive capacity from treatment6) ",
         "or base_capacity directly.")
  }

  capacity <- round(base_capacity * capacity_multiplier)
  n_random <- round(pct_random       * capacity)
  n_algo   <- round((1 - pct_random) * capacity)

  message(sprintf(
    "[%s] Capacity: %d | Random draw: %d | Algorithmic draw: %d",
    training_name, capacity, n_random, n_algo
  ))

  # ── Step 2: Selection ────────────────────────────────────────────────────────

  # 2a. Random draw from full population (without replacement)
  idx_random  <- sample(nrow(data), size = n_random, replace = FALSE)
  random_draw <- data[idx_random, ]
  random_draw$source <- "random"

  # 2b. Remaining population after random draw
  remaining <- data[-idx_random, ]

  # 2c. Algorithmic eligibility: score >= threshold, from remaining population
  eligible  <- remaining[remaining$risk_score_log >= score_threshold, ]
  not_elig  <- remaining[remaining$risk_score_log <  score_threshold, ]

  if (nrow(eligible) < n_algo) {
    warning(sprintf(
      "[%s] Only %d eligible rows but %d algorithmic slots — using all eligible.",
      training_name, nrow(eligible), n_algo
    ))
    n_algo <- nrow(eligible)
  }

  # 2d. Drawing random sample based on capacity from eligible set
  idx_algo           <- sample(nrow(eligible), size = n_algo, replace = FALSE)
  algo_draw          <- eligible[ idx_algo, ]
  algo_draw$source   <- "score"
  not_selected_above <- eligible[-idx_algo, ]   # eligible but not drawn

  # 2e. Combine selected and not-selected
  selected     <- rbind(random_draw, algo_draw)
  not_selected <- data[!(rownames(data) %in% rownames(selected)), ]

  # ── Verification ─────────────────────────────────────────────────────────────
  combined_check <- rbind(not_elig, not_selected_above)
  if (!setequal(rownames(not_selected), rownames(combined_check))) {
    warning("[", training_name, "] not_selected row-name check FAILED.")
  }

  # ── Step 3: Human reviewer decisions ────────────────────────────────────────
  # decision_fn receives the full "selected" data frame so it can use y_exit12 
  # for the outcome-oracle reviewer andrisk_score_log for the score-based 
  # reviewers. From the decision function approval probs per individuals as 
  # vector are registered.
  approval_probs          <- decision_fn(selected)
  selected$human_decision <- as.integer(runif(nrow(selected)) < approval_probs)

  # ── Step 4: Assign observed outcomes ────────────────────────────────────────
  # selected: if treated, draw from treatment IAPO; if not approved, draw 
  # from no-program IAPO
  selected$outcome <- rbinom(
    nrow(selected),
    size = 1, # each individual is a Bernoulli process like a coin toss
    prob = ifelse(
      selected$human_decision == 1, 
      selected[[iapo_col]],
      selected$iapo_no_program
    )
  )
  
  # not_selected: no review, binary draw from no-program IAPO
  not_selected$outcome <- rbinom(
    nrow(not_selected),
    size = 1,
    prob = not_selected$iapo_no_program
  )
  
  # not_selected_above
  not_selected_above <- not_selected[
    rownames(not_selected) %in% rownames(not_selected_above), ]
  not_selected_above$human_decision <- 0L
  not_selected_above$source         <- NA_character_

  # ── Step 5: Summary ─────────────────────────────────────────────────────────
  message(sprintf(
    "[%s] Selected: %d (treated: %d, not treated: %d) | Not selected: %d | Total: %d",
    training_name,
    nrow(selected),
    sum(selected$human_decision == 1),
    sum(selected$human_decision == 0),
    nrow(not_selected),
    nrow(selected) + nrow(not_selected)
  ))
  
  message(sprintf(
    "[%s] Remaining — above threshold (R >= %.2f): %d | below threshold (R < %.2f): %d",
    training_name,
    score_threshold,
    sum(remaining$risk_score_log >= score_threshold),
    score_threshold,
    sum(remaining$risk_score_log <  score_threshold)
  ))

  if (nrow(selected) + nrow(not_selected) != nrow(data)) {
    warning("[", training_name, "] Total row count mismatch — check selection logic.")
  }

  # ── Return ───────────────────────────────────────────────────────────────────
  list(
    name                         = training_name,
    selected                     = selected,
    not_selected                 = not_selected,
    not_selected_above_threshold = not_selected_above,
    remaining                    = remaining,
    params = list(
      iapo_col         = iapo_col,
      treatment_label  = treatment_label,
      base_capacity    = base_capacity,
      capacity_multiplier = capacity_multiplier,
      capacity         = capacity,
      n_random         = n_random,
      n_algo           = n_algo,
      pct_random       = pct_random,
      score_threshold  = score_threshold,
      decision_fn_name = deparse(substitute(decision_fn)),
      seed             = seed
    )
  )
}
