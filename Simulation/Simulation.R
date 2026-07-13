# =============================================================================
# Simulation.R  —  parameter-sweep driver for the four experiments
# =============================================================================
# Steps:
#   1. Load the raw ALMP dataset.
#   2. Source helpers, decision functions, and the simulation core.
#   3. Build the set of runs from whichever experiments are switched on.
#   4. For each run: simulate -> save -> run all measurement scenarios.
#      Every result is written to a folder tagged by ALL four parameters:
#        Results/<scenario>/<decision>/cap_<C>/ran_<R>/thres_<T>/[<training>/]
#      so no run overwrites another and the parameters are recoverable.
#   5. Consolidate everything into tidy master tables (Results/master/).
#
# To read results you never open folders: load the relevant master CSV and
# filter on run_decision / run_capacity / run_pct_random / run_threshold.
# =============================================================================

library(here)

# -----------------------------------------------------------------------------
# Load raw data (preprocessed based on Zezulka & Genin)
# -----------------------------------------------------------------------------
almp <- read.csv(here("Datasets", "ALMP_riskscore_treatment_clean.csv"))

# Global risk bins, shared by every script
GLOBAL_RISK_BREAKS <- seq(
  floor(min(almp$risk_score_log,   na.rm = TRUE) / 0.1) * 0.1,
  ceiling(max(almp$risk_score_log, na.rm = TRUE) / 0.1) * 0.1,
  by = 0.1
)

# -----------------------------------------------------------------------------
# Source functions
# -----------------------------------------------------------------------------
source(here("Simulation", "Helpers.R"))
source(here("Simulation", "decision_functions.R"))
source(here("Simulation", "simulation_core.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================
# Here we define the base values that are not changed across an experiment. We
# determine also the "sweeps", which are the values that change
BASE_DECISION    <- "decision_oracle"   
BASE_CAPACITIES    <- c(1.1,20)                   
BASE_PCT_RANDOM  <- 0.3
BASE_THRESHOLD   <- 0.4

# --- Sweep vectors, one per experiment ---------------------------------------
SWEEP_CAPACITY   <- c(1.1, 1.5, 2, 5, 10, 15, 20, 30)             # experiment 1
SWEEP_DECISION   <- c("decision_noisy", "decision_oracle",         # experiment 2
                 "decision_logistic",
                      "decision_no_decision", "decision_oracle_p75")
SWEEP_PCT_RANDOM <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)               # experiment 3
SWEEP_THRESHOLD  <- c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)       # experiment 4

# -----------------------------------------------------------------------------
# Run one setting
# -----------------------------------------------------------------------------
run_one <- function(decision_name, capacity, pct_random, threshold) {
  fn  <- get(decision_name, mode = "function")
  sim <- run_simulation(
    data                = almp,
    training_name       = "vocational_training",
    iapo_col            = "iapo_vocational",
    treatment_label     = "vocational",
    capacity_multiplier = capacity,
    pct_random          = pct_random,
    score_threshold     = threshold,
    decision_fn         = fn,
    seed                = 123
  )
  sim$params$decision_fn_name <- decision_name   
  sim
}

run_measurements <- function() {
  source(here("Measurements", "scenario_1.R"),   local = TRUE)
  source(here("Measurements", "scenario_2.R"),   local = TRUE)
  source(here("Measurements", "scenario_4.R"),   local = TRUE)
  source(here("Measurements", "scenario_5.R"),   local = TRUE)
  source(here("Measurements", "prioritarian.R"), local = TRUE)
}

# =============================================================================
# BUILD THE RUN LIST
# =============================================================================
# Each experiment varies one axis off the base. Everything else stays at base.
mk <- function(decision = BASE_DECISION, capacity,
               pct_random = BASE_PCT_RANDOM, threshold = BASE_THRESHOLD) {
  data.frame(decision = decision, capacity = capacity,
             pct_random = pct_random, threshold = threshold,
             stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Helper: repeat a sweep block once per base capacity, then stack
# -----------------------------------------------------------------------------
mk_across_capacities <- function(...) {
  do.call(rbind, lapply(BASE_CAPACITIES, function(cap) mk(capacity = cap, ...)))
}

RUN_EXPERIMENTS <- c("capacity", "decision", "pct_random", "threshold")

runs <- rbind(
  mk(capacity = SWEEP_CAPACITY),                        # experiment 1
  mk_across_capacities(decision   = SWEEP_DECISION),    # experiment 2
  mk_across_capacities(pct_random = SWEEP_PCT_RANDOM),  # experiment 3
  mk_across_capacities(threshold  = SWEEP_THRESHOLD)    # experiment 4
)

# Deduplicate: the base point is shared by every experiment, so run it once.
runs <- runs[!duplicated(runs[c("decision", "capacity", "pct_random", "threshold")]), ]
rownames(runs) <- NULL

cat(sprintf("\nExperiments %s -> %d unique settings to run.\n",
            paste(RUN_EXPERIMENTS, collapse = ", "), nrow(runs)))
print(runs)

# =============================================================================
# RUN THE SWEEP
# =============================================================================
# A failure in one setting is logged and skipped so the rest of the sweep still
# completes; failed settings are listed at the end.
failures <- list()

for (i in seq_len(nrow(runs))) {
  r   <- runs[i, ]
  tag <- sprintf("decision=%s | cap=%s | ran=%s | thres=%s",
                 r$decision, r$capacity, r$pct_random, r$threshold)
  message("\n========== RUN ", i, "/", nrow(runs), ": ", tag, " ==========")

  ok <- tryCatch({
    vt_sim <- run_one(r$decision, r$capacity, r$pct_random, r$threshold)
    save_sim(vt_sim)     # overwritten each pass, consumed immediately below
    run_measurements()   # reads the sim just saved; writes to tagged dirs
    TRUE
  }, error = function(e) {
    message("  !! RUN FAILED: ", conditionMessage(e))
    FALSE
  })

  if (!ok) failures[[length(failures) + 1]] <-
      data.frame(run = i, tag = tag, stringsAsFactors = FALSE)
}

if (length(failures)) {
  cat("\n=== SETTINGS THAT FAILED (inspect these) ===\n")
  print(do.call(rbind, failures), row.names = FALSE)
} else {
  cat("\nAll settings completed.\n")
}

# =============================================================================
# CONSOLIDATE: build Results/master/*.csv from every run on disk
# =============================================================================
source(here("Measurements", "aggregate_results.R"))
