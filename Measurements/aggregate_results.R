# =============================================================================
# aggregate_results.R  —  consolidate per-run results into master tables
# =============================================================================
# After a parameter sweep, every measurement script writes its CSVs into a
# run-specific folder:
#   Results/<scenario>/<decision>/cap_<C>/ran_<R>/thres_<T>/[<training>/]<file>.csv
#
# This script walks that tree, reads each headline CSV, attaches the run
# parameters (parsed from the folder tags) as explicit columns, and stacks them
# into one tidy "master" table per result type under Results/master/. 
# =============================================================================

library(here)
library(dplyr)

source(here("Simulation", "Helpers.R"))   # provides RESULTS_DIR

MASTER_DIR <- file.path(RESULTS_DIR, "master")
dir.create(MASTER_DIR, showWarnings = FALSE, recursive = TRUE)

# Every result file currently on disk (master folder excluded).
sep       <- .Platform$file.sep
all_files <- list.files(RESULTS_DIR, recursive = TRUE, full.names = TRUE)
all_files <- all_files[!grepl(paste0(sep, "master", sep), all_files, fixed = TRUE)]

# Pull a "<key>_<number>" tag out of a path (e.g. grab(path, "cap") -> 4).
grab <- function(path, key) {
  m <- regmatches(path, regexpr(paste0(key, "_[0-9.]+"), path))
  if (length(m) == 0) NA_real_ else as.numeric(sub(paste0(key, "_"), "", m))
}

# Decision = the path segment directly after the scenario/prioritarian folder.
grab_decision <- function(path) {
  parts <- strsplit(path, sep, fixed = TRUE)[[1]]
  idx   <- which(grepl("^(scenario_[0-9]+|prioritarian)$", parts))
  if (!length(idx)) return(NA_character_)
  nxt <- parts[idx[1] + 1]
  if (length(nxt) == 0 || is.na(nxt)) NA_character_ else nxt
}

# Read every file with a given exact basename, tag it, and stack the rows.
collect <- function(filename) {
  files <- all_files[basename(all_files) == filename]
  if (!length(files)) return(NULL)
  bind_rows(lapply(files, function(f) {
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$run_decision   <- grab_decision(f)
    df$run_capacity   <- grab(f, "cap")
    df$run_pct_random <- grab(f, "ran")
    df$run_threshold  <- grab(f, "thres")
    df$source_file    <- f
    df
  }))
}

# Build one master CSV from all files sharing an exact basename.
write_master <- function(filename, out_name) {
  df <- collect(filename)
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("  %-38s (no files found)", out_name))
    return(invisible(NULL))
  }
  run_cols <- c("run_decision", "run_capacity", "run_pct_random", "run_threshold")
  df <- df[, c(run_cols,
               setdiff(names(df), c(run_cols, "source_file")),
               "source_file")]
  df <- df[do.call(order, df[run_cols]), , drop = FALSE]
  write.csv(df, file.path(MASTER_DIR, out_name), row.names = FALSE)
  message(sprintf("  %-38s %d rows / %d runs",
                  out_name, nrow(df), length(unique(df$source_file))))
}

cat("\n=== AGGREGATING RESULTS INTO", MASTER_DIR, "===\n")

# headline file on disk                       master output
write_master("scenario1_rdd_results.csv",      "scenario_1_rdd_master.csv")
write_master("scenario_2_results.csv",         "scenario_2_results_master.csv")
write_master("scenario_4_overall.csv",         "scenario_4_overall_master.csv")
write_master("scenario_4_itt_by_bin.csv",      "scenario_4_itt_by_bin_master.csv")
write_master("scenario_4_late_by_bin.csv",     "scenario_4_late_by_bin_master.csv")
write_master("scenario_4_overlap_by_bin.csv",  "scenario_4_overlap_by_bin_master.csv")
write_master("scenario_5_overall.csv",         "scenario_5_overall_master.csv")
write_master("scenario_5_agas.csv",            "scenario_5_agas_master.csv")
write_master("scenario_5_ade_by_bin.csv",      "scenario_5_ade_by_bin_master.csv")
write_master("scenario_5_overlap_by_bin.csv",  "scenario_5_overlap_by_bin_master.csv")
write_master("representativeness.csv",         "prioritarian_representativeness_master.csv")
write_master("auc_risk_vs_Y0.csv",             "prioritarian_auc_master.csv")
write_master("variation_Y0.csv",               "prioritarian_variation_Y0_master.csv")

message("\nMaster tables written to: ", MASTER_DIR)
