# =============================================================================
# 04_data_exploration_and_cleaning.R  —  Dataset exploration and cleaning
# =============================================================================
# Script:
#   1. Loads and inspects the raw ALMP dataset.
#   2. Summarises treatment effect statistics per programme (ATE, variance,
#      correlation with risk score, share benefiting).
#   3. Derives observed programme counts from treatment6.
#   4. Selects the relevant columns and saves the cleaned dataset to
#      Datasets/ALMP_riskscore_treatment_clean.csv
# =============================================================================

library(here)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load raw data
# -----------------------------------------------------------------------------
raw_path <- here("Datasets", "1203_ALMP_effects_riskFemale001.csv")
df <- read.csv(raw_path)

cat("Rows:", nrow(df), "| Columns:", ncol(df), "\n")
cat("Column names:\n")
print(names(df))

# -----------------------------------------------------------------------------
# 2. Treatment effect analysis per programme
# -----------------------------------------------------------------------------
# IAPO = individual average potential outcome.
# Treatment effect for individual i = iapo_programme - iapo_no_program.
# Negative TE = programme reduces unemployment probability (beneficial).

risk_col       <- "risk_score_log"
no_program_col <- "iapo_no_program"
iapo_cols      <- grep("^iapo_", names(df), value = TRUE)
treatment_cols <- setdiff(iapo_cols, no_program_col)

te_results <- lapply(treatment_cols, function(col) {
  te <- df[[col]] - df[[no_program_col]]
  data.frame(
    Treatment     = sub("iapo_", "", col),
    ATE           = mean(te),
    Variance_TE   = var(te),
    Corr_Risk_TE  = cor(df[[risk_col]], te),
    Share_Benefit = mean(te < 0)   # share for whom programme is beneficial
  )
})

te_df <- do.call(rbind, te_results)
te_df <- te_df[order(-te_df$Variance_TE), ]   # sort by heterogeneity desc

cat("\n=== TREATMENT EFFECT SUMMARY (sorted by variance) ===\n")
print(te_df, row.names = FALSE, digits = 4)

# -----------------------------------------------------------------------------
# 3. Observed programme counts  (from treatment6 column)
# -----------------------------------------------------------------------------
cat("\n=== OBSERVED PROGRAMME COUNTS (from treatment6) ===\n")
programme_counts <- as.data.frame(table(Programme = df$treatment6))
programme_counts <- programme_counts[order(-programme_counts$Freq), ]
names(programme_counts)[2] <- "Count"
print(programme_counts, row.names = FALSE)

cat("\nUse these counts as base_capacity in run_simulation():\n")
for (i in seq_len(nrow(programme_counts))) {
  cat(sprintf("  %-15s  %d\n",
              programme_counts$Programme[i],
              programme_counts$Count[i]))
}

# -----------------------------------------------------------------------------
# 4. Select columns and save cleaned dataset
# -----------------------------------------------------------------------------
individual_cols <- c(
  # Demographics
  "ID", "age", "female", "swiss", "foreigner", "foreigner_married",
  "married", "other_mother_tongue",

  # Location
  "city", "city_big", "city_medium", "city_no",
  "canton_french", "canton_german", "canton_italian", "canton_moth_tongue",
  "unemp_rate", "gdp_pc",

  # Education and qualifications
  "qual", "qual_degree", "qual_semiskilled", "qual_unskilled", "qual_wo_degree",

  # Employment history
  "emp_share_last_2yrs", "emp_spells_5yrs", "ue_spells_last_2yrs",
  "past_income", "employability",

  # Previous job characteristics
  "prev_job", "prev_job_manager", "prev_job_skilled", "prev_job_unskilled",
  "prev_job_self", "prev_job_sec1", "prev_job_sec2", "prev_job_sec3",
  "prev_job_sec_cat", "prev_job_sec_mis",

  # Training indicator
  "training",

  # Treatment and outcomes
  "treatment6", "elap", "y_emp", "y_unemp", "y_exit12",

  # Selected IAPOs
  "iapo_no_program", "iapo_vocational", "iapo_job_search",

  # Risk score
  "risk_score_log"
)

# Keep only columns that actually exist
missing_cols <- setdiff(individual_cols, names(df))
if (length(missing_cols) > 0) {
  warning("These columns were not found and will be skipped: ",
          paste(missing_cols, collapse = ", "))
}
cols_to_keep <- intersect(individual_cols, names(df))
df_clean <- df[, cols_to_keep]

cat("\n=== CLEANED DATASET ===\n")
cat("Shape:", nrow(df_clean), "rows x", ncol(df_clean), "columns\n")
print(head(df_clean))

# Save
out_path <- here("Datasets", "ALMP_riskscore_treatment_clean.csv")
write.csv(df_clean, out_path, row.names = FALSE)
cat("\nSaved cleaned dataset to:", out_path, "\n")
