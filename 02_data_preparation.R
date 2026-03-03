# ==============================================================================
# 02_data_preparation.R — Shared Data Loading & Cleaning
# ==============================================================================
# Description: Loads all preprocessed .rds files, merges FEV1/ACQ/Neutrophil
#              time series, handles missingness, applies ARMCD corrections,
#              and filters to complete-case subjects.
#
# Input:  .rds files created by 01_data_processing.R
# Output: `filtered_df` — clean longitudinal dataframe ready for analysis
#
# This script is sourced by analysis scripts (03–08) to avoid code duplication.
# ==============================================================================

source("00_config.R")

# --- Load Preprocessed Data ---------------------------------------------------

FEV1_6         <- readRDS("FEV1_6months.rds")
FEV1           <- readRDS("spirometry.rds")
SevExa         <- readRDS("SevExRate12Mnths.rds")
AQLQ           <- readRDS("AQLQ.rds")
ACQ            <- readRDS("ACQ.rds")
NeutrophilCount <- readRDS("bloods.rds")
clinical       <- readRDS("clinical_covars.rds")

cat("All .rds files loaded.\n")

# --- Build FEV1 Time Series ---------------------------------------------------

FEV1_filtered <- FEV1[FEV1$PARAMCD == "FEV1" & FEV1$ATPT == "PRE",
                       c("USUBJID", "AVISIT", "AVAL", "ARMCD")]
colnames(FEV1_filtered)[3] <- "FEV1"

final_df <- FEV1_filtered

# Add 6-month FEV1 from separate source
FEV1_6_filtered <- FEV1_6[, c("USUBJID", "VISIT", "ZKORRES", "ARMCD")]
colnames(FEV1_6_filtered) <- c("USUBJID", "AVISIT", "FEV1", "ARMCD")
final_df <- rbind(final_df, FEV1_6_filtered)

# --- Add ACQ Total Score ------------------------------------------------------

ACQ <- ACQ %>%
  filter(AVISIT %in% TIME_POINTS_SPIROMETRY, PARAMCD == "ACQ_TOT") %>%
  select(AVAL, USUBJID, AVISIT) %>%
  rename(ACQ_TOTAL = AVAL)

final_df <- merge(final_df, ACQ,
                  by.x = c("USUBJID", "AVISIT"),
                  by.y = c("USUBJID", "AVISIT"),
                  all.x = TRUE, all.y = TRUE)

# --- Add Neutrophil Count -----------------------------------------------------
# Use NEUTLE (neutrophil/leukocyte ratio) from blood panel
# Add baseline neutrophil from clinical covariates

NeutrophilCount <- NeutrophilCount %>%
  filter(
    LBTESTCD == "NEUTLE" &
      AVISIT %in% c("BASELINE", "1 MONTH", "2 MONTHS", "3 MONTHS",
                     "4 MONTHS", "5 MONTHS", "6 MONTHS") &
      (AVISIT != "1 MONTH" | (AVISIT == "1 MONTH" & LBTPT == "PRE-DOSE"))
  ) %>%
  select(AVAL, USUBJID, AVISIT) %>%
  rename(NEUT = AVAL)

# Add baseline neutrophil index from clinical covariates
BaselineNeut <- clinical %>%
  select(Subject_ID, NeutrophilIndex) %>%
  rename(USUBJID = Subject_ID, NEUT = NeutrophilIndex) %>%
  mutate(AVISIT = "BASELINE")

NeutrophilCount <- rbind(NeutrophilCount, BaselineNeut)

final_df <- merge(final_df, NeutrophilCount,
                  by.x = c("USUBJID", "AVISIT"),
                  by.y = c("USUBJID", "AVISIT"),
                  all.x = TRUE)

# --- Remove 2-Week Visit & Clean Factors -------------------------------------

final_df <- final_df[final_df$AVISIT != "2 WEEKS", ]
final_df$AVISIT <- droplevels(final_df$AVISIT)

# --- Handle Missingness -------------------------------------------------------

# Calculate per-person missingness
final_df$missing_percentage <- rowSums(is.na(final_df)) / ncol(final_df) * 100

percentage_missing_per_person <- final_df %>%
  group_by(USUBJID) %>%
  summarise(percent_missing = mean(missing_percentage))

# Remove subjects with >10% missing data
missing_people <- percentage_missing_per_person[
  percentage_missing_per_person$percent_missing > MISSING_THRESHOLD, 1
]
final_df <- final_df[!final_df$USUBJID %in% missing_people$USUBJID, ]

cat("Removed", nrow(missing_people), "subjects with >",
    MISSING_THRESHOLD, "% missing data.\n")

# --- Apply ARMCD Corrections --------------------------------------------------

final_df <- fix_armcd(final_df)

# --- Filter to Complete Cases with All Required Visits -------------------------

visit_levels <- REQUIRED_VISITS
final_df$AVISIT <- factor(final_df$AVISIT, levels = visit_levels)
final_df <- final_df[complete.cases(final_df), ]

complete_ids <- final_df %>%
  filter(AVISIT %in% REQUIRED_VISITS) %>%
  group_by(ARMCD, USUBJID) %>%
  summarise(All_Visits = all(REQUIRED_VISITS %in% AVISIT), .groups = "drop") %>%
  filter(All_Visits) %>%
  pull(USUBJID)

filtered_df <- final_df %>%
  filter(USUBJID %in% complete_ids)

# --- Type Conversions ---------------------------------------------------------

filtered_df$FEV1 <- as.numeric(filtered_df$FEV1)
filtered_df <- filtered_df[, 1:6]  # Keep: USUBJID, AVISIT, FEV1, ARMCD, ACQ_TOTAL, NEUT

cat("Data preparation complete.", length(unique(filtered_df$USUBJID)),
    "subjects with complete data across", length(REQUIRED_VISITS), "visits.\n")
