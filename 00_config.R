# ==============================================================================
# 00_config.R — Shared Configuration & Setup
# ==============================================================================
# Description: Common setup sourced by all analysis scripts.
#              Sets working directory, loads libraries, defines constants.
#
# Usage: source("00_config.R") at the top of every analysis script
# ==============================================================================

# --- Environment Setup --------------------------------------------------------
rm(list = ls())

# Set working directory (works in RStudio and from command line)
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  project_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  # When sourced from Rscript, use the script's own location
  project_path <- getSrcDirectory(function(x) x)
  if (length(project_path) == 0 || project_path == "") {
    project_path <- getwd()
  }
}
setwd(project_path)

set.seed(784)

# --- Libraries ----------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)

# --- Constants ----------------------------------------------------------------

# Spirometry visit time points (ordered)
TIME_POINTS_SPIROMETRY <- c(
  "BASELINE", "1 MONTH", "2 MONTHS", "3 MONTHS", "4 MONTHS", "6 MONTHS"
)

# All required visits for complete-case analysis
REQUIRED_VISITS <- c(
  "BASELINE", "1 MONTH", "2 MONTHS", "3 MONTHS", "4 MONTHS", "6 MONTHS"
)

# Missingness threshold (%) — subjects above this are excluded
MISSING_THRESHOLD <- 10

# --- ARMCD Corrections --------------------------------------------------------
# Four subjects have incorrect treatment arm codes at specific visits.
# This function centralises the corrections applied across all analyses.
# NOTE: Subject IDs have been anonymised for public release.

# Define anonymised subject IDs (replace with actual IDs from trial data)
SUBJ_A <- "SUBJECT_A"  # Actual ID redacted — ARMCD should be AZD15 at 3 MONTHS
SUBJ_B <- "SUBJECT_B"  # Actual ID redacted — ARMCD should be AZD45 at 3 MONTHS
SUBJ_C <- "SUBJECT_C"  # Actual ID redacted — ARMCD should be AZD45 at 6 MONTHS
SUBJ_D <- "SUBJECT_D"  # Actual ID redacted — ARMCD should be PLACEBO at 4 MONTHS

fix_armcd <- function(df) {
  df[df$USUBJID == SUBJ_A & df$AVISIT == "3 MONTHS", "ARMCD"] <- "AZD15"
  df[df$USUBJID == SUBJ_B & df$AVISIT == "3 MONTHS", "ARMCD"] <- "AZD45"
  df[df$USUBJID == SUBJ_C & df$AVISIT == "6 MONTHS", "ARMCD"] <- "AZD45"
  df[df$USUBJID == SUBJ_D & df$AVISIT == "4 MONTHS", "ARMCD"] <- "PLACEBO"

  return(df)
}

cat("Config loaded.\n")
