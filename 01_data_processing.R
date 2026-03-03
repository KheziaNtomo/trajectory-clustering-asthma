# ==============================================================================
# 01_data_processing.R — Raw Data Ingestion
# ==============================================================================
# Description: Reads raw SAS (.sas7bdat) files from the Data/data directory
#              and extracts key clinical datasets, saving each as .rds files.
#
# Input:  Data/data/*.sas7bdat (raw clinical trial data)
# Output: Multiple .rds files used by downstream scripts
#
# This script only needs to be run ONCE (or whenever source data changes).
# ==============================================================================

source("00_config.R")
library(haven)

# --- Read Raw SAS Data --------------------------------------------------------

file_path <- "Data/data"
file_list <- list.files(path = file_path, pattern = "\\.sas7bdat$", full.names = TRUE)
data_list <- lapply(file_list, read_sas)

cat("Loaded", length(file_list), "SAS files.\n")

# --- Extract Spirometry Data (FEV1) -------------------------------------------

d1 <- data_list[[1]]
fev1 <- d1[d1$PARAMCD == "FEV1", ]
saveRDS(fev1, "FEV1.rds")

# Extract spirometry data (all parameters, split by category)
spirometry_split <- d1 %>%
  unite("PARCAT", PARCAT1:PARCAT5, na.rm = TRUE, sep = "") %>%
  group_split(PARCAT)

# Save full spirometry for downstream use
saveRDS(d1, "spirometry.rds")

cat("Spirometry data saved.\n")

# --- Extract Exacerbation & Hospital Endpoints --------------------------------

d2 <- data_list[[2]]

# Metrics to extract individually
metrics <- c(
  "Severe Exacerbation Rate over 12 Months",
  "Time to First Sev. Ex. in 6mths (Days)",
  "No. Hosp/ICU/ER Ad. in 12 Months",
  "No. Mild/Mod. Ex. with OCS <3 days in 12mths",
  "No. of Severe Exacerbations in 12 Months"
)
sub_names <- c(
  "SevExRate12Mnths", "Time to first exac", "HospAdm",
  "MildModOCS12Mon", "SevExa12Mnth"
)

for (i in seq_along(metrics)) {
  subset_data <- d2 %>% filter(PARAM == metrics[[i]])
  saveRDS(subset_data, paste0(sub_names[[i]], ".rds"))
}

# 6-month hospital/ICU/ER and mean duration endpoints
round2_met <- d2[d2$PARAM %in% c(
  "Mean Duration Sev. Ex. in 6mths (Days)",
  "No. Hosp/ICU/ER Ad. in 6 Months",
  "No. of Severe Exacerbations in 6 Months "
), ]
saveRDS(round2_met, "important_hosp.rds")

# Time to first exacerbation
time2first <- d2[d2$PARAM %in% c("Time to First Sev. Ex. in 6mths (Days)"), ]
saveRDS(time2first, "time2firstexa.rds")

cat("Exacerbation endpoints saved.\n")

# --- Extract Clinical Covariates ----------------------------------------------

d3 <- data_list[[3]]
d3 <- d3 %>%
  filter(ARMCD != "SCRNFAIL") %>%     # Remove screening failures
  filter(COMPTRFL == "Y") %>%          # Completed treatment
  filter(FASFL == "Y") %>%             # Full analysis set
  filter(COMPLFL == "Y") %>%           # Study completers
  select(USUBJID, SEX, ACTARMCD, AGE, BLBMI, SMOKSTAT,
         GROUP01, GROUP04, COVAR2, COVAR3) %>%
  rename(
    Subject_ID                   = USUBJID,
    Gender                       = SEX,
    Treatment_Group              = ACTARMCD,
    Baseline_BMI                 = BLBMI,
    Smoking_Status               = SMOKSTAT,
    OCS_USE                      = GROUP01,
    Baseline_OCS_USE             = GROUP04,
    NeutrophilIndex              = COVAR2,
    NeutrophilIndexCategorisation = COVAR3
  )
saveRDS(d3, "clinical_covars.rds")

cat("Clinical covariates saved.\n")

# --- Extract Controlled vs Uncontrolled Asthma --------------------------------

d4 <- data_list[[4]]
d4 <- d4[(d4$PARAMCD == "NOUNCW1" | d4$PARAMCD == "NOCONW1"), ]
saveRDS(d4, "controlled_vs_uncontrolled.rds")

# --- Extract Blood Biomarkers ------------------------------------------------

d18 <- data_list[[18]]

# Full blood panel (with lymphocytes)
bloods2 <- d18[d18$LBTESTCD %in% c("EOS", "EOSLE", "NEUT", "NEUTLE", "WBC", "IGE", "LYM"), ]
saveRDS(bloods2, "bloods2.rds")

# Core blood panel (without lymphocytes)
bloods <- d18[d18$LBTESTCD %in% c("EOS", "EOSLE", "NEUT", "NEUTLE", "WBC", "IGE"), ]
saveRDS(bloods, "bloods.rds")

cat("Blood biomarkers saved.\n")

# --- Extract 6-Month FEV1 from Final Dataset ----------------------------------

last_df <- data_list[[length(data_list)]]
last_df <- last_df[last_df$ZKCAT == "SPIROMETRY", ]
last_df <- last_df[last_df$ZKTESTCD == "FEV1", ]
last_df <- last_df[last_df$ZKTPT == "PRE", ]
saveRDS(last_df[last_df$VISIT == "6 MONTHS", ], "FEV1_6months.rds")

cat("6-month FEV1 saved.\n")

# --- Extract Full Blood Panel -------------------------------------------------
# (Used for imputation and fold-change analyses)

saveRDS(d18, "full_blood_panel.rds")

cat("\n=== Data processing complete. All .rds files saved. ===\n")
