# ==============================================================================
# 04_visualisations.R — Manuscript Figures & Supplementary Material
# ==============================================================================
# Description: Generates all figures and tables for the Brief Communication:
#
#   Figure 1:     Longitudinal cluster trajectories with CIs (FEV1, ACQ-5,
#                 NEUT, NEUTLE) stratified by treatment vs placebo
#   Figure 2:     Biomarker log-fold change heatmap with significance borders
#   Supp Table 1: Baseline demographics by treatment group
#   Supp Table 2: Baseline characteristics by trajectory cluster
#   Supp Figure 1: NLR change from baseline to 6 months
#   Supp Table 3: NLR paired test p-values
#
# Input:  .rds files + trajectory_consensus.rds from 03_trajectory_clustering.R
# Output: Figure1.pdf, Figure2.pdf, Excel tables
# ==============================================================================

source("00_config.R")
library(pheatmap)
library(openxlsx)
library(reshape2)
library(tableone)
library(lme4)

# ==============================================================================
# LOAD DATA
# ==============================================================================

trajectory   <- as.factor(readRDS("trajectory_consensus_no_baseline_blood_neut.rds"))
bloods_wide  <- readRDS("bloods_wide.rds")
bloods       <- readRDS("bloods.rds")
bloods2      <- readRDS("bloods2.rds")
clinical     <- readRDS("clinical_covars.rds")
all_data     <- readRDS("all_data.rds")

# Prepare trajectory assignments
outcomes          <- as.data.frame(trajectory)
outcomes$USUBJID  <- rownames(outcomes)

# Subset all_data columns
all_data <- all_data[, c("USUBJID", "ARMCD", "AVISIT", "FEV1", "ACQ_TOTAL")]
all_data <- merge(all_data, outcomes, by = "USUBJID")

# --- Merge blood biomarkers (EOS, NEUT, NEUTLE, EOSLE) -----------------------

bloods_alt <- bloods %>%
  filter(LBTESTCD %in% c("EOS", "EOSLE", "NEUT", "NEUTLE")) %>%
  group_by(USUBJID, AVISIT, LBTESTCD) %>%
  summarize(mean_value = mean(LBSTRESN, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = LBTESTCD, values_from = mean_value) %>%
  select(USUBJID, AVISIT, NEUTLE, NEUT, EOS, EOSLE)

bloods_alt[bloods_alt$AVISIT == "RANDOMISATION", "AVISIT"] <- "BASELINE"
all_data[all_data$AVISIT == "RANDOMISATION", "AVISIT"]     <- "BASELINE"

all_data <- merge(all_data, bloods_alt,
                  by.x = c("USUBJID", "AVISIT"),
                  by.y = c("USUBJID", "AVISIT"))

baseline <- all_data[all_data$AVISIT == "BASELINE", ]

# Order visit timepoints
all_data$AVISIT <- factor(all_data$AVISIT,
  levels = c("BASELINE", "1 MONTH", "2 MONTHS", "3 MONTHS", "4 MONTHS", "6 MONTHS")
)

# --- Define treatment groups --------------------------------------------------

all_data <- all_data %>%
  mutate(Treatment_Group = case_when(
    ARMCD == "PLACEBO"                       ~ "PLACEBO",
    ARMCD %in% c("AZD5", "AZD15", "AZD45")  ~ "TREATMENT",
    TRUE                                     ~ "OTHER"
  ))

# ==============================================================================
# FIGURE 1 — Longitudinal Cluster Trajectories
# ==============================================================================

cat("Generating Figure 1...\n")

# Reshape to long format
long_data <- all_data %>%
  pivot_longer(cols = c(FEV1, ACQ_TOTAL, NEUTLE, NEUT),
               names_to = "Variable", values_to = "Value")

long_data <- long_data %>%
  mutate(Treatment_Group = case_when(
    ARMCD == "PLACEBO"                       ~ "PLACEBO",
    ARMCD %in% c("AZD5", "AZD15", "AZD45")  ~ "TREATMENT",
    TRUE                                     ~ "OTHER"
  ))

# Calculate cluster means and 95% CIs
cluster_summary <- long_data %>%
  group_by(trajectory, AVISIT, Variable, Treatment_Group) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    lower = mean_value - 1.96 * sd(Value, na.rm = TRUE) / sqrt(n()),
    upper = mean_value + 1.96 * sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Use plotmath for formatted axis labels
cluster_summary <- cluster_summary %>%
  mutate(Variable = recode(Variable,
    "ACQ_TOTAL" = "bold('ACQ-5')",
    "FEV1"      = "bold('FEV'[1]*' (L)')",
    "NEUT"      = "bold('Neutrophil Count (\u00d710'^9*'/L)')",
    "NEUTLE"    = "bold('Neutrophil/Leukocyte (%)')"
  ))

fig1 <- ggplot(cluster_summary,
  aes(x = AVISIT, y = mean_value,
      color = factor(trajectory),
      group = interaction(trajectory, Treatment_Group))) +
  geom_line(aes(linetype = ifelse(Treatment_Group == "PLACEBO", "Placebo", "Treatment")),
            size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(trajectory)),
              alpha = 0.2, color = NA) +
  scale_linetype_manual(
    name   = "Treatment Group",
    values = c("Placebo" = "dashed", "Treatment" = "solid"),
    labels = c("Placebo (Dashed)", "Treatment (Solid)")
  ) +
  labs(title = "", x = "Visit", y = "Mean Value",
       color = "Cluster", fill = "Cluster") +
  facet_wrap(~ Variable, nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(size = 12)
  )

ggsave("Figure1.pdf", fig1, width = 12, height = 5, dpi = 300)
cat("Figure 1 saved.\n")

# ==============================================================================
# FIGURE 2 — Biomarker Log-Fold Change Heatmap
# ==============================================================================

cat("Generating Figure 2...\n")

# Merge bloods_wide with trajectory and clinical data
bloods_wide_full <- merge(baseline[, c("USUBJID", "trajectory")],
                          bloods_wide, by = "USUBJID")

clinical <- readRDS("clinical_covars.rds")
clinical <- merge(clinical, bloods_wide_full[, c("trajectory", "USUBJID")],
                  by.x = "Subject_ID", by.y = "USUBJID")

clinical_placebo <- clinical[clinical$Treatment_Group == "PLACEBO", ]

bloods_wide_full$placebo <- ifelse(
  bloods_wide_full$USUBJID %in% clinical_placebo$Subject_ID,
  "Placebo", "Treatment"
)

# Mean log-fold change by cluster × treatment
cluster_heatmap_data <- bloods_wide_full %>%
  group_by(trajectory, placebo) %>%
  summarise(across(c(ALB, ALP, ALT, AST, BASO, BASOLE, BILI, CA,
                     CREAT, CRP, EOS, EOSLE, HGB, K, LYM, LYMLE, MONO, MONOLE,
                     PLAT, RETI, SODIUM, WBC, NEUT, NEUTLE),
                   mean, na.rm = TRUE)) %>%
  ungroup()

rownames(cluster_heatmap_data) <- paste(cluster_heatmap_data$trajectory,
                                        cluster_heatmap_data$placebo, sep = "_")
heatmap_matrix <- as.matrix(cluster_heatmap_data[, -c(1, 2)])
rownames(heatmap_matrix) <- rownames(cluster_heatmap_data)

# --- Paired Wilcoxon tests for significance borders --------------------------

bloods_two_points <- readRDS("bloods_two_points.rds")
bloods_two_points <- merge(bloods_two_points, outcomes, by = "USUBJID")

df_sub <- bloods_two_points %>%
  filter(VISIT %in% c("RANDOMISATION", "6 MONTHS"))

df_wide <- df_sub %>%
  select(USUBJID, VISIT, LBTESTCD, LBSTRESN, trajectory) %>%
  pivot_wider(names_from = VISIT, values_from = LBSTRESN)

df_wide$placebo <- ifelse(
  df_wide$USUBJID %in% clinical_placebo$Subject_ID,
  "Placebo", "Treatment"
)

p_values <- df_wide %>%
  group_by(trajectory, LBTESTCD, placebo) %>%
  summarise(
    test = list(wilcox.test(`RANDOMISATION`, `6 MONTHS`, paired = TRUE, exact = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(p_value = sapply(test, function(x) x$p.value)) %>%
  select(trajectory, LBTESTCD, placebo, p_value)

result_long <- p_values %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Merge with heatmap data
heatmap_df           <- as.data.frame(heatmap_matrix)
heatmap_df$trajectory <- rownames(heatmap_df)
heatmap_df$trajectory <- as.factor(heatmap_df$trajectory)
heatmap_long         <- melt(heatmap_df, id.vars = "trajectory",
                              variable.name = "LBTESTCD", value.name = "value")

result_long <- result_long %>%
  unite("trajectory", trajectory, placebo, sep = "_", remove = FALSE)

heatmap_long <- merge(heatmap_long, result_long, by = c("trajectory", "LBTESTCD"))

# Create heatmap with significance borders
fig2 <- ggplot(heatmap_long, aes(x = LBTESTCD, y = trajectory, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0.15, name = "LogFC") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 4) +
  geom_tile(data = filter(heatmap_long, p_adj < 0.05),
            aes(x = LBTESTCD, y = trajectory),
            fill = NA, color = "black", size = 1) +
  labs(title = "", x = "Biomarker", y = "Trajectory") +
  theme_minimal()

ggsave("Figure2.pdf", fig2, width = 15, height = 5, dpi = 300)
cat("Figure 2 saved.\n")

# ==============================================================================
# MANN-WHITNEY TESTS — Treatment vs Placebo Within Clusters
# ==============================================================================

cat("Running Mann-Whitney tests...\n")

# Reload long_data (without plotmath labels)
long_data <- all_data %>%
  pivot_longer(cols = c(FEV1, ACQ_TOTAL, NEUTLE, NEUT),
               names_to = "Variable", values_to = "Value") %>%
  mutate(Treatment_Group = case_when(
    ARMCD == "PLACEBO"                       ~ "PLACEBO",
    ARMCD %in% c("AZD5", "AZD15", "AZD45")  ~ "TREATMENT",
    TRUE                                     ~ "OTHER"
  )) %>%
  filter(Treatment_Group %in% c("PLACEBO", "TREATMENT"))

mw_results <- long_data %>%
  group_by(trajectory, AVISIT, Variable) %>%
  summarise(
    test = list(wilcox.test(Value ~ Treatment_Group, exact = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(test, function(x) x$p.value),
    p_adj   = p.adjust(p_value, method = "BH")
  )

write.xlsx(mw_results[, c("trajectory", "AVISIT", "Variable", "p_value", "p_adj")],
           file = "mann_whitney_results_BH_traj.xlsx")
cat("Mann-Whitney results saved.\n")

# ==============================================================================
# PAIRED T-TESTS — Baseline vs 6 Months
# ==============================================================================

cat("Running paired t-tests (baseline vs 6 months)...\n")

baseline_6m <- all_data %>%
  filter(AVISIT %in% c("BASELINE", "6 MONTHS"))

long_data_b6 <- baseline_6m %>%
  pivot_longer(cols = c(FEV1, ACQ_TOTAL, NEUTLE, NEUT),
               names_to = "Variable", values_to = "Value") %>%
  mutate(Treatment_Group = case_when(
    ARMCD == "PLACEBO"                       ~ "PLACEBO",
    ARMCD %in% c("AZD5", "AZD15", "AZD45")  ~ "TREATMENT",
    TRUE                                     ~ "OTHER"
  )) %>%
  filter(Treatment_Group %in% c("PLACEBO", "TREATMENT"))

wide_data <- long_data_b6 %>%
  select(USUBJID, trajectory, Treatment_Group, AVISIT, Variable, Value) %>%
  pivot_wider(names_from = AVISIT, values_from = Value)

# Paired t-tests
ttest_results <- wide_data %>%
  group_by(trajectory, Treatment_Group, Variable) %>%
  summarise(
    test = list(t.test(`BASELINE`, `6 MONTHS`, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(test, function(x) x$p.value),
    p_adj   = p.adjust(p_value, method = "BH")
  )

# Formatted p-value table
formatted_pvals <- ttest_results %>%
  select(trajectory, Treatment_Group, Variable, p_adj) %>%
  pivot_wider(names_from = Variable, values_from = p_adj, names_prefix = "pval_")

write.xlsx(formatted_pvals, file = "formatted_baseline_vs_6months_pvalues.xlsx")
cat("Paired t-test results saved.\n")

# ==============================================================================
# SUPPLEMENTARY FIGURE 1 — Neutrophil-to-Lymphocyte Ratio (NLR)
# ==============================================================================

cat("Generating Supplementary Figure 1 (NLR)...\n")

# Calculate NLR from bloods2
neutlymph <- bloods2[bloods2$LBTESTCD == "NEUT", c("USUBJID", "AVISIT", "LBSTRESN")]
lymphlymph <- bloods2[bloods2$LBTESTCD == "LYM", c("USUBJID", "AVISIT", "LBSTRESN")]

# Baseline NLR
bl_neut  <- neutlymph[neutlymph$AVISIT == "ENROLMENT", ]
bl_lymph <- lymphlymph[lymphlymph$AVISIT == "ENROLMENT", ]
bl_ratio <- merge(bl_neut, bl_lymph, by = "USUBJID")
bl_ratio$ratio <- bl_ratio$LBSTRESN.x / bl_ratio$LBSTRESN.y

# 6-month NLR
m6_neut  <- neutlymph[neutlymph$AVISIT == "6 MONTHS", ]
m6_lymph <- lymphlymph[lymphlymph$AVISIT == "6 MONTHS", ]
m6_ratio <- merge(m6_neut, m6_lymph, by = "USUBJID")
m6_ratio$ratio <- m6_ratio$LBSTRESN.x / m6_ratio$LBSTRESN.y

# Merge baseline and 6-month ratios
ratio <- merge(bl_ratio[, c("USUBJID", "ratio")],
               m6_ratio[, c("USUBJID", "ratio")], by = "USUBJID")
ratio <- merge(ratio, outcomes, by = "USUBJID")
ratio$placebo <- ifelse(ratio$USUBJID %in% clinical_placebo$Subject_ID,
                        "Placebo", "Treatment")

# Reshape for plotting
df_long <- melt(ratio, id.vars = c("USUBJID", "trajectory", "placebo"),
                measure.vars = c("ratio.x", "ratio.y"))

df_long <- df_long %>%
  unite("trajectory", trajectory, placebo, sep = "_", remove = FALSE)

# Summary for line plot
summary_df <- df_long %>%
  mutate(Timepoint = recode(variable, ratio.x = "Baseline", ratio.y = "6 Months")) %>%
  group_by(trajectory, Timepoint) %>%
  summarise(mean_log = mean(value, na.rm = TRUE), .groups = "drop")

supp_fig1 <- ggplot(summary_df,
  aes(x = Timepoint, y = mean_log, group = trajectory, color = trajectory)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(x = "Timepoint", y = "Mean NLR",
       color = "Trajectory + Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.position = "right")

ggsave("Supplementary_Figure1.pdf", supp_fig1, width = 10, height = 6, dpi = 300)

# Supplementary Table 3: NLR paired tests
nlr_pvals <- ratio %>%
  group_by(trajectory, placebo) %>%
  summarise(p_value = wilcox.test(ratio.x, ratio.y, paired = TRUE)$p.value,
            .groups = "drop")

cat("NLR results:\n")
print(nlr_pvals)

cat("Supplementary Figure 1 saved.\n")

# ==============================================================================
# SUPPLEMENTARY TABLES 1 & 2 — Baseline Demographics
# ==============================================================================

cat("Generating Supplementary Tables...\n")

# Build clinical characteristics table
clinical <- readRDS("clinical_covars.rds")
clinical <- merge(clinical, bloods_wide_full[, c("trajectory", "USUBJID")],
                  by.x = "Subject_ID", by.y = "USUBJID")

# Supp Table 1: By treatment group
myvars_t1  <- c("Gender", "AGE", "Baseline_BMI", "Smoking_Status", "Baseline_OCS_USE")
catVars_t1 <- c("Gender", "AGE", "Baseline_BMI", "Smoking_Status", "Baseline_OCS_USE")

table1 <- CreateTableOne(vars = myvars_t1, strata = "Treatment_Group",
                         addOverall = TRUE, data = clinical,
                         factorVars = catVars_t1, includeNA = TRUE)
cat("\n=== Supplementary Table 1: Demographics by Treatment Group ===\n")
print(table1, showAllLevels = TRUE)

# Supp Table 2: By trajectory cluster
table2 <- CreateTableOne(vars = myvars_t1, strata = "trajectory",
                         addOverall = TRUE, data = clinical,
                         factorVars = catVars_t1, includeNA = TRUE)
cat("\n=== Supplementary Table 2: Demographics by Trajectory Cluster ===\n")
print(table2, showAllLevels = TRUE)

cat("\n=== All visualisations complete. ===\n")
