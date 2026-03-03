# ==============================================================================
# 03_trajectory_clustering.R — DTW-Based Trajectory Clustering
# ==============================================================================
# Description: Identifies patient response subgroups using Dynamic Time Warping
#              (DTW) with consensus clustering on multivariate time series
#              (FEV1, ACQ-5, Neutrophil count).
#
# Method:
#   1. Exclude baseline (cluster on raw changes from 1 month onward)
#   2. Scale features, create per-subject multivariate matrices
#   3. DTW-based partitional clustering via subsampling (stability analysis)
#   4. PAM applied to stability matrix; consensus score selects optimal k
#   5. t-SNE visualisation of cluster structure
#
# Input:  .rds files (via 02_data_preparation.R)
# Output: trajectory_consensus.rds — named vector of cluster assignments
#
# Reference: Methods section of the Brief Communication manuscript
# ==============================================================================

source("02_data_preparation.R")

library(dtwclust)
library(clValid)
library(cluster)
library(sharp)
library(dtw)
library(Rtsne)

set.seed(45653)

# ==============================================================================
# 1. PREPARE MULTIVARIATE TIME SERIES
# ==============================================================================

# Convert features to numeric
filtered_df$FEV1      <- as.numeric(filtered_df$FEV1)
filtered_df$ACQ_TOTAL <- as.numeric(filtered_df$ACQ_TOTAL)
filtered_df$NEUT      <- as.numeric(filtered_df$NEUT)

# Exclude baseline — cluster on post-treatment trajectories only
filtered_df <- filtered_df[filtered_df$AVISIT != "BASELINE", ]

# Scale numeric features
filtered_df <- filtered_df %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.))))

# Save for downstream use by visualisations
saveRDS(filtered_df, "data_by_armcd.rds")

# --- Create per-subject multivariate matrices ---------------------------------
# Each subject becomes a matrix: rows = timepoints, cols = FEV1/ACQ/NEUT

flat_list <- filtered_df %>%
  group_by(USUBJID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.$USUBJID))) %>%
  map(function(subject_data) {
    subject_data %>%
      slice(-1) %>%                               # Remove first post-baseline row
      arrange(AVISIT) %>%
      as.matrix()
  })

# --- Clean matrices: extract numeric columns only -----------------------------

flat_list_cleaned <- lapply(flat_list, function(mat) {
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  # Keep AVISIT (col 2), FEV1 (col 3), ACQ (col 5), NEUT (col 6)
  df <- df[, c(2, 3, 5, 6)]
  df[[1]] <- factor(df[[1]], levels = c("2 MONTHS", "3 MONTHS", "4 MONTHS", "6 MONTHS"))
  df[[2]] <- as.numeric(df[[2]])
  df[[3]] <- as.numeric(df[[3]])
  df[[4]] <- as.numeric(df[[4]])
  df <- df[2:4]  # Keep only numeric columns
})

# --- Extract neutrophil-only trajectories for univariate clustering -----------

the_names <- names(flat_list_cleaned)
just_neut <- list()
for (x in seq_along(flat_list_cleaned)) {
  just_neut[[x]] <- flat_list_cleaned[[x]][, 3]
}
flat_list_cleaned <- just_neut
names(flat_list_cleaned) <- the_names

# ==============================================================================
# 2. STABILITY ANALYSIS — Subsampling-based Cluster Validation
# ==============================================================================

calculate_stability_list <- function(data_list, k, n_subsamples = 10,
                                     sample_fraction = 0.8) {
  n <- length(data_list)
  co_cluster_matrix <- matrix(0, n, n)
  inclusion_matrix  <- matrix(0, n, n)

  for (i in 1:n_subsamples) {
    subsample_indices <- sample(1:n, size = floor(sample_fraction * n))
    subsample_data    <- data_list[subsample_indices]

    clustering_result <- tsclust(
      subsample_data,
      type     = "partitional",
      k        = k,
      distance = "dtw_basic",
      centroid = "pam",
      trace    = FALSE
    )

    cluster_assignments <- clustering_result@cluster

    for (x in seq_along(subsample_indices)) {
      for (y in seq_along(subsample_indices)) {
        idx_x <- subsample_indices[x]
        idx_y <- subsample_indices[y]
        inclusion_matrix[idx_x, idx_y] <- inclusion_matrix[idx_x, idx_y] + 1
        if (cluster_assignments[x] == cluster_assignments[y]) {
          co_cluster_matrix[idx_x, idx_y] <- co_cluster_matrix[idx_x, idx_y] + 1
        }
      }
    }
  }

  stability_matrix <- co_cluster_matrix / inclusion_matrix
  stability_matrix[is.nan(stability_matrix)] <- 0

  avg_stability <- mean(stability_matrix[upper.tri(stability_matrix)])

  return(list(
    stability_matrix = stability_matrix,
    avg_stability    = avg_stability,
    inclusion_matrix = inclusion_matrix
  ))
}

# --- Run stability analysis for k = 2..7 -------------------------------------

cat("Running stability analysis (100 subsamples, 50% fraction)...\n")

stabilities <- list()
stab_mat    <- list()
incl_mat    <- list()

for (i in 2:7) {
  result <- calculate_stability_list(
    flat_list_cleaned, k = i, n_subsamples = 100, sample_fraction = 0.5
  )

  stabilities[[i]] <- result$avg_stability
  stab_mat[[i]]    <- result$stability_matrix
  incl_mat[[i]]    <- result$inclusion_matrix

  cat("k =", i, "| Avg stability:", round(result$avg_stability, 4), "\n")
}

# ==============================================================================
# 3. CONSENSUS CLUSTERING — Select Optimal k
# ==============================================================================

cat("Computing consensus scores...\n")

scores <- list()
cluster_assignment <- list()

for (i in 2:7) {
  herh <- stab_mat[[i]]
  ho   <- incl_mat[[i]]

  pam_clust <- cluster::pam(as.dist(1 - herh), k = i)
  theta     <- CoMembership(groups = pam_clust$clustering)

  cluster_assignment[[i]] <- pam_clust$clustering

  silhouette_values <- silhouette(pam_clust$clustering, dist(1 - herh))
  plot(silhouette_values, main = paste("Silhouette plot for k =", i))

  scores[[i]] <- ConsensusScore(
    prop  = herh[upper.tri(herh)],
    K     = ho[upper.tri(ho)],
    theta = theta[upper.tri(theta)]
  )

  cat("k =", i, "| Consensus score:", round(scores[[i]], 4), "\n")
}

# ==============================================================================
# 4. t-SNE VISUALISATION
# ==============================================================================

cat("Computing DTW distance matrix for t-SNE...\n")

time_series_list <- flat_list_cleaned
n <- length(time_series_list)
dist_matrix <- matrix(0, n, n)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    dist_matrix[i, j] <- dist_matrix[j, i] <-
      dtw(time_series_list[[i]], time_series_list[[j]])$distance
  }
}

# Best k from consensus analysis (k=3 based on manuscript)
best_k <- 3

tsne_result <- Rtsne(as.dist(dist_matrix), dims = 2,
                     is_distance = TRUE, perplexity = 50)
tsne_df <- data.frame(
  X       = tsne_result$Y[, 1],
  Y       = tsne_result$Y[, 2],
  Cluster = as.factor(cluster_assignment[[best_k]])
)

ggplot(tsne_df, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Visualisation of Trajectory Clusters",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

# ==============================================================================
# 5. SAVE CLUSTER ASSIGNMENTS
# ==============================================================================

names(cluster_assignment[[best_k]]) <- names(flat_list_cleaned)
saveRDS(cluster_assignment[[best_k]], "trajectory_consensus.rds")

cat("\n=== Trajectory clustering complete. ===\n")
cat("Optimal k =", best_k, "\n")
cat("Cluster sizes:", table(cluster_assignment[[best_k]]), "\n")
cat("Saved to: trajectory_consensus.rds\n")
