# Environment and Dependencies ----

library(readxl)
library(psych)
library(boot)
library(here)
library(cluster)
library(NbClust)
library(e1071)
library(scatterplot3d)
library(ggplot2)
library(dplyr)
library(Matrix)
library(GPArotation)
library(DescTools)
library(irr)

source("krippendorff_alpha.R")
source("pca_helpers.R")
source("clustering_helpers.R")
source("visualisations.R")
source("simple_pca.R")
source("statistical_means.R")
source("helpers.R")

# Analysis Parameters and Configuration ----

set.seed(42)

n_boot              <- 1000
percentile          <- 95.0
seed                <- 42
manual_n_components <- 2 # Controls PCA dimensions (NULL = Horn's Parallel Analysis)
manual_n_clusters   <- 3 # Controls K-Means groups (NULL = NbClust consensus)

excel_file    <- here("..", "Data", "consolidated_results.xlsx")
sheet_name    <- "Whole"
sheets        <- c("CoC", "Au", "In", "Re", "La", "Fi", "Cq", "Rr", "Pc")
comp_names    <- c("Data Provenance", "Methodological Validity")
cluster_names <- c("Mechanism-Driven", "Provenance-Driven", "Legislation-Driven")
palette       <- c(
  "Legislation-Driven" = "#FF2A54",
  "Provenance-Driven"    = "#0066FF",
  "Mechanism-Driven"     = "#9B59B6"
)

# Column positions for Krippendorff's alpha across review rounds
round_cols <- list(
  Round_1 = 3:5,
  Round_2 = 6:8
)

# Preliminary Statistical Means ----

# =========================================================================
# Statistical Means
# =========================================================================

run_preliminary_stats(
  excel_file = excel_file,
  sheet_name = sheet_name
)

# =========================================================================
# Krippendorff's Alpha Reliability Pipeline
# =========================================================================

krippendorff_results <- compute_krippendorff_pipeline(
  excel_file = excel_file,
  sheets     = sheets,
  round_cols = round_cols,
  n_boot     = n_boot,
  seed       = seed
)

print(krippendorff_results)

# =========================================================================
# Ingest for PCA and Clustering
# =========================================================================

info <- extract_dataset(
  excel_file = excel_file,
  sheet_name = sheet_name
)

items_clean      <- info$items_clean
n_obs            <- info$n_obs
categories_clean <- info$categories_clean

# ========================================================================
# Run PCA
# ========================================================================

# To run the old, complex PCA that was used for small-size datasets (Redundant)
# pca_exact <- run_complex_pca(
#   items = items_clean,
#   n_boot = n_boot,
#   percentile = percentile,
#   manual_n_components = manual_n_components,
#   seed=seed
# )

# Run the new PCA tailored for larger datasets
info <- run_simple_pca(items = items_clean,
  n_obs = n_obs,
  n_boot = n_boot,
  percentile = percentile,
  manual_n_components = manual_n_components
)

pca_res <- info$pca_res
pa_res  <- info$pa_res

# =========================================================================
# Pre-Clustering Triangulation (Elbow, Silhouette, Gap Stat, Ward's)
# =========================================================================

scores_r3 <- project_pca_scores(items_clean, pca_res, comp_names)

print(head(scores_r3, 3))

triangulation_out <- triangulate_r3_clusters(scores_r3, k_range = 2:10)

# =========================================================================
# NbClust Multi-Criteria Cluster Validation (30 Indices Voting)
# =========================================================================

nb_consensus <- assess_nbclust_consensus(
  scores_matrix = scores_r3,
  min_k         = 2,
  max_k         = 10,
  method        = "kmeans"
)

# =========================================================================
# Project Publications, Cluster, and Decision Boundary Visualization
# =========================================================================

cluster_results <- characterize_pca_clusters(
  scores        = scores_r3,
  n_clusters    = manual_n_clusters,
  items         = items_clean,
  comp_names    = comp_names,
  cluster_names = cluster_names
)

cluster_results$assignments
cluster_results$score_centroids
cluster_results$item_profiles

# =========================================================================
# Publication Visualization
# =========================================================================


# For 3d plotting of 3 principal components
# s3d_plot <- plot_3d_images(
#     scores        = scores_r3,
#     clusters      = cluster_results$assignments$Cluster,
#     comp_names    = comp_names,
#     cluster_names = cluster_names,
#     palette       = palette
# )

# Calculate theoretical bounds for pca
# Compute correlation matrix of the raw items
R <- stats::cor(items_clean, use = "pairwise.complete.obs")
S <- as.matrix(pca_res$structure_matrix)

# Factor scoring weights: B = R^(-1) %*% S
B <- solve(R) %*% S

axis_limits <- get_factor_plot_limits(B, items_clean, min_val = 0, max_val = 2)
print(axis_limits)

# Include 5% boundary padding
x_bounds <- c(axis_limits$Min[1], axis_limits$Max[1]) * 1.05
y_bounds <- c(axis_limits$Min[2], axis_limits$Max[2]) * 1.05

s2d_plot <- plot_pca_clusters_2d(
  scores        = scores_r3,
  clusters      = cluster_results$assignments$Cluster,
  title         = "Projecting Publications onto Principal Components",
  comp_names    = comp_names,
  cluster_names = cluster_names,
  palette       = palette,
  xlim          = x_bounds,
  ylim          = y_bounds
)

plot_accompanying_evidence(
  cluster_results     = cluster_results,
  categories_clean    = categories_clean,
  cluster_names       = cluster_names,
  palette             = palette,
  pa_res              = pa_res,
  comp_names          = comp_names,
  manual_n_components = manual_n_components
)
