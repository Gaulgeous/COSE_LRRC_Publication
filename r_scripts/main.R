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

set.seed(42)

excel_file <- here("..", "Data", "consolidated_results.xlsx")
sheets     <- c("CoC", "Au", "In", "Re", "La", "Fi", "Cq", "Rr", "Pc")
n_boot     <- 1000
percentile  <- 95.0
seed       <- 42
manual_n_components <- 3   # Controls PCA dimensions (NULL = Horn's Parallel Analysis)
manual_n_clusters   <- 3   # Controls K-Means groups   (NULL = NbClust consensus)
comp_names <- c("Data Provenance", "Process Replication", "Authority & Expertise")
cluster_names <- c("Mechanism-Driven", "Admissibility-Driven", "Provenance-Driven")
palette <- c(

  "Mechanism-Driven"     = "#FF2A54",  # Red
  "Admissibility-Driven" = "#0066FF",  # Blue
  "Provenance-Driven"    = "#9B59B6"   # Purple
)

# =========================================================================
# Krippendorff's Alpha Reliability Pipeline
# =========================================================================

# Define column positions (assuming Col 1 is Publication)
round_cols <- list(
  Round_1 = 2:4,
  Round_2 = 5:7
)

final_results <- compute_krippendorff_pipeline(
  excel_file = excel_file,
  sheets     = sheets,
  round_cols = round_cols,
  n_boot     = n_boot,
  seed       = 42
)

print(final_results)

# =========================================================================
# Ingest for PCA and Clustering
# =========================================================================
df_raw    <- read_excel(file_path, sheet="Whole")

# Filter down to complete numeric/ordinal rating items
items <- na.omit(df_raw[, sapply(df_raw, is.numeric)])
n_obs <- nrow(items)

# Filter out the Frowis_2020 as an outlier
cleaned_raw <- df_raw[df_raw$Publication != "Frowis_2020", ]

# Extract strictly the 9 features (dropping the "Publications" column)
# Using rownames so the publication labels are preserved on the observations
rownames(cleaned_raw) <- cleaned_raw$Publication
items_clean <- cleaned_raw[, setdiff(colnames(cleaned_raw), "Publication")]

# Verify dimensions (should be 26 rows by 9 columns)
dim(items_clean)

cat("===================================================================\n")
cat("DATA SUMMARY\n")
cat(sprintf("Observations (Units): %d | Variables (Items): %d\n", n_obs, ncol(items_clean)))
cat("===================================================================\n\n")

# =========================================================================
# Polychoric and Kendall Correlation Matrix (2S-MLE, Continuity Correction, Smoothing)
# =========================================================================

kendall_matrix <- compute_psd_kendall(items_clean)

cat("--- Positive Semi-Definite Kendall Tau-b Correlation Matrix ---\n")
print(round(kendall_matrix, 3))

# Compute exact 2-step MLE without Taylor approximations or continuity distortion
#    global = FALSE -> Evaluates thresholds per pair
#    correct = 0    -> Prevents artificial off-diagonal probability inflation
#    smooth = FALSE -> Prevents cor.smooth from flattening strong signals
poly_mle <- polychoric(
  items, 
  global  = FALSE, 
  correct = 0, 
  smooth  = FALSE
)$rho

# Apply Higham nearest correlation projection
cor_polychoric <- project_near_psd(poly_mle)
dimnames(cor_polychoric) <- dimnames(poly_mle)

cat("--- Polychoric Matrix (Matching Python Implementation) ---\n")
print(round(cor_polychoric, 3))

# Lower triangle heatmap (standard publication format)
# plot_correlation_heatmap(
#   cor_mat = cor_polychoric,
#   type    = "lower",
#   title   = "Polychoric Correlation of Publication's Coverage of Scorecard Components"
# )

# Full matrix heatmap without text values
plot_correlation_heatmap(
  cor_mat     = cor_polychoric,
  type        = "full",
  title = "Polychoric Correlation of Publication's Coverage of Scorecard Components"
)

# =========================================================================
# Factorability Diagnostics (Bartlett's Test & Bootstrapped KMO)
# =========================================================================

assess_factorability(cor_polychoric, n_obs)

# =========================================================================
# Dimension Retention (Kaiser-Guttman & Horn's Parallel Analysis)
# =========================================================================

# Exact Parallel Analysis
pa_res <- perform_horn_parallel_analysis_exact(
  df             = items_clean, 
  n_observations = nrow(items_clean), # 26
  n_simulations  = n_boot, 
  percentile     = percentile
)

# =========================================================================
# Principal Component Analysis with Direct Oblimin Rotation
# =========================================================================

pca_exact <- perform_abels_pca_oblimin_exact(df = items_clean, n_factors = manual_n_components)

cat("--- Rotated Pattern Matrix (Table 3 Style) ---\n")
print(round(pca_exact$pattern_loadings, 3))

cat("\n--- Structure Matrix (Table 4 Style) ---\n")
print(round(pca_exact$structure_matrix, 3))

cat("\n--- Component Inter-Correlation Matrix (Phi) ---\n")
print(round(pca_exact$factor_correlations, 3))

cat("\n--- Item-to-Component Primary Assignments ---\n")
print(pca_exact$assignments_table, row.names = FALSE)

# =========================================================================
# Internal Reliability Diagnostics (Cronbach's Alpha & McDonald's Omega)
# =========================================================================

compute_internal_reliabilty(df = items_clean, corr_matrix = cor_polychoric)

# =========================================================================
# Pre-Clustering Triangulation (Elbow, Silhouette, Gap Stat, Ward's)
# =========================================================================

# Project observations (polarity is already inverted; semantic names assigned)
scores_r3 <- project_pca_scores(items_clean, pca_exact)

# Verify structure and dimensions
print(head(scores_r3, 3))

# Run triangulation directly on the 27x3 matrix
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

# Invocation using optimal k from assess_nbclust_consensus()
cluster_results <- characterize_pca_clusters(
  scores        = scores_r3,
  n_clusters    = manual_n_clusters,
  items         = items_clean,
  comp_names    = comp_names,
  cluster_names = cluster_names
)

# Access outputs programmatically
cluster_results$assignments
cluster_results$score_centroids
cluster_results$item_profiles

# =========================================================================
# Publication Visualization
# =========================================================================

s3d_plot <- plot_pca_clusters_3d(
  scores        = scores_r3,
  clusters      = cluster_results$assignments$Cluster,
  comp_names    = comp_names,
  cluster_names = cluster_names,
  palette       = palette,
  title         = "Projecting Publications onto Principal Components",
  angle         = 130 # Optimal cluster separation azimuth
)

# For finding the perfect angle for publication
# for (theta in seq(0, 360, by = 10)) {
#   s3d_plot <- plot_pca_clusters_3d(
#     scores        = scores_r3,
#     clusters      = cluster_results$assignments$Cluster,
#     comp_names    = comp_names,
#     cluster_names = cluster_names,
#     palette       = palette,
#     title         = sprintf("Projecting Publications onto Principal Components (Angle: %d°)", theta),
#     angle         = theta
#   )
  
#   # Pause briefly so frames render visually if viewing in an interactive plot window
#   Sys.sleep(0.1)
# }

# Optional: Add publication row labels directly to the points
# text(
#   s3d_plot$xyz.convert(scores_r3[, 1], scores_r3[, 2], scores_r3[, 3]),
#   labels = rownames(scores_r3),
#   pos    = 3,
#   cex    = 0.65
# )

# -------------------------------------------------------------------------
# Tally Cluster Counts and Compute Donut Chart Geometry
# -------------------------------------------------------------------------

# Seamless integration with previously computed clusters
donut_res <- plot_cluster_donut(
  clusters = cluster_results, 
  cluster_order = cluster_names, 
  palette  = palette,
  title    = "Publication Cluster Distribution"
)

# Access count summary table
print(donut_res$cluster_counts)

# -------------------------------------------------------------------------
# Prepare Eigenvalue Data (From Horn's Analysis)
# -------------------------------------------------------------------------

plot_eigenvalue_retention(
  eigenvalues = pa_res,  
  comp_names  = comp_names,
  top_n       = manual_n_components,
  palette     = c(
    "Data Provenance"       = "#0066FF",
    "Process Replication"   = "#FF2A54",
    "Authority & Expertise" = "#9B59B6"
  ),
  title      = "Horn's Parallel Analysis: Eigenvalue Retention"
)
