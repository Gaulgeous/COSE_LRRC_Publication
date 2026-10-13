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
# Statistical Means
# =========================================================================

# Read without default column names so both header rows (Rows 1 & 2) 
# are preserved intact as data rows for the parser
raw_data <- readxl::read_excel(
  path      = excel_file,
  sheet     = "Whole",
  col_names = FALSE
)

# Convert to standard data.frame for consistent matrix indexing
raw_df <- as.data.frame(raw_data)

# =========================================================================
# 2. Compute Category Profile Means
# =========================================================================

# Extracts Category from the Publication multi-column, coerces score
# columns to numeric, and calculates the mean for each category
category_means <- compute_category_means_simple(raw_df)

cat("--- Category Means Across Scorecard Dimensions ---\n")
print(round(category_means, 3))
cat("\n")

# =========================================================================
# 3. Pairwise Euclidean Distances & Heatmap
# =========================================================================

# Computes pairwise Euclidean distances between the 6 categories
# and renders the publication-ready heatmap
category_distances <- compute_category_euclidean_distance_simple(
  category_means = category_means,
  title          = "Euclidean Distance Between Publication Categories",
  show_values    = TRUE
)

cat("--- Pairwise Euclidean Distance Matrix ---\n")
print(round(category_distances, 3))
cat("\n")

# =========================================================================
# Krippendorff's Alpha Reliability Pipeline
# =========================================================================

# Define column positions (assuming Col 1 is Publication)
round_cols <- list(
  Round_1 = 3:5,
  Round_2 = 6:8
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

# Read starting from Row 2 so sub-headers ("Category", "Author", and the 9 items)
# become clean column names. This ensures numeric columns are parsed as numeric
# rather than being coerced to character by a text header row.
df_raw <- as.data.frame(readxl::read_excel(file_path, sheet = "Whole", skip = 1))

# Standardize column names (strip whitespace)
colnames(df_raw) <- trimws(colnames(df_raw))

# Identify which column holds the publication citation / key ("Author" or "Publication")
pub_col <- if ("Author" %in% colnames(df_raw)) {
  "Author"
} else if ("Publication" %in% colnames(df_raw)) {
  "Publication"
} else {
  colnames(df_raw)[1]
}

# Exclude Frowis_2020 outlier
cleaned_raw <- df_raw[!grepl("Frowis_2020", df_raw[[pub_col]], ignore.case = TRUE), ]

# Preserve publication keys as rownames (using standard data.frame)
rownames(cleaned_raw) <- make.unique(as.character(cleaned_raw[[pub_col]]))

# Exclude metadata columns ("Category", "Author", "Publication", etc.)
metadata_cols <- grep("^(category|author|publication)$", colnames(cleaned_raw), 
                      ignore.case = TRUE, value = TRUE)

feature_cols <- setdiff(colnames(cleaned_raw), metadata_cols)

# Extract and convert items to a numeric data frame, omitting incomplete cases
items_clean <- stats::na.omit(as.data.frame(lapply(cleaned_raw[, feature_cols, drop = FALSE], as.numeric)))
rownames(items_clean) <- rownames(cleaned_raw)[as.numeric(rownames(items_clean))]

# Calculate sample size after outlier removal (N = 89)
n_obs <- nrow(items_clean)

# =========================================================================
# 4. Verify Dimensions
# =========================================================================

# Expected: 89 rows (90 minus Frowis_2020) by 9 columns
cat(sprintf("Cleaned items dimensions: %d observations by %d scorecard items.\n", 
            nrow(items_clean), ncol(items_clean)))
print(dim(items_clean))

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
# Simplified, big dataset version
# =========================================================================
# =========================================================================
# 1. Polychoric Correlation Matrix (Exact 2-Step MLE)
# =========================================================================

# Compute exact pairwise MLE without Taylor approximations
poly_res       <- psych::polychoric(items_clean, global = FALSE, correct = 0, smooth = FALSE)
cor_polychoric <- poly_res$rho

# Conditional PSD check: project only if pairwise MLE yields a negative eigenvalue
if (min(eigen(cor_polychoric, symmetric = TRUE, only.values = TRUE)$values) < 0) {
  pd_obj         <- Matrix::nearPD(cor_polychoric, corr = TRUE, keepDiag = TRUE)
  cor_polychoric <- as.matrix(pd_obj$mat)
  dimnames(cor_polychoric) <- dimnames(poly_res$rho)
}

cat("--- Polychoric Correlation Matrix ---\n")
print(round(cor_polychoric, 3))

# Full matrix correlation heatmap
plot_correlation_heatmap(
  cor_mat = cor_polychoric,
  type    = "full",
  title   = "Polychoric Correlation of Publication's Coverage of Scorecard Components"
)

# =========================================================================
# 2. Factorability Diagnostics (Bartlett's Test & KMO)
# =========================================================================

factorability_res <- assess_factorability_simple(
  cor_mat = cor_polychoric, 
  n_obs   = n_obs, 
  verbose = TRUE
)

# =========================================================================
# 3. Dimension Retention (Scree & Horn's Parallel Analysis)
# =========================================================================

# Horn's Parallel Analysis (1,000 Monte Carlo runs with N = 90 noise)
pa_res <- perform_horn_parallel_analysis_simple(
  x               = cor_polychoric,
  n_obs           = n_obs,
  n_simulations   = 1000,
  percentile      = 95.0,
  plot_comparison = TRUE,
  verbose         = TRUE
)

# =========================================================================
# 4. PCA with Direct Oblimin Rotation & Anchor Alignment
# =========================================================================

# If manual_n_components is NULL, it will automatically fall back to pa_res$n_factors_parallel
n_comp_to_extract <- if (exists("manual_n_components") && !is.null(manual_n_components)) {
  manual_n_components
} else {
  pa_res$n_factors_parallel
}

pca_res <- perform_abels_pca_oblimin_simple(
  x             = cor_polychoric,
  n_factors     = n_comp_to_extract,
  n_obs         = n_obs,
  align_anchors = TRUE
)

cat("--- Rotated Pattern Matrix (Table 3 Style) ---\n")
print(pca_res$pattern_loadings)

cat("\n--- Structure Matrix (Table 4 Style) ---\n")
print(pca_res$structure_matrix)

cat("\n--- Component Inter-Correlation Matrix (Phi) ---\n")
print(pca_res$factor_correlations)

cat("\n--- Item-to-Component Primary Assignments ---\n")
print(pca_res$assignments_table, row.names = FALSE)

# =========================================================================
# 5. Internal Reliability (Raw Alpha, Ordinal Alpha, & McDonald's Omega)
# =========================================================================

# Self-contained: computes alpha, ordinal alpha, omega, and dropped diagnostics directly
rel_res <- compute_internal_reliability_simple(
  df      = items_clean, 
  verbose = TRUE
)



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
