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

set.seed(42)

# =========================================================================
# 1. Ingest Data
# =========================================================================
file_path <- here("..", "Data", "consolidated_results.xlsx")
df_raw    <- read_excel(file_path)

# Filter down to complete numeric/ordinal rating items
items <- na.omit(df_raw[, sapply(df_raw, is.numeric)])
n_obs <- nrow(items)

cat("===================================================================\n")
cat("DATA SUMMARY\n")
cat(sprintf("Observations (Units): %d | Variables (Items): %d\n", n_obs, ncol(items)))
cat("===================================================================\n\n")

# =========================================================================
# 2. Polychoric Correlation Matrix (2S-MLE, Continuity Correction, Smoothing)
# =========================================================================

spearman_ci <- cor.ci(items, method = "spearman", n.iter = 2000, alpha = 0.05, plot = FALSE)
kendall_ci  <- cor.ci(items, method = "kendall",  n.iter = 2000, alpha = 0.05, plot = FALSE)

cat("--- Spearman Correlation Matrix ---\n")
print(round(spearman_ci$rho, 3))
cat("\n")

cat("--- Kendall Correlation Matrix ---\n")
print(round(kendall_ci$rho, 3))
cat("\n")


poly_fit       <- polychoric(items, ML = FALSE, correct = 0.5, smooth = TRUE)
cor_polychoric <- poly_fit$rho

cat("--- Polychoric Correlation Matrix ---\n")
print(round(cor_polychoric, 3))
cat("\n")

# =========================================================================
# 3. Factorability Diagnostics (Bartlett's Test & Bootstrapped KMO)
# =========================================================================
cat("--- Bartlett's Test of Sphericity ---\n")
bartlett_res <- cortest.bartlett(cor_polychoric, n = n_obs)
cat(sprintf("Chi-Square: %.3f | df: %d | p-value: %.3e\n\n", 
            bartlett_res$chisq, bartlett_res$df, bartlett_res$p.value))

cat("--- Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy ---\n")
kmo_res <- KMO(cor_polychoric)

cat(sprintf("Overall MSA: %.3f\n\n", kmo_res$MSA))

kmo_table <- data.frame(
  Item = names(kmo_res$MSAi),
  MSA  = round(kmo_res$MSAi, 3),
  row.names = NULL
)

cat("Per-Item Measures of Sampling Adequacy (MSAi):\n")
print(kmo_table)
cat("\n")

# =========================================================================
# 4. Dimension Retention (Kaiser-Guttman & Horn's Parallel Analysis)
# =========================================================================
cat("--- Eigenvalue Decomposition & Kaiser-Guttman Rule ---\n")
ev       <- eigen(cor_polychoric, symmetric = TRUE, only.values = TRUE)$values
pct_var  <- (ev / sum(ev)) * 100
cum_var  <- cumsum(pct_var)
n_kaiser <- sum(ev >= 1.0)

variance_table <- data.frame(
  Component     = paste("Component", seq_along(ev)),
  Eigenvalue    = round(ev, 4),
  Pct_Variance  = round(pct_var, 2),
  Cum_Variance  = round(cum_var, 2),
  Retain_Kaiser = ifelse(ev >= 1.0, "Yes", "No")
)
print(variance_table, row.names = FALSE)
cat(sprintf("\nKaiser-Guttman suggests retaining %d component(s) (Cumulative Variance: %.2f%%).\n\n",
            n_kaiser, cum_var[n_kaiser]))

cat("--- Horn's Parallel Analysis (Monte Carlo Simulation) ---\n")
pa_out <- fa.parallel(
  cor_polychoric,
  n.obs  = n_obs,
  fa     = "pc",
  n.iter = 1000,
  quant  = 0.95,
  plot   = TRUE,
  main   = "Horn's Parallel Analysis vs. Scree Trace"
)

# Extract empirical and simulated eigenvalues safely
emp_ev <- pa_out$pc.values
sim_mu <- pa_out$pc.sim

# Safely locate the 95th percentile values (depends on psych internal structure)
sim_q <- if (!is.null(pa_out$sim.quantile$pc.quant)) {
  pa_out$sim.quantile$pc.quant
} else if (!is.null(pa_out$pc.sim.quantile)) {
  pa_out$pc.sim.quantile
} else {
  sim_mu # Fallback to mean simulation if quantile slot wasn't populated
}

# Ensure numeric vector
sim_q  <- as.numeric(sim_q)
sim_mu <- as.numeric(sim_mu)

pa_table <- data.frame(
  Component       = paste("Component", seq_along(emp_ev)),
  Empirical_EV    = round(emp_ev, 4),
  Synthetic_95Pct = round(sim_q, 4),
  Synthetic_Mean  = round(sim_mu, 4),
  Retain          = ifelse(emp_ev > sim_q, "Yes", "No")
)

print(pa_table, row.names = FALSE)
cat(sprintf("\nHorn's Parallel Analysis suggests retaining %d component(s).\n\n", pa_out$ncomp))

# =========================================================================
# 5. Principal Component Analysis with Direct Oblimin Rotation
# =========================================================================
MANUAL_N_COMPONENTS <- 3   # Controls PCA dimensions (NULL = Horn's Parallel Analysis)
MANUAL_N_CLUSTERS   <- 2   # Controls K-Means groups   (NULL = NbClust consensus)

n_comp <- if (!is.null(MANUAL_N_COMPONENTS)) MANUAL_N_COMPONENTS else pa_out$ncomp

cat(sprintf("--- Principal Component Analysis (Oblimin Rotation, nfactors = %d) ---\n", n_comp))

pca_fit <- principal(
  r        = cor.smooth(cor_polychoric),
  nfactors = n_comp,
  rotate   = "oblimin",
  covar    = FALSE
)

# Extract and format Pattern Matrix
pattern_mat <- unclass(pca_fit$loadings)
colnames(pattern_mat) <- paste0("RC", seq_len(n_comp))

cat("\nPattern Matrix (Direct Oblimin Loadings):\n")
print(round(pattern_mat, 3))

# Extract and format Structure Matrix
if (!is.null(pca_fit$Structure)) {
  structure_mat <- unclass(pca_fit$Structure)
  colnames(structure_mat) <- paste0("RC", seq_len(n_comp))
  cat("\nStructure Matrix (Pattern * Phi):\n")
  print(round(structure_mat, 3))
}

# Extract and format Factor Correlation Matrix (Phi)
if (!is.null(pca_fit$Phi)) {
  phi_mat <- pca_fit$Phi
  dimnames(phi_mat) <- list(colnames(pattern_mat), colnames(pattern_mat))
  cat("\nComponent Inter-Correlation Matrix (Phi):\n")
  print(round(phi_mat, 3))
}

# Assign items to primary component
primary_comp <- colnames(pattern_mat)[max.col(abs(pattern_mat))]
assignments  <- data.frame(
  Item               = rownames(pattern_mat),
  Primary_Component  = primary_comp,
  Primary_Loading    = round(pattern_mat[cbind(seq_along(primary_comp), match(primary_comp, colnames(pattern_mat)))], 3)
)

cat("\nItem-to-Component Primary Assignments:\n")
print(assignments, row.names = FALSE)

# =========================================================================
# 6. Psychometric Reliability Analysis (Alpha, Omega, Item-Dropped)
# =========================================================================
cat("===================================================================\n")
cat("PSYCHOMETRIC RELIABILITY ANALYSIS\n")
cat("===================================================================\n\n")

# A. Standard Cronbach's Alpha (Raw Data)
std_alpha <- psych::alpha(items, check.keys = FALSE)

# B. Ordinal Cronbach's Alpha (Polychoric Correlation Matrix)
ord_alpha <- psych::alpha(cor_polychoric, check.keys = FALSE)

# C. McDonald's Omega (1-Factor MINRES Solution)
om_fit <- suppressWarnings(
  omega(cor_polychoric, nfactors = 1, fm = "minres", plot = FALSE)
)

# Print Global Summary Metrics
cat("--- Scale Reliability Summary ---\n")
cat(sprintf("Standard Cronbach's Alpha (Raw) : %.3f\n", std_alpha$total$raw_alpha))
cat(sprintf("Ordinal Cronbach's Alpha (Poly) : %.3f\n", ord_alpha$total$raw_alpha))
cat(sprintf("McDonald's Omega Total (ω_t)     : %.3f\n\n", om_fit$omega.tot))

# -------------------------------------------------------------------------
# Single-Factor Loadings and Error Variances
# -------------------------------------------------------------------------
cat("--- Single-Factor Standardized Loadings (MINRES) ---\n")
loadings_table <- data.frame(
  Item             = colnames(items),
  Std_Loading_g    = round(as.vector(om_fit$schmid$sl[, "g"]), 3),
  Error_Variance_u = round(as.vector(om_fit$schmid$sl[, "u2"]), 3)
)
print(loadings_table, row.names = FALSE)
cat("\n")

# -------------------------------------------------------------------------
# Item-Deleted Diagnostics Table
# -------------------------------------------------------------------------
# Calculate Omega if an item is dropped across sub-matrices
omega_dropped <- sapply(seq_len(ncol(cor_polychoric)), function(i) {
  sub_r <- cor_polychoric[-i, -i]
  fit_sub <- suppressWarnings(omega(sub_r, nfactors = 1, fm = "minres", plot = FALSE))
  fit_sub$omega.tot
})

diagnostics_table <- data.frame(
  Item                     = colnames(items),
  Std_Alpha_If_Dropped     = round(std_alpha$alpha.drop[, "raw_alpha"], 3),
  Ordinal_Alpha_If_Dropped = round(ord_alpha$alpha.drop[, "raw_alpha"], 3),
  Delta_Ordinal_Alpha      = round(ord_alpha$alpha.drop[, "raw_alpha"] - ord_alpha$total$raw_alpha, 3),
  Omega_If_Dropped         = round(omega_dropped, 3),
  Delta_Omega              = round(omega_dropped - om_fit$omega.tot, 3),
  row.names                = NULL
)

cat("--- Item-Deleted Reliability Diagnostics ---\n")
print(diagnostics_table, row.names = FALSE)
cat("\n")

# =========================================================================
# 7. Pre-Clustering Triangulation (Elbow, Silhouette, Gap Stat, Ward's)
# =========================================================================
cat("===================================================================\n")
cat("PRE-CLUSTERING TRIANGULATION ACROSS ITEMS\n")
cat("===================================================================\n\n")

# 1. Angular correlation distance metric: D = sqrt(2 * (1 - R))
r_clipped   <- pmin(pmax(cor_polychoric, -1.0), 1.0)
dist_mat    <- as.dist(sqrt(pmax(2.0 * (1.0 - r_clipped), 0.0)))
dist_matrix <- as.matrix(dist_mat)
p_vars      <- ncol(cor_polychoric)

# Restrict k to a reasonable range strictly below variable count
k_candidates <- 2:min(5, p_vars - 1)
elbow_k_list <- 1:min(5, p_vars - 1)

# -------------------------------------------------------------------------
# A. Elbow Method (Total Within-Cluster Sum of Squares / Inertia)
# -------------------------------------------------------------------------
wss <- sapply(elbow_k_list, function(k) {
  kmeans(dist_mat, centers = k, nstart = 25)$tot.withinss
})

# Identify knee point via distance to chord line connecting (k_1, y_1) to (k_max, y_max)
x_pts <- elbow_k_list
y_pts <- wss
chord_vec  <- c(x_pts[length(x_pts)] - x_pts[1], y_pts[length(y_pts)] - y_pts[1])
chord_unit <- chord_vec / sqrt(sum(chord_vec^2))

diff_pts   <- cbind(x_pts - x_pts[1], y_pts - y_pts[1])
proj_len   <- diff_pts %*% chord_unit
proj_pts   <- proj_len %*% t(chord_unit)
dist_chord <- sqrt(rowSums((diff_pts - proj_pts)^2))
opt_k_elbow <- x_pts[which.max(dist_chord)]

# -------------------------------------------------------------------------
# B. Silhouette Analysis (Average Silhouette Width via cluster::silhouette)
# -------------------------------------------------------------------------
sil_scores <- sapply(k_candidates, function(k) {
  # Average linkage clustering on angular distance
  clusters <- cutree(hclust(dist_mat, method = "average"), k = k)
  ss <- silhouette(clusters, dist_mat)
  mean(ss[, "sil_width"])
})
opt_k_sil <- k_candidates[which.max(sil_scores)]

# -------------------------------------------------------------------------
# C. Tibshirani's Gap Statistic (cluster::clusGap with 1-SE rule)
# -------------------------------------------------------------------------
# clusGap handles uniform bounding-box Monte Carlo simulations natively
gap_res <- clusGap(
  dist_matrix,
  FUNcluster = kmeans,
  K.max      = max(k_candidates),
  B          = 500,
  nstart     = 25
)
gap_tab <- gap_res$Tab

# Extract Gap statistics and Standard Errors
gap_vals <- gap_tab[k_candidates, "gap"]
gap_se   <- gap_tab[k_candidates, "SE.sim"]

# 1-SE Rule: smallest k such that Gap(k) >= Gap(k+1) - SE(k+1)
opt_k_gap <- maxSE(gap_tab[, "gap"], gap_tab[, "SE.sim"], method = "firstSEmax")

# -------------------------------------------------------------------------
# D. Triangulation Summary Table
# -------------------------------------------------------------------------
triangulation_summary <- data.frame(
  Clusters         = paste0("k = ", k_candidates),
  WSS_Inertia      = round(wss[k_candidates], 3),
  Avg_Silhouette   = round(sil_scores, 3),
  Gap_Statistic    = round(gap_vals, 3),
  Gap_Std_Error    = round(gap_se, 3)
)

print(triangulation_summary, row.names = FALSE)
cat(sprintf("\nTriangulated Decisions: Elbow k = %d | Silhouette k = %d | Gap (1-SE) k = %d\n\n",
            opt_k_elbow, opt_k_sil, opt_k_gap))

# -------------------------------------------------------------------------
# E. Visual Dashboard (4 Diagnostic Panels)
# -------------------------------------------------------------------------
# Ward's Minimum Variance Hierarchical Linkage
ward_fit <- hclust(dist_mat, method = "ward.D2")

par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3, 1))

# Panel 1: Elbow Method
plot(elbow_k_list, wss, type = "b", pch = 19, col = "#8c564b", lwd = 2,
     xlab = "Number of Clusters (k)", ylab = "Total WSS (Inertia)",
     main = "Elbow Method (WCSS)")
abline(v = opt_k_elbow, col = "firebrick", lty = 2, lwd = 1.5)
grid()

# Panel 2: Silhouette Width
plot(k_candidates, sil_scores, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
     xlab = "Number of Clusters (k)", ylab = "Average Silhouette Width",
     main = "Silhouette Analysis")
abline(v = opt_k_sil, col = "firebrick", lty = 2, lwd = 1.5)
grid()

# Panel 3: Gap Statistic Curve with 1-SE Error Bars
plot(k_candidates, gap_vals, type = "b", pch = 15, col = "#2ca02c", lwd = 2,
     ylim = range(c(gap_vals - gap_se, gap_vals + gap_se)),
     xlab = "Number of Clusters (k)", ylab = "Gap(k)",
     main = "Gap Statistic (B=500)")
arrows(k_candidates, gap_vals - gap_se, k_candidates, gap_vals + gap_se,
       length = 0.05, angle = 90, code = 3, col = "#2ca02c")
abline(v = opt_k_gap, col = "firebrick", lty = 2, lwd = 1.5)
grid()

# Panel 4: Ward's Dendrogram
plot(as.dendrogram(ward_fit), horiz = TRUE,
     main = "Ward's Hierarchical Tree",
     xlab = "Ward Dissimilarity Distance")

# Reset plotting grid back to default
par(mfrow = c(1, 1))

# =========================================================================
# 8. NbClust Multi-Criteria Cluster Validation (30 Indices Voting)
# =========================================================================
cat("===================================================================\n")
cat("NBCLUST CONSENSUS CLUSTER VALIDATION (30 INDICES)\n")
cat("===================================================================\n\n")

# Project the 9 items from angular distance into an exact 3D Euclidean space
# (9 items x 3 dimensions: no singular scatter matrices, fully compatible with all 30 indices)
item_coords <- cmdscale(dist_mat, k = 3)
rownames(item_coords) <- colnames(cor_polychoric)

# Run NbClust with index = "all" across the coordinate space
nb_res <- suppressWarnings(
  NbClust(
    data     = item_coords,
    distance = "euclidean",
    min.nc   = 2,
    max.nc   = 5,
    method   = "ward.D2",
    index    = "all"
  )
)

# 1. Tally votes across all available indices
votes_raw   <- as.numeric(nb_res$Best.nc[1, ])
valid_votes <- votes_raw[!is.na(votes_raw) & is.finite(votes_raw) & votes_raw >= 2]

votes_table <- as.data.frame(table(Clusters_k = valid_votes))
colnames(votes_table) <- c("Clusters_k", "Metric_Votes")

cat("--- Majority Voting Consensus Across All Computed Indices ---\n")
print(votes_table, row.names = FALSE)

opt_k_nbclust <- as.numeric(as.character(votes_table$Clusters_k[which.max(votes_table$Metric_Votes)]))
cat(sprintf("\nNbClust Consensus Optimal Clusters: k = %d\n\n", opt_k_nbclust))

# 2. Item Partition Assignments
cat("--- Item Partition Assignments (Optimal k) ---\n")
partition_df <- data.frame(
  Item    = rownames(item_coords),
  Cluster = as.numeric(nb_res$Best.partition),
  row.names = NULL
)
print(partition_df)
cat("\n")

# =========================================================================
# 9. Project Publications, Cluster, and Decision Boundary Visualization
# =========================================================================
cat("===================================================================\n")
cat("PUBLICATION PROJECTION & CLUSTERING VISUALIZATION\n")
cat("===================================================================\n\n")

# 1. Standardize original metric ratings (Z-scores)
z_scores <- scale(items)

# 2. Project Component Scores via Ordinary Least Squares
# Pattern loadings from the Oblimin PCA fit:
L <- unclass(pca_fit$loadings)
n_components <- ncol(L)

# Least-squares projection: Scores = Z %*% (L(L'L)^-1)
proj_scores <- z_scores %*% (L %*% solve(t(L) %*% L))

comp_labels <- c("Data Provenance", "Process Replication", "Authority & Expertise")
colnames(proj_scores) <- comp_labels[seq_len(n_components)]

# 3. K-Means Clustering on Projected Component Scores

n_clusters <- if (!is.null(MANUAL_N_CLUSTERS)) {
  MANUAL_N_CLUSTERS
} else if (exists("opt_k_nbclust") && !is.na(opt_k_nbclust)) {
  opt_k_nbclust
} else {
  2
}

km_fit     <- kmeans(proj_scores, centers = n_clusters, nstart = 25)

# Cluster naming convention
cluster_names <- if (n_clusters == 2) {
  c("Process-Driven", "Authority-Driven")
} else {
  paste("Cluster", seq_len(n_clusters))
}
assigned_clusters <- factor(cluster_names[km_fit$cluster], levels = cluster_names)

# 4. Summaries & Metric Profiling
df_proj <- data.frame(
  Assigned_Cluster = assigned_clusters,
  round(proj_scores, 3)
)

cluster_centroids <- aggregate(proj_scores ~ Assigned_Cluster, data = df_proj, FUN = mean)
metric_profiles <- aggregate(as.data.frame(items), by = list(Assigned_Cluster = assigned_clusters), FUN = mean)
metric_profiles_clean <- metric_profiles
metric_profiles_clean[, -1] <- round(metric_profiles_clean[, -1], 3)

cat("--- Cluster Centroids (Component Space) ---\n")
print(cluster_centroids)
cat("\n")

cat("--- Original Item Profiles by Cluster ---\n")
print(metric_profiles_clean)
cat("\n")

# -------------------------------------------------------------------------
# 5. Publication-Ready Decision Boundary Plot (2D or 3D)
# -------------------------------------------------------------------------
palette <- c("#FF2A54", "#0066FF", "#2E6F40", "#C84B31")
pt_colors <- palette[as.numeric(assigned_clusters)]

if (n_components >= 3) {
  # --- 3D Scatter with Fitted SVM Hyperplane ---
  s3d <- scatterplot3d(
    proj_scores[, 1], proj_scores[, 2], proj_scores[, 3],
    color       = pt_colors,
    pch         = 19,
    cex.symbols = 1.3,
    xlab        = colnames(proj_scores)[1],
    ylab        = colnames(proj_scores)[2],
    zlab        = colnames(proj_scores)[3],
    main        = "3D PCA Component Projection & Clusters",
    angle       = 45,
    grid        = TRUE,
    box         = TRUE
  )
  
  # Fit linear SVM to extract decision hyperplane
  svm_mod <- svm(proj_scores[, 1:3], assigned_clusters, kernel = "linear", cost = 1)
  w <- t(svm_mod$coefs) %*% svm_mod$SV
  b <- -svm_mod$rho
  
  # Overlay plane: w1*x + w2*y + w3*z + b = 0 => z = (-w1*x - w2*y - b) / w3
  if (abs(w[3]) > 1e-4) {
    s3d$plane3d(
      Intercept = -b / w[3],
      x.coef    = -w[1] / w[3],
      y.coef    = -w[2] / w[3],
      col       = "#708090",
      lty       = "dotted"
    )
  }
  
  legend(
    "topright",
    legend = levels(assigned_clusters),
    col    = palette[seq_along(levels(assigned_clusters))],
    pch    = 19,
    title  = "Cluster Category",
    bty    = "n"
  )
  
} else {
  # --- 2D Scatter with Linear SVM Decision Boundary ---
  plot(
    proj_scores[, 1], proj_scores[, 2],
    col  = pt_colors,
    pch  = 19,
    cex  = 1.4,
    xlab = colnames(proj_scores)[1],
    ylab = colnames(proj_scores)[2],
    main = "2D Component Score Projections & Decision Boundary"
  )
  grid(lty = 3, col = "#cccccc")
  
  # Fit 2D linear SVM
  svm_mod <- svm(proj_scores[, 1:2], assigned_clusters, kernel = "linear", cost = 1)
  w <- t(svm_mod$coefs) %*% svm_mod$SV
  b <- -svm_mod$rho
  
  # Decision line: w1*x + w2*y + b = 0 => y = (-w1/w2)*x - (b/w2)
  if (abs(w[2]) > 1e-4) {
    abline(a = -b / w[2], b = -w[1] / w[2], col = "#333333", lty = 2, lwd = 1.5)
  }
  
  legend(
    "topright",
    legend = levels(assigned_clusters),
    col    = palette[seq_along(levels(assigned_clusters))],
    pch    = 19,
    title  = "Cluster Category",
    bty    = "n"
  )
}

# -------------------------------------------------------------------------
# 10. Tally Cluster Counts and Compute Donut Chart Geometry
# -------------------------------------------------------------------------
cluster_counts <- df_proj %>%
  count(Assigned_Cluster, name = "count") %>%
  mutate(
    pct      = count / sum(count) * 100,
    label    = sprintf("%.1f%%", pct),
    # Midpoint of stacked bars in polar space for inner percentage labels
    y_pos    = cumsum(count) - (count / 2)
  )

# Color mapping matching your Python script
color_map <- c(
  "Authority-Driven" = "#0066FF",
  "Process-Driven"   = "#FF2A54"
)

# -------------------------------------------------------------------------
# 11. Render Donut Chart
# -------------------------------------------------------------------------
ggplot(cluster_counts, aes(x = 2, y = count, fill = Assigned_Cluster)) +
  # Draw ring slices with white borders
  geom_col(color = "white", linewidth = 1.2, width = 0.8) +
  # Convert Cartesian bar chart into polar donut
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  # Hole radius (x < 1.6 creates donut center, x = 2 is ring)
  scale_fill_manual(values = color_map) +
  # Inner bold percentage labels
  geom_text(
    aes(y = y_pos, label = label),
    color = "white",
    size = 4.5,
    fontface = "bold"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 12, face = "bold"),
    plot.margin     = margin(10, 10, 10, 10)
  )

# -------------------------------------------------------------------------
# 13. Prepare Eigenvalue Data (First 3 Components)
# -------------------------------------------------------------------------
# Uses emp_ev extracted earlier from Horn's Parallel Analysis
df_eigen <- data.frame(
  Component  = factor(
    c("Data Provenance", "Process Replication", "Authority and Expertise"),
    levels = c("Data Provenance", "Process Replication", "Authority and Expertise")
  ),
  Eigenvalue = emp_ev[1:3],
  Color      = c("#FF2A54", "#9B27B0", "#0066FF")
)

# -------------------------------------------------------------------------
# 14. Render Publication-Grade Bar Plot
# -------------------------------------------------------------------------
ggplot(df_eigen, aes(x = Component, y = Eigenvalue, fill = Component)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  scale_fill_manual(values = setNames(df_eigen$Color, df_eigen$Component)) +
  geom_text(
    aes(label = sprintf("%.3f", Eigenvalue)),
    vjust    = -0.6,
    size     = 4,
    fontface = "bold"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(
    x = NULL,
    y = "Empirical Eigenvalue"
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text.x  = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10, color = "black"),
    plot.margin  = margin(15, 15, 15, 15)
  )

