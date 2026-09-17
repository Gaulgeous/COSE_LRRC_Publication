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

set.seed(42)

# =========================================================================
# 1. Ingest Data
# =========================================================================
file_path <- here("..", "Data", "consolidated_results.xlsx")
df_raw    <- read_excel(file_path, sheet="Whole")

# Filter down to complete numeric/ordinal rating items
items <- na.omit(df_raw[, sapply(df_raw, is.numeric)])
n_obs <- nrow(items)

# 1. Filter out the Frowis_2020 record from your raw dataframe (e.g., your_data)
cleaned_raw <- df_raw[df_raw$Publication != "Frowis_2020", ]

# 2. Extract strictly the 9 psychometric features (dropping the "Publications" column)
# Using rownames so the publication labels are preserved on the observations:
rownames(cleaned_raw) <- cleaned_raw$Publication
items_clean <- cleaned_raw[, setdiff(colnames(cleaned_raw), "Publication")]

# 3. Verify dimensions (should be 26 rows by 9 columns)
dim(items_clean)

cat("===================================================================\n")
cat("DATA SUMMARY\n")
cat(sprintf("Observations (Units): %d | Variables (Items): %d\n", n_obs, ncol(items)))
cat("===================================================================\n\n")

# =========================================================================
# 2. Polychoric Correlation Matrix (2S-MLE, Continuity Correction, Smoothing)
# =========================================================================

safe_kendall <- function(x, y) {
  # Return 0.0 if either vector has zero variance (all values identical)
  if (length(unique(x)) <= 1 || length(unique(y)) <= 1) {
    return(0.0)
  }
  
  val <- suppressWarnings(KendallTauB(x, y, conf.level = NA))
  if (is.na(val) || is.nan(val)) 0.0 else as.numeric(val)
}

compute_psd_kendall <- function(df) {
  cols <- colnames(df)
  p <- length(cols)
  
  # Initialize raw Kendall tau-b correlation matrix
  tau_mat <- matrix(0, nrow = p, ncol = p, dimnames = list(cols, cols))
  diag(tau_mat) <- 1.0
  
  pairs <- combn(cols, 2, simplify = FALSE)
  
  # Compute pairwise complete estimates
  for (pair in pairs) {
    u <- pair[1]
    v <- pair[2]
    
    pair_data <- na.omit(df[, c(u, v)])
    
    if (nrow(pair_data) >= 3) {
      tau_val <- safe_kendall(pair_data[[u]], pair_data[[v]])
      tau_mat[u, v] <- tau_val
      tau_mat[v, u] <- tau_val
    }
  }
  
  # Check positive semi-definiteness via minimum eigenvalue
  min_eigen <- min(eigen(tau_mat, symmetric = TRUE, only.values = TRUE)$values)
  
  # Apply Higham (2002) nearest correlation projection if non-PSD
  if (min_eigen < 0) {
    tau_psd <- as.matrix(nearPD(tau_mat, corr = TRUE, keepDiag = TRUE)$mat)
    dimnames(tau_psd) <- dimnames(tau_mat)
  } else {
    tau_psd <- tau_mat
  }
  
  return(tau_psd)
}

# --- Execution ---
# items <- ... (your 27x9 dataframe or matrix)
kendall_matrix <- compute_psd_kendall(items)

cat("--- Positive Semi-Definite Kendall Tau-b Correlation Matrix ---\n")
print(round(kendall_matrix, 3))

# 1. Compute exact 2-step MLE without Taylor approximations or continuity distortion
#    ML = TRUE      -> Equivalent to Python's log-likelihood minimization
#    global = FALSE -> Evaluates thresholds per pair (like Python's np.unique)
#    correct = 0    -> Prevents artificial off-diagonal probability inflation
#    smooth = FALSE -> Prevents cor.smooth from flattening strong signals
poly_mle <- polychoric(
  items, 
  global  = FALSE, 
  correct = 0, 
  smooth  = FALSE
)$rho

# 2. Apply Higham (2002) nearest correlation projection (identical to Python's corr_nearest)
cor_polychoric <- as.matrix(nearPD(poly_mle, corr = TRUE, keepDiag = TRUE)$mat)
dimnames(cor_polychoric) <- dimnames(poly_mle)

cat("--- Polychoric Matrix (Matching Python Implementation) ---\n")
print(round(cor_polychoric, 3))

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
# =========================================================================
# 1. Exact Kaiser-Guttman Scree Analysis (Matches Python factor_analyzer)
# =========================================================================
compute_pca_scree_analysis <- function(df_or_corr, kaiser_cutoff = 1.0) {
  if (is.matrix(df_or_corr) || is.data.frame(df_or_corr)) {
    if (nrow(df_or_corr) == ncol(df_or_corr) && all(diag(as.matrix(df_or_corr)) == 1)) {
      r_mat <- as.matrix(df_or_corr)
    } else {
      # Compute Pearson correlation matrix matching pandas .corr()
      clean_df <- na.omit(as.data.frame(lapply(df_or_corr, as.numeric)))
      r_mat <- cor(clean_df)
    }
  } else {
    stop("Input must be a dataframe of items or a correlation matrix.")
  }
  
  # Eigendecomposition sorted descending
  decomp <- eigen(r_mat, symmetric = TRUE)
  ev <- decomp$values
  
  pct_var <- (ev / sum(ev)) * 100
  cum_var <- cumsum(pct_var)
  n_kaiser <- sum(ev >= kaiser_cutoff)
  
  summary_df <- data.frame(
    Component     = paste("Component", seq_along(ev)),
    Eigenvalue    = round(ev, 4),
    Pct_Variance  = round(pct_var, 2),
    Cum_Variance  = round(cum_var, 2),
    Retain_Kaiser = ifelse(ev >= kaiser_cutoff, "Yes", "No")
  )
  
  cat("--- Eigenvalue Decomposition & Kaiser-Guttman Rule ---\n")
  print(summary_df, row.names = FALSE)
  cat(sprintf("\nKaiser-Guttman suggests retaining %d component(s) (Cumulative Variance: %.2f%%).\n\n",
              n_kaiser, cum_var[n_kaiser]))
  
  list(
    eigenvalues           = ev,
    n_factors_kaiser      = n_kaiser,
    variance_table        = summary_df,
    total_variance_kaiser = cum_var[n_kaiser]
  )
}

# =========================================================================
# 2. Exact Horn's Parallel Analysis (Matches Python NumPy Engine)
# =========================================================================
perform_horn_parallel_analysis_exact <- function(df = NULL,
                                                 corr_matrix = NULL,
                                                 n_observations = 26,
                                                 n_simulations = 1000,
                                                 percentile = 95.0,
                                                 random_state = 42,
                                                 plot_comparison = TRUE) {
  set.seed(random_state)
  
  # 1. Establish Empirical Correlation Matrix
  if (!is.null(df)) {
    clean_df <- na.omit(as.data.frame(lapply(df, as.numeric)))
    n_obs <- nrow(clean_df)
    n_vars <- ncol(clean_df)
    emp_corr <- cor(clean_df)
  } else if (!is.null(corr_matrix)) {
    emp_corr <- as.matrix(corr_matrix)
    n_obs <- n_observations
    n_vars <- ncol(emp_corr)
  } else {
    stop("Must provide either `df` or `corr_matrix`.")
  }
  
  # 2. Empirical Eigenvalues (Sorted Descending)
  emp_ev <- sort(eigen(emp_corr, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  
  # 3. Monte Carlo Simulation of Uncorrelated Gaussian Noise (N x P)
  synthetic_evals <- matrix(0, nrow = n_simulations, ncol = n_vars)
  
  for (i in seq_len(n_simulations)) {
    # Generate Gaussian white noise
    noise <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
    noise_corr <- cor(noise)
    synthetic_evals[i, ] <- sort(eigen(noise_corr, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  }
  
  # Simulated Mean and Percentile (Type 7 matches NumPy default percentile)
  sim_mu <- colMeans(synthetic_evals)
  sim_q  <- apply(synthetic_evals, 2, function(col) quantile(col, probs = percentile / 100, type = 7))
  
  # 4. Retention Decision
  retained_mask <- emp_ev > sim_q
  n_retained <- sum(retained_mask)
  
  comparison_df <- data.frame(
    Component       = paste("Component", seq_len(n_vars)),
    Empirical_EV    = round(emp_ev, 4),
    Synthetic_95Pct = round(sim_q, 4),
    Synthetic_Mean  = round(sim_mu, 4),
    Retain          = ifelse(retained_mask, "Yes", "No")
  )
  
  cat("--- Horn's Parallel Analysis (Exact Monte Carlo Simulation) ---\n")
  print(comparison_df, row.names = FALSE)
  cat(sprintf("\nHorn's Parallel Analysis suggests retaining %d component(s).\n\n", n_retained))
  
  # 5. Scree & Parallel Analysis Plot
  if (plot_comparison) {
    comp_indices <- seq_len(n_vars)
    y_limits <- range(c(emp_ev, sim_q, 1.0))
    
    plot(comp_indices, emp_ev, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
         ylim = y_limits, xlab = "Component Number", ylab = "Eigenvalue",
         main = "Horn's Parallel Analysis vs. Scree Trace", xaxt = "n")
    axis(1, at = comp_indices)
    lines(comp_indices, sim_q, type = "b", pch = 15, col = "#d62728", lty = 2, lwd = 1.75)
    lines(comp_indices, sim_mu, type = "b", pch = 17, col = "#7f7f7f", lty = 3, lwd = 1.25)
    abline(h = 1.0, col = "black", lty = 3)
    grid()
    
    legend("topright",
           legend = c("Empirical Eigenvalues",
                      sprintf("Simulated %dth %%ile", as.integer(percentile)),
                      "Simulated Mean",
                      "Kaiser (lambda = 1.0)"),
           col = c("#1f77b4", "#d62728", "#7f7f7f", "black"),
           pch = c(19, 15, 17, NA),
           lty = c(1, 2, 3, 3),
           lwd = c(2, 1.75, 1.25, 1),
           bty = "n")
  }
  
  list(
    n_factors_parallel   = n_retained,
    empirical_evals      = emp_ev,
    synthetic_percentile = sim_q,
    synthetic_mean       = sim_mu,
    comparison_table     = comparison_df
  )
}

# Step 1: Scree Analysis via Pearson Correlation
scree_res <- compute_pca_scree_analysis(df_or_corr = items_clean, kaiser_cutoff = 1.0)

# Step 2: Exact Parallel Analysis matching Python
pa_res <- perform_horn_parallel_analysis_exact(
  df             = items_clean, 
  n_observations = nrow(items_clean), # 26
  n_simulations  = 1000, 
  percentile     = 95.0, 
  random_state   = 42
)

# =========================================================================
# 5. Principal Component Analysis with Direct Oblimin Rotation
# =========================================================================
MANUAL_N_COMPONENTS <- 3   # Controls PCA dimensions (NULL = Horn's Parallel Analysis)
MANUAL_N_CLUSTERS   <- 3   # Controls K-Means groups   (NULL = NbClust consensus)

perform_abels_pca_oblimin_exact <- function(df = NULL, 
                                            corr_matrix = NULL, 
                                            n_factors = 3) {
  # 1. Establish Correlation Matrix (Matching Python clean_df.corr())
  if (!is.null(corr_matrix)) {
    r_mat <- as.matrix(corr_matrix)
    cols <- colnames(corr_matrix)
  } else if (!is.null(df)) {
    clean_df <- na.omit(as.data.frame(lapply(df, as.numeric)))
    cols <- colnames(clean_df)
    r_mat <- cor(clean_df)
  } else {
    stop("Must provide either `df` or `corr_matrix`.")
  }
  
  # Clean column names (strip non-breaking spaces and trim whitespace)
  cols <- trimws(gsub("\u00a0", " ", cols))
  colnames(r_mat) <- cols
  rownames(r_mat) <- cols
  n_vars <- ncol(r_mat)
  
  # 2. Positive Semi-Definiteness Guard (threshold <= 1e-7)
  evals_check <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(evals_check <= 1e-7)) {
    r_mat <- as.matrix(nearPD(r_mat, corr = TRUE, keepDiag = TRUE)$mat)
    dimnames(r_mat) <- list(cols, cols)
  }
  
  # 3. Eigendecomposition (Sorted Descending)
  decomp <- eigen(r_mat, symmetric = TRUE)
  sorted_evals <- decomp$values
  sorted_evecs <- decomp$vectors
  
  # 4. Unrotated Component Loadings: Lambda = V * sqrt(D)
  top_evals <- sorted_evals[seq_len(n_factors)]
  top_evecs <- sorted_evecs[, seq_len(n_factors), drop = FALSE]
  unrotated_loadings <- top_evecs %*% diag(sqrt(pmax(top_evals, 1e-12)), nrow = n_factors)
  
  comp_names <- paste0("Component ", seq_len(n_factors))
  rownames(unrotated_loadings) <- cols
  colnames(unrotated_loadings) <- comp_names
  
  # 5. Direct Oblimin Multi-Start (Quartimin, normalize = FALSE)
  sign_combinations <- list(
    c( 1,  1,  1),
    c(-1,  1,  1),
    c( 1, -1,  1),
    c( 1,  1, -1),
    c(-1, -1,  1),
    c(-1,  1, -1),
    c( 1, -1, -1),
    c(-1, -1, -1)
  )
  
  best_rot <- NULL
  best_loss <- Inf
  
  for (signs in sign_combinations) {
    L_start <- unrotated_loadings %*% diag(signs, nrow = n_factors)
    rot_try <- tryCatch(
      oblimin(L_start, gam = 0, normalize = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(rot_try)) {
      L_rot <- unclass(rot_try$loadings)
      loss <- 0
      for (j in 1:(n_factors - 1)) {
        for (k in (j + 1):n_factors) {
          loss <- loss + sum((L_rot[, j]^2) * (L_rot[, k]^2))
        }
      }
      if (loss < best_loss) {
        best_loss <- loss
        best_rot <- rot_try
      }
    }
  }
  
  # Extract and explicitly restore rownames stripped by GPArotation
  rotated_loadings <- unclass(best_rot$loadings)
  rownames(rotated_loadings) <- cols
  colnames(rotated_loadings) <- comp_names
  phi_matrix <- best_rot$Phi
  
  # 6. Robust Anchor Matching via Regex
  find_idx <- function(pattern, target_names) {
    m <- grep(pattern, target_names, ignore.case = TRUE)
    if (length(m) > 0) m[1] else NA_integer_
  }
  
  auth_row <- find_idx("Authentic", cols)
  law_row  <- find_idx("Lawful", cols)
  rel_row  <- find_idx("Reliab", cols)
  
  # Align columns if anchor items exist; otherwise retain original order
  if (!is.na(auth_row) && !is.na(law_row)) {
    c1_idx <- which.max(abs(rotated_loadings[auth_row, ]))
    c2_idx <- which.max(abs(rotated_loadings[law_row, ]))
    c3_idx <- setdiff(seq_len(3), c(c1_idx, c2_idx))
    
    col_order <- c(c1_idx, c2_idx, c3_idx)
    rotated_loadings <- rotated_loadings[, col_order]
    phi_matrix <- phi_matrix[col_order, col_order]
    
    # Align signs to match Python convention (C1 negative, C2 positive, C3 positive)
    target_signs <- c(
      if (rotated_loadings[auth_row, 1] > 0) -1 else 1,
      if (rotated_loadings[law_row, 2] < 0) -1 else 1,
      if (!is.na(rel_row) && rotated_loadings[rel_row, 3] < 0) -1 else 1
    )
    
    rotated_loadings <- rotated_loadings %*% diag(target_signs)
    phi_matrix <- diag(target_signs) %*% phi_matrix %*% diag(target_signs)
  }
  
  dimnames(rotated_loadings) <- list(cols, comp_names)
  dimnames(phi_matrix) <- list(comp_names, comp_names)
  
  # Pattern & Structure Matrices
  pattern_loadings <- as.data.frame(rotated_loadings)
  structure_matrix <- as.data.frame(as.matrix(pattern_loadings) %*% phi_matrix)
  colnames(structure_matrix) <- comp_names
  rownames(structure_matrix) <- cols
  factor_correlations <- as.data.frame(phi_matrix)
  
  # Item Clustering
  abs_pattern <- abs(as.matrix(pattern_loadings))
  primary_comp_idx <- max.col(abs_pattern)
  primary_assignments <- comp_names[primary_comp_idx]
  
  assignments_df <- data.frame(
    Item              = cols,
    Primary_Component = primary_assignments,
    Primary_Loading   = round(as.matrix(pattern_loadings)[cbind(seq_along(cols), primary_comp_idx)], 3),
    stringsAsFactors  = FALSE
  )
  
  list(
    eigenvalues         = sorted_evals,
    pattern_loadings    = pattern_loadings,
    structure_matrix    = structure_matrix,
    factor_correlations = factor_correlations,
    assignments_table   = assignments_df
  )
}

# --- Execution ---
# Pass raw `items` directly (or Pearson correlation matrix)
pca_exact <- perform_abels_pca_oblimin_exact(df = items_clean, n_factors = 3)

cat("--- Rotated Pattern Matrix (Table 3 Style) ---\n")
print(round(pca_exact$pattern_loadings, 3))

cat("\n--- Structure Matrix (Table 4 Style) ---\n")
print(round(pca_exact$structure_matrix, 3))

cat("\n--- Component Inter-Correlation Matrix (Phi) ---\n")
print(round(pca_exact$factor_correlations, 3))

cat("\n--- Item-to-Component Primary Assignments ---\n")
print(pca_exact$assignments_table, row.names = FALSE)


# =========================================================================
# 6. Psychometric Reliability Analysis (Alpha, Omega, Item-Dropped)
# =========================================================================
cat("===================================================================\n")
cat("PSYCHOMETRIC RELIABILITY ANALYSIS\n")
cat("===================================================================\n\n")

# A. Standard Cronbach's Alpha (Raw Data)
std_alpha <- psych::alpha(items_clean, check.keys = FALSE)

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
  Item             = colnames(items_clean),
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
  Item                     = colnames(items_clean),
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

# =========================================================================
# 1. Project Input Cases into R^3 Using Previous PCA Output
# =========================================================================
project_pca_scores <- function(df, pca_res) {
  # Clean and standardize raw observations
  clean_df <- na.omit(as.data.frame(lapply(df, as.numeric)))
  cols <- colnames(clean_df)
  cols_clean <- trimws(gsub("\u00a0", " ", cols))
  colnames(clean_df) <- cols_clean
  Z <- scale(clean_df)
  
  # Extract correlation matrix from standardized observations
  r_mat <- cor(Z)
  
  # Extract pattern loadings and factor correlations (Phi)
  pattern_mat <- as.matrix(pca_res$pattern_loadings)
  phi_mat     <- as.matrix(pca_res$factor_correlations)
  
  # Structure matrix: S = Pattern %*% Phi
  structure_mat <- pattern_mat %*% phi_mat
  
  # Standard regression score weights: W = R^(-1) %*% S
  W <- solve(r_mat) %*% structure_mat
  
  # Project observations into R^3: Scores = Z %*% W
  scores_r3 <- Z %*% W
  colnames(scores_r3) <- colnames(pattern_mat)
  rownames(scores_r3) <- paste0("Obs_", seq_len(nrow(scores_r3)))
  
  scores_r3[, 1] <- -1 * scores_r3[, 1]
  
  return(scores_r3)
}

# =========================================================================
# 2. Triangulation Engine in R^3 Coordinate Space
# =========================================================================
triangulate_r3_clusters <- function(scores_r3, k_range = 2:5, seed = 42) {
  set.seed(seed)
  n_obs <- nrow(scores_r3)
  k_candidates <- k_range[k_range < n_obs]
  elbow_k_list <- 1:max(k_candidates)
  
  # Euclidean distance in R^3
  dist_r3 <- dist(scores_r3, method = "euclidean")
  
  # A. Elbow Method (WCSS)
  wss <- sapply(elbow_k_list, function(k) {
    kmeans(scores_r3, centers = k, nstart = 25)$tot.withinss
  })
  
  # Knee detection via chord line
  x_pts <- elbow_k_list
  y_pts <- wss
  chord_vec  <- c(x_pts[length(x_pts)] - x_pts[1], y_pts[length(y_pts)] - y_pts[1])
  chord_unit <- chord_vec / sqrt(sum(chord_vec^2))
  
  diff_pts   <- cbind(x_pts - x_pts[1], y_pts - y_pts[1])
  proj_pts   <- (diff_pts %*% chord_unit) %*% t(chord_unit)
  dist_chord <- sqrt(rowSums((diff_pts - proj_pts)^2))
  opt_k_elbow <- x_pts[which.max(dist_chord)]
  
  # B. Silhouette Analysis (Average Linkage in R^3)
  sil_scores <- sapply(k_candidates, function(k) {
    clusters <- cutree(hclust(dist_r3, method = "average"), k = k)
    ss <- silhouette(clusters, dist_r3)
    mean(ss[, "sil_width"])
  })
  opt_k_sil <- k_candidates[which.max(sil_scores)]
  
  # C. Tibshirani's Gap Statistic (spaceH0 = 'original' matching Python)
  gap_res <- clusGap(
    scores_r3,
    FUNcluster = function(x, k) kmeans(x, centers = k, nstart = 25),
    K.max      = max(k_candidates),
    B          = 500,
    spaceH0    = "original"
  )
  gap_tab  <- gap_res$Tab
  gap_vals <- gap_tab[k_candidates, "gap"]
  gap_se   <- gap_tab[k_candidates, "SE.sim"]
  opt_k_gap <- maxSE(gap_tab[, "gap"], gap_tab[, "SE.sim"], method = "firstSEmax")
  
  # D. Summary Table
  triangulation_summary <- data.frame(
    Clusters       = paste0("k = ", k_candidates),
    WSS_Inertia    = round(wss[k_candidates], 3),
    Avg_Silhouette = round(sil_scores, 3),
    Gap_Statistic  = round(gap_vals, 3),
    Gap_Std_Error  = round(gap_se, 3)
  )
  
  cat("--- Triangulation Summary Table (Operating in R^3 PCA Space) ---\n")
  print(triangulation_summary, row.names = FALSE)
  cat(sprintf("\nTriangulated Decisions: Elbow k = %d | Silhouette k = %d | Gap (1-SE) k = %d\n\n",
              opt_k_elbow, opt_k_sil, opt_k_gap))
  
  # E. 4-Panel Diagnostic Dashboard
  ward_fit <- hclust(dist_r3, method = "ward.D2")
  
  par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3, 1))
  
  plot(elbow_k_list, wss, type = "b", pch = 19, col = "#8c564b", lwd = 2,
       xlab = "Number of Clusters (k)", ylab = "Total WSS (Inertia)",
       main = "Elbow Method (R^3)")
  abline(v = opt_k_elbow, col = "firebrick", lty = 2, lwd = 1.5)
  grid()
  
  plot(k_candidates, sil_scores, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
       xlab = "Number of Clusters (k)", ylab = "Average Silhouette Width",
       main = "Silhouette (R^3)")
  abline(v = opt_k_sil, col = "firebrick", lty = 2, lwd = 1.5)
  grid()
  
  plot(k_candidates, gap_vals, type = "b", pch = 15, col = "#2ca02c", lwd = 2,
       ylim = range(c(gap_vals - gap_se, gap_vals + gap_se)),
       xlab = "Number of Clusters (k)", ylab = "Gap(k)",
       main = "Gap Statistic (B=500, R^3)")
  arrows(k_candidates, gap_vals - gap_se, k_candidates, gap_vals + gap_se,
         length = 0.05, angle = 90, code = 3, col = "#2ca02c")
  abline(v = opt_k_gap, col = "firebrick", lty = 2, lwd = 1.5)
  grid()
  
  plot(as.dendrogram(ward_fit), horiz = TRUE,
       main = "Ward's Tree (Obs in R^3)",
       xlab = "Ward Dissimilarity")
  
  par(mfrow = c(1, 1))
  
  list(
    scores = scores_r3,
    summary = triangulation_summary,
    optimal_k = list(elbow = opt_k_elbow, silhouette = opt_k_sil, gap = opt_k_gap)
  )
}

# =========================================================================
# Execution Workflow
# =========================================================================
# Step 1: Reuse the existing PCA output
pca_res <- perform_abels_pca_oblimin_exact(df = items_clean, n_factors = 3)

# Step 2: Project observations into R^3
scores_r3 <- project_pca_scores(items_clean, pca_res)

# Step 3: Run triangulation directly on the 27x3 matrix
triangulation_out <- triangulate_r3_clusters(scores_r3, k_range = 2:5)

# =========================================================================
# 8. NbClust Multi-Criteria Cluster Validation (30 Indices Voting)
# =========================================================================
cat("===================================================================\n")
cat("NBCLUST CONSENSUS CLUSTER VALIDATION (30 INDICES)\n")
cat("===================================================================\n\n")

# 1. Ensure input is the 26x3 observation score matrix from project_pca_scores()
# Dim: 26 observations x 3 PCA component dimensions
obs_coords <- as.matrix(scores_r3)

# 2. Run NbClust across observations
# Set index = "alllong" or "all"
# Use method = "kmeans" (or "ward.D2" if prioritizing hierarchical linkage)
nb_res <- suppressWarnings(
  NbClust(
    data     = obs_coords,
    distance = "euclidean",
    min.nc   = 2,
    max.nc   = 5,
    method   = "kmeans",
    index    = "all"
  )
)

# 3. Extract the voting distribution
# Row 1 of Best.nc contains the recommended number of clusters for each index
votes_raw   <- as.numeric(nb_res$Best.nc[1, ])
valid_votes <- votes_raw[!is.na(votes_raw) & is.finite(votes_raw) & votes_raw >= 2]

votes_table <- as.data.frame(table(Clusters_k = valid_votes))
colnames(votes_table) <- c("Clusters_k", "Metric_Votes")

cat("--- Majority Voting Consensus Across All Computed Indices ---\n")
print(votes_table, row.names = FALSE)

# Optimal k by majority vote
opt_k_nbclust <- as.numeric(as.character(votes_table$Clusters_k[which.max(votes_table$Metric_Votes)]))
cat(sprintf("\nNbClust Consensus Optimal Clusters: k = %d\n\n", opt_k_nbclust))

# 4. Observation Partition Assignments
cat("--- Publication / Case Partition Assignments (Optimal k) ---\n")
partition_df <- data.frame(
  Publication = rownames(obs_coords),
  Cluster     = as.numeric(nb_res$Best.partition),
  stringsAsFactors = FALSE
)
print(partition_df)
cat("\n")

# =========================================================================
# 9. Project Publications, Cluster, and Decision Boundary Visualization
# =========================================================================
cat("===================================================================\n")
cat("PUBLICATION PROJECTION & CLUSTERING VISUALIZATION\n")
cat("===================================================================\n\n")

# =========================================================================
# 1. Coordinate Extraction & Semantic Label Alignment
# =========================================================================
# proj_scores is derived from project_pca_scores()
proj_scores  <- as.matrix(scores_r3)
proj_scores[, "Data Provenance"] <- -1 * proj_scores[, "Data Provenance"]
n_obs        <- nrow(proj_scores)
n_components <- ncol(proj_scores)

# Align semantic dimensions to match rotated pattern loading domains:
# Component 1 -> Data Provenance (Chain of Custody, Authenticity, Integrity)
# Component 2 -> Authority & Expertise (Lawful Acquisition, Competence)
# Component 3 -> Process Replication (Reliability, Repeatability)
comp_labels <- c("Data Provenance", "Authority & Expertise", "Process Replication")
colnames(proj_scores) <- comp_labels[seq_len(n_components)]

# =========================================================================
# 2. Strict Observation Alignment (Guarantees Matching Row Dimensions)
# =========================================================================
# Ensure item dataframe matches the exact N rows present in proj_scores
if (!is.null(rownames(proj_scores)) && all(rownames(proj_scores) %in% rownames(items_clean))) {
  items_matched <- items_clean[rownames(proj_scores), , drop = FALSE]
} else if (nrow(items_clean) == n_obs) {
  items_matched <- items_clean
} else {
  # Fallback subsetting if names do not match but row count was filtered
  items_matched <- head(items_clean, n_obs)
}

# =========================================================================
# 3. K-Means Clustering on Component Space Coordinates
# =========================================================================
n_clusters <- if (exists("MANUAL_N_CLUSTERS") && !is.null(MANUAL_N_CLUSTERS)) {
  MANUAL_N_CLUSTERS
} else if (exists("opt_k_nbclust") && !is.na(opt_k_nbclust)) {
  opt_k_nbclust
} else {
  2
}

set.seed(42)
km_fit <- kmeans(proj_scores, centers = n_clusters, nstart = 25)

cluster_names <- if (n_clusters == 2) {
  c("Process-Driven", "Authority-Driven")
} else {
  paste("Cluster", seq_len(n_clusters))
}

assigned_clusters <- factor(cluster_names[km_fit$cluster], levels = cluster_names)

# =========================================================================
# 4. Metric Profiling & Centroid Summaries
# =========================================================================
# Centroids in R^3 component score space
cluster_centroids <- aggregate(
  proj_scores,
  by  = list(Assigned_Cluster = assigned_clusters),
  FUN = mean
)
cluster_centroids[, -1] <- round(cluster_centroids[, -1], 3)

# Original Likert item means by cluster
metric_profiles <- aggregate(
  as.data.frame(items_matched),
  by  = list(Assigned_Cluster = assigned_clusters),
  FUN = mean
)
metric_profiles[, -1] <- round(metric_profiles[, -1], 3)

cat("--- Cluster Centroids (Component Space) ---\n")
print(cluster_centroids)
cat("\n")

cat("--- Original Item Profiles by Cluster ---\n")
print(metric_profiles)
cat("\n")

# =========================================================================
# 5. Publication-Grade Decision Boundary Visualization
# =========================================================================
palette   <- c("#FF2A54", "#0066FF", "#9B59B6")
library(scatterplot3d)
library(e1071)

# =========================================================================
# 1. Component Extraction and PC Naming
# =========================================================================
proj_scores  <- as.matrix(scores_r3)
n_obs        <- nrow(proj_scores)
n_components <- ncol(proj_scores)

# Principal component dimension names
comp_labels <- c("Data Provenance", "Authority & Expertise", "Process Replication")
colnames(proj_scores) <- comp_labels[seq_len(n_components)]

# Ensure item rows align with score projections
if (!is.null(rownames(proj_scores)) && all(rownames(proj_scores) %in% rownames(items))) {
  items_matched <- items[rownames(proj_scores), , drop = FALSE]
} else if (nrow(items) == n_obs) {
  items_matched <- items
} else {
  items_matched <- head(items, n_obs)
}

# =========================================================================
# 2. K-Means Clustering Mapped to Principal Component Names
# =========================================================================

# Create a 2-column canvas: 80% plot, 20% legend
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1.2))

# --- Panel 1: 3D Scatterplot ---
par(mar = c(5, 4, 4, 1))

s3d <- scatterplot3d(
  proj_scores[, 1], proj_scores[, 2], proj_scores[, 3],
  color       = pt_colors,
  pch         = 19,
  cex.symbols = 1.3,
  xlab        = colnames(proj_scores)[1],
  ylab        = colnames(proj_scores)[2],
  zlab        = colnames(proj_scores)[3],
  main        = "Projection of Forensics Publications onto Principal Components",
  angle       = 60,
  grid        = TRUE,
  box         = TRUE
)

# Add centroids
centroids_3d <- km_fit$centers
s3d$points3d(
  centroids_3d[, 1], centroids_3d[, 2], centroids_3d[, 3],
  col = palette[seq_len(n_clusters)],
  pch = 8,
  cex = 2,
  lwd = 2
)

# --- Panel 2: Dedicated Legend Canvas ---
par(mar = c(5, 0, 4, 1))
plot.new() # Blank canvas

legend(
  "left",
  legend  = c("Mechanism-Driven", "Admissibility-Driven", "Provenance-Driven"),
  col     = palette[seq_along(levels(assigned_clusters))],
  pch     = 19,
  title   = "Cluster Category",
  bty     = "o",
  box.col = "#444444",
  box.lwd = 1.1,
  bg      = "white",
  cex     = 0.9
)

# Reset layout to default
layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))
# -------------------------------------------------------------------------
# 10. Tally Cluster Counts and Compute Donut Chart Geometry (3 Groups)
# -------------------------------------------------------------------------
set.seed(42)
km_fit_3 <- kmeans(proj_scores, centers = 3, nstart = 50)

# Check cluster sizes to confirm all 3 groups are populated
print(table(km_fit_3$cluster))

# Assign semantic cluster names
# Cluster 1, 2, 3 mapped to your 3 analytical profiles
cluster_names <- c("Mechanism-Driven", "Admissibility-Driven", "Provenance-Driven")
assigned_clusters_3 <- factor(cluster_names[km_fit_3$cluster], levels = cluster_names)

# Update df_proj with the 3-cluster partition
df_proj <- data.frame(
  Publication      = rownames(proj_scores),
  Assigned_Cluster = assigned_clusters_3,
  round(proj_scores, 3),
  stringsAsFactors = FALSE
)

# =========================================================================
# 2. Tally Counts & Calculate Polar Midpoints
# =========================================================================
cluster_counts <- df_proj %>%
  count(Assigned_Cluster, name = "count") %>%
  mutate(
    pct   = count / sum(count) * 100,
    label = sprintf("%d\n(%.1f%%)", count, pct)
  ) %>%
  # Arrange descending so the polar coordinates match ggplot's stack order
  arrange(desc(Assigned_Cluster)) %>%
  mutate(
    y_pos = cumsum(count) - (count / 2)
  )

print(cluster_counts)

# =========================================================================
# 3. Render 3-Cluster Donut Chart
# =========================================================================
# 1. Enforce factor levels so colors, slices, and legend order align exactly
cluster_order <- c("Mechanism-Driven", "Admissibility-Driven", "Provenance-Driven")

cluster_counts <- cluster_counts %>%
  mutate(Assigned_Cluster = factor(Assigned_Cluster, levels = cluster_order)) %>%
  arrange(desc(Assigned_Cluster)) %>%
  mutate(y_pos = cumsum(count) - (count / 2))

# 2. Color mapping matching the 3D scatterplot palette
color_map <- c(
  "Admissibility-Driven" = "#0066FF",  # High Provenance, High Authority
  "Provenance-Driven"    = "#9B59B6",  # High Provenance, Low Authority
  "Mechanism-Driven"     = "#FF2A54"   # Low Provenance, Low Authority
)

# 3. Render Donut Plot with formatted title and matching colors
ggplot(cluster_counts, aes(x = 2, y = count, fill = Assigned_Cluster)) +
  geom_col(color = "white", linewidth = 1.2, width = 0.8) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = color_map) +
  geom_text(
    aes(y = y_pos, label = label),
    color      = "white",
    size       = 4.5,
    fontface   = "bold",
    lineheight = 0.9
  ) +
  labs(
    title    = "Publication Cluster Distribution",
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 4)),
    plot.subtitle   = element_text(hjust = 0.5, size = 11, color = "#555555", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 11, face = "bold"),
    plot.margin     = margin(12, 12, 12, 12)
  )

# -------------------------------------------------------------------------
# 13. Prepare Eigenvalue Data (First 3 Components)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 13. Prepare Eigenvalue Data (First 3 Components from Horn's Analysis)
# -------------------------------------------------------------------------
# Pull empirical eigenvalues directly from earlier Horn's parallel analysis
# (Supports pa_res$empirical_evals from perform_horn_parallel_analysis_exact 
#  or pa_out$pc.values / emp_ev from psych::fa.parallel)
empirical_values <- if (exists("pa_res") && !is.null(pa_res$empirical_evals)) {
  pa_res$empirical_evals[1:3]
} else if (exists("emp_ev")) {
  emp_ev[1:3]
} else if (exists("pa_out") && !is.null(pa_out$pc.values)) {
  pa_out$pc.values[1:3]
} else {
  stop("Empirical eigenvalues not found. Please run Horn's Parallel Analysis first.")
}

# Align component names and colors strictly in index order:
# PC1 = Data Provenance, PC2 = Authority & Expertise, PC3 = Process Replication
comp_names <- c("Data Provenance", "Authority & Expertise", "Process Replication")
comp_colors <- c(
  "Data Provenance"       = "#0066FF",
  "Authority & Expertise" = "#9B59B6",
  "Process Replication"   = "#FF2A54"
)

df_eigen <- data.frame(
  Component  = factor(comp_names, levels = comp_names),
  Eigenvalue = empirical_values,
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# 14. Render Publication-Grade Bar Plot
# -------------------------------------------------------------------------
ggplot(df_eigen, aes(x = Component, y = Eigenvalue, fill = Component)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  # Kaiser-Guttman threshold reference line
  geom_hline(
    yintercept = 1.0, 
    linetype = "dashed", 
    color = "#555555", 
    linewidth = 0.8
  ) +
  annotate(
    "text", 
    x = 0.7, 
    y = 1.05, 
    label = "Kaiser Cutoff (1.0)", 
    size = 3.5, 
    fontface = "italic", 
    color = "#555555"
  ) +
  scale_fill_manual(values = comp_colors) +
  geom_text(
    aes(label = sprintf("%.3f", Eigenvalue)),
    vjust    = -0.6,
    size     = 4.2,
    fontface = "bold"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = NULL,
    y = "Empirical Eigenvalue",
    title = "Empirical Factor Retention",
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle= element_text(size = 10.5, color = "#555555", hjust = 0.5, margin = margin(b = 12)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text.x  = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10, color = "black"),
    plot.margin  = margin(15, 15, 15, 15)
  )

