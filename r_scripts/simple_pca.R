source("helpers.R")

# Vectorized Kendall correlation with conditional PSD projection
compute_psd_kendall_simple <- function(df, use = "pairwise.complete.obs") {
  # Coerce to numeric matrix in a single step
  mat <- data.matrix(df)
  
  # Compute Kendall's tau-b across all 9 variables in compiled C
  tau_mat <- stats::cor(mat, method = "kendall", use = use)
  
  # Check if indefinite (negative eigenvalues)
  min_eigen <- min(eigen(tau_mat, symmetric = TRUE, only.values = TRUE)$values)
  
  # Project only if necessary
  if (min_eigen < 0) {
    pd_obj  <- Matrix::nearPD(tau_mat, corr = TRUE, keepDiag = TRUE)
    tau_mat <- as.matrix(pd_obj$mat)
    dimnames(tau_mat) <- list(colnames(df), colnames(df))
  }
  
  tau_mat
}

assess_factorability_simple <- function(cor_mat, n_obs, verbose = TRUE) {
  # 1. Bartlett's Test of Sphericity
  bartlett_res <- psych::cortest.bartlett(cor_mat, n = n_obs)
  
  # 2. Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy
  kmo_res <- psych::KMO(cor_mat)
  
  # Format per-item table cleanly
  item_msa <- data.frame(
    Variable = names(kmo_res$MSAi),
    MSA      = round(unname(kmo_res$MSAi), 3),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("--- Bartlett's Test of Sphericity ---\n")
    cat(sprintf("Chi-Square: %.2f | df: %d | p-value: %.3e\n\n", 
                bartlett_res$chisq, bartlett_res$df, bartlett_res$p.value))
    
    cat("--- KMO Sampling Adequacy ---\n")
    cat(sprintf("Overall MSA: %.3f\n\n", kmo_res$MSA))
    
    cat("Per-Item MSA:\n")
    print(item_msa, row.names = FALSE)
  }
  
  # Return structured results for paper tables or markdown reporting
  invisible(list(
    bartlett    = bartlett_res,
    kmo_overall = kmo_res$MSA,
    kmo_items   = item_msa
  ))
}

compute_pca_scree_analysis_simple <- function(x, kaiser_cutoff = 1.0, verbose = TRUE) {
  # 1. Obtain correlation matrix cleanly
  if (is.matrix(x) && nrow(x) == ncol(x) && all(abs(diag(x) - 1) < 1e-5)) {
    r_mat <- x
  } else {
    r_mat <- stats::cor(data.matrix(x), use = "pairwise.complete.obs")
  }
  
  # 2. Spectral decomposition (LAPACK returns eigenvalues pre-sorted descending)
  ev <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values
  
  # Guard against negative eigenvalues from rank correlations
  ev_clean <- pmax(ev, 0)
  
  pct_var  <- (ev_clean / sum(ev_clean)) * 100
  cum_var  <- cumsum(pct_var)
  n_kaiser <- sum(ev_clean >= kaiser_cutoff)
  
  # 3. Formatted variance summary table
  summary_df <- data.frame(
    Component     = paste0("PC", seq_along(ev_clean)),
    Eigenvalue    = round(ev_clean, 4),
    Pct_Variance  = round(pct_var, 2),
    Cum_Variance  = round(cum_var, 2),
    Retain_Kaiser = ifelse(ev_clean >= kaiser_cutoff, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("--- Eigenvalue Decomposition & Variance Explained ---\n")
    print(summary_df, row.names = FALSE)
    
    retained_pct <- if (n_kaiser > 0) cum_var[n_kaiser] else 0.0
    cat(sprintf("\nKaiser-Guttman (eigenvalue >= %.2f) suggests retaining %d component(s) (%.2f%% cumulative variance).\n\n",
                kaiser_cutoff, n_kaiser, retained_pct))
  }
  
  list(
    cor_matrix            = r_mat,
    eigenvalues           = ev_clean,
    pct_variance          = pct_var,
    cum_variance          = cum_var,
    n_factors_kaiser      = n_kaiser,
    variance_table        = summary_df,
    total_variance_kaiser = if (n_kaiser > 0) cum_var[n_kaiser] else 0.0
  )
}

perform_horn_parallel_analysis_simple <- function(x,
                                                  n_obs = NULL,
                                                  n_simulations = 1000,
                                                  percentile = 95.0,
                                                  kaiser_cutoff = 1.0,
                                                  plot_comparison = TRUE,
                                                  verbose = TRUE) {
  # 1. Determine correlation matrix and sample size (N=90)
  if (is.matrix(x) && nrow(x) == ncol(x) && all(abs(diag(x) - 1) < 1e-5)) {
    r_mat <- x
    if (is.null(n_obs)) {
      stop("When providing a correlation matrix, `n_obs` (e.g., 90) must be specified.")
    }
  } else {
    mat <- data.matrix(x)
    n_obs <- nrow(stats::na.omit(mat))
    r_mat <- stats::cor(mat, use = "pairwise.complete.obs")
  }
  
  n_vars <- ncol(r_mat)
  
  # 2. Empirical eigenvalues (LAPACK returns sorted descending)
  emp_ev <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values
  emp_ev <- pmax(emp_ev, 0)
  pct_var <- (emp_ev / sum(emp_ev)) * 100
  cum_var <- cumsum(pct_var)
  
  # 3. Monte Carlo Simulation of Uncorrelated Gaussian Noise
  # Matrix pre-allocation: 1000 x 9
  synthetic_evals <- matrix(0, nrow = n_simulations, ncol = n_vars)
  for (i in seq_len(n_simulations)) {
    noise <- matrix(stats::rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
    # LAPACK guarantees eigenvalues are descending; no sort() needed
    synthetic_evals[i, ] <- eigen(stats::cor(noise), symmetric = TRUE, only.values = TRUE)$values
  }
  
  sim_mu <- colMeans(synthetic_evals)
  sim_q  <- apply(synthetic_evals, 2, stats::quantile, probs = percentile / 100, type = 7)
  
  retained_mask <- emp_ev > sim_q
  n_retained    <- sum(retained_mask)
  n_kaiser      <- sum(emp_ev >= kaiser_cutoff)
  
  # 4. Summary Table for Reporting
  comparison_df <- data.frame(
    Component       = paste0("PC", seq_len(n_vars)),
    Empirical_EV    = round(emp_ev, 4),
    Pct_Variance    = round(pct_var, 2),
    Cum_Variance    = round(cum_var, 2),
    Synthetic_95Pct = round(sim_q, 4),
    Synthetic_Mean  = round(sim_mu, 4),
    Retain_Parallel = ifelse(retained_mask, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("--- Horn's Parallel Analysis (Monte Carlo Simulation) ---\n")
    cat(sprintf("Sample Size (N): %d | Variables (P): %d | Simulations: %d\n\n", n_obs, n_vars, n_simulations))
    print(comparison_df, row.names = FALSE)
    
    retained_cum_var <- if (n_retained > 0) cum_var[n_retained] else 0.0
    cat(sprintf("\nParallel Analysis suggests retaining %d component(s) (%.2f%% cumulative variance).\n", 
                n_retained, retained_cum_var))
    cat(sprintf("Kaiser-Guttman rule suggests retaining %d component(s).\n\n", n_kaiser))
  }
  
  # 5. Scree / Parallel Analysis Plot
  if (plot_comparison) {
    comp_idx <- seq_len(n_vars)
    y_limits <- range(c(emp_ev, sim_q, kaiser_cutoff))
    
    plot(comp_idx, emp_ev, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
         ylim = y_limits, xlab = "Component Number", ylab = "Eigenvalue",
         main = "Parallel Analysis vs. Empirical Scree Plot", xaxt = "n")
    axis(1, at = comp_idx)
    lines(comp_idx, sim_q, type = "b", pch = 15, col = "#d62728", lty = 2, lwd = 1.75)
    lines(comp_idx, sim_mu, type = "b", pch = 17, col = "#7f7f7f", lty = 3, lwd = 1.25)
    abline(h = kaiser_cutoff, col = "black", lty = 3)
    grid()
    
    legend("topright",
           legend = c("Empirical Eigenvalues",
                      sprintf("Simulated %dth %%ile", as.integer(percentile)),
                      "Simulated Mean",
                      sprintf("Kaiser Rule (lambda = %.1f)", kaiser_cutoff)),
           col = c("#1f77b4", "#d62728", "#7f7f7f", "black"),
           pch = c(19, 15, 17, NA),
           lty = c(1, 2, 3, 3),
           lwd = c(2, 1.75, 1.25, 1),
           bty = "n")
  }
  
  invisible(list(
    n_factors_parallel   = n_retained,
    n_factors_kaiser     = n_kaiser,
    empirical_evals      = emp_ev,
    pct_variance         = pct_var,
    cum_variance         = cum_var,
    synthetic_percentile = sim_q,
    synthetic_mean       = sim_mu,
    comparison_table     = comparison_df
  ))
}

perform_abels_pca_oblimin_simple <- function(x, n_factors, n_obs, align_anchors = TRUE) {
  # 1. Input checks & matrix conversion
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square correlation matrix.")
  }
  
  # 2. Extract full eigenvalues from the correlation matrix
  raw_evals <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  sorted_evals <- sort(raw_evals, decreasing = TRUE)
  
  # 3. PCA with Direct Oblimin rotation via psych::principal
  # rotate = "oblimin" performs direct oblimin (GPArotation)
  pca_fit <- psych::principal(
    r      = x,
    nfactors = n_factors,
    rotate = "oblimin",
    n.obs  = n_obs,
    covar  = FALSE
  )
  
  # Pattern matrix (P) and Component correlation matrix (Phi)
  pattern_loadings <- unclass(pca_fit$loadings)
  phi_matrix       <- pca_fit$Phi
  
  if (is.null(phi_matrix)) {
    # If n_factors == 1 or orthogonal fallback occurs
    phi_matrix <- diag(n_factors)
    colnames(phi_matrix) <- rownames(phi_matrix) <- colnames(pattern_loadings)
  }
  
  # 4. Optional anchor alignment (reflection of signs)
  # Ensures the dominant variable on each component has a positive loading
  if (align_anchors) {
    for (j in seq_len(n_factors)) {
      max_idx <- which.max(abs(pattern_loadings[, j]))
      if (pattern_loadings[max_idx, j] < 0) {
        pattern_loadings[, j] <- -pattern_loadings[, j]
        phi_matrix[j, ]       <- -phi_matrix[j, ]
        phi_matrix[, j]       <- -phi_matrix[, j]
      }
    }
  }
  
  # 5. Structure Matrix (S = P %*% Phi)
  # Relates variables to components accounting for factor correlations
  structure_matrix <- pattern_loadings %*% phi_matrix
  dimnames(structure_matrix) <- dimnames(pattern_loadings)
  
  # 6. Assignments table: Map each variable to its primary component
  var_names <- rownames(pattern_loadings)
  if (is.null(var_names)) {
    var_names <- paste0("Var_", seq_len(nrow(pattern_loadings)))
    rownames(pattern_loadings) <- var_names
  }
  
  comp_names <- colnames(pattern_loadings)
  max_comp_idx <- apply(abs(pattern_loadings), 1, which.max)
  
  assignments_df <- data.frame(
    Variable         = var_names,
    Assigned_Factor  = comp_names[max_comp_idx],
    Max_Loading      = round(mapply(function(r, c) pattern_loadings[r, c], 
                                    seq_along(max_comp_idx), 
                                    max_comp_idx), 4),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  
  # 7. Return formatted output list
  list(
    eigenvalues         = sorted_evals,
    n_factors           = n_factors,
    pattern_loadings    = round(as.data.frame(pattern_loadings), 4),
    structure_matrix    = round(as.data.frame(structure_matrix), 4),
    factor_correlations = round(as.data.frame(phi_matrix), 4),
    assignments_table   = assignments_df
  )
}

compute_internal_reliability_simple <- function(df, verbose = TRUE) {
  # 1. Clean raw data: drop missing cases and ensure numeric matrix/data.frame
  items_clean <- stats::na.omit(as.data.frame(df))
  n_obs       <- nrow(items_clean)
  n_items     <- ncol(items_clean)
  item_names  <- colnames(items_clean)
  
  # 2. Standard Cronbach's Alpha (Raw Data)
  std_alpha <- psych::alpha(items_clean, check.keys = FALSE, warnings = FALSE)
  
  # 3. Polychoric Correlation Matrix & Ordinal Alpha
  # Fast 2-step polychoric estimation for ordinal/Likert scales
  poly_obj       <- suppressWarnings(psych::polychoric(items_clean))
  cor_polychoric <- poly_obj$rho
  ord_alpha      <- psych::alpha(cor_polychoric, check.keys = FALSE, warnings = FALSE)
  
  # 4. McDonald's Omega (1-Factor MINRES Solution)
  om_fit <- suppressWarnings(
    psych::omega(cor_polychoric, nfactors = 1, fm = "minres", plot = FALSE)
  )
  
  # 5. Item-Dropped Omega Calculation
  omega_dropped <- sapply(seq_len(n_items), function(i) {
    sub_r   <- cor_polychoric[-i, -i]
    fit_sub <- suppressWarnings(
      psych::omega(sub_r, nfactors = 1, fm = "minres", plot = FALSE)
    )
    fit_sub$omega.tot
  })
  
  # 6. Structured Tables for Reporting
  summary_metrics <- data.frame(
    Metric = c("Standard Cronbach's Alpha (Raw)", 
               "Ordinal Cronbach's Alpha (Polychoric)", 
               "McDonald's Omega Total (omega_t)"),
    Value  = round(c(std_alpha$total$raw_alpha, 
                     ord_alpha$total$raw_alpha, 
                     om_fit$omega.tot), 3),
    stringsAsFactors = FALSE
  )
  
  loadings_table <- data.frame(
    Item             = item_names,
    Std_Loading_g    = round(as.vector(om_fit$schmid$sl[, "g"]), 3),
    Error_Variance_u = round(as.vector(om_fit$schmid$sl[, "u2"]), 3),
    stringsAsFactors = FALSE
  )
  
  diagnostics_table <- data.frame(
    Item                     = item_names,
    Std_Alpha_If_Dropped     = round(std_alpha$alpha.drop[, "raw_alpha"], 3),
    Ordinal_Alpha_If_Dropped = round(ord_alpha$alpha.drop[, "raw_alpha"], 3),
    Delta_Ordinal_Alpha      = round(ord_alpha$alpha.drop[, "raw_alpha"] - ord_alpha$total$raw_alpha, 3),
    Omega_If_Dropped         = round(omega_dropped, 3),
    Delta_Omega              = round(omega_dropped - om_fit$omega.tot, 3),
    stringsAsFactors         = FALSE
  )
  
  if (verbose) {
    cat("--- Scale Reliability Summary (N =", n_obs, ") ---\n")
    print(summary_metrics, row.names = FALSE)
    cat("\n--- Single-Factor Standardized Loadings (MINRES) ---\n")
    print(loadings_table, row.names = FALSE)
    cat("\n--- Item-Deleted Reliability Diagnostics ---\n")
    print(diagnostics_table, row.names = FALSE)
    cat("\n")
  }
  
  invisible(list(
    summary_metrics   = summary_metrics,
    loadings_table    = loadings_table,
    diagnostics_table = diagnostics_table,
    cor_polychoric    = cor_polychoric,
    raw_alpha_obj     = std_alpha,
    poly_alpha_obj    = ord_alpha,
    omega_obj         = om_fit
  ))
}

run_simple_pca <- function(items,
                           n_obs,
                           n_boot = 1000,
                           percentile = 95.0,
                           manual_n_components = 2,
                           ...) {
  # =========================================================================
  # 1. Polychoric Correlation Matrix (Exact 2-Step MLE)
  # =========================================================================
  
  kendall_matrix <- compute_psd_kendall_simple(items_clean)
  
  cat("--- Positive Semi-Definite Kendall Tau-b Correlation Matrix ---\n")
  print(round(kendall_matrix, 3))
  
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
    n_simulations   = n_boot,
    percentile      = percentile,
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
  
  info <- list(
    pca_res = pca_res,
    pa_res = pa_res
  )
  
  invisible(info)
}