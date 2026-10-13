# Dependencies:
# - readxl
# - psych
# - boot
# - here
# - cluster
# - NbClust
# - e1071
# - scatterplot3d
# - ggplot2
# - dplyr
# - Matrix
# - GPArotation
# - DescTools

source("helpers.R")

safe_kendall <- function(x, y) {
  # Return 0.0 if either vector has zero variance (all values identical)
  if (length(unique(x)) <= 1 || length(unique(y)) <= 1) {
    return(0.0)
  }
  
  val <- suppressWarnings(KendallTauB(x, y, conf.level = NA))
  if (is.na(val) || is.nan(val)) 0.0 else as.numeric(val)
}

# Conditional Higham projection only when non-PSD
project_near_psd <- function(mat, tol = 1e-8, corr = TRUE) {
  # Quick eigenvalue check on 9x9 matrix (near instantaneous)
  evals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  
  # If already positive semi-definite, return unchanged to preserve raw correlations
  if (min(evals) >= -tol) {
    return(mat)
  }
  
  # Project only if indefinite
  pd_result <- Matrix::nearPD(mat, corr = corr, keepDiag = TRUE)
  out <- as.matrix(pd_result$mat)
  dimnames(out) <- dimnames(mat)
  out
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
  
  # Apply Higham nearest correlation projection if non-PSD
  if (min_eigen < 0) {
    tau_psd <- project_near_psd(tau_mat)
  } else {
    tau_psd <- tau_mat
  }
  
  return(tau_psd)
}

assess_factorability <- function(cor_mat, n_obs) {

  cat("--- Bartlett's Test of Sphericity ---\n")
  bartlett_res <- cortest.bartlett(cor_mat, n = n_obs)
  cat(sprintf("Chi-Square: %.3f | df: %d | p-value: %.3e\n\n", 
              bartlett_res$chisq, bartlett_res$df, bartlett_res$p.value))

  cat("--- Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy ---\n")
  kmo_res <- KMO(cor_mat)

  cat(sprintf("Overall MSA: %.3f\n\n", kmo_res$MSA))

  kmo_table <- data.frame(
    Item = names(kmo_res$MSAi),
    MSA  = round(kmo_res$MSAi, 3),
    row.names = NULL
  )

  cat("Per-Item Measures of Sampling Adequacy (MSAi):\n")
  print(kmo_table)
  cat("\n")
}


# Core Scree & Eigenvalue Analysis Function
compute_pca_scree_analysis <- function(df_or_corr, kaiser_cutoff = 1.0, verbose = TRUE) {
  
  parsed <- resolve_corr_matrix(df_or_corr = df_or_corr)
  r_mat  <- ensure_psd(parsed$cor_matrix)
  
  # Eigendecomposition sorted descending
  decomp <- eigen(r_mat, symmetric = TRUE)
  ev <- sort(decomp$values, decreasing = TRUE)
  
  pct_var <- (ev / sum(ev)) * 100
  cum_var <- cumsum(pct_var)
  n_kaiser <- sum(ev >= kaiser_cutoff)
  
  summary_df <- data.frame(
    Component     = paste("Component", seq_along(ev)),
    Eigenvalue    = round(ev, 4),
    Pct_Variance  = round(pct_var, 2),
    Cum_Variance  = round(cum_var, 2),
    Retain_Kaiser = ifelse(ev >= kaiser_cutoff, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("--- Eigenvalue Decomposition & Kaiser-Guttman Rule ---\n")
    print(summary_df, row.names = FALSE)
    cat(sprintf("\nKaiser-Guttman suggests retaining %d component(s) (Cumulative Variance: %.2f%%).\n\n",
                n_kaiser, cum_var[n_kaiser]))
  }
  
  list(
    cor_matrix            = r_mat,
    eigenvalues           = ev,
    pct_variance          = pct_var,
    cum_variance          = cum_var,
    n_factors_kaiser      = n_kaiser,
    kaiser_cutoff         = kaiser_cutoff,
    variance_table        = summary_df,
    total_variance_kaiser = cum_var[n_kaiser]
  )
}

# Horn's Parallel Analysis with Integrated Variance Statistics
perform_horn_parallel_analysis_exact <- function(df = NULL,
                                                 corr_matrix = NULL,
                                                 df_or_corr = NULL,
                                                 n_observations = 26,
                                                 n_simulations = 1000,
                                                 percentile = 95.0,
                                                 kaiser_cutoff = 1.0,
                                                 plot_comparison = TRUE,
                                                 verbose = TRUE) {
  
  # Resolve input matrix / dataframe seamlessly
  parsed <- resolve_corr_matrix(df = df, corr_matrix = corr_matrix, df_or_corr = df_or_corr)
  n_obs  <- if (!is.na(parsed$n_obs)) parsed$n_obs else n_observations
  r_mat  <- ensure_psd(parsed$cor_matrix)
  
  # Delegate eigenvalue & variance decomposition to scree engine
  scree_res <- compute_pca_scree_analysis(r_mat, kaiser_cutoff = kaiser_cutoff, verbose = FALSE)
  emp_ev    <- scree_res$eigenvalues
  pct_var   <- scree_res$pct_variance
  cum_var   <- scree_res$cum_variance
  n_vars    <- length(emp_ev)
  
  # Sample size determination
  if (is.data.frame(r_mat) && !(nrow(r_mat) == ncol(r_mat) && all(diag(as.matrix(r_mat)) == 1))) {
    n_obs <- nrow(na.omit(r_mat))
  } else {
    n_obs <- n_observations
  }
  
  # Monte Carlo Simulation
  synthetic_evals <- matrix(0, nrow = n_simulations, ncol = n_vars)
  for (i in seq_len(n_simulations)) {
    noise <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
    synthetic_evals[i, ] <- sort(eigen(cor(noise), symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  }
  
  sim_mu <- colMeans(synthetic_evals)
  sim_q  <- apply(synthetic_evals, 2, function(col) quantile(col, probs = percentile / 100, type = 7))
  
  retained_mask <- emp_ev > sim_q
  n_retained    <- sum(retained_mask)
  
  # Merged comparison table
  comparison_df <- data.frame(
    Component       = paste("Component", seq_len(n_vars)),
    Empirical_EV    = round(emp_ev, 4),
    Pct_Variance    = round(pct_var, 2),
    Cum_Variance    = round(cum_var, 2),
    Synthetic_95Pct = round(sim_q, 4),
    Synthetic_Mean  = round(sim_mu, 4),
    Retain          = ifelse(retained_mask, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("--- Horn's Parallel Analysis (Exact Monte Carlo Simulation) ---\n")
    print(comparison_df, row.names = FALSE)
    cat(sprintf("\nHorn's Parallel Analysis suggests retaining %d component(s)", n_retained))
    if (n_retained > 0) {
      cat(sprintf(" (Cumulative Variance: %.2f%%).\n\n", cum_var[n_retained]))
    } else {
      cat(".\n\n")
    }
  }
  
  if (plot_comparison) {
    comp_indices <- seq_len(n_vars)
    y_limits     <- range(c(emp_ev, sim_q, kaiser_cutoff))
    
    plot(comp_indices, emp_ev, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
         ylim = y_limits, xlab = "Component Number", ylab = "Eigenvalue",
         main = "Horn's Parallel Analysis vs. Scree Trace", xaxt = "n")
    axis(1, at = comp_indices)
    lines(comp_indices, sim_q, type = "b", pch = 15, col = "#d62728", lty = 2, lwd = 1.75)
    lines(comp_indices, sim_mu, type = "b", pch = 17, col = "#7f7f7f", lty = 3, lwd = 1.25)
    abline(h = kaiser_cutoff, col = "black", lty = 3)
    grid()
    
    legend("topright",
           legend = c("Empirical Eigenvalues",
                      sprintf("Simulated %dth %%ile", as.integer(percentile)),
                      "Simulated Mean",
                      sprintf("Kaiser (lambda = %.1f)", kaiser_cutoff)),
           col = c("#1f77b4", "#d62728", "#7f7f7f", "black"),
           pch = c(19, 15, 17, NA),
           lty = c(1, 2, 3, 3),
           lwd = c(2, 1.75, 1.25, 1),
           bty = "n")
  }
  
  invisible(list(
    n_factors_parallel   = n_retained,
    n_factors_kaiser     = scree_res$n_factors_kaiser,
    empirical_evals      = emp_ev,
    pct_variance         = pct_var,
    cum_variance         = cum_var,
    synthetic_percentile = sim_q,
    synthetic_mean       = sim_mu,
    comparison_table     = comparison_df,
    scree_summary        = scree_res
  ))
}

perform_abels_pca_oblimin_exact <- function(df = NULL, 
                                            corr_matrix = NULL, 
                                            df_or_corr = NULL,
                                            n_factors = 3,
                                            n_observations = 26) {
  # Resolve Correlation Matrix & Clean Labels
  parsed <- resolve_corr_matrix(df = df, corr_matrix = corr_matrix, df_or_corr = df_or_corr)
  r_mat  <- ensure_psd(parsed$cor_matrix)
  cols   <- parsed$labels
  n_obs  <- if (!is.na(parsed$n_obs)) parsed$n_obs else n_observations
  
  # Dynamic Factor Retention (invokes Horn's Parallel Analysis if n_factors is NULL)
  if (is.null(n_factors)) {
    pa_res    <- perform_horn_parallel_analysis_exact(r_mat, n_observations = n_obs, plot_comparison = FALSE, verbose = FALSE)
    n_factors <- max(1, pa_res$n_factors_parallel)
  }
  
  # 3. Eigendecomposition
  decomp       <- eigen(r_mat, symmetric = TRUE)
  sorted_evals <- decomp$values
  sorted_evecs <- decomp$vectors
  
  # Unrotated Component Loadings: Lambda = V * sqrt(D)
  top_evals          <- sorted_evals[seq_len(n_factors)]
  top_evecs          <- sorted_evecs[, seq_len(n_factors), drop = FALSE]
  unrotated_loadings <- top_evecs %*% diag(sqrt(pmax(top_evals, 1e-12)), nrow = n_factors)
  
  comp_names <- paste0("Component ", seq_len(n_factors))
  rownames(unrotated_loadings) <- cols
  colnames(unrotated_loadings) <- comp_names
  
  # Direct Oblimin Multi-Start (Quartimin, normalize = FALSE)
  # Generate exhaustive sign combinations for the retained components
  sign_combinations <- as.matrix(expand.grid(rep(list(c(1, -1)), n_factors)))
  
  best_rot  <- NULL
  best_loss <- Inf
  
  for (row_idx in seq_len(nrow(sign_combinations))) {
    signs   <- sign_combinations[row_idx, ]
    L_start <- unrotated_loadings %*% diag(signs, nrow = n_factors)
    
    rot_try <- tryCatch(
      GPArotation::oblimin(L_start, gam = 0, normalize = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(rot_try)) {
      L_rot <- unclass(rot_try$loadings)
      loss  <- 0
      for (j in 1:(n_factors - 1)) {
        for (k in (j + 1):n_factors) {
          loss <- loss + sum((L_rot[, j]^2) * (L_rot[, k]^2))
        }
      }
      if (loss < best_loss) {
        best_loss <- loss
        best_rot  <- rot_try
      }
    }
  }
  
  rotated_loadings <- unclass(best_rot$loadings)
  rownames(rotated_loadings) <- cols
  colnames(rotated_loadings) <- comp_names
  phi_matrix <- best_rot$Phi
  
  # Anchor Matching & Alignment (Applied if 3 components are extracted)
  find_idx <- function(pattern, target_names) {
    m <- grep(pattern, target_names, ignore.case = TRUE)
    if (length(m) > 0) m[1] else NA_integer_
  }
  
  auth_row <- find_idx("Authentic", cols)
  rel_row  <- find_idx("Reliab", cols)
  law_row  <- find_idx("Lawful", cols)
  
  if (n_factors == 3 && !is.na(auth_row) && !is.na(rel_row)) {
    c1_idx <- which.max(abs(rotated_loadings[auth_row, ]))
    c2_idx <- which.max(abs(rotated_loadings[rel_row, ]))
    c3_idx <- setdiff(seq_len(3), c(c1_idx, c2_idx))
    
    col_order        <- c(c1_idx, c2_idx, c3_idx)
    rotated_loadings <- rotated_loadings[, col_order]
    phi_matrix       <- phi_matrix[col_order, col_order]
    
    target_signs <- c(
      if (rotated_loadings[auth_row, 1] > 0) -1 else 1,
      if (rotated_loadings[rel_row, 2] < 0)  -1 else 1,
      if (!is.na(law_row) && rotated_loadings[law_row, 3] < 0) -1 else 1
    )
    
    rotated_loadings <- rotated_loadings %*% diag(target_signs)
    phi_matrix       <- diag(target_signs) %*% phi_matrix %*% diag(target_signs)
  }
  
  # 6. Final Matrices & Assignment Table
  pattern_loadings <- as.data.frame(rotated_loadings)
  structure_matrix <- as.data.frame(as.matrix(pattern_loadings) %*% phi_matrix)
  colnames(structure_matrix) <- comp_names
  rownames(structure_matrix) <- cols
  
  abs_pattern         <- abs(as.matrix(pattern_loadings))
  primary_comp_idx    <- max.col(abs_pattern)
  primary_assignments <- comp_names[primary_comp_idx]
  
  assignments_df <- data.frame(
    Item              = cols,
    Primary_Component = primary_assignments,
    Primary_Loading   = round(as.matrix(pattern_loadings)[cbind(seq_along(cols), primary_comp_idx)], 3),
    stringsAsFactors  = FALSE
  )
  
  list(
    eigenvalues         = sorted_evals,
    n_factors           = n_factors,
    pattern_loadings    = pattern_loadings,
    structure_matrix    = structure_matrix,
    factor_correlations = as.data.frame(phi_matrix),
    assignments_table   = assignments_df
  )
}

compute_internal_reliabilty <- function(df = NULL, corr_matrix = NULL) {

  # Standard Cronbach's Alpha (Raw Data)
  std_alpha <- psych::alpha(items_clean, check.keys = FALSE)

  # Ordinal Cronbach's Alpha (Polychoric Correlation Matrix)
  ord_alpha <- psych::alpha(cor_polychoric, check.keys = FALSE)

  # McDonald's Omega (1-Factor MINRES Solution)
  om_fit <- suppressWarnings(
    omega(cor_polychoric, nfactors = 1, fm = "minres", plot = FALSE)
  )

  # Print Global Summary Metrics
  cat("--- Scale Reliability Summary ---\n")
  cat(sprintf("Standard Cronbach's Alpha (Raw) : %.3f\n", std_alpha$total$raw_alpha))
  cat(sprintf("Ordinal Cronbach's Alpha (Poly) : %.3f\n", ord_alpha$total$raw_alpha))
  cat(sprintf("McDonald's Omega Total (ω_t)     : %.3f\n\n", om_fit$omega.tot))

  # Single-Factor Loadings and Error Variances
  cat("--- Single-Factor Standardized Loadings (MINRES) ---\n")
  loadings_table <- data.frame(
    Item             = colnames(items_clean),
    Std_Loading_g    = round(as.vector(om_fit$schmid$sl[, "g"]), 3),
    Error_Variance_u = round(as.vector(om_fit$schmid$sl[, "u2"]), 3)
  )
  print(loadings_table, row.names = FALSE)
  cat("\n")

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
}