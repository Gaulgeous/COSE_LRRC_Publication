source("helpers.R")

# Pairwise and Matrix Correlation Utilities ----

#' Safely compute Kendall rank correlation coefficient
#'
#' Computes Kendall's tau-b between two vectors, returning zero if either vector
#' possesses zero variance or if the calculation yields non-finite values.
#'
#' @param x Numeric or ordinal vector.
#' @param y Numeric or ordinal vector of length matching `x`.
#'
#' @return A numeric scalar representing Kendall's tau-b correlation, or `0.0`
#'   if variance is degenerate or non-estimable.
#' @noRd
#'
#' @examples
#' safe_kendall(c(1, 2, 3), c(2, 3, 4))
#' safe_kendall(c(1, 1, 1), c(2, 3, 4))
safe_kendall <- function(x, y) {
  # Return 0.0 if either vector has zero variance (all values identical)
  if (length(unique(x)) <= 1 || length(unique(y)) <= 1) {
    return(0.0)
  }

  val <- suppressWarnings(DescTools::KendallTauB(x, y, conf.level = NA))
  if (is.na(val) || is.nan(val)) 0.0 else as.numeric(val)
}


#' Project indefinite matrix to nearest positive semi-definite matrix
#'
#' Inspects eigenvalues of a symmetric matrix and conditionally computes the Higham
#' nearest positive semi-definite (PSD) correlation projection using [Matrix::nearPD()]
#' only when negative eigenvalues exceed the tolerance.
#'
#' @param mat Square numeric matrix.
#' @param tol Numeric. Non-negative tolerance for negative eigenvalue detection (default 1e-8).
#' @param corr Logical. Whether to enforce unit diagonal elements in the projection.
#'
#' @return A positive semi-definite matrix of class `matrix` with identical dimension names.
#' @export
#'
#' @examples
#' set.seed(42)
#' m <- matrix(c(1, 0.9, 0.9, 0.9, 1, 0.9, 0.9, 0.9, 1), nrow = 3)
#' project_near_psd(m)
project_near_psd <- function(mat, tol = 1e-8, corr = TRUE) {
  evals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values

  # If already positive semi-definite, return unchanged to preserve raw correlations
  if (min(evals) >= -tol) {
    return(mat)
  }

  pd_result <- Matrix::nearPD(mat, corr = corr, keepDiag = TRUE)
  out <- as.matrix(pd_result$mat)
  dimnames(out) <- dimnames(mat)
  out
}


#' Compute positive semi-definite Kendall correlation matrix via pairwise evaluation
#'
#' Evaluates pairwise Kendall's tau-b correlations across columns using complete
#' cases, setting non-estimable pairs to zero, and projects the resulting matrix
#' to the nearest positive semi-definite form if necessary.
#'
#' @param df Data frame or matrix containing numeric or ordinal columns.
#'
#' @return A symmetric numeric matrix of positive semi-definite Kendall rank correlations.
#' @export
#'
#' @examples
#' set.seed(42)
#' sample_df <- data.frame(
#'   Item1 = sample(0:2, 20, replace = TRUE),
#'   Item2 = sample(0:2, 20, replace = TRUE),
#'   Item3 = sample(0:2, 20, replace = TRUE)
#' )
#' compute_psd_kendall(sample_df)
compute_psd_kendall <- function(df) {
  cols <- colnames(df)
  p <- length(cols)

  tau_mat <- matrix(0, nrow = p, ncol = p, dimnames = list(cols, cols))
  diag(tau_mat) <- 1.0

  pairs <- utils::combn(cols, 2, simplify = FALSE)

  for (pair in pairs) {
    u <- pair[1]
    v <- pair[2]

    pair_data <- stats::na.omit(df[, c(u, v)])

    if (nrow(pair_data) >= 3) {
      tau_val <- safe_kendall(pair_data[[u]], pair_data[[v]])
      tau_mat[u, v] <- tau_val
      tau_mat[v, u] <- tau_val
    }
  }

  min_eigen <- min(eigen(tau_mat, symmetric = TRUE, only.values = TRUE)$values)

  if (min_eigen < 0) {
    tau_psd <- project_near_psd(tau_mat)
  } else {
    tau_psd <- tau_mat
  }

  return(tau_psd)
}


# Sampling Adequacy and Scree Analysis ----

#' Assess sampling adequacy and factorability via Bartlett and KMO tests
#'
#' Executes Bartlett's Test of Sphericity and the Kaiser-Meyer-Olkin (KMO) Measure
#' of Sampling Adequacy, displaying overall and individual item MSA diagnostics.
#'
#' @param cor_mat Square numeric correlation matrix.
#' @param n_obs Integer. Number of observations used to estimate `cor_mat`.
#'
#' @return No return value, called for printing diagnostic results.
#' @export
#'
#' @examples
#' data(mtcars)
#' c_mat <- cor(mtcars[, 1:4])
#' assess_factorability(c_mat, n_obs = nrow(mtcars))
assess_factorability <- function(cor_mat, n_obs) {
  cat("--- Bartlett's Test of Sphericity ---\n")
  bartlett_res <- psych::cortest.bartlett(cor_mat, n = n_obs)
  cat(sprintf("Chi-Square: %.3f | df: %d | p-value: %.3e\n\n",
              bartlett_res$chisq, bartlett_res$df, bartlett_res$p.value))

  cat("--- Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy ---\n")
  kmo_res <- psych::KMO(cor_mat)

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


#' Compute PCA scree analysis and Kaiser-Guttman factor retention
#'
#' Resolves input data to a positive semi-definite correlation matrix, performs
#' spectral decomposition, and evaluates cumulative variance alongside Kaiser-Guttman
#' retention criteria.
#'
#' @param df_or_corr Data frame, matrix, or precomputed correlation matrix.
#' @param kaiser_cutoff Numeric. Minimum eigenvalue benchmark for retention (default 1.0).
#' @param verbose Logical. If `TRUE`, prints variance decomposition table.
#'
#' @return A named list containing:
#'   \item{cor_matrix}{The analyzed positive semi-definite correlation matrix.}
#'   \item{eigenvalues}{Numeric vector of sorted eigenvalues.}
#'   \item{pct_variance}{Percentage of total variance explained per component.}
#'   \item{cum_variance}{Cumulative percentage of variance explained.}
#'   \item{n_factors_kaiser}{Count of components meeting `kaiser_cutoff`.}
#'   \item{kaiser_cutoff}{Cutoff value applied.}
#'   \item{variance_table}{Data frame summarizing component variance and retention.}
#'   \item{total_variance_kaiser}{Cumulative variance accounted for by retained components.}
#' @export
#'
#' @examples
#' data(mtcars)
#' compute_pca_scree_analysis(mtcars[, 1:4], verbose = FALSE)
compute_pca_scree_analysis <- function(df_or_corr, kaiser_cutoff = 1.0, verbose = TRUE) {
  parsed <- resolve_corr_matrix(df_or_corr = df_or_corr)
  r_mat  <- ensure_psd(parsed$cor_matrix)

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


# Parallel Analysis and Simulation ----

#' Perform Horn parallel analysis using exact Monte Carlo simulation
#'
#' Generates uncorrelated Gaussian noise matrices matching sample observations and
#' variable counts, comparing empirical eigenvalues against simulation quantiles.
#'
#' @param df Optional raw data frame or matrix of observations.
#' @param corr_matrix Optional precomputed correlation matrix.
#' @param df_or_corr Generic parameter accepting either a data frame or correlation matrix.
#' @param n_observations Integer. Number of observations underlying the matrix (default 26).
#' @param n_simulations Integer. Number of Monte Carlo iterations (default 1000).
#' @param percentile Numeric. Retention quantile cutoff for synthetic noise (default 95.0).
#' @param kaiser_cutoff Numeric. Kaiser retention benchmark threshold (default 1.0).
#' @param plot_comparison Logical. If `TRUE`, renders base R scree comparison plot.
#' @param verbose Logical. If `TRUE`, prints summary retention comparisons.
#'
#' @return Invisibly returns a named list with components:
#'   \item{n_factors_parallel}{Integer count of components exceeding the simulation quantile.}
#'   \item{n_factors_kaiser}{Integer count of components exceeding `kaiser_cutoff`.}
#'   \item{empirical_evals}{Numeric vector of observed empirical eigenvalues.}
#'   \item{pct_variance}{Percentage of variance explained per component.}
#'   \item{cum_variance}{Cumulative percentage of variance explained.}
#'   \item{synthetic_percentile}{Quantile values from Monte Carlo simulations.}
#'   \item{synthetic_mean}{Mean eigenvalue across Monte Carlo simulations.}
#'   \item{comparison_table}{Data frame comparing empirical and simulated eigenvalue trajectories.}
#'   \item{scree_summary}{Output object from [compute_pca_scree_analysis()].}
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_df <- data.frame(matrix(rnorm(60), nrow = 20, ncol = 3))
#' perform_horn_parallel_analysis_exact(df = sim_df, n_simulations = 50,
#'                                      plot_comparison = FALSE, verbose = FALSE)
perform_horn_parallel_analysis_exact <- function(df = NULL,
                                                 corr_matrix = NULL,
                                                 df_or_corr = NULL,
                                                 n_observations = 26,
                                                 n_simulations = 1000,
                                                 percentile = 95.0,
                                                 kaiser_cutoff = 1.0,
                                                 plot_comparison = TRUE,
                                                 verbose = TRUE) {
  parsed <- resolve_corr_matrix(df = df, corr_matrix = corr_matrix, df_or_corr = df_or_corr)
  n_obs  <- if (!is.na(parsed$n_obs)) parsed$n_obs else n_observations
  r_mat  <- ensure_psd(parsed$cor_matrix)

  scree_res <- compute_pca_scree_analysis(r_mat, kaiser_cutoff = kaiser_cutoff, verbose = FALSE)
  emp_ev    <- scree_res$eigenvalues
  pct_var   <- scree_res$pct_variance
  cum_var   <- scree_res$cum_variance
  n_vars    <- length(emp_ev)

  if (is.data.frame(r_mat) && !(nrow(r_mat) == ncol(r_mat) && all(diag(as.matrix(r_mat)) == 1))) {
    n_obs <- nrow(stats::na.omit(r_mat))
  } else {
    n_obs <- n_observations
  }

  synthetic_evals <- matrix(0, nrow = n_simulations, ncol = n_vars)
  for (i in seq_len(n_simulations)) {
    noise <- matrix(stats::rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
    synthetic_evals[i, ] <- sort(eigen(stats::cor(noise), symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  }

  sim_mu <- colMeans(synthetic_evals)
  sim_q  <- apply(synthetic_evals, 2, function(col) stats::quantile(col, probs = percentile / 100, type = 7))

  retained_mask <- emp_ev > sim_q
  n_retained    <- sum(retained_mask)

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


# Direct Oblimin PCA and Rotation Multi-Start ----

#' Perform PCA with multi-start direct oblimin rotation and anchor alignment
#'
#' Extracts unrotated principal components, executes a multi-start direct oblimin
#' rotation across reflection sign combinations to evade local optima, permutes
#' component order to maximize anchor item mappings, and orients factor signs.
#'
#' @param df Optional raw data frame or matrix.
#' @param corr_matrix Optional precomputed correlation matrix.
#' @param df_or_corr Generic parameter accepting either a data frame or correlation matrix.
#' @param n_factors Integer or `NULL`. Number of factors to extract. If `NULL`,
#'   dynamically determines dimensions via parallel analysis.
#' @param n_observations Integer. Number of observations underlying correlation estimates (default 26).
#'
#' @return A named list containing:
#'   \item{eigenvalues}{Numeric vector of sorted eigenvalues.}
#'   \item{n_factors}{Integer count of extracted components.}
#'   \item{pattern_loadings}{Data frame of rotated pattern coefficients.}
#'   \item{structure_matrix}{Data frame of structure matrix coefficients.}
#'   \item{factor_correlations}{Data frame of component inter-correlations (Phi).}
#'   \item{assignments_table}{Data frame assigning items to primary components by max absolute loading.}
#' @export
#'
#' @examples
#' data(mtcars)
#' c_mat <- cor(mtcars[, 1:6])
#' perform_abels_pca_oblimin_exact(corr_matrix = c_mat, n_factors = 2, n_observations = nrow(mtcars))
perform_abels_pca_oblimin_exact <- function(df = NULL,
                                            corr_matrix = NULL,
                                            df_or_corr = NULL,
                                            n_factors = 3,
                                            n_observations = 26) {
  parsed <- resolve_corr_matrix(df = df, corr_matrix = corr_matrix, df_or_corr = df_or_corr)
  r_mat  <- ensure_psd(parsed$cor_matrix)
  cols   <- parsed$labels
  n_obs  <- if (!is.na(parsed$n_obs)) parsed$n_obs else n_observations

  # Fall back to Horn's parallel analysis if component count is unassigned
  if (is.null(n_factors)) {
    pa_res    <- perform_horn_parallel_analysis_exact(r_mat, n_observations = n_obs,
                                                      plot_comparison = FALSE, verbose = FALSE)
    n_factors <- max(1, pa_res$n_factors_parallel)
  }

  decomp       <- eigen(r_mat, symmetric = TRUE)
  sorted_evals <- decomp$values
  sorted_evecs <- decomp$vectors

  top_evals          <- sorted_evals[seq_len(n_factors)]
  top_evecs          <- sorted_evecs[, seq_len(n_factors), drop = FALSE]
  unrotated_loadings <- top_evecs %*% diag(sqrt(pmax(top_evals, 1e-12)), nrow = n_factors)

  comp_names <- paste0("Component ", seq_len(n_factors))
  dimnames(unrotated_loadings) <- list(cols, comp_names)

  # Multi-start oblimin rotation over 2^k sign reflections to avoid local minima
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
  phi_matrix       <- best_rot$Phi
  dimnames(rotated_loadings) <- list(cols, comp_names)
  dimnames(phi_matrix)       <- list(comp_names, comp_names)

  find_idx <- function(pattern, target_names) {
    m <- grep(pattern, target_names, ignore.case = TRUE)
    if (length(m) > 0) m[1] else NA_integer_
  }

  auth_row <- find_idx("Authentic", cols)
  rel_row  <- find_idx("Reliab", cols)
  law_row  <- find_idx("Lawful", cols)

  # Match 3-factor models to authentic, reliable, and lawful anchor items
  if (n_factors == 3 && !is.na(auth_row) && !is.na(rel_row)) {
    # Evaluate all 6 permutations to guarantee a valid 1-to-1 column mapping
    all_perms <- list(
      c(1, 2, 3), c(1, 3, 2),
      c(2, 1, 3), c(2, 3, 1),
      c(3, 1, 2), c(3, 2, 1)
    )

    # Select the permutation maximizing total absolute loading across anchors
    perm_scores <- sapply(all_perms, function(p) {
      abs(rotated_loadings[auth_row, p[1]]) +
        abs(rotated_loadings[rel_row,  p[2]]) +
        (if (!is.na(law_row)) abs(rotated_loadings[law_row, p[3]]) else 0)
    })

    col_order <- all_perms[[which.max(perm_scores)]]

    rotated_loadings <- rotated_loadings[, col_order, drop = FALSE]
    phi_matrix       <- phi_matrix[col_order, col_order, drop = FALSE]

    # Orient signs so anchor variables point in the positive direction
    target_signs <- c(
      if (rotated_loadings[auth_row, 1] < 0) -1 else 1,
      if (rotated_loadings[rel_row,  2] < 0) -1 else 1,
      if (!is.na(law_row) && rotated_loadings[law_row, 3] < 0) -1 else 1
    )

    sign_diag        <- diag(target_signs, nrow = 3, ncol = 3)
    rotated_loadings <- rotated_loadings %*% sign_diag
    phi_matrix       <- sign_diag %*% phi_matrix %*% sign_diag

    dimnames(rotated_loadings) <- list(cols, comp_names)
    dimnames(phi_matrix)       <- list(comp_names, comp_names)
  }

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


# Internal Reliability Diagnostics ----

#' Compute scale internal consistency metrics
#'
#' Evaluates standard Cronbach's alpha on raw scores, ordinal Cronbach's alpha
#' from polychoric correlation matrices, McDonald's omega total using single-factor
#' MINRES estimation, and item-deleted reliability diagnostics.
#'
#' @param df Data frame or matrix of raw item responses (defaults to `items_clean`).
#' @param corr_matrix Matrix of polychoric correlations (defaults to `cor_polychoric`).
#'
#' @return No return value, called for printing reliability summaries.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_internal_reliabilty(df = items_clean, corr_matrix = cor_polychoric)
#' }
compute_internal_reliabilty <- function(df = NULL, corr_matrix = NULL) {
  # Standard Cronbach's Alpha (Raw Data)
  std_alpha <- psych::alpha(items_clean, check.keys = FALSE)

  # Ordinal Cronbach's Alpha (Polychoric Correlation Matrix)
  ord_alpha <- psych::alpha(cor_polychoric, check.keys = FALSE)

  # McDonald's Omega (1-Factor MINRES Solution)
  om_fit <- suppressWarnings(
    psych::omega(cor_polychoric, nfactors = 1, fm = "minres", plot = FALSE)
  )

  cat("--- Scale Reliability Summary ---\n")
  cat(sprintf("Standard Cronbach's Alpha (Raw) : %.3f\n", std_alpha$total$raw_alpha))
  cat(sprintf("Ordinal Cronbach's Alpha (Poly) : %.3f\n", ord_alpha$total$raw_alpha))
  cat(sprintf("McDonald's Omega Total (\u03c9_t)     : %.3f\n\n", om_fit$omega.tot))

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
    fit_sub <- suppressWarnings(psych::omega(sub_r, nfactors = 1, fm = "minres", plot = FALSE))
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


# Pipeline Orchestration ----

#' Execute complete complex PCA workflow
#'
#' Orchestrates Kendall correlation estimation, 2-step MLE polychoric correlation
#' derivation with Higham projection, factorability diagnostics, exact Monte Carlo
#' parallel analysis, multi-start direct oblimin PCA, and internal scale reliability.
#'
#' @param items Matrix or data frame of observed item responses.
#' @param n_boot Integer. Number of Monte Carlo simulation runs for parallel analysis (default 1000).
#' @param percentile Numeric. Retention quantile cutoff for parallel analysis (default 95.0).
#' @param manual_n_components Integer. Number of components to extract (default 2).
#' @param seed Integer. Random seed for reproducible simulation generation (default 42).
#' @param ... Additional arguments passed to downstream analysis functions.
#'
#' @return Invisibly returns the PCA result list from [perform_abels_pca_oblimin_exact()].
#' @export
#'
#' @examples
#' \dontrun{
#' pca_exact <- run_complex_pca(
#'   items               = items_clean,
#'   n_boot              = 1000,
#'   percentile          = 95.0,
#'   manual_n_components = 2,
#'   seed                = 42
#' )
#' }
run_complex_pca <- function(items,
                            n_boot = 1000,
                            percentile = 95.0,
                            manual_n_components = 2,
                            seed = 42,
                            ...) {
  kendall_matrix <- compute_psd_kendall(items_clean)

  cat("--- Positive Semi-Definite Kendall Tau-b Correlation Matrix ---\n")
  print(round(kendall_matrix, 3))

  # Compute exact 2-step MLE without Taylor approximations or continuity distortion
  #    global = FALSE -> Evaluates thresholds per pair
  #    correct = 0    -> Prevents artificial off-diagonal probability inflation
  #    smooth = FALSE -> Prevents cor.smooth from flattening strong signals
  poly_mle <- psych::polychoric(
    items_clean,
    global  = FALSE,
    correct = 0,
    smooth  = FALSE
  )$rho

  # Apply Higham nearest correlation projection
  cor_polychoric <- project_near_psd(poly_mle)
  dimnames(cor_polychoric) <- dimnames(poly_mle)

  cat("--- Polychoric Matrix (Matching Python Implementation) ---\n")
  print(round(cor_polychoric, 3))

  # Full matrix heatmap without text values
  plot_correlation_heatmap(
    cor_mat = cor_polychoric,
    type    = "full",
    title   = "Polychoric Correlation of Publication's Coverage of Scorecard Components"
  )

  assess_factorability(cor_polychoric, n_obs)

  pa_res <- perform_horn_parallel_analysis_exact(
    df             = items_clean,
    n_observations = nrow(items_clean),
    n_simulations  = n_boot,
    percentile     = percentile
  )

  pca_exact <- perform_abels_pca_oblimin_exact(df = items_clean, n_factors = manual_n_components)

  cat("--- Rotated Pattern Matrix (Table 3 Style) ---\n")
  print(round(pca_exact$pattern_loadings, 3))

  cat("\n--- Structure Matrix (Table 4 Style) ---\n")
  print(round(pca_exact$structure_matrix, 3))

  cat("\n--- Component Inter-Correlation Matrix (Phi) ---\n")
  print(round(pca_exact$factor_correlations, 3))

  cat("\n--- Item-to-Component Primary Assignments ---\n")
  print(pca_exact$assignments_table, row.names = FALSE)

  compute_internal_reliabilty(df = items_clean, corr_matrix = cor_polychoric)

  invisible(pca_exact)
}