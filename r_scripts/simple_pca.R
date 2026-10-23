source("helpers.R")

# Correlation and Factorability Diagnostics ----

#' Compute positive semi-definite Kendall correlation matrix
#'
#' Computes pairwise Kendall rank correlation coefficients (tau-b) across numeric
#' columns and conditionally projects indefinite matrices to the nearest positive
#' semi-definite (PSD) space using [Matrix::nearPD()].
#'
#' @param df Data frame or matrix containing numeric item scores.
#' @param use Character string specifying missing data handling, passed to [stats::cor()].
#'   Defaults to `"pairwise.complete.obs"`.
#'
#' @return A symmetric numeric matrix containing positive semi-definite Kendall tau-b
#'   correlations with row and column names preserved.
#' @export
#'
#' @examples
#' set.seed(42)
#' sample_data <- matrix(sample(0:2, 60, replace = TRUE), ncol = 3,
#'                       dimnames = list(NULL, c("Item1", "Item2", "Item3")))
#' compute_psd_kendall_simple(sample_data)
compute_psd_kendall_simple <- function(df, use = "pairwise.complete.obs") {
  mat <- data.matrix(df)

  tau_mat <- stats::cor(mat, method = "kendall", use = use)

  # Check if indefinite (negative eigenvalues)
  min_eigen <- min(eigen(tau_mat, symmetric = TRUE, only.values = TRUE)$values)

  # Project to nearest positive semi-definite matrix only if negative eigenvalues exist
  if (min_eigen < 0) {
    pd_obj  <- Matrix::nearPD(tau_mat, corr = TRUE, keepDiag = TRUE)
    tau_mat <- as.matrix(pd_obj$mat)
    dimnames(tau_mat) <- list(colnames(df), colnames(df))
  }

  tau_mat
}


#' Assess sampling adequacy and factorability of a correlation matrix
#'
#' Evaluates whether a correlation matrix is suitable for factor extraction using
#' Bartlett's Test of Sphericity and the Kaiser-Meyer-Olkin (KMO) Measure of
#' Sampling Adequacy.
#'
#' @param cor_mat Square numeric correlation matrix.
#' @param n_obs Integer. Number of observations used to estimate `cor_mat`.
#' @param verbose Logical. If `TRUE`, prints formatted summary diagnostics to the console.
#'
#' @return Invisibly returns a named list with components:
#'   \item{bartlett}{Results object from [psych::cortest.bartlett()].}
#'   \item{kmo_overall}{Numeric overall Kaiser-Meyer-Olkin Measure of Sampling Adequacy.}
#'   \item{kmo_items}{Data frame of variable names and individual MSA values.}
#' @export
#'
#' @examples
#' data(mtcars)
#' c_mat <- cor(mtcars[, 1:4])
#' assess_factorability_simple(c_mat, n_obs = nrow(mtcars), verbose = FALSE)
assess_factorability_simple <- function(cor_mat, n_obs, verbose = TRUE) {
  bartlett_res <- psych::cortest.bartlett(cor_mat, n = n_obs)

  kmo_res <- psych::KMO(cor_mat)

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

  invisible(list(
    bartlett    = bartlett_res,
    kmo_overall = kmo_res$MSA,
    kmo_items   = item_msa
  ))
}


# Dimension Retention and Parallel Analysis ----

#' Perform PCA scree analysis and Kaiser-Guttman evaluation
#'
#' Computes eigenvalues and variance explained across components from a raw data
#' matrix or precomputed correlation matrix, applying the Kaiser-Guttman retention rule.
#'
#' @param x Square correlation matrix or raw numeric data frame/matrix.
#' @param kaiser_cutoff Numeric. Minimum eigenvalue threshold for retention (default 1.0).
#' @param verbose Logical. If `TRUE`, prints summary variance and retention output.
#'
#' @return A named list with components:
#'   \item{cor_matrix}{The correlation matrix analyzed.}
#'   \item{eigenvalues}{Numeric vector of non-negative eigenvalues in descending order.}
#'   \item{pct_variance}{Numeric vector of percentage variance explained by each component.}
#'   \item{cum_variance}{Numeric vector of cumulative percentage variance explained.}
#'   \item{n_factors_kaiser}{Integer count of components meeting `kaiser_cutoff`.}
#'   \item{variance_table}{Data frame summarizing eigenvalues and retention decisions.}
#'   \item{total_variance_kaiser}{Numeric cumulative variance accounted for by retained components.}
#' @export
#'
#' @examples
#' data(mtcars)
#' compute_pca_scree_analysis_simple(mtcars[, 1:5], verbose = FALSE)
compute_pca_scree_analysis_simple <- function(x, kaiser_cutoff = 1.0, verbose = TRUE) {
  if (is.matrix(x) && nrow(x) == ncol(x) && all(abs(diag(x) - 1) < 1e-5)) {
    r_mat <- x
  } else {
    r_mat <- stats::cor(data.matrix(x), use = "pairwise.complete.obs")
  }

  ev <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values

  # Guard against negative eigenvalues from rank correlations
  ev_clean <- pmax(ev, 0)

  pct_var  <- (ev_clean / sum(ev_clean)) * 100
  cum_var  <- cumsum(pct_var)
  n_kaiser <- sum(ev_clean >= kaiser_cutoff)

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


#' Perform Horn parallel analysis using Monte Carlo simulation
#'
#' Implements Horn's parallel analysis by generating synthetic Gaussian noise matrices
#' matching empirical sample size and variable dimensions, comparing empirical
#' eigenvalues against simulation percentiles.
#'
#' @param x Square correlation matrix or raw data frame/matrix.
#' @param n_obs Integer. Sample size. Required if `x` is a correlation matrix.
#' @param n_simulations Integer. Number of Monte Carlo simulation iterations (default 1000).
#' @param percentile Numeric. Percentile of simulated eigenvalues used as cutoff (default 95.0).
#' @param kaiser_cutoff Numeric. Kaiser benchmark cutoff line (default 1.0).
#' @param plot_comparison Logical. If `TRUE`, draws base R scree plot comparing empirical
#'   and simulated eigenvalue trajectories.
#' @param verbose Logical. If `TRUE`, prints tabular summary diagnostics to the console.
#'
#' @return Invisibly returns a named list with components:
#'   \item{n_factors_parallel}{Integer count of components exceeding the simulation quantile.}
#'   \item{n_factors_kaiser}{Integer count of components exceeding `kaiser_cutoff`.}
#'   \item{empirical_evals}{Numeric vector of observed empirical eigenvalues.}
#'   \item{pct_variance}{Numeric vector of percentage variance explained per component.}
#'   \item{cum_variance}{Numeric vector of cumulative percentage variance explained.}
#'   \item{synthetic_percentile}{Numeric vector of simulated eigenvalue quantile thresholds.}
#'   \item{synthetic_mean}{Numeric vector of mean simulated eigenvalues.}
#'   \item{comparison_table}{Data frame cross-tabulating empirical and simulated values.}
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' perform_horn_parallel_analysis_simple(sim_mat, n_simulations = 50, plot_comparison = FALSE, verbose = FALSE)
perform_horn_parallel_analysis_simple <- function(x,
                                                  n_obs = NULL,
                                                  n_simulations = 1000,
                                                  percentile = 95.0,
                                                  kaiser_cutoff = 1.0,
                                                  plot_comparison = TRUE,
                                                  verbose = TRUE) {
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

  emp_ev <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values
  emp_ev <- pmax(emp_ev, 0)
  pct_var <- (emp_ev / sum(emp_ev)) * 100
  cum_var <- cumsum(pct_var)

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

  if (plot_comparison) {
    comp_idx <- seq_len(n_vars)
    y_limits <- range(c(emp_ev, sim_q, kaiser_cutoff))

    plot(comp_idx, emp_ev, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
         ylim = y_limits, xlab = "Component Number", ylab = "Eigenvalue",
         main = "Comparison of Horn's Parallel Analysis & Scree Plot", xaxt = "n")
    axis(1, at = comp_idx)
    lines(comp_idx, sim_q, type = "b", pch = 15, col = "#d62728", lty = 2, lwd = 1.75)
    lines(comp_idx, sim_mu, type = "b", pch = 17, col = "#7f7f7f", lty = 3, lwd = 1.25)
    abline(h = kaiser_cutoff, col = "black", lty = 3)
    grid()

    legend("topright",
           legend = c("Empirical Eigenvalues",
                      sprintf("Simulated %dth Percentile", as.integer(percentile)),
                      "Simulated Mean",
                      sprintf("Kaiser-Gutmann Cut-Off", kaiser_cutoff)),
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


# Principal Component Analysis and Rotations ----

#' Perform PCA with direct oblimin rotation and anchor alignment
#'
#' Conducts principal component analysis on a correlation matrix using oblique
#' direct oblimin rotation, computes pattern and structure matrices, and reflects
#' signs so primary indicator variables have positive loadings.
#'
#' @param x Square numeric correlation matrix or data frame.
#' @param n_factors Integer. Number of principal components to extract.
#' @param n_obs Integer. Number of observations underlying the correlation matrix.
#' @param align_anchors Logical. If `TRUE`, inverts column signs so that the variable
#'   with the maximum absolute loading has a positive coefficient.
#'
#' @return A named list containing:
#'   \item{eigenvalues}{Numeric vector of sorted eigenvalues.}
#'   \item{n_factors}{Integer count of extracted components.}
#'   \item{pattern_loadings}{Data frame of rotated pattern loadings (P).}
#'   \item{structure_matrix}{Data frame of structure coefficients (S = P %*% Phi).}
#'   \item{factor_correlations}{Data frame of component inter-correlations (Phi).}
#'   \item{assignments_table}{Data frame assigning variables to primary components by maximum absolute loading.}
#' @export
#'
#' @examples
#' data(mtcars)
#' c_mat <- cor(mtcars[, 1:6])
#' perform_abels_pca_oblimin_simple(c_mat, n_factors = 2, n_obs = nrow(mtcars))
perform_abels_pca_oblimin_simple <- function(x, n_factors, n_obs, align_anchors = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square correlation matrix.")
  }

  raw_evals <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  sorted_evals <- sort(raw_evals, decreasing = TRUE)

  pca_fit <- psych::principal(
    r        = x,
    nfactors = n_factors,
    rotate   = "oblimin",
    n.obs    = n_obs,
    covar    = FALSE
  )

  pattern_loadings <- unclass(pca_fit$loadings)
  phi_matrix       <- pca_fit$Phi

  if (is.null(phi_matrix)) {
    # If n_factors == 1 or orthogonal fallback occurs
    phi_matrix <- diag(n_factors)
    colnames(phi_matrix) <- rownames(phi_matrix) <- colnames(pattern_loadings)
  }

  # Ensure dominant indicator on each component loads positively
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

  structure_matrix <- pattern_loadings %*% phi_matrix
  dimnames(structure_matrix) <- dimnames(pattern_loadings)

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

  list(
    eigenvalues         = sorted_evals,
    n_factors           = n_factors,
    pattern_loadings    = round(as.data.frame(pattern_loadings), 4),
    structure_matrix    = round(as.data.frame(structure_matrix), 4),
    factor_correlations = round(as.data.frame(phi_matrix), 4),
    assignments_table   = assignments_df
  )
}


# Scale Reliability and Boundary Analysis ----

#' Compute raw, ordinal, and omega internal consistency reliability
#'
#' Evaluates scale consistency by computing standard Cronbach's alpha from raw
#' scores, ordinal Cronbach's alpha from polychoric correlations, McDonald's total
#' omega (omega_t) using a single-factor MINRES solution, and item-deleted diagnostics.
#'
#' @param df Data frame or matrix of ordinal or Likert item responses.
#' @param verbose Logical. If `TRUE`, prints summary diagnostics and item-deleted metrics.
#'
#' @return Invisibly returns a named list with components:
#'   \item{summary_metrics}{Data frame comparing standard alpha, ordinal alpha, and McDonald's omega.}
#'   \item{loadings_table}{Data frame of single-factor standardized loadings and uniqueness values.}
#'   \item{diagnostics_table}{Data frame of item-deleted alpha and omega diagnostics.}
#'   \item{cor_polychoric}{Estimated polychoric correlation matrix.}
#'   \item{raw_alpha_obj}{Output object from [psych::alpha()] on raw data.}
#'   \item{poly_alpha_obj}{Output object from [psych::alpha()] on polychoric matrix.}
#'   \item{omega_obj}{Output object from [psych::omega()].}
#' @export
#'
#' @examples
#' set.seed(42)
#' df_items <- data.frame(matrix(sample(1:5, 100, replace = TRUE), ncol = 5))
#' compute_internal_reliability_simple(df_items, verbose = FALSE)
compute_internal_reliability_simple <- function(df, verbose = TRUE) {
  items_clean <- stats::na.omit(as.data.frame(df))
  n_obs       <- nrow(items_clean)
  n_items     <- ncol(items_clean)
  item_names  <- colnames(items_clean)

  std_alpha <- psych::alpha(items_clean, check.keys = FALSE, warnings = FALSE)

  # Two-step polychoric estimation for ordinal scales
  poly_obj       <- suppressWarnings(psych::polychoric(items_clean))
  cor_polychoric <- poly_obj$rho
  ord_alpha      <- psych::alpha(cor_polychoric, check.keys = FALSE, warnings = FALSE)

  om_fit <- suppressWarnings(
    psych::omega(cor_polychoric, nfactors = 1, fm = "minres", plot = FALSE)
  )

  omega_dropped <- sapply(seq_len(n_items), function(i) {
    sub_r   <- cor_polychoric[-i, -i]
    fit_sub <- suppressWarnings(
      psych::omega(sub_r, nfactors = 1, fm = "minres", plot = FALSE)
    )
    fit_sub$omega.tot
  })

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


#' Compute theoretical component score coordinate limits
#'
#' Computes the minimum and maximum possible component scores given linear scoring
#' weights and item rating boundaries, identifying extrema via weight signs.
#'
#' @param B Numeric scoring weight matrix (items in rows, components in columns).
#' @param items_data Matrix or data frame of observed item responses used to calculate
#'   standardization means and standard deviations.
#' @param min_val Numeric. Minimum possible rating on the original item response scale (default 0).
#' @param max_val Numeric. Maximum possible rating on the original item response scale (default 2).
#'
#' @return A data frame containing:
#'   \item{Factor}{Factor or component identifier string.}
#'   \item{Min}{Theoretical lower score bound.}
#'   \item{Max}{Theoretical upper score bound.}
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_items <- matrix(sample(0:2, 50, replace = TRUE), ncol = 5)
#' B_weights <- matrix(c(0.2, 0.4, -0.1, 0.3, 0.2), ncol = 1,
#'                     dimnames = list(NULL, "Dim1"))
#' get_factor_plot_limits(B_weights, sim_items, min_val = 0, max_val = 2)
get_factor_plot_limits <- function(B, items_data, min_val = 0, max_val = 2) {
  means <- colMeans(items_data, na.rm = TRUE)
  sds   <- apply(items_data, 2, stats::sd, na.rm = TRUE)

  # Standardized item boundaries
  z_min <- (min_val - means) / sds
  z_max <- (max_val - means) / sds

  n_factors <- ncol(B)
  factor_names <- colnames(B)
  if (is.null(factor_names)) factor_names <- paste0("Factor_", seq_len(n_factors))

  limits <- data.frame(
    Factor = factor_names,
    Min    = numeric(n_factors),
    Max    = numeric(n_factors)
  )

  for (k in seq_len(n_factors)) {
    w <- B[, k]

    # Positive weights take z_max to maximize and z_min to minimize; negative weights invert
    x_max <- ifelse(w >= 0, z_max, z_min)
    x_min <- ifelse(w >= 0, z_min, z_max)

    limits$Max[k] <- sum(w * x_max)
    limits$Min[k] <- sum(w * x_min)
  }

  return(limits)
}


# Pipeline Orchestration ----

#' Execute complete psychometric PCA pipeline
#'
#' Orchestrates the end-to-end PCA workflow: computes Kendall and polychoric
#' correlations, performs factorability diagnostics, tests dimension retention
#' via parallel analysis, extracts rotated oblimin components, and computes
#' internal consistency reliability metrics.
#'
#' @param items Matrix or data frame of numeric item responses.
#' @param n_obs Integer. Number of observations.
#' @param n_boot Integer. Number of Monte Carlo simulation runs for parallel analysis (default 1000).
#' @param percentile Numeric. Retention quantile cutoff for parallel analysis (default 95.0).
#' @param manual_n_components Integer or `NULL`. Specified number of components to extract;
#'   if `NULL`, automatically adopts the number suggested by parallel analysis.
#' @param ... Additional arguments forwarded to downstream routines.
#'
#' @return Invisibly returns a named list with components:
#'   \item{pca_res}{Results list from [perform_abels_pca_oblimin_simple()].}
#'   \item{pa_res}{Results list from [perform_horn_parallel_analysis_simple()].}
#' @export
#'
#' @examples
#' \dontrun{
#' info <- run_simple_pca(
#'   items               = items_clean,
#'   n_obs               = 90,
#'   n_boot              = 500,
#'   percentile          = 95.0,
#'   manual_n_components = 2
#' )
#' }
run_simple_pca <- function(items,
                           n_obs,
                           n_boot = 1000,
                           percentile = 95.0,
                           manual_n_components = 2,
                           ...) {
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

  plot_correlation_heatmap(
    cor_mat = cor_polychoric,
    type    = "full",
    title   = "Polychoric Correlation of Publication's Coverage of Scorecard Components"
  )

  factorability_res <- assess_factorability_simple(
    cor_mat = cor_polychoric,
    n_obs   = n_obs,
    verbose = TRUE
  )

  pa_res <- perform_horn_parallel_analysis_simple(
    x               = cor_polychoric,
    n_obs           = n_obs,
    n_simulations   = n_boot,
    percentile      = percentile,
    plot_comparison = TRUE,
    verbose         = TRUE
  )

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

