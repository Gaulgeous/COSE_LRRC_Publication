# String & Label Normalizer
clean_labels <- function(x, to_lower = FALSE, strip_nbsp = TRUE) {
  if (is.null(x)) return(NULL)
  
  if (is.character(x)) {
    res <- if (strip_nbsp) gsub("\u00a0", " ", x) else x
    res <- trimws(res)
    return(if (to_lower) tolower(res) else res)
  }
  
  # Recursively clean matrices or data frames
  if (is.matrix(x) || is.data.frame(x)) {
    if (!is.null(colnames(x))) colnames(x) <- clean_labels(colnames(x), to_lower = to_lower, strip_nbsp = strip_nbsp)
    if (!is.null(rownames(x))) rownames(x) <- clean_labels(rownames(x), to_lower = to_lower, strip_nbsp = strip_nbsp)
    return(x)
  }
  
  x
}

# Input Parser: Disentangle Raw Dataframe vs. Precomputed Correlation Matrix
resolve_corr_matrix <- function(df = NULL, corr_matrix = NULL, df_or_corr = NULL) {
  input <- if (!is.null(corr_matrix)) {
    corr_matrix
  } else if (!is.null(df)) {
    df
  } else if (!is.null(df_or_corr)) {
    df_or_corr
  } else {
    stop("Must provide either a dataframe (`df`) or correlation matrix (`corr_matrix`).")
  }
  
  if (is.matrix(input) || is.data.frame(input)) {
    mat <- as.matrix(input)
    is_square <- nrow(mat) == ncol(mat)
    is_corr   <- is_square && all(abs(diag(mat) - 1) < 1e-6, na.rm = TRUE)
    
    if (is_corr) {
      r_mat <- mat
      n_obs <- NA_integer_
    } else {
      clean_df <- na.omit(as.data.frame(lapply(as.data.frame(input), as.numeric)))
      n_obs    <- nrow(clean_df)
      r_mat    <- stats::cor(clean_df)
    }
  } else {
    stop("Input must be a numeric dataframe, tibble, or correlation matrix.")
  }
  
  # Standardize labels immediately
  r_mat <- clean_labels(r_mat)
  
  list(
    cor_matrix = r_mat,
    n_obs      = n_obs,
    n_vars     = ncol(r_mat),
    labels     = colnames(r_mat)
  )
}

# Positive Semi-Definite (PSD) Guard
ensure_psd <- function(mat, tol = 1e-7, corr = TRUE, keepDiag = TRUE, ...) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  
  evals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  
  # Only project via nearPD if non-positive definite
  if (any(evals <= tol)) {
    orig_names <- dimnames(mat)
    pd_obj     <- Matrix::nearPD(mat, corr = corr, keepDiag = keepDiag, ...)
    mat        <- as.matrix(pd_obj$mat)
    dimnames(mat) <- orig_names
  }
  
  mat
}