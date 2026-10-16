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

extract_dataset <- function(excel_file,
                          sheet_name,
                          ...) {
  
  # Skipping row 1 makes Row 2 ("Category", "Author", and the 9 items) the column names
  df_raw <- as.data.frame(readxl::read_excel(excel_file, sheet = sheet_name, skip = 1))
  
  # Clean column headers (trim leading/trailing whitespace)
  colnames(df_raw) <- trimws(colnames(df_raw))
  
  # Identify metadata columns
  cat_col <- grep("^category$", colnames(df_raw), ignore.case = TRUE, value = TRUE)[1]
  pub_col <- grep("^(author|publication)$", colnames(df_raw), ignore.case = TRUE, value = TRUE)[1]
  
  if (is.na(pub_col)) pub_col <- colnames(df_raw)[2]
  if (is.na(cat_col)) cat_col <- colnames(df_raw)[1]
  
  # Forward-fill merged Category cells so every publication has its category label
  raw_cat_vec <- trimws(as.character(df_raw[[cat_col]]))
  filled_cat  <- character(length(raw_cat_vec))
  curr_cat    <- NA_character_
  
  for (i in seq_along(raw_cat_vec)) {
    val <- raw_cat_vec[i]
    if (!is.na(val) && val != "" && val != "NA") {
      curr_cat <- val
    }
    filled_cat[i] <- curr_cat
  }
  df_raw[[cat_col]] <- filled_cat
  
  # Filter out empty/trailing rows (must have a valid Author/Publication entry)
  valid_rows <- !is.na(df_raw[[pub_col]]) & 
    trimws(df_raw[[pub_col]]) != "" & 
    trimws(df_raw[[pub_col]]) != "NA"
  
  cleaned_raw <- df_raw[valid_rows, , drop = FALSE]
  
  # Assign unique publication names as row identifiers
  rownames(cleaned_raw) <- make.unique(trimws(as.character(cleaned_raw[[pub_col]])))
  
  # Exclude metadata columns
  metadata_cols <- c(cat_col, pub_col)
  feature_cols  <- setdiff(colnames(cleaned_raw), metadata_cols)
  
  # Coerce scorecard items to numeric matrix and wrap into data.frame
  # Row names (author labels) are preserved automatically
  feature_mat <- apply(cleaned_raw[, feature_cols, drop = FALSE], 2, function(x) {
    as.numeric(trimws(as.character(x)))
  })
  rownames(feature_mat) <- rownames(cleaned_raw)
  
  # Keep complete cases
  items_clean <- stats::na.omit(as.data.frame(feature_mat))
  
  # Retain matched category vector for downstream category analyses
  categories_clean <- cleaned_raw[rownames(items_clean), cat_col]
  
  # Sample size verification (N = 89)
  n_obs <- nrow(items_clean)
  
  cat(sprintf("Successfully imported %d publications across %d scorecard items.\n", 
              n_obs, ncol(items_clean)))
  print(dim(items_clean))
  
  # Expected: 89 rows (90 minus Frowis_2020) by 9 columns
  cat(sprintf("Cleaned items dimensions: %d observations by %d scorecard items.\n", 
              nrow(items_clean), ncol(items_clean)))
  print(dim(items_clean))
  
  cat("===================================================================\n")
  cat("DATA SUMMARY\n")
  cat(sprintf("Observations (Units): %d | Variables (Items): %d\n", n_obs, ncol(items_clean)))
  cat("===================================================================\n\n")
  
  
  return(list(
    items_clean = items_clean,
    n_obs = n_obs,
    categories_clean = categories_clean
  ))
}

run_preliminary_stats <- function(excel_file,
                                  sheet_name,
                                  ...) {
  # Read without default column names so both header rows (Rows 1 & 2) 
  # are preserved intact as data rows for the parser
  raw_data <- readxl::read_excel(
    path      = excel_file,
    sheet     = sheet_name,
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
    title          = "Euclidean Distance Between Publication Category Centroids",
    show_values    = TRUE
  )
  
  cat("--- Pairwise Euclidean Distance Matrix ---\n")
  print(round(category_distances, 3))
  cat("\n")
}