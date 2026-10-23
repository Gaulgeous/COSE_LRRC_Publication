# Label Normalization and Matrix Helpers ----

#' Clean and normalize character labels
#'
#' Strips non-breaking spaces and whitespace, optionally converts to lowercase,
#' and recursively standardizes dimension names of matrices and data frames.
#'
#' @param x Character vector, matrix, data frame, or `NULL`.
#' @param to_lower Logical. If `TRUE`, converts characters to lowercase.
#' @param strip_nbsp Logical. If `TRUE`, replaces non-breaking space characters
#'   (`\u00a0`) with standard spaces.
#'
#' @return An object of the same class as `x` with trimmed, standardized character
#'   labels or dimension names.
#' @export
#'
#' @examples
#' clean_labels(c("  Item 1\u00a0 ", "Item 2  "))
#' df <- data.frame(a = 1, row.names = " Row 1\u00a0 ")
#' clean_labels(df)
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


#' Resolve correlation matrix from raw or precomputed inputs
#'
#' Disentangles raw rectangular observation frames from square precomputed correlation
#' matrices, computes Pearson correlation when raw data is provided, cleans dimension
#' labels, and extracts sample size.
#'
#' @param df Optional raw numeric data frame or tibble.
#' @param corr_matrix Optional square numeric correlation matrix.
#' @param df_or_corr Generic parameter accepting either a data frame or a correlation matrix.
#'
#' @return A named list containing:
#'   \item{cor_matrix}{Symmetric numeric correlation matrix with normalized labels.}
#'   \item{n_obs}{Integer observation count, or `NA_integer_` if a precomputed matrix was supplied.}
#'   \item{n_vars}{Integer count of variables (columns).}
#'   \item{labels}{Character vector of variable names.}
#' @export
#'
#' @examples
#' data(mtcars)
#' res <- resolve_corr_matrix(df = mtcars[, 1:4])
#' res$n_vars
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

  r_mat <- clean_labels(r_mat)

  list(
    cor_matrix = r_mat,
    n_obs      = n_obs,
    n_vars     = ncol(r_mat),
    labels     = colnames(r_mat)
  )
}


#' Ensure positive semi-definite matrix properties
#'
#' Verifies eigenvalue non-negativity of a symmetric matrix and conditionally projects
#' non-positive semi-definite (PSD) matrices to the nearest PSD space using [Matrix::nearPD()].
#'
#' @param mat Square numeric matrix.
#' @param tol Numeric. Tolerance threshold for detecting non-positive eigenvalues (default 1e-7).
#' @param corr Logical. Whether to enforce unit diagonal elements in the projected matrix.
#' @param keepDiag Logical. Whether to retain original diagonal elements if `corr = FALSE`.
#' @param ... Additional arguments passed to [Matrix::nearPD()].
#'
#' @return A positive semi-definite numeric matrix of class `matrix` with preserved dimension names.
#' @export
#'
#' @examples
#' set.seed(42)
#' m <- matrix(c(1, 0.95, 0.95, 0.95, 1, 0.95, 0.95, 0.95, 1), nrow = 3)
#' ensure_psd(m)
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


# Dataset Ingestion and Preprocessing ----

#' Extract and preprocess publication scorecard dataset from Excel
#'
#' Imports publication scorecard records from a target Excel sheet, forward-fills
#' merged category headers, cleans and assigns author names as unique row identifiers,
#' extracts numeric scorecard features, and filters complete cases.
#'
#' @param excel_file Character string. File path to the Excel workbook.
#' @param sheet_name Character string or integer. Name or index of sheet to import.
#' @param ... Additional arguments ignored for compatibility.
#'
#' @return A named list containing:
#'   \item{items_clean}{Data frame of complete numeric scorecard observations with publication row names.}
#'   \item{n_obs}{Integer count of complete observations.}
#'   \item{categories_clean}{Character vector of forward-filled category labels matching `items_clean`.}
#' @export
#'
#' @examples
#' \dontrun{
#' info <- extract_dataset("scorecards.xlsx", sheet_name = "Primary_Review")
#' }
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
    items_clean      = items_clean,
    n_obs            = n_obs,
    categories_clean = categories_clean
  ))
}


#' Run preliminary scorecard summary statistics and distance profiles
#'
#' Reads an Excel scorecard sheet preserving two-tier header structures, calculates
#' overall column dimension means, computes category-level decimal profiles, and
#' generates a pairwise Euclidean centroid distance heatmap.
#'
#' @param excel_file Character string. File path to the raw Excel workbook.
#' @param sheet_name Character string or integer. Sheet name or position index.
#' @param ... Additional arguments passed to downstream plotting functions.
#'
#' @return No return value, called for summary printing and heatmap visualization.
#' @export
#'
#' @examples
#' \dontrun{
#' run_preliminary_stats("scorecards.xlsx", sheet_name = "Primary_Review")
#' }
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

  col_names <- as.character(unlist(raw_df[2, 3:ncol(raw_df)]))
  col_names <- trimws(col_names)

  score_data <- raw_df[3:nrow(raw_df), 3:ncol(raw_df), drop = FALSE]

  score_data_num <- as.data.frame(
    lapply(score_data, function(col) as.numeric(as.character(col)))
  )

  col_averages <- colMeans(score_data_num, na.rm = TRUE)
  names(col_averages) <- col_names

  cat("--- Column Averages ---\n")
  print(round(col_averages, 3))
  cat("\n")

  # Extracts Category from the Publication multi-column, coerces score
  # columns to numeric, and calculates the mean for each category
  category_means <- compute_category_means_simple(raw_df)

  cat("--- Category Means Across Scorecard Dimensions ---\n")
  print(round(category_means, 3))
  cat("\n")

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
