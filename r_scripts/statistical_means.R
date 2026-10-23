# Category Profile Means and Distance Analysis ----

#' Compute category profile means from scorecard tables
#'
#' Ingests raw publication scorecard data in single-row or two-tier merged header
#' format, forward-fills merged category labels across publications, discards
#' empty rows, and computes author counts alongside mean scorecard dimensions per
#' category.
#'
#' @param df Data frame or matrix containing publication scorecard metrics. Accepts
#'   raw imported sheets with either single-tier headers or two-tier merged headers
#'   where metadata is located in rows 1 and 2.
#' @param na.rm Logical. Whether missing numeric scores should be excluded when
#'   computing category dimension means.
#'
#' @return A numeric matrix where rows correspond to unique categories and columns
#'   contain the author count followed by the rounded decimal mean for each
#'   scorecard dimension.
#' @export
#'
#' @examples
#' sample_df <- data.frame(
#'   Category = c("Empirical", "", "Theoretical", ""),
#'   Author   = c("Smith (2020)", "Jones (2021)", "Taylor (2019)", "Brown (2022)"),
#'   Metric1  = c(3, 4, 2, 5),
#'   Metric2  = c(1, 2, 5, 4),
#'   stringsAsFactors = FALSE
#' )
#' compute_category_means_simple(sample_df)
compute_category_means_simple <- function(df, na.rm = TRUE) {
  df_mat <- as.matrix(df)

  # Check if Row 1 & Row 2 contain header metadata (e.g., read_excel(col_names = FALSE))
  row1_txt <- tolower(as.character(df_mat[1, ]))
  row2_txt <- tolower(as.character(df_mat[2, ]))

  is_2tier <- any(grepl("publication", row1_txt)) &&
    any(grepl("category|author", row2_txt))

  if (is_2tier) {
    # Data starts on row 3
    data_rows <- df_mat[-c(1, 2), , drop = FALSE]

    cat_col_idx    <- which(grepl("category", row2_txt))[1]
    author_col_idx <- which(grepl("author", row2_txt))[1]

    score_col_idx  <- setdiff(seq_len(ncol(df_mat)), c(cat_col_idx, author_col_idx))
    score_names    <- as.character(df_mat[2, score_col_idx])
  } else {
    # Single-row header (e.g., read_excel(skip = 1))
    data_rows      <- df_mat
    all_cols_lower <- tolower(colnames(df))

    cat_col_idx    <- which(grepl("category", all_cols_lower))[1]
    author_col_idx <- which(grepl("author|publication", all_cols_lower))[1]

    score_col_idx  <- setdiff(seq_len(ncol(df)), c(cat_col_idx, author_col_idx))
    score_names    <- colnames(df)[score_col_idx]
  }

  raw_cat <- trimws(as.character(data_rows[, cat_col_idx]))
  authors <- trimws(as.character(data_rows[, author_col_idx]))

  # Carry category labels down across all observations in merged blocks
  cat_filled  <- character(length(raw_cat))
  current_cat <- NA_character_

  for (i in seq_along(raw_cat)) {
    val <- raw_cat[i]
    if (!is.na(val) && val != "" && val != "NA") {
      current_cat <- val
    }
    cat_filled[i] <- current_cat
  }

  score_mat <- apply(data_rows[, score_col_idx, drop = FALSE], 2, function(col) {
    suppressWarnings(as.numeric(trimws(col)))
  })
  colnames(score_mat) <- trimws(score_names)

  # Filter out empty or unpopulated rows
  valid_rows <- !is.na(cat_filled) & cat_filled != "" &
    !is.na(authors) & authors != "" & authors != "NA" &
    apply(score_mat, 1, function(r) any(!is.na(r)))

  cat_clean   <- cat_filled[valid_rows]
  score_clean <- score_mat[valid_rows, , drop = FALSE]

  unique_cats <- unique(cat_clean)

  out_rows <- lapply(unique_cats, function(cat_name) {
    row_idx <- which(cat_clean == cat_name)

    author_count <- length(row_idx)
    item_means <- colMeans(score_clean[row_idx, , drop = FALSE], na.rm = na.rm)

    c(Author = author_count, round(item_means, 3))
  })

  result_matrix <- do.call(rbind, out_rows)
  rownames(result_matrix) <- unique_cats

  return(result_matrix)
}


#' Compute pairwise Euclidean distances and plot category heatmap
#'
#' Calculates pairwise Euclidean distances between category mean profiles and
#' produces a base R heatmap matrix with dynamic margin adjustment, contrasting
#' numerical text overlays, and orientation from top-left to bottom-right.
#'
#' @param category_means Numeric matrix or data frame of category summary means,
#'   where rows represent distinct categories and columns represent dimensions.
#' @param title Character string. Main plot header.
#' @param show_values Logical. Whether to overlay numeric distance coefficients on
#'   each heatmap tile.
#' @param palette Character vector of hex colors defining the gradient from low
#'   distance (similar) to high distance (dissimilar).
#'
#' @return Invisibly returns the symmetric pairwise Euclidean distance matrix
#'   of class `matrix`.
#' @export
#'
#' @examples
#' set.seed(42)
#' cat_profiles <- matrix(
#'   runif(12, 1, 5),
#'   nrow = 4,
#'   dimnames = list(paste("Category", LETTERS[1:4]), paste("Dim", 1:3))
#' )
#' compute_category_euclidean_distance_simple(cat_profiles)
compute_category_euclidean_distance_simple <- function(category_means,
                                                       title = "Euclidean Distance Between Publication Categories",
                                                       show_values = TRUE,
                                                       palette = c("#2166ac", "#f7f7f7", "#b2182b")) {
  if (!is.matrix(category_means) && !is.data.frame(category_means)) {
    stop("Input `category_means` must be a numeric matrix or data frame.")
  }

  mean_mat <- data.matrix(category_means)
  cats     <- rownames(mean_mat)
  n_cats   <- nrow(mean_mat)

  if (is.null(cats) || any(cats == "")) {
    cats <- paste("Category", seq_len(n_cats))
    rownames(mean_mat) <- cats
  }

  dist_obj <- stats::dist(mean_mat, method = "euclidean")
  dist_mat <- as.matrix(dist_obj)
  dimnames(dist_mat) <- list(cats, cats)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Convert maximum label string length in inches to margin lines
  max_label_len <- max(graphics::strwidth(cats, units = "inches", cex = 0.85))
  margin_needed <- (max_label_len / par("cin")[2]) + 1.5

  bottom_mar <- max(5, margin_needed)
  left_mar   <- max(5, margin_needed)

  par(mar = c(bottom_mar, left_mar, 4, 3))

  color_func <- grDevices::colorRampPalette(palette)
  n_colors   <- 100
  colors     <- color_func(n_colors)

  # Transpose and reverse rows so graphics::image() renders [1, 1] at top-left
  t_mat <- t(dist_mat[n_cats:1, ])

  graphics::image(
    x = seq_len(n_cats),
    y = seq_len(n_cats),
    z = t_mat,
    col = colors,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title
  )

  graphics::axis(1, at = seq_len(n_cats), labels = cats, las = 2, cex.axis = 0.85)
  graphics::axis(2, at = seq_len(n_cats), labels = rev(cats), las = 2, cex.axis = 0.85)
  graphics::box()

  if (show_values) {
    grid_coords <- expand.grid(x = seq_len(n_cats), y = seq_len(n_cats))
    for (idx in seq_len(nrow(grid_coords))) {
      x_pos <- grid_coords$x[idx]
      y_pos <- grid_coords$y[idx]
      val   <- t_mat[x_pos, y_pos]

      # Shift text color for contrast against saturated tiles
      text_col <- if (val > (max(dist_mat) * 0.75)) "white" else "black"
      graphics::text(x_pos, y_pos, labels = sprintf("%.2f", val), col = text_col, cex = 0.85)
    }
  }

  invisible(dist_mat)
}