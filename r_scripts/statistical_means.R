compute_category_means_simple <- function(df, na.rm = TRUE) {
  # -------------------------------------------------------------------------
  # 1. Parse Headers & Structure (Handles 2-Tier Merged Headers)
  # -------------------------------------------------------------------------
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
  
  # -------------------------------------------------------------------------
  # 2. Forward-Fill Merged Category Cells Across Papers
  # -------------------------------------------------------------------------
  raw_cat <- trimws(as.character(data_rows[, cat_col_idx]))
  authors <- trimws(as.character(data_rows[, author_col_idx]))
  
  # Carry the category label down across all authors in that merged group
  cat_filled  <- character(length(raw_cat))
  current_cat <- NA_character_
  
  for (i in seq_along(raw_cat)) {
    val <- raw_cat[i]
    if (!is.na(val) && val != "" && val != "NA") {
      current_cat <- val
    }
    cat_filled[i] <- current_cat
  }
  
  # -------------------------------------------------------------------------
  # 3. Numeric Scorecard Extraction & Row Validation
  # -------------------------------------------------------------------------
  score_mat <- apply(data_rows[, score_col_idx, drop = FALSE], 2, function(col) {
    suppressWarnings(as.numeric(trimws(col)))
  })
  colnames(score_mat) <- trimws(score_names)
  
  # Filter out empty rows (keeps rows that have an author and valid scores)
  valid_rows <- !is.na(cat_filled) & cat_filled != "" & 
    !is.na(authors) & authors != "" & authors != "NA" &
    apply(score_mat, 1, function(r) any(!is.na(r)))
  
  cat_clean   <- cat_filled[valid_rows]
  score_clean <- score_mat[valid_rows, , drop = FALSE]
  
  # -------------------------------------------------------------------------
  # 4. Compute Author Counts and Decimals Means per Category
  # -------------------------------------------------------------------------
  unique_cats <- unique(cat_clean)
  
  out_rows <- lapply(unique_cats, function(cat_name) {
    row_idx <- which(cat_clean == cat_name)
    
    # Total count of papers/authors in this category
    author_count <- length(row_idx)
    
    # Decimal means across all papers in this category
    item_means <- colMeans(score_clean[row_idx, , drop = FALSE], na.rm = na.rm)
    
    c(Author = author_count, round(item_means, 3))
  })
  
  result_matrix <- do.call(rbind, out_rows)
  rownames(result_matrix) <- unique_cats
  
  return(result_matrix)
}

compute_category_euclidean_distance_simple <- function(category_means, 
                                                       title = "Euclidean Distance Between Publication Categories",
                                                       show_values = TRUE,
                                                       palette = c("#2166ac", "#f7f7f7", "#b2182b")) {
  # 1. Input Validation
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
  
  # 2. Pairwise Euclidean Distance Calculation
  # dist() operates across rows (categories)
  dist_obj <- stats::dist(mean_mat, method = "euclidean")
  dist_mat <- as.matrix(dist_obj)
  dimnames(dist_mat) <- list(cats, cats)
  
  # 3. Render Publication-Quality Heatmap (Base R graphics, no external dependencies)
  # Set up margins for category labels
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Expand bottom and left margins to accommodate category labels
  par(mar = c(6, 6, 4, 3))
  
  # Color palette: low distance (similar) to high distance (dissimilar)
  color_func <- grDevices::colorRampPalette(palette)
  n_colors   <- 100
  colors     <- color_func(n_colors)
  
  # Base R image() maps [x, y] onto [0, 1] grid, displaying bottom-to-top
  # We transpose and flip rows so [1,1] aligns with top-left
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
  
  # Add Category Labels
  graphics::axis(1, at = seq_len(n_cats), labels = cats, las = 2, cex.axis = 0.85)
  graphics::axis(2, at = seq_len(n_cats), labels = rev(cats), las = 2, cex.axis = 0.85)
  graphics::box()
  
  # 4. Optional: Print Numeric Distance Values in Heatmap Cells
  if (show_values) {
    grid_coords <- expand.grid(x = seq_len(n_cats), y = seq_len(n_cats))
    for (idx in seq_len(nrow(grid_coords))) {
      x_pos <- grid_coords$x[idx]
      y_pos <- grid_coords$y[idx]
      val   <- t_mat[x_pos, y_pos]
      
      # Flip text color for high-contrast readability against dark backgrounds
      text_col <- if (val > (max(dist_mat) * 0.75)) "white" else "black"
      graphics::text(x_pos, y_pos, labels = sprintf("%.2f", val), col = text_col, cex = 0.85)
    }
  }
  
  # Return distance matrix for downstream tables or tests
  invisible(dist_mat)
}