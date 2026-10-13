compute_category_means_simple <- function(df, 
                                          category_col = NULL, 
                                          na.rm = TRUE) {
  # -------------------------------------------------------------------------
  # 1. Parse Multi-Column Structure
  # -------------------------------------------------------------------------
  
  # Check if df has a raw 2-tier header in rows 1 and 2 (e.g. read.csv(header=FALSE))
  has_2row_header <- any(grepl("publication", as.character(df[1, ]), ignore.case = TRUE)) &&
                     any(grepl("category", as.character(df[2, ]), ignore.case = TRUE))
  
  if (has_2row_header) {
    top_header <- as.character(df[1, ])
    sub_header <- as.character(df[2, ])
    
    # Identify publication vs non-publication columns
    pub_mask <- grepl("publication", top_header, ignore.case = TRUE)
    cat_idx  <- which(pub_mask & grepl("category", sub_header, ignore.case = TRUE))[1]
    
    cat_vector <- as.character(df[-c(1, 2), cat_idx])
    
    # Non-publication data subset (coerced to numeric)
    score_df <- df[-c(1, 2), !pub_mask, drop = FALSE]
    score_mat <- data.matrix(score_df)
    colnames(score_mat) <- sub_header[!pub_mask]
    
  } else if (is.data.frame(df$Publication) || is.matrix(df$Publication)) {
    # Nested data.frame / matrix column (e.g., df$Publication$Category)
    pub_df     <- as.data.frame(df$Publication)
    cat_col_id <- if (!is.null(category_col)) category_col else grep("category", colnames(pub_df), ignore.case = TRUE, value = TRUE)[1]
    cat_vector <- as.character(pub_df[[cat_col_id]])
    
    # All columns other than df$Publication
    non_pub_df <- df[, setdiff(names(df), "Publication"), drop = FALSE]
    score_mat  <- data.matrix(non_pub_df)
    
  } else {
    # Delimited/prefixed column names (e.g., Publication.Category or Publication_Category)
    all_cols <- colnames(df)
    pub_mask <- grepl("^publication[._/ ]?", all_cols, ignore.case = TRUE)
    
    if (!any(pub_mask)) {
      # Fallback: look for a column literally named 'Category'
      pub_mask <- grepl("^(category|author)$", all_cols, ignore.case = TRUE)
    }
    
    if (!is.null(category_col)) {
      cat_vector <- as.character(df[[category_col]])
    } else {
      cat_name   <- grep("category", all_cols[pub_mask], ignore.case = TRUE, value = TRUE)[1]
      cat_vector <- as.character(df[[cat_name]])
    }
    
    # Non-publication score columns
    non_pub_df <- df[, !pub_mask, drop = FALSE]
    score_mat  <- data.matrix(non_pub_df)
  }
  
  # -------------------------------------------------------------------------
  # 2. Iterate Over Categories & Compute Column Means
  # -------------------------------------------------------------------------
  
  unique_categories <- unique(stats::na.omit(cat_vector))
  n_cats            <- length(unique_categories)
  n_cols            <- ncol(score_mat)
  col_names         <- colnames(score_mat)
  
  # Pre-allocate output matrix: rows = categories, cols = scorecard items
  mean_mat <- matrix(
    NA_real_, 
    nrow     = n_cats, 
    ncol     = n_cols, 
    dimnames = list(unique_categories, col_names)
  )
  
  for (i in seq_along(unique_categories)) {
    cat_val <- unique_categories[i]
    row_idx <- which(cat_vector == cat_val)
    
    # Calculate column means for this category subset
    cat_subset <- score_mat[row_idx, , drop = FALSE]
    mean_mat[i, ] <- colMeans(cat_subset, na.rm = na.rm)
  }
  
  return(mean_mat)
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