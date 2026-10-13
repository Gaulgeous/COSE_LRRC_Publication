# Plot observations in 3D principal component space partitioned by cluster
plot_pca_clusters_3d <- function(scores,
                                 clusters,
                                 comp_names = NULL,
                                 cluster_names = NULL,
                                 palette = NULL,
                                 title = "Projection of Observations onto Principal Components",
                                 angle = 60,
                                 layout_widths = c(3.8, 1.2),
                                 pch = 19,
                                 cex_symbols = 1.3,
                                 ...) {
  scores <- as.matrix(scores)
  if (ncol(scores) < 3) {
    stop("`scores` must contain at least 3 dimensions for 3D projection.")
  }
  
  # Component column labeling
  if (!is.null(comp_names)) {
    colnames(scores)[1:3] <- comp_names[1:3]
  } else if (is.null(colnames(scores))) {
    colnames(scores)[1:3] <- paste("Component", 1:3)
  }
  
  # Standardize cluster factor levels and labels
  if (is.factor(clusters)) {
    raw_clusters <- as.character(clusters)
  } else {
    raw_clusters <- clusters
  }
  
  # Determine distinct cluster names if not explicitly supplied
  if (is.null(cluster_names)) {
    unique_vals <- unique(raw_clusters)
    cluster_names <- if (is.numeric(unique_vals)) {
      paste("Cluster", sort(unique_vals))
    } else {
      as.character(sort(unique_vals))
    }
  }
  
  # Map numeric IDs or existing strings to cluster_names
  if (is.numeric(raw_clusters)) {
    assigned_clusters <- factor(cluster_names[raw_clusters], levels = cluster_names)
  } else if (!all(raw_clusters %in% cluster_names)) {
    # Extract embedded numbers if strings like "Cluster 1" were passed
    c_ids <- as.numeric(gsub("[^0-9]", "", as.character(raw_clusters)))
    assigned_clusters <- factor(cluster_names[c_ids], levels = cluster_names)
  } else {
    assigned_clusters <- factor(raw_clusters, levels = cluster_names)
  }
  
  # Configure color palette
  if (is.null(palette)) {
    default_colors <- c("#FF2A54", "#0066FF", "#9B59B6", "#2ECC71", "#E67E22", "#1ABC9C")
    palette <- default_colors[seq_along(cluster_names)]
    names(palette) <- cluster_names
  } else if (is.null(names(palette))) {
    names(palette) <- cluster_names[seq_along(palette)]
  }
  
  pt_colors <- palette[as.character(assigned_clusters)]
  
  # Manage graphic device state safely via on.exit
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  })
  
  # Split canvas for plot and dedicated legend panel
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = layout_widths)
  
  # Panel 1: 3D Scatterplot
  graphics::par(mar = c(5, 4, 4, 1))
  s3d <- scatterplot3d::scatterplot3d(
    scores[, 1], scores[, 2], scores[, 3],
    color       = pt_colors,
    pch         = pch,
    cex.symbols = cex_symbols,
    xlab        = colnames(scores)[1],
    ylab        = colnames(scores)[2],
    zlab        = colnames(scores)[3],
    main        = title,
    angle       = angle,
    grid        = TRUE,
    box         = TRUE,
    ...
  )
  
  # Panel 2: Dedicated Legend Canvas
  graphics::par(mar = c(5, 0, 4, 1))
  graphics::plot.new()
  
  graphics::legend(
    "left",
    legend    = names(palette),
    col       = palette,
    pch       = pch,
    title     = "Cluster Category",
    bty       = "o",
    box.col   = "#444444",
    box.lwd   = 1.1,
    bg        = "white",
    cex       = 0.9,
    y.intersp = 1.3
  )
  
  # Return scatterplot3d object invisibly for downstream coordinate conversions/overlays
  invisible(s3d)
}



plot_cluster_donut <- function(clusters,
                               cluster_order = NULL,
                               palette = NULL,
                               title = "Publication Cluster Distribution",
                               label_size = 4.5,
                               hole_ratio = 0.7) {
  # Extract vector of cluster labels
  if (is.list(clusters) && !is.data.frame(clusters)) {
    if ("assignments" %in% names(clusters)) {
      clusters <- clusters$assignments
    }
  }
  
  if (is.data.frame(clusters)) {
    if ("Cluster_Name" %in% names(clusters)) {
      raw_vec <- clusters$Cluster_Name
    } else if ("Cluster" %in% names(clusters)) {
      raw_vec <- clusters$Cluster
    } else {
      raw_vec <- clusters[[1]]
    }
  } else {
    raw_vec <- clusters
  }
  
  # Establish level order
  if (is.null(cluster_order)) {
    cluster_order <- if (is.factor(raw_vec)) levels(raw_vec) else unique(as.character(raw_vec))
  }
  
  # Tabulate occurrences across specified order
  counts_table <- table(factor(as.character(raw_vec), levels = cluster_order))
  
  plot_df <- data.frame(
    Cluster = factor(cluster_order, levels = cluster_order),
    Count   = as.numeric(counts_table[cluster_order]),
    stringsAsFactors = FALSE
  )
  
  # Compute exact angular boundaries and slice midpoints
  total_obs          <- sum(plot_df$Count)
  plot_df$Fraction   <- if (total_obs > 0) plot_df$Count / total_obs else 0
  plot_df$Percentage <- plot_df$Fraction * 100
  plot_df$ymax       <- cumsum(plot_df$Fraction)
  plot_df$ymin       <- c(0, head(plot_df$ymax, -1))
  plot_df$ymid       <- (plot_df$ymin + plot_df$ymax) / 2
  plot_df$Label      <- ifelse(
    plot_df$Count > 0,
    sprintf("%d\n(%.1f%%)", plot_df$Count, plot_df$Percentage),
    ""
  )
  
  # Radial dimensions
  r_outer <- 3
  r_inner <- r_outer * hole_ratio
  r_mid   <- (r_inner + r_outer) / 2
  
  # Standardize color palette
  if (is.null(palette)) {
    default_cols <- c("#9B59B6", "#FF2A54", "#0066FF", "#2ECC71", "#F39C12")
    palette <- stats::setNames(default_cols[seq_along(cluster_order)], cluster_order)
  } else if (is.null(names(palette))) {
    names(palette) <- cluster_order[seq_along(palette)]
  }
  
  # Construct plot using explicit geometric rectangles
  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        ymin = ymin,
        ymax = ymax,
        xmin = r_inner,
        xmax = r_outer,
        fill = Cluster
      ),
      color     = "white",
      linewidth = 1.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x     = r_mid,
        y     = ymid,
        label = Label
      ),
      color      = "white",
      size       = label_size,
      fontface   = "bold",
      lineheight = 0.9
    ) +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::xlim(c(0, r_outer + 0.5)) +
    ggplot2::scale_fill_manual(
      values = palette,
      breaks = cluster_order,
      drop   = FALSE
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold", margin = ggplot2::margin(b = 10)),
      legend.position = "bottom",
      legend.title    = ggplot2::element_blank(),
      legend.text     = ggplot2::element_text(size = 11, face = "bold"),
      plot.margin     = ggplot2::margin(12, 12, 12, 12)
    )
  
  print(p)
  
  invisible(list(
    plot           = p,
    cluster_counts = plot_df[, c("Cluster", "Count", "Percentage")]
  ))
}


# Render empirical factor retention bar plot against a Kaiser-Guttman cutoff
plot_eigenvalue_retention <- function(eigenvalues,
                                      comp_names = NULL,
                                      top_n = 3,
                                      cutoff = 1.0,
                                      palette = NULL,
                                      title = "Empirical Factor Retention",
                                      subtitle = NULL,
                                      show_cutoff = TRUE) {

  # Ingest eigenvalue inputs directly from vectors or upstream function outputs
  vals <- if (is.list(eigenvalues)) {
    if (!is.null(eigenvalues$empirical_evals)) {
      eigenvalues$empirical_evals            # From perform_horn_parallel_analysis_exact()
    } else if (!is.null(eigenvalues$eigenvalues)) {
      eigenvalues$eigenvalues                # From compute_pca_scree_analysis() or PCA
    } else if (!is.null(eigenvalues$pc.values)) {
      eigenvalues$pc.values                  # From psych::fa.parallel()
    } else {
      stop("Could not locate eigenvalue vector in the supplied list.")
    }
  } else if (is.numeric(eigenvalues)) {
    eigenvalues
  } else {
    stop("`eigenvalues` must be a numeric vector or a results object from scree/parallel analysis.")
  }
  
  # Subset to top_n components if specified
  if (!is.null(top_n)) {
    vals <- vals[seq_len(min(top_n, length(vals)))]
  }
  n_comp <- length(vals)
  
  # Label component names
  if (is.null(comp_names)) {
    comp_names <- if (!is.null(names(vals))) names(vals) else paste("Component", seq_len(n_comp))
  } else {
    comp_names <- comp_names[seq_len(n_comp)]
  }
  
  # Configure dynamic color palette
  if (is.null(palette)) {
    default_colors <- c("#0066FF", "#9B59B6", "#FF2A54", "#2ECC71", "#E67E22", "#1ABC9C")
    palette <- if (n_comp <= length(default_colors)) {
      default_colors[seq_len(n_comp)]
    } else {
      grDevices::hcl.colors(n_comp, "Dark 3")
    }
    names(palette) <- comp_names
  } else if (is.null(names(palette))) {
    names(palette) <- comp_names[seq_along(palette)]
  }
  
  # Build plotting dataframe
  df_eigen <- data.frame(
    Component  = factor(comp_names, levels = comp_names),
    Eigenvalue = as.numeric(vals),
    stringsAsFactors = FALSE
  )
  
  # Render ggplot bar chart
  p <- ggplot2::ggplot(df_eigen, ggplot2::aes(x = Component, y = Eigenvalue, fill = Component)) +
    ggplot2::geom_col(width = 0.55, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", Eigenvalue)),
      vjust    = -0.6,
      size     = 4.2,
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x        = NULL,
      y        = "Empirical Eigenvalue",
      title    = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = 10.5, color = "#555555", hjust = 0.5, margin = ggplot2::margin(b = 12)),
      axis.title.y  = ggplot2::element_text(size = 12, face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text.x   = ggplot2::element_text(size = 11, face = "bold", color = "black"),
      axis.text.y   = ggplot2::element_text(size = 10, color = "black"),
      plot.margin   = ggplot2::margin(15, 15, 15, 15)
    )
  
  # Dynamic reference cutoff line and label positioning
  if (show_cutoff && !is.null(cutoff)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = cutoff,
        linetype   = "dashed",
        color      = "#555555",
        linewidth  = 0.8
      ) +
      ggplot2::annotate(
        "text",
        x        = 0.55,
        y        = cutoff + (max(df_eigen$Eigenvalue) * 0.03),
        label    = sprintf("Kaiser Cutoff (%.1f)", cutoff),
        size     = 3.5,
        fontface = "italic",
        color    = "#555555",
        hjust    = 0
      )
  }
  
  print(p)
  
  invisible(list(
    plot = p,
    data = df_eigen
  ))
}

plot_correlation_heatmap <- function(cor_mat,
                                     type = c("lower", "full", "upper"),
                                     palette = c("#0066FF", "#FFFFFF", "#FF2A54"),
                                     show_values = TRUE,
                                     digits = 2,
                                     title = "Correlation Heatmap") {
  type <- match.arg(type)
  cor_mat <- as.matrix(cor_mat)
  
  # Mask upper or lower triangle if requested
  plot_mat <- cor_mat
  if (type == "lower") {
    plot_mat[upper.tri(plot_mat)] <- NA
  } else if (type == "upper") {
    plot_mat[lower.tri(plot_mat)] <- NA
  }
  
  # Reshape into tidy format
  df <- as.data.frame(as.table(plot_mat))
  colnames(df) <- c("Var1", "Var2", "Correlation")
  df <- df[!is.na(df$Correlation), ]
  
  # Orient axes to display diagonal from top-left to bottom-right
  df$Var1 <- factor(df$Var1, levels = rev(rownames(cor_mat)))
  df$Var2 <- factor(df$Var2, levels = colnames(cor_mat))
  
  # Dynamic text contrast for readability over dark tiles
  df$text_color <- ifelse(abs(df$Correlation) > 0.6, "white", "black")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradient2(
      low      = palette[1],
      mid      = palette[2],
      high     = palette[3],
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Pearson\nCorr"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 10)),
      axis.text.x     = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 10),
      axis.text.y     = ggplot2::element_text(face = "bold", size = 10),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right"
    )
  
  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf(paste0("%.", digits, "f"), Correlation),
        color = text_color
      ),
      size     = 3.8,
      fontface = "bold",
      show.legend = FALSE
    ) +
      ggplot2::scale_color_identity()
  }
  
  return(p)
}