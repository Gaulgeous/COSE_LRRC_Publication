# 3D and 2D Cluster Visualizations ----

#' Plot observations in 3D principal component space partitioned by cluster
#'
#' Projects multivariate observations across their first three principal component
#' dimensions using a 3D scatterplot and places a dedicated categorical legend
#' in an adjacent canvas panel.
#'
#' @param scores Matrix or data frame containing at least 3 numeric columns representing
#'   component projections.
#' @param clusters Vector or factor indicating cluster membership for each observation.
#' @param comp_names Character vector of at least length 3 specifying axis labels for the
#'   first three components. Defaults to existing column names or "Component 1:3".
#' @param cluster_names Character vector defining ordered labels for unique clusters.
#'   Inferred from `clusters` if omitted.
#' @param palette Named or unnamed character vector of hex colors for cluster categories.
#' @param title Character string. Main plot header.
#' @param angle Numeric. Azimuthal viewing angle in degrees for the 3D perspective.
#' @param layout_widths Numeric vector of length 2 defining the relative layout widths
#'   of the scatterplot and legend panels.
#' @param pch Point symbol code passed to [scatterplot3d::scatterplot3d()].
#' @param cex_symbols Point size expansion factor.
#' @param ... Additional graphical arguments forwarded to [scatterplot3d::scatterplot3d()].
#'
#' @return Invisibly returns the `scatterplot3d` plotting object, enabling downstream
#'   coordinate conversions (e.g., via `$xyz.convert()`) and point labeling.
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_scores <- matrix(rnorm(60), ncol = 3, dimnames = list(NULL, c("PC1", "PC2", "PC3")))
#' sim_clusters <- sample(c("Cluster A", "Cluster B", "Cluster C"), 20, replace = TRUE)
#' plot_pca_clusters_3d(sim_scores, sim_clusters, angle = 45)
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

  if (!is.null(comp_names)) {
    colnames(scores)[1:3] <- comp_names[1:3]
  } else if (is.null(colnames(scores))) {
    colnames(scores)[1:3] <- paste("Component", 1:3)
  }

  if (is.factor(clusters)) {
    raw_clusters <- as.character(clusters)
  } else {
    raw_clusters <- clusters
  }

  if (is.null(cluster_names)) {
    unique_vals <- unique(raw_clusters)
    cluster_names <- if (is.numeric(unique_vals)) {
      paste("Cluster", sort(unique_vals))
    } else {
      as.character(sort(unique_vals))
    }
  }

  if (is.numeric(raw_clusters)) {
    assigned_clusters <- factor(cluster_names[raw_clusters], levels = cluster_names)
  } else if (!all(raw_clusters %in% cluster_names)) {
    # Extract embedded numbers if strings like "Cluster 1" were passed
    c_ids <- as.numeric(gsub("[^0-9]", "", as.character(raw_clusters)))
    assigned_clusters <- factor(cluster_names[c_ids], levels = cluster_names)
  } else {
    assigned_clusters <- factor(raw_clusters, levels = cluster_names)
  }

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

  invisible(s3d)
}


#' Plot observations on a 2D PCA space grouped by cluster
#'
#' Renders a 2D scatterplot of observations across their leading two principal component
#' dimensions, accompanied by grid reference guides and a dedicated legend side panel.
#'
#' @param scores Matrix or data frame containing observation scores. First 2 columns
#'   are projected onto the Cartesian axes.
#' @param clusters Vector or factor indicating cluster membership for each observation.
#' @param comp_names Character vector of length 2 specifying axis labels.
#' @param cluster_names Named character vector or vector mapping cluster levels to human-readable labels.
#' @param palette Character vector of custom hex colors matching cluster levels.
#' @param title Character string. Main plot header.
#' @param layout_widths Numeric vector of length 2 defining the relative widths of the
#'   plot canvas and the legend panel.
#' @param pch Point symbol code passed to [graphics::plot()].
#' @param cex_symbols Point size expansion factor.
#' @param cex_lab Character expansion factor for axis title text.
#' @param cex_axis Character expansion factor for axis tick labels.
#' @param cex_main Character expansion factor for the main plot title.
#' @param cex_legend Character expansion factor for legend text entries.
#' @param ... Additional graphical arguments forwarded to [graphics::plot()].
#'
#' @return Invisibly returns a list with components:
#'   \item{scores}{Matrix of coordinates for the first two components.}
#'   \item{clusters}{Factor vector of assigned cluster memberships.}
#'   \item{palette}{Named character vector of hex colors mapped to cluster levels.}
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_scores <- matrix(rnorm(40), ncol = 2, dimnames = list(NULL, c("Dim1", "Dim2")))
#' sim_clusters <- factor(sample(c("Group 1", "Group 2"), 20, replace = TRUE))
#' plot_pca_clusters_2d(sim_scores, sim_clusters)
plot_pca_clusters_2d <- function(scores,
                                 clusters,
                                 comp_names = NULL,
                                 cluster_names = NULL,
                                 palette = NULL,
                                 title = "Projecting Publications onto Principal Components",
                                 layout_widths = c(4.8, 1.2),
                                 pch = 19,
                                 cex_symbols = 1.3,
                                 cex_lab = 1.2,
                                 cex_axis = 1.05,
                                 cex_main = 1.25,
                                 cex_legend = 0.95,
                                 ...) {
  scores <- as.matrix(scores)
  if (ncol(scores) < 2) {
    stop("`scores` must contain at least 2 dimensions for 2D projection.")
  }
  
  # Component column labeling
  if (!is.null(comp_names)) {
    colnames(scores)[1:2] <- comp_names[1:2]
  } else if (is.null(colnames(scores))) {
    colnames(scores)[1:2] <- paste("Component", 1:2)
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

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  })

  graphics::layout(matrix(c(1, 2), nrow = 1), widths = layout_widths)

  graphics::par(mar = c(5, 5, 4, 1), font.lab = 1, font.main = 2)
  graphics::plot(
    scores[, 1], scores[, 2],
    col      = pt_colors,
    pch      = pch,
    cex      = cex_symbols,
    xlab     = colnames(scores)[1],
    ylab     = colnames(scores)[2],
    main     = title,
    cex.lab  = cex_lab,
    cex.axis = cex_axis,
    cex.main = cex_main,
    ...
  )
  graphics::grid(col = "grey85", lty = "dotted")
  graphics::abline(h = 0, v = 0, lty = 2, col = "grey60", lwd = 0.8)

  # Redraw points over the grid lines for crisp visibility
  graphics::points(
    scores[, 1], scores[, 2],
    col = pt_colors,
    pch = pch,
    cex = cex_symbols
  )

  # Dedicated legend canvas with tight margins
  graphics::par(mar = c(5, 0, 4, 0.5))
  graphics::plot.new()

  graphics::legend(
    "center",
    legend    = names(palette),
    col       = palette,
    pch       = pch,
    pt.cex    = cex_symbols * 0.9,
    title     = "Cluster Category",
    bty       = "o",
    box.col   = "#444444",
    box.lwd   = 1.0,
    bg        = "white",
    cex       = cex_legend,
    y.intersp = 1.35,
    inset     = c(0, 0)
  )

  invisible(list(scores = scores[, 1:2], clusters = assigned_clusters, palette = palette))
}


# Cluster Distribution Donut Charts ----

#' Plot cluster composition as a donut chart
#'
#' Generates a publication-ready polar donut chart depicting cluster frequencies,
#' slice percentages, and counts using `ggplot2`.
#'
#' @param clusters Vector, factor, data frame, or results list containing cluster
#'   assignments. If a list or data frame is supplied, attempts to extract assignments
#'   from `$assignments`, `$Cluster_Name`, or `$Cluster`.
#' @param cluster_order Character vector defining the display factor level ordering.
#'   Defaults to factor levels or unique values.
#' @param palette Named or unnamed character vector of hex colors matching `cluster_order`.
#' @param title Character string. Main plot header.
#' @param label_size Numeric. Size of the text labels inside the donut slices.
#' @param hole_ratio Numeric. Ratio of the inner hole radius to the outer ring radius
#'   (between 0 and 1).
#'
#' @return Invisibly returns a named list with components:
#'   \item{plot}{The generated [ggplot2::ggplot()] object.}
#'   \item{cluster_counts}{Data frame containing counts and percentages per cluster.}
#' @export
#'
#' @examples
#' clus <- factor(c("Cluster 1", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 2"))
#' plot_cluster_donut(clus, title = "Sample Distribution")
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

  if (is.null(cluster_order)) {
    cluster_order <- if (is.factor(raw_vec)) levels(raw_vec) else unique(as.character(raw_vec))
  }

  counts_table <- table(factor(as.character(raw_vec), levels = cluster_order))

  plot_df <- data.frame(
    Cluster = factor(cluster_order, levels = cluster_order),
    Count   = as.numeric(counts_table[cluster_order]),
    stringsAsFactors = FALSE
  )

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

  r_outer <- 3
  r_inner <- r_outer * hole_ratio
  r_mid   <- (r_inner + r_outer) / 2

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


# Dimension Retention and Factor Analysis Plots ----

#' Plot empirical factor retention against Kaiser-Guttman cutoff
#'
#' Visualizes empirical eigenvalues across principal components or factors using a
#' bar chart, highlighting eigenvalues retained above a specified threshold.
#'
#' @param eigenvalues Numeric vector of eigenvalues, or an analysis list containing
#'   components `$empirical_evals`, `$eigenvalues`, or `$pc.values`.
#' @param comp_names Character vector specifying component labels for the x-axis.
#' @param top_n Integer. Maximum number of leading components to display.
#' @param cutoff Numeric. Kaiser-Guttman retention benchmark line (defaults to 1.0).
#' @param palette Named or unnamed character vector of hex colors for the bars.
#' @param title Character string. Main plot header.
#' @param subtitle Optional character string. Secondary header text.
#' @param show_cutoff Logical. If `TRUE`, draws a horizontal dashed cutoff benchmark line.
#'
#' @return Invisibly returns a named list with components:
#'   \item{plot}{The generated [ggplot2::ggplot()] bar plot.}
#'   \item{data}{Data frame containing the components and mapped eigenvalue metrics.}
#' @export
#'
#' @examples
#' sample_evals <- c(2.84, 1.42, 0.91, 0.65, 0.38)
#' plot_eigenvalue_retention(sample_evals, top_n = 4, cutoff = 1.0)
plot_eigenvalue_retention <- function(eigenvalues,
                                      comp_names = NULL,
                                      top_n = 3,
                                      cutoff = 1.0,
                                      palette = NULL,
                                      title = "Empirical Factor Retention",
                                      subtitle = NULL,
                                      show_cutoff = TRUE) {
  vals <- if (is.list(eigenvalues)) {
    if (!is.null(eigenvalues$empirical_evals)) {
      eigenvalues$empirical_evals
    } else if (!is.null(eigenvalues$eigenvalues)) {
      eigenvalues$eigenvalues
    } else if (!is.null(eigenvalues$pc.values)) {
      eigenvalues$pc.values
    } else {
      stop("Could not locate eigenvalue vector in the supplied list.")
    }
  } else if (is.numeric(eigenvalues)) {
    eigenvalues
  } else {
    stop("`eigenvalues` must be a numeric vector or a results object from scree/parallel analysis.")
  }

  if (!is.null(top_n)) {
    vals <- vals[seq_len(min(top_n, length(vals)))]
  }
  n_comp <- length(vals)

  if (is.null(comp_names)) {
    comp_names <- if (!is.null(names(vals))) names(vals) else paste("Component", seq_len(n_comp))
  } else {
    comp_names <- comp_names[seq_len(n_comp)]
  }

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

  df_eigen <- data.frame(
    Component  = factor(comp_names, levels = comp_names),
    Eigenvalue = as.numeric(vals),
    stringsAsFactors = FALSE
  )

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


# Correlation Heatmaps ----

#' Plot a triangular or full correlation heatmap
#'
#' Produces a formatted correlation tile matrix with divergent gradient colors,
#' optional numerical labels, and adaptive label contrast for dark background tiles.
#'
#' @param cor_mat Square numeric correlation matrix or object coercible to matrix.
#' @param type Character string indicating whether to display `"lower"`, `"upper"`,
#'   or `"full"` matrix portions.
#' @param palette Character vector of 3 hex colors defining the low, mid, and high
#'   color points of the gradient.
#' @param show_values Logical. Whether to overlay numeric correlation coefficients.
#' @param digits Integer. Number of decimal places when formatting correlation numbers.
#' @param title Character string. Main plot header.
#'
#' @return A [ggplot2::ggplot()] object representing the heatmap.
#' @export
#'
#' @examples
#' data(mtcars)
#' c_mat <- cor(mtcars[, 1:4])
#' plot_correlation_heatmap(c_mat, type = "lower")
plot_correlation_heatmap <- function(cor_mat,
                                     type = c("lower", "full", "upper"),
                                     palette = c("#0066FF", "#FFFFFF", "#FF2A54"),
                                     show_values = TRUE,
                                     digits = 2,
                                     title = "Correlation Heatmap") {
  type <- match.arg(type)
  cor_mat <- as.matrix(cor_mat)

  plot_mat <- cor_mat
  if (type == "lower") {
    plot_mat[upper.tri(plot_mat)] <- NA
  } else if (type == "upper") {
    plot_mat[lower.tri(plot_mat)] <- NA
  }

  df <- as.data.frame(as.table(plot_mat))
  colnames(df) <- c("Var1", "Var2", "Correlation")
  df <- df[!is.na(df$Correlation), ]

  # Orient axes so diagonal runs top-left to bottom-right
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


# Publication and Category Integration ----

#' Map cluster assignments to raw multi-row publication table
#'
#' Parses a multi-row rectangular table containing publication categories and authors,
#' forward-fills merged category headers, validates length against cluster assignments,
#' and returns a standardized data frame.
#'
#' @param raw_data Data frame or matrix matching the raw table structure where column 1
#'   is Category and column 2 is Author or Publication.
#' @param cluster_assignments Vector or factor of cluster labels corresponding to the
#'   extracted publications.
#'
#' @return A data frame containing:
#'   \item{Category}{Character vector of forward-filled category labels.}
#'   \item{Author}{Character vector of author/publication names.}
#'   \item{Assigned_Cluster}{Factor of standardized cluster memberships.}
#' @export
#'
#' @examples
#' sample_tbl <- data.frame(
#'   Col1 = c("Methodology", "", "Applications"),
#'   Col2 = c("Smith (2020)", "Jones (2021)", "Taylor (2022)"),
#'   stringsAsFactors = FALSE
#' )
#' assign_clusters_to_publications(sample_tbl, c("1", "2", "1"))
assign_clusters_to_publications <- function(raw_data, cluster_assignments) {
  raw_mat <- as.matrix(raw_data)

  start_row <- 1
  for (r in seq_len(min(5, nrow(raw_mat)))) {
    val1 <- tolower(trimws(as.character(raw_mat[r, 1])))
    val2 <- tolower(trimws(as.character(raw_mat[r, 2])))
    if (val1 %in% c("publication", "category") || val2 %in% c("author", "publication")) {
      start_row <- r + 1
    }
  }

  clean_mat <- raw_mat[start_row:nrow(raw_mat), 1:2, drop = FALSE]

  cats  <- clean_mat[, 1]
  auths <- trimws(as.character(clean_mat[, 2]))

  # Discard empty trailing or padding rows
  valid_rows <- which(!is.na(auths) & auths != "" & auths != "NA")
  cats  <- cats[valid_rows]
  auths <- auths[valid_rows]

  # Forward-fill merged category cells
  filled_cats <- character(length(cats))
  current_cat <- NA_character_

  for (i in seq_along(cats)) {
    val <- cats[i]
    if (!is.na(val) && trimws(val) != "" && trimws(val) != "NA") {
      current_cat <- trimws(val)
    }
    filled_cats[i] <- current_cat
  }

  if (length(cluster_assignments) != length(auths)) {
    stop(sprintf(
      "Length mismatch: Data table has %d publications, but %d cluster assignments were passed.",
      length(auths), length(cluster_assignments)
    ))
  }

  cluster_vec <- as.character(cluster_assignments)
  if (!all(grepl("^Cluster", cluster_vec, ignore.case = TRUE))) {
    cluster_vec <- paste("Cluster", cluster_vec)
  }
  cluster_factor <- factor(cluster_vec, levels = sort(unique(cluster_vec)))

  data.frame(
    Category         = filled_cats,
    Author           = auths,
    Assigned_Cluster = cluster_factor,
    stringsAsFactors = FALSE
  )
}


#' Append cleaned categories to cluster assignments data frame
#'
#' Merges cleaned category labels with an existing cluster assignment table and produces
#' category-by-cluster contingency counts.
#'
#' @param assignments_df Data frame containing observation cluster assignments.
#' @param categories Character vector or factor of category labels with length matching
#'   `nrow(assignments_df)`.
#'
#' @return A named list containing:
#'   \item{augmented_assignments}{Data frame with the `Category` factor column appended.}
#'   \item{unique_categories}{Character vector of unique sorted category levels.}
#'   \item{contingency_table}{Contingency [table()] cross-tabulating clusters and categories.}
#' @export
#'
#' @examples
#' assign_df <- data.frame(Case = 1:3, Cluster_Name = c("C1", "C2", "C1"))
#' cats <- c("Empirical", "Theoretical", "Empirical")
#' prepare_cluster_category_data(assign_df, cats)
prepare_cluster_category_data <- function(assignments_df, categories) {
  if (!is.data.frame(assignments_df)) {
    assignments_df <- as.data.frame(assignments_df, stringsAsFactors = FALSE)
  }

  categories <- as.character(categories)

  if (nrow(assignments_df) != length(categories)) {
    stop(sprintf(
      "Length mismatch: assignments has %d rows, but categories has length %d.",
      nrow(assignments_df), length(categories)
    ))
  }

  augmented_df <- assignments_df
  augmented_df$Category <- factor(trimws(categories))

  unique_labels <- levels(augmented_df$Category)

  clust_col <- if ("Cluster_Name" %in% colnames(augmented_df)) {
    "Cluster_Name"
  } else if ("Cluster" %in% colnames(augmented_df)) {
    "Cluster"
  } else {
    "Cluster_Num"
  }

  summary_tbl <- table(
    Cluster  = augmented_df[[clust_col]],
    Category = augmented_df$Category
  )

  list(
    augmented_assignments = augmented_df,
    unique_categories     = unique_labels,
    contingency_table     = summary_tbl
  )
}


#' Plot faceted donut charts of category distribution across clusters
#'
#' Produces multi-panel donut charts showing the proportional breakdown of categories
#' within each cluster using `ggplot2` facetting.
#'
#' @param assignments_df Data frame containing a `Category` column alongside cluster
#'   grouping identifiers.
#' @param cluster_col Character string specifying the cluster grouping column name.
#'   Defaults to `"Cluster_Name"` (falling back to `"Cluster"` or `"Cluster_Num"`).
#' @param category_palette Optional named character vector of hex colors mapped to category levels.
#' @param hole_ratio Numeric. Relative size of the central hole (default 0.55).
#' @param title Character string. Main plot header.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#'
#' @examples
#' df_test <- data.frame(
#'   Case = 1:6,
#'   Cluster_Name = rep(c("Cluster 1", "Cluster 2"), each = 3),
#'   Category = c("Cat A", "Cat B", "Cat A", "Cat B", "Cat B", "Cat C")
#' )
#' plot_cluster_category_donuts(df_test)
plot_cluster_category_donuts <- function(assignments_df,
                                         cluster_col = "Cluster_Name",
                                         category_palette = NULL,
                                         hole_ratio = 0.55,
                                         title = "Publication Categories Distributed Across Clusters") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required. Install it using install.packages('ggplot2').")
  }

  if (!cluster_col %in% colnames(assignments_df)) {
    if ("Cluster" %in% colnames(assignments_df)) {
      cluster_col <- "Cluster"
    } else {
      cluster_col <- "Cluster_Num"
    }
  }

  if (!"Category" %in% colnames(assignments_df)) {
    stop("`assignments_df` must contain a 'Category' column. Run `prepare_cluster_category_data` first.")
  }

  df <- assignments_df
  df$Target_Cluster <- as.factor(df[[cluster_col]])
  df$Category       <- as.factor(df$Category)

  comp_agg <- stats::aggregate(
    Case ~ Target_Cluster + Category,
    data = df,
    FUN = length
  )
  names(comp_agg)[names(comp_agg) == "Case"] <- "Count"

  comp_agg <- do.call(rbind, lapply(split(comp_agg, comp_agg$Target_Cluster), function(sub) {
    sub <- sub[order(sub$Count, decreasing = TRUE), ]
    total <- sum(sub$Count)
    sub$Fraction   <- sub$Count / total
    sub$Percentage <- round(sub$Fraction * 100, 1)
    sub$Ymax       <- cumsum(sub$Fraction)
    sub$Ymin       <- c(0, head(sub$Ymax, n = -1))
    sub$Ymid       <- (sub$Ymin + sub$Ymax) / 2
    sub$Slice_Label <- paste0(sub$Count, "\n(", sub$Percentage, "%)")
    sub
  }))

  unique_cats <- levels(df$Category)
  if (is.null(category_palette)) {
    base_palette <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                      "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                      "#9C755F", "#BAB0AC", "#D37295", "#FABFD2")
    category_palette <- rep(base_palette, length.out = length(unique_cats))
    names(category_palette) <- unique_cats
  }

  p <- ggplot2::ggplot(comp_agg, ggplot2::aes(
    ymin = Ymin, ymax = Ymax,
    xmin = hole_ratio + 1, xmax = 2 + 1,
    fill = Category
  )) +
    ggplot2::geom_rect(color = "white", linewidth = 0.8) +
    # Display count and percentage on slices >= 4% of cluster total
    ggplot2::geom_text(
      data = subset(comp_agg, Fraction >= 0.04),
      ggplot2::aes(x = (hole_ratio + 1 + 3) / 2, y = Ymid, label = Slice_Label),
      color = "white",
      fontface = "bold",
      size = 3.2,
      lineheight = 0.85
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(0.5, 3.2)) +
    ggplot2::facet_wrap(~ Target_Cluster, nrow = 1) +
    ggplot2::scale_fill_manual(values = category_palette, drop = FALSE) +
    ggplot2::labs(
      title    = title,
      fill     = "Publication Category"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5, size = 15, margin = ggplot2::margin(b = 4)),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey40", margin = ggplot2::margin(b = 15)),
      strip.text       = ggplot2::element_text(face = "bold", size = 12, margin = ggplot2::margin(b = 8)),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", size = 10),
      legend.text      = ggplot2::element_text(size = 9),
      plot.margin      = ggplot2::margin(15, 15, 15, 15)
    )

  print(p)
  return(p)
}


# Pipeline and Batch Visualization Helpers ----

#' Animate and evaluate 3D principal component perspectives
#'
#' Iterates through viewing angles to inspect cluster separation in 3D principal component space
#' and labels individual observations at their converted spatial coordinates.
#'
#' @param scores Matrix or data frame containing at least 3 principal component columns.
#' @param clusters Vector or factor of cluster memberships.
#' @param comp_names Character vector specifying axis names for the components.
#' @param cluster_names Character vector defining ordered labels for clusters.
#' @param palette Character vector of colors corresponding to cluster levels.
#' @param ... Additional arguments passed to [plot_pca_clusters_3d()].
#'
#' @return Invisibly returns the final `scatterplot3d` plot object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_3d_images(
#'   scores        = scores_r3,
#'   clusters      = cluster_results$assignments$Cluster,
#'   comp_names    = c("PC1", "PC2", "PC3"),
#'   cluster_names = c("Cluster 1", "Cluster 2"),
#'   palette       = c("#0066FF", "#FF2A54")
#' )
#' }
plot_3d_images <- function(scores,
                           clusters,
                           comp_names,
                           cluster_names,
                           palette,
                           ...) {
  s3d_plot <- plot_pca_clusters_3d(
    scores        = scores_r3,
    clusters      = cluster_results$assignments$Cluster,
    comp_names    = comp_names,
    cluster_names = cluster_names,
    palette       = palette,
    title         = "Projecting Publications onto Principal Components",
    angle         = 130 # Optimal cluster separation azimuth
  )

  # Rotate through angles to assess visual separation across perspectives
  for (theta in seq(0, 360, by = 10)) {
    s3d_plot <- plot_pca_clusters_3d(
      scores        = scores_r3,
      clusters      = cluster_results$assignments$Cluster,
      comp_names    = comp_names,
      cluster_names = cluster_names,
      palette       = palette,
      title         = sprintf("Projecting Publications onto Principal Components (Angle: %d°)", theta),
      angle         = theta
    )

    Sys.sleep(0.1)
  }

  text(
    s3d_plot$xyz.convert(scores_r3[, 1], scores_r3[, 2], scores_r3[, 3]),
    labels = rownames(scores_r3),
    pos    = 3,
    cex    = 0.65
  )

  invisible(s3d_plot)
}


#' Generate supporting diagnostic plots for cluster and dimension retention
#'
#' Executes a sequence of visualization diagnostics: augments cluster assignments
#' with category metadata, renders category-by-cluster faceted donut charts,
#' plots overall cluster counts, and displays empirical eigenvalue retention bars.
#'
#' @param cluster_results A list containing clustering results, including an
#'   `assignments` data frame.
#' @param categories_clean Character vector or factor of cleaned publication categories
#'   matching observation rows.
#' @param cluster_names Character vector of cluster level labels in desired ordering.
#' @param palette Named or unnamed character vector of hex colors mapped to clusters.
#' @param pa_res Results list or numeric vector from parallel/scree analysis passed
#'   to [plot_eigenvalue_retention()].
#' @param comp_names Character vector specifying component labels.
#' @param manual_n_components Integer. Number of leading components to display.
#' @param ... Additional arguments passed to downstream plotting routines.
#'
#' @return No return value, called for plotting side effects.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_accompanying_evidence(
#'   cluster_results     = cluster_results,
#'   categories_clean    = categories_clean,
#'   cluster_names       = cluster_names,
#'   palette             = palette,
#'   pa_res              = pa_res,
#'   comp_names          = comp_names,
#'   manual_n_components = 2
#' )
#' }
plot_accompanying_evidence <- function(cluster_results,
                                       categories_clean,
                                       cluster_names,
                                       palette,
                                       pa_res,
                                       comp_names,
                                       manual_n_components,
                                       ...) {
  prepared_data <- prepare_cluster_category_data(
    assignments_df = cluster_results$assignments,
    categories     = categories_clean
  )

  cluster_results$assignments <- prepared_data$augmented_assignments

  plot_cluster_category_donuts(
    assignments_df = cluster_results$assignments,
    cluster_col    = "Cluster_Name"
  )

  donut_res <- plot_cluster_donut(
    clusters      = cluster_results,
    cluster_order = cluster_names,
    palette       = palette,
    title         = "Publication Cluster Distribution"
  )

  print(donut_res$cluster_counts)

  plot_eigenvalue_retention(
    eigenvalues = pa_res,
    comp_names  = comp_names,
    top_n       = manual_n_components,
    title       = "Horn's Parallel Analysis: Eigenvalue Retention"
  )
}
