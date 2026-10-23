source("helpers.R")

# PCA Score Projections ----

#' Project observations onto rotated principal component scores
#'
#' Projects standardized observation items onto rotated principal component dimensions
#' using regression factor scoring weights derived from the empirical correlation
#' matrix and PCA structure coefficients: \eqn{W = R^{-1} S}.
#'
#' @param df Data frame or matrix of raw numeric observations.
#' @param pca_res Results list from a PCA model containing a `structure_matrix` component.
#' @param comp_names Character vector of labels assigned to component score columns.
#'
#' @return A numeric matrix of projected component scores with preserved observation
#'   row names and assigned component column names.
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_df <- data.frame(
#'   Item1 = rnorm(20),
#'   Item2 = rnorm(20),
#'   Item3 = rnorm(20),
#'   row.names = paste0("Obs_", 1:20)
#' )
#' fake_pca <- list(
#'   structure_matrix = matrix(
#'     c(0.8, 0.2, 0.1, 0.1, 0.7, 0.6),
#'     nrow = 3, ncol = 2,
#'     dimnames = list(c("Item1", "Item2", "Item3"), c("Dim1", "Dim2"))
#'   )
#' )
#' project_pca_scores(sim_df, fake_pca, c("Dimension 1", "Dimension 2"))
project_pca_scores <- function(df, pca_res, comp_names) {
  clean_df <- clean_labels(df)
  clean_df <- na.omit(as.data.frame(lapply(clean_df, as.numeric)))
  rownames(clean_df) <- rownames(na.omit(df))

  Z <- scale(clean_df)

  r_mat <- ensure_psd(cor(Z))

  structure_mat <- as.matrix(pca_res$structure_matrix)

  # Standard regression score weights: W = R^(-1) %*% S
  W <- solve(r_mat) %*% structure_mat

  scores_r3 <- Z %*% W
  rownames(scores_r3) <- rownames(clean_df)

  colnames(scores_r3) <- comp_names

  # Double check that this column shouldn't be inverted
  # scores_r3[, "Data Provenance"] <- -1 * scores_r3[, "Data Provenance"]

  scores_r3
}


# Cluster Triangulation and Diagnostics ----

#' Triangulate optimal cluster count across multiple diagnostic metrics
#'
#' Evaluates candidate cluster numbers across an observation coordinate matrix using
#' three complementary criteria: normalized elbow curvature (scaled knee detection),
#' average silhouette width, and Tibshirani's gap statistic with Monte Carlo reference
#' distributions. Selects consensus cluster count by majority vote.
#'
#' @param scores_r3 Numeric matrix or data frame of spatial observation coordinates.
#' @param k_range Integer vector of candidate cluster counts to evaluate (default `2:5`).
#' @param b_gap Integer. Number of bootstrap reference datasets for the gap statistic (default 500).
#' @param method Character string specifying clustering algorithm (`"kmeans"`, `"ward.D2"`,
#'   or `"average"`).
#' @param seed Optional integer for reproducible random starts and bootstrap generation.
#' @param plot Logical. If `TRUE`, renders diagnostic multi-panel plots for elbow,
#'   silhouette, gap, and hierarchical dendrograms.
#'
#' @return Invisibly returns a named list with components:
#'   \item{scores}{Matrix of coordinates analyzed.}
#'   \item{clusters}{Named integer vector of final cluster assignments.}
#'   \item{consensus_k}{Consensus optimal cluster number determined by vote.}
#'   \item{summary_table}{Data frame summarizing WSS inertia, silhouette width, and gap values.}
#'   \item{optimal_k}{Named list of individually optimal \eqn{k} values per metric.}
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_coords <- matrix(rnorm(60), ncol = 3, dimnames = list(paste0("Case_", 1:20), NULL))
#' triangulate_r3_clusters(sim_coords, k_range = 2:4, b_gap = 20, plot = FALSE)
triangulate_r3_clusters <- function(scores_r3,
                                   k_range = 2:5,
                                   b_gap = 500,
                                   method = c("kmeans", "ward.D2", "average"),
                                   seed = 42,
                                   plot = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)

  scores_r3    <- as.matrix(scores_r3)
  n_obs        <- nrow(scores_r3)
  k_candidates <- k_range[k_range < n_obs]
  elbow_k_list <- 1:max(k_candidates)
  dist_r3      <- stats::dist(scores_r3, method = "euclidean")

  get_clusters <- function(k) {
    if (method == "kmeans") {
      stats::kmeans(scores_r3, centers = k, nstart = 25)$cluster
    } else {
      stats::cutree(stats::hclust(dist_r3, method = method), k = k)
    }
  }

  wss <- sapply(elbow_k_list, function(k) {
    stats::kmeans(scores_r3, centers = k, nstart = 25)$tot.withinss
  })

  # Min-Max Normalization to eliminate axis-scale distortion in knee detection
  norm_x <- (elbow_k_list - min(elbow_k_list)) / (max(elbow_k_list) - min(elbow_k_list))
  norm_y <- (wss - min(wss)) / (max(wss) - min(wss))

  chord_vec   <- c(norm_x[length(norm_x)] - norm_x[1], norm_y[length(norm_y)] - norm_y[1])
  chord_unit  <- chord_vec / sqrt(sum(chord_vec^2))
  diff_pts    <- cbind(norm_x - norm_x[1], norm_y - norm_y[1])
  proj_pts    <- (diff_pts %*% chord_unit) %*% t(chord_unit)
  dist_chord  <- sqrt(rowSums((diff_pts - proj_pts)^2))
  opt_k_elbow <- elbow_k_list[which.max(dist_chord)]

  sil_scores <- sapply(k_candidates, function(k) {
    cl <- get_clusters(k)
    ss <- cluster::silhouette(cl, dist_r3)
    mean(ss[, "sil_width"])
  })
  opt_k_sil <- k_candidates[which.max(sil_scores)]

  gap_cluster_fun <- if (method == "kmeans") {
    function(x, k) stats::kmeans(x, centers = k, nstart = 25)
  } else {
    function(x, k) list(cluster = stats::cutree(stats::hclust(stats::dist(x), method = method), k = k))
  }

  gap_res <- cluster::clusGap(
    scores_r3,
    FUNcluster = gap_cluster_fun,
    K.max      = max(k_candidates),
    B          = b_gap,
    spaceH0    = "original"
  )
  gap_tab   <- gap_res$Tab
  gap_vals  <- gap_tab[k_candidates, "gap"]
  gap_se    <- gap_tab[k_candidates, "SE.sim"]
  opt_k_gap <- cluster::maxSE(gap_tab[, "gap"], gap_tab[, "SE.sim"], method = "firstSEmax")

  # Consensus K (Majority Vote, defaulting to Silhouette on tie)
  k_votes     <- c(opt_k_elbow, opt_k_sil, opt_k_gap)
  consensus_k <- as.integer(names(which.max(table(k_votes))))
  if (length(unique(k_votes)) == 3) consensus_k <- opt_k_sil

  final_assignments <- get_clusters(consensus_k)
  names(final_assignments) <- rownames(scores_r3)

  triangulation_summary <- data.frame(
    Clusters       = paste0("k = ", k_candidates),
    WSS_Inertia    = round(wss[k_candidates], 3),
    Avg_Silhouette = round(sil_scores, 3),
    Gap_Statistic  = round(gap_vals, 3),
    Gap_Std_Error  = round(gap_se, 3),
    stringsAsFactors = FALSE
  )

  cat("--- Triangulation Summary Table (Operating in R^3 PCA Space) ---\n")
  print(triangulation_summary, row.names = FALSE)
  cat(sprintf("\nDecisions: Elbow k = %d | Silhouette k = %d | Gap (1-SE) k = %d => Consensus k = %d\n\n",
              opt_k_elbow, opt_k_sil, opt_k_gap, consensus_k))

  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    ward_tree <- stats::hclust(dist_r3, method = "ward.D2")
    par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3, 1))

    plot(elbow_k_list, wss, type = "b", pch = 19, col = "#8c564b", lwd = 2,
         xlab = "k", ylab = "WSS", main = "Elbow (Scaled Knee)")
    abline(v = opt_k_elbow, col = "firebrick", lty = 2, lwd = 1.5)
    grid()

    plot(k_candidates, sil_scores, type = "b", pch = 19, col = "#1f77b4", lwd = 2,
         xlab = "k", ylab = "Avg Silhouette", main = sprintf("Silhouette (%s)", method))
    abline(v = opt_k_sil, col = "firebrick", lty = 2, lwd = 1.5)
    grid()

    plot(k_candidates, gap_vals, type = "b", pch = 15, col = "#2ca02c", lwd = 2,
         ylim = range(c(gap_vals - gap_se, gap_vals + gap_se)),
         xlab = "k", ylab = "Gap", main = "Gap Statistic")
    arrows(k_candidates, gap_vals - gap_se, k_candidates, gap_vals + gap_se,
           length = 0.04, angle = 90, code = 3, col = "#2ca02c")
    abline(v = opt_k_gap, col = "firebrick", lty = 2, lwd = 1.5)
    grid()

    plot(stats::as.dendrogram(ward_tree), horiz = TRUE, main = "Ward Tree", xlab = "Dissimilarity")
  }

  invisible(list(
    scores        = scores_r3,
    clusters      = final_assignments,
    consensus_k   = consensus_k,
    summary_table = triangulation_summary,
    optimal_k     = list(elbow = opt_k_elbow, silhouette = opt_k_sil, gap = opt_k_gap)
  ))
}


#' Assess optimal cluster count via NbClust majority voting
#'
#' Runs the comprehensive [NbClust::NbClust()] battery across multiple clustering
#' indices, tabulates metric recommendations to determine the consensus number
#' of clusters, and extracts observation partitions.
#'
#' @param scores_matrix Numeric matrix or data frame of observation coordinates.
#' @param min_k Integer. Minimum cluster count to evaluate (default 2).
#' @param max_k Integer. Maximum cluster count to evaluate (default 10).
#' @param method Character string specifying clustering algorithm (default `"kmeans"`).
#' @param distance Character string specifying distance metric (default `"euclidean"`).
#' @param index Character string specifying clustering indices to evaluate (default `"all"`).
#' @param verbose Logical. Whether to print tabulated votes and partition allocations.
#'
#' @return Invisibly returns a named list with components:
#'   \item{optimal_k}{Consensus optimal number of clusters selected by majority vote.}
#'   \item{votes_table}{Data frame cross-tabulating candidate cluster counts and metric votes.}
#'   \item{partitions}{Data frame mapping each case to its assigned cluster number.}
#'   \item{nb_object}{Underlying output object from [NbClust::NbClust()].}
#' @export
#'
#' @examples
#' set.seed(42)
#' sim_mat <- matrix(rnorm(60), ncol = 3, dimnames = list(paste0("Case_", 1:20), NULL))
#' assess_nbclust_consensus(sim_mat, min_k = 2, max_k = 4, verbose = FALSE)
assess_nbclust_consensus <- function(scores_matrix,
                                     min_k = 2,
                                     max_k = 10,
                                     method = "kmeans",
                                     distance = "euclidean",
                                     index = "all",
                                     verbose = TRUE) {
  obs_coords <- as.matrix(scores_matrix)

  if (is.null(rownames(obs_coords))) {
    rownames(obs_coords) <- paste0("Obs_", seq_len(nrow(obs_coords)))
  }

  nb_res <- suppressWarnings(
    NbClust::NbClust(
      data     = obs_coords,
      distance = distance,
      min.nc   = min_k,
      max.nc   = max_k,
      method   = method,
      index    = index
    )
  )

  # Extract majority voting distribution from recommended cluster count per index
  votes_raw   <- as.numeric(nb_res$Best.nc[1, ])
  valid_votes <- votes_raw[!is.na(votes_raw) & is.finite(votes_raw) & votes_raw >= min_k]

  votes_table <- as.data.frame(table(Clusters_k = valid_votes))
  colnames(votes_table) <- c("Clusters_k", "Metric_Votes")

  opt_k <- as.numeric(as.character(votes_table$Clusters_k[which.max(votes_table$Metric_Votes)]))

  partition_df <- data.frame(
    Case             = rownames(obs_coords),
    Cluster          = as.numeric(nb_res$Best.partition),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  if (verbose) {
    cat(sprintf("--- NbClust Majority Voting (Method: %s, Distance: %s) ---\n", method, distance))
    print(votes_table, row.names = FALSE)
    cat(sprintf("\nNbClust Consensus Optimal Clusters: k = %d\n\n", opt_k))

    cat("--- Case Partition Assignments (Optimal k) ---\n")
    print(partition_df)
    cat("\n")
  }

  invisible(list(
    optimal_k   = opt_k,
    votes_table = votes_table,
    partitions  = partition_df,
    nb_object   = nb_res
  ))
}


# Cluster Profiling and Characterization ----

#' Characterize cluster profiles and calculate centroids
#'
#' Computes component score centroids and original indicator variable profiles across
#' cluster groupings, supporting precomputed partitions (e.g., from `NbClust`) or
#' freshly fitted k-means partitions with semantic label mapping.
#'
#' @param scores Matrix or data frame of observation component scores.
#' @param clusters Optional vector, factor, or results list containing cluster partitions.
#' @param n_clusters Integer. Target cluster count for k-means when `clusters` is omitted (default 3).
#' @param items Optional data frame or matrix of original unreduced observation indicators.
#' @param cluster_names Optional character vector or named vector mapping cluster numbers to semantic titles.
#' @param comp_names Optional character vector of names for component dimensions.
#' @param nstart Integer. Number of random initializations for k-means (default 25).
#' @param seed Optional integer for reproducible k-means clustering.
#' @param verbose Logical. If `TRUE`, prints centroid and item summary profiles.
#'
#' @return Invisibly returns a named list with components:
#'   \item{assignments}{Data frame mapping cases to numeric IDs and semantic cluster names.}
#'   \item{score_centroids}{Data frame of mean component scores per cluster.}
#'   \item{item_profiles}{Data frame of mean original indicator values per cluster (if `items` provided).}
#'   \item{kmeans_fit}{Fitted `kmeans` object, or `NULL` if external partitions were supplied.}
#' @export
#'
#' @examples
#' set.seed(42)
#' coords <- matrix(rnorm(60), ncol = 3, dimnames = list(paste0("Paper_", 1:20), paste0("PC", 1:3)))
#' characterize_pca_clusters(coords, n_clusters = 2, verbose = FALSE)
characterize_pca_clusters <- function(scores,
                                      clusters = NULL,
                                      n_clusters = 3,
                                      items = NULL,
                                      cluster_names = NULL,
                                      comp_names = NULL,
                                      nstart = 25,
                                      seed = 42,
                                      verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  scores <- as.matrix(scores)
  n_obs  <- nrow(scores)

  if (!is.null(comp_names)) {
    colnames(scores)[seq_along(comp_names)] <- comp_names
  }

  # Use exact partition if provided, otherwise compute k-means
  if (!is.null(clusters)) {
    if (is.list(clusters) && "partitions" %in% names(clusters)) {
      c_ids <- as.numeric(clusters$partitions$Cluster)
    } else {
      c_ids <- as.numeric(clusters)
    }
    km_fit <- NULL
  } else {
    km_fit <- stats::kmeans(scores, centers = n_clusters, nstart = nstart)
    c_ids  <- km_fit$cluster
  }

  # Map cluster IDs to semantic names via positional or named vector
  if (is.null(cluster_names)) {
    cluster_names  <- paste("Cluster", sort(unique(c_ids)))
    assigned_names <- cluster_names[c_ids]
  } else if (!is.null(names(cluster_names)) && all(as.character(unique(c_ids)) %in% names(cluster_names))) {
    assigned_names <- cluster_names[as.character(c_ids)]
  } else {
    assigned_names <- cluster_names[c_ids]
  }
  assigned_factor <- factor(assigned_names, levels = unique(cluster_names))

  # Calculate centroids grouped by both Cluster_Num and Cluster_Name
  score_centroids <- stats::aggregate(
    scores,
    by  = list(Cluster_Num = c_ids, Cluster_Name = assigned_factor),
    FUN = mean
  )
  score_cols <- setdiff(names(score_centroids), c("Cluster_Num", "Cluster_Name", "Cluster"))
  score_centroids[score_cols] <- round(score_centroids[score_cols], 3)

  # Item profiling grouped by both Cluster_Num and Cluster_Name
  metric_profiles <- NULL
  if (!is.null(items)) {
    items_matched <- as.data.frame(items)
    if (!is.null(rownames(scores)) && all(rownames(scores) %in% rownames(items_matched))) {
      items_matched <- items_matched[rownames(scores), , drop = FALSE]
    }
    metric_profiles <- stats::aggregate(
      items_matched,
      by  = list(Cluster_Num = c_ids, Cluster_Name = assigned_factor),
      FUN = mean
    )
    metric_cols <- setdiff(names(metric_profiles), c("Cluster_Num", "Cluster_Name", "Cluster"))
    metric_profiles[metric_cols] <- round(metric_profiles[metric_cols], 3)
  }

  assignments_df <- data.frame(
    Case         = if (!is.null(rownames(scores))) rownames(scores) else paste0("Obs_", seq_len(n_obs)),
    Cluster_Num  = c_ids,
    Cluster_Name = assigned_factor,
    Cluster      = assigned_factor,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(sprintf("--- Cluster Centroids in Component Space (k = %d) ---\n", length(unique(c_ids))))
    print(score_centroids)
    cat("\n")
    if (!is.null(metric_profiles)) {
      cat("--- Original Item Profiles by Cluster ---\n")
      print(metric_profiles)
      cat("\n")
    }
  }

  invisible(list(
    assignments     = assignments_df,
    score_centroids = score_centroids,
    item_profiles   = metric_profiles,
    kmeans_fit      = km_fit
  ))
}
