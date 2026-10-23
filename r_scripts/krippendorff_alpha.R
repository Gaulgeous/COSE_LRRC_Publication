# Inter-Rater Reliability: Krippendorff's Alpha ----

#' Compute Krippendorff alpha across sheets with bootstrap confidence intervals
#'
#' Ingests multiple Excel sheets, transposes rating matrices for designated review
#' rounds, computes ordinal Krippendorff's alpha point estimates, and derives
#' non-parametric bootstrap confidence intervals across publications.
#'
#' @param excel_file Character string. File path to the multi-sheet Excel workbook.
#' @param sheets Character vector of sheet names to iterate through.
#' @param round_cols Named list of column index vectors indicating rater columns for
#'   each round (defaults to `Round_1 = 3:5` and `Round_2 = 6:8`).
#' @param n_boot Integer. Number of bootstrap resamples across units (default 1000).
#' @param ci_probs Numeric vector of length 2 specifying lower and upper percentile bounds
#'   for confidence intervals (default `c(0.025, 0.975)`).
#' @param skip Integer. Number of header rows to skip when reading each sheet (default 2).
#' @param seed Optional integer for reproducible bootstrap sampling.
#'
#' @return A data frame with one row per sheet and columns containing the category name,
#'   point estimates (`Alpha_*`), and bootstrap bounds (`Alpha_*_CI_Low`, `Alpha_*_CI_High`).
#' @export
#'
#' @examples
#' \dontrun{
#' res <- compute_krippendorff_pipeline(
#'   excel_file = "ratings.xlsx",
#'   sheets     = c("Category_A", "Category_B"),
#'   round_cols = list(Round_1 = 3:5, Round_2 = 6:8),
#'   n_boot     = 500,
#'   seed       = 42
#' )
#' }
compute_krippendorff_pipeline <- function(
  excel_file,
  sheets,
  round_cols = list(Round_1 = 3:5, Round_2 = 6:8),
  n_boot = 1000,
  ci_probs = c(0.025, 0.975),
  skip = 2,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Fallback to perfect agreement (1.0) if all ratings across units are identical
  calc_alpha <- function(mat) {
    if (length(unique(stats::na.omit(as.vector(mat)))) <= 1) return(1.0)
    suppressWarnings(irr::kripp.alpha(mat, method = "ordinal")$value)
  }

  results <- lapply(sheets, function(sheet) {
    raw <- readxl::read_excel(excel_file, sheet = sheet, skip = 2, col_names = FALSE)

    row_out <- list(Category = sheet)

    for (rnd in names(round_cols)) {
      # Transpose to [raters x units] expected by irr::kripp.alpha
      mat <- t(as.matrix(sapply(raw[, round_cols[[rnd]]], as.numeric)))
      n_units <- ncol(mat)

      row_out[[paste0("Alpha_", rnd)]] <- calc_alpha(mat)

      # Resample publications (units) with replacement
      boot_alphas <- replicate(n_boot, {
        boot_idx <- sample(n_units, size = n_units, replace = TRUE)
        calc_alpha(mat[, boot_idx])
      })

      ci <- stats::quantile(boot_alphas, probs = c(0.025, 0.975), na.rm = TRUE)
      row_out[[paste0("Alpha_", rnd, "_CI_Low")]]  <- unname(ci[1])
      row_out[[paste0("Alpha_", rnd, "_CI_High")]] <- unname(ci[2])
    }

    as.data.frame(row_out)
  })

  final_results <- do.call(rbind, results)

}


#' Compute Krippendorff alpha for dual-round rating tables
#'
#' Evaluates inter-rater agreement for two designated rating rounds from an in-memory
#' data frame, filtering incomplete cases, transposing matrices to rater-by-unit orientation,
#' and reporting against standard benchmark thresholds.
#'
#' @param df Data frame containing rater columns across rating rounds.
#' @param r1_cols Column indices or names for Round 1 raters (default `3:5`).
#' @param r2_cols Column indices or names for Round 2 raters (default `6:8`).
#' @param method Character string specifying measurement scale (`"ordinal"`, `"interval"`,
#'   `"nominal"`, or `"ratio"`). Defaults to `"ordinal"`.
#' @param verbose Logical. Whether to print formatted summary table and interpretation guides.
#'
#' @return Invisibly returns a named list with components:
#'   \item{summary_table}{Data frame summarizing rater counts, evaluated units, skipped rows, and alpha values.}
#'   \item{round_1_fit}{Underlying result object from [irr::kripp.alpha()] for Round 1.}
#'   \item{round_2_fit}{Underlying result object from [irr::kripp.alpha()] for Round 2.}
#' @export
#'
#' @examples
#' set.seed(42)
#' sample_ratings <- data.frame(
#'   ID = 1:10,
#'   Item = paste0("Paper_", 1:10),
#'   R1_A = sample(1:3, 10, replace = TRUE),
#'   R1_B = sample(1:3, 10, replace = TRUE),
#'   R1_C = sample(1:3, 10, replace = TRUE),
#'   R2_A = sample(1:3, 10, replace = TRUE),
#'   R2_B = sample(1:3, 10, replace = TRUE),
#'   R2_C = sample(1:3, 10, replace = TRUE)
#' )
#' compute_krippendorff_alpha_simple(sample_ratings, r1_cols = 3:5, r2_cols = 6:8, verbose = FALSE)
compute_krippendorff_alpha_simple <- function(df,
                                              r1_cols = 3:5,
                                              r2_cols = 6:8,
                                              method = c("ordinal", "interval", "nominal", "ratio"),
                                              verbose = TRUE) {
  method <- match.arg(method)

  run_round_alpha <- function(sub_df, round_name) {
    mat <- data.matrix(sub_df)

    complete_mask <- stats::complete.cases(mat)
    n_total       <- nrow(mat)
    n_valid       <- sum(complete_mask)
    n_skipped     <- n_total - n_valid

    if (n_valid == 0) {
      warning(sprintf("%s has 0 rows with complete ratings across all 3 annotators.", round_name))
      return(list(
        alpha     = NA_real_,
        n_units   = 0,
        n_skipped = n_skipped,
        n_raters  = ncol(mat),
        raw_fit   = NULL
      ))
    }

    # irr::kripp.alpha requires a (raters x units) matrix
    rater_by_unit_matrix <- t(mat[complete_mask, , drop = FALSE])

    fit <- irr::kripp.alpha(rater_by_unit_matrix, method = method)

    list(
      alpha     = fit$value,
      n_units   = fit$subjects,
      n_skipped = n_skipped,
      n_raters  = fit$raters,
      raw_fit   = fit
    )
  }

  res_r1 <- run_round_alpha(df[, r1_cols, drop = FALSE], "Round 1")
  res_r2 <- run_round_alpha(df[, r2_cols, drop = FALSE], "Round 2")

  summary_df <- data.frame(
    Round               = c("Round 1", "Round 2"),
    Annotators          = c(res_r1$n_raters, res_r2$n_raters),
    Evaluated_Units     = c(res_r1$n_units, res_r2$n_units),
    Skipped_Rows        = c(res_r1$n_skipped, res_r2$n_skipped),
    Measurement_Level   = rep(method, 2),
    Krippendorff_Alpha  = round(c(res_r1$alpha, res_r2$alpha), 3),
    stringsAsFactors    = FALSE
  )

  if (verbose) {
    cat("--- Inter-Rater Reliability: Krippendorff's Alpha ---\n")
    print(summary_df, row.names = FALSE)
    cat("\nBenchmark interpretation (Krippendorff, 2004):\n")
    cat("  alpha >= 0.800 : High reliability (customary threshold for conclusive data)\n")
    cat("  0.667 <= alpha < 0.800 : Tentative conclusions permissible\n")
    cat("  alpha < 0.667  : Low reliability (inadequate agreement)\n\n")
  }

  invisible(list(
    summary_table = summary_df,
    round_1_fit   = res_r1$raw_fit,
    round_2_fit   = res_r2$raw_fit
  ))
}

