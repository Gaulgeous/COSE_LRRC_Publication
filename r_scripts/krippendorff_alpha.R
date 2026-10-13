# Dependencies:
# - readxl
# - irr
# - here

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

  # Fast alpha calculation expecting a [3 raters x N units] matrix
  calc_alpha <- function(mat) {
    if (length(unique(na.omit(as.vector(mat)))) <= 1) return(1.0)
    suppressWarnings(kripp.alpha(mat, method = "ordinal")$value)
  }

  results <- lapply(sheets, function(sheet) {
    # Skip the two header rows and pull raw data directly
    raw <- read_excel(excel_file, sheet = sheet, skip = 2, col_names = FALSE)
    
    row_out <- list(Category = sheet)
    
    for (rnd in names(round_cols)) {
      # Extract columns, ensure numeric, and transpose to [raters x units]
      mat <- t(as.matrix(sapply(raw[, round_cols[[rnd]]], as.numeric)))
      n_units <- ncol(mat)
      
      # Point estimate
      row_out[[paste0("Alpha_", rnd)]] <- calc_alpha(mat)
      
      # Resample publications (columns) with replacement
      boot_alphas <- replicate(n_boot, {
        boot_idx <- sample(n_units, size = n_units, replace = TRUE)
        calc_alpha(mat[, boot_idx])
      })
      
      ci <- quantile(boot_alphas, probs = c(0.025, 0.975), na.rm = TRUE)
      row_out[[paste0("Alpha_", rnd, "_CI_Low")]]  <- unname(ci[1])
      row_out[[paste0("Alpha_", rnd, "_CI_High")]] <- unname(ci[2])
    }
    
    as.data.frame(row_out)
  })

  final_results <- do.call(rbind, results)

}


compute_krippendorff_alpha_simple <- function(df, 
                                              r1_cols = 3:5, 
                                              r2_cols = 6:8, 
                                              method = c("ordinal", "interval", "nominal", "ratio"),
                                              verbose = TRUE) {
  method <- match.arg(method)
  
  # Helper to process a single round
  run_round_alpha <- function(sub_df, round_name) {
    # Coerce annotator ratings to numeric
    mat <- data.matrix(sub_df)
    
    # Retain strictly complete rows (where annotator 1, 2, AND 3 are present)
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
    # mat[complete_mask, ] is (units x raters) -> transpose to (3 x valid_units)
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
  
  # 1. Compute for Round 1 (Cols 3-5) and Round 2 (Cols 6-8)
  res_r1 <- run_round_alpha(df[, r1_cols, drop = FALSE], "Round 1")
  res_r2 <- run_round_alpha(df[, r2_cols, drop = FALSE], "Round 2")
  
  # 2. Formatted Summary Table
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