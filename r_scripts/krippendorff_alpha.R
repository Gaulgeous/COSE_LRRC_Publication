# Dependencies:
# - readxl
# - irr
# - here

compute_krippendorff_pipeline <- function(
  excel_file,
  sheets,
  round_cols = list(Round_1 = 2:4, Round_2 = 5:7),
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
