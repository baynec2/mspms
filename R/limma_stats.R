#' limma_stats
#'
#' Calculates statistics for each condition relative to time 0 using limma for
#' differential analysis. Results are then formatted to be consistent with
#' results produced by other statistic approaches used in the mspms package
#' (log2fc_t_test).
#'
#' @param processed_qf mspms data in a QFeatures object.
#'
#' @returns a tibble containing statistics
#' @export
#'
#' @examples
#' mspms_limma_results <- limma_stats(mspms::processed_qf)
limma_stats <- function(processed_qf) {
  # Preparing the data for limma analysis
  norm_data <- SummarizedExperiment::assay(
    processed_qf,
    "peptides_log_impute_norm"
  )
  colData <- SummarizedExperiment::colData(processed_qf)
  # Ensure time and condition columns are factors
  colData$time <- factor(colData$time, levels = unique(colData$time))
  colData$condition <- factor(colData$condition)
  # Create the design matrix with a valid naming format using periods
  design_mat <- calc_limma_design_matrix(colData, norm_data)
  # Fit the linear model
  fit <- limma::lmFit(norm_data, design_mat)
  # Define contrasts for each condition's time points relative to T0
  contrast_matrix <- calc_limma_contrasts(colData, design_mat)
  # Fit contrasts and apply empirical Bayes
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit3 <- limma::eBayes(fit2)
  # Retrieve results for each contrast
  results <- lapply(colnames(contrast_matrix), function(contrast_name) {
    coef_index <- which(colnames(fit3$coefficients) == contrast_name)
    if (length(coef_index) > 0) {
      result <- limma::topTable(fit3, coef = coef_index, number = Inf)
      result$contrast <- contrast_name
      return(result)
    } else {
      warning(paste(
        "Contrast", contrast_name,
        "not found in fit coefficients. Skipping."
      ))
      return(NULL)
    }
  })

  # Combine results into a single tibble
  results_combined <- do.call(rbind, results) %>%
    tibble::rownames_to_column("peptide") %>%
    # duplicate peptides got a number in process, removing here
    dplyr::mutate(peptide = gsub("\\d*", "", peptide)) %>%
    dplyr::mutate(
      # Split the contrast at " - " to separate the two parts
      group1 = sapply(
        strsplit(as.character(contrast), " - "),
        function(x) x[1]
      ),
      # Extract the condition and time for the first part (before " - ")
      condition = sapply(strsplit(group1, "\\."), function(x) x[2]),
      time = sapply(strsplit(group1, "\\."), function(x) as.numeric(x[4])),
      # Now do the same for the second part (after " - ")
      group2 = sapply(strsplit(as.character(contrast), " - "), function(x) x[2])
    ) %>%
    rstatix::add_significance("adj.P.Val", "p.adj.signif")
  # Add cleavage data back to table
  results_combined <- add_peptide_data(results_combined, processed_qf)
  # Adding control mean to results
  t0 <- mspms::mspms_tidy(processed_qf, "peptides_log_impute_norm") %>%
    dplyr::filter(time == 0) %>%
    dplyr::group_by(.data$peptide, .data$condition) %>%
    dplyr::summarise(control_mean = mean(.data$peptides_log_impute_norm))

  # Convert results to a consistent format.
  results_formatted <- results_combined %>%
    dplyr::inner_join(t0, by = c("peptide", "condition")) %>%
    dplyr::select("condition", "peptide", "library_id", "peptide_type", "control_mean",
      "time",
      "sample_mean" = "AveExpr", "comparison" = "contrast",
      "log2fc" = "logFC", "cleavage_seq", "cleavage_pos", "group1",
      "group2", "p" = "P.Value", "p.adj" = "adj.P.Val", "p.adj.signif"
    )
  return(results_formatted)
}
