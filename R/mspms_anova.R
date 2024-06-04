#' mspms_anova
#'
#' Performs anova on data prepared with prepare_for_stats().
#' Tests the effect of time for each condition. FDR corrected.
#'
#' @param mspms_data = this is the data that has been run through the mspms
#' id pipeline.
#'
#' @return a data frame with anova results for each peptide.
#' @export
#'
#' @examples
#' anova_results <- mspms_anova(mspms::mspms_data)
mspms_anova <- function(mspms_data) {
  anova <- mspms_data %>%
    dplyr::group_by(.data$condition, .data$Peptide,.data$cleavage_pos) %>%
    rstatix::anova_test(value ~ time) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    tibble::as_tibble()

  return(anova)
}
