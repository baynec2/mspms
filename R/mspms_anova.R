#' mspms_anova
#'
#' Performs anova on data prepared with prepare_for_stats(). Tests the effect of time for each condition. FDR corrected.
#'
#' @param mspms_data = this is the data that has been run through the mspms pipeline.
#'
#' @return a data frame with anova results for each peptide.
#' @export
#'
#' @examples
#' anova_results = mspms_anova(mspms::mspms_data)
mspms_anova = function(mspms_data){

  anova = prepared_for_stats %>%
    dplyr::group_by(.data$condition,.data$Peptide) %>%
    rstatix::anova_test(value ~ time) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    tibble::as_tibble()

  return(anova)

}
