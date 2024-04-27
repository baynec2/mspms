#' mspms_anova
#'
#' Performs anova on data prepared with prepare_for_stats(). Tests the effect of time for each condition. FDR corrected.
#'
#' @param prepared_for_stats = this is the data that has been prepared using prepare_for_stats()
#'
#' @return a data frame with anova results for each peptide.
#' @export
#'
#' @examples
#' anova_results = mspms_anova(mspms::prepare_for_stats)
mspms_anova = function(prepared_for_stats){

  anova = prepared_for_stats %>%
    dplyr::group_by(condition,Peptide) %>%
    rstatix::anova_test(value ~ time) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    tibble::as_tibble()

  return(anova)

}
