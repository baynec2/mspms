#' log2fc_t_test
#'
#' Gives log2fc and statistics from t tests across conditions. Gives differences relative to time 0 within each condition.
#'
#' @param prepared_for_stats = this is data that has been prepared using prepare_for_stats().
#'
#' @return a tibblecontaining log2fc and t test statistics across conditions
#' @export
#'
#' @examples
#'
#' log2fc_and_t_test = log2fc_t_test(mspms::prepared_for_stats)
#'
#'
log2fc_t_test = function(prepared_for_stats){

  # Calculating the log2fc
  log2fc = mspms::mspms_log2fc(prepared_for_stats)

  # Calculating the stats
  stat = mspms::mspms_t_tests(prepared_for_stats)

  #Combining them together
  log2fc_stats = dplyr::inner_join(log2fc,stat,by = c("Peptide","comparison")) %>%
    dplyr::mutate(time = .data$group2,
           condition = .data$condition.x) %>%
    dplyr::select(-.data$condition.x,-.data$condition.y) %>%
    tibble::as_tibble()

  return(log2fc_stats)

}
