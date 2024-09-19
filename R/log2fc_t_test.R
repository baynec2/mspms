#' log2fc_t_test
#'
#' Calculates the log2 fold change and t-test statistics given a user specified
#' reference variable and value. 
#'
#' @param processed_qf mspms data in a QFeatures object.
#' @param reference_variable the colData variable to use as reference
#' @param reference_value the value of the colData variable to use as reference
#'
#' @return a tibble containing log2fc and t test statistics 
#' @export
#'
#' @examples
#' log2fc_and_t_test <- log2fc_t_test(mspms::processed_qf)
log2fc_t_test <- function(processed_qf,
                          reference_variable = "time",
                          reference_value = 0) {
  
  # Calculating the log2fc
  log2fc <- mspms_log2fc(
    processed_qf,
    reference_variable,
    reference_value
  )

  # Calculating the stats
  stat <- mspms_t_tests(processed_qf,
    reference_variable,
    reference_value = as.character(reference_value)
  )

  # Combining them together
  log2fc_stats <- dplyr::inner_join(log2fc, stat, by = c(
    "peptide",
    "comparison"
  )) %>%
    tibble::as_tibble() %>% 
    dplyr::rename("condition" = "condition.x") %>% 
    dplyr::select(-"condition.y")

  return(log2fc_stats)
}
