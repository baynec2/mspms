#' log2fct_condition
#'
#' preforms log2fc and t test analysis at a user specified time point for
#' two specified conditions.
#'
#' @param mspms_data processed mspms data
#' @param at_time the time you would like to use for comparison
#' @param ref_condition the reference (denominator) condition
#' @param comparison_condition the comparison (numerator) condition
#'
#' @return a tibble
#' @export
#'
#' @examples
#' log2t <- log2fct_condition(mspms::mspms_data,
#'   at_time = 60,
#'   ref_condition = "DMSO",
#'   comparison_condition = "MZB"
#' )
log2fct_condition <- function(mspms_data,
                              at_time,
                              ref_condition,
                              comparison_condition) {
  # Filtering the data set to include the relevant conditions and time point.
  f <- mspms_data %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$condition %in% c(ref_condition, comparison_condition)) %>%
    dplyr::filter(.data$time == at_time)

  # Performing the statistics
  stat <- f %>%
    dplyr::group_by(.data$Peptide, .data$cleavage_seq, .data$cleavage_pos) %>%
    rstatix::t_test(value ~ condition, ref.group = ref_condition) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    tibble::as_tibble()

  # Extracting the control data
  reference_data <- f %>%
    dplyr::filter(.data$condition == ref_condition) %>%
    dplyr::group_by(.data$Peptide) %>%
    dplyr::summarise(reference_mean = mean(.data$value, na.rm = TRUE))

  # Extracting the experimental data
  comparison_data <- f %>%
    dplyr::filter(.data$condition == comparison_condition) %>%
    dplyr::group_by(.data$Peptide) %>%
    dplyr::summarise(comparison_mean = mean(.data$value, na.rm = TRUE))

  # Calculating the log2fc
  log2fc <- dplyr::inner_join(comparison_data, reference_data,
    by = c("Peptide")
  ) %>%
    dplyr::mutate(log2fc = log2(.data$comparison_mean /
      .data$reference_mean)) %>%
    dplyr::mutate(comparison = paste0(
      comparison_condition,
      "/", ref_condition, " at ", at_time
    )) %>%
    tibble::as_tibble()


  # Combining log2fc with statistics
  out <- dplyr::inner_join(log2fc, stat, by = "Peptide")

  return(out)
}
