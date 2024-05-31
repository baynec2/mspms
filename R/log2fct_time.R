#' log2fct_time
#'
#' preforms log2fc and t test analysis for a user specified condition comparing
#'  two specified times.
#'
#' @param mspms_data processed mspms data
#' @param within_condition the condition you would like to compare within
#' @param ref_time the reference (denominator) time
#' @param comparison_time the comparison (numerator) time
#'
#' @return a tibble
#' @export
#'
#' @examples
#'
#' log2t <- log2fct_time(mspms::mspms_data,
#'   within_condition = "DMSO",
#'   ref_time = "0",
#'   comparison_time = "60"
#' )
#'
log2fct_time <- function(mspms_data,
                         within_condition,
                         ref_time,
                         comparison_time) {
  # Filtering the data set to only include the relevant conditions and time
  point.
  f <- mspms_data %>%
    dplyr::filter(
      time %in% c(ref_time, comparison_time),
      condition == within_condition
    )

  # Performing the statistics
  stat <- f %>%
    dplyr::group_by(Peptide, cleavage_seq) %>%
    rstatix::t_test(value ~ time, ref.group = ref_time) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    tibble::as_tibble()

  # Extracting the reference means
  reference_data <- f %>%
    dplyr::filter(condition == ref_time) %>%
    dplyr::group_by(.data$Peptide) %>%
    dplyr::summarise(reference_mean = mean(.data$value, na.rm = TRUE))

  # Extracting the comparison means
  comparison_data <- f %>%
    dplyr::filter(condition == comparison_time) %>%
    dplyr::group_by(.data$Peptide) %>%
    dplyr::summarise(comparison_mean = mean(.data$value, na.rm = TRUE))

  # Calculating the log2fc
  log2fc <- dplyr::inner_join(comparison_data, reference_data,
    by = c("Peptide")
  ) %>%
    dplyr::mutate(log2fc = log2(.data$comparison_mean /
      .data$reference_mean)) %>%
    dplyr::mutate(comparison = paste0(
      comparison_time,
      "/", ref_time, " within ", within_condition
    )) %>%
    tibble::as_tibble()


  # Combining log2fc with statistics
  out <- dplyr::inner_join(log2fc, stat, by = "Peptide")

  return(out)
}
