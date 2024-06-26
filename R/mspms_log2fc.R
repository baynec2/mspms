#' mspms_log2fc
#'
#' Calculates the log2fc for each timepoint within each condition relative to
#' time 0.
#'
#' @param mspms_data = this is the data that has been run through the mspms
#'  pipeline
#'
#' @return a data frame with the log2fc for each timepoint within each condition
#' relative to time 0
#' @export
#'
#' @examples
#'
#' log2fc <- mspms_log2fc(mspms::prepared_for_stats)
mspms_log2fc <- function(mspms_data) {
  # Extracting just the control data
  control_data <- mspms_data %>%
    dplyr::filter(.data$time == 0) %>%
    dplyr::group_by(.data$condition, .data$Peptide) %>%
    dplyr::summarise(control_mean = mean(.data$value, na.rm = TRUE))


  # Extracting just the reference data
  reference_data <- mspms_data %>%
    dplyr::group_by(.data$Peptide, .data$time, .data$condition) %>%
    dplyr::summarise(reference_mean = mean(.data$value, na.rm = TRUE))

  # Calculating the log2fc

  log2fc <- dplyr::inner_join(control_data, reference_data, by = c(
    "Peptide",
    "condition"
  )) %>%
    dplyr::mutate(
      comparison = paste0(
        .data$condition, ".T", .data$time, "/",
        .data$condition, ".T0"
      ),
      log2fc = log2(.data$reference_mean / .data$control_mean)
    ) %>%
    tibble::as_tibble()


  return(log2fc)
}
