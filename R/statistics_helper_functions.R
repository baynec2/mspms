#' remaining_cd_names
#'
#' determine what the remaining colData names are when removing the reference
#' variable.
#'
#' @param processed_qf a QFeatures object
#' @param reference_variable name of reference variable
#'
#' @return a vector of the remaining names in the colData
#' @keywords internal
remaining_cd_names <- function(processed_qf, reference_variable) {
  cd_names <- SummarizedExperiment::colData(processed_qf) %>%
    as.data.frame() %>%
    names()

  ex_col_names <- c(reference_variable, "group", "quantCols")
  remaining_names <- cd_names[cd_names %!in% ex_col_names]

  if (length(remaining_names) != 1) {
    stop("expected there to be one remaining name in colData")
  }
  return(remaining_names)
}

#' mspms_log2fc
#'
#' calculates the log2fc for each time point within each condition relative to
#' a specified value for a specified reference variable.
#' @param processed_qf a QFeatures object with a SummarizedExperiment named
#' "peptides_norm".
#' @param reference_variable the variable to used as a reference (denominator
#' of log2 fold change).
#' @param reference_value  the value of the reference variable to use as
#' the reference
#' @return a tibble with the t test statistics for each peptide within each
#' group with the supplied value at the supplied variable as reference.
#' @keywords internal
mspms_log2fc <- function(processed_qf,
                         reference_variable = "time",
                         reference_value = 0) {
  # converting reference_variable to symbol
  sym <- dplyr::sym(reference_variable)
  # creating enquosure to determine what column to filter by
  filter_reference_by <- dplyr::enquo(sym)
  # Determining the remaining variables
  remaining_variables <- remaining_cd_names(processed_qf, reference_variable)
  remaining_variable_syms <- dplyr::sym(remaining_variables)
  group_by_cols <- dplyr::enquo(remaining_variable_syms)
  # Converting the data in the qf to tibble
  mspms_data <- mspms_tidy(processed_qf)
  # Extracting just the control data
  control_data <- mspms_data %>%
    dplyr::filter(!!filter_reference_by == reference_value) %>%
    dplyr::group_by(!!group_by_cols, .data$peptide) %>%
    dplyr::summarise(control_mean = mean(.data$peptides_norm, na.rm = TRUE))
  # Extracting just the reference data
  reference_data <- mspms_data %>%
    dplyr::group_by(.data$peptide, !!filter_reference_by, !!group_by_cols) %>%
    dplyr::summarise(
      sample_mean = mean(.data$peptides_norm, na.rm = TRUE),
      sample_sd = stats::sd(.data$peptides_norm, na.rm = TRUE)
    )
  # Calculating the log2fc
  log2fc <- dplyr::inner_join(control_data, reference_data, by = c(
    "peptide",
    remaining_variables
  )) %>%
    dplyr::mutate(
      comparison = paste0(
        !!remaining_variable_syms, ".", !!filter_reference_by, "/",
        !!remaining_variable_syms, ".", reference_value
      ),
      log2fc = log2(.data$sample_mean / .data$control_mean),
      difference = .data$sample_mean - .data$control_mean
    ) %>%
    tibble::as_tibble()
  return(log2fc)
}

#' mspms_t_tests
#'
#' performs t-tests for each peptide within each group for the user specified.
#' FDR adjustment is performed.
#'
#' @param processed_qf a QFeatures object with a SummarizedExperiment named
#' "peptides_norm".
#' @param reference_variable the variable to used as a reference.
#' @param reference_value  the value of the reference variable to use as
#' the reference
#' @return a tibble with the t test statistics for each peptide within each
#' group with the supplied value at the supplied variable as reference.
#' @keywords internal
mspms_t_tests <- function(processed_qf,
                          reference_variable = "time",
                          reference_value = "0") {
  if (!is.character(reference_value)) {
    reference_value <- as.character(reference_value)
  }
  # converting reference_variable to symbol
  sym <- dplyr::sym(reference_variable)
  # creating enquosure to determine what column to filter by
  filter_reference_by <- dplyr::enquo(sym)
  # Determining the remaining variables
  remaining_variables <- remaining_cd_names(processed_qf, reference_variable)
  remaining_variable_syms <- dplyr::sym(remaining_variables)
  group_by_cols <- dplyr::enquo(remaining_variable_syms)
  mspms_data <- mspms_tidy(processed_qf)
  # Defining formula to use in t_test
  formula <- stats::as.formula(paste("peptides_norm ~", reference_variable))
  # Doing T tests
  stat <- mspms_data %>%
    dplyr::group_by(
      .data$cleavage_seq,
      !!group_by_cols,
      .data$peptide,
      .data$cleavage_pos
    ) %>%
    rstatix::t_test(
      formula = formula,
      ref.group = reference_value
    ) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    dplyr::mutate(
      comparison = paste0(
        !!remaining_variable_syms, ".", .data$group2, "/",
        !!remaining_variable_syms, ".", .data$group1
      )
    ) %>%
    tibble::as_tibble()
  return(stat)
}
