#' remaining_cd_names
#'
#' @param processed_qf 
#' @param col_name 
#'
#' @return
#' @keywords internal
remaining_cd_names = function(processed_qf,reference_variable){
  
  cd_names = SummarizedExperiment::colData(processed_qf)%>% 
  as.data.frame() %>% 
  names()
  
  ex_col_names = c(reference_variable,"group","quantCols")
  remaining_names = cd_names[cd_names %!in% ex_col_names]
 
  if(length(remaining_names) != 1){
    stop("expected there to be one remaining name in colData")
  }
  return(remaining_names)
}

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
#' @keywords internal
#'
#' @examples
#'
#' log2fc <- mspms_log2fc(mspms::prepared_for_stats)
mspms_log2fc <- function(processed_qf,
                         reference_variable = "time",
                         reference_value = 0) {
  
  # converting reference_variable to symbol
  sym = dplyr::sym(reference_variable)
  # creating enquosure to determine what column to filter by
  filter_reference_by = dplyr::enquo(sym)
  
  # Determining the remaining variables
  remaining_variables <- remaining_cd_names(processed_qf,reference_variable) 
  
  remaining_variable_syms = dplyr::sym(remaining_variables)
  group_by_cols = dplyr::enquo(remaining_variable_syms)
  
  # Converting the data in the qf to tibble
  mspms_data = mspms_tidy(processed_qf)
  
  # Extracting just the control data
  control_data <- mspms_data %>%
    dplyr::filter(!!filter_reference_by == reference_value) %>%
    dplyr::group_by(!!group_by_cols, .data$peptide) %>%
    dplyr::summarise(control_mean = mean(.data$peptides_norm, na.rm = TRUE))
  
  # Extracting just the reference data
  reference_data <- mspms_data %>%
    dplyr::group_by(.data$peptide, !!filter_reference_by,!!group_by_cols) %>%
    dplyr::summarise(reference_mean = mean(.data$peptides_norm, na.rm = TRUE))
  
  # Calculating the log2fc
  log2fc <- dplyr::inner_join(control_data, reference_data, by = c(
    "peptide",
    remaining_variables
  )) %>%
    dplyr::mutate(
      comparison = paste0(
        !!remaining_variable_syms, ".",!!filter_reference_by, "/",
        !!remaining_variable_syms, ".",reference_value),
      log2fc = log2(.data$reference_mean / .data$control_mean)
    ) %>%
    tibble::as_tibble()
  
  
  return(log2fc)
}

#' mspms_t_tests
#'
#' Performs ttests for each peptide within each group with t0 as reference.
#' P value adjustment performed using FDR correction
#'
#' @param mspms_data = this is the data that has been run through the mspms
#'  workflow.
#'
#' @return a tibble with the t test statistics for each peptide within each
#' group with t0 as reference
#' @keywords internal
#'
#' @examples
#' t_tests <- mspms_t_tests(mspms::mspms_data)
mspms_t_tests <- function(processed_qf,
                          reference_variable = "time",
                          reference_value = "0") {
  
  # converting reference_variable to symbol
  sym = dplyr::sym(reference_variable)
  # creating enquosure to determine what column to filter by
  filter_reference_by = dplyr::enquo(sym)
  
  # Determining the remaining variables
  remaining_variables = remaining_cd_names(processed_qf,reference_variable) 
  
  remaining_variable_syms = dplyr::sym(remaining_variables)
  group_by_cols = dplyr::enquo(remaining_variable_syms)
  
  mspms_data = mspms_tidy(processed_qf)
  
  #Defining formula to use in t_test
  formula = as.formula(paste("peptides_norm ~",reference_variable))
  # Doing T tests
  stat <- mspms_data %>%
    dplyr::group_by(.data$cleavage_seq,
                    !!group_by_cols,
                    .data$peptide,
                    .data$cleavage_pos) %>%
    rstatix::t_test(formula = formula, ref.group = reference_value) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    dplyr::mutate(
      comparison = paste0(
        !!remaining_variable_syms, ".",.data$group2, "/",
        !!remaining_variable_syms, ".",.data$group1)) %>%
    tibble::as_tibble()
  
  return(stat)
}

