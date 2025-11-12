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
    dplyr::summarise(control_mean = mean(.data$peptides_norm,
      na.rm = TRUE
    ))
  # Extracting just the reference data
  reference_data <- mspms_data %>%
    dplyr::group_by(
      .data$peptide, !!filter_reference_by,
      !!group_by_cols
    ) %>%
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
  mspms_data <- mspms_tidy(processed_qf, "peptides_log_impute_norm")
  # Defining formula to use in t_test
  formula <- stats::as.formula(paste("peptides_log_impute_norm ~", reference_variable))
  # Doing T tests
  stat <- mspms_data %>%
    dplyr::group_by(
      .data$cleavage_seq,
      !!group_by_cols,
      .data$peptide,
      .data$cleavage_pos,
      .data$peptide_type
    ) %>%
    rstatix::t_test(
      formula = formula,
      ref.group = reference_value
    ) %>%
    dplyr::select(-"p.adj.signif") %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    rstatix::add_significance(p.col = "p.adj") %>%
    dplyr::mutate(
      comparison = paste0(
        !!remaining_variable_syms, ".", .data$group2, "/",
        !!remaining_variable_syms, ".", .data$group1
      )
    ) %>%
    tibble::as_tibble()
  return(stat)
}

#' add_peptide_data
#'
#' adds peptide information for every peptide in the data.
#'
#' @param tibble tibble you would like to add peptide info to. Must have column
#' named peptide
#' @param qf a QFeatures object with rowData for peptides. cleavage_seq,
#' cleavage_pos, and cleavage_type.
#' @return a tibble with column named peptide.
#' @keywords internal
add_peptide_data <- function(tibble, qf) {
  # extract peptide data
  pd <- SummarizedExperiment::rowData(qf[["peptides_log_impute_norm"]]) %>%
    tibble::as_tibble() %>%
    dplyr::select("peptide", "library_id", "peptide_type", "cleavage_seq", "cleavage_pos")
  # innerjoin
  out <- dplyr::inner_join(pd, tibble, by = "peptide")
  return(out)
}


#' calc_limma_design_matrix
#'
#' Calculates a limma compatible design matrix for mspms data.
#'
#' @param colData colData with condition and time variables as factors
#' @param norm_data normalized data from QFeatures object to use
#' @returns a model matrix
#' @keywords internal
calc_limma_design_matrix <- function(colData,
                                     norm_data) {
  # Create the design matrix with a valid naming format using periods
  if (length(unique(colData$condition)) == 1) {
    design <- stats::model.matrix(~ 0 + time, data = colData)
  } else {
    design <- stats::model.matrix(~ 0 + condition:time, data = colData)
  }
  colnames(design) <- apply(
    expand.grid(levels(colData$condition), levels(colData$time)), 1,
    function(x) paste0("condition.", x[1], ".time.", x[2])
  )
  rownames(design) <- colnames(norm_data) # Match sample names
  return(design)
}

#' calc_limma_contrasts
#'
#' Calculates limma contrasts given colData. The contrasts returned are pairwise
#' relative to T0 for each timepoint assayed.
#'
#' @param colData colData from mspms experiment
#' @param design_mat design_mat as returned by calc_limma_design_matrix
#' @returns a contrast matrix
#' @keywords internal
calc_limma_contrasts <- function(colData, design_mat) {
  contrast_list <- list()
  conditions <- levels(colData$condition)

  for (cond in conditions) {
    # Get all time points available for the condition
    valid_time_points <- unique(colData$time[colData$condition == cond])

    for (t in valid_time_points[-1]) { # Skip time 0 for each condition
      # Construct contrast names to match the standardized format with periods
      base_name <- paste0("condition.", cond, ".time.", t)
      ref_name <- paste0("condition.", cond, ".time.0")

      # Contrast name and formula
      contrast_name <- paste0(base_name, " - ", ref_name)
      contrast_formula <- paste0(base_name, " - ", ref_name)
      contrast_list[[contrast_name]] <- contrast_formula
    }
  }
  # Create contrast matrix
  contrast_matrix <- limma::makeContrasts(
    contrasts = unlist(contrast_list),
    levels = design_mat
  )
  return(contrast_matrix)
}
