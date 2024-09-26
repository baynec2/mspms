#' calc_per_samples_library_nd
#' Calculate the percentage of samples each library_id peptide was not detected
#' in.
#' @param processed_qf a QFeatures object with a SummarizedExperiment named
#' "peptides". Intended to be prepared by one of the pre-processing
#' prepare_x_data functions of the mspms R package.
#' @param peptide_library_ids a character vector containing the names of the
#' library_ids
#' @return a tibble containing percentage of samples each library id was
#' detected in,  both as full length, and as cleavage products.
#' @keywords internal

calc_per_samples_library_nd <- function(processed_qf,
                                        peptide_library_ids =
                                          mspms::peptide_library$library_id) {
  # tidying the data so we can work with it
  mspms_data <- mspms_tidy(processed_qf, "peptides")
  # calculating n of samples
  n_samples <- length(unique(mspms_data$quantCols))
  # converting to wide format, easier to reason with this way
  wide <- mspms_data %>%
    dplyr::select(
      "peptide", "library_id", "peptides",
      "peptide_type", "quantCols"
    )
  # considering full length peptides only
  full_length <- wide %>%
    dplyr::filter(.data$peptide_type == "full_length") %>%
    tidyr::pivot_wider(
      names_from = "quantCols",
      values_from = "peptides"
    )
  # What peptides are not detected at all
  nd_full <- peptide_library_ids[peptide_library_ids %!in% unique(
    full_length$library_id
  )]
  # Building a tibble of the peptides not detected
  nd_full <- tibble::tibble(
    library_id = nd_full,
    n_samples = n_samples,
    n_missing = n_samples
  )
  # Figuring the number of full length library ids missing per sample
  n_missing_full_length <- tibble::tibble(
    library_id = full_length$library_id,
    n_samples = n_samples,
    n_missing = rowSums(is.na(
      full_length
    ))
  )
  # Combining the completely missing with partially missing data
  full <- dplyr::bind_rows(n_missing_full_length, nd_full) %>%
    dplyr::mutate(peptide_type = "full_length", .after = "library_id")
  # Now considering cleavage products, conting the number of non nas, per sample
  cleavage_product <- wide %>%
    dplyr::filter(.data$peptide_type == "cleavage_product") %>%
    dplyr::select(-"peptide") %>%
    tidyr::pivot_wider(
      names_from = "quantCols",
      values_from = "peptides",
      values_fn = ~ sum(!is.na(.))
    )
  # If there are 0 non NAs, there must only be NA values for that sample
  cleavage_product[cleavage_product == 0] <- NA
  # Creating a tibble with all the data
  cp_row_sums <- tibble::tibble(
    library_id = cleavage_product$library_id,
    n_samples = n_samples,
    n_missing = rowSums(is.na(cleavage_product))
  )
  nd_cleavage <- peptide_library_ids[peptide_library_ids %!in% unique(
    cleavage_product$library_id
  )]
  nd_cleavage <- tibble::tibble(
    library_id = nd_cleavage,
    n_samples = n_samples,
    n_missing = n_samples
  )
  cp_final <- dplyr::bind_rows(
    cp_row_sums,
    nd_cleavage
  ) %>%
    dplyr::mutate(peptide_type = "cleavage_product", .after = "library_id")
  # Combining all data
  out <- dplyr::bind_rows(full, cp_final) %>%
    dplyr::mutate(per_samples_undetected = .data$n_missing /
      .data$n_samples * 100)
  return(out)
}
#' prepare_qc_check
#' Run simple quality control checks on the data. This checks to see how many
#' peptides belonging to the library were identified in the data in
#' each sample. Computes full length, and cleavage products independantly.
#' @param processed_qf a QFeatures object with a SummarizedExperiment named
#' "peptides". Intended to be prepared by one of the pre-processing
#' prepare_x_data functions of the mspms R package.
#' @param peptide_library_ids a character vector containing the names of the
#' library_ids
#' @return a tibble containing percentage of library_ids detected per sample,
#' both as full length, and as cleavage products.
#'
#' @keywords internal

prepare_qc_check_data <- function(processed_qf,
                                  peptide_library_ids =
                                    mspms::peptide_library$library_id) {
  # Filtering to only include detected peptides
  long <- processed_qf %>%
    mspms_tidy(se_name = "peptides") %>%
    dplyr::filter(!is.na(.data$peptides))
  # Counting # of peptides in peptide library
  library_num <- length(peptide_library_ids)
  # checking to see what percentage of the library is detected in each sample.
  check <- long %>%
    dplyr::group_by(
      .data$quantCols,
      .data$peptide_type,
      .data$condition,
      .data$time,
      .data$group
    ) %>%
    dplyr::summarise(
      n_detected = dplyr::n_distinct(.data$library_id),
      n_total = library_num,
      per_library_id_detected = round(.data$n_detected / .data$n_total
        * 100, 2),
      per_library_id_undetected = 100 -
        .data$per_library_id_detected
    )
  return(check)
}
#' icelogo_col_scheme
#' Defining a color scheme for our iceLogos
#'
#' @return a ggseqlogo color scheme function
#' @keywords internal
#'
icelogo_col_scheme <- function() {
  col_scheme <- ggseqlogo::make_col_scheme(
    chars = c(
      "G", "S", "T", "Y", "C", "N", "Q", "K", "R", "H", "D", "E", "P", "A", "W",
      "F", "L", "I", "M", "V", "n", "X"
    ),
    group = c(
      rep("Polar", 5), rep("Neutral", 2), rep("Basic", 3), rep("Acidic", 2),
      rep("Hydrophobic", 9),
      "Past Terminus"
    ),
    col = c(
      rep("#058644", 5), rep("#720091", 2), rep("#0046C5", 3),
      rep("#C5003E", 2),
      rep("#2E2E2E", 9), "#808080"
    ),
    name = "chemistry"
  )
  return(col_scheme)
}
#' count_cleavages_per_pos
#'
#' Count the number of cleavages per position
#'
#' @param data a tibble containing columns named peptide,cleavage_pos,condition,
#' and time. Other column names can be included.
#' @return a ggplot2 object
#' @keywords internal

count_cleavages_per_pos <- function(data,
                                    peptide_library = mspms::peptide_library) {
  peptide_length <- unique(nchar(peptide_library$library_match_sequence))

  positions <- seq_len(peptide_length - 1)
  # Counting
  count <- data %>%
    dplyr::group_by(.data$condition, .data$time, .data$cleavage_pos) %>%
    dplyr::summarise(n = dplyr::n())
  # Adding missing positions
  final <- tibble::tibble()
  for (i in unique(data$condition)) {
    f <- count %>% dplyr::filter(.data$condition == i)
    for (t in unique(f$time)) {
      count_f <- dplyr::filter(
        f,
        .data$condition == i,
        .data$time == t
      )
      missing <- positions[!positions %in% count_f$cleavage_pos]
      missing_t <- tibble::tibble(
        condition = i, time = t,
        cleavage_pos = missing, n = rep(
          0,
          length(missing)
        )
      )
      out <- dplyr::bind_rows(count_f, missing_t)
      final <- dplyr::bind_rows(final, out)
    }
  }
  final$time <- as.character(final$time)
  return(final)
}

#' prepare_for_PCA()
#'
#' prepare QFeatures object for PCA analysis
#'
#' @param mspms_tidy_data tidy mspms data (prepared from QFeatures object
#'  by mspms_tidy())
#' @param value_colname the name of the column containing values.
#'
#' @return a tibble
#' @keywords internal
prepare_for_PCA <- function(mspms_tidy_data, value_colname = "peptides_norm") {
  value_colname <- dplyr::sym(value_colname)
  # dealing with no visible binding for global variable ‘.’ NOTE
  . <- NULL
  PCA_data <- mspms_tidy_data %>%
    dplyr::select(
      "quantCols", "peptide", "group", "condition",
      "time", !!value_colname
    ) %>%
    tidyr::pivot_wider(
      names_from = "peptide", values_from = !!value_colname,
      values_fn = NULL
    ) %>%
    # if a peptide has an na remove it
    dplyr::select_if(~ !any(is.na(.)))
  dat <- PCA_data %>%
    dplyr::select(5:length(.))
  # prcomp
  prcomp <- prcomp(dat, center = TRUE, scale = TRUE)
  # Getting PCA values
  PCA_df <- stats::predict(prcomp, dat) %>%
    tibble::as_tibble()
  # Extracting the metadata.
  md <- PCA_data %>%
    dplyr::select(seq_len(4))
  Prop_of_var <- data.frame(summary(prcomp)$importance)[2, ]
  all <- dplyr::bind_cols(md, PCA_df) %>%
    dplyr::mutate(time = as.factor(.data$time))
  list <- list(all, Prop_of_var)
  return(list)
}
