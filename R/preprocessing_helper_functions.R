# Adding Global Variables for use with NSE
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "peptide",
    "contrast",
    "group1",
    "group2",
    "group3",
    "time",
    "peptide_length",
    "quantCols",
    "Protein.Group",
    "Stripped.Sequnce",
    "library_match_sequence",
    "cleavage_pos",
    "prev_aa_i",
    "next_aa_i",
    "library_id",
    "library_real_sequence",
    "Prev AA",
    "Peptide Sequence",
    "Next AA",
    "Protein",
    "q_value",
    "Stripped.Sequence",
    "proteins",
    "terminus",
    "terminus_cleavage_pos",
    ":="
  ))
}

# Dummy reference to avoid R CMD check note
.onLoad <- function(libname, pkgname) {
  if (FALSE) {
    imputeLCMD::impute.MAR()
  }
}


# Negate
`%!in%` <- Negate(`%in%`)

#' Generic preparation function for MSP-MS input files
#'
#' @param filepath path to the input file
#' @param colData_filepath path to colData CSV
#' @param peptide_library peptide library used in experiment
#' @param n_residues number of residues to include around cleavage site
#' @param read_fun function to read the file (e.g., read_tsv, read_csv, read_delim)
#' @param validate_fun function to validate the file (check_file_is_valid_*)
#' @param transform_fun function to transform file into standard peptide format
#'
#' @return a QFeatures object
#' @keywords internal
prepare_file <- function(filepath,
                         colData_filepath,
                         peptide_library = mspms::peptide_library,
                         n_residues = 4,
                         read_fun,
                         validate_fun,
                         transform_fun) {
  # Step 1: Read file
  raw_data <- read_fun(filepath)

  # Step 2: Validate
  validate_fun(raw_data)

  # Step 3: Transform into standard peptide format
  peptides <- transform_fun(raw_data, peptide_library)

  # Step 4: Load colData
  colData <- load_colData(colData_filepath)

  # Step 5: Convert to QFeatures
  qf <- prepared_to_qf(peptides, colData, peptide_library, n_residues)

  return(qf)
}


#' Transform PEAKS LFQ file into standard peptide format
#'
#' @param df PEAKS LFQ data read in with read_csv
#' @param peptide_library peptide library used in the experiment
#' @param quality_threshold minimum peptide quality to keep (default 0.3)
#'
#' @return a tibble with columns: peptide, library_id, sample intensities
#' @keywords internal
transform_peaks <- function(df, peptide_library, quality_threshold = 0.3) {
  # Filter out low quality peptides
  df <- df %>% dplyr::filter(.data$Quality > quality_threshold)

  # Report how many were removed
  n_removed <- sum(df$Quality <= quality_threshold, na.rm = TRUE)
  if (n_removed > 0) {
    percentage <- round(n_removed / nrow(df) * 100)
    message(
      n_removed, " peptides removed (quality < ", quality_threshold,
      ") = ", percentage, "% of total"
    )
  }

  # Identify relevant columns: from Avg. Area to Sample Profile (Ratio)
  start <- which(names(df) == "Avg. Area") + 1
  end <- which(names(df) == "Sample Profile (Ratio)") - 1

  # Create the standard peptide column
  df <- df %>%
    dplyr::mutate(Peptide = gsub("\\.", "_", .data$Peptide)) %>%
    dplyr::select(
      peptide = "Peptide",
      library_id = "Protein Accession",
      dplyr::any_of(start:end)
    )

  # Join with peptide library to ensure consistent metadata
  df <- dplyr::inner_join(peptide_library, df, by = "library_id")

  # Remove unnecessary library columns before returning
  df <- df %>%
    dplyr::select(peptide, library_id, dplyr::everything()) %>%
    dplyr::select(-library_match_sequence, -library_real_sequence)

  return(df)
}

#' Transform Proteome Discoverer PeptideGroups.txt into standard peptide format
#'
#' @param df Proteome Discoverer PeptideGroups.txt read with read_delim
#' @param peptide_library peptide library used in the experiment
#'
#' @return a tibble with columns: peptide, library_id, sample intensities
#' @keywords internal
transform_pd <- function(df, peptide_library) {
  # Rename to standard columns
  df <- df %>%
    dplyr::rename(
      peptide = "Annotated Sequence",
      library_id = "Master Protein Accessions"
    )

  # Convert peptide notation from PD to standard format
  df <- df %>%
    dplyr::mutate(
      peptide = gsub("\\[-\\]", "", .data$peptide),
      peptide = gsub("^\\.", "", .data$peptide),
      peptide = gsub("\\.$", "", .data$peptide),
      peptide = gsub("\\].", "_", .data$peptide),
      peptide = gsub("\\.\\[", "_", .data$peptide),
      peptide = gsub("\\]", "", .data$peptide),
      peptide = gsub("\\[", "", .data$peptide)
    )

  # Remove abundance normalization/scaling/count columns
  df <- df %>%
    dplyr::select(
      -dplyr::contains("Normalized"),
      -dplyr::contains("Scaled"),
      -dplyr::contains("Count")
    )

  # Keep only peptide, library_id, and abundance columns
  df <- df %>%
    dplyr::select("peptide", "library_id", dplyr::contains("Abundance")) %>%
    dplyr::rename_with(~ gsub("Abundance: ", "", .)) %>%
    dplyr::rename_with(~ gsub(": Sample", "", .))

  # Join with peptide library to ensure consistent library_id mapping
  df <- dplyr::inner_join(peptide_library, df, by = "library_id")

  # Remove extra library columns
  df <- df %>%
    dplyr::select(peptide, library_id, dplyr::everything()) %>%
    dplyr::select(-library_match_sequence, -library_real_sequence)

  return(df)
}


#' Transform DIA-NN report (pr_matrix.tsv) into standard peptide format
#'
#' @param df DIA-NN pr_matrix.tsv read with read_tsv
#' @param peptide_library peptide library used in the experiment
#'
#' @return a tibble with columns: peptide, library_id, sample intensities
#' @keywords internal
transform_diann <- function(df, peptide_library) {
  # start and and end indicies of quantification columns
  start <- which(names(df) == "Precursor.Id") + 1
  end <- length(names(df)) - 2
  nsamples <- end - start

  # Convert to peptides by summarizing precursor intensities (sum) per peptide
  peptides <- df %>%
    dplyr::group_by(
      library_id = `Protein.Group`,
      peptide = `Stripped.Sequence`
    ) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(names(df)[start:end]), \(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    )

  peptides <- dplyr::inner_join(peptide_library, peptides, by = "library_id") %>%
    dplyr::mutate(
      cleavage_pos = stringr::str_locate(library_match_sequence, peptide),
      prev_aa_i = cleavage_pos[, "start"] - 1,
      next_aa_i = cleavage_pos[, "end"] + 1,
      prev_aa = dplyr::if_else(prev_aa_i > 0, stringr::str_sub(library_match_sequence, prev_aa_i, prev_aa_i), ""),
      next_aa = dplyr::if_else(next_aa_i <= stringr::str_length(library_match_sequence), stringr::str_sub(library_match_sequence, next_aa_i, next_aa_i), "")
    ) %>%
    dplyr::mutate(peptide = dplyr::case_when(
      prev_aa != "" ~ paste0(prev_aa, "_", peptide),
      next_aa != "" ~ paste0(peptide, "_", next_aa),
      TRUE ~ peptide
    ))
  peptides <- peptides %>%
    dplyr::select(1:(length(names(peptides)) - 5)) %>%
    dplyr::select(peptide, library_id, dplyr::everything()) %>%
    dplyr::select(-library_match_sequence, -library_real_sequence)

  return(peptides)
}


#' Transform FragPipe combined_peptide.tsv into standard peptide format
#'
#' @param df FragPipe combined_peptide.tsv read with read_tsv
#' @param peptide_library peptide library used in the experiment
#'
#' @return a tibble with columns: peptide, library_id, sample intensities
#' @keywords internal
transform_fragpipe <- function(df, peptide_library) {
  # Replace NA Prev/Next AA with empty string
  df <- df %>%
    dplyr::mutate(
      `Prev AA` = dplyr::if_else(.data$`Prev AA` == "-", "", .data$`Prev AA`),
      `Next AA` = dplyr::if_else(.data$`Next AA` == "-", "", .data$`Next AA`)
    ) %>%
    # Build peptide with _ cleavage indicators
    dplyr::mutate(
      peptide = paste0(`Prev AA`, "_", `Peptide Sequence`, "_", `Next AA`),
      peptide = gsub("^_", "", peptide),
      peptide = gsub("_$", "", peptide)
    )

  # Select peptide, library_id, and abundance columns
  df <- df %>%
    dplyr::select(
      peptide,
      library_id = Protein,
      dplyr::contains("MaxLFQ Intensity")
    ) %>%
    dplyr::rename_with(~ gsub(" MaxLFQ Intensity", "", .))

  # Convert 0 to NA
  df[df == 0] <- NA

  # Join with peptide library to ensure consistency
  df <- dplyr::inner_join(peptide_library, df, by = "library_id") %>%
    dplyr::select(peptide, library_id, dplyr::everything()) %>%
    dplyr::select(-library_match_sequence, -library_real_sequence)

  return(df)
}


#' Transform Sage lfq.tsv into standard peptide format
#'
#' @param df sage lfq.tsv read with read_tsv
#' @param peptide_library peptide library used in the experiment
#'
#' @return a tibble with columns: peptide, library_id, sample intensities
#' @keywords internal
transform_sage <- function(df, peptide_library) {
  # start and and end indicies of quantification columns
  start <- which(names(df) == "spectral_angle") + 1
  end <- length(names(df))

  # Convert to peptides by summarizing precursor intensities (sum) per peptide
  peptides <- df %>%
    # Only including peptides confidently detected FDR < 0.01
    dplyr::filter(q_value <= 0.01) %>%
    dplyr::group_by(
      library_id = `proteins`,
      peptide
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(names(df)[start:end]),
        \(x) sum(x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  peptides <- dplyr::inner_join(peptide_library, peptides, by = "library_id") %>%
    dplyr::mutate(
      cleavage_pos = stringr::str_locate(library_match_sequence, peptide),
      prev_aa_i = cleavage_pos[, "start"] - 1,
      next_aa_i = cleavage_pos[, "end"] + 1,
      prev_aa = dplyr::if_else(prev_aa_i > 0, stringr::str_sub(library_match_sequence, prev_aa_i, prev_aa_i), ""),
      next_aa = dplyr::if_else(next_aa_i <= stringr::str_length(library_match_sequence), stringr::str_sub(library_match_sequence, next_aa_i, next_aa_i), "")
    ) %>%
    dplyr::mutate(peptide = dplyr::case_when(
      prev_aa != "" ~ paste0(prev_aa, "_", peptide),
      next_aa != "" ~ paste0(peptide, "_", next_aa),
      TRUE ~ peptide
    ))
  peptides <- peptides %>%
    dplyr::select(1:(length(names(peptides)) - 5)) %>%
    dplyr::select(peptide, library_id, dplyr::everything()) %>%
    dplyr::select(-library_match_sequence, -library_real_sequence)

  # Converting 0 to NA
  peptides[peptides == 0] <- NA

  # Removing .mzmL extension from filenames
  new_names <- gsub("\\.mzML", "", names(peptides))
  names(peptides) <- new_names
  return(peptides)
}


#' check_file_is_valid
#'
#' Validate that an input data frame contains all expected columns.
#'
#' @param data A data frame read into R.
#' @param expected_names A character vector of expected column names.
#' @param tool_name A short string identifying the originating software
#' (e.g., "PEAKS", "PD",etc).
#' @return Raises an error (`stop`) with an informative message if required
#' columns are missing; otherwise returns `invisible(NULL)`.
#' @keywords internal

check_file_is_valid <- function(data, expected_names, tool_name) {
  missing_names <- expected_names[expected_names %!in% names(data)]
  if (length(missing_names) > 0) {
    stop(
      paste0(
        "This doesn't look like the expected ", tool_name, " file. ",
        "The following columns are missing: ",
        paste(missing_names, collapse = ", ")
      )
    )
  }
  invisible(NULL)
}

#' check_file_is_valid_peaks
#' Check to make sure the input data looks like the expected PEAKS file.
#' @param peaks_data protein-peptides-lfq.csv file generated by PEAKS read
#' into R.
#' @return a stop command with a informative message if file looks unexpected.
#' otherwise, nothing.
#' @keywords internal
check_file_is_valid_peaks <- function(peaks_data) {
  # These are the names of the columns expected to be in every PEAKS file
  expected_names <- c(
    "Protein Group", "Protein ID", "Protein Accession", "Peptide",
    "Used", "Candidate", "Quality", "Significance", "Avg. ppm", "Avg. Area",
    "Sample Profile (Ratio)", "Group Profile (Ratio)", "Max Ratio", "#Vector",
    "Start", "End", "PTM"
  )
  check_file_is_valid(peaks_data, expected_names, "PEAKS")
}


#' check_file_is_valid_fragpipe
#' Check to make sure the input data looks like the expected FragPipe file.
#' @param fragpipe_data combined_peptide.tsv file generated by FragPipe read
#' into R.
#'
#' @return a stop command with a informative message if file looks unexpected.
#' otherwise, nothing.
check_file_is_valid_fragpipe <- function(fragpipe_data) {
  # These are the names of the columns expected to be in every fragpipe file
  expected_names <- c(
    "Peptide Sequence", "Prev AA", "Next AA", "Start", "End", "Peptide Length",
    "Charges", "Protein", "Protein ID", "Entry Name", "Gene",
    "Protein Description", "Mapped Genes", "Mapped Proteins"
  )

  check_file_is_valid(fragpipe_data, expected_names, "Fragpipe")
}

#' check_file_is_valid_pd
#' Check to make sure the input data looks like the expected ProteomeDiscoverer
#' file.
#' @param pd_data PeptideGroups.txt file generated by ProteomeDiscover and read
#' into R.
#'
#' @return a stop command with a informative message if file looks unexpected.
#' otherwise, nothing.
check_file_is_valid_pd <- function(pd_data) {
  # These are the names of the columns expected to be in every PEAKS file
  expected_names <- c(
    "Peptide Groups Peptide Group ID", "Checked", "Confidence",
    "Annotated Sequence", "Modifications", "Qvality PEP",
    "Qvality q-value", "# Protein Groups", "# Proteins",
    "# PSMs", "Master Protein Accessions",
    "Positions in Master Proteins",
    "Modifications in Master Proteins",
    "# Missed Cleavages"
  )
  check_file_is_valid(pd_data, expected_names, "Proteome Discoverer")
}
#' check_file_is_valid_diann
#' Check to make sure the input data looks like the expected DIA-NN output
#' file.
#' @param diann_data pg_matrix.tsv file generated by DIA-NN and read
#' into R.
#'
#' @return a stop command with a informative message if file looks unexpected.
#' otherwise, nothing.
check_file_is_valid_diann <- function(diann_data) {
  expected_names <- c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic",
    "Stripped.Sequence", "Modified.Sequence",
    "Precursor.Charge", "Precursor.Id", "All Mapped Proteins",
    "All Mapped Genes"
  )

  check_file_is_valid(diann_data, expected_names, "DIA-NN")
}


#' check_file_is_valid_sage
#' Check to make sure the input data looks like the expected PEAKS file.
#' @param sage_data read in lfq.tsv file output prodcued by Sage
#' into R.
#' @return a stop command with a informative message if file looks unexpected.
#' otherwise, nothing.
#' @keywords internal
check_file_is_valid_sage <- function(sage_data) {
  # These are the names of the columns expected to be in every PEAKS file
  expected_names <- c(
    "peptide", "charge", "proteins", "q_value", "score", "spectral_angle"
  )
  check_file_is_valid(sage_data, expected_names, "sage")
}

#' Generalized cleavage function with dynamic column names
#'
#' Finds cleavage sequences for either N-terminal or C-terminal cleavages
#' relative to a peptide library sequence.
#'
#' @param peptide_sequence Peptide sequence, single-letter code. "_" denotes cleavage site.
#' @param library_match_sequence Sequence matched by proteomics software (may differ from real sequence).
#' @param library_real_sequence True peptide sequence.
#' @param n_residues Number of residues to include on each side of the cleavage site.
#' @param terminus "nterm" or "cterm", specifying which terminus to analyze.
#'
#' @return tibble with peptide, cleavage sequence, and cleavage position. Column names are
#' dynamically named based on the terminus.
#' @keywords internal
cleavage <- function(peptide_sequence,
                     library_match_sequence,
                     library_real_sequence,
                     n_residues = 4,
                     terminus = c("nterm", "cterm")) {
  terminus <- match.arg(terminus)
  has_cleavage <- grepl("_", peptide_sequence)

  # Return NA if no cleavage detected
  if (!has_cleavage) {
    return(tibble::tibble(
      peptide = peptide_sequence,
      !!terminus := NA_character_,
      !!paste0(terminus, "_cleavage_pos") := NA_integer_
    ))
  }

  # Pad library sequences with Xs
  n_x <- paste0(rep("X", n_residues), collapse = "")
  x_mod_match <- paste0(n_x, library_match_sequence, n_x)
  x_mod_real <- paste0(n_x, library_real_sequence, n_x)

  if (terminus == "cterm" && grepl("_.$", peptide_sequence)) {
    temp <- gsub("_", "", substr(peptide_sequence, 1, nchar(peptide_sequence) - 1))
    left_ref_start <- regexpr(temp, x_mod_match)[[1]]
    cleavage_pos <- left_ref_start + nchar(temp) - 1
  } else if (terminus == "nterm" && regexpr("_", peptide_sequence)[[1]] == 2) {
    temp <- gsub("_", "", substr(peptide_sequence, 3, nchar(peptide_sequence)))
    cleavage_pos <- regexpr(temp, x_mod_match)[[1]] - 1
  } else {
    return(tibble::tibble(
      peptide = peptide_sequence,
      !!terminus := NA_character_,
      !!paste0(terminus, "_cleavage_pos") := NA_integer_
    ))
  }

  # Extract left and right sequences around cleavage
  left_seq <- substr(x_mod_real, cleavage_pos - n_residues + 1, cleavage_pos)
  right_seq <- substr(x_mod_real, cleavage_pos + 1, cleavage_pos + n_residues)

  tibble::tibble(
    peptide = peptide_sequence,
    !!terminus := paste0(left_seq, right_seq),
    !!paste0(terminus, "_cleavage_pos") := regexpr(temp, library_match_sequence)[[1]] +
      ifelse(terminus == "cterm", nchar(temp) - 1, -1)
  )
}

#' N-terminal cleavage sequence extraction
#'
#' Wrapper for `cleavage()` that extracts the N-terminal cleavage sequence
#' and cleavage position from a peptide relative to its library sequence.
#'
#' @inheritParams cleavage
#' @return A tibble with columns:
#'   - `peptide`: the input peptide sequence
#'   - `nterm`: the N-terminal cleavage sequence (n residues on each side)
#'   - `nterm_cleavage_pos`: the position of the N-terminal cleavage in the library sequence
#' @keywords internal
nterm_cleavage <- function(...) cleavage(..., terminus = "nterm")

#' C-terminal cleavage sequence extraction
#'
#' Wrapper for `cleavage()` that extracts the C-terminal cleavage sequence
#' and cleavage position from a peptide relative to its library sequence.
#'
#' @inheritParams cleavage
#' @return A tibble with columns:
#'   - `peptide`: the input peptide sequence
#'   - `cterm`: the C-terminal cleavage sequence (n residues on each side)
#'   - `cterm_cleavage_pos`: the position of the C-terminal cleavage in the library sequence
#' @keywords internal
cterm_cleavage <- function(...) cleavage(..., terminus = "cterm")
#' add_cleavages
#'
#' Adds cleavage information to a tibble by wraping the n_term_cleavage
#' and c_term_cleavage functions into a consolidated function.
#'
#' @param joined_with_library  a tibble containing columns named "peptide",
#' "library_match_sequence", and "library_real_sequence".
#' @param n_residues the number of residues to the left and right of the
#'  cleavage site to include in the output.
#' @return a tibble with cleavage information added.
#' @keywords internal
add_cleavages <- function(joined_with_library, n_residues = 4) {
  # Iterating through and applying nterm_clevage
  nterm <- purrr::pmap_df(
    list(
      joined_with_library$peptide,
      joined_with_library$library_match_sequence,
      joined_with_library$library_real_sequence,
      n_residues
    ),
    nterm_cleavage
  )
  # Iterating though and applying cterm_cleavage
  cterm <- purrr::pmap_df(
    list(
      joined_with_library$peptide,
      joined_with_library$library_match_sequence,
      joined_with_library$library_real_sequence,
      n_residues
    ),
    cterm_cleavage
  )
  # Combining nterm and cterm
  cleavages <- dplyr::bind_cols(nterm, cterm[, 2:3])
  joined_with_library <- dplyr::select(joined_with_library, -"peptide")
  # Building final data frame.
  output <- dplyr::bind_cols(cleavages, joined_with_library)
  return(output)
}

#' consolidate_cleavages
#'
#' Consolidate the n term and c term cleavage data. The nterm and cterm
#' cleavage information  are consolidated into a single column and rows
#  that have both nterm and cterm cleavage information are removed.
#'
#' @param cleavage_added_data a tibble where cleavage information has
#' been added by add_cleavages()
#'
#' @return a tibble with the cleavage information combined into a single column
#'  and rows with no cleavage information or double information removed.
#' @keywords internal
consolidate_cleavages <- function(cleavage_added_data) {
  out <- cleavage_added_data %>%
    # consolidating cleavage sequence
    dplyr::mutate(cleavage_seq = dplyr::case_when(
      !is.na(.data$nterm) & is.na(.data$cterm) ~ .data$nterm,
      !is.na(.data$cterm) & is.na(.data$nterm) ~ .data$cterm,
      TRUE ~ NA_character_
    ), .after = "cterm_cleavage_pos") %>%
    # Removing peptides with double cleavages
    dplyr::filter(!(!is.na(.data$cterm) & !is.na(.data$nterm))) %>%
    dplyr::mutate(cleavage_pos = dplyr::case_when(
      is.na(.data$cterm_cleavage_pos) ~ .data$nterm_cleavage_pos,
      TRUE ~ .data$cterm_cleavage_pos
    ), .after = "cleavage_seq") %>%
    # Adding character specifying what are cleaved vs not
    dplyr::mutate(
      peptide_type = dplyr::case_when(
        is.na(.data$cleavage_pos) ~ "full_length",
        TRUE ~ "cleavage_product"
      ),
      .after = "cleavage_pos"
    ) %>%
    dplyr::select(
      -"nterm", -"cterm", -"nterm_cleavage_pos",
      -"cterm_cleavage_pos"
    )
  return(out)
}
#' convert prepared data to a QFeatures object
#'
#' @param prepared_data data prepared within one of the prepare functions
#' @param colData  sample metadata
#' @param peptide_library the peptide library used.
#' @param n_residues the number of residues reported in the cleavage site
#'
#' @return a QFeatures object
#' @keywords internal
prepared_to_qf <- function(prepared_data,
                           colData,
                           peptide_library = mspms::peptide_library,
                           n_residues = 4) {
  # adding peptide length to prepared_data
  prepared_data <- prepared_data %>%
    dplyr::mutate(
      peptide_length = gsub("^._", "", peptide),
      peptide_length = nchar(gsub("_.$", "", peptide_length)),
      .after = peptide
    )
  # check peptide library for correct names.
  check_peptide_library(peptide_library)
  # Making sure that prepared_data and peptide library are consistent
  peptide_library_ids <- peptide_library$library_id
  peptide_library_ids_data <- unique(prepared_data$library_id)
  missing <- peptide_library_ids_data[peptide_library_ids_data %!in%
    peptide_library_ids]
  missing_mes <- paste0(missing, collapse = ",")
  if (length(missing > 1)) {
    stop("There are peptide library ids in your data that are not in your
         peptide library. Specificially ", missing_mes, "are missing from your
         peptide library.")
  }
  # combining peptide sequences
  combined <- dplyr::inner_join(peptide_library, prepared_data,
    by = "library_id"
  ) %>%
    add_cleavages(n_residues = n_residues) %>%
    consolidate_cleavages()
  name_in_prepared_data <- names(prepared_data)[3:length(prepared_data)]
  n_coldata_nin_prepared_data <- sum(colData$quantCols %!in%
    name_in_prepared_data)
  missing_name <- paste0(
    name_in_prepared_data[
      colData$quantCols %!in% name_in_prepared_data
    ],
    collapse = " "
  )
  if (n_coldata_nin_prepared_data > 0) {
    stop(
      "the quantCol names in your colData do not match those in your",
      " proteomics data. Specifically the column(s) ", missing_name,
      " are present in your proteomics data but not in your",
      " colData"
    )
  }
  # Need columns in the combined data frame to be in the same order as colData
  # Extract the ordered sample names (columns 9 onward in combined)
  order_in_combined <- names(combined)[9:ncol(combined)]

  # arrange colData to be in the same order
  colData <- colData |>
    dplyr::arrange(factor(quantCols, levels = order_in_combined))

  QF <- QFeatures::readQFeatures(combined,
    quantCols = 9:length(combined),
    fnames = "peptide",
    colData = colData,
    name = "peptides"
  )
  return(QF)
}

#' load_colData
#'
#' load a .csv file containing sample colData. Check for errors
#'
#' @param colData_filepath filepath to .csv file containing colData.
#'
#' @return a tibble
#' @keywords internal

load_colData <- function(colData_filepath) {
  colData <- readr::read_csv(colData_filepath)
  # Expected names
  `%!in%` <- Negate(`%in%`)
  expected_names <- c("quantCols", "group", "condition", "time")
  colData_names <- names(colData)
  unexpected_names <- expected_names %!in% colData_names
  n_unexpected_names <- sum(unexpected_names)
  if (n_unexpected_names > 0) {
    missing_names <- expected_names[unexpected_names]
    stop(
      paste0(c("colData must have columns named \"quantCols\",\"group\",
    \"condition\", and \"time\" ", "you are missing ", paste(missing_names,
        collapse = ", "
      ))),
      collapse = ""
    )
  }
  # Check if the column types match
  expected_types <- c("character", "character", "character", "numeric")
  colData_types <- sapply(colData, class)
  unexpected_types <- (colData_types != expected_types[1:4])
  n_unexpected_types <- sum(unexpected_types)

  if (n_unexpected_types > 0) {
    missing_types <- names(colData)[unexpected_types]
    stop(paste0(
      "colData columns must be as follows: quantCols = character,
    group = character, condition = character, time = numeric.",
      "You have the wrong type for the following column(s): ",
      paste(missing_types, collapse = ", ")
    ))
  }
  # Check to make sure none of the columns have spaces in them
  n_spaces <- sum(sapply(colData, function(x) {
    grepl(".* .*", x)
  }))
  if (n_spaces > 0) {
    stop("Your colData has spaces in it (ie condition = Cathepsin A). Please
    change it to not contain any spaces to be compatible with limma (ie
    condition  = cathepsin_A)")
  }
  return(colData)
}

#' check_peptide_library
#'
#' @param peptide_library
#'
#' @return an informative error if the column names of the peptide library are
#' unexpected. Otherwise nothing.
#' @keywords internal

check_peptide_library <- function(peptide_library) {
  pl_names <- names(peptide_library)
  if (!identical(pl_names, c(
    "library_id",
    "library_match_sequence",
    "library_real_sequence"
  ))) {
    stop("the first three columns of the peptide library .csv are not as
         expected. They must be library_id, library_match_sequence, and
         library_real_sequence")
  }
}


#' check_peptide_library
#'
#' @param peptide_library
#'
#' @return an informative error if the column names of the peptide library are
#' unexpected. Otherwise nothing.
#' @keywords internal
