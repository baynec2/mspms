# Negate
`%!in%` = Negate(`%in%`)

#' check_file_is_valid_peaks
#' check to make sure a loaded tibble is what a peaks file is expected to look
#' like.
#' @param peaks_data 
#'
#' @return
#' @keywords internal
check_file_is_valid_peaks = function(peaks_data){
  # These are the names of the columns expected to be in every PEAKS file
  expected_names <- c(
    "Protein Group", "Protein ID", "Protein Accession", "Peptide",
    "Used", "Candidate", "Quality", "Significance", "Avg. ppm", "Avg. Area",
    "Sample Profile (Ratio)", "Group Profile (Ratio)", "Max Ratio", "#Vector",
    "Start", "End", "PTM"
  )
  names_in_data = names(peaks_data)
  # What names are missing
  missing_names = expected_names[expected_names %!in% names_in_data]
  if(length(missing_names) > 1){
    stop(
      paste("This doesn't look like the expected PEAKS file. The following
          columns are missing:",missing_names)
    )
  }
  # PEAKS files with PTMS are not supported
  if(sum(!is.na(peaks_data$PTM)) > 1){
    stop("There are PTMS in your file, redo PEAKS search without
    PTMs. PTM analyis is not currently supported")
  }
  # Add something to check the colData to make sure it matches the files
}

check_file_is_valid_fragpipe = function(fragpipe_data){
  # These are the names of the columns expected to be in every fragpipe file
  expected_names <- c(
    "Peptide Sequence", "Prev AA", "Next AA", "Start", "End", "Peptide Length",
    "Charges", "Protein", "Protein ID", "Entry Name", "Gene",
    "Protein Description","Mapped Genes", "Mapped Proteins",
    "CatA_0000_1 Spectral Count"
  )
  
  names_in_data = names(fragpipe_data)
  
  # What names are missing
  missing_names = expected_names[expected_names %!in% names_in_data]
  if(length(missing_names) > 1){
    stop(
      paste("This doesn't look like the expected fragpipe file. The following
          columns are missing:",missing_names)
    )
  }
  # Add something that looks at the colData to make sure it matches 
  
}

#' cterm_cleavage
#'
#' Finding the cleavage sequences on the C terminus of a given peptide in
#' reference to the peptide library it was derived from
#'
#' @param peptide_sequence = this is the peptide sequence in single letter
#' AA code. _ denotes cleavage site.
#' @param library_match_sequence  = this is the sequence of the corresponding
#' peptide in the library to the peptide_sequence.
#' @param library_real_sequence = this is the library match sequence with the
#' Ls transformed to Ms.
#' @param n_residues = this is the number of residues to the left and right of
#' the cleavage site that you want to extract.

#' @return a tibble with the peptide sequence, cleavage sequences n number of
#' AA on the left and right of the c term cleavage, and the position of the
#' c-term cleavage in the library sequence.
#' @keywords internal
#' @examples
#' peptide_sequence <- "YWLSTHLAGKR_R"
#' library_match_sequence <- "YWLSTHLAGKRRDW"
#' library_real_sequence <- "YWMSTHLAGKRRDW"
#' cterm_cleavage(
#'   peptide_sequence, library_match_sequence,
#'   library_real_sequence
#' )
#'
cterm_cleavage <- function(peptide_sequence,
                           library_match_sequence,
                           library_real_sequence,
                           n_residues = 4) {
  # First, we determine if there is a cterm cleavage
  # We do this by looking to see if there is a _ (denotes a cleavage)
  # If there is, we check to see if it is in the second to last position
  # Or there might be 2 _s in the peptide sequence
  # (meaning there is a cterm and nterm cleavage)
  
  n_of_match <- length(gregexpr("_", peptide_sequence)[[1]])
  if ((grepl("_", peptide_sequence) == TRUE) &&
      (gregexpr("_", peptide_sequence)[[1]][n_of_match] == nchar(peptide_sequence) - 1)) {
    # if there is a c term cleavage, it is the last position - 1 e
    pos <- nchar(peptide_sequence) - 1
    
    # now we define the sequence on the left side of the cleavage
    temp <- substr(peptide_sequence, pos - n_residues, pos - 1)
    
    # Checking to see what part of the reference sequence this matches.
    left_reference_beginning <- regexpr(temp, library_match_sequence)[[1]][1]
    # it is four AA long, so the end is the beginning + 3
    left_reference_end <- left_reference_beginning + (n_residues - 1)
    
    # now we can determine the right side of the cleavage sequence
    right_reference_beginning <- left_reference_end + 1
    
    # the cterm cleavage position is the right reference beginning minus 1
    cterm_cleavage_pos <- right_reference_beginning - 1
    
    # it is four AA long, so the end is the beginning + 3
    right_reference_end <- right_reference_beginning + (n_residues - 1)
    
    # Extracting the sequences from the reference sequence
    right_sequence <- substr(
      library_real_sequence,
      right_reference_beginning,
      right_reference_end
    )
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    right_sequence <- paste0(
      right_sequence,
      paste0(
        rep(
          "X",
          n_residues - nchar(right_sequence)
        ),
        collapse = ""
      )
    )
    
    # Doing the same on the left side of the cleavage event
    
    left_sequence <- substr(
      library_real_sequence,
      left_reference_beginning,
      left_reference_end
    )
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    left_sequence <- paste0(paste0(rep("X", n_residues - nchar(left_sequence)),
                                   collapse = ""
    ), left_sequence)
    
    
    cterm <- paste(c(left_sequence, right_sequence), collapse = "")
  } else {
    cterm <- NA
    cterm_cleavage_pos <- NA
  }
  output <- tibble::tibble(
    peptide = peptide_sequence,
    cterm = cterm,
    cterm_cleavage_pos = cterm_cleavage_pos
  )
  return(output)
}



#' nterm_cleavage
#'
#' Finding the cleavage sequences on the N terminus of a given peptide in
#' reference to the peptide library it was derived from.
#'
#' @param peptide_sequence = this is the peptide sequence in single leter AA
#' code. _ denotes cleavage site..
#' @param library_match_sequence  = this is the sequence of the corresponding
#'  peptide in the library to the peptide_sequence.
#' @param library_real_sequence = this is the library match sequence with the
#' Ls transformed to Ms (This is what the legacy code did so it is kept this way
#'  in case there was a good reason for it)
#' @param n_residues = this is the number of residues to the left and right of
#' the cleavage site that you want to extract.

#' @return a tibble with the peptide sequence, cleavage sequences n specified
#'  number of AA on the left and right of the n term cleavage, and the position
#'   of the n term cleavage in the library sequence.
#' @keywords internal
#' @examples
#'
#' peptide_sequence <- "F_LAHWVGI"
#' library_real_sequence <- "GMYYKRFMAHWVGI"
#' library_match_sequence <- "GLYYKRFLAHWVGI"
#'
#'
#' nterm_cleavage(peptide_sequence, library_real_sequence, library_match_sequence)
nterm_cleavage <- function(peptide_sequence,
                           library_match_sequence,
                           library_real_sequence,
                           n_residues = 4) {
  # Determining if there is a n term cleavage
  # _ denotes a cleavage, and if it is the second position, it is on the n term!
  if ((grepl("_", peptide_sequence) == TRUE) &&
      (regexpr("_", peptide_sequence)[[1]][1] == 2)) {
    # if there is a n term cleavage, the first letter of the right side is the third letter our sequence according to the logic above.
    pos <- 2 + 1
    
    # taking the sequence from right after the _ to .
    temp <- substr(peptide_sequence, pos, pos + (n_residues - 1))
    
    # Checking to see what part of the reference sequence this matches.
    right_reference_start <- regexpr(temp, library_match_sequence)[[1]][1]
    right_reference_end <- right_reference_start + (n_residues - 1)
    
    # the position of the n term cleavage is right reference start - 1
    nterm_cleavage_pos <- right_reference_start - 1
    
    # Now determining the left side of the cleavage event.
    left_reference_end <- right_reference_start - 1
    left_reference_start <- left_reference_end - (n_residues - 1)
    
    # Extracting the sequences from the reference sequence
    right_sequence <- substr(
      library_real_sequence,
      right_reference_start,
      right_reference_end
    )
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    right_sequence <- paste0(
      right_sequence,
      paste0(rep("X", n_residues - nchar(right_sequence),
                 collapse = ""
      ))
    )
    
    left_sequence <- substr(
      library_real_sequence,
      left_reference_start,
      left_reference_end
    )
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    left_sequence <- paste0(
      paste0(rep("X", n_residues - nchar(left_sequence)),
             collapse = ""
      ),
      left_sequence
    )
    
    
    nterm <- paste(c(left_sequence, right_sequence), collapse = "")
  } else {
    nterm <- NA
    nterm_cleavage_pos <- NA
  }
  output <- tibble::tibble(
    peptide = peptide_sequence,
    nterm = nterm,
    nterm_cleavage_pos = nterm_cleavage_pos
  )
  return(output)
}

#' add_cleavages
#'
#' This function adds cleavage information to the tibble that has been
#' normalized, outlier removed, imputed, and joined with (peptide) library.
#' The cleavage information is represented as the 4 amino acids to the left and
#' right of the cleavage site.
#' If there is no amino acid in a position of original peptide in the library
#' that was cleaved, it is represented as an X.
#' This wraps the mspms::n_term_cleavage and mspms::c_term_cleavage functions
#' into a consolidated function.
#'
#' @param joined_with_library = this is the tibble that has been normalized,
#'  outlier removed, imputed, and joined with the library.
#' @param n_residues = the number of residues to the left and right of the
#'  cleavage site to include in the output.
#' @return a tibble with cleavage information added.
#' @keywords internal
#' @examples
#' # adding the clevages 5 AA to the left and right of the cleavage site.
#' add_cleavages(mspms::joined_with_library, n_residues = 5)
#'
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
  joined_with_library = dplyr::select(joined_with_library,-"peptide")
  # Building final data frame.
  output <- dplyr::bind_cols(cleavages,joined_with_library)
  return(output)
}

#' consolidate_cleavages
#'
#' This function is used to consolidatae the cleavage data. It combines the 
#' nterm and cterm cleavage information into a single column and removes rows 
#' that don't have any cleavage information or have both nterm and cterm
#' cleavage information.
#'
#' @param cleavage_added_data = this is a tibble where cleavage information has
#' been added by add_cleavages()
#'
#' @return a tibble with the cleavage information combined into a single column
#'  and rows with no cleavage information or double information removed.
#' @export
#' @examples
#' polished <- polish(mspms::cleavage_added_data)
consolidate_cleavages <- function(cleavage_added_data) {
  out <- cleavage_added_data %>%
    # consolidating cleavage sequence
    dplyr::mutate(cleavage_seq = dplyr::case_when(
      !is.na(nterm) & is.na(cterm) ~ nterm,
      !is.na(cterm) & is.na(nterm) ~ cterm,
      TRUE ~ NA
    ), .after = "cterm_cleavage_pos") %>%
    # Removing peptides with double cleavages
    dplyr::filter(!(!is.na(.data$cterm) & !is.na(.data$nterm))) %>%
    dplyr::mutate(cleavage_pos = dplyr::case_when(
      is.na(cterm_cleavage_pos) ~ nterm_cleavage_pos,
      TRUE ~ cterm_cleavage_pos
    ), .after = "cleavage_seq") %>%
    # Adding character specifying what are cleaved vs not
    dplyr::mutate(peptide_type = dplyr::case_when(is.na(cleavage_pos) ~ "full_length",
                                                  TRUE ~ "cleavage_product"),
                  .after = cleavage_pos) %>% 
    dplyr::select(-"nterm",-"cterm",-"nterm_cleavage_pos",-"cterm_cleavage_pos") 
  
  return(out)
}
#' convert prepared data to a QFeatures object
#'
#' @param prepared_data = data prepared by one of the prepare functions
#' @param colData = metadata
#' @param peptide_library = the peptide library used. 
#' @param n_residues = the number of residues reported in the cleavage site
#'
#' @return
#' @keywords internal
#'
#' @examples
prepared_to_qf = function(prepared_data,
                          colData,
                          peptide_library = mspms::peptide_library,
                          n_residues = 4){
  
  # combining peptide sequences
  combined = dplyr::inner_join(peptide_library,prepared_data,
                               by = "library_id") %>% 
    add_cleavages(n_residues = n_residues) %>% 
    consolidate_cleavages()
  
  # Creating a QFeatures object
  QF = QFeatures::readQFeatures(combined,
                                quantCols = 8:length(combined),
                                fnames = "peptide",
                                colData = colData,
                                name = "peptides"
  )
  
  return(QF)
}