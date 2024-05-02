#' prepare_pd
#'
#' prepare proteome discoverer output for mspms analysis.
#'
#' @param filepath = this is the filepath to the proteome discoverer excel formatted file.
#'
#' @return a tibble with the data formatted for use with normalyze
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#' prepared_proteome_discoverer = prepare_pd("tests/testdata/proteome_discoverer_output.xlsx")

prepare_pd = function(filepath){

  # Read in the file
  data = readxl::read_excel(filepath)


  # Throw error if this doesn't look like a proteome discoverer file

  if(!all(c("Annotated Sequence",
            "Master Protein Accessions",
            "Qvality PEP") %in% names(data))){
    stop("This does not look like a proteome discoverer file.
         The file should contain columns named: Annotated Sequence, Master Protein Accessions, Qvality PEP")
  }

  #Make the names consistent with what is used in the package.

  name_fixed = data %>%
    dplyr::rename(Peptide = .data$`Annotated Sequence`,
                  `Protein Accession` = .data$`Master Protein Accessions`)

  # Converting the peptide notation from how it is expressed in proteome discover to what we have used.
  peptide_fixed = name_fixed %>%
    dplyr::mutate(Peptide = gsub("\\[-\\]","",.data$Peptide)) %>%
    dplyr::mutate(Peptide = gsub("^\\.","",.data$Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\.$","",.data$Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\].","_",.data$Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\.\\[","_",.data$Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\]","",.data$Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\[","",.data$Peptide)) %>%
    #Only keep max quality peptide, for cases where there is more than one match.
    dplyr::group_by(.data$Peptide) %>%
    dplyr::filter(.data$`Qvality PEP` == max(.data$`Qvality PEP`)) %>%
    dplyr::ungroup()

  #Converting the column names to reasonable sample ID names
  sample_name_fixed = peptide_fixed %>%
    dplyr::rename_with(~gsub("Abundance: .*: ","",.))


  #Let's only keep the columns that contain the data we are interested in.

  sample_cols = which(grepl("Abundance.*",names(data)))

  out = sample_name_fixed %>%
    dplyr::select(.data$Peptide,.data$`Protein Accession`,dplyr::all_of(sample_cols))


 return(sample_name_fixed)

}

