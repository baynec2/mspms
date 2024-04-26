#' prepare_pd
#'
#' prepare proteome discoverer output for mspms analysis.
#'
#' @param filepath = this is the filepath to the proteome discoverer excel formatted file.
#'
#' @return
#' @export
#'
#' @examples
#'
#'
prepare_pd = function(filepath){

  # Read in the file
  data = readxl::read_excel(filepath)

  #Make the names consistent with what is used in the package.

  name_fixed = data %>%
    dplyr::rename(Peptide = `Annotated Sequence`,
                  `Protein Accession` = `Master Protein Accessions`)

  # Converting the peptide notation from how it is expressed in proteome discover to what we have used.
  peptide_fixed = name_fixed %>%
    dplyr::mutate(Peptide = gsub("\\[-\\]","",Peptide)) %>%
    dplyr::mutate(Peptide = gsub("^\\.","",Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\.$","",Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\].","_",Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\.\\[","_",Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\]","",Peptide)) %>%
    dplyr::mutate(Peptide = gsub("\\[","",Peptide))

  #Converting the column names to reasonable sample ID names
  sample_name_fixed = peptide_fixed %>%
    dplyr::rename_with(~gsub("Abundance: .*: ","",.))


  #Let's only keep the columns that contain the data we are interested in.

  sample_cols = which(grepl("Abundance.*",names(data)))

  out = sample_name_fixed %>%
    dplyr::select(Peptide,`Protein Accession`,dplyr::all_of(sample_cols))


 return(sample_name_fixed)

}
