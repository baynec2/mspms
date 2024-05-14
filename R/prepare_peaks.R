#' prepare_peaks
#' This function prepares the data out from peaks for downstream normalization via normalyzer.
#' Briefly, it combines the two files by peptide IDs and filters out peptides that have quality scores <= 0.3. 0s are also replaced with NAs
#' Note that there can be multiple matches between these, currently only keeps the values for the best quality scores. This could be a bug.
#'
#' @param lfq_filepath = this is the filepath to the first PEAKS output table
#' @param id_filepath = this is the filepath to the second PEAKS output table
#'
#' @return a tibble containing the combined columns from the lfq and id files, with quality scores > 0.3 and 0s replaced with NAs.
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' prepare_peaks("tests/testdata/protein-peptides-lfq.csv",
#'              "tests/testdata/protein-peptides-id.csv")
#'
prepare_peaks = function(lfq_filepath,
                         id_filepath){

  # Reading in the label free quantification data
  lfq = readr::read_csv(lfq_filepath,guess_max = 10,na = c("-",""))

  `%!in%` = Negate(`%in%`)

  # Making sure the lfq file has the correct headers
  if(sum(c("Protein Group",
           "Protein ID",
           "Protein Accession",
           "Peptide",
           "Used",
           "Candidate",
           "Quality",
           "Significance",
           "Avg. ppm",
           "Avg. Area") %!in% names(lfq))>0){
    stop(paste("The lfq file does not have the correct headers.
                The headers should contain:", "Protein Group","Protein ID","Protein Accession","Peptide",
                "Used","Candidate","Quality","Significance","Avg.ppm","Avg.Area"))

  }

  # Reading in the ids
  id = readr::read_csv(id_filepath,na = c("-",""))

  #Making sure id file has the correct headers
  if(sum(c("Protein Group",
           "Protein ID",
           "Protein Accession",
           "Peptide",
           "Unique",
           "-10lgP",
           "Mass",
           "Length",
           "ppm",
           "m/z",
           "z",
           "RT") %!in%
    names(id))>0){
    stop(paste("The id file does not have the correct headers.
                The headers should contain:","Protein Group","Protein ID","Protein Accession","Peptide",
                "Unique","-10lgP","Mass","Length","ppm","m/z","z","RT", "and columns corresponding to your sample names"))

  }

  id = id %>%
    #only keep the Peptide ID with the highest score in case there are more than one.
    # This is essentially what the current script does, is it intended? Seems like a bug to me
    dplyr::group_by(.data$Peptide) %>%
    dplyr::filter(.data$`-10lgP` == max(.data$`-10lgP`)) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$Peptide,.data$RT,6:12)

  # Combining data frame and filtering to only contain quality scores > 0.3
  output = lfq %>%
    # Sometimes a peptide is detected more than once at different retention times.
    # We will only keep the one with the highest quality score.
    dplyr::group_by(.data$Peptide) %>%
    dplyr::filter(.data$Quality == max(.data$Quality)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(id,by = c("Peptide"),multiple = "first") %>%
    dplyr::filter(.data$Quality > 0.3) %>%
    tibble::as_tibble()


  # finding the columns with our samples

  start = which(names(output) == "Avg. Area") + 1
  end = which(names(output) == "Sample Profile (Ratio)") - 1

  # Selecting only the columns that we care about
  output = output %>%
    dplyr::mutate(Peptide = gsub("\\.","_",.data$Peptide)) %>%
    dplyr::select(.data$Peptide,.data$RT,.data$`Protein Accession`,dplyr::any_of(start:end)) %>%
    #Dealing with the case where there are PTMs, removing the mod from peptide seq report.
    dplyr::mutate(Peptide = gsub("//(.*\\)","",Peptide))

  # Replacing 0 with NA
  output[output == 0] = NA_real_


  return(output)
}


