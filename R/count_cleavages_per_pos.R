#' count_cleavages_per_pos
#'
#' count the number of cleavages per position
#'
#' @param mspms_data = this is the data that has gone through the mspms pipeline.
#' @param sig_peptides = this is a list of peptides of interest.
#'
#' @return a tibble with the number of cleavages per position
#' @export
#'
#' @examples
#'# Extracting the significant peptides for DMSO
#'sig_peptides = mspms::mspms_anova(mspms::mspms_data) %>%
#'  dplyr::filter(p.adj <= 0.05,
#'                condition == "MZB") %>%
#'  dplyr::pull(Peptide)
#'
#'# Counting the number of cleavages per position
#'counts = count_cleavages_per_pos(mspms::mspms_data,sig_peptides)
#'
#'#Plotting the results
#'p1 = counts %>%
#'  ggplot2::ggplot(ggplot2::aes(x = cleavage_pos,y = n)) +
#'  ggplot2::geom_point()+
#'  ggplot2::geom_line()
#'p1
#'
count_cleavages_per_pos = function(mspms_data,sig_peptides){

  positions = seq_len(13)

  count = mspms_data %>%
    dplyr::select(Peptide,cleavage_pos) %>%
    dplyr::distinct() %>%
    dplyr::filter(Peptide %in% sig_peptides) %>%
    dplyr::group_by(cleavage_pos) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup()


  missing = positions[!positions %in% count$cleavage_pos]

  missing_t = tibble::tibble(cleavage_pos = missing, n = rep(0,length(missing)))

  out = dplyr::bind_rows(count,missing_t) %>%
    dplyr::arrange(n)

  return(out)

}


