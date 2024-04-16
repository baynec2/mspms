#' polish
#'
#' This function is used to polish the cleavage data. It combines the nterm and cterm cleavage information into a single column and removes rows that
#' don't have any cleavage information or have both nterm and cterm cleavage information.
#'
#' @param cleavage_added_data
#'
#' @return a data frame with the cleavage information combined into a single column and rows with no cleavage information or double information removed.
#' @export
#'
#' @examples
polish = function(cleavage_added_data){

  out = cleavage_added_data %>%
    dplyr::mutate(cleavage_seq = dplyr::case_when(!is.na(nterm) & is.na(cterm) ~ nterm,
                                                  !is.na(cterm) & is.na(nterm) ~ cterm,
                                                  TRUE ~NA),.after = "cterm_cleavage_pos") %>%
    dplyr::filter(!is.na(cleavage_seq)) %>%
    tibble::as.tibble()


  return(out)
}
