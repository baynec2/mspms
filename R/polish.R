#' polish
#'
#' This function is used to polish the cleavage data. It combines the nterm and
#' cterm cleavage information into a single column and removes rows that
#' don't have any cleavage information or have both nterm and cterm cleavage
#' information.
#'
#' @param cleavage_added_data = this is a tibble where cleavage information has
#' been added by add_cleavages()
#'
#' @return a tibble with the cleavage information combined into a single column
#'  and rows with no cleavage information or double information removed.
#' @export
#' @examples
#' polished <- polish(mspms::cleavage_added_data)
polish <- function(cleavage_added_data) {
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
    tibble::as_tibble()
  
  return(out)
}
