#' prepare_sig_pdif
#'
#' prepare significant percent difference data frame for ice logo
#'
#' @param percent_difference = data frame containing the percent differences
#' @param sig_zscores = matrix of significant AAs at each position based on z-scores
#'
#' @return a tibble
#' @export
#' @examples
prepare_sig_p_dif = function(percent_difference,sig_zscores){

  sig_final = percent_difference %>%
    tibble::rownames_to_column("AA") %>%
    tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "percent_difference") %>%
    dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position)) %>%
    dplyr::filter(.data$aa_position %in% sig_zscores$aa_position) %>%
    dplyr::select(-.data$aa_position) %>%
    tidyr::pivot_wider(
      names_from = .data$position,
      values_from = .data$percent_difference,
      names_sort = FALSE
    )

  # Dealing with the problem introduced by the pivot wider function
  missing_cols = ncol(percent_difference) - (ncol(sig_final)-1)
  missing_data = as.data.frame(matrix(nrow = nrow(sig_final), ncol = missing_cols, NA))
  `%!in%` = Negate(`%in%`)
  missing_names = names(percent_difference)[names(percent_difference) %!in% names(sig_final)]
  names(missing_data) = missing_names

  final = sig_final %>%
    dplyr::bind_cols(missing_data) %>%
    dplyr::relocate(paste0("P", c((ncol(percent_difference) / 2):1,
                                  paste0(
                                    1:(ncol(percent_difference) / 2), "'"
                                  )))) %>%
    tibble::column_to_rownames("AA") %>%
    as.matrix()

  return(final)

}
