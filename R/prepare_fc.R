#' prepare_fc
#'
#' prepare fold changes of AA by position for Icelogo visualization.
#'
#' @param fold_change = this is a matrix of the fold changes of the AA by
#' position.
#' @param sig_zscores = this is a tibble of the significant zscores.
#'
#' @return a matrix of the fold changes of the significant AAs at each position.
#' @export
#' @examples
prepare_fc <- function(fold_change, sig_zscores) {
  sig_final <- fold_change %>%
    tibble::rownames_to_column("AA") %>%
    tidyr::pivot_longer(2:length(.),
      names_to = "position",
      values_to = "fold_change"
    ) %>%
    dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position)) %>%
    dplyr::filter(.data$aa_position %in% sig_zscores$aa_position) %>%
    dplyr::select(-.data$aa_position) %>%
    tidyr::pivot_wider(
      names_from = .data$position,
      values_from = .data$fold_change,
      names_sort = FALSE
    )

  # Dealing with the problem introduced by the pivot wider function
  missing_cols <- ncol(percent_difference) - (ncol(sig_final) - 1)
  missing_data <- as.data.frame(matrix(
    nrow = nrow(sig_final),
    ncol = missing_cols, NA
  ))
  `%!in%` <- Negate(`%in%`)
  missing_names <- names(pd)[names(pd) %!in% names(sig_final)]
  names(missing_data) <- missing_names

  final <- sig_final %>%
    dplyr::bind_cols(missing_data) %>%
    dplyr::relocate(paste0("P", c(
      (ncol(percent_difference) / 2):1,
      paste0(
        seq_len(ncol(percent_difference) / 2), "'"
      )
    ))) %>%
    tibble::column_to_rownames("AA") %>%
    as.matrix()


  return(final)
}
