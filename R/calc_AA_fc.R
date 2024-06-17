#' calc_AA_fc
#'
#' calculate the fold change of each AA by position.
#'
#' @param experimental_prop_matrix this is a matrix of the experimental
#' proportions (from your vector of cleavage sequences) at each position.
#' @param background_prop_matrix this is a matrix of the background proportions
#'  of AAs at each position
#' @param sig_zscores this is a tibble of the significant zscores.
#' @return a matrix
#' @export
#' @examplesIf isTRUE(FALSE)
#' 
calc_AA_fc <- function(experimental_prop_matrix,
                       background_prop_matrix,
                       sig_zscores) {
  #calculating the fold change of each amino acid by position
  fold_change <- (experimental_prop_matrix * 100) /
    (background_prop_matrix * 100)
  
  # Filtering to only include significant z-scores
  prepared_fc = mspms::prepare_fc(as.data.frame(fold_change),
                                  sig_zscores)

  converted_fc <- prepared_fc
  # Here we run into a problem with infinite values.
  # If a value is 0 in the count, but greater than that in the background
  # Icelogo shows infinite value as taking up the entire scale.
  # We will mimic that here

  max_sum <- max(purrr::map_df(as.data.frame(converted_fc), sum,na.rm = TRUE),
                 na.rm = TRUE)
  final_converted_fc <- data.frame(row.names = rownames(converted_fc))

  for (i in seq_len(ncol(converted_fc))) {
    if (0 %in% converted_fc[, i]) {
      num_zero <- sum(converted_fc[, i] == 0, na.rm = TRUE)

      col <- replace(
        converted_fc[, i],
        converted_fc[, i] == 0,
        (max_sum / num_zero * -1)
      )

      final_converted_fc <- cbind(final_converted_fc, col)
    } else {
      final_converted_fc <- cbind(final_converted_fc, converted_fc[, i])
    }
  }
  
  colnames(final_converted_fc) <- colnames(converted_fc)

  # now converting values less than 1 to converted fold change. 
 final_converted_fc2 = final_converted_fc %>% 
   dplyr::mutate(
     dplyr::across(dplyr::all_of(names(final_converted_fc)),
          ~case_when(.x < 1 & .x > 0  ~ (1/.x) * -1,
                      .default = .x)
          )
   )
 
 final_converted_fc2 = as.matrix(final_converted_fc2)

  return(final_converted_fc2)
}
