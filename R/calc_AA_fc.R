#'calc_AA_fc
#'
#'calculate the fold change of each AA by position.
#'
#' @param experimental_prop_matrix this is a matrix of the experimental proportions (from your vector of cleavage sequences) at each position.
#' @param background_prop_matrix this is a matrix of the background proportions of AAs at each position
#' @return a data frame.
#' @export
#' @examples
#' experimental_prop_matrix = matrix(c(0.4,0.2),nrow = 1,dimnames = list("A"))
#' background_prop_matrix = matrix(c(0.2,0.4),nrow = 1,dimnames = list("A"))
#' calc_AA_fc(experimental_prop_matrix, background_prop_matrix)
#'
#'
#'
#'
calc_AA_fc= function(experimental_prop_matrix,
                     background_prop_matrix){

  fold_change = (experimental_prop_matrix * 100) / (background_prop_matrix * 100)

  converted_fc = fold_change

  # Here we run into a problem with infinite values. If a value is 0 in the count, but greater than that in the background - we have a problem.
  # To deal with this, ice logo just shows an infinite value as taking up the entire scale.
  # We could mimic that here by setting the columns with infinite values to the max value of the fold change/ the number of things in the columnns

  max_sum = max(purrr::map_df(as.data.frame(converted_fc), sum), na.rm = T)
  final_converted_fc  = data.frame(row.names = rownames(converted_fc))

  for (i in seq_len(ncol(converted_fc))) {
    if (0 %in% converted_fc[, i]) {
      num_zero = sum(converted_fc[, i] == 0, na.rm = TRUE)

      col = replace(converted_fc[, i],
                    converted_fc[, i] == 0,
                    (max_sum / num_zero * -1))

      final_converted_fc =  cbind(final_converted_fc, col)
    } else{
      final_converted_fc = cbind(final_converted_fc, converted_fc[, i])
    }
  }

  colnames(final_converted_fc) = colnames(converted_fc)

  return(final_converted_fc)
}
