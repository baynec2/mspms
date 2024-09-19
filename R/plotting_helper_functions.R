#' icelogo_col_scheme
#'
#' @return
#' @keywords internal
#'
#' @examples
icelogo_col_scheme = function(){
col_scheme =  ggseqlogo::make_col_scheme(
    chars = c(
      "G", "S", "T", "Y", "C", "N", "Q", "K", "R", "H", "D", "E", "P", "A", "W",
      "F", "L", "I", "M", "V", "n", "X"
    ),
    group = c(
      rep("Polar", 5), rep("Neutral", 2), rep("Basic", 3), rep("Acidic", 2),
      rep("Hydrophobic", 9),
      "Past Terminus"
    ),
    col = c(
      rep("#058644", 5), rep("#720091", 2), rep("#0046C5", 3), rep("#C5003E", 2),
      rep("#2E2E2E", 9), "#808080"
    ),
    name = "chemistry"
  )
return(col_scheme)
}

#' count_cleavages_per_pos
#'
#' count the number of cleavages per position
#'
#' @param data a tibble containing columns named Peptide,cleavage_pos,condition,
#' and time. Other column names can be included.
#' @return a ggplot2 object
#' @keywords internal
#'
#' @examples
#'data = data.frame(condition = rep("DMSO", 13),
#'                  time = rep(60,13),
#'                  cleavage_pos = rep(13,13)
#'                  )
#'count_cleavages_per_pos(data)

count_cleavages_per_pos <- function(data,
                                    peptide_library = mspms::peptide_library) {
  
  peptide_length = unique(nchar(peptide_library$library_match_sequence))
  
  positions <- seq_len(peptide_length-1)
  #Counting
  count <- data %>%
    dplyr::group_by(.data$condition,.data$time,.data$cleavage_pos) %>%
    dplyr::summarise(n = dplyr::n())
  #Adding missing positions
  final = tibble::tibble()
  for(i in unique(data$condition)){
    f = count %>% dplyr::filter(.data$condition == i)
    for(t in unique(f$time)){
      count_f <- dplyr::filter(f,
                               .data$condition == i,
                               .data$time == t)
      missing <- positions[!positions %in% count_f$cleavage_pos]
      missing_t <- tibble::tibble(condition = i, time = t,
                                  cleavage_pos = missing, n = rep(
                                    0,
                                    length(missing)
                                  ))
      out <- dplyr::bind_rows(count_f, missing_t)
      final <- dplyr::bind_rows(final,out) 
    }
  }
  final$time = as.character(final$time)
  
  return(final)
}



