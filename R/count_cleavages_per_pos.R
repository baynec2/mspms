#' count_cleavages_per_pos
#'
#' count the number of cleavages per position
#'
#' @param data a tibble containing columns named Peptide,cleavage_pos,condition,
#' and time. Other column names can be included.
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'data = data.frame(condition = rep("DMSO", 13),
#'                  time = rep(60,13),
#'                  cleavage_pos = rep(13,13)
#'                  )
#'count_cleavages_per_pos(data)

count_cleavages_per_pos <- function(data) {
  positions <- seq_len(13)
  #Counting
  count <- data %>%
    dplyr::group_by(.data$condition,.data$time,.data$cleavage_pos) %>%
    dplyr::summarise(n = dplyr::n())
  #Adding missing positions
  final = data.frame()
  for(i in unique(data$condition)){
    f = count %>% dplyr::filter(.data$condition == i)
    for(t in unique(f$time)){
    count_f <- dplyr::filter(f,
                            .data$condition == i,
                            .data$time == t)
    missing <- positions[!positions %in% count_f$cleavage_pos]
    missing_t <- tibble::tibble(condition = i, time = t, cleavage_pos = missing, n = rep(
      0,
      length(missing)
    ))
    out <- dplyr::bind_rows(count_f, missing_t)
    final <- dplyr::bind_rows(final,out) %>%
      dplyr::mutate(time = forcats::fct_inseq(as.character(.data$time)))
    }
  }

  return(final)
}



